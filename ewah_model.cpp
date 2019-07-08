#include <cstring> //memcpy
#include "ewah_model.h"
#include "compressors.h"

namespace djinn {

/*======   EWAH container   ======*/

djn_ewah_model_t::djn_ewah_model_t() :
    p(nullptr),
    p_len(0), u_len(0), p_cap(0), p_free(false),
    n_variants(0)
{
    
}

djn_ewah_model_t::~djn_ewah_model_t() {
    if (p_free) delete[] p;
}

void djn_ewah_model_t::reset() {
    pbwt.Reset();
    n_variants = 0;
}

int djn_ewah_model_t::StartEncoding(bool use_pbwt, bool reset) {
    if (p_cap == 0) { // initiate a buffer if there is none
        delete[] p;
        p = new uint8_t[10000000];
        p_len = 0;
        p_cap = 10000000;
        p_free = true;
    }
    if (reset) this->reset();
    p_len = 0;
    u_len = 0;

    return 1;
}

size_t djn_ewah_model_t::FinishEncoding(uint8_t* support_buffer, uint32_t support_cap, CompressionStrategy strat) {
    // int ZstdCompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
    GeneralCompressor comp;
    switch(strat) {
        case (CompressionStrategy::ZSTD): comp = &ZstdCompress; break;
        case (CompressionStrategy::LZ4):  comp = &Lz4Compress;  break;
    }
    u_len = p_len;
    int ret = (*comp)(p, p_len, support_buffer, support_cap, 1);
    memcpy(p, support_buffer, ret); // copy data back to p
    p_len = ret;
    // std::cerr << "[djn_ewah_model_t::FinishEncoding] debug=" << u_len << "->" << p_len << "->" << ret << std::endl;
    // return p_len;
    return ret;
}

int djn_ewah_model_t::StartDecoding(uint8_t* support_buffer, uint32_t support_cap, CompressionStrategy strat, bool use_pbwt, bool reset) {
    if (reset) this->reset();
    if (p == nullptr) return -2;

    GeneralDecompressor dec;
    switch(strat) {
        case (CompressionStrategy::ZSTD): dec = &ZstdDecompress; break;
        case (CompressionStrategy::LZ4):  dec = &Lz4Decompress;  break;
    }
    // int ZstdDecompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity)
    // std::cerr << "[model decode]" << p_len << std::endl;
    int ret = (*dec)(p, p_len, support_buffer, support_cap);
    // std::cerr << "[model decode]" << p_len << "->" << ret << " capacity=" << p_cap << std::endl;
    memcpy(p, support_buffer, ret); // copy data back to p
    p_len = 0;
    // std::cerr << "return model decode" << std::endl;

    return 1;
}

/*======   Variant EWAH model   ======*/

djinn_ewah_model::djinn_ewah_model() : 
    codec(CompressionStrategy::ZSTD),
    p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
    q(nullptr), q_len(0), q_alloc(0), q_free(true)
{
    
}

djinn_ewah_model::~djinn_ewah_model() { 
    if (p_free) delete[] p;
    if (q_free) delete[] q;
}

int djinn_ewah_model::EncodeBcf(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) {
    if (data == nullptr) return -2;
    if (len_data % ploidy != 0) return -3;
    // Currently limited to 14 alt alleles + missing + EOV marker (total of 16).
    assert(alt_alleles < 14);

    // Todo: check that ploidy_dict model is initiated!
    std::shared_ptr<djn_ewah_model_container_t> tgt_container;

    const uint64_t tuple = ((uint64_t)len_data << 32) | ploidy;
    auto search = ploidy_map.find(tuple);
    if (search != ploidy_map.end()) {
        tgt_container = ploidy_models[search->second];
        // ploidy_dict->EncodeSymbol(search->second);
        p[p_len++] = search->second;
    } else {
        // std::cerr << "Not found. Inserting: " << len_data << "," << ploidy << "(" << tuple << ") as " << ploidy_models.size() << std::endl;
        ploidy_map[tuple] = ploidy_models.size();
        // ploidy_dict->EncodeSymbol(ploidy_models.size());
        p[p_len++] = ploidy_models.size();

        ploidy_models.push_back(std::make_shared<djn_ewah_model_container_t>(len_data, ploidy, (bool)use_pbwt));
        tgt_container = ploidy_models[ploidy_models.size() - 1];
        tgt_container->StartEncoding(use_pbwt, init);
    }
    assert(tgt_container.get() != nullptr);

    // Compute allele counts.
    memset(hist_alts, 0, 256*sizeof(uint32_t));
    for (int i = 0; i < len_data; ++i) {
        ++hist_alts[BCF_UNPACK_GENOTYPE_GENERAL(data[i])];
    }

    // for (int i = 0; i < 256; ++i) {
    //     if (hist_alts[i]) std::cerr << i << ":" << hist_alts[i] << ",";
    // }
    // std::cerr << " data2m=" << tgt_container->model_2mc->range_coder->OutSize() << " dataNm=" << tgt_container->model_nm->range_coder->OutSize() << " permute=" << permute << std::endl;

    // Biallelic, no missing, and no special EOV symbols.
    if (alt_alleles <= 2 && hist_alts[14] == 0 && hist_alts[15] == 0) {
        // tgt_container->marchetype->EncodeSymbol(0); // add archtype as 2mc
        tgt_container->p[tgt_container->p_len++] = 0;

        int ret = -1;
        if (use_pbwt) {
            // std::cerr << "encoding with pbwt" << std::endl;
            // Todo: add lower limit to stored parameters during serialization
            if (hist_alts[1] < 10) { // dont update if < 10 alts
                for (int i = 0; i < len_data; ++i) {
                    tgt_container->model_2mc->pbwt.prev[i] = BCF_UNPACK_GENOTYPE(data[tgt_container->model_2mc->pbwt.ppa[i]]);
                }
            } else {
                tgt_container->model_2mc->pbwt.UpdateBcf(data, 1);
            }

            ret = (tgt_container->Encode2mc(tgt_container->model_2mc->pbwt.prev, len_data));
        } else {
            // std::cerr << "encoding nopbwt" << std::endl;
            ret = (tgt_container->Encode2mc(data, len_data, TWK_BCF_GT_UNPACK, 1));
        }

        if (ret > 0) ++n_variants;
        return ret;
    } else { // Otherwise.
        // tgt_container->marchetype->EncodeSymbol(1);
        tgt_container->p[tgt_container->p_len++] = 1;
        
        int ret = -1;
        if (use_pbwt) {
            tgt_container->model_nm->pbwt.UpdateBcfGeneral(data, 1); // otherwise
            ret = (tgt_container->EncodeNm(tgt_container->model_nm->pbwt.prev, len_data));
        } else {
            ret = (tgt_container->EncodeNm(data, len_data, TWK_BCF_GT_UNPACK_GENERAL, 1));
        }
        if (ret > 0) ++n_variants;
        return ret;
    }
}

int djinn_ewah_model::Encode(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) {
    if (data == nullptr) return -2;
    if (len_data % ploidy != 0) return -3;
    // Currently limited to 14 alt alleles + missing + EOV marker (total of 16).
    assert(alt_alleles < 14);

    // Todo: check that ploidy_dict model is initiated!
    std::shared_ptr<djn_ewah_model_container_t> tgt_container;

    const uint64_t tuple = ((uint64_t)len_data << 32) | ploidy;
    auto search = ploidy_map.find(tuple);
    if (search != ploidy_map.end()) {
        tgt_container = ploidy_models[search->second];
        // ploidy_dict->EncodeSymbol(search->second);
        p[p_len++] = search->second;
    } else {
        // std::cerr << "Not found. Inserting: " << len_data << "," << ploidy << "(" << tuple << ") as " << ploidy_models.size() << std::endl;
        ploidy_map[tuple] = ploidy_models.size();
        // ploidy_dict->EncodeSymbol(ploidy_models.size());
        p[p_len++] = ploidy_models.size();
        ploidy_models.push_back(std::make_shared<djn_ewah_model_container_t>(len_data, ploidy, (bool)use_pbwt));
        tgt_container = ploidy_models[ploidy_models.size() - 1];
        tgt_container->StartEncoding(use_pbwt, init);
    }
    assert(tgt_container.get() != nullptr);

    // Compute allele counts.
    memset(hist_alts, 0, 256*sizeof(uint32_t));
    for (int i = 0; i < len_data; ++i) {
        ++hist_alts[data[i]];
    }

    // for (int i = 0; i < 256; ++i) {
    //     if (hist_alts[i]) std::cerr << i << ":" << hist_alts[i] << ",";
    // }
    // std::cerr << " data2m=" << tgt_container->model_2mc->range_coder->OutSize() << " dataNm=" << tgt_container->model_nm->range_coder->OutSize() << " permute=" << permute << std::endl;

    // Biallelic, no missing, and no special EOV symbols.
    if (alt_alleles <= 2 && hist_alts[14] == 0 && hist_alts[15] == 0) {
        // tgt_container->marchetype->EncodeSymbol(0); // add archtype as 2mc
        tgt_container->p[tgt_container->p_len++] = 0;

        int ret = -1;
        if (use_pbwt) {
            // std::cerr << "encoding with pbwt" << std::endl;
            // Todo: add lower limit to stored parameters during serialization
            if (hist_alts[1] < 10) { // dont update if < 10 alts
                for (int i = 0; i < len_data; ++i) {
                    tgt_container->model_2mc->pbwt.prev[i] = data[tgt_container->model_2mc->pbwt.ppa[i]];
                }
            } else {
                tgt_container->model_2mc->pbwt.Update(data, 1);
            }

            ret = (tgt_container->Encode2mc(tgt_container->model_2mc->pbwt.prev, len_data));
        } else {
            // std::cerr << "encoding nopbwt" << std::endl;
            ret = (tgt_container->Encode2mc(data, len_data, DJN_MAP_NONE, 0));
        }

        if (ret > 0) ++n_variants;
        return ret;
    } else { // Otherwise.
        // tgt_container->marchetype->EncodeSymbol(1);
        tgt_container->p[tgt_container->p_len++] = 1;
        
        int ret = -1;
        if (use_pbwt) {
            tgt_container->model_nm->pbwt.Update(data, 1);
            ret = (tgt_container->EncodeNm(tgt_container->model_nm->pbwt.prev, len_data));
        } else {
            ret = (tgt_container->EncodeNm(data, len_data, DJN_MAP_NONE, 0));
        }
        if (ret > 0) ++n_variants;
        return ret;
    }
}

int djinn_ewah_model::DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len) {
    if (ewah_data == nullptr) return -1;
    if (ret_buffer == nullptr) return -1;

    // Decode stream archetype.
    uint8_t type = p[p_len];
    ++p_len;
    assert(type < ploidy_models.size());

    std::shared_ptr<djn_ewah_model_container_t> tgt_container = ploidy_models[type];

    return(tgt_container->DecodeNext(ewah_data,ret_ewah,ret_buffer,ret_len));
}

int djinn_ewah_model::DecodeNext(djinn_variant_t*& variant) {
     // Decode stream archetype.
    uint8_t type = p[p_len];
    ++p_len;
    std::shared_ptr<djn_ewah_model_container_t> tgt_container = ploidy_models[type];

    if (variant == nullptr) {
        variant = new djinn_variant_t;
        variant->data_alloc = tgt_container->n_samples + 65536;
        variant->data = new uint8_t[variant->data_alloc];
        variant->data_free = true;
    } else if (tgt_container->n_samples >= variant->data_alloc) {
        if (variant->data_free) delete[] variant->data;
        variant->data_alloc = tgt_container->n_samples + 65536;
        variant->data = new uint8_t[variant->data_alloc];
        variant->data_free = true;
    }

    if (q == nullptr) {
        q_alloc = tgt_container->n_samples + 65536;
        q = new uint8_t[q_alloc];
        q_len = 0;
        q_free = true;
    } else if (tgt_container->n_samples >= q_alloc) {
        if (q_free) delete[] q;
        q_alloc = tgt_container->n_samples + 65536;
        q = new uint8_t[q_alloc];
        q_len = 0;
        q_free = true;
    }
    q_len = 0; // Reset q_len for next iteration.

    variant->ploidy   = tgt_container->ploidy;
    variant->data_len = 0;
    variant->errcode  = 0;
    variant->unpacked = DJN_UN_IND;

    int ret = tgt_container->DecodeNext(q,q_len,variant->data,variant->data_len);
    
    // Compute number of alleles.
    uint32_t max_allele = 0;
    for (int i = 0; i < 256; ++i) {
        max_allele = tgt_container->hist_alts[i] != 0 ? i : max_allele;
        // if (tgt_container->hist_alts[i]) std::cerr << i << ":" << tgt_container->hist_alts[i] << ",";
    }
    // std::cerr << " max=" << max_allele << std::endl;
    variant->n_allele = max_allele + 1;

    return ret;
}

int djinn_ewah_model::DecodeNextRaw(djinn_variant_t*& variant) {
    // Decode stream archetype.
    uint8_t type = p[p_len];
    ++p_len;
    std::shared_ptr<djn_ewah_model_container_t> tgt_container = ploidy_models[type];
    return tgt_container->DecodeNextRaw(variant);
}

int djinn_ewah_model::DecodeNextRaw(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;
    // uint8_t type = ploidy_dict->DecodeSymbol();
    uint8_t type = p[p_len];
    ++p_len;

    std::shared_ptr<djn_ewah_model_container_t> tgt_container = ploidy_models[type];
    return(tgt_container->DecodeNextRaw(data, len));
}

void djinn_ewah_model::StartEncoding(bool use_pbwt, bool reset) {
    // Store parameters for decoding.
    this->use_pbwt = use_pbwt;
    this->init = reset;
    n_variants = 0;
    p_len = 0;

    // std::cerr << "[djinn_ewah_model::StartEncoding] models start encoding" << std::endl;
    for (int i = 0; i < ploidy_models.size(); ++i) {
        ploidy_models[i]->StartEncoding(use_pbwt, reset);
    }
}

size_t djinn_ewah_model::FinishEncoding() {
    if (q == nullptr) {
        q = new uint8_t[10000000];
        q_free = true;
        q_alloc = 10000000;
        q_len = 0;
    }
    q_len = 0;

    GeneralCompressor comp;
    switch(codec) {
        case (CompressionStrategy::ZSTD): comp = &ZstdCompress; break;
        case (CompressionStrategy::LZ4):  comp = &Lz4Compress;  break;
    }

    int ret = (*comp)(p, p_len, q, q_alloc, 1);
    memcpy(p, q, ret); // copy data back to p
    p_len = ret;
    // std::cerr << "[djinn_ewah_model::FinishEncoding] p_len=" << p_len << "->" << ret << std::endl;

    size_t s_models = 0;
    for (int i = 0; i < ploidy_models.size(); ++i) {
        int ret = ploidy_models[i]->FinishEncoding(q, q_alloc, codec); // use q as support buffer
        s_models += ret;
    }
    
    return ret + s_models;
}

int djinn_ewah_model::StartDecoding() {
    assert(p != nullptr);
    assert(ploidy_models.size() != 0);

    if (q == nullptr) {
        q = new uint8_t[10000000];
        q_free = true;
        q_alloc = 10000000;
        q_len = 0;
    }
    q_len = 0;

    GeneralDecompressor dec;
    switch(codec) {
        case (CompressionStrategy::ZSTD): dec = &ZstdDecompress; break;
        case (CompressionStrategy::LZ4):  dec = &Lz4Decompress;  break;
    }
    // int ZstdDecompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity)
    // std::cerr << "before decompress: " << (int)codec << std::endl;
    int ret = (*dec)(p, p_len, q, q_alloc);
    // std::cerr << "after decompress=" << p_len << "->" << ret << std::endl;
    memcpy(p, q, ret); // copy data back to p
    p_len = 0;
    // std::cerr << "p_len now=" << p_len << std::endl;
    // std::cerr << "q_alloc=" << q_alloc << std::endl;

    for (int i = 0; i <ploidy_models.size(); ++i) {
        ploidy_models[i]->StartDecoding(q,q_alloc,codec,use_pbwt, init);
    }

    return 1;
}

int djinn_ewah_model::Serialize(uint8_t* dst) const {
    // Serialize as (int,uint32_t,uint32_t,uint8_t*,ctx1,ctx2):
    // #models,p_len,p,[models...]
    uint32_t offset = 0;

    // Reserve space for offset
    offset += sizeof(uint32_t);

    // Add
    *((int*)&dst[offset]) = (int)codec; // codec used
    offset += sizeof(int);
    *((int*)&dst[offset]) = (int)ploidy_models.size(); // number of models
    offset += sizeof(int);
    *((uint32_t*)&dst[offset]) = n_variants; // number of variants
    offset += sizeof(uint32_t);

    // Serialize bit-packed controller.
    uint8_t pack = (use_pbwt << 7) | (init << 6) | (unused << 0);
    dst[offset] = pack;
    offset += sizeof(uint8_t);

    *((uint32_t*)&dst[offset]) = p_len; // data length
    offset += sizeof(uint32_t);
    
    // Store model selection data.
    memcpy(&dst[offset], p, p_len); // data
    offset += p_len;

    // Serialize each model.
    for (int i = 0; i < ploidy_models.size(); ++i) {
        offset += ploidy_models[i]->Serialize(&dst[offset]);    
    }
    *((uint32_t*)&dst[0]) = offset;
    
    return offset;
}

int djinn_ewah_model::Serialize(std::ostream& stream) const {
    // Serialize as (int,uint32_t,uint32_t,uint8_t*,ctx1,ctx2):
    // #models,p_len,p,[models...]
    uint32_t out_len = GetSerializedSize();
    stream.write((char*)&out_len, sizeof(uint32_t));

    int out_codec = (int)codec;
    stream.write((char*)&out_codec, sizeof(int));

    int n_models = ploidy_models.size();
    stream.write((char*)&n_models, sizeof(int));
    stream.write((char*)&n_variants, sizeof(uint32_t));

    // Serialize bit-packed controller.
    uint8_t pack = (use_pbwt << 7) | (init << 6) | (unused << 0);
    stream.write((char*)&pack, sizeof(uint8_t));
    stream.write((char*)&p_len, sizeof(uint32_t));
    stream.write((char*)p, p_len);

    // Serialize each model.
    for (int i = 0; i < ploidy_models.size(); ++i)
        ploidy_models[i]->Serialize(stream);

    return out_len;
}

int djinn_ewah_model::GetSerializedSize() const {
    int ret = sizeof(uint32_t) + 2*sizeof(int) + sizeof(uint32_t) + sizeof(uint8_t) + sizeof(uint32_t) + p_len;
    for (int i = 0; i < ploidy_models.size(); ++i) {
        ret += ploidy_models[i]->GetSerializedSize();
    }
    return ret;
}

int djinn_ewah_model::GetCurrentSize() const {
    int ret = p_len;
    for (int i = 0; i < ploidy_models.size(); ++i) {
        ret += ploidy_models[i]->GetCurrentSize();
    }
    return ret;
}

// Deserialize data from an external buffer.
int djinn_ewah_model::Deserialize(uint8_t* src) {
    // Reset.
    // ploidy_remap.clear();

    uint32_t offset = 0;

    // Read total offset.
    uint32_t tot_offset = *((uint32_t*)&src[offset]);
    offset += sizeof(uint32_t);

    // Codec
    codec = CompressionStrategy(*((int*)&src[offset]));
    offset += sizeof(int);

    // Read #models,#variants,packed bools
    int n_models = *((int*)&src[offset]);
    offset += sizeof(int);
    n_variants = *((uint32_t*)&src[offset]);
    offset += sizeof(uint32_t);
    uint8_t pack = src[offset];
    use_pbwt = (pack >> 7) & 1;
    init = (pack >> 6) & 1;
    unused = 0;
    offset += sizeof(uint8_t);

    // Read p_len,p
    p_len = *((uint32_t*)&src[offset]);
    offset += sizeof(uint32_t);

    if (p_cap == 0 || p == nullptr || p_len > p_cap) {
        // std::cerr << "[Deserialize] Limit. p_cap=" << p_cap << "," << "p is nullptr=" << (p == nullptr ? "yes" : "no") << ",p_len=" << p_len << "/" << p_cap << std::endl;
        if (p_free) delete[] p;
        p = new uint8_t[p_len + 65536];
        p_cap = p_len + 65536;
        p_free = true;
    }    
    
    // Store model selection data.
    memcpy(p, &src[offset], p_len);
    offset += p_len;
    
    // Read each model.
    for (int i = 0; i < n_models; ++i) {
        // First peek at their content to invoke the correct constructor as it
        // requires #samples, #ploidy, use_pbwt
        int pl = *((int*)&src[offset]);
        uint32_t n_s = *((uint32_t*)&src[offset+sizeof(int)]);

        // Update ploidy map and insert data in the correct position relative
        // the encoding order.
        const uint64_t tuple = ((uint64_t)n_s << 32) | pl;
        auto search = ploidy_map.find(tuple);
        if (search != ploidy_map.end()) {
            offset += ploidy_models[search->second]->Deserialize(&src[offset]);
        } else {
            ploidy_map[tuple] = ploidy_models.size();
            ploidy_models.push_back(std::make_shared<djinn::djn_ewah_model_container_t>(n_s, pl, (bool)use_pbwt));
            offset += ploidy_models.back()->Deserialize(&src[offset]);
        }
    }
    // std::cerr << "[Deserialize] Decoded=" << offset << "/" << tot_offset << std::endl;
    assert(offset == tot_offset);
    return offset;
}

// Deserialize data from a IO stream. This approach is more efficient
// as data does NOT have to be copied and the current model can be
// reused.
int djinn_ewah_model::Deserialize(std::istream& stream) {
    // uint32_t out_len = GetSerializedSize();
    uint32_t out_len = 0;
    stream.read((char*)&out_len, sizeof(uint32_t));

    // Codec
    int in_codec = 0;
    stream.read((char*)&in_codec, sizeof(int));
    codec = CompressionStrategy(in_codec);

    // std::cerr << "codec:" << (int)codec << std::endl;

    // n_models,n_variants
    int n_models = 0;
    stream.read((char*)&n_models, sizeof(int));
    stream.read((char*)&n_variants, sizeof(uint32_t));

    // Deserialize bit-packed controller.
    uint8_t pack = 0;
    stream.read((char*)&pack, sizeof(uint8_t));
    use_pbwt = (pack >> 7) & 1;
    init = (pack >> 6) & 1;
    unused = 0;

    stream.read((char*)&p_len, sizeof(uint32_t));
    if (p_cap == 0 || p == nullptr || p_len > p_cap) {
        // std::cerr << "[Deserialize] Limit. p_cap=" << p_cap << "," << "p is nullptr=" << (p == nullptr ? "yes" : "no") << ",p_len=" << p_len << "/" << p_cap << std::endl;
        if (p_free) delete[] p;
        p = new uint8_t[p_len + 65536];
        p_cap = p_len + 65536;
        p_free = true;
    }

    // std::cerr << "[Deserialize] Limit. p_cap=" << p_cap << "," << "p is nullptr=" << (p == nullptr ? "yes" : "no") << ",p_len=" << p_len << "/" << p_cap << std::endl;
    

    // std::cerr << "p_cap=" << p_cap << " p_len=" << p_len << std::endl;
    stream.read((char*)p, p_len);

    // Serialize each model.
    for (int i = 0; i < n_models; ++i) {
        // First peek at their content to invoke the correct constructor as it
        // requires #samples, #ploidy, use_pbwt
        int pl = 0; 
        uint32_t n_s = 0; 
        stream.read((char*)&pl,  sizeof(int));
        stream.read((char*)&n_s, sizeof(uint32_t));

        // std::cerr << "[Deserialize] #pl=" << pl << ", #n_s=" << n_s << std::endl;

        // Update ploidy map and insert data in the correct position relative
        // the encoding order.
        const uint64_t tuple = ((uint64_t)n_s << 32) | pl;
        auto search = ploidy_map.find(tuple);
        if (search != ploidy_map.end()) {
            ploidy_models[search->second]->Deserialize(stream);
        } else {
            // std::cerr << "[Deserialize] Adding map [" << tuple << "] for [" << n_s << "," << pl << "]" << std::endl; 
            ploidy_map[tuple] = ploidy_models.size();
            ploidy_models.push_back(std::make_shared<djinn::djn_ewah_model_container_t>(n_s, pl, (bool)use_pbwt));
            // std::cerr << "before deserialize: " << ploidy_models.size() << std::endl;
            ploidy_models.back()->Deserialize(stream);
        }
    }

    return stream.tellg();
}

/*======   Container   ======*/

constexpr uint32_t djn_ewah_model_container_t::nm_ref_bits[16];

djn_ewah_model_container_t::djn_ewah_model_container_t(int64_t n_s, int pl, bool use_pbwt) : 
    use_pbwt(use_pbwt),
    ploidy(pl), n_samples(n_s), n_variants(0),
    n_samples_wah(std::ceil((float)n_samples / 32) * 32), 
    n_samples_wah_nm(std::ceil((float)n_samples * 4/32) * 8),
    n_wah(n_samples_wah_nm / 8),
    wah_bitmaps(new uint32_t[n_wah]), 
    p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
    model_2mc(std::make_shared<djn_ewah_model_t>()),
    model_nm(std::make_shared<djn_ewah_model_t>())
{
    assert(n_s % pl == 0); // #samples/#ploidy must be divisible
}

djn_ewah_model_container_t::djn_ewah_model_container_t(int64_t n_s, int pl, bool use_pbwt, uint8_t* src, uint32_t src_len) : 
    use_pbwt(use_pbwt),
    ploidy(pl), n_samples(n_s), n_variants(0),
    n_samples_wah(std::ceil((float)n_samples / 32) * 32), 
    n_samples_wah_nm(std::ceil((float)n_samples * 4/32) * 8),
    n_wah(n_samples_wah_nm / 8),
    wah_bitmaps(new uint32_t[n_wah]), 
    p(src), p_len(src_len), p_cap(0), p_free(false),
    model_2mc(std::make_shared<djn_ewah_model_t>()),
    model_nm(std::make_shared<djn_ewah_model_t>())
{
    assert(n_s % pl == 0); // #samples/#ploidy must be divisible
}

djn_ewah_model_container_t::~djn_ewah_model_container_t() {
    if (p_free) delete[] p;
    delete[] wah_bitmaps;
}

void djn_ewah_model_container_t::StartEncoding(bool use_pbwt, bool reset) {
    if (model_2mc.get() == nullptr) return;
    if (model_nm.get() == nullptr)  return;

    // std::cerr << "in start encoding: " << model_2mc->pbwt.n_symbols << std::endl; 
    if (use_pbwt) {
        if (model_2mc->pbwt.n_symbols == 0) {
            // std::cerr << "init pbwt: " << n_samples << "," << 2 << std::endl;
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_2mc->pbwt.Initiate(n_samples, 2);
        }

        if (model_nm->pbwt.n_symbols == 0) {
            // std::cerr << "init pbwt: " << n_samples << "," << 16 << std::endl;
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_nm->pbwt.Initiate(n_samples, 16);
        }
    }

    if (reset) {
        // std::cerr << "[djn_ewah_model_container_t::StartEncoding] resetting" << std::endl;
        model_2mc->reset();
        model_nm->reset();
    } else {
        model_2mc->n_variants = 0;
        model_nm->n_variants = 0;
    }

    this->use_pbwt = use_pbwt;
    n_variants = 0;
    p_len = 0;

    model_2mc->StartEncoding(use_pbwt, reset);
    model_nm->StartEncoding(use_pbwt, reset);
}

size_t djn_ewah_model_container_t::FinishEncoding(uint8_t* support_buffer, uint32_t support_cap, CompressionStrategy strat) {
    if (model_2mc.get() == nullptr) return -1;
    if (model_nm.get() == nullptr) return -1;

    int s_2mc = model_2mc->FinishEncoding(support_buffer, support_cap, strat);
    int s_nm  = model_nm->FinishEncoding(support_buffer, support_cap, strat);

    GeneralCompressor comp;
    switch(strat) {
        case (CompressionStrategy::ZSTD): comp = &ZstdCompress; break;
        case (CompressionStrategy::LZ4):  comp = &Lz4Compress; break;
    }

    int ret = (*comp)(p, p_len, support_buffer, support_cap, 1);
    memcpy(p, support_buffer, ret); // copy data back to p
    p_len = ret;

    size_t s_rc  = ret;
    // std::cerr << "container finish=" << s_rc << " and " << s_2mc << " and " << s_nm << std::endl; 
    
    return s_rc + s_2mc + s_nm;
}

void djn_ewah_model_container_t::StartDecoding(uint8_t* support_buffer, uint32_t support_cap, CompressionStrategy strat, bool use_pbwt, bool reset) {
    if (model_2mc.get() == nullptr) return;
    if (model_nm.get() == nullptr) return;

    // std::cerr << "in start decoding: " << model_2mc->pbwt.n_symbols << std::endl; 
    if (use_pbwt) {
        if (model_2mc->pbwt.n_symbols == 0) {
            // std::cerr << "init pbwt: " << n_samples << "," << 2 << std::endl;
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_2mc->pbwt.Initiate(n_samples, 2);
        }

        if (model_nm->pbwt.n_symbols == 0) {
            // std::cerr << "init pbwt: " << n_samples << "," << 16 << std::endl;
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_nm->pbwt.Initiate(n_samples, 16);
        }
    }

    if (reset) {
        // std::cerr << "[djn_ewah_model_container_t::StartDecoding] resetting" << std::endl;
        model_2mc->reset();
        model_nm->reset();
    } else {
        model_2mc->n_variants = 0;
        model_nm->n_variants = 0;
    }
    // std::cerr << "[djn_ewah_model_container_t::StartDecoding] p_len=" << p_len << std::endl;

    this->use_pbwt = use_pbwt;

    GeneralDecompressor dec;
    switch(strat) {
        case (CompressionStrategy::ZSTD): dec = &ZstdDecompress; break;
        case (CompressionStrategy::LZ4):  dec = &Lz4Decompress;  break;
    }
    // int ZstdDecompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity)
    // std::cerr << "[cont]" << p_len << "," << support_cap << std::endl;
    int ret = (*dec)(p, p_len, support_buffer, support_cap);
    // std::cerr << "inflated=" << p_len << "->" << ret << std::endl;
    memcpy(p, support_buffer, ret); // copy data back to p
    p_len = 0;

    model_2mc->StartDecoding(support_buffer, support_cap, strat, use_pbwt, reset);
    model_nm->StartDecoding(support_buffer, support_cap, strat, use_pbwt, reset);
}

int djn_ewah_model_container_t::Encode2mc(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;
    
    if (wah_bitmaps == nullptr) {
        n_samples_wah = std::ceil((float)n_samples / 32) * 32;
        n_samples_wah_nm = std::ceil((float)n_samples * 4/32) * 8; // Can fit eight 4-bit entries in a 32-bit word
        n_wah = n_samples_wah_nm / 8; 
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    // uint32_t alts = 0;
    for (int i = 0; i < n_samples; ++i) {
        if (data[i]) {
            wah_bitmaps[i / 32] |= 1L << (i % 32);
            // ++alts;
        }
    }

    // if (alts > 30) {
    //     std::cerr << "alts: " << alts;
    //     for (int i = 0; i < n_wah/4; ++i) {
    //         std::cerr << " " << std::bitset<32>(wah_bitmaps[i]);
    //     }
    //     std::cerr << std::endl;
    // }
    // std::cerr << "encoding 2mc=" << n_wah/4 << "->" << n_wah/4*32 << " alts=" << alts << std::endl;

    // return EncodeWah(wah_bitmaps, n_wah/4); // 4 bytes in a 32-bit word
    return EncodeWah(wah_bitmaps, n_samples_wah >> 5); // n_samples_wah / 32
}

int djn_ewah_model_container_t::Encode2mc(uint8_t* data, uint32_t len, const uint8_t* map, const int shift) {
    if (data == nullptr) return -1;
    if (n_samples == 0)  return -2;
    if (map == nullptr)  return -3;

    if (wah_bitmaps == nullptr) {
        n_samples_wah = std::ceil((float)n_samples / 32) * 32;
        n_samples_wah_nm = std::ceil((float)n_samples * 4/32) * 8; // Can fit eight 4-bit entries in a 32-bit word
        n_wah = n_samples_wah_nm / 8; 
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    // uint32_t alts = 0;
    for (int i = 0; i < n_samples; ++i) {
        if (map[data[i] >> shift]) {
            wah_bitmaps[i / 32] |= 1L << (i % 32);
            // ++alts;
        }
    }
    // std::cerr << "encoding 2mc=" << n_wah/4 << "->" << n_wah/4*32 << " alts=" << alts << std::endl;

    // return EncodeWah(wah_bitmaps, n_wah/4); // 4 bytes in a 32-bit word
    return EncodeWah(wah_bitmaps, n_samples_wah >> 5); // n_samples_wah / 32
}

int djn_ewah_model_container_t::EncodeNm(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;
    
    if (wah_bitmaps == nullptr) {
        n_wah = std::ceil((float)n_samples / 32) * 4;
        n_samples_wah = (n_wah * 32) / 4;
        n_samples_wah_nm = std::ceil((float)n_samples * 4 / 32) * 8;
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    // Todo: move into ResetBitmaps() functions
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    for (int i = 0; i < n_samples; ++i) {
        wah_bitmaps[i / 8] |= (uint64_t)data[i] << (4*(i % 8));
    }

    // for (int i = 0; i < n_wah; ++i) std::cerr << std::bitset<64>(wah_bitmaps[i]) << " ";
    // std::cerr << std::endl;

    // return EncodeWahNm(wah_bitmaps, n_samples_wah/4);
    return EncodeWahNm(wah_bitmaps, n_samples_wah_nm >> 3); // n_samples_wah_nm / 8
    // return -1;
}

int djn_ewah_model_container_t::EncodeNm(uint8_t* data, uint32_t len, const uint8_t* map, const int shift) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;
    if (map == nullptr) return -3;

    if (wah_bitmaps == nullptr) {
        n_samples_wah = std::ceil((float)n_samples / 32) * 32;
        n_samples_wah_nm = std::ceil((float)n_samples * 4/32) * 8; // Can fit eight 4-bit entries in a 32-bit word
        n_wah = n_samples_wah_nm / 8; 
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    // Todo: move into ResetBitmaps() functions
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    for (int i = 0; i < n_samples; ++i) {
        wah_bitmaps[i / 8] |= (uint64_t)map[data[i] >> shift] << (4*(i % 8));
    }

    // for (int i = 0; i < n_wah; ++i) std::cerr << std::bitset<64>(wah_bitmaps[i]) << " ";
    // std::cerr << std::endl;

    // return EncodeWahNm(wah_bitmaps, n_samples_wah/4);
    return EncodeWahNm(wah_bitmaps, n_samples_wah_nm >> 3); // n_samples_wah_nm / 8
    // return -1;
}

int djn_ewah_model_container_t::DecodeRaw_nm(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;
    if (model_nm.get() == nullptr) return -2;

    int64_t n_samples_obs = 0;
    int objects = 0;

    // Compute als
    memset(hist_alts, 0, 256*sizeof(uint32_t));

    while(true) {
        // Emit empty EWAH marker.
        djinn_ewah_t* ewah = (djinn_ewah_t*)&data[len]; 
        ewah->reset();
        len += sizeof(djinn_ewah_t);
        
        djinn_ewah_t* e = (djinn_ewah_t*)&model_nm->p[model_nm->p_len]; 
        // std::cerr << "EWAH: ref=" << e->ref << ",clean=" << e->clean << ",dirty=" << e->dirty << " offset=" << model_nm->p_len << "/" << model_nm->p_cap << std::endl;
        assert(e->clean + e->dirty > 0);
        ewah->ref   = e->ref;
        ewah->clean = e->clean;
        ewah->dirty = e->dirty;
        model_nm->p_len += sizeof(djinn_ewah_t);
        n_samples_obs += ewah->clean*8;
        n_samples_obs += ewah->dirty*8;
        hist_alts[ewah->ref & 15] += 8*ewah->clean;

        for (int i = 0; i < ewah->dirty; ++i) {
            uint32_t r = *((const uint32_t*)&model_nm->p[model_nm->p_len]); // copy
            *((uint32_t*)&data[len]) = r;
            for (int j = 0; j < 8; ++j) {
                ++hist_alts[r & 15];
                r >>= 4;
            }
            len += sizeof(uint32_t);
            model_nm->p_len += sizeof(uint32_t);
        }
        ++objects;

        // std::cerr << "samples now=" << n_samples_obs << "/" << n_samples_wah_nm << std::endl;
        if (n_samples_obs == n_samples_wah_nm) {
            // std::cerr << "obs=" << n_samples_obs << "/" << n_samples_wah << std::endl;
            break;
        }

        if (n_samples_obs > n_samples_wah_nm) {
            std::cerr << "[djn_ewah_model_container_t::DecodeRaw_nm] Decompression corruption: " << n_samples_obs << "/" << n_samples_wah  << std::endl;
            exit(1);
        }
    }

    // std::cerr << "[decoded nm] " << hist_alts[0] << "," <<hist_alts[1] << std::endl;
    return objects;
}

int djn_ewah_model_container_t::EncodeWah(uint32_t* wah, uint32_t len) { // input WAH-encoded data
    if (wah == nullptr) return -1;
    if (model_2mc.get() == nullptr) return -1;

    // Resize if necessary.
    if (model_2mc->p_len + n_samples + 65536 > model_2mc->p_cap) {
        const uint32_t rc_size = model_2mc->p_len;
        // std::cerr << "[djn_ewah_model_container_t::EncodeWah][RESIZE] resizing from: " << rc_size << "->" << 2*rc_size << std::endl;
        uint8_t* prev = model_2mc->p; // old
        model_2mc->p_cap = model_2mc->p_len + 2*n_samples + 65536;
        model_2mc->p = new uint8_t[model_2mc->p_cap]; // double size. should rarely occur
        memcpy(model_2mc->p, prev, rc_size);
        if (model_2mc->p_free) delete[] prev;
        model_2mc->p_free = true;
    }

    // Debug
    uint32_t n_objs = 1;
    uint32_t n_obs  = 0;
    
    // djinn_ewah_t ewah; ewah.reset();
    djinn_ewah_t* ewah = (djinn_ewah_t*)&model_2mc->p[model_2mc->p_len];
    ewah->reset();
    model_2mc->p_len += sizeof(djinn_ewah_t);
    // uint32_t wah_ref = (wah[0] & 15);

    for (int i = 0; i < len; ++i) {
        // Is dirty
        if (wah[i] != 0 && wah[i] != std::numeric_limits<uint32_t>::max()) {
            // if (model_2mc->n_variants == 524) std::cerr << "dirty: " << std::bitset<32>(wah[i]) << std::endl;
            ++ewah->dirty;
            *((uint32_t*)&model_2mc->p[model_2mc->p_len]) = wah[i];
            model_2mc->p_len += sizeof(uint32_t);
            ++n_obs;
        } 
        // Is clean
        else {
            // Dirty have been set then make new EWAH
            if (ewah->dirty) {
                // Only make a new EWAH if anything is set
                n_obs += ewah->clean;
                assert(ewah->clean + ewah->dirty > 0);
                ewah = (djinn_ewah_t*)&model_2mc->p[model_2mc->p_len];
                ewah->reset();
                model_2mc->p_len += sizeof(djinn_ewah_t);
                ++ewah->clean;
                ewah->ref = (wah[i] & 1);
                ++n_objs;
            } 
            // No dirty have been set
            else {
                if (ewah->ref == (wah[i] & 1)) ++ewah->clean;
                else { // Different clean word
                    n_obs += ewah->clean;
                    if (ewah->clean) {
                        assert(ewah->clean + ewah->dirty > 0);
                        // if (model_2mc->n_variants == 524) std::cerr << "ewah: ref=" << ewah->ref << ",clean=" << ewah->clean << ",dirty=" << ewah->dirty << std::endl;
                        ewah = (djinn_ewah_t*)&model_2mc->p[model_2mc->p_len];
                        ewah->reset();
                        model_2mc->p_len += sizeof(djinn_ewah_t);
                    }
                    ++ewah->clean;
                    ewah->ref = (wah[i] & 1);
                    ++n_objs;
                }
            }
        }
    }
    // if (model_2mc->n_variants == 524 && ewah->clean) std::cerr << "ewah: ref=" << ewah->ref << ",clean=" << ewah->clean << ",dirty=" << ewah->dirty << std::endl;
    n_obs += ewah->clean;
    // std::cerr << "n_obs=" << n_obs << "/" << len << std::endl;
    assert(n_obs == len);

    ++model_2mc->n_variants;
    ++n_variants;
    return n_objs;
}

int djn_ewah_model_container_t::EncodeWahNm(uint32_t* wah, uint32_t len) { // input WAH-encoded data
    if (wah == nullptr) return -1;
    if (model_nm.get() == nullptr) return -1;

    // Resize if necessary.
    if (model_nm->p_len + n_samples + 65536 > model_nm->p_cap) {
        const uint32_t rc_size = model_nm->p_len;
        // std::cerr << "[djn_ewah_model_container_t::EncodeWah][RESIZE] resizing from: " << rc_size << "->" << 2*rc_size << std::endl;
        uint8_t* prev = model_nm->p; // old
        model_nm->p_cap = model_nm->p_len + 2*n_samples + 65536;
        model_nm->p = new uint8_t[model_nm->p_cap]; // double size. should rarely occur
        memcpy(model_nm->p, prev, rc_size);
        if (model_nm->p_free) delete[] prev;
        model_nm->p_free = true;
    }

    // Debug
    uint32_t n_objs = 1;
    uint32_t n_obs  = 0;
    
    // djinn_ewah_t ewah; ewah.reset();
    djinn_ewah_t* ewah = (djinn_ewah_t*)&model_nm->p[model_nm->p_len];
    ewah->reset();
    model_nm->p_len += sizeof(djinn_ewah_t);
    // uint32_t wah_ref = (wah[0] & 15);

    for (int i = 0; i < len; ++i) {
        // Is dirty
        if (wah[i] != djn_ewah_model_container_t::nm_ref_bits[wah[i] & 15]) {
            // if (model_nm->n_variants == 0) std::cerr << "dirty nm: " << std::bitset<32>(wah[i]) << std::endl;
            ++ewah->dirty;
            *((uint32_t*)&model_nm->p[model_nm->p_len]) = wah[i];
            model_nm->p_len += sizeof(uint32_t);
            ++n_obs;
        } 
        // Is clean
        else {
            // Dirty have been set then make new EWAH
            if (ewah->dirty) {
                // Only make a new EWAH if anything is set
                n_obs += ewah->clean;
                // if (model_nm->n_variants == 0) std::cerr << "ewah: ref=" << ewah->ref << ",clean=" << ewah->clean << ",dirty=" << ewah->dirty << std::endl;
                assert(ewah->clean + ewah->dirty > 0);
                ewah = (djinn_ewah_t*)&model_nm->p[model_nm->p_len];
                ewah->reset();
                model_nm->p_len += sizeof(djinn_ewah_t);
                ++ewah->clean;
                ewah->ref = (wah[i] & 15);
                ++n_objs;
            } 
            // No dirty have been set
            else {
                if (ewah->ref == (wah[i] & 15)) ++ewah->clean;
                else { // Different clean word
                    n_obs += ewah->clean;
                    if (ewah->clean) {
                        assert(ewah->clean + ewah->dirty > 0);
                        // if (model_nm->n_variants == 0) std::cerr << "ewah: ref=" << ewah->ref << ",clean=" << ewah->clean << ",dirty=" << ewah->dirty << std::endl;
                        ewah = (djinn_ewah_t*)&model_nm->p[model_nm->p_len];
                        ewah->reset();
                        model_nm->p_len += sizeof(djinn_ewah_t);
                    }
                    ++ewah->clean;
                    ewah->ref = (wah[i] & 15);
                    ++n_objs;
                }
            }
        }
    }
    // if (model_nm->n_variants == 0 && ewah->clean) std::cerr << "ewah: ref=" << ewah->ref << ",clean=" << ewah->clean << ",dirty=" << ewah->dirty << std::endl;
    n_obs += ewah->clean;
    // std::cerr << "n_obs=" << n_obs << "/" << len << std::endl;
    assert(n_obs == len);

    ++model_nm->n_variants;
    ++n_variants;
    return n_objs;
}

int djn_ewah_model_container_t::DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len) {
    if (ewah_data == nullptr)  return -1;
    if (ret_buffer == nullptr) return -1;
    if (model_2mc.get() == nullptr) return -1;
    if (model_nm.get() == nullptr)  return -1;

    // Decode stream archetype.
    uint8_t type = p[p_len];
    ++p_len;

    size_t ret_ewah_init = ret_ewah;
    int objs = 0;
    switch(type) {
    case 0: objs = DecodeRaw(ewah_data, ret_ewah);    break;
    case 1: objs = DecodeRaw_nm(ewah_data, ret_ewah); break;
    default: std::cerr << "[djn_ewah_model_container_t::DecodeNext] decoding error: " << (int)type << " (valid=[0,1])" << std::endl; return -1;
    }

    if (objs <= 0) return -1;

    if (use_pbwt) {
        if (type == 0) {
            if (hist_alts[1] >= 10) {
                // std::cerr << "[2nm usepbwt]" << std::endl;
                model_2mc->pbwt.ReverseUpdateEWAH(ewah_data, ret_ewah, ret_buffer); 
                ret_len = n_samples;
            } else {
                // std::cerr << "[2nm NOpbwt]" << std::endl;
                // Unpack EWAH into literals according to current PPA
                uint32_t local_offset = ret_ewah_init;
                uint32_t ret_pos = 0;
                for (int j = 0; j < objs; ++j) {
                    djinn_ewah_t* ewah = (djinn_ewah_t*)&ewah_data[local_offset];
                    local_offset += sizeof(djinn_ewah_t);

                    // Clean words.
                    uint32_t to = ret_pos + ewah->clean*32 > n_samples ? n_samples : ret_pos + ewah->clean*32;
                    for (int i = ret_pos; i < to; ++i) {
                        ret_buffer[model_2mc->pbwt.ppa[i]] = (ewah->ref & 1);
                    }
                    ret_pos = to;

                    for (int i = 0; i < ewah->dirty; ++i) {
                        to = ret_pos + 32 > n_samples ? n_samples : ret_pos + 32;
                        
                        uint32_t dirty = *((uint32_t*)(&ewah_data[local_offset])); // copy
                        for (int j = ret_pos; j < to; ++j) {
                            ret_buffer[model_2mc->pbwt.ppa[j]] = (dirty & 1);
                            dirty >>= 1;
                        }
                        local_offset += sizeof(uint32_t);
                        ret_pos = to;
                    }
                    assert(ret_pos <= n_samples);
                    assert(local_offset <= ret_ewah);
                }
                assert(ret_pos == n_samples);
                ret_len = n_samples;
            }
        } else { // is NM
            // std::cerr << "[nm usepbwt]" << std::endl;
            model_nm->pbwt.ReverseUpdateEWAHNm(ewah_data, ret_ewah, ret_buffer);
            ret_len = n_samples;
        }
    } else { // not using PBWT
        if (type == 0) {
            // Unpack EWAH into literals
            uint32_t local_offset = ret_ewah_init;
            uint32_t ret_pos = 0;
            for (int j = 0; j < objs; ++j) {
                djinn_ewah_t* ewah = (djinn_ewah_t*)&ewah_data[local_offset];
                local_offset += sizeof(djinn_ewah_t);
                
                // Clean words.
                uint32_t to = ret_pos + ewah->clean*32 > n_samples ? n_samples : ret_pos + ewah->clean*32;
                for (int i = ret_pos; i < to; ++i) {
                    ret_buffer[i] = (ewah->ref & 1);
                }
                ret_pos = to;

                for (int i = 0; i < ewah->dirty; ++i) {
                    to = ret_pos + 32 > n_samples ? n_samples : ret_pos + 32;
                    
                    uint32_t dirty = *((uint32_t*)(&ewah_data[local_offset])); // copy
                    for (int j = ret_pos; j < to; ++j) {
                        ret_buffer[j] = (dirty & 1);
                        dirty >>= 1;
                    }
                    local_offset += sizeof(uint32_t);
                    ret_pos = to;
                }
                assert(ret_pos <= n_samples);
                assert(local_offset <= ret_ewah);
            }
            assert(ret_pos == n_samples);
            ret_len = n_samples;
        } else {
            // Unpack EWAH into literals
            uint32_t local_offset = ret_ewah_init;
            uint32_t n_s_obs = 0;

            for (int j = 0; j < objs; ++j) {
                djinn_ewah_t* ewah = (djinn_ewah_t*)&ewah_data[local_offset];
                local_offset += sizeof(djinn_ewah_t);

                // Clean words.
                uint64_t to = n_s_obs + ewah->clean * 8 > n_samples ? n_samples : n_s_obs + ewah->clean * 8;
                for (int i = n_s_obs; i < to; ++i) {
                    ret_buffer[i] = (ewah->ref & 15); // update prev when non-zero
                }
                n_s_obs = to;

                // Loop over dirty bitmaps.
                for (int i = 0; i < ewah->dirty; ++i) {
                    to = n_s_obs + 8 > n_samples ? n_samples : n_s_obs + 8;
                    assert(n_s_obs < n_samples);
                    assert(to <= n_samples);
                    
                    uint32_t dirty = *((uint32_t*)(&ewah_data[local_offset])); // copy
                    for (int j = n_s_obs; j < to; ++j) {
                        ret_buffer[j] = (dirty & 15);
                        dirty >>= 4;
                    }
                    local_offset += sizeof(uint32_t);
                    n_s_obs = to;
                }
                assert(local_offset <= ret_ewah);
            }
            assert(n_s_obs == n_samples);
            ret_len = n_samples;
        }
    }

    return objs;
}

int djn_ewah_model_container_t::DecodeNextRaw(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;
    if (model_2mc.get() == nullptr) return -1;
    if (model_nm.get() == nullptr) return -1;
    
    // Decode stream archetype.
    uint8_t type = p[p_len];
    ++p_len;

    switch(type) {
    case 0: return(DecodeRaw(data, len)); break;
    case 1: return(DecodeRaw_nm(data, len)); break;
    default: std::cerr << "[djn_ewah_model_container_t::DecodeNextRaw] decoding error: " << (int)type << " (valid=[0,1])" << std::endl; return -1;
    }

    // Never reached.
    return -3;
}

// Return raw, potentially permuted, EWAH encoding
int djn_ewah_model_container_t::DecodeRaw(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;
    if (model_2mc.get() == nullptr) return -2;

    int64_t n_samples_obs = 0;
    int objects = 0;

    // Compute als
    memset(hist_alts, 0, 256*sizeof(uint32_t));

    while(true) {
        // Emit empty EWAH marker.
        djinn_ewah_t* ewah = (djinn_ewah_t*)&data[len]; 
        ewah->reset();
        len += sizeof(djinn_ewah_t);
        
        djinn_ewah_t* e = (djinn_ewah_t*)&model_2mc->p[model_2mc->p_len]; 
        // std::cerr << "EWAH: ref=" << e->ref << ",clean=" << e->clean << ",dirty=" << e->dirty << " offset=" << model_2mc->p_len << "/" << model_2mc->p_cap << std::endl;
        assert(e->clean + e->dirty > 0);
        ewah->ref   = e->ref;
        ewah->clean = e->clean;
        ewah->dirty = e->dirty;
        model_2mc->p_len += sizeof(djinn_ewah_t);
        n_samples_obs += ewah->clean*32;
        n_samples_obs += ewah->dirty*32;
        hist_alts[ewah->ref & 1] += 32*ewah->clean;

        for (int i = 0; i < ewah->dirty; ++i) {
            const uint32_t* r = (const uint32_t*)&model_2mc->p[model_2mc->p_len];
            *((uint32_t*)&data[len]) = *r;
            hist_alts[1] += __builtin_popcount(*r);
            hist_alts[0] += __builtin_popcount(~(*r));
            len += sizeof(uint32_t);
            model_2mc->p_len += sizeof(uint32_t);
        }
        ++objects;

        if (n_samples_obs == n_samples_wah) {
            // std::cerr << "obs=" << n_samples_obs << "/" << n_samples_wah << std::endl;
            break;
        }

        if (n_samples_obs > n_samples_wah) {
            std::cerr << "[djn_ewah_model_container_t::DecodeRaw] Decompression corruption: " << n_samples_obs << "/" << n_samples_wah  << std::endl;
            exit(1);
        }
    }

    std::cerr << "[decoded 2nm] " << hist_alts[0] << "," <<hist_alts[1] << std::endl;

    return objects;
}

int djn_ewah_model_container_t::DecodeNextRaw(djinn_variant_t*& variant) {
    if (variant == nullptr) {
        variant = new djinn_variant_t;
        variant->data_alloc = n_samples + 65536;
        variant->data = new uint8_t[variant->data_alloc];
        variant->data_free = true;
    } else if (n_samples >= variant->data_alloc) {
        if (variant->data_free) delete[] variant->data;
        variant->data_alloc = n_samples + 65536;
        variant->data = new uint8_t[variant->data_alloc];
        variant->data_free = true;
    }

    variant->ploidy   = ploidy;
    variant->data_len = 0;
    variant->errcode  = 0;
    variant->unpacked = DJN_UN_EWAH;

    // Decode stream archetype.
    uint8_t type = p[p_len++];

    int ret = 0;
    switch(type) {
    case 0: ret = DecodeRaw(variant->data, variant->data_len);    break;
    case 1: ret = DecodeRaw_nm(variant->data, variant->data_len); break;
    default: std::cerr << "[djn_ctx_model_container_t::DecodeNextRaw] decoding error: " << (int)type << " (valid=[0,1])" << std::endl; return -1;
    }

    if (ret <= 0) return ret;

    if (variant->d == nullptr) {
        variant->d = new djn_variant_dec_t;
    }

    if (ret >= variant->d->m_ewah) {
        delete[] variant->d->ewah;
        variant->d->ewah = new djinn_ewah_t*[ret + 512];
        variant->d->m_ewah = ret + 512;
    }

    if (variant->d->m_dirty < n_samples_wah_nm) {
        delete[] variant->d->dirty;
        variant->d->dirty = new uint32_t*[n_samples_wah_nm];
        variant->d->m_dirty = n_samples_wah_nm;
    }

    // Reset
    variant->d->n_ewah  = 0;
    variant->d->n_dirty = 0;
    // Set bitmap type.
    variant->d->dirty_type = type;

    // Construct EWAH mapping
    uint32_t local_offset = 0;
    uint32_t ret_pos = 0;

    for (int j = 0; j < ret; ++j) {
        djinn_ewah_t* ewah = (djinn_ewah_t*)&variant->data[local_offset];
        variant->d->ewah[variant->d->n_ewah++] = (djinn_ewah_t*)&variant->data[local_offset];
        local_offset += sizeof(djinn_ewah_t);
    
        // Clean words.
        uint32_t to = ret_pos + ewah->clean*32 > n_samples ? n_samples : ret_pos + ewah->clean*32;
        ret_pos = to;

        for (int i = 0; i < ewah->dirty; ++i) {
            variant->d->dirty[variant->d->n_dirty++] = (uint32_t*)(&variant->data[local_offset]);
            to = ret_pos + 32 > n_samples ? n_samples : ret_pos + 32;
            local_offset += sizeof(uint32_t);
            ret_pos = to;
        }
        assert(ret_pos <= n_samples);
        assert(local_offset <= variant->data_len);
    }
    assert(ret_pos == n_samples);

    // Compute number of alleles.
    uint32_t max_allele = 0;
    for (int i = 0; i < 256; ++i) {
        max_allele = hist_alts[i] != 0 ? i : max_allele;
        // if (tgt_container->hist_alts[i]) std::cerr << i << ":" << tgt_container->hist_alts[i] << ",";
    }
    // std::cerr << " max=" << max_allele << std::endl;
    variant->n_allele = max_allele + 1;

    return 1;
}

int djn_ewah_model_container_t::Serialize(uint8_t* dst) const {
    // Serialize as (int,uint32_t,uint32_t,uint8_t*,ctx1,ctx2):
    // ploidy,n_samples,n_variants,p_len,p,model_2mc,model_nm
    uint32_t offset = 0;
    *((int*)&dst[offset]) = ploidy; // ploidy
    offset += sizeof(int);
    *((uint32_t*)&dst[offset]) = n_samples; // number of samples
    offset += sizeof(uint32_t);
    *((uint32_t*)&dst[offset]) = n_variants; // number of variants
    offset += sizeof(uint32_t);
    *((uint32_t*)&dst[offset]) = p_len; // data length
    offset += sizeof(uint32_t);
    memcpy(&dst[offset], p, p_len); // data
    offset += p_len;
    offset += model_2mc->Serialize(&dst[offset]);
    offset += model_nm->Serialize(&dst[offset]);
    return offset;
}

int djn_ewah_model_container_t::Serialize(std::ostream& stream) const {
    // Serialize as (int,uint32_t,uint32_t,uint8_t*,ctx1,ctx2):
    // ploidy,n_samples,n_variants,p_len,p,model_2mc,model_nm
    stream.write((char*)&ploidy, sizeof(int));
    stream.write((char*)&n_samples, sizeof(uint32_t));
    stream.write((char*)&n_variants, sizeof(uint32_t));
    stream.write((char*)&p_len, sizeof(uint32_t));
    stream.write((char*)p, p_len);
    model_2mc->Serialize(stream);
    model_nm->Serialize(stream);
    return stream.tellp();
}

int djn_ewah_model_container_t::GetSerializedSize() const {
    int ret = sizeof(int) + 3*sizeof(uint32_t) + p_len + model_2mc->GetSerializedSize() + model_nm->GetSerializedSize();
    return ret;
}

int djn_ewah_model_container_t::GetCurrentSize() const {
    // int ret = range_coder->OutSize();
    int ret = p_len;
    ret += model_2mc->p_len;
    ret += model_nm->p_len;
    return ret;
}

int djn_ewah_model_container_t::Deserialize(uint8_t* dst) {
    uint32_t offset = 0;
    int pl = *((int*)&dst[offset]);
    // std::cerr << "pl=" << pl << " ploidy=" << ploidy << std::endl;
    assert(pl == ploidy);
    offset += sizeof(int);
    uint32_t n_s = *((uint32_t*)&dst[offset]);
    assert(n_s == n_samples);
    offset += sizeof(uint32_t);
    n_variants = *((uint32_t*)&dst[offset]);
    offset += sizeof(uint32_t);
    p_len = *((uint32_t*)&dst[offset]);
    offset += sizeof(uint32_t);

    // initiate a buffer if there is none or it's too small
    if (p_cap == 0 || p == nullptr || p_len > p_cap) {
        // std::cerr << "[Deserialize] Limit. p_cap=" << p_cap << "," << "p is nullptr=" << (p == nullptr ? "yes" : "no") << ",p_len=" << p_len << "/" << p_cap << std::endl;
        if (p_free) delete[] p;
        p = new uint8_t[p_len + 65536];
        p_cap = p_len + 65536;
        p_free = true;
    }

    memcpy(p, &dst[offset], p_len); // data
    offset += p_len;
    // Todo objects
    offset += model_2mc->Deserialize(&dst[offset]);
    offset += model_nm->Deserialize(&dst[offset]);

    return(offset);
}

int djn_ewah_model_container_t::Deserialize(std::istream& stream) {
    // #pl and #n_s read outside of this function in Deserialize() for
    // the parent.
    //
    // stream.read((char*)&ploidy, sizeof(int));
    // stream.read((char*)&n_samples, sizeof(uint32_t));
    stream.read((char*)&n_variants, sizeof(uint32_t));
    stream.read((char*)&p_len, sizeof(uint32_t));
    
    // initiate a buffer if there is none or it's too small
    if (p_cap == 0 || p == nullptr || p_len > p_cap) {
        // std::cerr << "[Deserialize] Limit. p_cap=" << p_cap << "," << "p is nullptr=" << (p == nullptr ? "yes" : "no") << ",p_len=" << p_len << "/" << p_cap << std::endl;
        if (p_free) delete[] p;
        p = new uint8_t[p_len + 65536];
        p_cap = p_len + 65536;
        p_free = true;
    }

    stream.read((char*)p, p_len);
    model_2mc->Deserialize(stream);
    model_nm->Deserialize(stream);
    return stream.tellg();
}


}