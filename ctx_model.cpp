#include <cstring> //memcpy
#include "ctx_model.h"

namespace djinn {

/*======   Context container   ======*/

djn_ctx_model_t::djn_ctx_model_t() :
    p(nullptr),
    p_len(0), p_cap(0), p_free(false),
    n_variants(0)
{
    
}

djn_ctx_model_t::~djn_ctx_model_t() {
    if (p_free) delete[] p;
}

void djn_ctx_model_t::Initiate2mc() {
    range_coder = std::make_shared<RangeCoder>();
    mref = std::make_shared<GeneralModel>(2, 512, 18, 32, range_coder);
    mlog_rle = std::make_shared<GeneralModel>(16, 32, 18, 16, range_coder);  // 2^32 max -> log2() = 5
    mrle = std::make_shared<GeneralModel>(256, 64, 18, 32, range_coder); // 1 bit ref, 5 bit alt
    mrle2_1 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    mrle2_2 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    mrle4_1 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    mrle4_2 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    mrle4_3 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    mrle4_4 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    dirty_wah = std::make_shared<GeneralModel>(256, 256, 18, 16, range_coder);
    mtype = std::make_shared<GeneralModel>(2, 512, 18, 1, range_coder);
}

void djn_ctx_model_t::InitiateNm() {
    range_coder = std::make_shared<RangeCoder>();
    mref = std::make_shared<GeneralModel>(2, 512, 18, 32, range_coder);
    mlog_rle = std::make_shared<GeneralModel>(64, 32, 18, 16, range_coder);  // 4 bit ref, log2(32) = 5 bit -> 2^9 = 512
    mrle = std::make_shared<GeneralModel>(256, 512, 24, 1, range_coder); // 4 bits ref + 5 bits log2(run) -> 2^9
    mrle2_1 = std::make_shared<GeneralModel>(256, 512, range_coder);
    mrle2_2 = std::make_shared<GeneralModel>(256, 512, range_coder);
    mrle4_1 = std::make_shared<GeneralModel>(256, 512, range_coder);
    mrle4_2 = std::make_shared<GeneralModel>(256, 512, range_coder);
    mrle4_3 = std::make_shared<GeneralModel>(256, 512, range_coder);
    mrle4_4 = std::make_shared<GeneralModel>(256, 512, range_coder);
    dirty_wah = std::make_shared<GeneralModel>(256, 256, 18, 32, range_coder);
    mtype = std::make_shared<GeneralModel>(2, 512, 18, 1, range_coder);
}

void djn_ctx_model_t::reset() {
    if (mref.get() != nullptr) mref->Reset();
    if (mlog_rle.get() != nullptr) mlog_rle->Reset();
    if (mrle.get() != nullptr) mrle->Reset();
    if (mrle2_1.get() != nullptr) mrle2_1->Reset(); 
    if (mrle2_2.get() != nullptr) mrle2_2->Reset();
    if (mrle4_1.get() != nullptr) mrle4_1->Reset(); 
    if (mrle4_2.get() != nullptr) mrle4_2->Reset(); 
    if (mrle4_3.get() != nullptr) mrle4_3->Reset(); 
    if (mrle4_4.get() != nullptr) mrle4_4->Reset();
    if (dirty_wah.get() != nullptr) dirty_wah->Reset();
    if (mtype.get() != nullptr) mtype->Reset();
    if (mref.get() != nullptr) mref->Reset();
    pbwt.Reset();
    n_variants = 0; p_len = 0;
}

int djn_ctx_model_t::StartEncoding(bool use_pbwt, bool reset) {
    if (range_coder.get() == nullptr) return -1;

    if (p_cap == 0) { // initiate a buffer if there is none
        delete[] p;
        p = new uint8_t[10000000];
        p_len = 0;
        p_cap = 10000000;
        p_free = true;
    }
    if (reset) this->reset();
    p_len = 0;

    range_coder->SetOutput(p);
    range_coder->StartEncode();
    return 1;
}

size_t djn_ctx_model_t::FinishEncoding() {
    if (range_coder.get() == nullptr) return -1;
    range_coder->FinishEncode();
    p_len = range_coder->OutSize();
    // std::cerr << "[djn_ctx_model_t::FinishEncoding] p_len=" << p_len << std::endl;
    return range_coder->OutSize();
}

int djn_ctx_model_t::StartDecoding(bool use_pbwt, bool reset) {
    if (range_coder.get() == nullptr) return -1;
    if (reset) this->reset();
    if (p == nullptr) return -2; // or result in corruption as range coder immediately loads data

    range_coder->SetInput(p);
    range_coder->StartDecode();
    return 1;
}

/*======   Varaint context model   ======*/

djinn_ctx_model::djinn_ctx_model() : 
    p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
    q(nullptr), q_len(0), q_alloc(0), q_free(true),
    range_coder(std::make_shared<RangeCoder>()), 
    ploidy_dict(std::make_shared<GeneralModel>(256, 256, range_coder))
{
    
}

djinn_ctx_model::~djinn_ctx_model() { 
    if (p_free) delete[] p;
    if (q_free) delete[] q;
}

int djinn_ctx_model::EncodeBcf(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) {
    if (data == nullptr) return -2;
    if (len_data % ploidy != 0) return -3;
    // Currently limited to 14 alt alleles + missing + EOV marker (total of 16).
    assert(alt_alleles < 14);

    // Todo: check that ploidy_dict model is initiated!
    std::shared_ptr<djn_ctx_model_container_t> tgt_container;

    const uint64_t tuple = ((uint64_t)len_data << 32) | ploidy;
    auto search = ploidy_map.find(tuple);
    if (search != ploidy_map.end()) {
        tgt_container = ploidy_models[search->second];
        ploidy_dict->EncodeSymbol(search->second);
    } else {
        // std::cerr << "Not found. Inserting: " << len_data << "," << ploidy << "(" << tuple << ") as " << ploidy_models.size() << std::endl;
        ploidy_map[tuple] = ploidy_models.size();
        ploidy_dict->EncodeSymbol(ploidy_models.size());
        ploidy_models.push_back(std::make_shared<djn_ctx_model_container_t>(len_data, ploidy, (bool)use_pbwt));
        tgt_container = ploidy_models[ploidy_models.size() - 1];
        tgt_container->StartEncoding(use_pbwt, init);
    }
    assert(tgt_container.get() != nullptr);

    // Compute allele counts.
    memset(hist_alts, 0, 256*sizeof(uint32_t));
    for (int i = 0; i < len_data; ++i) {
        ++hist_alts[BCF_UNPACK_GENOTYPE_GENERAL(data[i])];
    }
    // Check.
    assert(tgt_container->marchetype.get() != nullptr);

    // for (int i = 0; i < 256; ++i) {
        // if (hist_alts[i]) std::cerr << i << ":" << hist_alts[i] << ",";
    // }
    // std::cerr << " n_variants=" << tgt_container->n_variants << " data2m=" << tgt_container->model_2mc->range_coder->OutSize() << " dataNm=" << tgt_container->model_nm->range_coder->OutSize() << " permute=" << (int)use_pbwt << std::endl;

    // Biallelic, no missing, and no special EOV symbols.
    if (alt_alleles <= 2 && hist_alts[14] == 0 && hist_alts[15] == 0) {
        tgt_container->marchetype->EncodeSymbol(0); // add archtype as 2mc

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

            // std::cerr << "here 2mc: " << hist_alts[1] << std::endl;
            ret = (tgt_container->Encode2mc(tgt_container->model_2mc->pbwt.prev, len_data));
        } else {
            // std::cerr << "encoding nopbwt" << std::endl;
            ret = (tgt_container->Encode2mc(data, len_data, TWK_BCF_GT_UNPACK, 1));
        }

        if (ret > 0) ++n_variants;
        return ret;
    } else { // Otherwise.
        tgt_container->marchetype->EncodeSymbol(1);
        
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

int djinn_ctx_model::Encode(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) {
    if (data == nullptr) return -2;
    if (len_data % ploidy != 0) return -3;
    // Currently limited to 14 alt alleles + missing + EOV marker (total of 16).
    assert(alt_alleles < 14);

    // Todo: check that ploidy_dict model is initiated!
    std::shared_ptr<djn_ctx_model_container_t> tgt_container;

    const uint64_t tuple = ((uint64_t)len_data << 32) | ploidy;
    auto search = ploidy_map.find(tuple);
    if (search != ploidy_map.end()) {
        tgt_container = ploidy_models[search->second];
        ploidy_dict->EncodeSymbol(search->second);
    } else {
        // std::cerr << "Not found. Inserting: " << len_data << "," << ploidy << "(" << tuple << ") as " << ploidy_models.size() << std::endl;
        ploidy_map[tuple] = ploidy_models.size();
        ploidy_dict->EncodeSymbol(ploidy_models.size());
        ploidy_models.push_back(std::make_shared<djn_ctx_model_container_t>(len_data, ploidy, (bool)use_pbwt));
        tgt_container = ploidy_models[ploidy_models.size() - 1];
        tgt_container->StartEncoding(use_pbwt, init); // Todo: fix me
    }
    assert(tgt_container.get() != nullptr);

    // Compute allele counts.
    memset(hist_alts, 0, 256*sizeof(uint32_t));
    for (int i = 0; i < len_data; ++i) {
        ++hist_alts[data[i]];
    }
    // Check.
    assert(tgt_container->marchetype.get() != nullptr);

    // for (int i = 0; i < 256; ++i) {
    //     if (hist_alts[i]) std::cerr << i << ":" << hist_alts[i] << ",";
    // }
    // std::cerr << " data2m=" << tgt_container->model_2mc->range_coder->OutSize() << " dataNm=" << tgt_container->model_nm->range_coder->OutSize() << " permute=" << (bool)use_pbwt << std::endl;

    // Biallelic, no missing, and no special EOV symbols.
    if (alt_alleles <= 2 && hist_alts[14] == 0 && hist_alts[15] == 0) {
        tgt_container->marchetype->EncodeSymbol(0); // add archtype as 2mc

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
        tgt_container->marchetype->EncodeSymbol(1);
        
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

int djinn_ctx_model::DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len) {
    if (ewah_data == nullptr) return -1;
    if (ret_buffer == nullptr) return -1;

    // Decode stream archetype.
    uint8_t type = ploidy_dict->DecodeSymbol();
    // std::cerr << "[Decode model] Stream=" << (int)type << std::endl;
    std::shared_ptr<djn_ctx_model_container_t> tgt_container = ploidy_models[type];

    return(tgt_container->DecodeNext(ewah_data,ret_ewah,ret_buffer,ret_len));
}

int djinn_ctx_model::DecodeNext(djinn_variant_t*& variant) {
    // Decode stream archetype.
    uint8_t type = ploidy_dict->DecodeSymbol();
    // std::cerr << "[Decode model] Stream=" << (int)type << std::endl;
    std::shared_ptr<djn_ctx_model_container_t> tgt_container = ploidy_models[type];

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

    variant->ploidy = tgt_container->ploidy;
    variant->data_len = 0;
    variant->errcode = 0;

    return(tgt_container->DecodeNext(q,q_len,variant->data,variant->data_len));
}

int djinn_ctx_model::DecodeNextRaw(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;
    uint8_t type = ploidy_dict->DecodeSymbol();

    std::shared_ptr<djn_ctx_model_container_t> tgt_container = ploidy_models[type];
    return(tgt_container->DecodeNextRaw(data, len));
}

void djinn_ctx_model::StartEncoding(bool use_pbwt, bool reset) {
    // Store parameters for decoding.
    this->use_pbwt = use_pbwt;
    this->init = reset;
    
    // If resetting the model_2mc->
    if (reset) {
        // std::cerr << "[djinn_ctx_model::StartEncoding] resetting" << std::endl;
        for (int i = 0; i < ploidy_models.size(); ++i) {
            ploidy_models[i]->model_2mc->reset();
            ploidy_models[i]->model_nm->reset();
        }
    } else {
        for (int i = 0; i < ploidy_models.size(); ++i) {
            ploidy_models[i]->model_2mc->n_variants = 0;
            ploidy_models[i]->model_nm->n_variants  = 0;
        }
    }

    n_variants = 0;
    p_len = 0;

    // Local range coder
    range_coder->SetOutput(p);
    range_coder->StartEncode();

    // std::cerr << "[djinn_ctx_model::StartEncoding] models start encoding" << std::endl;
    for (int i = 0; i < ploidy_models.size(); ++i) {
        ploidy_models[i]->StartEncoding(use_pbwt, reset);
    }
}

size_t djinn_ctx_model::FinishEncoding() {
    if (range_coder.get() == nullptr) return -1;
    range_coder->FinishEncode();
    p_len = range_coder->OutSize();

    // std::cerr << "[djinn_ctx_model::FinishEncoding] p_len=" << p_len << std::endl;

    size_t s_models = 0;
    for (int i = 0; i < ploidy_models.size(); ++i) {
        s_models += ploidy_models[i]->FinishEncoding();
        // std::cerr << "Model-" << i << ": " << ploidy_models[i]->range_coder->OutSize() << std::endl;
    }

    size_t s_rc  = range_coder->OutSize();
    // std::cerr << "Model sizes=" << s_rc << " and models=" << s_models << std::endl;
    
    return s_rc + s_models;
}

int djinn_ctx_model::StartDecoding() {
    assert(p != nullptr);
    assert(ploidy_models.size() != 0);

    for (int i = 0; i <ploidy_models.size(); ++i) {
        ploidy_models[i]->StartDecoding(use_pbwt, init);
    }

    range_coder->SetInput(p);
    range_coder->StartDecode();
    return 1;
}

// Read/write
int djinn_ctx_model::Serialize(uint8_t* dst) const {
    // Serialize as (int,uint32_t,uint32_t,uint8_t*,ctx1,ctx2):
    // #models,p_len,p,[models...]
    uint32_t offset = 0;

    // std::cerr << "[djinn_ctx_model::Serialize] p_len=" << p_len << std::endl;

    // Reserve space for offset
    offset += sizeof(uint32_t);

    // Add
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

    // std::cerr << "[Serialize] Offset here=" << offset << std::endl;

    // Serialize each model.
    for (int i = 0; i < ploidy_models.size(); ++i) {
        offset += ploidy_models[i]->Serialize(&dst[offset]);    
    }
    
    *((uint32_t*)&dst[0]) = offset;
    // std::cerr << "test=" << GetSerializedSize() << "/" << offset << std::endl;
    // assert(GetSerializedSize() == offset);
    
    return offset;
}

int djinn_ctx_model::Serialize(std::ostream& stream) const {
    // Serialize as (int,uint32_t,uint32_t,uint8_t*,ctx1,ctx2):
    // #models,p_len,p,[models...]
    uint32_t out_len = GetSerializedSize();
    stream.write((char*)&out_len, sizeof(uint32_t));

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

int djinn_ctx_model::GetSerializedSize() const {
    int ret = sizeof(uint32_t) + sizeof(int) + sizeof(uint32_t) + sizeof(uint8_t) + sizeof(uint32_t) + p_len;
    for (int i = 0; i < ploidy_models.size(); ++i) {
        ret += ploidy_models[i]->GetSerializedSize();
    }
    return ret;
}

int djinn_ctx_model::GetCurrentSize() const {
    int ret = range_coder->OutSize();
    for (int i = 0; i < ploidy_models.size(); ++i) {
        ret += ploidy_models[i]->GetCurrentSize();
    }
    return ret;
}

// Deserialize data from an external buffer.
int djinn_ctx_model::Deserialize(uint8_t* src) {
    // Reset.
    ploidy_remap.clear();

    uint32_t offset = 0;

    // Read total offset.
    uint32_t tot_offset = *((uint32_t*)&src[offset]);
    offset += sizeof(uint32_t);

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
    // Store model selection data.
    memcpy(p, &src[offset], p_len);
    offset += p_len;
    // Deserialize each model.
    // std::cerr << "[Deserialize] Offset here=" << offset << " #models=" << n_models << std::endl;
    // std::cerr << "[Deserialize] tot_offset=" << tot_offset << ",n_variants=" << n_variants << ",use_pbwt=" << (int)use_pbwt << ",init=" << (bool)init << ",p_len=" << p_len << std::endl;
    
    // Read each model.
    for (int i = 0; i < n_models; ++i) {
        // First peek at their content to invoke the correct constructor as it
        // requires #samples, #ploidy, use_pbwt
        int pl = *((int*)&src[offset]);
        uint32_t n_s = *((uint32_t*)&src[offset+sizeof(int)]);
        // std::cerr << "[Deserialize] #pl=" << pl << ", #n_s=" << n_s << std::endl;
        
        // Update ploidy map and insert data in the correct position relative
        // the encoding order.
        const uint64_t tuple = ((uint64_t)n_s << 32) | pl;
        auto search = ploidy_map.find(tuple);
        if (search != ploidy_map.end()) {
            // std::cerr << "Illegal! Cannot exist! Corruption..." << std::endl;
            // exit(1);
            // std::cerr << "Already exist: " << search->second << std::endl;
            // ploidy_remap[search->second] = search->second;
            offset += ploidy_models[search->second]->Deserialize(&src[offset]);
        } else {
            // std::cerr << "[Deserialize] Adding map [" << tuple << "] for [" << n_s << "," << pl << "]" << std::endl; 
            ploidy_map[tuple] = ploidy_models.size();
            ploidy_models.push_back(std::make_shared<djinn::djn_ctx_model_container_t>(n_s, pl, (bool)use_pbwt));
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
int djinn_ctx_model::Deserialize(std::istream& stream) {
    // uint32_t out_len = GetSerializedSize();
    uint32_t out_len = 0;
    stream.read((char*)&out_len, sizeof(uint32_t));

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
    stream.read((char*)p, p_len);

    // Serialize each model.
    for (int i = 0; i < n_models; ++i) {
        // First peek at their content to invoke the correct constructor as it
        // requires #samples, #ploidy, use_pbwt
        int pl = 0; 
        uint32_t n_s = 0; 
        stream.read((char*)&pl, sizeof(int));
        stream.read((char*)&n_s, sizeof(uint32_t));

        // std::cerr << "[Deserialize] #pl=" << pl << ", #n_s=" << n_s << std::endl;

        // Update ploidy map and insert data in the correct position relative
        // the encoding order.
        const uint64_t tuple = ((uint64_t)n_s << 32) | pl;
        auto search = ploidy_map.find(tuple);
        if (search != ploidy_map.end()) {
            // std::cerr << "Illegal! Cannot exist! Corruption..." << std::endl;
            // exit(1);
            // std::cerr << "Already exist: " << search->second << std::endl;
            // ploidy_remap[search->second] = search->second;
            ploidy_models[search->second]->Deserialize(stream);
        } else {
            // std::cerr << "[Deserialize] Adding map [" << tuple << "] for [" << n_s << "," << pl << "]" << std::endl; 
            ploidy_map[tuple] = ploidy_models.size();
            ploidy_models.push_back(std::make_shared<djinn::djn_ctx_model_container_t>(n_s, pl, (bool)use_pbwt));
            ploidy_models.back()->Deserialize(stream);
        }
    }
    // std::cerr << "[Deserialize] Decoded=" << offset << "/" << tot_offset << std::endl;

    return stream.tellg();
}

/*======   Container   ======*/

djn_ctx_model_container_t::djn_ctx_model_container_t(int64_t n_s, int pl, bool use_pbwt) : 
    use_pbwt(use_pbwt),
    ploidy(pl), n_samples(n_s), n_variants(0),
    n_samples_wah(std::ceil((float)n_samples / 32) * 32), 
    n_samples_wah_nm(std::ceil((float)n_samples * 4/32) * 8),
    n_wah(n_samples_wah_nm / 8),
    wah_bitmaps(new uint32_t[n_wah]), 
    p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
    range_coder(std::make_shared<RangeCoder>()), 
    marchetype(std::make_shared<GeneralModel>(2, 1024, range_coder)),
    model_2mc(std::make_shared<djn_ctx_model_t>()),
    model_nm(std::make_shared<djn_ctx_model_t>())
{
    assert(n_s % pl == 0); // #samples/#ploidy must be divisible
    model_2mc->Initiate2mc();
    model_nm->InitiateNm();
}

djn_ctx_model_container_t::djn_ctx_model_container_t(int64_t n_s, int pl, bool use_pbwt, uint8_t* src, uint32_t src_len) : 
    use_pbwt(use_pbwt),
    ploidy(pl), n_samples(n_s), n_variants(0),
    n_samples_wah(std::ceil((float)n_samples / 32) * 32), 
    n_samples_wah_nm(std::ceil((float)n_samples * 4/32) * 8),
    n_wah(n_samples_wah_nm / 8),
    wah_bitmaps(new uint32_t[n_wah]), 
    p(src), p_len(src_len), p_cap(0), p_free(false),
    range_coder(std::make_shared<RangeCoder>()), 
    marchetype(std::make_shared<GeneralModel>(2, 1024, range_coder)),
    model_2mc(std::make_shared<djn_ctx_model_t>()),
    model_nm(std::make_shared<djn_ctx_model_t>())
{
    assert(n_s % pl == 0); // #samples/#ploidy must be divisible
    model_2mc->Initiate2mc();
    model_nm->InitiateNm();
}

djn_ctx_model_container_t::~djn_ctx_model_container_t() {
    if (p_free) delete[] p;
    delete[] wah_bitmaps;
}

void djn_ctx_model_container_t::StartEncoding(bool use_pbwt, bool reset) {
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
        // std::cerr << "[djn_ctx_model_container_t::StartEncoding] resetting" << std::endl;
        model_2mc->reset();
        model_nm->reset();
    } else {
        model_2mc->n_variants = 0;
        model_nm->n_variants = 0;
    }

    this->use_pbwt = use_pbwt;
    n_variants = 0;

    // Local range coder
    if (range_coder.get() == nullptr)
        range_coder = std::make_shared<RangeCoder>();
    range_coder->SetOutput(p);
    range_coder->StartEncode();

    model_2mc->StartEncoding(use_pbwt, reset);
    model_nm->StartEncoding(use_pbwt, reset);
}

size_t djn_ctx_model_container_t::FinishEncoding() {
    if (range_coder.get() == nullptr) return -1;
    range_coder->FinishEncode();
    model_2mc->FinishEncoding();
    model_nm->FinishEncoding();
    p_len = range_coder->OutSize();

    size_t s_rc  = range_coder->OutSize();
    size_t s_2mc = model_2mc->range_coder->OutSize();
    size_t s_nm  = model_nm->range_coder->OutSize();
    std::cerr << "container finish. type=" << s_rc << " 2mc=" << s_2mc << " nm=" << s_nm << std::endl; 
    
    return s_rc + s_2mc + s_nm;
}

void djn_ctx_model_container_t::StartDecoding(bool use_pbwt, bool reset) {
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
        // std::cerr << "[djn_ctx_model_container_t::StartDecoding] resetting" << std::endl;
        model_2mc->reset();
        model_nm->reset();
    } else {
        model_2mc->n_variants = 0;
        model_nm->n_variants = 0;
    }
    // std::cerr << "[djn_ctx_model_container_t::StartDecoding] p_len=" << p_len << std::endl;

    this->use_pbwt = use_pbwt;

    // Local range coder
    if (range_coder.get() == nullptr)
        range_coder = std::make_shared<RangeCoder>();
    range_coder->SetInput(p);
    range_coder->StartDecode();

    model_2mc->StartDecoding(use_pbwt, reset);
    model_nm->StartDecoding(use_pbwt, reset);
}

int djn_ctx_model_container_t::Encode2mc(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;

    if (wah_bitmaps == nullptr) {
        n_samples_wah = std::ceil((float)n_samples / 32) * 32;
        n_samples_wah_nm = std::ceil((float)n_samples * 4/32) * 8; // Can fit eight 4-bit entries in a 32-bit word
        n_wah = n_samples_wah_nm / 8; 
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    for (int i = 0; i < n_samples; ++i) {
        if (data[i]) {
            wah_bitmaps[i / 32] |= 1L << (i % 32);
        }
    }

    return EncodeWah(wah_bitmaps, n_samples_wah >> 5); // n_samples_wah / 32
}

int djn_ctx_model_container_t::Encode2mc(uint8_t* data, uint32_t len, const uint8_t* map, const int shift) {
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
    // std::cerr << "encoding 2mc=" << n_wah/4 << "->" << n_wah/4*32 << std::endl;

    return EncodeWah(wah_bitmaps, n_samples_wah >> 5); // n_samples_wah / 32
}

int djn_ctx_model_container_t::EncodeNm(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;
    
    if (wah_bitmaps == nullptr) {
        n_samples_wah = std::ceil((float)n_samples / 32) * 32;
        n_samples_wah_nm = std::ceil((float)n_samples * 4/32) * 8; // Can fit eight 4-bit entries in a 32-bit word
        n_wah = n_samples_wah_nm / 8; 
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    // Todo: move into ResetBitmaps() functions
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    for (int i = 0; i < n_samples; ++i) {
        wah_bitmaps[i / 8] |= (uint64_t)data[i] << (4*(i % 8));
    }

    for (int i = 0; i < (n_samples_wah_nm >> 3); ++i) { 
        if (wah_bitmaps[i] != 0) std::cerr << i << ":" << std::bitset<32>(wah_bitmaps[i]) << std::endl;
    }
    // std::cerr << std::endl;

    std::cerr << "[djn_ctx_model_container_t::EncodeNm] n_wah=" << n_wah << " and input=" << (n_samples_wah_nm >> 3) << std::endl;
    // exit(1);

    return EncodeWahNm(wah_bitmaps, n_samples_wah_nm >> 3); // n_samples_wah_nm / 8
}

int djn_ctx_model_container_t::EncodeNm(uint8_t* data, uint32_t len, const uint8_t* map, const int shift) {
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

    return EncodeWahNm(wah_bitmaps, n_samples_wah_nm >> 3); // n_samples_wah_nm / 8
}

int djn_ctx_model_container_t::DecodeRaw_nm(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;

    int64_t n_samples_obs = 0;
    int objects = 0;

    // std::cerr << "[djn_ctx_model_container_t::DecodeRaw_nm] start pos=" << len << std::endl;

    // Emit empty EWAH marker.
    djinn_ewah_t* ewah = (djinn_ewah_t*)&data[len]; 
    ewah->reset();
    len += sizeof(djinn_ewah_t);

    while(true) {
        uint8_t type = model_nm->mtype->DecodeSymbol();
        // data[len++] = type; // store archetype

        // std::cerr << "Type=" << (int)type << std::endl;
        if (type == 0) { // bitmaps
            ++ewah->dirty;
            // uint64_t wah = 0;
            // std::cerr << "Bitmap=";
            uint32_t* t = (uint32_t*)&data[len];
            for (int i = 0; i < 4; ++i) {
                data[len++] = model_nm->dirty_wah->DecodeSymbol();
                // wah |= model_nm->dirty_wah->DecodeSymbol();
                // std::cerr << std::bitset<8>(wah&255);
                // wah <<= 8;
            }
            // std::cerr << std::endl;
            // std::cerr << "[djn_ctx_model_container_t::DecodeRaw_nm] dirty " << std::bitset<32>(*t) << std::endl;
            n_samples_obs += 8;

        } else { // is RLE
            // Emit new EWAH marker
            if (ewah->clean > 0 || ewah->dirty > 0) {
                std::cerr << "[djn_ctx_model_container_t::DecodeRaw_nm] clean=" << ewah->clean << ", dirty=" << ewah->dirty << ", ref=" << ewah->ref << std::endl;
                ewah = (djinn_ewah_t*)&data[len];
                ewah->reset();
                len += sizeof(djinn_ewah_t);
                ++objects;
            }

            // Decode an RLE
            uint32_t ref = 1; uint32_t len = 1;
            DecodeWahRLE_nm(ref, len, model_nm);
            // std::cerr << "decoded RLE=" << ref << " len=" << len << std::endl;
            ewah->ref = ref & 15;
            ewah->clean = len;

            n_samples_obs += ewah->clean*8;
        }

        if (n_samples_obs == n_samples_wah_nm) {
            // std::cerr << "[djn_ctx_model_container_t::DecodeRaw_nm] obs=" << n_samples_obs << "/" << n_samples_wah_nm << std::endl;
            break;
        }
        if (n_samples_obs > n_samples_wah_nm) {
            std::cerr << "[djn_ctx_model_container_t::DecodeRaw_nm] Decompression corruption: " << n_samples_obs << "/" << n_samples_wah_nm  << std::endl;
            exit(1);
        }
    }

    if (ewah->clean > 0 || ewah->dirty > 0) {
        ++objects;
        std::cerr << "[djn_ctx_model_container_t::DecodeRaw_nm] clean=" << ewah->clean << ", dirty=" << ewah->dirty << ", ref=" << ewah->ref << std::endl;
        // std::cerr << "final=" << std::bitset<64>(*ewah) << " ref=" << ewah_ref << ",run=" << ewah_run << ",dirty=" << ewah_dirty << std::endl;
    }
    // std::cerr << "n_objs=" << objects << std::endl;

    return 1;
}

int djn_ctx_model_container_t::EncodeWah(uint32_t* wah, uint32_t len) { // input WAH-encoded data
    if (wah == nullptr) return -1;
    ++model_2mc->n_variants;

    uint32_t wah_ref = wah[0];
    uint32_t wah_run = 1;

    // Resize if necessary.
    // std::cerr << "[ENCODE] " << model_2mc->range_coder->OutSize() << "/" << model_2mc->p_cap << std::endl;
    if (model_2mc->range_coder->OutSize() + n_samples > model_2mc->p_cap) {
        const uint32_t rc_size = model_2mc->range_coder->OutSize();
        std::cerr << "[djn_ctx_model_container_t::EncodeWah][RESIZE] resizing from: " << rc_size << "->" << (rc_size + 2*n_samples + 65536) << std::endl;
        uint8_t* prev = model_2mc->p; // old
        model_2mc->p_cap = rc_size + 2*n_samples + 65536;
        model_2mc->p = new uint8_t[model_2mc->p_cap];
        memcpy(model_2mc->p, prev, rc_size);
        if (model_2mc->p_free) delete[] prev;
        model_2mc->p_free = true;
        model_2mc->range_coder->out_buf = model_2mc->p + rc_size; 
        model_2mc->range_coder->in_buf  = model_2mc->p;
    }

    for (int i = 1; i < len; ++i) {
        if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint32_t>::max()) || (wah_ref != wah[i])) {
            if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint32_t>::max()) || wah_run == 1) {
                model_2mc->mtype->EncodeSymbol(0);
                
                // std::cerr << "Bitmap=" << std::bitset<32>(wah_ref) << ": " << __builtin_popcount(wah_ref) << std::endl;
                for (int i = 0; i < 4; ++i) {
                    model_2mc->dirty_wah->EncodeSymbol(wah_ref & 255);
                    // std::cerr << "bitmap-" << i << " " << std::bitset<32>(wah_ref * 255) << " ";
                    // ref_alt[1] += __builtin_popcount(wah_ref & 255);
                    wah_ref >>= 8;
                }
                // std::cerr << std::endl;
            } else {
                model_2mc->mtype->EncodeSymbol(1);
                // std::cerr << "wah_ref: " << (int)(wah_ref&1) << " wah_run=" << wah_run << std::endl;
                // ref_alt[wah_ref & 1] += wah_run * 32;
                EncodeWahRLE(wah_ref, wah_run, model_2mc);
            }
            wah_run = 0;
            wah_ref = wah[i];
            // wah_ref_pos = i;
        }
        ++wah_run;
    }

    if (wah_run) {
        // std::cerr << "wah run left=" << wah_run << std::endl;
        if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint32_t>::max()) || wah_run == 1) {
            model_2mc->mtype->EncodeSymbol(0);
            
            // std::cerr << "Bitmap=" << std::bitset<32>(wah_ref) << ": " << __builtin_popcount(wah_ref) << std::endl;
            for (int i = 0; i < 4; ++i) {
                model_2mc->dirty_wah->EncodeSymbol(wah_ref & 255);
                // ref_alt[1] += __builtin_popcount(wah_ref & 255);
                wah_ref >>= 8;
            }
            
        } else {
            model_2mc->mtype->EncodeSymbol(1);
            // ref_alt[wah_ref & 1] += wah_run * 32;
            // std::cerr << "wah_ref: " << (int)(wah_ref&1) << " wah_run=" << wah_run << std::endl;
            EncodeWahRLE(wah_ref, wah_run, model_2mc);
        }
    }
    // std::cerr << "encoding wah ref_alt=" << ref_alt[0] << "," << ref_alt[1] << std::endl;

    ++n_variants;
    return 1;
}

int djn_ctx_model_container_t::EncodeWahNm(uint32_t* wah, uint32_t len) { // input WAH-encoded data
    if (wah == nullptr) return -1;

    ++model_nm->n_variants;

    // std::cerr << "refnm=" << std::bitset<32>(wah[0]) << std::endl;

    // uint32_t wah_ref = wah[0];
    // uint64_t wah_run = 1;

    // Resize if necessary.
    // std::cerr << "[ENCODE] " << model_nm->range_coder->OutSize() << "/" << model_nm->p_cap << std::endl;
    if (model_nm->range_coder->OutSize() + n_samples > model_nm->p_cap) {
        const uint32_t rc_size = model_nm->range_coder->OutSize();
        std::cerr << "[djn_ctx_model_container_t::EncodeWahNm][RESIZE] resizing from: " << rc_size << "->" << (rc_size + 2*n_samples + 65536) << std::endl;
        uint8_t* prev = model_nm->p; // old
        model_nm->p_cap = rc_size + 2*n_samples + 65536;
        model_nm->p = new uint8_t[model_nm->p_cap];
        memcpy(model_nm->p, prev, rc_size);
        if (model_nm->p_free) delete[] prev;
        model_nm->p_free = true;
        model_nm->range_coder->out_buf = model_nm->p + rc_size; 
        model_nm->range_coder->in_buf  = model_nm->p;
    }
    
    djinn_ewah_t ewah; ewah.reset();
    
    // Build reference by broadcasting lower 4 bits
    // to all 8 four-bit positions in a 32-bit word.
    uint32_t wah_ref = (wah[0] & 15);
    uint32_t wah_cmp = wah_ref;
    for (int i = 1; i < 8; ++i) wah_cmp |= wah_ref << (i*4);

    // Word is dirty if the word is different from the expected clean word.
    ewah.dirty += (wah_ref != wah_cmp);
    if (ewah.dirty) {
        std::cerr << "First word is dirty!" << std::endl;
        model_nm->mtype->EncodeSymbol(0);
        for (int i = 0; i < 4; ++i) {
            model_nm->dirty_wah->EncodeSymbol(wah_ref & 255);
            wah_ref >>= 8;
        }
    }

    // Word is clean: pattern is repeated.
    ewah.clean += (wah_ref == wah_cmp);
    if (ewah.clean) ewah.ref = (wah[0] & 15);
    assert(ewah.dirty + ewah.clean == 1);

    uint32_t n_objs = 1;
    uint32_t n_obs = 0;

    for (int i = 1; i < len; ++i) {
        // Build reference.
        wah_ref = (wah[i] & 15);
        wah_cmp = wah_ref;
        for (int j = 1; j < 8; ++j) wah_cmp |= wah_ref << (j*4);

        // if (i == 136 || i == 338 || i == 386 || i == 522) {
        //     std::cerr << "Debug@" << i << ": " << wah_ref << "," << wah_cmp << "," << std::bitset<32>(wah[i]) << std::endl;
        //     std::cerr << "Equality test=" << (wah[i] == wah_cmp) << std::endl;
        // }

        // Is dirty
        if (wah[i] != wah_cmp) {
            if (ewah.clean) {
                model_nm->mtype->EncodeSymbol(1);
                EncodeWahRLE_nm(ewah.ref, ewah.clean, model_nm);
                n_obs += ewah.clean;
                ewah.reset();
            }

            ++ewah.dirty;
            model_nm->mtype->EncodeSymbol(0);
            
            wah_ref = wah[i];
            // std::cerr << "[djn_ctx_model_container_t::EncodeWahNm] Dirty: " << i << "/" << len << " " << std::bitset<32>(wah_ref) << std::endl;
            for (int i = 0; i < 4; ++i) {
                model_nm->dirty_wah->EncodeSymbol(wah_ref & 255);
                wah_ref >>= 8;
            }
            ++n_obs;
        } 
        // Is clean
        else {
            // Dirty have been set then make new EWAH
            if (ewah.dirty) {
                if (ewah.clean) {
                    model_nm->mtype->EncodeSymbol(1);
                    EncodeWahRLE_nm(ewah.ref, ewah.clean, model_nm);
                    n_obs += ewah.clean;
                }
                // n_obs += ewah.run;

                ewah.reset();
                ++ewah.clean;
                ewah.ref = wah_ref;
                ++n_objs;
            } 
            // No dirty have been set
            else {
                if (ewah.clean == 0) ewah.ref = wah_ref;
                if (ewah.ref == wah_ref) ++ewah.clean;
                else { // Different clean word
                    if (ewah.clean) {
                        model_nm->mtype->EncodeSymbol(1);
                        EncodeWahRLE_nm(ewah.ref, ewah.clean, model_nm);
                        n_obs += ewah.clean;
                        // n_obs += wah_run;
                    }

                    ewah.reset();
                    ++ewah.clean;
                    ewah.ref = wah_ref;
                    ++n_objs;
                }
            }
        }
    }

    if (ewah.dirty || ewah.clean) {
        // std::cerr << "[djn_ctx_model_container_t::EncodeWahNm] Residuals = " << ewah.clean << "," << ewah.dirty << std::endl;    
        if (ewah.clean) {
            model_nm->mtype->EncodeSymbol(1);
            EncodeWahRLE_nm(ewah.ref, ewah.clean, model_nm);
            n_obs += ewah.clean;
            ++n_objs;
        }
    }

    std::cerr << "[djn_ctx_model_container_t::EncodeWahNm] Done. Total = " << n_obs << "/" << len << std::endl;
    assert(n_obs == len);

    ++n_variants;
    return 1;
}

int djn_ctx_model_container_t::EncodeWahRLE(uint32_t ref, uint32_t len, std::shared_ptr<djn_ctx_model_t> model) {
    // Length cannot be 0 so remove 1
    // len -= 1;

    model->mref->EncodeSymbol(ref&1);
    uint32_t log_length = round_log2(len);
    model->mlog_rle->EncodeSymbol(log_length);

    // std::cerr << "RLE=" << (ref&1) << ",len=" << len << ",log=" << log_length << std::endl;

    if (log_length < 2) {
        // std::cerr << "single=" << len << "," << (ref&1) << std::endl;
    }
    else if (log_length <= 8) {
        assert(len < 256);
        model->mrle->model_context <<= 1;
        model->mrle->model_context |= (ref & 1);
        model->mrle->model_context <<= 5;
        model->mrle->model_context |= log_length;
        model->mrle->model_context &= model->mrle->model_ctx_mask;
        model->mrle->EncodeSymbolNoUpdate(len & 255);

    } else if (log_length <= 16) {
        assert(len < 65536);
        model->mrle2_1->model_context <<= 1;
        model->mrle2_1->model_context |= (ref & 1);
        model->mrle2_1->model_context <<= 5;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        model->mrle2_1->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle2_2->model_context <<= 1;
        model->mrle2_2->model_context |= (ref & 1);
        model->mrle2_2->model_context <<= 5;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        model->mrle2_2->EncodeSymbolNoUpdate(len & 255);
    } else {
        model->mrle4_1->model_context <<= 1;
        model->mrle4_1->model_context |= (ref & 1);
        model->mrle4_1->model_context <<= 5;
        model->mrle4_1->model_context |= log_length;
        model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
        model->mrle4_1->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_2->model_context <<= 1;
        model->mrle4_2->model_context |= (ref & 1);
        model->mrle4_2->model_context <<= 5;
        model->mrle4_2->model_context |= log_length;
        model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
        model->mrle4_2->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_3->model_context <<= 1;
        model->mrle4_3->model_context |= (ref & 1);
        model->mrle4_3->model_context <<= 5;
        model->mrle4_3->model_context |= log_length;
        model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
        model->mrle4_3->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_4->model_context <<= 1;
        model->mrle4_4->model_context |= (ref & 1);
        model->mrle4_4->model_context <<= 5;
        model->mrle4_4->model_context |= log_length;
        model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
        model->mrle4_4->EncodeSymbolNoUpdate(len & 255);
    }

    return 1;
}

int djn_ctx_model_container_t::EncodeWahRLE_nm(uint32_t ref, uint32_t len, std::shared_ptr<djn_ctx_model_t> model) {
    model->mref->EncodeSymbol(ref & 15);
    uint32_t log_length = round_log2(len);
    model->mlog_rle->EncodeSymbol(log_length);

    // std::cerr << "nmRLE=" << (ref&15) << ",len=" << len << ",log=" << log_length << std::endl;

    if (log_length < 2) {
        std::cerr << "single=" << len << "," << (ref & 15) << std::endl;
    }
    else if (log_length <= 8) {
        assert(len < 256);
        model->mrle->model_context = (ref & 15);
        model->mrle->model_context <<= 5;
        model->mrle->model_context |= log_length;
        model->mrle->model_context &= model->mrle->model_ctx_mask;
        model->mrle->EncodeSymbolNoUpdate(len & 255);

    } else if (log_length <= 16) {
        assert(len < 65536);
        model->mrle2_1->model_context = (ref & 15);
        model->mrle2_1->model_context <<= 5;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        model->mrle2_1->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle2_2->model_context = (ref & 15);
        model->mrle2_2->model_context <<= 5;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        model->mrle2_2->EncodeSymbolNoUpdate(len & 255);
    } else {
        model->mrle4_1->model_context = (ref & 15);
        model->mrle4_1->model_context <<= 5;
        model->mrle4_1->model_context |= log_length;
        model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
        model->mrle4_1->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_2->model_context = (ref & 15);
        model->mrle4_2->model_context <<= 5;
        model->mrle4_2->model_context |= log_length;
        model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
        model->mrle4_2->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_3->model_context = (ref & 15);
        model->mrle4_3->model_context <<= 5;
        model->mrle4_3->model_context |= log_length;
        model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
        model->mrle4_3->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_4->model_context = (ref & 15);
        model->mrle4_4->model_context <<= 5;
        model->mrle4_4->model_context |= log_length;
        model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
        model->mrle4_4->EncodeSymbolNoUpdate(len & 255);
    }

    return 1;
}

int djn_ctx_model_container_t::DecodeWahRLE(uint32_t& ref, uint32_t& len, std::shared_ptr<djn_ctx_model_t> model) {
    ref = model->mref->DecodeSymbol();
    uint32_t log_length = model->mlog_rle->DecodeSymbol();
    // std::cerr << "ref=" << ref << " log=" << log_length << std::endl;

    if (log_length < 2) {
        std::cerr << "single=" << len << "," << (ref&1) << std::endl;
    }
    else if (log_length <= 8) {
        model->mrle->model_context <<= 1;
        model->mrle->model_context |= (ref & 1);
        model->mrle->model_context <<= 5;
        model->mrle->model_context |= log_length;
        model->mrle->model_context &= model->mrle->model_ctx_mask;
        len = model->mrle->DecodeSymbolNoUpdate();
    } else if (log_length <= 16) {
        model->mrle2_1->model_context <<= 1;
        model->mrle2_1->model_context |= (ref & 1);
        model->mrle2_1->model_context <<= 5;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        len = model->mrle2_1->DecodeSymbolNoUpdate();
        model->mrle2_2->model_context <<= 1;
        model->mrle2_2->model_context |= (ref & 1);
        model->mrle2_2->model_context <<= 5;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        len |= (uint32_t)model->mrle2_2->DecodeSymbolNoUpdate() << 8;
    } else {
        model->mrle4_1->model_context <<= 1;
        model->mrle4_1->model_context |= (ref & 1);
        model->mrle4_1->model_context <<= 5;
        model->mrle4_1->model_context |= log_length;
        model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
        len = model->mrle4_1->DecodeSymbolNoUpdate();
        model->mrle4_2->model_context <<= 1;
        model->mrle4_2->model_context |= (ref & 1);
        model->mrle4_2->model_context <<= 5;
        model->mrle4_2->model_context |= log_length;
        model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
        len |= (uint32_t)model->mrle4_2->DecodeSymbolNoUpdate() << 8;
        model->mrle4_3->model_context <<= 1;
        model->mrle4_3->model_context |= (ref & 1);
        model->mrle4_3->model_context <<= 5;
        model->mrle4_3->model_context |= log_length;
        model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
        len |= (uint32_t)model->mrle4_3->DecodeSymbolNoUpdate() << 16;
        model->mrle4_4->model_context <<= 1;
        model->mrle4_4->model_context |= (ref & 1);
        model->mrle4_4->model_context <<= 5;
        model->mrle4_4->model_context |= log_length;
        model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
        len |= (uint32_t)model->mrle4_4->DecodeSymbolNoUpdate() << 24;
    }

    return 1;
}

int djn_ctx_model_container_t::DecodeWahRLE_nm(uint32_t& ref, uint32_t& len, std::shared_ptr<djn_ctx_model_t> model) {
    ref = model->mref->DecodeSymbol();
    uint32_t log_length = model->mlog_rle->DecodeSymbol();
    // std::cerr << "ref=" << ref << " log=" << log_length << std::endl;

    if (log_length < 2) {
        std::cerr << "single=" << len << "," << (ref&15) << std::endl;
    }
    else if (log_length <= 8) {
        model->mrle->model_context = (ref & 15);
        model->mrle->model_context <<= 5;
        model->mrle->model_context |= log_length;
        model->mrle->model_context &= model->mrle->model_ctx_mask;
        len = model->mrle->DecodeSymbolNoUpdate();
    } else if (log_length <= 16) {
        model->mrle2_1->model_context = (ref & 15);
        model->mrle2_1->model_context <<= 5;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        len = model->mrle2_1->DecodeSymbolNoUpdate();
        model->mrle2_2->model_context = (ref & 15);
        model->mrle2_2->model_context <<= 5;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        len |= (uint32_t)model->mrle2_2->DecodeSymbolNoUpdate() << 8;
    } else {
        model->mrle4_1->model_context = (ref & 15);
        model->mrle4_1->model_context <<= 5;
        model->mrle4_1->model_context |= log_length;
        model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
        len = model->mrle4_1->DecodeSymbolNoUpdate();
        model->mrle4_2->model_context = (ref & 15);
        model->mrle4_2->model_context <<= 5;
        model->mrle4_2->model_context |= log_length;
        model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
        len |= (uint32_t)model->mrle4_2->DecodeSymbolNoUpdate() << 8;
        model->mrle4_3->model_context = (ref & 15);
        model->mrle4_3->model_context <<= 5;
        model->mrle4_3->model_context |= log_length;
        model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
        len |= (uint32_t)model->mrle4_3->DecodeSymbolNoUpdate() << 16;
        model->mrle4_4->model_context = (ref & 15);
        model->mrle4_4->model_context <<= 5;
        model->mrle4_4->model_context |= log_length;
        model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
        len |= (uint32_t)model->mrle4_4->DecodeSymbolNoUpdate() << 24;
    }

    return 1;
}

int djn_ctx_model_container_t::DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len) {
    if (ewah_data == nullptr) return -1;
    if (ret_buffer == nullptr) return -2;

    // Decode stream archetype.
    uint8_t type = marchetype->DecodeSymbol();
    // std::cerr << "[djn_ctx_model_container_t::DecodeNext] Type=" << (int)type << std::endl;

    size_t ret_ewah_init = ret_ewah;
    int objs = 0;
    switch(type) {
    case 0: objs = DecodeRaw(ewah_data, ret_ewah); break;
    case 1: objs = DecodeRaw_nm(ewah_data, ret_ewah); break;
    default: std::cerr << "[djn_ctx_model_container_t::DecodeNext] decoding error: " << (int)type << " (valid=[0,1])" << std::endl; return -1;
    }

    if (objs <= 0) return -1;

    // std::cerr << "alts=" << hist_alts[0] << "," << hist_alts[1] << std::endl;

    if (use_pbwt) {
        if (type == 0) {
            if (hist_alts[1] >= 10) {
                model_2mc->pbwt.ReverseUpdateEWAH(ewah_data, ret_ewah, ret_buffer); 
                ret_len = n_samples;
            } else {
                // Unpack EWAH into literals according to current PPA
                uint32_t local_offset = ret_ewah_init;
                uint32_t ret_pos = 0;
                for (int j = 0; j < objs; ++j) {
                    djinn_ewah_t* ewah = (djinn_ewah_t*)&ewah_data[local_offset];
                    local_offset += sizeof(djinn_ewah_t);
                    // std::cerr << "ewah=" << ewah->ref << "," << ewah->clean << "," << ewah->dirty << std::endl;
                    
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

                    // local_offset += ewah->dirty * sizeof(uint32_t);
                    // vals += ewah->dirty + ewah->clean;

                    // std::cerr << "local=" << local_offset << "/" << ret_ewah << std::endl;
                    assert(ret_pos <= n_samples);
                    assert(local_offset <= ret_ewah);
                }
                assert(ret_pos == n_samples);
                ret_len = n_samples;
            }
        } else { // is NM
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
                // std::cerr << "ewah=" << ewah->ref << "," << ewah->clean << "," << ewah->dirty << std::endl;
                
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

                // local_offset += ewah->dirty * sizeof(uint32_t);
                // vals += ewah->dirty + ewah->clean;

                // std::cerr << "local=" << local_offset << "/" << ret_ewah << std::endl;
                assert(ret_pos <= n_samples);
                assert(local_offset <= ret_ewah);
            }
            assert(ret_pos == n_samples);
            ret_len = n_samples;
        } else {
            // Todo:
            std::cerr << "[NM PBWT] not avail" << std::endl;
            exit(1);
        }
    }

    return objs;
}

int djn_ctx_model_container_t::DecodeNextRaw(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;
    uint8_t type = marchetype->DecodeSymbol();
    
    switch(type) {
    case 0: return DecodeRaw(data, len);
    case 1: return DecodeRaw_nm(data, len);
    default: std::cerr << "[djn_ctx_model_container_t::DecodeNext] " << type << std::endl; return -1;
    }

    // Never reached.
    return -3;
}

// Return raw, potentially permuted, EWAH encoding
int djn_ctx_model_container_t::DecodeRaw(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;

    int64_t n_samples_obs = 0;
    int objects = 0;

    // Emit empty EWAH marker.
    djinn_ewah_t* ewah = (djinn_ewah_t*)&data[len]; 
    ewah->reset();
    len += sizeof(djinn_ewah_t);

    // Compute als
    memset(hist_alts, 0, 256*sizeof(uint32_t));

    while(true) {
        uint8_t type = model_2mc->mtype->DecodeSymbol();

        // std::cerr << "Type=" << (int)type << std::endl;
        if (type == 0) { // bitmaps
            ++ewah->dirty;

            // std::cerr << "Bitmap=";
            // std::cerr << "bitmap" << std::endl;
            for (int i = 0; i < 4; ++i) {
                data[len] = model_2mc->dirty_wah->DecodeSymbol();
                hist_alts[1] += __builtin_popcount(data[len]);
                ++len;
                // wah |= model_2mc->dirty_wah->DecodeSymbol();
                // std::cerr << std::bitset<8>(wah&255);
                // wah <<= 8;
            }
            // std::cerr << std::endl;
            n_samples_obs += 32;

        } else { // is RLE
            // Emit new EWAH marker
            if (ewah->clean > 0 || ewah->dirty > 0) {
                ewah = (djinn_ewah_t*)&data[len];
                ewah->reset();
                len += sizeof(djinn_ewah_t);
                ++objects;
            }

            // Decode an RLE
            uint32_t ref = 1; uint32_t len = 1;
            DecodeWahRLE(ref, len, model_2mc);
            // std::cerr << "decoded RLE=" << ref << " len=" << len << std::endl;
            ewah->ref = ref & 1; // todo: transform missing and EOV to symbols 14,15
            ewah->clean = len;

            // ref_alt[ewah->ref&1] += ewah->clean*32;
            hist_alts[ewah->ref] += ewah->clean * 32;

            n_samples_obs += len*32;
        }

        if (n_samples_obs == n_samples_wah) {
            // std::cerr << "obs=" << n_samples_obs << "/" << n_samples_wah << std::endl;
            break;
        }
        if (n_samples_obs > n_samples_wah) {
            std::cerr << "[djn_ctx_model_container_t::DecodeNextRaw] Decompression corruption: " << n_samples_obs << "/" << n_samples_wah  << std::endl;
            exit(1);
        }
    }

    if (ewah->clean > 0 || ewah->dirty > 0) {
        ++objects;
        // std::cerr << "final=" << std::bitset<64>(*ewah) << " ref=" << ewah_ref << ",run=" << ewah_run << ",dirty=" << ewah_dirty << std::endl;
        // ref_alt[ewah->ref&1] += ewah->clean*32;
        hist_alts[ewah->ref] += ewah->clean * 32;
        // ++len;
        // Todo: dirty
    }

    // std::cerr << "ewah decode ref_alt=" << ref_alt[0] << "," << ref_alt[1] << std::endl;

    return objects;
}

}