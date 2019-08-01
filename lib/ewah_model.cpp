#include <cstring> //memcpy

#include "djinn.h"
#include "pbwt.h"
#include "compressors.h"
#include "fastdelta.h"

namespace djinn {

/*======   EWAH container   ======*/

djn_ewah_model_t::djn_ewah_model_t() :
    pbwt(std::make_shared<PBWT>()),
    n_variants(0)
{
    
}

djn_ewah_model_t::~djn_ewah_model_t() {

}

void djn_ewah_model_t::reset() {
    if (pbwt.get() != nullptr) pbwt->Reset();
    n_variants = 0;
}

int djn_ewah_model_t::StartEncoding(bool use_pbwt, bool reset, bool store_offset) {
    if (reset)
        this->reset();
    
    data.p_len = 0;
    data.u_len = 0;
    offsets.p_len = 0;
    offsets.u_len = 0;
    internal_offsets.p_len = 0;
    internal_offsets.u_len = 0;

    return 1;
}

size_t djn_ewah_model_t::FinishEncoding(uint8_t*& support_buffer, uint32_t& support_cap, CompressionStrategy strat, int c_level) {
    GeneralCompressor comp;
    switch(strat) {
        case (CompressionStrategy::ZSTD): comp = &ZstdCompress; break;
        case (CompressionStrategy::LZ4):  comp = &Lz4Compress;  break;
    }
    data.u_len = data.p_len;
    if (data.p_len != 0) {
        int ret = (*comp)(data.p, data.p_len, support_buffer, support_cap, c_level);
        memcpy(data.p, support_buffer, ret); // copy data back to p
        data.p_len = ret;
    }

    if (offsets.p_len) {
        // std::cerr << "q_Len=" << q_len << std::endl;
        offsets.u_len = offsets.p_len;
        // delta
        compute_deltas_inplace(offsets.p, offsets.p_len, 0);
        int ret = (*comp)((uint8_t*)offsets.p, offsets.p_len*sizeof(uint32_t), support_buffer, support_cap, c_level);
        memcpy(offsets.p, support_buffer, ret); // copy data back to p
        offsets.p_len = ret;
        std::cerr << "[djn_ewah_model_t::FinishEncoding Q] debug=" << offsets.u_len << "->" << offsets.p_len << "->" << ret << std::endl;
    }

    if (internal_offsets.p_len) {
        // std::cerr << "q_Len=" << q_len << std::endl;
        internal_offsets.u_len = internal_offsets.p_len;
        // delta
        compute_deltas_inplace(internal_offsets.p, internal_offsets.p_len, 0);
        int ret = (*comp)((uint8_t*)internal_offsets.p, internal_offsets.p_len*sizeof(uint32_t), support_buffer, support_cap, c_level);
        memcpy(internal_offsets.p, support_buffer, ret); // copy data back to p
        internal_offsets.p_len = ret;
        std::cerr << "[djn_ewah_model_t::FinishEncoding OFF] debug=" << internal_offsets.u_len << "->" << internal_offsets.p_len << "->" << ret << std::endl;
    }

    // std::cerr << "[djn_ewah_model_t::FinishEncoding] debug=" << u_len << "->" << p_len << "->" << ret << std::endl;
    // std::cerr << "returning: " << p_len << " and " << q_len << "->" << p_len + q_len << std::endl;
    return data.p_len + offsets.p_len + internal_offsets.p_len;
}

int djn_ewah_model_t::StartDecoding(uint8_t*& support_buffer, uint32_t& support_cap, CompressionStrategy strat, bool use_pbwt, bool reset) {
    if (reset) this->reset();
    if (data.p == nullptr) return -2;

    GeneralDecompressor dec;
    switch(strat) {
        case (CompressionStrategy::ZSTD): dec = &ZstdDecompress; break;
        case (CompressionStrategy::LZ4):  dec = &Lz4Decompress;  break;
    }

    if (support_cap < data.u_len) {
        uint8_t* old = support_buffer;
        support_buffer = new uint8_t[data.u_len + 65536];
        memcpy(support_buffer, old, support_cap);
        delete[] old;
        support_cap = data.u_len + 65536;
    }

    if (data.p_len != 0) {
        int ret = (*dec)(data.p, data.p_len, support_buffer, support_cap);
        memcpy(data.p, support_buffer, ret); // copy data back to p
    }
    data.p_len = 0;

    if (offsets.p_len != 0) {
        int ret = (*dec)((uint8_t*)offsets.p, offsets.p_len, support_buffer, support_cap);
        memcpy((uint8_t*)offsets.p, support_buffer, ret); // copy data back to p
        offsets.p_len = offsets.u_len / sizeof(uint32_t);
        assert(ret == offsets.u_len * sizeof(uint32_t));
        compute_prefix_sum_inplace(offsets.p, offsets.u_len, 0);
        // for (int i = 0; i < offsets.p_len; ++i) {
        //     std::cerr << "," << offsets.p[i];
        // }
        // std::cerr << std::endl;
        offsets.p_len = 0;
    }

    if (internal_offsets.p_len != 0) {
        int ret = (*dec)((uint8_t*)internal_offsets.p, internal_offsets.p_len, support_buffer, support_cap);
        memcpy((uint8_t*)internal_offsets.p, support_buffer, ret); // copy data back to p
        internal_offsets.p_len = internal_offsets.u_len / sizeof(uint32_t);
        assert(ret == internal_offsets.u_len * sizeof(uint32_t));
        compute_prefix_sum_inplace(internal_offsets.p, internal_offsets.u_len, 0);
        for (int i = 0; i < internal_offsets.p_len; ++i) {
            std::cerr << "," << internal_offsets.p[i];
        }
        std::cerr << std::endl;
        internal_offsets.p_len = 0;
    }

    return 1;
}

int djn_ewah_model_t::Serialize(uint8_t* dst) const {
    // Serialize as (uint32_t,uint32_t,uint8_t*):
    // p_len, n_variants, p
    uint32_t offset = 0;
    *((uint32_t*)&dst[offset]) = data.p_len; // data length
    offset += sizeof(uint32_t);
    *((uint32_t*)&dst[offset]) = data.u_len; // data length
    offset += sizeof(uint32_t);
    *((uint32_t*)&dst[offset]) = n_variants; // number of variants
    offset += sizeof(uint32_t);
    // q data
    *((uint32_t*)&dst[offset]) = offsets.p_len;
    offset += sizeof(uint32_t);
    *((uint32_t*)&dst[offset]) = offsets.u_len;
    offset += sizeof(uint32_t);

    // q data
    *((uint32_t*)&dst[offset]) = internal_offsets.p_len;
    offset += sizeof(uint32_t);
    *((uint32_t*)&dst[offset]) = internal_offsets.u_len;
    offset += sizeof(uint32_t);
    
    memcpy(&dst[offset], data.p, data.p_len); // data
    offset += data.p_len;
    
    if (offsets.p_len) {
        memcpy(&dst[offset], offsets.p, offsets.p_len); // data
        offset += offsets.p_len;
    }

    if (internal_offsets.p_len) {
        memcpy(&dst[offset], internal_offsets.p, internal_offsets.p_len); // data
        offset += internal_offsets.p_len;
    }

    return offset;
}

int djn_ewah_model_t::Serialize(std::ostream& stream) const {
    // Serialize as (uint32_t,uint32_t,uint8_t*):
    // p_len, n_variants, p
    stream.write((char*)&data.p_len, sizeof(uint32_t));
    stream.write((char*)&data.u_len, sizeof(uint32_t));
    stream.write((char*)&n_variants, sizeof(uint32_t));
    stream.write((char*)&offsets.p_len, sizeof(uint32_t));
    stream.write((char*)&offsets.u_len, sizeof(uint32_t));
    stream.write((char*)&internal_offsets.p_len, sizeof(uint32_t));
    stream.write((char*)&internal_offsets.u_len, sizeof(uint32_t));   
    stream.write((char*)data.p, data.p_len);
    
    if (offsets.p_len) {
        stream.write((char*)offsets.p, offsets.p_len);
    }

    if (internal_offsets.p_len) {
        stream.write((char*)internal_offsets.p, internal_offsets.p_len);
    }
    return stream.tellp();
}

int djn_ewah_model_t::GetSerializedSize() const {
    int ret = 7*sizeof(uint32_t) + data.p_len;
    if (offsets.p_len) ret += offsets.p_len;
    if (internal_offsets.p_len) ret += internal_offsets.p_len;
    return ret;
}

int djn_ewah_model_t::Deserialize(uint8_t* dst) {
    uint32_t offset = 0;
    data.p_len = *((uint32_t*)&dst[offset]);
    offset += sizeof(uint32_t);
    data.u_len = *((uint32_t*)&dst[offset]);
    offset += sizeof(uint32_t);
    n_variants = *((uint32_t*)&dst[offset]);
    offset += sizeof(uint32_t);

    offsets.p_len = *((uint32_t*)&dst[offset]);
    offset += sizeof(uint32_t);
    offsets.u_len = *((uint32_t*)&dst[offset]);
    offset += sizeof(uint32_t);

    internal_offsets.p_len = *((uint32_t*)&dst[offset]);
    offset += sizeof(uint32_t);
    internal_offsets.u_len = *((uint32_t*)&dst[offset]);
    offset += sizeof(uint32_t);

    // initiate a buffer if there is none or it's too small
    if (data.p_cap == 0 || data.p == nullptr || data.u_len > data.p_cap) {
        // std::cerr << "[Deserialize] Limit. p_cap=" << p_cap << "," << "p is nullptr=" << (p == nullptr ? "yes" : "no") << ",p_len=" << p_len << "/" << p_cap << std::endl;
        if (data.p_free) delete[] data.p;
        data.p_cap = (data.u_len > data.p_len ? data.u_len : data.p_len) + 65536;
        data.p = new uint8_t[data.p_cap + 65536];
        data.p_free = true;
    }

    memcpy(data.p, &dst[offset], data.p_len); // data
    offset += data.p_len;

    if (offsets.p_len > offsets.p_cap || offsets.p_len && offsets.p == nullptr) {
        if (offsets.p_free) delete[] offsets.p;
        offsets.p_cap = (offsets.u_len > offsets.p_len ? offsets.u_len : offsets.p_len) + 65536;
        offsets.p = new uint32_t[offsets.p_cap + 65536];
        offsets.p_free = true;
    }

    if (offsets.p_len) {
        memcpy((char*)offsets.p, &dst[offset], offsets.p_len); // data
        offset += offsets.p_len;
    }

    if (internal_offsets.p_len > internal_offsets.p_cap || internal_offsets.p_len && internal_offsets.p == nullptr) {
        if (internal_offsets.p_free) delete[] internal_offsets.p;
        internal_offsets.p_cap = (internal_offsets.u_len > internal_offsets.p_len ? internal_offsets.u_len : internal_offsets.p_len) + 65536;
        internal_offsets.p = new uint32_t[internal_offsets.p_cap + 65536];
        internal_offsets.p_free = true;
    }

    if (internal_offsets.p_len) {
        memcpy((char*)internal_offsets.p, &dst[offset], internal_offsets.p_len); // data
        offset += internal_offsets.p_len;
    }
    
    return offset;
}

int djn_ewah_model_t::Deserialize(std::istream& stream) {
    // Serialize as (uint32_t,uint32_t,uint8_t*):
    // p_len, n_variants, p
    stream.read((char*)&data.p_len, sizeof(uint32_t));
    stream.read((char*)&data.u_len, sizeof(uint32_t));
    stream.read((char*)&n_variants, sizeof(uint32_t));
    stream.read((char*)&offsets.p_len, sizeof(uint32_t));
    stream.read((char*)&offsets.u_len, sizeof(uint32_t));
    stream.read((char*)&internal_offsets.p_len, sizeof(uint32_t));
    stream.read((char*)&internal_offsets.u_len, sizeof(uint32_t));

    // initiate a buffer if there is none or it's too small
    // std::cerr << "[Deserialize] " << p_len << "," << u_len << std::endl;
    if (data.p_cap == 0 || data.p == nullptr || data.u_len > data.p_cap) {
        // std::cerr << "[Deserialize] Limit. p_cap=" << p_cap << "," << "p is nullptr=" << (p == nullptr ? "yes" : "no") << ",p_len=" << p_len << "/" << p_cap << std::endl;
        if (data.p_free) delete[] data.p;
        data.p_cap = (data.u_len > data.p_len ? data.u_len : data.p_len) + 65536;
        data.p = new uint8_t[data.p_cap];
        data.p_free = true;
    }

    stream.read((char*)data.p, data.p_len);

    if (offsets.p_len > offsets.p_cap || offsets.p_len && offsets.p == nullptr) {
        if (offsets.p_free) delete[] offsets.p;
        offsets.p_cap = (offsets.u_len > offsets.p_len ? offsets.u_len : offsets.p_len) + 65536;
        offsets.p = new uint32_t[offsets.p_cap + 65536];
        offsets.p_free = true;
    }

    if (offsets.p_len) {
        stream.read((char*)offsets.p, offsets.p_len);
    }

    if (internal_offsets.p_len > internal_offsets.p_cap || internal_offsets.p_len && internal_offsets.p == nullptr) {
        if (internal_offsets.p_free) delete[] internal_offsets.p;
        internal_offsets.p_cap = (internal_offsets.u_len > internal_offsets.p_len ? internal_offsets.u_len : internal_offsets.p_len) + 65536;
        internal_offsets.p = new uint32_t[internal_offsets.p_cap + 65536];
        internal_offsets.p_free = true;
    }

    if (internal_offsets.p_len) {
        stream.read((char*)internal_offsets.p, internal_offsets.p_len);
    }

    return stream.tellg();
}

/*======   Variant EWAH model   ======*/

djinn_ewah_model::djinn_ewah_model() : 
    codec(CompressionStrategy::ZSTD), compression_level(DJINN_CLEVEL_DEFAULT),
    p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
    q(nullptr), q_len(0), q_alloc(0), q_free(true)
{
}

djinn_ewah_model::djinn_ewah_model(CompressionStrategy codec, int c_level) : 
    codec(codec), compression_level(c_level),
    p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
    q(nullptr), q_len(0), q_alloc(0), q_free(true)
{
    
    if (c_level <= 0) c_level = 1;
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
        ploidy_models.back()->store_offset = store_offset;
        tgt_container = ploidy_models[ploidy_models.size() - 1];
        tgt_container->StartEncoding(use_pbwt, init);
    }
    assert(tgt_container.get() != nullptr);

    // Compute allele counts.
    memset(hist_alts, 0, 256*sizeof(uint32_t));
    for (int i = 0; i < len_data; ++i) {
        ++hist_alts[DJN_BCF_UNPACK_GENOTYPE_GENERAL(data[i])];
    }

    // for (int i = 0; i < 256; ++i) {
    //     if (hist_alts[i]) std::cerr << i << ":" << hist_alts[i] << ",";
    // }
    // std::cerr << " data2m=" << tgt_container->model_2mc->range_coder->OutSize() << " dataNm=" << tgt_container->model_nm->range_coder->OutSize() << " permute=" << permute << std::endl;

    // Biallelic, no missing, and no special EOV symbols.
    if (alt_alleles <= 2 && hist_alts[14] == 0 && hist_alts[15] == 0) {
        tgt_container->p[tgt_container->p_len++] = 0;

        int ret = -1;
        if (use_pbwt) {
            // Todo: add lower limit to stored parameters during serialization
            if (hist_alts[1] < 10) { // dont update if < 10 alts
                for (int i = 0; i < len_data; ++i) {
                    tgt_container->model_2mc->pbwt->prev[i] = DJN_BCF_UNPACK_GENOTYPE(data[tgt_container->model_2mc->pbwt->ppa[i]]);
                }
            } else {
                tgt_container->model_2mc->pbwt->UpdateBcf(data, 1);
            }

            ret = (tgt_container->Encode2mc(tgt_container->model_2mc->pbwt->prev, len_data));
        } else {
            // std::cerr << "encoding nopbwt" << std::endl;
            ret = (tgt_container->Encode2mc(data, len_data, DJN_BCF_GT_UNPACK, 1));
        }

        if (ret > 0) ++n_variants;
        return ret;
    } else { // Otherwise.
        tgt_container->p[tgt_container->p_len++] = 1;
        
        int ret = -1;
        if (use_pbwt) {
            tgt_container->model_nm->pbwt->UpdateBcfGeneral(data, 1); // otherwise
            ret = (tgt_container->EncodeNm(tgt_container->model_nm->pbwt->prev, len_data));
        } else {
            ret = (tgt_container->EncodeNm(data, len_data, DJN_BCF_GT_UNPACK_GENERAL, 1));
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

    std::shared_ptr<djn_ewah_model_container_t> tgt_container;

    const uint64_t tuple = ((uint64_t)len_data << 32) | ploidy;
    auto search = ploidy_map.find(tuple);
    if (search != ploidy_map.end()) {
        tgt_container = ploidy_models[search->second];
        p[p_len++] = search->second;
    } else {
        // std::cerr << "Not found. Inserting: " << len_data << "," << ploidy << "(" << tuple << ") as " << ploidy_models.size() << std::endl;
        ploidy_map[tuple] = ploidy_models.size();
        p[p_len++] = ploidy_models.size();
        ploidy_models.push_back(std::make_shared<djn_ewah_model_container_t>(len_data, ploidy, (bool)use_pbwt));
        ploidy_models.back()->store_offset = store_offset;
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
        tgt_container->p[tgt_container->p_len++] = 0;// add archtype as 2mc

        int ret = -1;
        if (use_pbwt) {
            // Todo: add lower limit to stored parameters during serialization
            if (hist_alts[1] < 10) { // dont update if < 10 alts
                for (int i = 0; i < len_data; ++i) {
                    tgt_container->model_2mc->pbwt->prev[i] = data[tgt_container->model_2mc->pbwt->ppa[i]];
                }
            } else {
                tgt_container->model_2mc->pbwt->Update(data, 1);
            }

            ret = (tgt_container->Encode2mc(tgt_container->model_2mc->pbwt->prev, len_data));
        } else {
            ret = (tgt_container->Encode2mc(data, len_data, DJN_MAP_NONE, 0));
        }

        if (ret > 0) ++n_variants;
        return ret;
    } else { // Otherwise.
        tgt_container->p[tgt_container->p_len++] = 1;
        
        int ret = -1;
        if (use_pbwt) {
            tgt_container->model_nm->pbwt->Update(data, 1);
            ret = (tgt_container->EncodeNm(tgt_container->model_nm->pbwt->prev, len_data));
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
    variant->n_allele = 0;

    int ret = tgt_container->DecodeNext(q,q_len,variant->data,variant->data_len);
    if (ret <= 0) {
        variant->errcode = 1;
        return ret;
    }
    
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
    uint8_t type = p[p_len++];
    std::shared_ptr<djn_ewah_model_container_t> tgt_container = ploidy_models[type];
    return tgt_container->DecodeNextRaw(variant);
}

int djinn_ewah_model::DecodeNextRaw(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;
    uint8_t type = p[p_len++];

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
        ploidy_models[i]->store_offset = store_offset; // tell model to store offsets or not
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

    if (p_len != 0) {
        int ret = (*comp)(p, p_len, q, q_alloc, compression_level);
        memcpy(p, q, ret); // copy data back to p
        p_len = ret;
    }

    size_t s_models = 0;
    for (int i = 0; i < ploidy_models.size(); ++i) {
        uint32_t qa = q_alloc; // workaround for not being able to pass bit-field by reference
        int ret = ploidy_models[i]->FinishEncoding(q, qa, codec, compression_level); // use q as support buffer
        q_alloc = qa;
        s_models += ret;
    }
    
    return p_len + s_models;
}

int djinn_ewah_model::StartDecoding() {
    assert(p != nullptr);
    assert(ploidy_models.size() != 0);

    if (q == nullptr) {
        q_alloc = 10000000;
        q = new uint8_t[q_alloc];
        q_free = true;
        q_len  = 0;
    }
    q_len = 0;

    GeneralDecompressor dec;
    switch(codec) {
        case (CompressionStrategy::ZSTD): dec = &ZstdDecompress; break;
        case (CompressionStrategy::LZ4):  dec = &Lz4Decompress;  break;
    }
    if (p_len != 0) {
        int ret = (*dec)(p, p_len, q, q_alloc);
        memcpy(p, q, ret); // copy data back to p
    }
    p_len = 0;

    for (int i = 0; i <ploidy_models.size(); ++i) {
        uint32_t qa = q_alloc; // workaround for not being able to pass bit-field by reference
        ploidy_models[i]->StartDecoding(q,qa,codec,use_pbwt,init);
        q_alloc = qa;
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
    uint8_t pack = (use_pbwt << 7) | (init << 6) | (store_offset << 5) | (unused << 0);
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
    uint8_t pack = (use_pbwt << 7) | (init << 6) | (store_offset << 5) | (unused << 0);
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
    store_offset = (pack >> 5) & 1;
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
            ploidy_models.back()->store_offset = store_offset;
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
    uint32_t out_len = 0;
    stream.read((char*)&out_len, sizeof(uint32_t));

    // Codec
    int in_codec = 0;
    stream.read((char*)&in_codec, sizeof(int));
    codec = CompressionStrategy(in_codec);

    // n_models,n_variants
    int n_models = 0;
    stream.read((char*)&n_models,   sizeof(int));
    stream.read((char*)&n_variants, sizeof(uint32_t));

    // Deserialize bit-packed controller.
    uint8_t pack = 0;
    stream.read((char*)&pack, sizeof(uint8_t));
    use_pbwt = (pack >> 7) & 1;
    init = (pack >> 6) & 1;
    store_offset = (pack >> 5) & 1;
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
            ploidy_models.back()->store_offset = store_offset;
            ploidy_models.back()->Deserialize(stream);
        }
    }

    return stream.tellg();
}

/*======   Container   ======*/

djn_ewah_model_container_t::djn_ewah_model_container_t(int64_t n_s, int pl, bool use_pbwt) : 
    use_pbwt(use_pbwt),
    store_offset(false),
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
    store_offset(false),
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

    if (use_pbwt) {
        if (model_2mc->pbwt->n_symbols == 0) {
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_2mc->pbwt->Initiate(n_samples, 2);
        }

        if (model_nm->pbwt->n_symbols == 0) {
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_nm->pbwt->Initiate(n_samples, 16);
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

    model_2mc->StartEncoding(use_pbwt, reset, store_offset);
    model_nm->StartEncoding(use_pbwt, reset, store_offset);
}

size_t djn_ewah_model_container_t::FinishEncoding(uint8_t*& support_buffer, uint32_t& support_cap, CompressionStrategy strat, int c_level) {
    if (model_2mc.get() == nullptr) return -1;
    if (model_nm.get() == nullptr) return -1;

    int s_2mc = model_2mc->FinishEncoding(support_buffer, support_cap, strat, c_level);
    int s_nm  = model_nm->FinishEncoding(support_buffer, support_cap, strat, c_level);

    GeneralCompressor comp;
    switch(strat) {
        case (CompressionStrategy::ZSTD): comp = &ZstdCompress; break;
        case (CompressionStrategy::LZ4):  comp = &Lz4Compress; break;
    }

    if (p_len != 0) {
        int ret = (*comp)(p, p_len, support_buffer, support_cap, c_level);
        memcpy(p, support_buffer, ret); // copy data back to p
        p_len = ret;
    }

    size_t s_rc  = p_len;
    return s_rc + s_2mc + s_nm;
}

void djn_ewah_model_container_t::StartDecoding(uint8_t*& support_buffer, uint32_t& support_cap, CompressionStrategy strat, bool use_pbwt, bool reset) {
    if (model_2mc.get() == nullptr) return;
    if (model_nm.get() == nullptr) return;

    if (use_pbwt) {
        if (model_2mc->pbwt->n_symbols == 0) {
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_2mc->pbwt->Initiate(n_samples, 2);
        }

        if (model_nm->pbwt->n_symbols == 0) {
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_nm->pbwt->Initiate(n_samples, 16);
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
    if (p_len != 0) {    
        int ret = (*dec)(p, p_len, support_buffer, support_cap);
        memcpy(p, support_buffer, ret); // copy data back to p
    }
    p_len = 0;

    model_2mc->StartDecoding(support_buffer, support_cap, strat, use_pbwt, reset);
    model_nm->StartDecoding(support_buffer, support_cap, strat, use_pbwt, reset);
}

int djn_ewah_model_container_t::Encode2mc(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0)  return -2;
    
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

    for (int i = 0; i < n_samples; ++i) {
        if (map[data[i] >> shift]) {
            wah_bitmaps[i / 32] |= 1L << (i % 32);
        }
    }
    return EncodeWah(wah_bitmaps, n_samples_wah >> 5); // n_samples_wah / 32
}

int djn_ewah_model_container_t::EncodeNm(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0)  return -2;
    
    if (wah_bitmaps == nullptr) {
        n_wah = std::ceil((float)n_samples / 32) * 4;
        n_samples_wah = (n_wah * 32) / 4;
        n_samples_wah_nm = std::ceil((float)n_samples * 4 / 32) * 8;
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    for (int i = 0; i < n_samples; ++i) {
        wah_bitmaps[i / 8] |= (uint64_t)data[i] << (4*(i % 8));
    }

    return EncodeWahNm(wah_bitmaps, n_samples_wah_nm >> 3); // n_samples_wah_nm / 8
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
    
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    for (int i = 0; i < n_samples; ++i) {
        wah_bitmaps[i / 8] |= (uint64_t)map[data[i] >> shift] << (4*(i % 8));
    }

    return EncodeWahNm(wah_bitmaps, n_samples_wah_nm >> 3); // n_samples_wah_nm / 8
}

int djn_ewah_model_container_t::DecodeRaw_nm(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;
    if (model_nm.get() == nullptr) return -2;

    int64_t n_samples_obs = 0;
    int objects = 0;

    // Reset alts
    memset(hist_alts, 0, 256*sizeof(uint32_t));

    while(true) {
        // Emit empty EWAH marker.
        djinn_ewah_t* ewah = (djinn_ewah_t*)&data[len]; 
        ewah->reset();
        len += sizeof(djinn_ewah_t);
        
        djinn_ewah_t* e = (djinn_ewah_t*)&model_nm->data.p[model_nm->data.p_len]; 
        assert(e->clean + e->dirty > 0);
        ewah->ref   = e->ref;
        ewah->clean = e->clean;
        ewah->dirty = e->dirty;
        model_nm->data.p_len += sizeof(djinn_ewah_t);
        n_samples_obs += ewah->clean*8;
        n_samples_obs += ewah->dirty*8;
        hist_alts[ewah->ref & 15] += 8*ewah->clean;

        for (int i = 0; i < ewah->dirty; ++i) {
            uint32_t r = *((const uint32_t*)&model_nm->data.p[model_nm->data.p_len]); // copy
            *((uint32_t*)&data[len]) = r;
            for (int j = 0; j < 8; ++j) {
                ++hist_alts[r & 15];
                r >>= 4;
            }
            len += sizeof(uint32_t);
            model_nm->data.p_len += sizeof(uint32_t);
        }
        ++objects;

        if (n_samples_obs == n_samples_wah_nm) {
            break;
        }

        if (n_samples_obs > n_samples_wah_nm) {
            std::cerr << "[djn_ewah_model_container_t::DecodeRaw_nm] Decompression corruption: " << n_samples_obs << "/" << n_samples_wah  << std::endl;
            exit(1);
        }
    }

    uint32_t n_alts_obs = 0;
    for (int i = 0; i < 256; ++i) {
        n_alts_obs += hist_alts[i];
    }
    assert(n_alts_obs == n_samples_wah_nm);

    return objects;
}

int djn_ewah_model_container_t::EncodeWah(uint32_t* wah, uint32_t len) { // input WAH-encoded data
    if (wah == nullptr) return -1;
    if (model_2mc.get() == nullptr) return -1;

    // Resize if necessary.
    if (model_2mc->data.p_len + n_samples + 65536 > model_2mc->data.p_cap) {
        const uint32_t rc_size = model_2mc->data.p_len;
        // std::cerr << "[djn_ewah_model_container_t::EncodeWah][RESIZE] resizing from: " << rc_size << "->" << 2*rc_size << std::endl;
        uint8_t* prev = model_2mc->data.p; // old
        model_2mc->data.p_cap = model_2mc->data.p_len + 2*n_samples + 65536;
        model_2mc->data.p = new uint8_t[model_2mc->data.p_cap]; // double size. should rarely occur
        memcpy(model_2mc->data.p, prev, rc_size);
        if (model_2mc->data.p_free) delete[] prev;
        model_2mc->data.p_free = true;
    }

    if (store_offset && model_2mc->offsets.p_len + n_samples + 65536 > model_2mc->offsets.p_cap) {
        // std::cerr << "resizing qlen" << std::endl;
        const uint32_t rc_size = model_2mc->offsets.p_len;
        // std::cerr << "[djn_ewah_model_container_t::EncodeWah][RESIZE] resizing from: " << rc_size << "->" << 2*rc_size << std::endl;
        uint32_t* prev = model_2mc->offsets.p; // old
        model_2mc->offsets.p_cap = model_2mc->offsets.p_len + 2*n_samples + 65536;
        model_2mc->offsets.p = new uint32_t[model_2mc->offsets.p_cap]; // double size. should rarely occur
        memcpy(model_2mc->offsets.p, prev, rc_size*sizeof(uint32_t));
        if (model_2mc->offsets.p_free) delete[] prev;
        model_2mc->offsets.p_free = true;
    }

    if (store_offset && model_2mc->internal_offsets.p_len + n_samples + 65536 > model_2mc->internal_offsets.p_cap) {
        // std::cerr << "resizing qlen" << std::endl;
        const uint32_t rc_size = model_2mc->internal_offsets.p_len;
        // std::cerr << "[djn_ewah_model_container_t::EncodeWah][RESIZE] resizing from: " << rc_size << "->" << 2*rc_size << std::endl;
        uint32_t* prev = model_2mc->internal_offsets.p; // old
        model_2mc->internal_offsets.p_cap = model_2mc->internal_offsets.p_len + 2*n_samples + 65536;
        model_2mc->internal_offsets.p = new uint32_t[model_2mc->internal_offsets.p_cap]; // double size. should rarely occur
        memcpy(model_2mc->internal_offsets.p, prev, rc_size*sizeof(uint32_t));
        if (model_2mc->internal_offsets.p_free) delete[] prev;
        model_2mc->internal_offsets.p_free = true;
    }

    uint32_t n_objs = 1;
    uint32_t n_obs  = 0;
    
    djinn_ewah_t* ewah = (djinn_ewah_t*)&model_2mc->data.p[model_2mc->data.p_len];
    ewah->reset();
    if (store_offset) {
        model_2mc->offsets.p[model_2mc->offsets.p_len++] = model_2mc->data.p_len; // each EWAH position
        model_2mc->internal_offsets.p[model_2mc->internal_offsets.p_len++] = model_2mc->data.p_len; //
        // std::cerr << model_2mc->data.p_len << std::endl;
    }
    model_2mc->data.p_len += sizeof(djinn_ewah_t);

    for (int i = 0; i < len; ++i) {
        // Is dirty
        if (wah[i] != 0 && wah[i] != std::numeric_limits<uint32_t>::max()) {
            ++ewah->dirty;
            *((uint32_t*)&model_2mc->data.p[model_2mc->data.p_len]) = wah[i];
            model_2mc->data.p_len += sizeof(uint32_t);
            ++n_obs;
        } 
        // Is clean
        else {
            // Dirty have been set then make new EWAH
            if (ewah->dirty) {
                // Only make a new EWAH if anything is set
                n_obs += ewah->clean;
                assert(ewah->clean + ewah->dirty > 0);
                ewah = (djinn_ewah_t*)&model_2mc->data.p[model_2mc->data.p_len];
                ewah->reset();
                if (store_offset) {
                    model_2mc->offsets.p[model_2mc->offsets.p_len++] = model_2mc->data.p_len;
                    // model_2mc->internal_offsets.p[model_2mc->internal_offsets.p_len++] = model_2mc->data.p_len;
                    // std::cerr << model_2mc->data.p_len << std::endl;
                }
                model_2mc->data.p_len += sizeof(djinn_ewah_t);
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
                        ewah = (djinn_ewah_t*)&model_2mc->data.p[model_2mc->data.p_len];
                        ewah->reset();
                        if (store_offset) {
                            model_2mc->offsets.p[model_2mc->offsets.p_len++] = model_2mc->data.p_len;
                            // model_2mc->internal_offsets.p[model_2mc->internal_offsets.p_len++] = model_2mc->data.p_len;
                            // std::cerr << model_2mc->data.p_len << std::endl;
                        }
                        model_2mc->data.p_len += sizeof(djinn_ewah_t);
                    }
                    ++ewah->clean;
                    ewah->ref = (wah[i] & 1);
                    ++n_objs;
                }
            }
        }
    }
    n_obs += ewah->clean;
    assert(n_obs == len);

    ++model_2mc->n_variants;
    ++n_variants;
    return n_objs;
}

int djn_ewah_model_container_t::EncodeWahNm(uint32_t* wah, uint32_t len) { // input WAH-encoded data
    if (wah == nullptr) return -1;
    if (model_nm.get() == nullptr) return -1;

    // Resize if necessary.
    if (model_nm->data.p_len + n_samples + 65536 > model_nm->data.p_cap) {
        const uint32_t rc_size = model_nm->data.p_len;
        uint8_t* prev = model_nm->data.p; // old
        model_nm->data.p_cap = model_nm->data.p_len + 2*n_samples + 65536;
        model_nm->data.p = new uint8_t[model_nm->data.p_cap]; // double size. should rarely occur
        memcpy(model_nm->data.p, prev, rc_size);
        if (model_nm->data.p_free) delete[] prev;
        model_nm->data.p_free = true;
    }

    if (store_offset && model_nm->offsets.p_len + n_samples + 65536 > model_nm->offsets.p_cap) {
        const uint32_t rc_size = model_nm->offsets.p_len;
        // std::cerr << "[djn_ewah_model_container_t::EncodeWah][RESIZE] resizing from: " << rc_size << "->" << 2*rc_size << std::endl;
        uint32_t* prev = model_nm->offsets.p; // old
        model_nm->offsets.p_cap = model_nm->offsets.p_len + 2*n_samples + 65536;
        model_nm->offsets.p = new uint32_t[model_nm->offsets.p_cap]; // double size. should rarely occur
        memcpy(model_nm->offsets.p, prev, rc_size*sizeof(uint32_t));
        if (model_nm->offsets.p_free) delete[] prev;
        model_nm->offsets.p_free = true;
    }

    if (store_offset && model_nm->internal_offsets.p_len + n_samples + 65536 > model_nm->internal_offsets.p_cap) {
        // std::cerr << "resizing qlen" << std::endl;
        const uint32_t rc_size = model_nm->internal_offsets.p_len;
        // std::cerr << "[djn_ewah_model_container_t::EncodeWah][RESIZE] resizing from: " << rc_size << "->" << 2*rc_size << std::endl;
        uint32_t* prev = model_nm->internal_offsets.p; // old
        model_nm->internal_offsets.p_cap = model_nm->internal_offsets.p_len + 2*n_samples + 65536;
        model_nm->internal_offsets.p = new uint32_t[model_nm->internal_offsets.p_cap]; // double size. should rarely occur
        memcpy(model_nm->internal_offsets.p, prev, rc_size*sizeof(uint32_t));
        if (model_nm->internal_offsets.p_free) delete[] prev;
        model_nm->internal_offsets.p_free = true;
    }

    // Debug
    uint32_t n_objs = 1;
    uint32_t n_obs  = 0;
    
    djinn_ewah_t* ewah = (djinn_ewah_t*)&model_nm->data.p[model_nm->data.p_len];
    ewah->reset();
    if (store_offset) {
        model_nm->offsets.p[model_nm->offsets.p_len++] = model_nm->data.p_len;
        model_nm->internal_offsets.p[model_nm->internal_offsets.p_len++] = model_nm->data.p_len;
    }
    model_nm->data.p_len += sizeof(djinn_ewah_t);

    for (int i = 0; i < len; ++i) {
        // Is dirty
        if (wah[i] != DJN_NM_REF_BITS[wah[i] & 15]) {
            ++ewah->dirty;
            *((uint32_t*)&model_nm->data.p[model_nm->data.p_len]) = wah[i];
            model_nm->data.p_len += sizeof(uint32_t);
            ++n_obs;
        } 
        // Is clean
        else {
            // Dirty have been set then make new EWAH
            if (ewah->dirty) {
                n_obs += ewah->clean;
                assert(ewah->clean + ewah->dirty > 0);
                ewah = (djinn_ewah_t*)&model_nm->data.p[model_nm->data.p_len];
                ewah->reset();
                if (store_offset) {
                    model_nm->offsets.p[model_nm->offsets.p_len++] = model_nm->data.p_len;
                    model_nm->internal_offsets.p[model_nm->internal_offsets.p_len++] = model_nm->data.p_len;
                }
                model_nm->data.p_len += sizeof(djinn_ewah_t);
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
                        ewah = (djinn_ewah_t*)&model_nm->data.p[model_nm->data.p_len];
                        ewah->reset();
                        if (store_offset) {
                            model_nm->offsets.p[model_nm->offsets.p_len++] = model_nm->data.p_len;
                            model_nm->internal_offsets.p[model_nm->internal_offsets.p_len++] = model_nm->data.p_len;
                        }
                        model_nm->data.p_len += sizeof(djinn_ewah_t);
                    }
                    ++ewah->clean;
                    ewah->ref = (wah[i] & 15);
                    ++n_objs;
                }
            }
        }
    }
    n_obs += ewah->clean;
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
    uint8_t type = p[p_len++];

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
                 model_2mc->pbwt->ReverseUpdateEWAH(ewah_data, ret_ewah, ret_buffer); 
                ret_len = n_samples;
            } else {
                // Unpack EWAH into literals according to current PPA
                uint32_t local_offset = ret_ewah_init;
                uint32_t ret_pos = 0;
                for (int j = 0; j < objs; ++j) {
                    djinn_ewah_t* ewah = (djinn_ewah_t*)&ewah_data[local_offset];
                    local_offset += sizeof(djinn_ewah_t);

                    // Clean words.
                    uint32_t to = ret_pos + ewah->clean*32 > n_samples ? n_samples : ret_pos + ewah->clean*32;
                    for (int i = ret_pos; i < to; ++i) {
                        ret_buffer[model_2mc->pbwt->ppa[i]] = (ewah->ref & 1);
                    }
                    ret_pos = to;

                    for (int i = 0; i < ewah->dirty; ++i) {
                        to = ret_pos + 32 > n_samples ? n_samples : ret_pos + 32;
                        
                        uint32_t dirty = *((uint32_t*)(&ewah_data[local_offset])); // copy
                        for (int j = ret_pos; j < to; ++j) {
                            ret_buffer[model_2mc->pbwt->ppa[j]] = (dirty & 1);
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
            model_nm->pbwt->ReverseUpdateEWAHNm(ewah_data, ret_ewah, ret_buffer);
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
    uint8_t type = p[p_len++];

    switch(type) {
    case 0: return(DecodeRaw(data, len)); break;
    case 1: return(DecodeRaw_nm(data, len)); break;
    default: std::cerr << "[djn_ewah_model_container_t::DecodeNextRaw] decoding error: " << (int)type << " (valid=[0,1])" << std::endl; return -1;
    }

    // Never reached.
    return -3;
}

int djn_ewah_model_container_t::DecodeRawRandomAccess(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;
    if (model_2mc.get() == nullptr) return -1;
    if (model_nm.get() == nullptr) return -1;
    
    uint32_t start_offset = model_2mc->internal_offsets.p[model_2mc->internal_offsets.p_len++];
    // uint32_t start_offset = model_2mc->internal_offsets.p[model_2mc->internal_offsets.p_len++];


    return -1;
}

int djn_ewah_model_container_t::DecodeRawRandomAccess(djinn_variant_t*& variant) {
    return -1;
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
        
        djinn_ewah_t* e = (djinn_ewah_t*)&model_2mc->data.p[model_2mc->data.p_len]; 
        // std::cerr << "EWAH: ref=" << e->ref << ",clean=" << e->clean << ",dirty=" << e->dirty << " offset=" << model_2mc->data.p_len << "/" << model_2mc->data.p_cap << std::endl;
        assert(e->clean + e->dirty > 0);
        ewah->ref   = e->ref;
        ewah->clean = e->clean;
        ewah->dirty = e->dirty;
        model_2mc->data.p_len += sizeof(djinn_ewah_t);
        n_samples_obs += ewah->clean*32;
        n_samples_obs += ewah->dirty*32;
        hist_alts[ewah->ref & 1] += 32*ewah->clean;

        for (int i = 0; i < ewah->dirty; ++i) {
            const uint32_t* r = (const uint32_t*)&model_2mc->data.p[model_2mc->data.p_len];
            *((uint32_t*)&data[len]) = *r;
            hist_alts[1] += __builtin_popcount(*r);
            hist_alts[0] += __builtin_popcount(~(*r));
            len += sizeof(uint32_t);
            model_2mc->data.p_len += sizeof(uint32_t);
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

    // std::cerr << "[decoded 2nm] " << hist_alts[0] << "," <<hist_alts[1] << std::endl;

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
    variant->n_allele = 0;

    // Decode stream archetype.
    uint8_t type = p[p_len++];

    int ret = 0;
    switch(type) {
    case 0: ret = DecodeRaw(variant->data, variant->data_len);    break;
    case 1: ret = DecodeRaw_nm(variant->data, variant->data_len); break;
    default: std::cerr << "[djn_ctx_model_container_t::DecodeNextRaw] decoding error: " << (int)type << " (valid=[0,1])" << std::endl; return -1;
    }

    if (ret <= 0) {
        variant->errcode = 1;
        return ret;
    }

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
    variant->d->n_samples  = n_samples;

    // Construct EWAH mapping
    uint32_t local_offset = 0;
    uint32_t ret_pos = 0;

    // std::cerr << "here" << std::endl;
    for (int j = 0; j < ret; ++j) {
        djinn_ewah_t* ewah = (djinn_ewah_t*)&variant->data[local_offset];
        variant->d->ewah[variant->d->n_ewah++] = (djinn_ewah_t*)&variant->data[local_offset];
        local_offset += sizeof(djinn_ewah_t);
    
        // Clean words.
        uint32_t to = ret_pos + ewah->clean*32 > n_samples ? n_samples : ret_pos + ewah->clean*32;
        ret_pos = to;

        // Store first dirty pointer only.
        variant->d->dirty[variant->d->n_dirty++] = (uint32_t*)(&variant->data[local_offset]);

        for (int i = 0; i < ewah->dirty; ++i) {
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
    }
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
    int ret = p_len;
    ret += model_2mc->data.p_len;
    ret += model_nm->data.p_len;
    return ret;
}

int djn_ewah_model_container_t::Deserialize(uint8_t* dst) {
    uint32_t offset = 0;
    int pl = *((int*)&dst[offset]);
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
    offset += model_2mc->Deserialize(&dst[offset]);
    offset += model_nm->Deserialize(&dst[offset]);

    return(offset);
}

int djn_ewah_model_container_t::Deserialize(std::istream& stream) {
    // #pl and #n_s read outside of this function in Deserialize() for
    // the parent.
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