#include <cstring> //memcpy
#include "ewah_model.h"
#include "compressors.h"

namespace djinn {

/*======   Context container   ======*/

djn_ewah_model_t::djn_ewah_model_t() :
    p(nullptr),
    p_len(0), p_cap(0), p_free(false),
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

    return 1;
}

size_t djn_ewah_model_t::FinishEncoding() {
    return p_len;
}

int djn_ewah_model_t::StartDecoding(bool use_pbwt, bool reset) {
    if (reset) this->reset();
    if (p == nullptr) return -2; // or result in corruption as range coder immediately loads data

    return 1;
}

/*======   Varaint context model   ======*/

djinn_ewah_model::djinn_ewah_model() : 
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
        // tgt_container = ploidy_models[search->second];
        ploidy_dict->EncodeSymbol(search->second);
        p[p_len++] = search->second;
    } else {
        // std::cerr << "Not found. Inserting: " << len_data << "," << ploidy << "(" << tuple << ") as " << ploidy_models.size() << std::endl;
        ploidy_map[tuple] = ploidy_models.size();
        // ploidy_dict->EncodeSymbol(ploidy_models.size());
        p[p_len++] = ploidy_models.size();
        ploidy_models.push_back(std::make_shared<djn_ewah_model_container_t>(len_data, ploidy, (bool)use_pbwt));
        tgt_container = ploidy_models[ploidy_models.size() - 1];
        tgt_container->StartEncoding(use_pbwt, init); // Todo: fix me
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
    uint8_t type = ploidy_dict->DecodeSymbol();
    // std::cerr << "[Decode model] Stream=" << (int)type << std::endl;
    std::shared_ptr<djn_ewah_model_container_t> tgt_container = ploidy_models[type];

    return(tgt_container->DecodeNext(ewah_data,ret_ewah,ret_buffer,ret_len));
}

int djinn_ewah_model::DecodeNext(djinn_variant_t*& variant) {
    // Decode stream archetype.
    uint8_t type = ploidy_dict->DecodeSymbol();
    // std::cerr << "[Decode model] Stream=" << (int)type << std::endl;
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

    variant->ploidy = tgt_container->ploidy;
    variant->data_len = 0;
    variant->errcode = 0;

    return(tgt_container->DecodeNext(q,q_len,variant->data,variant->data_len));
}

int djinn_ewah_model::DecodeNextRaw(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;
    uint8_t type = ploidy_dict->DecodeSymbol();

    std::shared_ptr<djn_ewah_model_container_t> tgt_container = ploidy_models[type];
    return(tgt_container->DecodeNextRaw(data, len));
}

void djinn_ewah_model::StartEncoding(bool use_pbwt, bool reset) {
    // Store parameters for decoding.
    this->use_pbwt = use_pbwt;
    this->init = reset;
    n_variants = 0;

    std::cerr << "[djinn_ewah_model::StartEncoding] models start encoding" << std::endl;
    for (int i = 0; i < ploidy_models.size(); ++i) {
        ploidy_models[i]->StartEncoding(use_pbwt, reset);
    }
}

size_t djinn_ewah_model::FinishEncoding() {
    std::cerr << "[djinn_ewah_model::FinishEncoding] p_len=" << p_len << std::endl;

    size_t s_models = 0;
    for (int i = 0; i < ploidy_models.size(); ++i) {
        s_models += ploidy_models[i]->FinishEncoding();
    }

    std::cerr << "Model sizes=" << p_len << " and models=" << s_models << std::endl;
    
    return p_len + s_models;
}

int djinn_ewah_model::StartDecoding() {
    assert(p != nullptr);
    assert(ploidy_models.size() != 0);

    for (int i = 0; i <ploidy_models.size(); ++i) {
        ploidy_models[i]->StartDecoding(use_pbwt, init);
    }

    return 1;
}

/*======   Container   ======*/

djn_ewah_model_container_t::djn_ewah_model_container_t(int64_t n_s, int pl, bool use_pbwt) : 
    use_pbwt(use_pbwt),
    ploidy(pl), n_samples(n_s), n_variants(0),
    n_wah(std::ceil((float)n_samples / 32) * 4), 
    n_samples_wah((n_wah * 32) / 4), 
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
    n_wah(std::ceil((float)n_samples / 32) * 4), 
    n_samples_wah((n_wah * 32) / 4), 
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
    std::cerr << "in start encoding: " << model_2mc->pbwt.n_symbols << std::endl; 
    if (use_pbwt) {
        if (model_2mc->pbwt.n_symbols == 0) {
            std::cerr << "init pbwt: " << n_samples << "," << 2 << std::endl;
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_2mc->pbwt.Initiate(n_samples, 2);
        }

        if (model_nm->pbwt.n_symbols == 0) {
            std::cerr << "init pbwt: " << n_samples << "," << 16 << std::endl;
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_nm->pbwt.Initiate(n_samples, 16);
        }
    }

    if (reset) {
        std::cerr << "[djn_ewah_model_container_t::StartEncoding] resetting" << std::endl;
        model_2mc->reset();
        model_nm->reset();
    } else {
        model_2mc->n_variants = 0;
        model_nm->n_variants = 0;
    }

    this->use_pbwt = use_pbwt;
    n_variants = 0;

    model_2mc->StartEncoding(use_pbwt, reset);
    model_nm->StartEncoding(use_pbwt, reset);
}

size_t djn_ewah_model_container_t::FinishEncoding() {
    model_2mc->FinishEncoding();
    model_nm->FinishEncoding();

    size_t s_rc  = p_len;
    size_t s_2mc = model_2mc->p_len;
    size_t s_nm  = model_nm->p_len;
    std::cerr << "container finish=" << s_rc << " and " << s_2mc << " and " << s_nm << std::endl; 
    
    return s_rc + s_2mc + s_nm;
}

void djn_ewah_model_container_t::StartDecoding(bool use_pbwt, bool reset) {
    // std::cerr << "in start decoding: " << model_2mc->pbwt.n_symbols << std::endl; 
    if (use_pbwt) {
        if (model_2mc->pbwt.n_symbols == 0) {
            std::cerr << "init pbwt: " << n_samples << "," << 2 << std::endl;
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_2mc->pbwt.Initiate(n_samples, 2);
        }

        if (model_nm->pbwt.n_symbols == 0) {
            std::cerr << "init pbwt: " << n_samples << "," << 16 << std::endl;
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_nm->pbwt.Initiate(n_samples, 16);
        }
    }

    if (reset) {
        std::cerr << "[djn_ewah_model_container_t::StartDecoding] resetting" << std::endl;
        model_2mc->reset();
        model_nm->reset();
    } else {
        model_2mc->n_variants = 0;
        model_nm->n_variants = 0;
    }
    std::cerr << "[djn_ewah_model_container_t::StartDecoding] p_len=" << p_len << std::endl;

    this->use_pbwt = use_pbwt;

    model_2mc->StartDecoding(use_pbwt, reset);
    model_nm->StartDecoding(use_pbwt, reset);
}

int djn_ewah_model_container_t::Encode2mc(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;
    
    if (wah_bitmaps == nullptr) {
        n_wah = std::ceil((float)n_samples / 32) * 4;
        n_samples_wah = (n_wah * 32) / 4;
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    for (int i = 0; i < n_samples; ++i) {
        if (data[i]) {
            wah_bitmaps[i / 32] |= 1L << (i % 32);
        }
    }

    return EncodeWah(wah_bitmaps, n_wah/4); // 4 bytes in a 32-bit word
}

int djn_ewah_model_container_t::Encode2mc(uint8_t* data, uint32_t len, const uint8_t* map, const int shift) {
    if (data == nullptr) return -1;
    if (n_samples == 0)  return -2;
    if (map == nullptr)  return -3;
    
    if (wah_bitmaps == nullptr) {
        n_wah = std::ceil((float)n_samples / 32) * 4;
        n_samples_wah = (n_wah * 32) / 4;
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

    return EncodeWah(wah_bitmaps, n_wah/4); // 4 bytes in a 32-bit word
}

int djn_ewah_model_container_t::EncodeNm(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;
    
    if (wah_bitmaps == nullptr) {
        n_wah = std::ceil((float)n_samples / 32) * 4;
        n_samples_wah = (n_wah * 32) / 4;
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    // Todo: move into ResetBitmaps() functions
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    for (int i = 0; i < n_samples; ++i) {
        wah_bitmaps[i / 8] |= (uint64_t)data[i] << (4*(i % 8));
    }

    // for (int i = 0; i < n_wah; ++i) std::cerr << std::bitset<64>(wah_bitmaps[i]) << " ";
    // std::cerr << std::endl;

    return EncodeWahNm(wah_bitmaps, n_wah);
    // return -1;
}

int djn_ewah_model_container_t::EncodeNm(uint8_t* data, uint32_t len, const uint8_t* map, const int shift) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;
    if (map == nullptr) return -3;

    if (wah_bitmaps == nullptr) {
        n_wah = std::ceil((float)n_samples / 32) * 4;
        n_samples_wah = (n_wah * 32) / 4;
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    // Todo: move into ResetBitmaps() functions
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    for (int i = 0; i < n_samples; ++i) {
        wah_bitmaps[i / 8] |= (uint64_t)map[data[i] >> shift] << (4*(i % 8));
    }

    // for (int i = 0; i < n_wah; ++i) std::cerr << std::bitset<64>(wah_bitmaps[i]) << " ";
    // std::cerr << std::endl;

    return EncodeWahNm(wah_bitmaps, n_wah);
    // return -1;
}

int djn_ewah_model_container_t::DecodeRaw_nm(uint8_t* data, uint32_t& len) {
    return 1;
}

int djn_ewah_model_container_t::EncodeWah(uint32_t* wah, uint32_t len) { // input WAH-encoded data
    if (wah == nullptr) return -1;
    ++model_2mc->n_variants;

    // Resize if necessary.
    if (model_2mc->p_len + 2*n_samples > model_2mc->p_cap) {
        const uint32_t rc_size = model_2mc->p_len;
        std::cerr << "[djn_ewah_model_container_t::EncodeWah][RESIZE] resizing from: " << rc_size << "->" << 2*rc_size << std::endl;
        uint8_t* prev = model_2mc->p; // old
        model_2mc->p = new uint8_t[rc_size*2]; // double size. should rarely occur
        model_2mc->p_cap = rc_size*2;
        memcpy(model_2mc->p, prev, rc_size);
        if (model_2mc->p_free) delete[] prev;
        model_2mc->p_free = true;
        std::cerr << "done resizing" << std::endl;
    }

    djinn_ewah_t* ewah = (djinn_ewah_t*)&model_2mc->p[model_2mc->p_len];
    ewah->reset();
    model_2mc->p_len += sizeof(djinn_ewah_t);
    
    // bool prev_dirty = (wah[0] != 0 && wah[0] != std::numeric_limits<uint32_t>::max());
    ewah->dirty +=  (wah[0] != 0 && wah[0] != std::numeric_limits<uint32_t>::max());
    if (ewah->dirty) {
        *((uint32_t*)&model_2mc->p[model_2mc->p_len]) = wah[0];
        model_2mc->p_len += sizeof(uint32_t);
    }

    ewah->clean += !(wah[0] != 0 && wah[0] != std::numeric_limits<uint32_t>::max());
    if (ewah->clean) {
        ewah->ref = (wah[0] & 1);
    }
    assert(ewah->dirty + ewah->clean == 1);

    for (int i = 1; i < len; ++i) {
        // If WAH != 0 AND WAH != 111..1 then dirty
        // If previous word was dirty (count > 0) keep extending
        // If previous word was clean then add (count == 0)
        
        // Is dirty
        if ((wah[i] != 0 && wah[i] != std::numeric_limits<uint32_t>::max())) {
            ++ewah->dirty;
            *((uint32_t*)&model_2mc->p[model_2mc->p_len]) = wah[i];
            model_2mc->p_len += sizeof(uint32_t);
        } 
        // Is clean
        else {
            // Dirty have been set then make new EWAH
            if (ewah->dirty) {
                ewah = (djinn_ewah_t*)&model_2mc->p[model_2mc->p_len];
                ewah->reset();
                model_2mc->p_len += sizeof(djinn_ewah_t);
                ++ewah->clean;
                ewah->ref = (wah[i] & 1);
            } 
            // No dirty have been set
            else {
                if (ewah->ref == (wah[i] & 1)) ++ewah->clean;
                else { // Different clean word
                    ewah = (djinn_ewah_t*)&model_2mc->p[model_2mc->p_len];
                    ewah->reset();
                    model_2mc->p_len += sizeof(djinn_ewah_t);
                    ++ewah->clean;
                    ewah->ref = (wah[i] & 1);
                }
            }
        }
    }

    ++n_variants;
    return 1;
}

int djn_ewah_model_container_t::EncodeWahNm(uint32_t* wah, uint32_t len) { // input WAH-encoded data
    if (wah == nullptr) return -1;
    ++model_nm->n_variants;

    ++n_variants;
    return 1;
}

int djn_ewah_model_container_t::DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len) {
    if (ewah_data == nullptr) return -1;
    if (ret_buffer == nullptr) return -1;

    return 0;
}

int djn_ewah_model_container_t::DecodeNextRaw(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;
    
    return 1;
}

// Return raw, potentially permuted, EWAH encoding
int djn_ewah_model_container_t::DecodeRaw(uint8_t* data, uint32_t& len) {
    if (data == nullptr) return -1;

    return 1;
}

}