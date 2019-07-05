/*
* Copyright (c) 2019 Marcus D. R. Klarqvist
* Author(s): Marcus D. R. Klarqvist
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*   http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing,
* software distributed under the License is distributed on an
* "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
* KIND, either express or implied.  See the License for the
* specific language governing permissions and limitations
* under the License.
*/
#ifndef DJINN_CTX_MODEL_H_
#define DJINN_CTX_MODEL_H_

#include <unordered_map> //std::unordered_map
#include <limits> //std::numeric_limit

#include "djinn.h" // core components
#include "frequency_model.h" // RangeCoder and FrequencyModel
#include "pbwt.h" // PBWT algorithms

namespace djinn {

struct djn_ctx_model_t {
public:
    djn_ctx_model_t();
    ~djn_ctx_model_t();

    // Initiate this model as accepting either biallelic complete
    // input data (Initiate2mc) or any other (InitiateNm). Running
    // either of these functions are REQUIRED before starting encoding
    // or decoding.
    void Initiate2mc();
    void InitiateNm();

    int StartEncoding(bool use_pbwt, bool reset = false);
    size_t FinishEncoding();
    int StartDecoding(bool use_pbwt, bool reset = false);
    size_t FinishDecoding() { return 0; } // no effect
    
    void reset();

    // Read/write
    int Serialize(uint8_t* dst) const {
        // Serialize as (uint32_t,uint32_t,uint8_t*):
        // p_len, n_variants, p
        uint32_t offset = 0;
        *((uint32_t*)&dst[offset]) = p_len; // data length
        offset += sizeof(uint32_t);
        *((uint32_t*)&dst[offset]) = n_variants; // number of variants
        offset += sizeof(uint32_t);
        memcpy(&dst[offset], p, p_len); // data
        offset += p_len;

        // assert(range_coder->OutSize() == p_len);
        return offset;
    }

    int Serialize(std::ostream& stream) const {
        // Serialize as (uint32_t,uint32_t,uint8_t*):
        // p_len, n_variants, p
        stream.write((char*)&p_len, sizeof(uint32_t));
        stream.write((char*)&n_variants, sizeof(uint32_t));
        stream.write((char*)p, p_len);
        return stream.tellp();
    }
    
    int GetSerializedSize() const {
        int ret = sizeof(uint32_t) + sizeof(uint32_t) + p_len;
        return ret;
    }

    int GetCurrentSize() const {
        if (range_coder.get() == nullptr) return -1;
        return range_coder->OutSize();
    }
    
    int Deserialize(uint8_t* dst) {
        uint32_t offset = 0;
        p_len = *((uint32_t*)&dst[offset]);
        offset += sizeof(uint32_t);
        n_variants = *((uint32_t*)&dst[offset]);
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
        
        return offset;
    }

    int Deserialize(std::istream& stream) {
        // Serialize as (uint32_t,uint32_t,uint8_t*):
        // p_len, n_variants, p
        stream.read((char*)&p_len, sizeof(uint32_t));
        stream.read((char*)&n_variants, sizeof(uint32_t));

        // initiate a buffer if there is none or it's too small
        if (p_cap == 0 || p == nullptr || p_len > p_cap) {
            // std::cerr << "[Deserialize] Limit. p_cap=" << p_cap << "," << "p is nullptr=" << (p == nullptr ? "yes" : "no") << ",p_len=" << p_len << "/" << p_cap << std::endl;
            if (p_free) delete[] p;
            p = new uint8_t[p_len + 65536];
            p_cap = p_len + 65536;
            p_free = true;
        }

        stream.read((char*)p, p_len);
        return stream.tellg();
    }

public:
    PBWT pbwt;
    std::shared_ptr<RangeCoder>   range_coder; // shared range coder
    std::shared_ptr<GeneralModel> mref; // reference symbol for RLE encoding
    std::shared_ptr<GeneralModel> mlog_rle; // log2(run length)
    std::shared_ptr<GeneralModel> mrle, mrle2_1, mrle2_2, mrle4_1, mrle4_2, mrle4_3, mrle4_4; // rle models for 1-byte, 2-byte, or 4-byte run lengths
    std::shared_ptr<GeneralModel> dirty_wah; // Dirty bitmap words
    std::shared_ptr<GeneralModel> mtype; // Archetype encoding: either bitmap or RLE
    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1; // capacity (memory allocated), flag for data ownership
    uint32_t n_variants; // number of variants encoded
};

/*======   Haplotype context model   ======*/

// Container for haplotype context models
struct djn_ctx_model_container_t {
private:
    // 0000,0001,0010,0011
    // 0100,0101,0110,0111
    // etc.
    static constexpr const uint32_t ref_bits[16] = 
        {0,         16843009,  33686018,  50529027, 
         67372036,  84215045,  101058054, 117901063, 
         134744072, 151587081, 168430090, 185273099, 
         202116108, 218959117, 235802126, 252645135};

public:
    djn_ctx_model_container_t(int64_t n_s, int pl, bool use_pbwt);
    djn_ctx_model_container_t(int64_t n_s, int pl, bool use_pbwt, uint8_t* src, uint32_t src_len);
    ~djn_ctx_model_container_t();

    // Delete move and copy ctors
    djn_ctx_model_container_t(const djn_ctx_model_container_t& other) = delete;
    djn_ctx_model_container_t(djn_ctx_model_container_t&& other) = delete;
    djn_ctx_model_container_t& operator=(const djn_ctx_model_container_t& other) = delete;
    djn_ctx_model_container_t& operator=(djn_ctx_model_container_t&& other) = delete;
    
    void StartEncoding(bool use_pbwt, bool reset = false);
    size_t FinishEncoding();
    void StartDecoding(bool use_pbwt, bool reset = false);

    inline void ResetBitmaps() { memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t)); }

    int DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len);
    int DecodeNextRaw(uint8_t* data, uint32_t& len);

    // Read/write
    int Serialize(uint8_t* dst) const {
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

    int Serialize(std::ostream& stream) const {
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

    int GetSerializedSize() const {
        int ret = sizeof(int) + 3*sizeof(uint32_t) + p_len + model_2mc->GetSerializedSize() + model_nm->GetSerializedSize();
        return ret;
    }

    int GetCurrentSize() const {
        int ret = range_coder->OutSize();
        ret += model_2mc->GetCurrentSize();
        ret += model_nm->GetCurrentSize();
        return ret;
    }
    
    int Deserialize(uint8_t* dst) {
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

    int Deserialize(std::istream& stream) {
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

public:
    int Encode2mc(uint8_t* data, uint32_t len);
    int Encode2mc(uint8_t* data, uint32_t len, const uint8_t* map, const int shift = 1);
    int EncodeNm(uint8_t* data, uint32_t len);
    int EncodeNm(uint8_t* data, uint32_t len, const uint8_t* map, const int shift = 1);
    int EncodeWah(uint32_t* wah, uint32_t len);
    int EncodeWahNm(uint32_t* wah, uint32_t len);
    int EncodeWahRLE(uint32_t ref, uint32_t len, std::shared_ptr<djn_ctx_model_t> model);
    int EncodeWahRLE_nm(uint32_t ref, uint32_t len, std::shared_ptr<djn_ctx_model_t> model);

public:
    int DecodeRaw(uint8_t* data, uint32_t& len);
    int DecodeRaw_nm(uint8_t* data, uint32_t& len);
    int DecodeWahRLE(uint32_t& ref, uint32_t& len, std::shared_ptr<djn_ctx_model_t> model);
    int DecodeWahRLE_nm(uint32_t& ref, uint32_t& len, std::shared_ptr<djn_ctx_model_t> model);

public:
    bool use_pbwt;

    int ploidy;
    int64_t n_samples; // number of "samples" = haplotypes
    int64_t n_variants;

    // Fields for constructing EWAH encodings from input data.
    int64_t n_samples_wah; // Number of samples rounded up to closes 32-bit boundary
    int64_t n_samples_wah_nm;
    uint64_t n_wah; // Number of allocated 32-bit bitmaps
    uint32_t* wah_bitmaps; // Bitmaps

    uint8_t* p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1; // allocated data length, ownership of data flag

    // Shared range coder: All context models share this range coder and emit
    // encodings to a shared buffer.
    std::shared_ptr<RangeCoder> range_coder;
    // Context model encoding the model archetype: either zero (0) for 2MC
    // or one (1) for 2M, or two (2) for everything else.
    // This information is required to differentiate
    std::shared_ptr<GeneralModel> marchetype; // 0 for 2MC, 2 else
    std::shared_ptr<djn_ctx_model_t> model_2mc;
    std::shared_ptr<djn_ctx_model_t> model_nm;

    // Supportive array for computing allele counts to determine the presence
    // of missing values and/or end-of-vector symbols (in Bcf-encodings).
    uint32_t hist_alts[256];
};

class djinn_ctx_model : public djinn_model {
public:
    djinn_ctx_model();
    ~djinn_ctx_model();

    /**
     * Compress and encode data provided as a Bcf-encoded input vector of bytes.
     * The length of the data must be divisible by the ploidy to ascertain
     * correctness. Ploidy number is equal to the stride size in the Bcf-format.
     * PBWT-permutation is enabled by default.
     * 
     * If the provided (data length, ploidy)-tuple has not been observed before
     * then we will treat this as a new compression group. Using this tuple instead
     * of data length alone allows for us to distinguish between data lengths with
     * different ploidy, e.g. (120, 4) and (120, 2) and (120, 1). Although this
     * behaviour is illegal in Vcf/Bcf (different sample numbers per site) it is 
     * perfectly legal in Djinn. There is currently no provided subroutines for
     * maintaining this sample-site relationships so this needs to be implemented by
     * the end-user, if desired.
     * 
     * @param data 
     * @param len_data 
     * @param ploidy 
     * @param alt_alleles 
     * @return int 
     */
    int EncodeBcf(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) override;
    int Encode(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) override;

    /**
     * Builds all the neccessary objects before passing data for encoding. 
     * Calling this function is REQUIRED prior to calling EncodeBcf or Encode
     * on your target data.
     * 
     * @param use_pbwt 
     * @param reset 
     */
    void StartEncoding(bool use_pbwt, bool reset = false) override;
    size_t FinishEncoding() override;
    int StartDecoding() override;

    // Read/write
    int Serialize(uint8_t* dst) const override;
    int Serialize(std::ostream& stream) const override;
    int GetSerializedSize() const;
    int GetCurrentSize() const;
    
    // Deserialize data from an external buffer.
    int Deserialize(uint8_t* src) override;

    // Deserialize data from a IO stream. This approach is more efficient
    // as data does NOT have to be copied and the current model can be
    // reused.
    int Deserialize(std::istream& stream) override;

public:
    int DecodeNext(djinn_variant_t*& variant) override;
    int DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len) override;
    int DecodeNextRaw(uint8_t* data, uint32_t& len) override;

public:
    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1; // allocated data length, ownership of data flag

    // Support buffer. Currently only used for decoding EWAH.
    uint8_t *q;     // data
    uint32_t q_len; // data length
    uint32_t q_alloc:31, q_free:1; // allocated data length, ownership of data flag

    // Shared range coder: All context models share this range coder and emit
    // encodings to a shared buffer.
    std::shared_ptr<RangeCoder> range_coder;
    std::shared_ptr<GeneralModel> ploidy_dict; // 0 for first dict encoding, 1 for second etc.
    std::unordered_map<uint64_t, uint32_t> ploidy_map; // hash table that maps (data length, ploidy) packed into a 64-bit word to model offsets
    std::unordered_map<uint64_t, uint32_t> ploidy_remap; // use to remap in the case when reset is set to false
    std::vector< std::shared_ptr<djn_ctx_model_container_t> > ploidy_models;
};

}
#endif