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
#ifndef DJINN_EWAH_MODEL_H_
#define DJINN_EWAH_MODEL_H_
#include <unordered_map> //std::unordered_map
#include <limits> //std::numeric_limit

#include "djinn.h" // core components
#include "pbwt.h" // PBWT algorithms

namespace djinn {

// Compression strategy used.
enum class CompressionStrategy : uint32_t { ZSTD = 0, LZ4 = 1 };

struct djn_ewah_model_t {
public:
    djn_ewah_model_t();
    ~djn_ewah_model_t();

    int StartEncoding(bool use_pbwt, bool reset = false);
    size_t FinishEncoding(uint8_t* support_buffer, uint32_t support_cap, CompressionStrategy strat);
    int StartDecoding(bool use_pbwt, bool reset = false);
    size_t FinishDecoding() { return 0; } // no effect
    
    void reset();

    // Read/write
    int Serialize(uint8_t* dst) const;
    int Serialize(std::ostream& stream) const;
    int Deserialize(uint8_t* dst);
    int Deserialize(std::istream& stream);

public:
    PBWT pbwt;
    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1; // capacity (memory allocated), flag for data ownership
    uint32_t n_variants; // number of variants encoded
};

struct djn_ewah_model_container_t {
public:
    djn_ewah_model_container_t(int64_t n_s, int pl, bool use_pbwt);
    djn_ewah_model_container_t(int64_t n_s, int pl, bool use_pbwt, uint8_t* src, uint32_t src_len);
    ~djn_ewah_model_container_t();

    void StartEncoding(bool use_pbwt, bool reset = false);
    size_t FinishEncoding(uint8_t* support_buffer, uint32_t support_cap, CompressionStrategy strat);
    void StartDecoding(bool use_pbwt, bool reset = false);

    inline void ResetBitmaps() { memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t)); }

    int DecodeNext(uint8_t* data, uint32_t& len);
    int DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len);
    int DecodeNextRaw(uint8_t* data, uint32_t& len);

public:
    int Encode2mc(uint8_t* data, uint32_t len);
    int Encode2mc(uint8_t* data, uint32_t len, const uint8_t* map, const int shift = 1);
    int EncodeNm(uint8_t* data, uint32_t len);
    int EncodeNm(uint8_t* data, uint32_t len, const uint8_t* map, const int shift = 1);
    int EncodeWah(uint32_t* wah, uint32_t len);
    int EncodeWahNm(uint32_t* wah, uint32_t len);

public:
    int DecodeRaw(uint8_t* data, uint32_t& len);
    int DecodeRaw_nm(uint8_t* data, uint32_t& len);

public:
    bool use_pbwt;
    int ploidy;
    int64_t n_samples; // number of "samples" = haplotypes
    int64_t n_variants;

    // Fields for constructing EWAH encodings from input data.
    uint64_t n_wah; // Number of allocated 32-bit bitmaps
    int64_t n_samples_wah; // Number of samples rounded up to closes 32-bit boundary
    int64_t n_samples_wah_nm;
    uint32_t* wah_bitmaps; // Bitmaps

    uint8_t* p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1; // allocated data length, ownership of data flag
    
    // std::shared_ptr<GeneralModel> marchetype; // 0 for 2MC, 2 else
    std::shared_ptr<djn_ewah_model_t> model_2mc;
    std::shared_ptr<djn_ewah_model_t> model_2m; // unused
    std::shared_ptr<djn_ewah_model_t> model_nm;

    // Supportive array for computing allele counts to determine the presence
    // of missing values and/or end-of-vector symbols (in Bcf-encodings).
    uint32_t hist_alts[256];
};

class djinn_ewah_model : public djinn_model {
public:
    djinn_ewah_model();
    ~djinn_ewah_model();

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
    int Serialize(uint8_t* dst) const override { return -1; };
    int Serialize(std::ostream& stream) const override { return -1; };
    int Deserialize(uint8_t* src) override { return -1; };
    // Deserialize data from a IO stream. This approach is more efficient
    // as data does NOT have to be copied and the current model can be
    // reused.
    int Deserialize(std::istream& stream) override { return -1; };

public:
    int DecodeNext(djinn_variant_t*& variant) override;
    int DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len) override;
    int DecodeNextRaw(uint8_t* data, uint32_t& len) override;

public:
    CompressionStrategy codec; // Either ZSTD or LZ4 at the moment.

    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1; // allocated data length, ownership of data flag

    // Support buffer. Currently only used for decoding EWAH.
    uint8_t *q;     // data
    uint32_t q_len; // data length
    uint32_t q_alloc:31, q_free:1; // allocated data length, ownership of data flag

    std::shared_ptr<GeneralModel> ploidy_dict; // 0 for first dict encoding, 1 for second etc.
    std::unordered_map<uint64_t, uint32_t> ploidy_map; // maps (data length, ploidy) packed into a 64-bit word to model offsets
    std::vector< std::shared_ptr<djn_ewah_model_container_t> > ploidy_models;
};

}
#endif