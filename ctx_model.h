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
#include <unordered_map>//std::unordered_map
#include <limits>//std::numeric_limit

#include "djinn.h"
#include "pbwt.h"

namespace djinn {

static 
uint32_t round_log2(uint32_t x) {
	uint32_t r = 0;
	for (/**/; x; ++r) x >>= 1;
	return r;
}

/*======   Context model container   ======*/

class GeneralModel {
public:
    GeneralModel() noexcept;
    GeneralModel(int n_symbols, int model_size);
    GeneralModel(int n_symbols, int model_size, std::shared_ptr<RangeCoder> rc);
    GeneralModel(int n_symbols, int model_size, int shift, int step);
    GeneralModel(int n_symbols, int model_size, int shift, int step, std::shared_ptr<RangeCoder> rc);
    ~GeneralModel();

    int Initiate(int n_symbols, int model_size);
    int Initiate(int n_symbols, int model_size, std::shared_ptr<RangeCoder> rc);
    int Initiate(int n_symbols, int model_size, int shift, int step);
    int Initiate(int n_symbols, int model_size, int shift, int step, std::shared_ptr<RangeCoder> rc);

    int FinishEncoding();
    int FinishDecoding();
    void StartEncoding();
    void StartDecoding(uint8_t* data);

    void EncodeSymbol(const uint16_t symbol);
    void EncodeSymbolNoUpdate(const uint16_t symbol);
    
    uint16_t DecodeSymbol();
    uint16_t DecodeSymbolNoUpdate();

    void ResetModels();
    void ResetContext();
    void Reset();

public:
    int max_model_symbols;
    int model_context_shift;
    uint32_t model_context, model_ctx_mask;
    std::shared_ptr<RangeCoder> range_coder;
    std::vector < std::shared_ptr<FrequencyModel> > models;
    size_t n_additions; // number of updates performed
    size_t n_buffer; // buffer size
    uint8_t* buffer; // buffer. todo: fixme
};

/*
Desired interface: 

EncodeBcf(); // input data formatted according to bcf spec
Encode() // any input data properly formatted
Decode() // get a proper variant back
DecodeRaw() // get uncompressed but encoded and unpermuted data back

StartEncoding()
StopEncoding()
StartDecoding()
StopDecoding() // does nothing
*/

/*======   GT context variant   ======*/

// Todo
struct djn_ctx_variant_t {
    uint8_t* types;
    uint64_t* data;
    uint32_t n_types, n_data;
    uint32_t m_types: 31, m_types_free: 1; 
    uint32_t m_data: 31, m_data_free: 1;
};

/*======   GT context model type   ======*/

struct djinn_ctx_model_t {
public:
    djinn_ctx_model_t();
    ~djinn_ctx_model_t();

    // Initiate this model as accepting either biallelic complete
    // input data (Initiate2mc) or any other (InitiateNm). Running
    // either of these functions are REQUIRED before starting encoding
    // or decoding.
    void Initiate2mc();
    void InitiateNm();

    int StartEncoding(bool use_pbwt, bool reset = false);
    size_t FinishEncoding();
    int StartDecoding(uint8_t* data, bool reset = false);
    size_t FinishDecoding() { return 0; } // no effect
    
    void reset();

public:
    uint8_t use_pbwt: 1, init: 1, unused: 6;
    PBWT pbwt;
    std::shared_ptr<RangeCoder> range_coder; // shared range coder
    std::shared_ptr<GeneralModel> mref; // reference symbol for RLE encoding
    std::shared_ptr<GeneralModel> mlog_rle;  // log2(run length)
    std::shared_ptr<GeneralModel> mrle, mrle2_1, mrle2_2, mrle4_1, mrle4_2, mrle4_3, mrle4_4; // rle models for 1-byte, 2-byte, or 4-byte run lengths
    std::shared_ptr<GeneralModel> dirty_wah; // Dirty bitmap words
    std::shared_ptr<GeneralModel> mtype; // Archetype encoding: either bitmap or RLE
    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1; // capacity (memory allocated), flag for data ownership
    uint32_t n_variants; // number of variants encoded
};

/*======   Haplotype base model   ======*/

class djinn_model {
public:
    djinn_model() : n_samples(0), n_variants(0) {}
    djinn_model(uint64_t n_s) : n_samples(n_s), n_variants(0) {}
    ~djinn_model() {}

    virtual void SetSamples(int64_t n_s) { n_samples = n_s; }

    virtual int EncodeBcf(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles, const bool permute) =0;

    virtual void StartEncoding(bool use_pbwt, bool reset = false) =0;
    virtual size_t FinishEncoding() =0;
    virtual int StartDecoding(djinn_block_t* block) =0;

    // Decoding requires:
    // One buffer for decoding into EWAH values.
    // One buffer of size no smaller than ploidy*samples for storing the return vector of alleles.
    virtual int DecodeNext(uint8_t* data, size_t& len) =0;
    virtual int DecodeNext(uint8_t* ewah_data, size_t& ret_ewah, uint8_t* ret_buffer, size_t& ret_len) =0;
    virtual int DecodeNextRaw(uint8_t* data, size_t& len) =0;
    
    // Retrieval.
    virtual int GetBlockReference(djinn_block_t*& block) =0;
    
    // Todo
    // virtual int GetBlock(djinn_block_t*& block) =0;

public:
    int64_t n_samples; // number of samples or haplotypes
    uint32_t n_variants; // number of encoded variants
    // Supportive array for computing allele counts to determine the presence
    // of missing values and/or end-of-vector symbols (in Bcf-encodings).
    uint32_t hist_alts[256];
};

/*======   Haplotype context model   ======*/

// Container for haplotype context models
struct djinn_ctx_model_container_t {
public:
    djinn_ctx_model_container_t(int64_t n_s, int pl, bool use_pbwt) : 
        ploidy(pl), n_samples(n_s),
        n_wah(std::ceil((float)n_samples / 32) * 4), 
        n_samples_wah((n_wah * 32) / 4), 
        wah_bitmaps(new uint32_t[n_wah]), 
        p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
        range_coder(std::make_shared<RangeCoder>()), 
        marchetype(std::make_shared<GeneralModel>(2, 1024, range_coder)),
        model_2mc(std::make_shared<djinn_ctx_model_t>()),
        model_2m(std::make_shared<djinn_ctx_model_t>()),
        model_nm(std::make_shared<djinn_ctx_model_t>())
    {
        assert(n_s % pl == 0); // #samples/#ploidy must be divisible
        model_2mc->Initiate2mc();
        model_nm->InitiateNm();
        if (use_pbwt) {
            model_2mc->pbwt.Initiate(2, n_s);
            model_nm->pbwt.Initiate(16, n_s);
        }
    }

    ~djinn_ctx_model_container_t() {
        if (p_free) {
            delete[] p;
        }
        delete[] wah_bitmaps;
    }

public:
    int ploidy;
    int64_t n_samples; // number of "samples" = haplotypes

    // Fields for constructing EWAH encodings from input data.
    uint64_t n_wah; // Number of allocated 32-bit bitmaps
    int64_t n_samples_wah; // Number of samples rounded up to closes 32-bit boundary
    uint32_t* wah_bitmaps; // Bitmaps

    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1; // allocated data length, ownership of data flag

    // Shared range coder: All context models share this range coder and emit
    // encodings to a shared buffer.
    std::shared_ptr<RangeCoder> range_coder;
    // Context model encoding the model archetype: either zero (0) for 2MC
    // or one (1) for 2M, or two (2) for everything else.
    // This information is required to differentiate
    std::shared_ptr<GeneralModel> marchetype; // 0 for 2MC, 1 for 2M, 2 else
    std::shared_ptr<djinn_ctx_model_t> model_2mc;
    std::shared_ptr<djinn_ctx_model_t> model_2m;
    std::shared_ptr<djinn_ctx_model_t> model_nm;
};

class djinn_ctx_model : public djinn_model {
private:
    // 00000000,00000001,00000010,00000011
    // 00000100,00000101,00000110,00000111
    // etc.
    static constexpr const uint32_t ref_bits[16] = 
        {0,         16843009,  33686018,  50529027, 
         67372036,  84215045,  101058054, 117901063, 
         134744072, 151587081, 168430090, 185273099, 
         202116108, 218959117, 235802126, 252645135};

public:
    djinn_ctx_model();
    djinn_ctx_model(uint64_t n_s);
    ~djinn_ctx_model();

    // Todo: delete this
    void SetSamples(int64_t n_s) override;

    /**
     * Compress and encode data provided as a Bcf-encoded input vector of bytes.
     * The length of the data must be divisible by the ploidy to ascertain
     * correctness. Ploidy number is equal to the stride size in the Bcf-format.
     * PBWT-permutation is enabled by default. Override the "permute" parameter
     * by setting it to false.
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
     * @param permute 
     * @return int 
     */
    int EncodeBcf(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles, const bool permute = true) override;

    void StartEncoding(bool use_pbwt, bool reset = false) override;
    size_t FinishEncoding() override;
    // Implement me
    int StartDecoding(djinn_block_t* block) override;
    
    // Retrieval.
    int GetBlockReference(djinn_block_t*& block) override;
    
    // Todo
    // int GetBlock(djinn_block_t*& block) override;
    
// Internal functions
public:
    int Encode2mc(uint8_t* data, uint32_t len);
    int EncodeNm(uint8_t* data, uint32_t len);
    int EncodeWah(uint32_t* wah, uint32_t len);
    int EncodeWahNm(uint32_t* wah, uint32_t len);
    int EncodeWahRLE(uint32_t ref, uint32_t len, std::shared_ptr<djinn_ctx_model_t> model);
    int EncodeWahRLE_nm(uint32_t ref, uint32_t len, std::shared_ptr<djinn_ctx_model_t> model);

public:
    int DecodeNext(uint8_t* data, size_t& len) override;
    int DecodeNext(uint8_t* ewah_data, size_t& ret_ewah, uint8_t* ret_buffer, size_t& ret_len) override;
    int DecodeNextRaw(uint8_t* data, size_t& len) override;

    // int DecodeNextRaw(uint8_t* data);
    int DecodeRaw(uint8_t* data, size_t& len);
    int DecodeRaw_nm(uint8_t* data, size_t& len);
    int DecodeWahRLE(uint32_t& ref, uint32_t& len, std::shared_ptr<djinn_ctx_model_t> model);
    int DecodeWahRLE_nm(uint32_t& ref, uint32_t& len, std::shared_ptr<djinn_ctx_model_t> model);

public:
    // Fields for constructing EWAH encodings from input data.
    uint64_t n_wah; // Number of allocated 32-bit bitmaps
    int64_t n_samples_wah; // Number of samples rounded up to closes 32-bit boundary
    uint32_t* wah_bitmaps; // Bitmaps

    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1; // allocated data length, ownership of data flag

    // Shared range coder: All context models share this range coder and emit
    // encodings to a shared buffer.
    std::shared_ptr<RangeCoder> range_coder;
    // Context model encoding the model archetype: either zero (0) for 2MC
    // or one (1) for 2M, or two (2) for everything else.
    // This information is required to differentiate
    std::shared_ptr<GeneralModel> marchetype; // 0 for 2MC, 1 for 2M, 2 else
    std::shared_ptr<djinn_ctx_model_t> model_2mc;
    std::shared_ptr<djinn_ctx_model_t> model_2m;
    std::shared_ptr<djinn_ctx_model_t> model_nm;

    // Todo: these objects should be placed in a container
    // to serve arbitrary N-ploidy
    // Todo: map from offset to ploidy (ploidy here is interpreted as number of haplotypes)
    std::shared_ptr<GeneralModel> ploidy_dict; // 0 for first dict encoding, 1 for second etc.
    std::unordered_map<uint64_t, uint32_t> ploidy_map; // maps (data length, ploidy) packed into a 64-bit word to model offsets
    std::vector< std::shared_ptr<djinn_ctx_model_container_t> > ploidy_models;
};

}