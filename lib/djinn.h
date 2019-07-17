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
#ifndef DJINN_H_
#define DJINN_H_

#include <cstddef>//size_t
#include <cstdint>//uint
#include <cmath>//ceil,floor
#include <cstring>//memcpy
#include <cassert>//assert

#include <iostream>//debug
#include <bitset>//debug

#include <unordered_map> //std::unordered_map
#include <limits> //std::numeric_limit
#include <vector> //std::Vector
#include <memory> //std::shared_ptr

namespace djinn {

/**
 * Naming convention: djinn_* structures/classes are considered front-end
 * and djn_* objects are considered back-end.
 */

/*------   Version   ------*/
const int32_t DJINN_VERSION_MAJOR = 0;
const int32_t DJINN_VERSION_MINOR = 1;
const int32_t DJINN_VERSION_PATCH = 0;
const int32_t DJINN_VERSION_NUMBER  = (DJINN_VERSION_MAJOR *100*100 + DJINN_VERSION_MINOR *100 + DJINN_VERSION_PATCH);
const std::string DJINN_LIB_VERSION = std::to_string(DJINN_VERSION_MAJOR) + '.' + std::to_string(DJINN_VERSION_MINOR) + '.' + std::to_string(DJINN_VERSION_PATCH);

/* *************************************
 *  Default constant
 ***************************************/
#ifndef DJINN_CLEVEL_DEFAULT
#  define DJINN_CLEVEL_DEFAULT 3
#endif

/* *************************************
 *  Genotype maps
 ***************************************/

// Map HtsLib-styled Bcf-encoding to our internal format. This map is used for
// biallelic variants with support for missing values and the special end-of-vector
// (EOV) symbol.
// Maps missing to 2, 1->0, 2->1, and EOV -> 3.
constexpr const uint8_t DJN_BCF_GT_UNPACK[65] = 
{2,0,1,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3};

// Map HtsLib-styled Bcf-encoding to our internal format. This map is used for
// multi-allelic variants with support for missing values and the special 
// end-of-vector (EOV) symbol.
// MISSING -> 14, EOV -> 15, other values as normal
// Note that currently we can only store up to 14 ALT alleles in addition to 
// missing values and the EOV symbol.
constexpr const uint8_t DJN_BCF_GT_UNPACK_GENERAL[256] = 
{14, // missing value maps to 14
0,1,2,3,4,5,6,7,8,9,10,11,12,13,16,17,18,19,20,21,22,23,24,25,
26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,
47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,
15, // EOV value (64) maps to 15
66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,
87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,
106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,
121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,
137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,
153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,
169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,
185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,
201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,
217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,
233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,
249,250,251,252,253,254,255};

// Identity map. Maps 0->0, 1->1, 2->2, 3->3, etc.
constexpr const uint8_t DJN_MAP_NONE[256] = 
{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,
81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,
101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,
116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,
131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,
146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,
161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,
191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,
206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,
221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,
236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,
251,252,253,254,255};

// Packed 4-bit symbols into 32-bit words. Used for checking for symbol
// uniformity in a bitvector when the symbol (stride) size is > 1.
// S  Integer    Symbol repeat (4-bit)
// 0  0          00000000000000000000000000000000
// 1  286331153  00010001000100010001000100010001
// 2  572662306  00100010001000100010001000100010
// 3  858993459  00110011001100110011001100110011
// 4  1145324612 01000100010001000100010001000100
// 5  1431655765 01010101010101010101010101010101
// 6  1717986918 01100110011001100110011001100110
// 7  2004318071 01110111011101110111011101110111
// 8  2290649224 10001000100010001000100010001000
// 9  2576980377 10011001100110011001100110011001
// 10 2863311530 10101010101010101010101010101010
// 11 3149642683 10111011101110111011101110111011
// 12 3435973836 11001100110011001100110011001100
// 13 3722304989 11011101110111011101110111011101
// 14 4008636142 11101110111011101110111011101110
// 15 4294967295 11111111111111111111111111111111
static constexpr uint32_t DJN_NM_REF_BITS[16] = 
    {0,          286331153,  572662306,  858993459, 
    1145324612, 1431655765, 1717986918, 2004318071, 
    2290649224, 2576980377, 2863311530, 3149642683, 
    3435973836, 3722304989, 4008636142, 4294967295};

#define DJN_BCF_UNPACK_GENOTYPE(A) DJN_BCF_GT_UNPACK[(A) >> 1]
#define DJN_BCF_UNPACK_GENOTYPE_GENERAL(A) DJN_BCF_GT_UNPACK_GENERAL[(A) >> 1]

#define DJN_UN_NONE 0 // No unpacking
#define DJN_UN_EWAH 1 // EWAH level
#define DJN_UN_IND  2 // Individual level

#define DJN_DIRTY_2MC 0 // 1-bit or
#define DJN_DIRTY_NM  1 // 4-bit encoding in dirty bitmaps

// EWAH structure
#pragma pack(push, 1)
struct djinn_ewah_t {
    djinn_ewah_t() : ref(0), clean(0), dirty(0){}
    void reset() { ref = clean = dirty = 0; }

    uint64_t ref: 4, clean: 30, dirty: 30;
};
#pragma pack(pop)

struct djn_variant_dec_t {
    djn_variant_dec_t() : m_ewah(0), m_dirty(0), n_ewah(0), n_dirty(0), n_samples(0), dirty_type(0), ewah(nullptr), dirty(nullptr) {}
    ~djn_variant_dec_t() {
        delete[] ewah;
        delete[] dirty;
    }

    /**
     * Allocates 'n' elements worth of data for both ewah and dirty pointers.
     * Resizing of 'ewah' and 'dirty' may occur in some functions however this
     * subroutine MUST be called prior to usage. This is normally handled
     * internally.
     * 
     * @param n Allocate this many elements worth of memory.
     */
    void Allocate(uint32_t n) {
        delete[] ewah; delete[] dirty;
        m_ewah = n;
        m_dirty = n;
        n_ewah = 0, n_dirty = 0;
        n_samples = 0;
        ewah = new djinn_ewah_t*[m_ewah];
        dirty = new uint32_t*[m_dirty];
    }

    int m_ewah, m_dirty; // allocated bytes for ewah (m_ewah) or dirty (m_dirty) pointers
    int n_ewah, n_dirty; // number of elements for ewah (n_ewah) or dirty (n_dirty)
    int n_samples; // number of samples (alleles)
    int dirty_type; // one of DJN_DIRTY_*
    djinn_ewah_t **ewah; // pointers to ewah structs in parent djinn_variant_t::data
    uint32_t **dirty; // pointers to dirty uint32_t words in parent djinn_variant_t::data
};

// Variant record
struct djinn_variant_t {
    djinn_variant_t();
    ~djinn_variant_t();

    // Todo: implement me
    // int ToBcf(char* out, const char phasing = '|') const;

    /**
     * Silly function that adds a trailing "\t\n". This output is used in 
     * debugging checksums when comparing against the output of
     * bcftools query -f "[%GT\t]\n".
     *  
     * @param out     Output buffer.
     * @param phasing Character to use as phasing separator if ploidy > 1.
     * @return int    Returns the number of bytes used.
     */
    int ToVcfDebug(char* out, const char phasing = '|') const;
    
    /**
     * Convert data into a valid VCF string irrespective of wether the internally
     * stored data is EWAH compressed or unpacked into byte literals. Note that
     * if raw EWAH data is used (for example as retrieved from the GetNextRaw function)
     * it is possible that the output VCF order is permuted!
     * 
     * No overflow checks are made for the output buffer and the buffer must be
     * pre-allocated by the user.
     * 
     * @param out     Output buffer.
     * @param phasing Character to use as phasing separator if ploidy > 1.
     * @return int    Returns the number of bytes used.
     */
    int ToVcf(char* out, const char phasing = '|') const;

    // Todo: implement me
    int Merge(const djinn_variant_t& other) {
        if (data == nullptr && other.data == nullptr) return 0;
        if (errcode || other.errcode) return -1*errcode;
        if (unpacked != other.unpacked) return -1;

        // Resize if required
        if (data_len + other.data_len > data_alloc) {
            uint8_t* old = data;
            data_alloc = data_len + other.data_len + 65536;
            data = new uint8_t[data_alloc];
            memcpy(data, old, data_len);
            if (data_free) delete[] old;
            data_free = true;
        }

        return 1;
    }

    int ploidy, n_allele; // ploidy: data stride size for unpacked data, n_allele: number of alleles including ref
    uint8_t* data;
    uint32_t data_len;
    uint32_t data_alloc: 31, data_free: 1;
    int errcode; // error code when something goes wrong
    int unpacked; // one of DJN_UN_*
    djn_variant_dec_t* d;
};

/*======  Helper functions  ======*/

/**
 * Supportive function for converting an input variable 
 * x to floor(log2(x)).
 * 
 * @param x 
 * @return uint32_t 
 */
static uint32_t round_log2(uint32_t x) {
	uint32_t r = 0;
	for (/**/; x; ++r) x >>= 1;
	return r;
}

/*======   Base interface for Djinn   ======*/

/***************************************
*  Simple API
***************************************/
class djinn_model {
public:
    djinn_model() : use_pbwt(true), init(true), unused(0), n_variants(0) {}
    virtual ~djinn_model() {}

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
    virtual int EncodeBcf(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) =0;
    /**
     * Same as EncodeBcf but accepts [0,N-1]-encoded byte vectors.
     */
    virtual int Encode(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) =0;

    /**
     * Builds all the neccessary objects before passing data for encoding. 
     * Calling this function is REQUIRED prior to calling EncodeBcf or Encode
     * on your target data.
     * 
     * @param use_pbwt 
     * @param reset 
     */
    virtual void StartEncoding(bool use_pbwt, bool reset) =0;
    
    /**
     * Finishes the current encoding. Calling this function is REQUIRED prior to
     * calling Serialize or starting to Decode data as additional information is
     * written.
     * 
     * @return size_t 
     */
    virtual size_t FinishEncoding() =0;

    /**
     * Starts encoding the internal data that MUST have been provided prior to calling
     * this function. Easiest way to correctly read data is to Deserialize() from a
     * IO-stream or a preloaded buffer.
     * 
     * @return int 
     */
    virtual int StartDecoding() =0;

    // Decoding requires:
    // One buffer for decoding into EWAH values.
    // One buffer of size no smaller than ploidy*samples for storing the return vector of alleles.
    virtual int DecodeNext(djinn_variant_t*& variant) =0;
    virtual int DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len) =0;
    
    virtual int DecodeNextRaw(djinn_variant_t*& variant) =0;
    virtual int DecodeNextRaw(uint8_t* data, uint32_t& len) =0;

    // Read write
    virtual int Serialize(uint8_t* dst) const =0;
    virtual int Serialize(std::ostream& stream) const =0;
    
    /**
     * Deserialize a djinn-model-derived object from a given buffer. This
     * function involves a copy of the provided data internally.
     * 
     * @param dst 
     * @return int 
     */
    virtual int Deserialize(uint8_t* dst) =0;
    
    /**
     * Deserialize a djinn_model-derived object from an IO-stream. This approach
     * is generally much more efficient as no copies are involved and the current
     * model can be reused each iteration.
     * 
     * @param stream 
     * @return int 
     */
    virtual int Deserialize(std::istream& stream) =0;

    /**
     * Support function that computes the size of this object when serialized.
     * Used internally when calling the Serialize() functions. Is also useful
     * for knowing how much memory to allocate if serializing to an external 
     * buffer.
     * 
     * @return int 
     */
    virtual int GetSerializedSize() const =0;
    
    /**
     * Returns the current size of *encodings*.
     */
    virtual int GetCurrentSize() const =0;

    // Todo:
    // virtual int Merge(djinn_model* b1, djinn_model* b2);

public:
    uint8_t use_pbwt: 1, // PBWT pre-processor is used
            init: 1,     // Models should be reset
            unused: 6;   // Reserved space
    uint32_t n_variants; // Number of encoded variants

    // Supportive array for computing allele counts to determine the presence
    // of missing values and/or end-of-vector symbols (in Bcf-encodings).
    uint32_t hist_alts[256];
};

/***************************************
*  Context model
***************************************/
// Forward declarations.
class RangeCoder; // Range Coder
class GeneralModel;// General Model for Context Modelling
class PBWT; // PBWT

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
    int Serialize(uint8_t* dst) const;
    int Serialize(std::ostream& stream) const;
    int GetSerializedSize() const;
    int GetCurrentSize() const;
    int Deserialize(uint8_t* dst);
    int Deserialize(std::istream& stream);

public:
    std::shared_ptr<PBWT> pbwt;
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
    int DecodeNextRaw(djinn_variant_t*& variant);

    // Read/write
    int Serialize(uint8_t* dst) const;
    int Serialize(std::ostream& stream) const;
    int GetSerializedSize() const;
    int GetCurrentSize() const;
    int Deserialize(uint8_t* dst);
    int Deserialize(std::istream& stream);

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

    int EncodeBcf(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) override;
    int Encode(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) override;

    void StartEncoding(bool use_pbwt, bool reset = false) override;
    size_t FinishEncoding() override;
    int StartDecoding() override;

    // Read/write
    int Serialize(uint8_t* dst) const override;
    int Serialize(std::ostream& stream) const override;
    int Deserialize(uint8_t* src) override;
    int Deserialize(std::istream& stream) override;
    int GetSerializedSize() const override;
    int GetCurrentSize() const override;

public:
    int DecodeNext(djinn_variant_t*& variant) override;
    int DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len) override;
    int DecodeNextRaw(uint8_t* data, uint32_t& len) override;
    int DecodeNextRaw(djinn_variant_t*& variant) override;

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

/***************************************
*  EWAH model
***************************************/
// Compression strategy used.
enum class CompressionStrategy : uint32_t { ZSTD = 0, LZ4 = 1 };

struct djn_ewah_model_t {
public:
    djn_ewah_model_t();
    ~djn_ewah_model_t();

    int StartEncoding(bool use_pbwt, bool reset = false);
    size_t FinishEncoding(uint8_t*& support_buffer, uint32_t& support_cap, CompressionStrategy strat, int c_level);
    int StartDecoding(uint8_t*& support_buffer, uint32_t& support_cap, CompressionStrategy strat, bool use_pbwt, bool reset = false);
    size_t FinishDecoding() { return 0; } // no effect
    
    void reset();

    // Read/write
    int Serialize(uint8_t* dst) const;
    int Serialize(std::ostream& stream) const;
    int GetSerializedSize() const;
    int Deserialize(uint8_t* dst);
    int Deserialize(std::istream& stream);

public:
    std::shared_ptr<PBWT> pbwt;
    uint8_t *p;     // data
    uint32_t p_len, u_len; // data length
    uint32_t p_cap:31, p_free:1; // capacity (memory allocated), flag for data ownership
    uint32_t n_variants; // number of variants encoded
};

struct djn_ewah_model_container_t {
public:
    djn_ewah_model_container_t(int64_t n_s, int pl, bool use_pbwt);
    djn_ewah_model_container_t(int64_t n_s, int pl, bool use_pbwt, uint8_t* src, uint32_t src_len);
    ~djn_ewah_model_container_t();

    void StartEncoding(bool use_pbwt, bool reset = false);
    size_t FinishEncoding(uint8_t*& support_buffer, uint32_t& support_cap, CompressionStrategy strat, int c_level);
    void StartDecoding(uint8_t*& support_buffer, uint32_t& support_cap, CompressionStrategy strat, bool use_pbwt, bool reset = false);

    inline void ResetBitmaps() { memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t)); }

    int DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len);
    int DecodeNextRaw(uint8_t* data, uint32_t& len);
    int DecodeNextRaw(djinn_variant_t*& variant);

    // Read/write
    int Serialize(uint8_t* dst) const;
    int Serialize(std::ostream& stream) const;
    int GetSerializedSize() const;
    int GetCurrentSize() const;
    int Deserialize(uint8_t* dst);
    int Deserialize(std::istream& stream);

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
    int64_t n_samples_wah; // Number of samples rounded up to closes 32-bit boundary
    int64_t n_samples_wah_nm;
    uint64_t n_wah; // Number of allocated 32-bit bitmaps
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
    djinn_ewah_model(CompressionStrategy codec, int c_level = 1);
    ~djinn_ewah_model();

    int EncodeBcf(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) override;
    int Encode(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) override;

    void StartEncoding(bool use_pbwt, bool reset = false) override;
    size_t FinishEncoding() override;
    int StartDecoding() override;

    // Read/write
    int Serialize(uint8_t* dst) const override;
    int Serialize(std::ostream& stream) const override;
    int Deserialize(uint8_t* src) override;
    int Deserialize(std::istream& stream) override;
    int GetSerializedSize() const override;
    int GetCurrentSize() const override;

    // Todo: Merge EWAH data pairwise.
    // If the data is PBWT-permuted then first unpermute and add
    // raw data together.
    int Merge(djinn_ewah_model& model1, djinn_ewah_model& model2) {
        std::unordered_map<uint32_t, uint64_t> merge_map; // map ploidy -> new length
        std::unordered_map<uint32_t, uint64_t> model_map; // map ploidy -> offset in tgt_containers
        // Vector of vector of model containers to merge.
        std::vector< std::vector< std::shared_ptr<djn_ewah_model_container_t> > > tgt_containers;

        djinn_ewah_model* models[2] = {&model1, &model2};
        
        // Cycle over pairs of ploidy models in either container.
        // Record the ploidy maps where ploidy in the (ploidy,len)-tuple are the same
        // and add a new container with (ploidy,lenA + lenB).
        for (int i = 0; i < 2; ++i) {
            std::unordered_map<uint64_t, uint32_t>& target_map = models[i]->ploidy_map;
            for (auto it = target_map.begin(); it != target_map.end(); ++it) {
                // it->first  = (ploidy,len)-tuple
                // it->second = offset into ploidy_models
                // const uint64_t tuple = ((uint64_t)len_data << 32) | ploidy;
                uint32_t ploidy = it->first & ((1L << 32) - 1);
                uint32_t len = (it->first >> 32) & ((1L << 32) - 1);
                std::cerr << it->first << "(" << ploidy << "," << len << ") -> " << it->second << std::endl;

                auto search = merge_map.find(ploidy);
                if (search != merge_map.end()) {
                    // Already in map: increment by current length.
                    search->second += len;
                    tgt_containers.back().push_back(models[i]->ploidy_models[it->second]);
                } else {
                    // std::cerr << "Not found. Inserting: " << len_data << "," << ploidy << "(" << tuple << ") as " << ploidy_models.size() << std::endl;
                    merge_map[ploidy] = len;
                    model_map[ploidy] = tgt_containers.size();
                    tgt_containers.push_back(std::vector< std::shared_ptr<djn_ewah_model_container_t> >());
                    tgt_containers.back().push_back(models[i]->ploidy_models[it->second]);
                }
            }
        }

        for (auto it = merge_map.begin(); it != merge_map.end(); ++it) {
            // it->first  = ploidy
            // it->second = new len
            const uint64_t tuple = ((uint64_t)it->second << 32) | it->first;
            ploidy_map[tuple] = ploidy_models.size();

            // Add new containers.
            ploidy_models.push_back(std::make_shared<djn_ewah_model_container_t>(it->second, it->first, (bool)use_pbwt));
            ploidy_models.back()->StartEncoding(use_pbwt, init);
            // tgt_container = ploidy_models[ploidy_models.size() - 1];
            // tgt_container->StartEncoding(use_pbwt, init);
        }

        return 1;
    }

public:
    int DecodeNext(djinn_variant_t*& variant) override;
    int DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len) override;
    int DecodeNextRaw(uint8_t* data, uint32_t& len) override;
    int DecodeNextRaw(djinn_variant_t*& variant) override;

public:
    CompressionStrategy codec; // Either ZSTD or LZ4 at the moment.
    int compression_level;

    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1; // allocated data length, ownership of data flag

    // Support buffer. Currently only used for decoding EWAH.
    uint8_t *q;     // data
    uint32_t q_len; // data length
    uint32_t q_alloc:31, q_free:1; // allocated data length, ownership of data flag

    std::unordered_map<uint64_t, uint32_t> ploidy_map; // maps (data length, ploidy) packed into a 64-bit word to model offsets
    std::vector< std::shared_ptr<djn_ewah_model_container_t> > ploidy_models;
};

/***************************************
*  Sample-centric bitmap model
***************************************/
struct djn_bitmap_model_t {
public:
    djn_bitmap_model_t() : p(nullptr), p_len(0), u_len(0), p_cap(0), p_free(true) {}
    ~djn_bitmap_model_t() { 
        if (p_free) delete[] p;
    }

    // Read/write
    int Serialize(uint8_t* dst) const;
    int Serialize(std::ostream& stream) const;
    int GetSerializedSize() const;
    int Deserialize(uint8_t* dst);
    int Deserialize(std::istream& stream);

public:
    // Data must be aligned to the largest memory boundary available to 
    // ascertain good vectorization performance.
    uint64_t *p;     // bitmaps
    uint32_t p_len, u_len; // data length
    uint32_t p_cap:31, p_free:1; // capacity (memory allocated), flag for data ownership
};

class djinn_bitmap_model {
public:
    djinn_bitmap_model() : 
        codec(CompressionStrategy::LZ4), compression_level(1), 
        q(new uint8_t[10000000]), q_len(0), 
        q_alloc(10000000), q_free(true), n_samples(0), n_samples_bitmap(0), n_variants(0) { }

    djinn_bitmap_model(CompressionStrategy codec, int c_level = 1) : 
        codec(codec), compression_level(c_level), 
        q(new uint8_t[10000000]), q_len(0), 
        q_alloc(10000000), q_free(true), n_samples(0), n_samples_bitmap(0), n_variants(0) { }

    ~djinn_bitmap_model() {
        if (q_free) delete[] q;
    }

    // Only biallelic haploid or diploid samples are supported in this model.
    // This is a concious design decision to maximize computational throughput
    // for the most common use-case.
    int EncodeBcf(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles);
    
    int Encode(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) {
        if (data == nullptr) return -1;
        if (len_data == 0)   return  0;
        if (alt_alleles > 4) return -2;
        if (n_variants == n_samples) return -3;

        // Compute allele counts.
        memset(hist_alts, 0, 256*sizeof(uint32_t));
        for (int i = 0; i < len_data; ++i) {
            ++hist_alts[data[i]];
        }

        if (ploidy == 1) {
            for (int i = 0; i < len_data; ++i) {
                const uint8_t& gt = DJN_BCF_UNPACK_GENOTYPE(data[i]);
                bitmaps[i]->p[n_variants / 64] |= (uint64_t)gt << (n_variants % 64);
            }
        } else if (ploidy == 2) {
                for (int i = 0; i < len_data; ++i) {
                const uint8_t& gt = DJN_BCF_UNPACK_GENOTYPE(data[i]);
                bitmaps[i]->p[n_variants / 32] |= (uint64_t)gt << (2*(n_variants % 32)+0);
                bitmaps[i]->p[n_variants / 32] |= (uint64_t)gt << (2*(n_variants % 32)+1);
            }
        }
        ++n_variants;

        return 1;
    }

    // Allele counts must be either 2 (biallelic, no missing) or 4
    // biallelic (including missing).
    int StartEncoding(uint32_t n_samples, uint32_t n_alleles) {
        if (n_samples == 0) return -1;
        if (this->n_samples != n_samples) {
            this->n_samples = n_samples;
            n_samples_bitmap = std::ceil(n_samples/32.0f)*32;
            bitmaps.resize(n_samples);
            n_variants = 0;
        } else {
            for (int i = 0; i < bitmaps.size(); ++i) 
                bitmaps[i]->p_len = 0;
        }
        return 1;
    }

    int StartDecoding();

    // Read/write
    int Serialize(uint8_t* dst) const;
    int Serialize(std::ostream& stream) const;
    int Deserialize(uint8_t* src);
    int Deserialize(std::istream& stream);
    int GetSerializedSize() const;
    int GetCurrentSize() const;

    // Accessors.
    std::shared_ptr<djn_bitmap_model_t> operator[](const uint32_t p) { return bitmaps[p]; }
    std::shared_ptr<djn_bitmap_model_t> at(const uint32_t p) { return bitmaps[p]; }

public:
    CompressionStrategy codec; // Either ZSTD or LZ4 at the moment.
    int compression_level;

    // Support buffer. Currently only used for using additional compressors.
    uint8_t *q;     // data
    uint32_t q_len; // data length
    uint32_t q_alloc:31, q_free:1; // allocated data length, ownership of data flag

    // Vector of bitmaps.
    uint32_t n_samples, n_samples_bitmap, n_variants;
    std::vector< std::shared_ptr<djn_bitmap_model_t> > bitmaps;

    // Supportive array for computing allele counts to determine the presence
    // of missing values and/or end-of-vector symbols (in Bcf-encodings).
    uint32_t hist_alts[256];
};

}

#endif