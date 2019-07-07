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

#include <cstdlib>//malloc
#include <cstring>//memcpy
#include <cassert>//assert

#include <iostream>//debug
#include <bitset>//debug

namespace djinn {

// Map missing to 2, 1->0, 2->1, and EOV -> 3.
constexpr const uint8_t TWK_BCF_GT_UNPACK[65] = 
{2,0,1,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3};

constexpr const char DJN_UNPACK_2M_VCF[4] = {'0', '1', '.', 'x'}; // position 3 is the EOV marker
// const char DJN_UNPACK_VCF[16] = {}

// MISSING -> 14, EOV -> 15, other values as normal
// Note that currently we can only store up to 14 ALT
// alleles in addition to missing values and EOV symbols.
constexpr const uint8_t TWK_BCF_GT_UNPACK_GENERAL[256] = 
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

constexpr const uint8_t TWK_BCF_GT_UNPACK_GENERAL_REV[131] = 
{0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,15,16,17,18,19,20,21,22,
23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,
43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,
63,15,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,
83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,
103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,
118,119,120,121,122,123,124,125,126,127,128,129};

constexpr const uint8_t DJN_MAP_NONE[256] = {
0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
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
251,252,253,254,255
};


#define BCF_UNPACK_GENOTYPE(A) TWK_BCF_GT_UNPACK[(A) >> 1]
#define BCF_UNPACK_GENOTYPE_GENERAL(A) TWK_BCF_GT_UNPACK_GENERAL[(A) >> 1]

//
#define DJN_UN_NONE 0 // up to nothing
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
    djn_variant_dec_t() : m_ewah(0), m_dirty(0), n_ewah(0), n_dirty(0), dirty_type(0), ewah(nullptr), dirty(nullptr) {}
    ~djn_variant_dec_t() {
        delete[] ewah;
        delete[] dirty;
    }

    void Allocate(uint32_t n) {
        delete[] ewah; delete[] dirty;
        m_ewah = n;
        m_dirty = n;
        n_ewah = 0, n_dirty = 0;
        ewah = new djinn_ewah_t*[m_ewah];
        dirty = new uint32_t*[m_dirty];
    }

    int m_ewah, m_dirty;
    int n_ewah, n_dirty;
    int dirty_type; // one of DJN_DIRTY_*
    djinn_ewah_t **ewah;
    uint32_t **dirty;
};

// Variant record
struct djinn_variant_t {
    djinn_variant_t() : ploidy(0), data(nullptr), data_len(0), data_alloc(0), data_free(false), errcode(0), unpacked(0), d(nullptr) {}
    ~djinn_variant_t() {
        if (data_free) delete[] data;
        delete d;
    }

    int ploidy;
    uint8_t* data;
    uint32_t data_len;
    uint32_t data_alloc: 31, data_free: 1;
    int errcode;
    int unpacked; // one of DJN_UN_*
    djn_variant_dec_t* d;
};

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

/**
 * Naming convention: djinn_* structures/classes are considered front-end
 * and djn_*_t* are considered back-end.
 * 
 * GeneralModel
 * djinn_model               base model for context modelling and EWAH
 * 
 * djn_block_t
 * djn_ctx_store_t
 * djn_ctx_store_tt          st
 * djn_ctx_model_t           context models for storing haplotypes: models for 2MC and NM
 * djn_ctx_model_container_t model container for each (#samples,ploidy)-tuple context model
 * djinn_ctx_model           root interface extending djinn_model
 */

class djinn_model {
public:
    djinn_model() : use_pbwt(true), init(true), unused(0), n_variants(0) {}
    virtual ~djinn_model() {}

    virtual int EncodeBcf(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) =0;
    virtual int Encode(uint8_t* data, size_t len_data, int ploidy, uint8_t alt_alleles) =0;

    virtual void StartEncoding(bool use_pbwt, bool reset = false) =0;
    virtual size_t FinishEncoding() =0;
    virtual int StartDecoding() =0;

    // Decoding requires:
    // One buffer for decoding into EWAH values.
    // One buffer of size no smaller than ploidy*samples for storing the return vector of alleles.
    virtual int DecodeNext(uint8_t* ewah_data, uint32_t& ret_ewah, uint8_t* ret_buffer, uint32_t& ret_len) =0;
    virtual int DecodeNext(djinn_variant_t*& variant) =0;
    virtual int DecodeNextRaw(uint8_t* data, uint32_t& len) =0;
    virtual int DecodeNextRaw(djinn_variant_t*& variant) =0;

    // Read write
    virtual int Serialize(uint8_t* dst) const =0;
    virtual int Serialize(std::ostream& stream) const =0;
    virtual int Deserialize(uint8_t* dst) =0;
    virtual int Deserialize(std::istream& stream) =0;

    // Todo:
    // virtual int Merge(djinn_model* b1, djinn_model* b2);

public:
    uint8_t use_pbwt: 1, 
            init: 1, 
            unused: 6;
    uint32_t n_variants; // number of encoded variants

    // Supportive array for computing allele counts to determine the presence
    // of missing values and/or end-of-vector symbols (in Bcf-encodings).
    uint32_t hist_alts[256];
};

}

#endif