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

#include <limits>

#include "djinn.h"
#include "pbwt.h"

namespace djinn {

static 
uint32_t round_log2(uint32_t x) {
	uint32_t r = 0;
	for (/**/; x; ++r) x >>= 1;
	return r;
}

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
    void Initiate2mc();
    void InitiateNm();

    void reset();
    int StartEncoding(bool use_pbwt, bool reset = false);
    size_t FinishEncoding();
    int StartDecoding(uint8_t* data, bool reset = false);

public:
    uint8_t use_pbwt: 1, init: 1, unused: 6;
    PBWT pbwt;
    std::shared_ptr<RangeCoder> range_coder; // shared range coder
    std::shared_ptr<GeneralModel> mref, mlog_rle, mrle, mrle2_1, mrle2_2, mrle4_1, mrle4_2, mrle4_3, mrle4_4; // rle models
    std::shared_ptr<GeneralModel> dirty_wah; // dirty word model
    std::shared_ptr<GeneralModel> mtype; // mtype model
    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1;
    uint32_t n_variants;
};

/*======   GT context model   ======*/

class djinn_ctx_model {
private:
    static constexpr uint64_t ref_bits[16] = 
        {0,                     72340172838076672ULL,  144680345676153344ULL,  217020518514230016ULL, 
         289360691352306688ULL, 361700864190383360ULL, 434041037028460032ULL,  506381209866536704ULL, 
         578721382704613376ULL, 651061555542690048ULL, 723401728380766720ULL,  795741901218843392ULL, 
         868082074056920064ULL, 940422246894996736ULL, 1012762419733073408ULL, 1085102592571150080ULL};

public:
    djinn_ctx_model();
    djinn_ctx_model(uint64_t n_s);
    ~djinn_ctx_model();

    void SetSamples(int64_t n_s);
    int EncodeBcf(uint8_t* data, uint8_t alt_alleles);

    void StartEncoding(bool use_pbwt, bool reset = false);
    size_t FinishEncoding();
    int StartDecoding(djinn_block_t* block, bool reset = false);
    int GetBlockReference(djinn_block_t*& block);
    // Todo
    int GetBlock(djinn_block_t*& block);
    
// Internal functions
public:
    int Encode2mc(uint8_t* data, uint32_t len);
    int EncodeNm(uint8_t* data, uint32_t len);
    int EncodeWah(uint64_t* wah, uint32_t len);
    int EncodeWahNm(uint64_t* wah, uint32_t len);
    int EncodeWahRLE(uint64_t ref, uint32_t len, djinn_ctx_model_t* model);
    int EncodeWahRLE_nm(uint64_t ref, uint32_t len, djinn_ctx_model_t* model);
    
    int DecodeNext(uint8_t* data);
    // int DecodeNextRaw(uint8_t* data);

    int DecodeRaw(uint8_t* data);
    int DecodeRaw_nm(uint8_t* data);
    int DecodeWahRLE(uint64_t& ref, uint32_t& len, djinn_ctx_model_t* model);
    int DecodeWahRLE_nm(uint64_t& ref, uint32_t& len, djinn_ctx_model_t* model);

public:
    int64_t n_samples;
    
    uint64_t n_wah;
    int64_t n_samples_wah;
    uint64_t* wah_bitmaps;

    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1;
    uint32_t n_variants;

    uint32_t hist_alts[256];
    
    std::shared_ptr<RangeCoder> range_coder;
    std::shared_ptr<GeneralModel> marchetype; // 0 for 2MC, 1 for 2M, 2 else
    djinn_ctx_model_t model_2mc;
    djinn_ctx_model_t model_2m;
    djinn_ctx_model_t model_nm;
};

}