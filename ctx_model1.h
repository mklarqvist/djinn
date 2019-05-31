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

template<typename T>
uint32_t round_log2(T x)
{
	uint32_t r = 0;

	for (; x; ++r)
		x >>= 1;

	return r;
}

class djinn_ctx_model {
public:
    djinn_ctx_model() : own_buffer(false), buffer(nullptr), n_samples(0), n_wah(0), wah_bitmaps(nullptr) {}
    djinn_ctx_model(uint64_t n_s) : own_buffer(false), buffer(nullptr), n_samples(n_s), n_wah(std::ceil((float)2*n_samples / 64)), wah_bitmaps(new uint64_t[n_wah])
    {

    }
    
    ~djinn_ctx_model() { 
        if(own_buffer) delete[] buffer; 
        delete[] wah_bitmaps;
    }

    void SetSamples(int64_t n_s) { n_samples = n_s; } 

    int Encode(uint8_t* data, uint32_t len) {
        if (data == nullptr) return -1;
        if (n_samples == 0) return -2;
        if (n_wah == 0) {
            n_wah = std::ceil((float)n_samples / 64);
            wah_bitmaps = new uint64_t[n_wah];
        }
        
        memset(wah_bitmaps, 0, n_wah*sizeof(uint64_t));

        for (int i = 0; i < n_samples; ++i) {
            if (data[i]) {
                wah_bitmaps[i / 64] |= 1L << (i % 64);
            }
        }

        return EncodeWah(wah_bitmaps, n_wah);
    }

    int EncodeWah(uint64_t* wah, uint32_t len) { // input WAH-encoded data
        if (wah == nullptr) return -1;

        uint64_t wah_ref = wah[0];
        uint64_t wah_run = 1;

        for (int i = 1; i < len; ++i) {
            if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint64_t>::max()) || (wah_ref != wah[i])) {
                if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint64_t>::max()) || wah_run == 1) {
                    // std::cerr << "Dirty: " << std::bitset<64>(wah_ref) << " " << __builtin_popcountll(wah_ref) << "->" << observed_alts + __builtin_popcountll(wah_ref) << std::endl;
                    mtype->EncodeSymbol(0);
                    // observed_alts += __builtin_popcountll(wah_ref);
                    
                    // uint8_t wah_bin = 0;
                    
                    for (int i = 0; i < 8; ++i) {
                        // wah_bin |= ((wah_ref & 255) != 0) << i;
                        dirty_wah->EncodeSymbol(wah_ref & 255);
                        wah_ref >>= 8;
                    }
                    // dirty_partition->EncodeSymbol(wah_bin);

                } else {
                    // std::cerr << "Run: " << wah_run << "|" << (wah_ref&1) << " " << (wah_run * 64 * (wah_ref&1)) << "->" << observed_alts + (wah_run * 64 * (wah_ref&1)) << std::endl;
                    // observed_alts += wah_run * 64 * (wah_ref&1);
                    mtype->EncodeSymbol(1);
                    EncodeWahRLE(wah_ref, wah_run);
                }
                wah_run = 0;
                wah_ref = wah[i];
            }
            ++wah_run;
        }

        if (wah_run) {
            // std::cerr << "Final=" << wah_run << std::endl;
            if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint64_t>::max()) || wah_run == 1) {
                // std::cerr << "F Dirty: " << std::bitset<64>(wah_ref) << " " << __builtin_popcountll(wah_ref) << "->" << observed_alts + __builtin_popcountll(wah_ref) << std::endl;
                mtype->EncodeSymbol(0);
                // observed_alts += __builtin_popcountll(wah_ref);
                
                // uint8_t wah_bin = 0;
                    
                for (int i = 0; i < 8; ++i) {
                    // wah_bin |= ((wah_ref & 255) != 0) << i;
                    dirty_wah->EncodeSymbol(wah_ref & 255);
                    wah_ref >>= 8;
                }
                // dirty_partition->EncodeSymbol(wah_bin);
                
            } else {
                // std::cerr << "F Run: " << wah_run << "|" << (wah_ref&1) << " " << (wah_run * 64 * (wah_ref&1)) << "->" << observed_alts + (wah_run * 64 * (wah_ref&1)) << std::endl;
                // observed_alts += wah_run * 64 * (wah_ref&1);
                mtype->EncodeSymbol(1);
                EncodeWahRLE(wah_ref, wah_run);
            }
        }

        return 1;
    }

    int EncodeWahRLE(uint64_t ref, uint32_t len) {
        mref->EncodeSymbol(ref&1);
        uint32_t log_length = round_log2(len);
        mlog_rle->EncodeSymbol(log_length);

        if (log_length < 2) {
            std::cerr << "single=" << len << "," << (ref&1) << std::endl;
        }
        else if (log_length <= 8) {
            assert(len < 256);

            // uint32_t max_value_for_this_prefix = 1u << (log_length - 1);
            // assert(len - max_value_for_this_prefix < 256);

            mrle->model_context <<= 1;
            mrle->model_context |= (ref & 1);
            mrle->model_context <<= 4;
            mrle->model_context |= log_length;
            mrle->model_context &= mrle->model_ctx_mask;
            mrle->EncodeSymbolNoUpdate(len & 255);

        } else if (log_length <= 16) {
            assert(len < 65536);

            // uint32_t log_length2 = round_log2(len >> 8);
            // std::cerr << "second=" << (len>>8) << "->" << log_length2 << std::endl;
            // mlog_rle->EncodeSymbol(log_length2);

            // uint32_t max_value_for_this_prefix = 1u << (log_length - 8 - 1);
            // uint32_t max_value_for_this_prefix2 = 1u << (8 - 1);
            // std::cerr << len << "=" << log_length << "-> " << max_value_for_this_prefix << " and " << max_value_for_this_prefix2 << std::endl;
            // std::cerr << (len & 255) - max_value_for_this_prefix2 << " and " << (len >> 8) - max_value_for_this_prefix << std::endl;

            mrle2_1->model_context <<= 1;
            mrle2_1->model_context |= (ref & 1);
            mrle2_1->model_context <<= 4;
            mrle2_1->model_context |= log_length;
            mrle2_1->model_context &= mrle2_1->model_ctx_mask;
            mrle2_1->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            mrle2_2->model_context <<= 1;
            mrle2_2->model_context |= (ref & 1);
            mrle2_2->model_context <<= 4;
            mrle2_2->model_context |= log_length;
            mrle2_2->model_context &= mrle2_2->model_ctx_mask;
            mrle2_2->EncodeSymbolNoUpdate(len & 255);
        } else {
            mrle4_1->model_context <<= 1;
            mrle4_1->model_context |= (ref & 1);
            mrle4_1->model_context <<= 4;
            mrle4_1->model_context |= log_length;
            mrle4_1->model_context &= mrle4_1->model_ctx_mask;
            mrle4_1->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            mrle4_2->model_context <<= 1;
            mrle4_2->model_context |= (ref & 1);
            mrle4_2->model_context <<= 4;
            mrle4_2->model_context |= log_length;
            mrle4_2->model_context &= mrle4_2->model_ctx_mask;
            mrle4_2->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            mrle4_3->model_context <<= 1;
            mrle4_3->model_context |= (ref & 1);
            mrle4_3->model_context <<= 4;
            mrle4_3->model_context |= log_length;
            mrle4_3->model_context &= mrle4_3->model_ctx_mask;
            mrle4_3->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            mrle4_4->model_context <<= 1;
            mrle4_4->model_context |= (ref & 1);
            mrle4_4->model_context <<= 4;
            mrle4_4->model_context |= log_length;
            mrle4_4->model_context &= mrle4_4->model_ctx_mask;
            mrle4_4->EncodeSymbolNoUpdate(len & 255);
        }

        return 1;
    }

    int Compress(uint8_t* out); // compress ctx_block to buffer
    int Compress(djinn_ctx_block_t*& block); // compress ctx_block to given block

    void InitiateEncoding() {
        if (buffer == nullptr) {
             buffer = new uint8_t[20000000];
             own_buffer = true;
        }

        Initiate();
    }

    void Initiate() {
        if (range_coder.get() == nullptr) range_coder = std::make_shared<RangeCoder>();
        if (mref.get() == nullptr) mref = std::make_shared<GeneralModel>(2, 4096, range_coder);
        if (mlog_rle.get() == nullptr) mlog_rle = std::make_shared<GeneralModel>(32, 4096, range_coder); // 2^(4+1) order-3
        if (mrle.get() == nullptr) mrle = std::make_shared<GeneralModel>(256, 1024, 24, 1, range_coder); // 2^10 order-2
        if (mrle2_1.get() == nullptr) mrle2_1 = std::make_shared<GeneralModel>(256, 1024, range_coder); // 2^10 order-2
        if (mrle2_2.get() == nullptr) mrle2_2 = std::make_shared<GeneralModel>(256, 1024, range_coder);
        if (mrle4_1.get() == nullptr) mrle4_1 = std::make_shared<GeneralModel>(256, 1024, range_coder);
        if (mrle4_2.get() == nullptr) mrle4_2 = std::make_shared<GeneralModel>(256, 1024, range_coder);
        if (mrle4_3.get() == nullptr) mrle4_3 = std::make_shared<GeneralModel>(256, 1024, range_coder);
        if (mrle4_4.get() == nullptr) mrle4_4 = std::make_shared<GeneralModel>(256, 1024, range_coder);
        if (dirty_wah.get() == nullptr) dirty_wah = std::make_shared<GeneralModel>(256, 65536, 24, 1, range_coder);
        if (mtype.get() == nullptr) mtype = std::make_shared<GeneralModel>(2, 512, range_coder);
        if (mref.get() == nullptr) mref = std::make_shared<GeneralModel>(2, 1024, range_coder);
    }

    void StartEncoding(bool reset = false) {
        InitiateEncoding();

        // If resetting the model.
        if (reset) {
            mref->Reset();
            mlog_rle->Reset();
            mrle->Reset();
            mrle2_1->Reset(); mrle2_2->Reset();
            mrle4_1->Reset(); mrle4_2->Reset(); mrle4_3->Reset(); mrle4_4->Reset();
            dirty_wah->Reset();
            mtype->Reset();
            mref->Reset();
        }

        range_coder->SetOutput(buffer);
        range_coder->StartEncode();
    }

    size_t FinishEncoding() {
        if (range_coder.get() == nullptr) return -1;
        range_coder->FinishEncode();
        return range_coder->OutSize();
    }

    void StartDecoding(uint8_t* data) {
        if (own_buffer) delete[] buffer;
        buffer = data;
        own_buffer = false;

        range_coder->SetInput(buffer);
        range_coder->StartDecode();
    }

public:
    bool use_pbwt, own_buffer;
    int64_t n_samples;
    uint64_t n_wah;
    uint64_t* wah_bitmaps;

    std::shared_ptr<RangeCoder> range_coder;
    std::shared_ptr<PBWT> pbwt_haploid, pbwt_diploid;
    std::shared_ptr<GeneralModel> mref, mlog_rle, mrle, mrle2_1, mrle2_2, mrle4_1, mrle4_2, mrle4_3, mrle4_4;
    std::shared_ptr<GeneralModel> dirty_wah, mtype;
    uint8_t* buffer;
};

}