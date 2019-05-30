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
    int Encode(uint64_t* wah, uint32_t len) { // input WAH-encoded data
        uint64_t wah_ref = wah[0];
        uint64_t wah_run = 1;
        uint64_t observed_alts = 0;
        for (int i = 1; i < len; ++i) {
            if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint64_t>::max()) || (wah_ref != wah[i])) {
                if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint64_t>::max()) || wah_run == 1) {
                    // std::cerr << "Dirty: " << std::bitset<64>(wah_ref) << " " << __builtin_popcountll(wah_ref) << "->" << observed_alts + __builtin_popcountll(wah_ref) << std::endl;
                    mtype->EncodeSymbol(0);
                    observed_alts += __builtin_popcountll(wah_ref);
                    
                    // uint8_t wah_bin = 0;
                    
                    for (int i = 0; i < 8; ++i) {
                        // wah_bin |= ((wah_ref & 255) != 0) << i;
                        dirty_wah->EncodeSymbol(wah_ref & 255);
                        wah_ref >>= 8;
                    }
                    // dirty_partition->EncodeSymbol(wah_bin);

                } else {
                    // std::cerr << "Run: " << wah_run << "|" << (wah_ref&1) << " " << (wah_run * 64 * (wah_ref&1)) << "->" << observed_alts + (wah_run * 64 * (wah_ref&1)) << std::endl;
                    observed_alts += wah_run * 64 * (wah_ref&1);
                    mtype->EncodeSymbol(1);
                    EncodeRLE(wah_ref, wah_run);
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
                observed_alts += __builtin_popcountll(wah_ref);
                
                // uint8_t wah_bin = 0;
                    
                for (int i = 0; i < 8; ++i) {
                    // wah_bin |= ((wah_ref & 255) != 0) << i;
                    dirty_wah->EncodeSymbol(wah_ref & 255);
                    wah_ref >>= 8;
                }
                // dirty_partition->EncodeSymbol(wah_bin);
                
            } else {
                // std::cerr << "F Run: " << wah_run << "|" << (wah_ref&1) << " " << (wah_run * 64 * (wah_ref&1)) << "->" << observed_alts + (wah_run * 64 * (wah_ref&1)) << std::endl;
                observed_alts += wah_run * 64 * (wah_ref&1);
                mtype->EncodeSymbol(1);
                EncodeRLE(wah_ref, wah_run);
            }
        }
    }

    int EncodeRLE(uint64_t ref, uint32_t len) {
        mref->EncodeSymbol(ref&1);
        uint32_t log_length = round_log2(len);
        // if (log_length < 2) mlog_rle->EncodeSymbol(0);
        // else if (log_length < 8) mlog_rle->EncodeSymbol(1);
        // else if (log_length < 16) mlog_rle->EncodeSymbol(2);
        // else mlog_rle->EncodeSymbol(3);

        mlog_rle->model_context  = (ref & 1) << 4;
        mlog_rle->model_context |= log_length;
        // std::cerr << std::bitset<32>(mlog_rle->model_context) << std::endl;
        mlog_rle->EncodeSymbolNoUpdate(log_length);

        // uint32_t max_value_prefix = 1u << (log_length);
        // int32_t  add = max_value_prefix - len;

        if (log_length < 2) {
            // std::cerr << "single=" << n_run << "," << add << std::endl;
        }
        else 
        if (log_length <= 8) {
            // std::cerr << "length=" << wah_run << "->" << log_length << std::endl;
            assert(len < 256);
            mrle->model_context = (ref & 1);
            mrle->model_context <<= 4;
            mrle->model_context |= log_length;
            mrle->model_context &= mrle->model_ctx_mask;
            mrle->EncodeSymbolNoUpdate(len & 255);

        } else if (log_length <= 16) {
            // std::cerr << "Log length= " << log_length << " for " << len << std::endl;
            assert(len < 65536);
            // std::cerr << "length=" << wah_run << "->" << log_length << std::endl;
            // mrle2_1->model_context = 0;
            // mrle2_1->model_context <<= 1;
            mrle2_1->model_context = (ref & 1);
            mrle2_1->model_context <<= 4;
            mrle2_1->model_context |= log_length;
            mrle2_1->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            // mrle2_2->model_context = 0;
            // mrle2_2->model_context <<= 1;
            mrle2_2->model_context = (ref & 1);
            mrle2_2->model_context <<= 4;
            mrle2_2->model_context |= log_length;
            mrle2_2->EncodeSymbolNoUpdate(len & 255);
        } else {
            // std::cerr << "here=" << n_run << std::endl;
            // mrle4_1->model_context = 0;
            // mrle4_1->model_context <<= 1;
            mrle4_1->model_context = (ref & 1);
            mrle4_1->model_context <<= 4;
            mrle4_1->model_context |= log_length;
            mrle4_1->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            // mrle4_2->model_context = 0;
            // mrle4_2->model_context <<= 1;
            mrle4_2->model_context = (ref & 1);
            mrle4_2->model_context <<= 4;
            mrle4_2->model_context |= log_length;
            mrle4_2->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            // mrle4_3->model_context = 0;
            // mrle4_3->model_context <<= 1;
            mrle4_3->model_context = (ref & 1);
            mrle4_3->model_context <<= 4;
            mrle4_3->model_context |= log_length;
            mrle4_3->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            // mrle4_4->model_context = 0;
            // mrle4_4->model_context <<= 1;
            mrle4_4->model_context = (ref & 1);
            mrle4_4->model_context <<= 4;
            mrle4_4->model_context |= log_length;
            mrle4_4->EncodeSymbolNoUpdate(len & 255);
        }

        return 1;
    }

    int Compress(uint8_t* out); // compress ctx_block to buffer
    int Compress(djinn_ctx_block_t*& block); // compress ctx_block to given block

    std::shared_ptr<GeneralModel> mref, mlog_rle, mrle, mrle2_1, mrle2_2, mrle4_1, mrle4_2, mrle4_3, mrle4_4;
    std::shared_ptr<GeneralModel> dirty_wah, mtype;
};

}