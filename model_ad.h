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
#ifndef MODEL_AD_H_
#define MODEL_AD_H_

#include <vector>

#include <zstd.h>
#include <zstd_errors.h>

#include "frequency_model.h"
#include "range_coder.h"
#include "vcf_reader.h"

namespace djinn {

static int NearestPowerOfTwo(int n) {
    int v = n; 

    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++; // next power of 2

    int x = v >> 1; // previous power of 2

    return (v - n) > (n - x) ? x : v;
}

int LZstdCompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
    int ret = ZSTD_compress(out, out_capacity, in, n_in, c_level);
    return(ret);
}

int LZstdDecompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
    int ret = ZSTD_decompress(out, out_capacity, in, n_in);
    return(ret);
}

class FormatAlelleDepth {
private:
    struct ad_object {
        ad_object() : n(0), prev(0), vals(nullptr){}
        ~ad_object(){ delete[] vals;}

        uint32_t n;
        uint8_t  prev;
        uint8_t* vals;
    };

    static inline int LossWrapper(const int in, const bool lossy = false) {
        return lossy ? NearestPowerOfTwo(in) : in;
    }

public:
    // FormatAlelleDepth();
    FormatAlelleDepth(const uint32_t block_size, const uint32_t n_s, const bool lossy = false) : 
        lossy_encoding(lossy), 
        block_size(block_size),
        n_samples(n_s),
        tot_in(0), tot_out(0),
        buffer1(new uint8_t[block_size+65536]),
        buffer2(new uint8_t[block_size+65536]),
        model8(new pil::FrequencyModel[32768]),
        entries(new ad_object[n_samples]),
        entries2(new ad_object[n_samples])
    {
        for (int i = 0; i < n_samples; ++i) 
            entries[i].vals = new uint8_t[block_size];

        for (int i = 0; i < n_samples; ++i) 
            entries2[i].vals = new uint8_t[block_size];
    }

    ~FormatAlelleDepth() { 
        delete[] entries; 
        delete[] entries2;
        delete[] buffer1;
        delete[] buffer2;
        delete[] model8;
    }

    int AddRecord(const bcf_fmt_t* fmt) {
        if (fmt == NULL) return 0;

        // Todo: fix this
        if(fmt->p_len != 2*n_samples) {
            return 0;
        }

        if (entries[0].n == block_size) {
            Compress();
        }

        // Temp
        assert(fmt->type == BCF_BT_INT8);

        if (entries[0].n == 0) {
            for (int i = 0, j = 0; i + 2 <= fmt->p_len; i += 2, ++j) { // vector length = number of samples
                if (fmt->p[i+1] == 0) {
                    entries[j].vals[0] = 0;
                    entries[j].prev = 0;
                } else {
                    entries[j].vals[0] = LossWrapper(fmt->p[i], lossy_encoding);
                    entries[j].prev    = LossWrapper(fmt->p[i], lossy_encoding);
                }
                ++entries[j].n;
            }

            for (int i = 1, j = 0; i + 2 <= fmt->p_len; i += 2, ++j) {
                if (fmt->p[i+1] == 0) {
                    entries2[j].vals[0] = 0;
                    entries2[j].prev = 0;
                } else {
                    entries2[j].vals[0] = LossWrapper(fmt->p[i], lossy_encoding);
                    entries2[j].prev    = LossWrapper(fmt->p[i], lossy_encoding);
                }
                ++entries[j].n;
            }
        }
        else
        {
            for (int i = 0, j = 0; i + 2 <= fmt->p_len; i += 2, ++j) { // vector length = number of samples
                if (fmt->p[i+1] == 0) {
                    entries[j].vals[entries[j].n] = 0;
                    entries[j].prev = 0;
                } else {
                    entries[j].vals[entries[j].n] = LossWrapper(fmt->p[i], lossy_encoding) - entries[j].prev;
                    entries[j].prev = LossWrapper(fmt->p[i], lossy_encoding);
                    entries[j].vals[entries[j].n] = (entries[j].vals[entries[j].n] << 1) ^ ((entries[j].vals[entries[j].n] >> 7)); // zig-zag
                }
                ++entries[j].n;
            }

            for (int i = 1, j = 0; i + 2 <= fmt->p_len; i += 2, ++j) {
                if (fmt->p[i] == 0) {
                    entries2[j].vals[entries2[j].n] = 0;
                    entries2[j].prev = 0;
                } else {
                    entries2[j].vals[entries2[j].n] = LossWrapper(fmt->p[i], lossy_encoding) - entries2[j].prev;
                    entries2[j].prev = LossWrapper(fmt->p[i], lossy_encoding);
                    entries2[j].vals[entries2[j].n] = (entries2[j].vals[entries2[j].n] << 1) ^ ((entries2[j].vals[entries2[j].n] >> 7)); // zig-zag
                }
                ++entries2[j].n;
            }
        }
    }

    int Compress() {
        uint32_t n_out = 0;
        uint32_t n_out2 = 0;
        
        pil::RangeCoder rc_ad;
        rc_ad.SetOutput(buffer1);
        rc_ad.StartEncode();

        uint64_t ad_total_in = 0, ad_total_out = 0;

        for (int i = 0; i < n_samples; ++i) {
            uint32_t state = 0;
            rc_ad.SetOutput(buffer1);
            rc_ad.StartEncode();

            for (int i = 0; i < 32768; ++i) model8[i].Initiate(8,8);

            ad_total_in += entries[i].n;
            for (int j = 0; j < entries[i].n; ++j) {
                if (entries[i].vals[j] < 8) {
                    model8[state].EncodeSymbol(&rc_ad, entries[i].vals[j]);

                    ++n_out;
                    state <<= 3;
                    state |= entries[i].vals[j];
                    state &= 32767;
                } else {
                    buffer2[n_out2++] = entries[i].vals[j] - 8;
                    state <<= 3;
                    state &= 32767;
                }
                assert(state < 32768);
                assert(rc_ad.OutSize() < block_size + 65536);
            }
            rc_ad.FinishEncode();
            int ret = LZstdCompress((uint8_t*)buffer2, n_out2, buffer1, block_size+65536);
            ad_total_out += ret + rc_ad.OutSize();
            n_out2 = 0;
            n_out = 0;
            entries[i].n = 0;
        }
        std::cerr << "AD1 total=" << ad_total_in << "->" << ad_total_out << " (" << (double)ad_total_in/ad_total_out << "-fold)" << std::endl;
        tot_in += ad_total_in;
        tot_out += ad_total_out;

        ad_total_in = 0, ad_total_out = 0;
        n_out2 = 0, n_out = 0;
        for (int i = 0; i < n_samples; ++i) {
            uint32_t state = 0;
            rc_ad.SetOutput(buffer1);
            rc_ad.StartEncode();

            for (int i = 0; i < 32768; ++i) model8[i].Initiate(8,8);

            ad_total_in += entries2[i].n;
            for (int j = 0; j < entries2[i].n; ++j) {
                if (entries2[i].vals[j] < 8) {
                    model8[state].EncodeSymbol(&rc_ad, entries2[i].vals[j]);
                    
                    ++n_out;
                    state <<= 3;
                    state |= entries2[i].vals[j];
                    state &= 32767;
                } else {
                    buffer2[n_out2++] = entries2[i].vals[j] - 8;
                    state <<= 3;
                    state &= 32767;
                }
                assert(state < 32768);
                assert(rc_ad.OutSize() < block_size + 65536);
            }
            rc_ad.FinishEncode();
            int ret = LZstdCompress((uint8_t*)buffer2, n_out2, buffer1, block_size+65536);
            ad_total_out += ret + rc_ad.OutSize();
            n_out2 = 0;
            n_out = 0;
            entries2[i].n = 0;
        }
        std::cerr << "AD2 total=" << ad_total_in << "->" << ad_total_out << " (" << (double)ad_total_in/ad_total_out << "-fold)" << std::endl;

        tot_in += ad_total_in;
        tot_out += ad_total_out;
    }

public:
    bool lossy_encoding; // should this field be encoded using lossy transformations?
    uint32_t block_size;
    uint32_t n_samples;
    uint64_t tot_in, tot_out;
    uint8_t  *buffer1, *buffer2;
    pil::FrequencyModel* model8;
    ad_object* entries;
    ad_object* entries2;
};

}

#endif