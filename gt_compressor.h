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
#ifndef GT_COMPRESSOR_H_
#define GT_COMPRESSOR_H_

#include "vcf_reader.h"
#include "pbwt.h"

struct Buffer {
    Buffer() noexcept : len(0), cap(0), data(nullptr){}
    Buffer(const size_t size) noexcept : len(0), cap(size), data(new uint8_t[size]){}
    ~Buffer() { delete[] data; }

    const size_t& size() const { return len; }
    const size_t& capacity() const { return cap; }

    int resize() {
        uint8_t* old = data;
        size_t new_cap = len * 1.2 - len < 65536 ? 65536 : len * 1.2;
        data = new uint8_t[new_cap];
        cap = new_cap;
        memcpy(data,old,len);
        delete[] old;
        
        return 1;
    }

    int resize(const size_t desired_size) {
        if (desired_size < cap) {
            len = desired_size < len ? desired_size : len;
            return 1;
        }

        uint8_t* old = data;
        data = new uint8_t[desired_size];
        cap = desired_size;
        memcpy(data,old,len);
        delete[] old;
        
        return(0);
    }

    void reset() { len = 0; }

    size_t len, cap;
    uint8_t* data;
};

class GenotypeCompressor {
public:
    GenotypeCompressor(int64_t n_s) :
        n_samples(n_s),
        block_size(8192),
        processed_lines(0),
        processed_lines_local(0),
        bytes_in(0), bytes_out(0),
        models(nullptr)
    {
        base_models[0].Construct(n_samples, 2);
        base_models[1].Construct(n_samples, 2);
        base_models[2].Construct(n_samples, 3);
        base_models[3].Construct(n_samples, 3);
        base_models_complex[0].Construct(n_samples, 16);
        base_models_complex[1].Construct(n_samples, 16);

        buf_general[0].resize(10000000);
        buf_general[1].resize(10000000);
        buf_raw.resize(10000000);

        base_models[0].StartEncoding();
        base_models[1].StartEncoding();
        base_models[2].StartEncoding();
        base_models[3].StartEncoding();
        base_models_complex[0].StartEncoding();
        base_models_complex[1].StartEncoding();
    }

    ~GenotypeCompressor() {
        delete[] models;
    }

    inline bool CheckLimit() const {
        return (processed_lines_local == 8196 || 
                base_models[0].range_coder->OutSize() > 9000000 || 
                base_models[1].range_coder->OutSize() > 9000000 || 
                base_models[2].range_coder->OutSize() > 9000000 || 
                base_models[3].range_coder->OutSize() > 9000000 ||
                base_models_complex[0].range_coder->OutSize() > 9000000 ||
                base_models_complex[1].range_coder->OutSize() > 9000000 ||
                buf_raw.len > 9000000);
    }

    //
    int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) {
        if (bcf == NULL) return 0;
        const bcf_fmt_t* fmt = bcf_get_fmt(hdr, bcf, "GT");
        if (fmt == NULL) return 0;

        bytes_in += bcf->d.fmt[0].p_len;
        
        return(Encode2N(bcf));
    }

private:
    // Diploid wrapper.
    int Encode2N(const bcf1_t* bcf) {
        if (bcf->n_allele == 2) return(Encode2N2M(bcf->d.fmt));
        else if (bcf->n_allele < 4) return(Encode2NXM(bcf->d.fmt)); 
        else {
            // std::cerr << "alleles=" << bcf->n_allele << std::endl;
            uint8_t* gts = bcf->d.fmt[0].p;
            for (int i = 0; i < 2*n_samples; ++i) {
                buf_raw.data[buf_raw.len++] = gts[i];
            }

            ++processed_lines_local;
            ++processed_lines;
        }
    }

    // Wrapper for 2N2M
    int Encode2N2M(const bcf_fmt_t* fmt) {
        if (CheckLimit()) {
            Compress();
        }

        // Todo: assert genotypes are set for this variant.
        const uint8_t* gts = fmt[0].p; // data pointer
        
        // Todo: Assert that total number of alleles < 15.
        uint32_t alts[256] = {0};
        for (int i = 0; i < 2*n_samples; ++i) {
            ++alts[BCF_UNPACK_GENOTYPE_GENERAL(gts[i])];
        }

        if (alts[15] == 0) { // No missing values.
            return Encode2N2MC(fmt);

        } else { // Having missing values.
            // std::cerr << "using extended model" << std::endl;
            return Encode2N2MM(fmt);
        }

        return 1;
    }

    // 2N2M complete
    int Encode2N2MC(const bcf_fmt_t* fmt) {
        const uint8_t* gts = fmt[0].p; // data pointer
        
        base_models[0].pbwt->Update(&gts[0], 2);
        base_models[1].pbwt->Update(&gts[1], 2);

        base_models[0].ResetContext();
        for (int j = 0; j < n_samples; ++j) {
            assert(base_models[0].pbwt->prev[j] < 2);
            base_models[0].EncodeSymbol(base_models[0].pbwt->prev[j]);
        }
    
        base_models[1].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            assert(base_models[1].pbwt->prev[i] < 2);
            base_models[1].EncodeSymbol(base_models[1].pbwt->prev[i]);   
        }

        ++processed_lines_local;
        ++processed_lines;
        return 1;
    }

    // 2N2M with missing
    int Encode2N2MM(const bcf_fmt_t* fmt) {
        const uint8_t* gts = fmt[0].p; // data pointer
        
        base_models[3].pbwt->Update(&gts[0], 2);
        base_models[4].pbwt->Update(&gts[1], 2);

        base_models[3].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            assert(base_models[3].pbwt->prev[i] < 3);
            base_models[3].EncodeSymbol(base_models[3].pbwt->prev[i]);
        }

        base_models[4].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            assert(base_models[4].pbwt->prev[i] < 3);
            base_models[4].EncodeSymbol(base_models[4].pbwt->prev[i]);
        }

        ++processed_lines_local;
        ++processed_lines;
        return 1;
    }

    // 2N any M (up to 16)
    int Encode2NXM(const bcf_fmt_t* fmt) {
        if (CheckLimit()) {
            Compress();
        }

        // Todo: assert genotypes are set for this variant.
        const uint8_t* gts = fmt[0].p; // data pointer
        base_models_complex[0].pbwt->UpdateGeneral(&gts[0], 2);
        
        for (int i = 0; i < n_samples; ++i) {
            assert(base_models_complex[0].pbwt->prev[i] < 16);
            base_models_complex[0].EncodeSymbol(base_models_complex[0].pbwt->prev[i]);
        }

        base_models_complex[1].pbwt->UpdateGeneral(&gts[1], 2);

        for (int i = 0; i < n_samples; ++i) {
            assert(base_models_complex[1].pbwt->prev[i] < 16);
            base_models_complex[1].EncodeSymbol(base_models_complex[1].pbwt->prev[i]);
        }

        ++processed_lines_local;
        ++processed_lines;
        return 1;
    }

public:
    int Compress() {
        // flush: temp
        int p1  = base_models[0].FinishEncoding();
        int p2  = base_models[1].FinishEncoding();
        int p1E = base_models[2].FinishEncoding();
        int p2E = base_models[3].FinishEncoding();
        int p2X = base_models_complex[0].FinishEncoding();
        int p2X2 = base_models_complex[1].FinishEncoding();
        int praw = ZstdCompress(buf_raw.data, buf_raw.len,
                                buf_general[0].data, buf_general[0].capacity(),
                                10);

        base_models[0].Reset();
        base_models[0].StartEncoding();
        base_models[1].Reset();
        base_models[1].StartEncoding();
        base_models[2].Reset();
        base_models[2].StartEncoding();
        base_models[3].Reset();
        base_models[3].StartEncoding();
        base_models_complex[0].Reset();
        base_models_complex[0].StartEncoding();
        base_models_complex[1].Reset();
        base_models_complex[1].StartEncoding();

        processed_lines_local = 0;
        std::cerr << "Flushed: " << p1 << "," << p2 << "," << p1E << "," << p2E << "," << p2X << "," << p2X2 << "," << praw << std::endl;
        buf_raw.reset();
        bytes_out += p1+p2+p1E+p2E+p2X+p2X2+praw;
        std::cerr << "[PROGRESS] " << bytes_in << "->" << bytes_out << " (" << (double)bytes_in/bytes_out << "-fold)" << std::endl;
    }

public:
    int64_t n_samples;
    uint32_t block_size;
    uint32_t processed_lines;
    uint32_t processed_lines_local;
    uint64_t bytes_in, bytes_out;
    GeneralPBWTModel base_models[4];
    GeneralPBWTModel base_models_complex[2];
    Buffer buf_general[2];
    Buffer buf_raw;
    GeneralPBWTModel* models;
};

#endif