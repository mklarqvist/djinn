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

// temp
#include <bitset>

class MixedCompressor {
public:

    int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) {

        // Compression scheme:
        // If run-length > X then use RLE.
        // Otherwise use 64-bit bitmap.
        // Compress with Zstd or LZ4.
        if (bcf == NULL) return 0;
        const bcf_fmt_t* fmt = bcf_get_fmt(hdr, bcf, "GT");
        if (fmt == NULL) return 0;

        uint8_t* gts = fmt->p;
        for (int i = 0; i < 2504; ++i) {
            
        }
    }

    uint32_t lbuf;
    uint8_t* buf;
};

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

static uint8_t temp_unpack[3] = {2, 0, 1};

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

    virtual ~GenotypeCompressor() { delete[] models; }
    
    // Encode data using literals.
    // virtual int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) =0;
    virtual int Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles) =0;
    virtual int Encode2N2M(uint8_t* data, const int32_t n_data) =0;
    virtual int Encode2N2MC(uint8_t* data, const int32_t n_data) =0;
    virtual int Encode2N2MM(uint8_t* data, const int32_t n_data) =0;
    virtual int Encode2NXM(uint8_t* data, const int32_t n_data) =0;

    // Encode data using htslib.
    virtual 
    int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) {
        if (bcf == NULL) return 0;
        if (hdr == NULL) return 0;

        const bcf_fmt_t* fmt = bcf_get_fmt(hdr, bcf, "GT");
        if (fmt == NULL) return 0;
        if (fmt->p_len / n_samples != 2) {
            std::cerr << "input is not divisible by 2" << std::endl;
            return 0;
        }

        bytes_in += fmt->p_len;

        // Todo: extend beyond 2N
        return(Encode2N(fmt->p, fmt->p_len, bcf->n_allele));
    }
    virtual int Encode2N(bcf1_t* bcf, const bcf_hdr_t* hdr) =0;

    //
    int32_t RemapGenotypeEOV(uint8_t* data, const uint32_t len) {
        if (data == NULL) {
            std::cerr << "data==NULL" << std::endl;
            return -1;
        }

        memset(alts, 0, 256*sizeof(uint32_t));
        
        uint8_t max_val = 0;
        for (int i = 0; i < 2*n_samples; ++i) {
            uint8_t val = BCF_UNPACK_GENOTYPE_GENERAL(data[i]);
            ++alts[val];
            val = (val == 64 ? 0 : val);
            max_val = val > max_val ? val : max_val;
        }
        
        int32_t n_replaced = 0;
        if (alts[64]) {
            std::cerr << "replacing 64 with " << (int32_t)max_val << std::endl;
            for (int i = 0; i < len; ++i) {
                if (BCF_UNPACK_GENOTYPE_GENERAL(data[i]) == 64) {
                    data[i] = max_val;
                    ++n_replaced;
                }
            }
            std::cerr << "replaced=" << n_replaced << std::endl;
        }

        return n_replaced;
    }
    
public:
    int64_t n_samples;
    uint32_t block_size;
    uint32_t processed_lines;
    uint32_t processed_lines_local;
    
    uint64_t bytes_in, bytes_out;
    
    Buffer buf_general[2];
    Buffer buf_raw;

    

    uint32_t alts[256];

    GeneralPBWTModel base_models[4];
    GeneralPBWTModel base_models_complex[2];
    GeneralPBWTModel* models;
};

class GenotypeCompressorModelling : public GenotypeCompressor {
public:
    GenotypeCompressorModelling(int64_t n_s) : GenotypeCompressor(n_s)
    {
        nonsense[0].Construct(1, 2);
        nonsense[1].Construct(1, 2);
        nonsense[0].StartEncoding();
        nonsense[1].StartEncoding();
    }

    ~GenotypeCompressorModelling() {
        
    }

    bool CheckLimit() const {
        return (processed_lines_local == block_size || 
                base_models[0].range_coder->OutSize() > 9000000 || 
                base_models[1].range_coder->OutSize() > 9000000 || 
                base_models[2].range_coder->OutSize() > 9000000 || 
                base_models[3].range_coder->OutSize() > 9000000 ||
                base_models_complex[0].range_coder->OutSize() > 9000000 ||
                base_models_complex[1].range_coder->OutSize() > 9000000 ||
                buf_raw.len > 9000000);
    }

    //
    int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) override {
        if (bcf == NULL) return 0;
        if (hdr == NULL) return 0;
        
        const bcf_fmt_t* fmt = bcf_get_fmt(hdr, bcf, "GT");
        if (fmt == NULL) return 0;
        if (fmt->p_len / n_samples != 2) {
            std::cerr << "input is not divisible by 2" << std::endl;
            return 0;
        }

        bytes_in += bcf->d.fmt->p_len;

        // Todo: extend beyond 2N
        return(Encode2N(bcf,hdr));
    }

private:
    // Diploid wrapper.
    int Encode2N(bcf1_t* bcf, const bcf_hdr_t* hdr) override {
        if (bcf == NULL) return 0;
        if (hdr == NULL) return 0;
        
        const bcf_fmt_t* fmt = bcf_get_fmt(hdr, bcf, "GT");
        if (fmt == NULL) return 0;
        if (fmt->p_len / n_samples != 2) {
            std::cerr << "input is not divisible by 2" << std::endl;
            return 0;
        }

        if (bcf->n_allele == 2) return(Encode2N2M(fmt->p, fmt->p_len)); // biallelic
        else if (bcf->n_allele < 4) {

            int replaced = RemapGenotypeEOV(fmt->p, fmt->p_len);
            
            // std::cerr << "second path taken for " << bcf->n_allele << std::endl; 
            return(Encode2NXM(fmt->p, fmt->p_len));  // #alleles < 4
        }
        else {
            // std::cerr << "alleles=" << bcf->n_allele << std::endl;
            uint8_t* gts = fmt->p;
            for (int i = 0; i < 2*n_samples; ++i) {
                buf_raw.data[buf_raw.len++] = gts[i];
            }

            ++processed_lines_local;
            ++processed_lines;
        }
        return 1;
    }

    int Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles) {
        if (n_alleles == 2) return(Encode2N2M(data, n_data)); // biallelic
        else if (n_alleles < 4) return(Encode2NXM(data, n_data));  // #alleles < 4
        else {
            // std::cerr << "alleles=" << n_alleles << std::endl;
            const uint8_t* gts = data;
            for (int i = 0; i < 2*n_samples; ++i) {
                buf_raw.data[buf_raw.len++] = gts[i];
            }

            ++processed_lines_local;
            ++processed_lines;
        }
        return 1;
    }

    // Wrapper for 2N2M
    int Encode2N2M(uint8_t* data, const int32_t n_data) override {
        if (CheckLimit()) {
            Compress();
        }

        // Todo: assert genotypes are set for this variant.
        int replaced = RemapGenotypeEOV(data, n_data);
        if (replaced) {
            std::cerr << "replaced: divert to missing with " << replaced << std::endl;
            return Encode2NXM(data ,n_data);
        }

        if (alts[130] == 0) { // No missing values.
            return Encode2N2MC(data, n_data);
        } else { // Having missing values.
            // std::cerr << "using extended model" << std::endl;
            return Encode2N2MM(data, n_data);
        }

        return 1;
    }

    // 2N2M complete
    int Encode2N2MC(uint8_t* data, const int32_t n_data) {
        if (CheckLimit()) {
            Compress();
        }
        
        const uint8_t* gts = data; // data pointer
        base_models[0].pbwt->Update(&gts[0], 2);
        base_models[1].pbwt->Update(&gts[1], 2);

        // EncodeRLEBitmap(0);
        // EncodeRLEBitmap(1);

#if 1
        int n_steps = 16;
        uint32_t step_size = std::ceil((float)n_samples / n_steps);
        uint64_t bins1 = 0;
        for (int i = 0; i < n_samples; ++i) {
            if (base_models[0].pbwt->prev[i]) {
                bins1 |= (1 << (i/step_size));
            }
        }

        base_models[0].ResetContext();
        nonsense[0].ResetContext();
        uint32_t offset = 0, offset_end = 0;
        for (int i = 0; i < n_steps; ++i) {
            if (bins1 & (1 << i)) {
                nonsense[0].EncodeSymbol(1);
                offset_end = offset + step_size < n_samples ? offset + step_size : n_samples;
                for (int j = offset; j < offset_end; ++j) {
                    base_models[0].EncodeSymbol(base_models[0].pbwt->prev[j]);
                }
            } else nonsense[0].EncodeSymbol(0);
            offset += step_size;
        }

        uint64_t bins2 = 0;
        for (int i = 0; i < n_samples; ++i) {
            if (base_models[1].pbwt->prev[i]) {
                bins2 |= (1 << (i/step_size));
            }
        }
   
        base_models[1].ResetContext();
        nonsense[1].ResetContext();
        offset = 0, offset_end = 0;
        for (int i = 0; i < n_steps; ++i) {
            if (bins2 & (1 << i)) {
                nonsense[1].EncodeSymbol(1);
                offset_end = offset + step_size < n_samples ? offset + step_size : n_samples;
                for (int j = offset; j < offset_end; ++j) {
                    base_models[1].EncodeSymbol(base_models[1].pbwt->prev[j]);
                }
            } else nonsense[1].EncodeSymbol(0);
            offset += step_size;
        }
#endif

#if 0
        base_models[0].ResetContext();
        for (int j = 0; j < n_samples; ++j) {
            // assert(base_models[0].pbwt->prev[j] < 2);
            base_models[0].EncodeSymbol(base_models[0].pbwt->prev[j]);
        }
    
        base_models[1].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            // assert(base_models[1].pbwt->prev[i] < 2);
            base_models[1].EncodeSymbol(base_models[1].pbwt->prev[i]);
        }
#endif

        ++processed_lines_local;
        ++processed_lines;
        return 1;
    }

    // 2N2M with missing
    int Encode2N2MM(uint8_t* data, const int32_t n_data) {
        const uint8_t* gts = data; // data pointer
        
        base_models[2].pbwt->Update(&gts[0], 2);
        base_models[3].pbwt->Update(&gts[1], 2);

        base_models[2].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            assert(base_models[2].pbwt->prev[i] < 3);
            base_models[2].EncodeSymbol(base_models[2].pbwt->prev[i]);
        }

        base_models[3].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            assert(base_models[3].pbwt->prev[i] < 3);
            base_models[3].EncodeSymbol(base_models[3].pbwt->prev[i]);
        }

        ++processed_lines_local;
        ++processed_lines;
        return 1;
    }

    // 2N any M (up to 16)
    int Encode2NXM(uint8_t* data, const int32_t n_data) {
        if (CheckLimit()) {
            Compress();
        }

        const uint8_t* gts = data; // data pointer

        // Todo: assert genotypes are set for this variant.
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
                                20);

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

        int extra1 = nonsense[0].FinishEncoding();
        int extra2 = nonsense[1].FinishEncoding();
        std::cerr << processed_lines_local << "->" << extra1 << " and " << extra2 << std::endl;
        nonsense[0].Reset();
        nonsense[0].StartEncoding();
        nonsense[1].Reset();
        nonsense[1].StartEncoding();

        processed_lines_local = 0;
        std::cerr << "Flushed: " << p1 << "," << p2 << "," << p1E << "," << p2E << "," << p2X << "," << p2X2 << "," << praw << std::endl;
        buf_raw.reset();
        bytes_out += p1 + p2 + p1E + p2E + p2X + p2X2 + praw + extra1 + extra2;
        // bytes_out += p1 + p2 + p1E + p2E + extra1 + extra2;
        std::cerr << "[PROGRESS] " << bytes_in << "->" << bytes_out << " (" << (double)bytes_in/bytes_out << "-fold)" << std::endl;

        // Debug compression:
        // std::cerr << "[PROGRESS ZSTD] " << bytes_in1 << "->" << bytes_out_zstd1 << " (" << (double)bytes_in1/bytes_out_zstd1 << "-fold)" << std::endl;
        // std::cerr << "[PROGRESS LZ4] " << bytes_in1 << "->" << bytes_out_lz4 << " (" << (double)bytes_in1/bytes_out_lz4 << "-fold)" << std::endl;

        return 1;
    }

public:
    GeneralPBWTModel nonsense[2];
};

class GenotypeCompressorRLEBitmap : public GenotypeCompressor {
public:
    GenotypeCompressorRLEBitmap(int64_t n_s) : GenotypeCompressor(n_s),
        bytes_in1(0), bytes_out_zstd1(0), bytes_out_lz4(0)
    {
        gt_width.resize(n_s, 0);
        buf_wah[0].resize(10000000);
        buf_wah[1].resize(10000000);
    }

    int Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles) override;
    int Encode2N2M(uint8_t* data, const int32_t n_data) override;
    int Encode2N2MC(uint8_t* data, const int32_t n_data) override;
    int Encode2N2MM(uint8_t* data, const int32_t n_data) override;
    int Encode2NXM(uint8_t* data, const int32_t n_data) override;

    // Encode data using htslib.
    int Encode2N(bcf1_t* bcf, const bcf_hdr_t* hdr) override;

    int Compress() {
        // RLE-bitmap hybrid.
        if (buf_wah[0].len) {
            // Temporary debug before compressing.
            DebugRLEBitmap(buf_wah[0].data, buf_wah[0].len, n_samples);

            int praw = ZstdCompress(buf_wah[0].data, buf_wah[0].len,
                                    buf_general[0].data, buf_general[0].capacity(),
                                    20);

            bytes_out_zstd1 += praw;

            std::cerr << "Zstd->compressed: " << buf_wah[0].len << "->" << praw << std::endl;
            size_t plz4 = Lz4Compress(buf_wah[0].data, buf_wah[0].len,
                                      buf_general[0].data, buf_general[0].capacity(),
                                      9);
            
            // std::cerr << "writing=" << buf_wah[0].len << " and " << plz4 << std::endl; 
            
            // Temp: emit data to cout.
            std::cout.write((char*)&buf_wah[0].len, sizeof(size_t));
            std::cout.write((char*)&plz4, sizeof(size_t));
            std::cout.write((char*)buf_general[0].data, plz4);
            std::cout.flush();                      
            
            std::cerr << "Lz4->compressed: " << buf_wah[0].len << "->" << plz4 << std::endl;
            bytes_out_lz4 += plz4;

            buf_wah[0].len = 0;
        }

        base_models[0].Reset();
        base_models[1].Reset();
        base_models[2].Reset();
        base_models[3].Reset();
        base_models_complex[0].Reset();
        base_models_complex[1].Reset();

        processed_lines_local = 0;
        buf_raw.reset();
        // bytes_out += p1 + p2 + p1E + p2E + p2X + p2X2 + praw + extra1 + extra2;
        // bytes_out += p1 + p2 + p1E + p2E + extra1 + extra2;
        std::cerr << "[PROGRESS] " << bytes_in << "->" << bytes_out << " (" << (double)bytes_in/bytes_out << "-fold)" << std::endl;

        // Debug compression:
        std::cerr << "[PROGRESS ZSTD] " << bytes_in1 << "->" << bytes_out_zstd1 << " (" << (double)bytes_in1/bytes_out_zstd1 << "-fold)" << std::endl;
        std::cerr << "[PROGRESS LZ4] " << bytes_in1 << "->" << bytes_out_lz4 << " (" << (double)bytes_in1/bytes_out_lz4 << "-fold)" << std::endl;

    }

    int EncodeRLEBitmap(const int target) {
        assert(target == 0 || target == 1);

        bytes_in1 += n_samples;

        // RLE word -> 1 bit typing, 1 bit allele, word*8-2 bits for RLE
        const uint8_t* prev = base_models[target].pbwt->prev;
        uint8_t* data = buf_wah[0].data;
        size_t& data_len = buf_wah[0].len;

        // Setup.
        uint8_t ref = prev[0];
        uint32_t rle_len = 1;
        uint32_t start_rle = 0, end_rle = 1;

        // Debug.
        uint32_t rle_cost = 0;
        uint32_t n_objects = 0;
        uint32_t observed_alts = 0;
        uint32_t observed_length = 0;

        for (int i = 1; i < n_samples; ++i) {
            if (prev[i] != ref || rle_len == 16383) {
                end_rle = i;
                ++n_objects;

                // If the run length is < 30 then the word is considered "dirty"
                // and a 32-bit bitmap will be used instead.
                if (end_rle - start_rle < 30) {
                    // Increment position
                    i += 31 - (end_rle - start_rle);
                    // Make sure the end position is within range.
                    end_rle = start_rle + 31 < n_samples ? start_rle + 31 : n_samples;
                    // Update observed path.
                    observed_length += start_rle + 31 < n_samples ? 31 : n_samples - start_rle;
                    // Assertion.
                    assert(end_rle <= n_samples);
                    
                    // Debug:
                    // std::cerr << "Use Bitmap=" << (end_rle-start_rle) << "(" << start_rle << "," << end_rle << ")" << std::endl;
                    
                    // Construct 32-bit bitmap
                    uint32_t bitmap = 0;
                    for (int j = start_rle; j < end_rle; ++j) {
                        bitmap |= (prev[j] & 1);
                        bitmap <<= 1;
                    }
                    assert((bitmap & 1) == 0);
                    observed_alts += __builtin_popcount(bitmap);
                    
                    memcpy(&data[data_len], &bitmap, sizeof(uint32_t));
                    data_len += sizeof(uint32_t);

                    start_rle = i;
                    // Set new reference.
                    ref = prev[i];
                    // Set run-length. If this is the final object then set to 0
                    // for downstream filter.
                    rle_len = (end_rle == n_samples) ? 0 : 1;
                    // Update cost.
                    rle_cost += sizeof(uint32_t);
                    continue;
                   
                } else {
                    // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
                    
                    // Model:
                    // Lowest bit => 1 iff RLE, 0 otherwise
                    // If RLE => bit 2 is the reference
                    // Remainder => data (bitmap / run-length)
                    uint16_t rle_pack = (rle_len << 2) | ((ref & 1) << 1) | 1;
                    memcpy(&data[data_len], &rle_pack, sizeof(uint16_t));
                    data_len += sizeof(uint16_t);
                    
                    // Update observed alts.
                    observed_alts += (ref == 1) * rle_len; // predicate multiply
                    // Update observed path.
                    observed_length += rle_len;

                    // Reset length;
                    rle_len = 0;
                    // Set new reference.
                    ref = prev[i];
                    // Update cost.
                    rle_cost += sizeof(uint16_t);
                    // Update start position.
                    start_rle = i;
                }
            }
            ++rle_len; 
        }

        // Add only if last element is an unfinished RLE.
        if (rle_len) {
            // std::cerr << "Add final=" << rle_len << std::endl;
            ++n_objects;

            uint16_t rle_pack = (rle_len << 2) | ((ref & 1) << 1) | 1;
            memcpy(&data[data_len], &rle_pack, sizeof(uint16_t));
            data_len += sizeof(uint16_t);

            rle_cost += sizeof(uint16_t);
            observed_length += rle_len;
        }
        observed_alts += (ref == 1) * rle_len;

        // Ascertain correctness:
        assert(observed_alts == base_models[target].pbwt->n_queue[1]);
        assert(observed_length == n_samples);

        ++gt_width[n_objects];

        return 1;
    }

public:
    uint64_t bytes_in1;
    uint64_t bytes_out_zstd1, bytes_out_lz4;
    Buffer buf_wah[2];
    std::vector<uint32_t> gt_width;

};

#endif