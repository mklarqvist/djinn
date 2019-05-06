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

#include "compressors.h"

// temp
#include <bitset>

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
    GenotypeCompressor(int64_t n_s);
    virtual ~GenotypeCompressor();
    
    // Encode data using literals.
    // virtual int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) =0;
    virtual int Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles) =0;
    virtual int Encode2N2M(uint8_t* data, const int32_t n_data) =0;
    virtual int Encode2N2MC(uint8_t* data, const int32_t n_data) =0;
    virtual int Encode2N2MM(uint8_t* data, const int32_t n_data) =0;
    virtual int Encode2NXM(uint8_t* data, const int32_t n_data) =0;

    // Encode data using htslib.
    virtual int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr);
    virtual int Encode2N(bcf1_t* bcf, const bcf_hdr_t* hdr) =0;

    //
    int32_t RemapGenotypeEOV(uint8_t* data, const uint32_t len);
    
public:
    int64_t  n_samples;
    uint32_t block_size;
    uint32_t processed_lines;
    uint32_t processed_lines_local;
    
    uint64_t bytes_in, bytes_out;
    
    Buffer buf_compress;
    Buffer buf_raw;

    uint32_t alts[256];

    GeneralPBWTModel base_models[4];
    GeneralPBWTModel base_models_complex[2];
    GeneralPBWTModel* models;
};

class GenotypeCompressorModelling : public GenotypeCompressor {
public:
    GenotypeCompressorModelling(int64_t n_s);
    ~GenotypeCompressorModelling();

    bool CheckLimit() const;

    //
    int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) override;

private:
    // Diploid wrapper.
    int Encode2N(bcf1_t* bcf, const bcf_hdr_t* hdr) override;
    int Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles) override;
    int Encode2N2M(uint8_t* data, const int32_t n_data) override;
    int Encode2N2MC(uint8_t* data, const int32_t n_data) override;

    // 2N2M with missing
    int Encode2N2MM(uint8_t* data, const int32_t n_data) override;

    // 2N any M (up to 16)
    int Encode2NXM(uint8_t* data, const int32_t n_data) override;

    int Compress();

public:
    GeneralPBWTModel nonsense[2];
};

class GenotypeCompressorRLEBitmap : public GenotypeCompressor {
public:
    GenotypeCompressorRLEBitmap(int64_t n_s);
    ~GenotypeCompressorRLEBitmap();

    int Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles) override;

    int Encode2N2M(uint8_t* data, const int32_t n_data) override { return(-1); }
    int Encode2N2MC(uint8_t* data, const int32_t n_data) override { return(-1); }
    int Encode2N2MM(uint8_t* data, const int32_t n_data) override { return(-1); }
    int Encode2NXM(uint8_t* data, const int32_t n_data) override { return(-1); }

    // Encode data using htslib.
    int Encode2N(bcf1_t* bcf, const bcf_hdr_t* hdr) override { return(-1); }

    bool CheckLimit() const;

    int Compress();
    int EncodeRLEBitmap(const int target);

protected:
    // Compression strategy used.
    enum class CompressionStrategy : uint32_t { ZSTD = 0, LZ4 = 1 };

public:
    CompressionStrategy strategy;
    uint64_t bytes_in1;
    uint64_t bytes_out_zstd1, bytes_out_lz4;
    Buffer buf_wah[2];
    std::vector<uint32_t> gt_width;
    PBWT base_pbwt[2];
};

class GTCompressor {
public:
    // Compression strategy used.
    enum class CompressionStrategy : uint32_t { CONTEXT_MODEL = 0, RLE_BITMAP = 1 };

public:
    GTCompressor& SetStrategy(CompressionStrategy strategy, int64_t n_samples) {
        switch(strategy) {
        case(CompressionStrategy::CONTEXT_MODEL):
            instance = std::make_shared<GenotypeCompressorModelling>(n_samples);
            this->strategy = strategy;
            break;
        case(CompressionStrategy::RLE_BITMAP):
            instance = std::make_shared<GenotypeCompressorRLEBitmap>(n_samples);
            this->strategy = strategy;
            break;
        }
        return *this;
    }

public:
    inline int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) { this->instance->Encode(bcf, hdr); }
    // inline int Encode(uint8_t* data, const int32_t n_data, const int n_alleles) { this->instance->Encode(data, n_data, n_alleles); }

private:
    CompressionStrategy strategy;
    std::shared_ptr<GenotypeCompressor> instance;
};

#endif