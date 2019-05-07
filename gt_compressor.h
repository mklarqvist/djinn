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
#include "gt_decompressor.h"


#include <openssl/sha.h>

// #define DEBUG_PBWT 1
#define DEBUG_WAH 1

struct DataDigest {
    DataDigest() : len(0), capac(0), buffer(nullptr), has_initialized(false) {}
    DataDigest(uint32_t l) : len(0), capac(l), buffer(new uint8_t[l]), has_initialized(false) {}
    ~DataDigest() { delete[] buffer; }

    void resize(uint32_t l) {
        delete[] buffer;
        len = 0; capac = l;
        buffer = new uint8_t[l];
        has_initialized = false;
    }

    void reset() { len = 0; has_initialized = false; }

    inline bool InitDigest() {
		if (!SHA512_Init(&context)) return false;
        has_initialized = true;
		return true;
	}

    bool UpdateDigest(const uint8_t* data, const uint32_t data_len) {
		if (!has_initialized) {
            bool init_passed = this->InitDigest();
            if (!init_passed) return false;
        }
        // std::cerr << "Adding: " << len << "->" << len+data_len << "(" << data_len << ")" << "/" << capac << std::endl;
        assert(len+data_len < capac);
        memcpy(&buffer[len], data, data_len);
        len += data_len;

		if (!SHA512_Update(&context, data, data_len))
			return false;

		return true;
	}

    bool UpdateDigestStride(const uint8_t* data, const uint32_t data_len, const int stride = 1) {
		if (!has_initialized) {
            bool init_passed = this->InitDigest();
            if (!init_passed) return false;
        }
        // std::cerr << "Adding: " << len << "->" << len+data_len << "(" << data_len << ")" << "/" << capac << std::endl;
        // assert(len+data_len < capac);
        // memcpy(&buffer[len], data, data_len);
        // len += data_len;
        for (int i = 0; i + stride <= data_len; i += stride) {
            buffer[len++] = data[i];
        }

		if (!SHA512_Update(&context, data, data_len))
			return false;

		return true;
	}

    bool FinalizeDigest() {
		if (!has_initialized) return false;

		if (!SHA512_Final(&digest[0], &context))
			return false;

        has_initialized = false;
        len = 0;
		return true;
	}

    uint32_t len;
    uint32_t capac;
    uint8_t* buffer;

    bool has_initialized;
    SHA512_CTX context;
	uint8_t    digest[64];
};

// temp
// Ascertain that the binary output can be used to restore the input data.
// Computes allele counts as an example.
static
int64_t DebugRLEBitmap(const uint8_t* data, const uint32_t data_len, const uint32_t n_samples, const bool print = false) {
    // Data.
    // const uint8_t* data = buf_wah[target].data;
    // const uint32_t data_len = buf_wah[target].len;
    
    // Setup.
    uint32_t n_run = 0;
    uint32_t offset = 0;
    int64_t n_variants = 0;
    uint32_t n_alts = 0;

    while (true) {
        // Looking at most-significant byte for the target compression type.
        const uint8_t type = (data[offset] & 1);
        if (type == 0) {
            // assert((*reinterpret_cast<const uint32_t*>(&data[offset]) & 1) == 0);
            n_alts += __builtin_popcount(*reinterpret_cast<const uint32_t*>(&data[offset]) >> 1);
            offset += sizeof(uint32_t);
            n_run  += 31;
            if (n_run >= n_samples) {
                ++n_variants;
                n_run = 0;
                // std::cout << "AC=" << n_alts << '\n';
                if (print && (n_variants % 2) == 0) printf("AC=%d\n",n_alts);
                assert(n_alts <= 2*n_samples);
                if ((n_variants % 2) == 0) n_alts = 0;
            }
        }
        else if (type == 1) {
            uint16_t val = *reinterpret_cast<const uint16_t*>(&data[offset]);
            n_run  += (val >> 2);
            n_alts += ((val >> 1) & 1) * (val >> 2);

            // assert((*reinterpret_cast<const uint16_t*>(&data[offset]) & 1) == 1);
            offset += sizeof(uint16_t);
            if (n_run >= n_samples) {
                ++n_variants;
                n_run = 0;
                // std::cout << "AC=" << n_alts << '\n';
                if (print && (n_variants % 2) == 0) printf("AC=%d\n",n_alts);
                assert(n_alts <= 2*n_samples);
                if ((n_variants % 2) == 0) n_alts = 0;
            }
        }
        
        // Exit conditions.
        if (offset == data_len) {
            // std::cerr << "exit correct" << std::endl;
            break;
        }
        if (offset > data_len) {
            std::cerr << "overflow error: " << offset << "/" << data_len << std::endl;
            exit(1);
        }
    }

    // std::cerr << "Number of variants=" << (n_variants>>1) << std::endl;

    return n_variants;
}

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

    //
#if DEBUG_PBWT
    int32_t DebugPBWT() {
        std::shared_ptr<PBWT> pbwt1 = std::make_shared<PBWT>(n_samples, 2);  
        std::shared_ptr<PBWT> pbwt2 = std::make_shared<PBWT>(n_samples, 2);  

        uint32_t offset1 = 0, offset2 = 0;
        for (int i = 0; i < processed_lines_local; ++i) {
            pbwt1->ReverseUpdate(&debug_pbwt[0].buffer[offset1]);
            pbwt2->ReverseUpdate(&debug_pbwt[1].buffer[offset2]);
            offset1 += n_samples;
            offset2 += n_samples;
        }

        return 1;
    }
#endif

    //
    virtual bool CheckLimit() const =0;
    virtual int Compress() =0;
    
protected:
    int64_t  n_samples;
    uint32_t block_size;
    uint32_t processed_lines;
    uint32_t processed_lines_local;
    
    uint64_t bytes_in, bytes_out;
    
    Buffer buf_compress;
    Buffer buf_raw;

    uint32_t alts[256];

    GeneralPBWTModel base_models[4]; // 0-1: diploid biallelic no-missing; 2-3: diploid biallelic missing
    GeneralPBWTModel base_models_complex[2]; // 0-1: diploid n-allelic
    GeneralPBWTModel* models;
#if DEBUG_PBWT
    DataDigest debug_pbwt[2];
#endif
};

class GenotypeCompressorModelling : public GenotypeCompressor {
public:
    GenotypeCompressorModelling(int64_t n_s);
    ~GenotypeCompressorModelling();

    bool CheckLimit() const override;

    //
    // int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) override;

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

    int Compress() override;

private:
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

    bool CheckLimit() const override;

    int Compress() override;
    int EncodeRLEBitmap(const int target);

private:
    // Compression strategy used.
    enum class CompressionStrategy : uint32_t { ZSTD = 0, LZ4 = 1 };

private:
    CompressionStrategy strategy;
    uint64_t bytes_in1;//debug
    uint64_t bytes_out_zstd1, bytes_out_lz4;//debug
    Buffer buf_wah[2];
    std::vector<uint32_t> gt_width; // debug
    PBWT base_pbwt[2]; // diploid biallelic no-missing models
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