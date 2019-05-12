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

#include "djinn.h"
#include "vcf_reader.h"
#include "pbwt.h"
#include "compressors.h"

// temp
#include <bitset>
#include <chrono>

#include "gt_decompressor.h"
#include <openssl/sha.h>
#include <roaring/roaring.h>


#include "range_coder64.h"

namespace djinn {

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
        assert(len + data_len < capac);
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

        // if (len + (data_len/stride) >= capac) std::cerr << "Adding: " << len << "->" << len+(data_len/stride) << "(" << data_len/stride << ")" << "/" << capac << std::endl;
        assert(len + (data_len/stride) < capac);
        // memcpy(&buffer[len], data, data_len);
        // len += data_len;
        for (int i = 0; i + stride <= data_len; i += stride) {
            buffer[len++] = data[i];
        }
        // assert(len < capac);

		if (!SHA512_Update(&context, data, data_len))
			return false;

		return true;
	}

    bool FinalizeDigest() {
		if (!has_initialized) return false;

		if (!SHA512_Final(&digest[0], &context))
			return false;

        has_initialized = false;
        // len = 0;
		return true;
	}

    uint32_t len;
    uint32_t capac;
    uint8_t* buffer;

    bool has_initialized;
    SHA512_CTX context;
	uint8_t    digest[64];
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


// Forward declare for friendship
class GTCompressor;

class GenotypeCompressor {
public:
    GenotypeCompressor(int64_t n_s);
    virtual ~GenotypeCompressor();
    
    // Encode data using literals.
    virtual int Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles) =0;
    virtual int Encode2N2M(uint8_t* data, const int32_t n_data) =0;
    virtual int Encode2N2MC(uint8_t* data, const int32_t n_data) =0;
    virtual int Encode2N2MM(uint8_t* data, const int32_t n_data) =0;
    virtual int Encode2NXM(uint8_t* data, const int32_t n_data, const int32_t n_alleles) =0;

    // Encode data using htslib.
    virtual int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr);
    virtual int Encode(uint8_t* data, const int32_t n_data, const int32_t ploidy, const int32_t n_alleles);
    virtual int Encode2N(bcf1_t* bcf, const bcf_hdr_t* hdr) =0;
    //
    int32_t RemapGenotypeEOV(uint8_t* data, const uint32_t len);

    //
    virtual bool CheckLimit() const =0;
    virtual int Compress() =0;
    // virtual int Compress(djinn_block_t*& block) =0;

public:
    // Compression strategy used.
    friend GTCompressor; // Wrapper class friendship.
    enum class CompressionStrategy : uint32_t { ZSTD = 0, LZ4 = 1, CONTEXT = 2 };
    
protected:
    bool permute_pbwt;
    int64_t  n_samples;
    uint32_t block_size;
    uint32_t processed_lines;
    uint32_t processed_lines_local;
    
    uint64_t bytes_in, bytes_in_vcf, bytes_out;
    CompressionStrategy strategy;
    
    Buffer buf_compress;
    Buffer buf_raw;

    uint32_t alts[256];

    djinn_block_t* block;

#if DEBUG_PBWT
    DataDigest debug_pbwt[6]; // 2 complete, 2 missing, 2 complex
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
    int Encode2NXM(uint8_t* data, const int32_t n_data, const int32_t n_alleles) override;

    int Compress() override;

private:
    GeneralPBWTModel base_models[4]; // 0-1: diploid biallelic no-missing; 2-3: diploid biallelic missing
    GeneralPBWTModel base_models_complex[2]; // 0-1: diploid n-allelic
    GeneralPBWTModel* models;
    GeneralPBWTModel base_model_bitmaps[2];

    // SimpleModel64 s64[2];
};

class GenotypeCompressorRLEBitmap : public GenotypeCompressor {
public:
    GenotypeCompressorRLEBitmap(int64_t n_s);
    ~GenotypeCompressorRLEBitmap();

    int Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles) override;
    int Encode2N2M(uint8_t* data, const int32_t n_data) override;
    int Encode2N2MC(uint8_t* data, const int32_t n_data) override;
    int Encode2N2MM(uint8_t* data, const int32_t n_data) override;
    int Encode2NXM(uint8_t* data, const int32_t n_data, const int32_t n_alleles) override;

    // Encode data using htslib.
    int Encode2N(bcf1_t* bcf, const bcf_hdr_t* hdr) override { return(-1); }

    bool CheckLimit() const override;

    int Compress() override;

    // Compression routines.
    int EncodeRLEBitmap2N2MC(const int target);
    int EncodeRLEBitmap2N2MM(const int target);
    int EncodeRLEBitmap2NXM(const int target);
    int EncodeRLEBitmap2N2MC(const int target, const uint8_t* data, const int stride);
    int EncodeRLEBitmap2N2MM(const int target, const uint8_t* data, const int stride);
    int EncodeRLEBitmap2NXM(const int target, const uint8_t* data, const int stride);

private:
    uint64_t bytes_out_zstd1, bytes_out_lz4;//debug
    Buffer buf_wah[3];
    std::vector<uint32_t> gt_width; // debug
    PBWT base_pbwt[4]; // diploid biallelic no-missing models, missing + EOV
    PBWT complex_pbwt[2]; // diploid n-allelic
};


class HaplotypeCompressor {
public:
    HaplotypeCompressor(int64_t n_s) : n_samples(n_s), cum_pos(0), bitmaps(new roaring_bitmap_t*[n_samples])
    {
        for (int i = 0; i < n_samples; ++i) bitmaps[i] = roaring_bitmap_create();  
    }

    ~HaplotypeCompressor() {
        for (int i = 0; i < n_samples; ++i) roaring_bitmap_free(bitmaps[i]);
        delete[] bitmaps;
    }

    int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) {
        if (bcf == NULL) return 0;
        if (hdr == NULL) return 0;

        const bcf_fmt_t* fmt = bcf_get_fmt(hdr, bcf, "GT");
        if (fmt == NULL) return 0;
        if (fmt->p_len != n_samples) {
            std::cerr << "input is not divisible by 2" << std::endl;
            return 0;
        }

        if (bcf->n_allele != 2) return 0;
        assert(fmt->p_len == n_samples);

        const uint8_t* data = fmt->p;
        for (int i = 0; i < n_samples; ++i) {
            uint8_t add = BCF_UNPACK_GENOTYPE(data[i]);
            if (add == 1)
                roaring_bitmap_add(bitmaps[i], cum_pos);
        }
        ++cum_pos;

        return 1;
    }
    
    void PrintSizes() const {
        uint64_t total_cost = 0;
        uint32_t max_size = 0;

        for (int i = 0; i < n_samples; ++i) {
            uint32_t cost = roaring_bitmap_size_in_bytes(bitmaps[i]);
            max_size = std::max(max_size, cost);
            std::cerr << "Haplotype-" << i << ": " << cost << std::endl;
            total_cost += cost;
        }
        std::cerr << "Total cost=" << total_cost << std::endl;

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        for (int j = 0; j < n_samples; ++j) {
            std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
            uint64_t sim = 0;
            for (int i = j+1; i < n_samples; ++i) {
                sim += roaring_bitmap_and_cardinality(bitmaps[i], bitmaps[j]);
            }
            std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
            auto time_span_inner = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);
            std::cerr << j << "->" << sim << " time=" << time_span_inner.count() << "ms" << std::endl;
        }
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
        std::cerr << "Total time=" << time_span.count() << "s" << std::endl;
        

    }

public:
    int64_t n_samples;
    int64_t cum_pos;
    roaring_bitmap_t** bitmaps;
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

    GTCompressor& SetStrategy(GenotypeCompressor::CompressionStrategy strategy) {
        if (instance.get() == nullptr) return(*this);

        switch(strategy) {
            case(GenotypeCompressor::CompressionStrategy::CONTEXT): instance->strategy = GenotypeCompressor::CompressionStrategy::CONTEXT; break;
            case(GenotypeCompressor::CompressionStrategy::LZ4): instance->strategy = GenotypeCompressor::CompressionStrategy::LZ4; break;
            case(GenotypeCompressor::CompressionStrategy::ZSTD): instance->strategy = GenotypeCompressor::CompressionStrategy::ZSTD; break;
        }
        return *this;
    }

public:
    inline int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) { this->instance->Encode(bcf, hdr); }
    inline int Encode(uint8_t* data, const int32_t n_data, const int n_ploidy, const int n_alleles) { this->instance->Encode(data, n_data, n_ploidy, n_alleles); }
    inline void Compress() { this->instance->Compress(); }

private:
    CompressionStrategy strategy;
    std::shared_ptr<GenotypeCompressor> instance;
};

}

#endif