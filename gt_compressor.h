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

#include "ctx_model.h"

// temp
#include <bitset>
#include <chrono>

#include "gt_decompressor.h"

#if DEBUG_PBWT
#include <openssl/sha.h>
#endif

namespace djinn {

template<typename T>
uint32_t ilog2(T x)
{
	uint32_t r = 0;

	for (; x; ++r)
		x >>= 1;

	return r;
}

#if DEBUG_PBWT
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

    bool UpdateDigest(uint8_t data) {
		if (!has_initialized) {
            bool init_passed = this->InitDigest();
            if (!init_passed) return false;
        }

        // std::cerr << "Adding: " << len << "->" << len+data_len << "(" << data_len << ")" << "/" << capac << std::endl;
        // assert(len + data_len < capac);
        // memcpy(&buffer[len], data, data_len);
        // len += data_len;
        buffer[len++] = data;

		if (!SHA512_Update(&context, &data, 1))
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
#endif

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

    void SetPermutePbwt(const bool yes) { permute_pbwt = yes; }
    
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
    virtual int Compress(djinn_block_t*& block) =0;

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

#if DEBUG_PBWT
    DataDigest debug_pbwt[6]; // 2 complete, 2 missing, 2 complex
#endif
};

class GenotypeCompressorModelling : public GenotypeCompressor {
public:
    GenotypeCompressorModelling(int64_t n_s);
    ~GenotypeCompressorModelling();

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

    int Compress(djinn_block_t*& block) override;

#if DEBUG_CONTEXT
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
    typedef int(GenotypeDecompressorContext::*context_debug_decode)(uint8_t*);
    int DebugContext(uint8_t* in, size_t len_in, uint8_t* in_part, size_t len_part, uint8_t* ref_data, size_t n_cycles, int pbwt_sym, context_debug_decode decode_fn, const uint8_t* lookup_fn);
    int DebugContext(uint8_t* in, size_t len_in, uint8_t* ref_data, size_t n_cycles, int pbwt_sym, context_debug_decode decode_fn, const uint8_t* lookup_fn);
#endif

public:
   
    uint64_t bytes_out2, bytes_out3, bytes_out4;

    djinn_ctx_model djn_ctx, djn_ctx_decode;
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

    int Compress(djinn_block_t*& block) override;

    // Compression routines.
    int EncodeRLEBitmap2N2MC(const int target);
    int EncodeRLEBitmap2N2MM(const int target);
    int EncodeRLEBitmap2NXM(const int target);
    int EncodeRLEBitmap2N2MC(const int target, const uint8_t* data, const int stride);
    int EncodeRLEBitmap2N2MM(const int target, const uint8_t* data, const int stride);
    int EncodeRLEBitmap2NXM(const int target, const uint8_t* data, const int stride);

#if DEBUG_WAH
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
    typedef int(GenotypeDecompressorRLEBitmap::*bitmap_debug_decode)(uint8_t*, const uint32_t);
    int DebugWAH(uint8_t* in, size_t len_in, uint8_t* ref_data, size_t n_cycles, int pbwt_sym, bitmap_debug_decode decode_fn, const uint8_t* lookup_fn);
#endif

private:
    uint64_t bytes_out_zstd1, bytes_out_lz4;//debug
    size_t n_variants[6];
    Buffer buf_wah[6];
    // std::vector<uint32_t> gt_width; // debug
    PBWT base_pbwt[4]; // diploid biallelic no-missing models, missing + EOV
    PBWT complex_pbwt[2]; // diploid n-allelic
};

class HaplotypeCompressor {
public:
    HaplotypeCompressor(int64_t n_s) {
        n_samples = 2*n_s;
        offset = 0;
        n_in = 0;
        bitmaps = new uint64_t*[n_samples];
        for (int i = 0; i < n_samples; ++i) {
            bitmaps[i] = new uint64_t[8192];
            memset(bitmaps[i], 0, 8192*sizeof(uint64_t));
        }
    }

    ~HaplotypeCompressor() {
        if (bitmaps != nullptr) {
            for (int i = 0; i < n_samples; ++i) delete[] bitmaps[i];
            delete[] bitmaps;
        }
    }

    int EncodeBitmap(bcf1_t* bcf, const bcf_hdr_t* hdr) {
        if (bcf == NULL) return 0;
        if (hdr == NULL) return 0;

        const bcf_fmt_t* fmt = bcf_get_fmt(hdr, bcf, "GT");
        if (fmt == NULL) return 0;
        if (fmt->p_len != n_samples) {
            std::cerr << "input is not equal" << std::endl;
            return 0;
        }

        if (bcf->n_allele != 2) return 0;
        n_in += fmt->p_len;

        for (int i = 0; i < n_samples; ++i) {
            // 1 bit reserved for archetyping.
            bitmaps[i][offset / 63] |= (uint64_t)BCF_UNPACK_GENOTYPE(fmt->p[i]) << (offset % 63);
        }
        ++offset;

        return 1;
    }

    int Compress() {
        uint8_t* buf = new uint8_t[2000000];
        uint8_t* buf2 = new uint8_t[2000000];
        uint32_t buf2_off = 0;

        uint32_t n_bitmaps = offset/63;
        PBWT test(n_bitmaps*63,2);
        uint8_t* unpack = new uint8_t[n_bitmaps*63];
        for (int i = 0; i < n_samples; ++i) {
            uint32_t of = 0;
            for (int j = 0; j < n_bitmaps; ++j) {
                uint64_t ref = bitmaps[i][j];
                for (int k = 0; k < 63; ++k, ++of) {
                    unpack[of] = ref & 1;
                    ref >>= 1;
                }
            }
            test.Update(unpack, 1);
        }

        std::cerr << test.ToPrettyString() << std::endl;

        uint64_t* pbwt_bitmap = new uint64_t[n_bitmaps];
        for (int i = 0; i < n_samples; ++i) {
            memset(pbwt_bitmap,0,sizeof(uint64_t)*n_bitmaps);
            uint32_t of = 0;
            for (int j = 0; j < n_bitmaps; ++j) {
                uint64_t ref = bitmaps[i][j];
                for (int k = 0; k < 63; ++k, ++of) {
                    unpack[of] = ref & 1;
                    ref >>= 1;
                }
            }
            
            // assert(test.)
            // std::cerr << "repacking" << std::endl;
            for (int j = 0; j < n_bitmaps*63; ++j) {
                pbwt_bitmap[j / 63] |= (uint64_t)unpack[test.ppa[j]] << (j % 63);
            }

            // std::cerr << "wah+compress" << std::endl;

            // compress repacked
            uint64_t bref = pbwt_bitmap[0];
            uint32_t n_run = 1;
            uint32_t n_len = 0;
            uint32_t cost = 0;
            buf2_off = 0;
        
            uint32_t types[2] = {0};

            for (int j = 1; j < n_bitmaps; ++j) {
                // std::cerr << j << "/" << n_bitmaps << std::endl;
                if (pbwt_bitmap[j] != bref) {
                    if (n_run == 1) {
                        // std::cerr << std::bitset<32>(bref) << " ";
                        n_len += 63;
                        
                        *reinterpret_cast<uint64_t*>(&buf2[buf2_off]) = bref;
                        buf2_off += sizeof(uint64_t);

                        cost += sizeof(uint64_t);
                        ++types[0];
                    }
                    else {
                        // std::cerr << n_run << ":" << std::bitset<32>(bref) << " ";
                        cost += sizeof(uint32_t);
                        n_len += 63*n_run;
                        ++types[1];

                        uint32_t prun = ((bref & 1) << 31) | n_run;
                        *reinterpret_cast<uint32_t*>(&buf2[buf2_off]) = prun;
                        buf2_off += sizeof(uint32_t);
                    }

                    n_run = 0;
                    bref = pbwt_bitmap[j];
                }
                ++n_run;
            }

            if (n_run) {
                if (n_run == 1) {
                    // std::cerr << std::bitset<32>(bref) << " ";
                    n_len += 63;

                    *reinterpret_cast<uint64_t*>(&buf2[buf2_off]) = bref;
                    buf2_off += sizeof(uint64_t);

                    cost += sizeof(uint64_t);
                    ++types[0];
                }
                else {
                    // std::cerr << n_run << ":" << std::bitset<32>(bref) << " ";
                    cost += sizeof(uint32_t);
                    n_len += 63*n_run;
                    ++types[1];

                    uint32_t prun = ((bref & 1) << 31) | n_run;
                    *reinterpret_cast<uint32_t*>(&buf2[buf2_off]) = prun;
                    buf2_off += sizeof(uint32_t);
                }
            }
            // std::cerr << "compressing=" << buf2_off << std::endl;
            size_t praw = Lz4Compress(buf2, buf2_off, buf, 2000000, 9);
            
            std::cerr << "Sample-" << i << ": " << offset << "->" << praw << " (" << (float)offset/praw << "-fold) " << offset << "->" << offset/63 << " typing=" << (float)types[0]/(types[0]+types[1]) << " and " << (float)types[1]/(types[0]+types[1]) << std::endl;
            

        }

        delete[] pbwt_bitmap;
        delete[] unpack;
        delete[] buf2;


        uint64_t total = 0;
        for (int i = 0; i < n_samples; ++i) {
            uint64_t bref = bitmaps[i][0];
            uint32_t n_run = 1;
            uint32_t n_len = 0;
            uint32_t cost = 0;
        
            uint32_t types[2] = {0};

            
            for (int j = 1; j < n_bitmaps; ++j) {
                if (bitmaps[i][j] != bref) {
                    if (n_run == 1) {
                        // std::cerr << std::bitset<32>(bref) << " ";
                        n_len += 63;
                        cost += sizeof(uint64_t);
                        ++types[0];
                    }
                    else {
                        // std::cerr << n_run << ":" << std::bitset<32>(bref) << " ";
                        cost += sizeof(uint32_t);
                        n_len += 63*n_run;
                        ++types[1];
                    }

                    n_run = 0;
                    bref = bitmaps[i][j];
                }
                ++n_run;
            }

            if (n_run) {
                if (n_run == 1) {
                    // std::cerr << std::bitset<32>(bref) << " ";
                    n_len += 63;
                    cost += sizeof(uint64_t);
                    ++types[0];
                }
                else {
                    // std::cerr << n_run << ":" << std::bitset<32>(bref) << " ";
                    cost += sizeof(uint32_t);
                    n_len += 63*n_run;
                    ++types[1];
                }
            }
            size_t praw = Lz4Compress((uint8_t*)&bitmaps[i][0], n_bitmaps*sizeof(uint64_t), buf, 2000000, 9);
            
            std::cerr << "Sample-" << i << ": " << cost << "->" << praw << " (" << (float)offset/praw << "-fold) " << offset << "->" << offset/63 << " typing=" << (float)types[0]/(types[0]+types[1]) << " and " << (float)types[1]/(types[0]+types[1]) << std::endl;
            total += cost;
        }
        offset = 0;
        std::cerr << "[hapWAH] " << n_in << "->" << total << " (" << (float)n_in/total << "-fold)" << std::endl;
        n_in = 0;

        for (int i = 0; i < n_samples; ++i) {
            memset(bitmaps[i], 0, 8192*sizeof(uint64_t));
        }

        delete[] buf;

        return total;
    }

    int64_t n_samples;
    int64_t n_in;
    // std::vector<std::shared_ptr<Buffer>> haps;
    uint32_t offset;
    uint64_t** bitmaps;
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
            case(GenotypeCompressor::CompressionStrategy::LZ4):     instance->strategy = GenotypeCompressor::CompressionStrategy::LZ4;     break;
            case(GenotypeCompressor::CompressionStrategy::ZSTD):    instance->strategy = GenotypeCompressor::CompressionStrategy::ZSTD;    break;
        }
        return *this;
    }

    inline void SetPermutePbwt(const bool yes) { this->instance->SetPermutePbwt(yes); }

    inline int Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) { return this->instance->Encode(bcf, hdr); }
    inline int Encode(uint8_t* data, const int32_t n_data, const int n_ploidy, const int n_alleles) { return this->instance->Encode(data, n_data, n_ploidy, n_alleles); }
    inline int Compress(djinn_block_t*& block) { return this->instance->Compress(block); }

public:
    CompressionStrategy strategy;
    std::shared_ptr<GenotypeCompressor> instance;
};

}

#endif