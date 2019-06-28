
#include "gt_compressor.h"

#include <cstdio>//printf: debug
#include "gt_decompressor.h" //debug

namespace djinn {

/*======   Base model   ======*/

GenotypeCompressor::GenotypeCompressor(int64_t n_s) :
    permute_pbwt(true),
    n_samples(n_s),
    block_size(8192),
    processed_lines(0),
    processed_lines_local(0),
    bytes_in(0), bytes_in_vcf(0), bytes_out(0),
    strategy(CompressionStrategy::LZ4)
{
    buf_compress.resize(10000000);
    buf_raw.resize(10000000);

#if DEBUG_PBWT
    for (int i = 0; i < 6; ++i) {
        debug_pbwt[i].resize(n_samples*block_size + 65536);
    }
#endif
}

GenotypeCompressor::~GenotypeCompressor() { }

// Encode data using htslib.
int GenotypeCompressor::Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) {
    if (bcf == NULL) return 0;
    if (hdr == NULL) return 0;

    const bcf_fmt_t* fmt = bcf_get_fmt(hdr, bcf, "GT");
    if (fmt == NULL) return 0;
    if (fmt->p_len / n_samples != 2) {
        std::cerr << "input is not divisible by 2: " << fmt->p_len << "/" << n_samples << std::endl;
        return 0;
    }

    // Keep track of input bytes.
    bytes_in += fmt->p_len;
    bytes_in_vcf += fmt->p_len * 2 - 1; // (char)(sep)(char)(tab)

    // This is not used: results in worse compression performance across
    // the board.
    bool permute = true;
    if (bcf->n_allele > 2) {

    } else {
        if (strlen(bcf->d.allele[0]) != 1) {
            // std::cerr << "not snv: " << bcf->d.allele[0] << std::endl;
            // permute = false;
        }
        if (strlen(bcf->d.allele[1]) != 1) {
            // std::cerr << "not snv: " << bcf->d.allele[1] << std::endl;
            // permute = false;
        }
    }

    // if (bcf->n_allele < 16 && bcf->n_allele > 2) {
    //     for (int i = 0; i < bcf->n_allele; ++i) {
    //         std::cerr << bcf->d.allele[i] << ",";
    //     }
    //     std::cerr << std::endl;
    // }

    // Todo: extend beyond 2N
    return(Encode2N(fmt->p, fmt->p_len, bcf->n_allele, permute));
}

int GenotypeCompressor::Encode(uint8_t* data, const int32_t n_data, const int32_t ploidy, const int32_t n_alleles, const bool permute) {
    if (data == nullptr) return 0;
    if (n_data == 0) return 0;
    if (n_alleles <= 0) return 0;

    if (ploidy != 2) {
        std::cerr << "only diploid data is currently supported" << std::endl;
        return 0;
    }

    bytes_in += n_data;
    bytes_in_vcf += n_data * 2 - 1;

    return(Encode2N(data, n_data, n_alleles, permute));
}

//
int32_t GenotypeCompressor::RemapGenotypeEOV(uint8_t* data, const uint32_t len) {
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
        // std::cerr << "replacing 64 with " << (int32_t)max_val << std::endl;
        for (int i = 0; i < len; ++i) {
            if (BCF_UNPACK_GENOTYPE_GENERAL(data[i]) == 64) {
                data[i] = max_val;
                ++n_replaced;
            }
        }
        // std::cerr << "replaced=" << n_replaced << std::endl;
    }

    return n_replaced;
}

/*======   Context-modelling approach   ======*/

GenotypeCompressorModelling::GenotypeCompressorModelling(int64_t n_s) : GenotypeCompressor(n_s)
{
    strategy = CompressionStrategy::CONTEXT;
    bytes_out2 = 0; bytes_out3 = 0;

#if DEBUG_PBWT
    for (int i = 0; i < 2; ++i) {
        debug_bins[i].resize(n_samples*block_size + 65536);
    }
#endif

    // model
    bytes_out4 = 0;
    // djn_ctx.SetSamples(2*n_samples);
    djn_ctx.StartEncoding(true, false);
    // djn_ctx_decode.SetSamples(2*n_samples);

    debug_buffer = new char[10000000];
    len_debug = 0;
}

GenotypeCompressorModelling::~GenotypeCompressorModelling() { }

int GenotypeCompressorModelling::Encode2N(bcf1_t* bcf, const bcf_hdr_t* hdr) {
    if (bcf == NULL) return 0;
    if (hdr == NULL) return 0;
    
    const bcf_fmt_t* fmt = bcf_get_fmt(hdr, bcf, "GT");
    if (fmt == NULL) return 0;
    if (fmt->p_len / n_samples != 2) {
        std::cerr << "input is not divisible by 2" << std::endl;
        return 0;
    }

    if (bcf->n_allele == 2) return(Encode2N2M(fmt->p, fmt->p_len, true)); // biallelic
    else if (bcf->n_allele < 14) {

        int replaced = RemapGenotypeEOV(fmt->p, fmt->p_len);
        if (replaced) {
            std::cerr << "replaced: skipping " << replaced << std::endl;
            return Encode2N2MM(fmt->p, fmt->p_len, true); // 2N2MM supports missing+EOV in RLE-bitmap mode.
        }
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

int GenotypeCompressorModelling::Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles, const bool permute) {
    if (n_alleles == 2) return(Encode2N2M(data, n_data, permute)); // biallelic
    else if (n_alleles < 14) return(Encode2NXM(data, n_data, n_alleles, permute));  // #alleles < 14
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
int GenotypeCompressorModelling::Encode2N2M(uint8_t* data, const int32_t n_data, const bool permute) {
    // Todo: assert genotypes are set for this variant.
    int replaced = RemapGenotypeEOV(data, n_data);
    if (replaced) {
        std::cerr << "replaced: divert to missing with " << replaced << std::endl;
        return Encode2NXM(data ,n_data, 4, permute); // todo: see if this is true
    }

    // for (int i = 0; i < 256; ++i) {
        // if (alts[i]) std::cerr << i << ":" << alts[i] << ",";
    // }
    // std::cerr << std::endl;

    if (alts[14] == 0 && alts[15] == 0) { // No missing values.
        // std::cerr << "divert 2N2MC" << std::endl;
        return Encode2N2MC(data, n_data, permute);
    } else { // Having missing values.
        // std::cerr << "using extended model" << std::endl;
        // std::cerr << "divert 2N2MM" << std::endl;
        return Encode2N2MM(data, n_data, permute);
    }

    return 1;
}

// 2N2M complete
int GenotypeCompressorModelling::Encode2N2MC(uint8_t* data, const int32_t n_data, const bool permute) {
#if DEBUG_PBWT
    assert(debug_pbwt[0].UpdateDigestStride(&data[0], n_data, 2));
    assert(debug_pbwt[1].UpdateDigestStride(&data[1], n_data, 2));
#endif

    djn_ctx.EncodeBcf(data, 2*n_samples, 2, 2);

    ++processed_lines_local;
    ++processed_lines;
    return 1;
    
}

// 2N2M with missing
int GenotypeCompressorModelling::Encode2N2MM(uint8_t* data, const int32_t n_data, const bool permute) {
#if DEBUG_PBWT
    assert(debug_pbwt[2].UpdateDigestStride(&data[0], n_data, 2));
    assert(debug_pbwt[3].UpdateDigestStride(&data[1], n_data, 2));
#endif

    djn_ctx.EncodeBcf(data, 2*n_samples, 2, 15);

    // std::cerr << "After Encode2N2MM: " << 15 << std::endl;

    ++processed_lines_local;
    ++processed_lines;
    return 1;
}

// 2N any M (up to 16)
int GenotypeCompressorModelling::Encode2NXM(uint8_t* data, const int32_t n_data, const int32_t n_alleles, const bool permute) {

#if DEBUG_PBWT
    assert(debug_pbwt[4].UpdateDigestStride(&data[0], n_data, 2));
    assert(debug_pbwt[5].UpdateDigestStride(&data[1], n_data, 2));
#endif

    djn_ctx.EncodeBcf(data, 2*n_samples, 2, n_alleles);

    // std::cerr << "After Ecndoe2nXM: " << n_alleles << std::endl;

    ++processed_lines_local;
    ++processed_lines;
    return 1;
}


#if DEBUG_CONTEXT
int GenotypeCompressorModelling::DebugContext(uint8_t* in, size_t len_in, uint8_t* in_part, size_t len_part, uint8_t* ref_data, size_t n_cycles, int pbwt_sym, context_debug_decode decode_fn, const uint8_t* lookup_fn) {   
    if (len_in > 8) { // The return size is 8 bytes when empty.
        // PBWT has to be created even though it is not used.
        std::shared_ptr<GenotypeDecompressorContext> debug1 = std::make_shared<GenotypeDecompressorContext>(in, len_in, in_part, len_part, n_cycles, n_samples, pbwt_sym);
            
        // if (permute_pbwt)
        //     debug1->InitPbwt(pbwt_sym);

        //TWK_BCF_GT_UNPACK[(A >> 1)]

        uint8_t* out = new uint8_t[n_samples];
        uint32_t n_debug = 0;
        uint32_t debug_offset = 0;
        if (permute_pbwt) {
            for (int i = 0; i < n_cycles; ++i) {
                int32_t alts1 = CALL_MEMBER_FN(*debug1,decode_fn)(out);
                debug1->model->pbwt->ReverseUpdate(out);
                // assert(alts1 == debug1->pbwt->n_queue[1]);
                for (int j = 0; j < n_samples; ++j) {
                    assert(lookup_fn[ref_data[debug_offset + j] >> 1] == debug1->model->pbwt->prev[j]);
                }

                debug_offset += n_samples;
                ++n_debug;
            }
        } else {
            for (int i = 0; i < n_cycles; ++i) {
                int32_t alts1 = CALL_MEMBER_FN(*debug1,decode_fn)(out);
                for (int j = 0; j < n_samples; ++j) {
                    // assert(BCF_UNPACK_GENOTYPE_GENERAL(ref_data[debug_offset + j]) == out[j]);
                    assert(lookup_fn[ref_data[debug_offset + j] >> 1] == out[j]);
                }

                debug_offset += n_samples;
                ++n_debug;
            }
        }
        
        delete[] out;
        return 1;
    }
    return 1;
}

int GenotypeCompressorModelling::DebugContext(uint8_t* in, size_t len_in, uint8_t* ref_data, size_t n_cycles, int pbwt_sym, context_debug_decode decode_fn, const uint8_t* lookup_fn) {    
    if (len_in > 8) { // The return size is 8 bytes when empty.
        // PBWT has to be created even though it is not used.
        std::shared_ptr<GenotypeDecompressorContext> debug1 = std::make_shared<GenotypeDecompressorContext>(in, len_in, n_cycles, n_samples, pbwt_sym);
            
        // if (permute_pbwt)
        //     debug1->InitPbwt(pbwt_sym);

        //TWK_BCF_GT_UNPACK[(A >> 1)]

        uint8_t* out = new uint8_t[n_samples];
        uint32_t n_debug = 0;
        uint32_t debug_offset = 0;
        if (permute_pbwt) {
            for (int i = 0; i < n_cycles; ++i) {
                int32_t alts1 = CALL_MEMBER_FN(*debug1,decode_fn)(out);
                debug1->model->pbwt->ReverseUpdate(out);
                // assert(alts1 == debug1->pbwt->n_queue[1]);
                for (int j = 0; j < n_samples; ++j) {
                    assert(lookup_fn[ref_data[debug_offset + j] >> 1] == debug1->model->pbwt->prev[j]);
                }

                debug_offset += n_samples;
                ++n_debug;
            }
        } else {
            for (int i = 0; i < n_cycles; ++i) {
                int32_t alts1 = CALL_MEMBER_FN(*debug1,decode_fn)(out);
                for (int j = 0; j < n_samples; ++j) {
                    // assert(BCF_UNPACK_GENOTYPE_GENERAL(ref_data[debug_offset + j]) == out[j]);
                    assert(lookup_fn[ref_data[debug_offset + j] >> 1] == out[j]);
                }

                debug_offset += n_samples;
                ++n_debug;
            }
        }
        
        delete[] out;
        return 1;
    }
    return 1;
}
#endif

/*======   RLE-bitmap approach   ======*/

GenotypeCompressorRLEBitmap::GenotypeCompressorRLEBitmap(int64_t n_s) : GenotypeCompressor(n_s),
    bytes_out_zstd1(0), bytes_out_lz4(0)
{
    strategy = CompressionStrategy::LZ4;

    // gt_width.resize(n_s, 0);
    for (int i = 0; i < 6; ++i)
        buf_wah[i].resize(10000000);
        
    for (int i = 0; i < 6; ++i) n_variants[i] = 0;
    
    base_pbwt[0].Initiate(2*n_s, 2);
    base_pbwt[1].Initiate(n_s, 2);
    base_pbwt[2].Initiate(2*n_s, 4);
    base_pbwt[3].Initiate(n_s, 4);

    complex_pbwt[0].Initiate(2*n_s, 16);
    complex_pbwt[1].Initiate(n_s, 16);
}

GenotypeCompressorRLEBitmap::~GenotypeCompressorRLEBitmap() { }

int GenotypeCompressorRLEBitmap::Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles, const bool permute) {
    if (data == nullptr) return -1;
    if (n_data == 0) return -2;
    if (n_alleles == 0) return -3;

    // std::cerr << "n_alleles=" << n_alleles << std::endl;

    if (n_alleles == 2) return(Encode2N2M(data, n_data, permute)); // biallelic
    else if (n_alleles < 14) return(Encode2NXM(data, n_data, n_alleles, permute));  // alleles < 14 and 2 reserved for MISSING and EOV
    else {
        std::cerr << "raw alleles=" << n_alleles << std::endl;
        const uint8_t* gts = data;
        for (int i = 0; i < 2*n_samples; ++i) {
            buf_raw.data[buf_raw.len++] = gts[i];
        }

        ++processed_lines_local;
        ++processed_lines;
        return 1;
    }
    return -4;
}

int GenotypeCompressorRLEBitmap::Encode2N2M(uint8_t* data, const int32_t n_data, const bool permute) {
    // std::cerr << "2N2M" << std::endl;

    // Todo: assert genotypes are set for this variant.
    int replaced = RemapGenotypeEOV(data, n_data);
    if (replaced) {
        std::cerr << "replaced: skipping " << replaced << std::endl;
        return Encode2N2MM(data, n_data, permute); // 2N2MM supports missing+EOV in RLE-bitmap mode.
    }

    // for (int i = 0; i < 256; ++i) {
    //     if (alts[i]) std::cerr << i << ":" << alts[i] << ",";
    // }
    // std::cerr << std::endl;

    if (alts[15] == 0 && alts[14] == 0) { // No missing values and no EOV values
        return Encode2N2MC(data, n_data, permute);
    } else { // Having missing values.
        // std::cerr << "using extended model" << std::endl;
        return Encode2N2MM(data, n_data, permute);
    }

    return -1;
}

int GenotypeCompressorRLEBitmap::Encode2N2MC(uint8_t* data, const int32_t n_data, const bool permute) {
    // std::cerr << "2N2MC" << std::endl;
#if DEBUG_PBWT
    assert(debug_pbwt[0].UpdateDigestStride(&data[0], n_data, 2));
    assert(debug_pbwt[1].UpdateDigestStride(&data[1], n_data, 2));
#endif

    ++n_variants[0];
    ++n_variants[1];

#if DEBUG_WAH
    uint32_t buf_wah_before = buf_wah[0].len;
#endif

    if (permute_pbwt) {
        // Update diploid, biallelic, no-missing PBWT for haplotype 1 and 2.
        base_pbwt[0].UpdateBcf(&data[0], 1);
        // base_pbwt[1].UpdateBcf(&data[1], 2);

        EncodeRLEBitmap2N2MC(0);
    } else {
        EncodeRLEBitmap2N2MC(0, &data[0], 2);
    }

#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    std::shared_ptr<GenotypeDecompressorRLEBitmap> debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[0].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts1 = debug1->Decode2N2MC(buf_compress.data);
    if (permute_pbwt) {
        assert(alts1 == base_pbwt[0].n_queue[1]);
        // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 1.
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)base_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(base_pbwt[0].prev[i] == buf_compress.data[i]);
        }
    } else {
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)base_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(BCF_UNPACK_GENOTYPE(data[0+i*2]) == buf_compress.data[i]);
        }
    }
    buf_wah_before = buf_wah[1].len;
#endif

    // if (permute_pbwt) {
    //     EncodeRLEBitmap2N2MC(1);
    // } else {
    //     EncodeRLEBitmap2N2MC(1, &data[1], 2);
    // }

#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[1].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts2 = debug1->Decode2N2MC(buf_compress.data);
    if (permute_pbwt) {
        assert(alts2 == base_pbwt[1].n_queue[1]);
        // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 2.
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)base_pbwt[1].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(base_pbwt[1].prev[i] == buf_compress.data[i]);
        }
    } else {
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)base_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(BCF_UNPACK_GENOTYPE(data[1+i*2]) == buf_compress.data[i]);
        }
    }
#endif

    ++processed_lines_local;
    ++processed_lines;
    return 1;
}

int GenotypeCompressorRLEBitmap::Encode2N2MM(uint8_t* data, const int32_t n_data, const bool permute) {
#if DEBUG_PBWT
    assert(debug_pbwt[2].UpdateDigestStride(&data[0], n_data, 2));
    assert(debug_pbwt[3].UpdateDigestStride(&data[1], n_data, 2));
#endif

    ++n_variants[2];
    ++n_variants[3];

    // std::cerr << "2N2MM" << std::endl;

#if DEBUG_WAH
        uint32_t buf_wah_before = buf_wah[2].len;
#endif

    if (permute_pbwt) {
        // Update diploid, biallelic, no-missing PBWT for haplotype 1 and 2.
        base_pbwt[2].UpdateBcf(&data[0], 1);
        // base_pbwt[3].UpdateBcf(&data[1], 2);

        int ret1 = EncodeRLEBitmap2N2MM(2);
    } else {
        int ret1 = EncodeRLEBitmap2N2MM(2, &data[0], 2);
    }

#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    std::shared_ptr<GenotypeDecompressorRLEBitmap> debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[2].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts1 = debug1->Decode2N2MM(buf_compress.data);
    // assert(alts1 == complex_pbwt[0].n_queue[1]);

    if (permute_pbwt) {
        // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 1.
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)base_pbwt[2].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(base_pbwt[2].prev[i] == buf_compress.data[i]);
        }
    } else {
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)base_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(BCF_UNPACK_GENOTYPE(data[0+i*2]) == buf_compress.data[i]);
        }
    }
    buf_wah_before = buf_wah[3].len;
#endif
    
    // if (permute_pbwt) {
    //     int ret2 = EncodeRLEBitmap2N2MM(3);
    // } else {
    //     int ret2 = EncodeRLEBitmap2N2MM(3, &data[1], 2);
    // }

#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[3].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts2 = debug1->Decode2N2MM(buf_compress.data);
    // assert(alts1 == complex_pbwt[0].n_queue[1]);

    if (permute_pbwt) {
        // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 1.
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)base_pbwt[2].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(base_pbwt[3].prev[i] == buf_compress.data[i]);
        }
    } else {
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)base_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(BCF_UNPACK_GENOTYPE(data[1+i*2]) == buf_compress.data[i]);
        }
    }
#endif

    ++processed_lines_local;
    ++processed_lines;

    return 0;
}

int GenotypeCompressorRLEBitmap::Encode2NXM(uint8_t* data, const int32_t n_data, const int32_t n_alleles, const bool permute) { 
#if DEBUG_PBWT
    assert(debug_pbwt[4].UpdateDigestStride(&data[0], n_data, 2));
    assert(debug_pbwt[5].UpdateDigestStride(&data[1], n_data, 2));
#endif

    ++n_variants[4];
    ++n_variants[5];

    // Todo: check for missing and EOV and ascertain total < 16
    if (n_alleles >= 14) {
        // check for missing and EOV
    }

    // memset(alts, 0, 256*sizeof(uint32_t));
    
    // for (int i = 0; i < 2*n_samples; ++i) {
    //     uint8_t val = BCF_UNPACK_GENOTYPE_GENERAL(data[i]);
    //     ++alts[val];
    // }

    // for (int i = 0; i < 256; ++i) {
    //     if (alts[i]) std::cerr << i << ":" << alts[i] << ",";
    // }
    // std::cerr << std::endl;

#if DEBUG_WAH
    uint32_t buf_wah_before = buf_wah[4].len;
#endif

    if (permute_pbwt) {
        complex_pbwt[0].UpdateBcfGeneral(&data[0], 1);
        // complex_pbwt[1].UpdateBcfGeneral(&data[1], 2);
        
        EncodeRLEBitmap2NXM(0);
    } else {
        EncodeRLEBitmap2NXM(0, &data[0], 2);
    }

#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    std::shared_ptr<GenotypeDecompressorRLEBitmap> debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[4].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts1 = debug1->Decode2NXM(buf_compress.data);
    // assert(alts1 == complex_pbwt[0].n_queue[1]);

    if (permute_pbwt) {
        // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 1.
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)complex_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(complex_pbwt[0].prev[i] == buf_compress.data[i]);
        }
    } else {
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)base_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(BCF_UNPACK_GENOTYPE_GENERAL(data[0+i*2]) == buf_compress.data[i]);
        }
    }
    buf_wah_before = buf_wah[5].len;
#endif

    // if (permute_pbwt) {
    //     EncodeRLEBitmap2NXM(1);
    // } else {
    //     EncodeRLEBitmap2NXM(1, &data[1], 2);
    // }

#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[5].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts2 = debug1->Decode2NXM(buf_compress.data);
    // assert(alts1 == complex_pbwt[0].n_queue[1]);

    if (permute_pbwt) {
        // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 1.
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)complex_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(complex_pbwt[1].prev[i] == buf_compress.data[i]);
        }
    } else {
        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << i << " -> " << (int)base_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
            assert(BCF_UNPACK_GENOTYPE_GENERAL(data[1+i*2]) == buf_compress.data[i]);
        }
    }
#endif

    ++processed_lines_local;
    ++processed_lines;

    return 1;
}

int GenotypeCompressorRLEBitmap::EncodeRLEBitmap2N2MC(const int target) {
    assert(target == 0 || target == 1);

    // RLE word -> 1 bit typing, 1 bit allele, word*8-2 bits for RLE
    const uint8_t* prev = base_pbwt[target].prev;
    uint8_t* dst = buf_wah[target].data;
    size_t& dst_len = buf_wah[target].len;

    // Setup.
    uint8_t ref = prev[0];
    uint32_t rle_len = 1;
    uint32_t start_rle = 0, end_rle = 1;

    // Debug.
    uint32_t rle_cost = 0;
    uint32_t n_objects = 0;
    uint32_t observed_alts = 0;
    uint32_t observed_length = 0;

    for (int i = 1; i < 2*n_samples; ++i) {
        if (prev[i] != ref || rle_len == 16383) { // 2^14-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 31 then the word is considered "dirty"
            // and a 32-bit bitmap will be used instead.
            if (end_rle - start_rle < 31) {
                // Increment position
                i += 31 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 31 < 2*n_samples ? start_rle + 31 : 2*n_samples;
                // Update observed path.
                observed_length += start_rle + 31 < 2*n_samples ? 31 : 2*n_samples - start_rle;
                // Assertion.
                assert(end_rle <= 2*n_samples);
                
                // Construct 32-bit bitmap
                uint32_t bitmap = 0;
                for (int j = start_rle; j < end_rle; ++j) {
                    bitmap |= (prev[j] & 1);
                    bitmap <<= 1;
                }
                assert((bitmap & 1) == 0);
                observed_alts += __builtin_popcount(bitmap);
                
                // Debug:
                // std::cerr << "Use Bitmap=" << (end_rle-start_rle) << "(" << start_rle << "," << end_rle << ") " << std::bitset<32>(bitmap) << std::endl;
                
                memcpy(&dst[dst_len], &bitmap, sizeof(uint32_t));
                dst_len += sizeof(uint32_t);

                start_rle = i;
                // Set new reference.
                ref = prev[i];
                // Set run-length. If this is the final object then set to 0
                // for downstream filter.
                rle_len = (end_rle == 2*n_samples) ? 0 : 1;
                // Update cost.
                rle_cost += sizeof(uint32_t);
                continue;
                
            } else {
                // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
                
                // Model:
                // Lowest bit => 1 iff RLE, 0 otherwise
                // If RLE => bit 2 is the reference
                // Remainder => dst (bitmap / run-length)
                uint16_t rle_pack = (rle_len << 2) | ((ref & 1) << 1) | 1;
                memcpy(&dst[dst_len], &rle_pack, sizeof(uint16_t));
                dst_len += sizeof(uint16_t);
                
                // Update observed alts.
                observed_alts += ((rle_pack >> 1) & 1) * rle_len; // predicate multiply
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
        // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
        ++n_objects;

        uint16_t rle_pack = (rle_len << 2) | ((ref & 1) << 1) | 1;
        memcpy(&dst[dst_len], &rle_pack, sizeof(uint16_t));
        dst_len += sizeof(uint16_t);

        rle_cost += sizeof(uint16_t);
        observed_length += rle_len;
    }
    observed_alts += (ref == 1) * rle_len;

    // Ascertain correctness:
    assert(observed_alts == base_pbwt[target].n_queue[1]);
    assert(observed_length == 2*n_samples);

    // ++gt_width[n_objects];

    return n_objects;
}

int GenotypeCompressorRLEBitmap::EncodeRLEBitmap2N2MC(const int target, const uint8_t* data, const int stride) {
    assert(target == 0 || target == 1);

    // RLE word -> 1 bit typing, 1 bit allele, word*8-2 bits for RLE
    uint8_t* dst = buf_wah[target].data;
    size_t& dst_len = buf_wah[target].len;

    // Setup.
    uint8_t ref = BCF_UNPACK_GENOTYPE(data[0]);
    uint32_t rle_len = 1;
    uint32_t start_rle = 0, end_rle = 1;

    // Debug.
    uint32_t rle_cost = 0;
    uint32_t n_objects = 0;
    uint32_t observed_alts = 0;
    uint32_t observed_length = 0;

    for (int i = 1; i < 2*n_samples; ++i) {
        if (BCF_UNPACK_GENOTYPE(data[i*stride]) != ref || rle_len == 16383) { // 2^14-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 31 then the word is considered "dirty"
            // and a 32-bit bitmap will be used instead.
            if (end_rle - start_rle < 31) {
                // Increment position
                i += 31 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 31 < 2*n_samples ? start_rle + 31 : 2*n_samples;
                // Update observed path.
                observed_length += start_rle + 31 < 2*n_samples ? 31 : 2*n_samples - start_rle;
                // Assertion.
                assert(end_rle <= 2*n_samples);
                
                // Construct 32-bit bitmap
                uint32_t bitmap = 0;
                for (int j = start_rle; j < end_rle; ++j) {
                    bitmap |= (BCF_UNPACK_GENOTYPE(data[j*stride]) & 1);
                    bitmap <<= 1;
                }
                assert((bitmap & 1) == 0);
                observed_alts += __builtin_popcount(bitmap);
                
                // Debug:
                // std::cerr << "Use Bitmap=" << (end_rle-start_rle) << "(" << start_rle << "," << end_rle << ") " << std::bitset<32>(bitmap) << std::endl;
                
                memcpy(&dst[dst_len], &bitmap, sizeof(uint32_t));
                dst_len += sizeof(uint32_t);

                start_rle = i;
                // Set new reference.
                ref = BCF_UNPACK_GENOTYPE(data[i*stride]);
                // Set run-length. If this is the final object then set to 0
                // for downstream filter.
                rle_len = (end_rle == 2*n_samples) ? 0 : 1;
                // Update cost.
                rle_cost += sizeof(uint32_t);
                continue;
                
            } else {
                // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
                
                // Model:
                // Lowest bit => 1 iff RLE, 0 otherwise
                // If RLE => bit 2 is the reference
                // Remainder => dst (bitmap / run-length)
                uint16_t rle_pack = (rle_len << 2) | ((ref & 1) << 1) | 1;
                memcpy(&dst[dst_len], &rle_pack, sizeof(uint16_t));
                dst_len += sizeof(uint16_t);
                
                // Update observed alts.
                observed_alts += ((rle_pack >> 1) & 1) * rle_len; // predicate multiply
                // Update observed path.
                observed_length += rle_len;

                // Reset length;
                rle_len = 0;
                // Set new reference.
                ref = BCF_UNPACK_GENOTYPE(data[i*stride]);
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
        // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
        ++n_objects;

        uint16_t rle_pack = (rle_len << 2) | ((ref & 1) << 1) | 1;
        memcpy(&dst[dst_len], &rle_pack, sizeof(uint16_t));
        dst_len += sizeof(uint16_t);

        rle_cost += sizeof(uint16_t);
        observed_length += rle_len;
    }
    observed_alts += (ref == 1) * rle_len;

    // Ascertain correctness:
    assert(observed_length == 2*n_samples);

    // ++gt_width[n_objects];

    return n_objects;
}

int GenotypeCompressorRLEBitmap::EncodeRLEBitmap2N2MM(const int target) {
    assert(target == 2 || target == 3);

    // RLE word -> 1 bit typing, 2 bit allele, word*8-3 bits for RLE
    const uint8_t* prev = base_pbwt[target].prev;
    uint8_t* dst = buf_wah[target].data;
    size_t& dst_len = buf_wah[target].len;

    // Setup.
    uint8_t ref = prev[0];
    uint32_t rle_len = 1;
    uint32_t start_rle = 0, end_rle = 1;

    // Debug.
    uint32_t rle_cost = 0;
    uint32_t n_objects = 0;
    // uint32_t observed_alts = 0;
    uint32_t observed_length = 0;

    for (int i = 1; i < 2*n_samples; ++i) {
        if (prev[i] != ref || rle_len == 8191) { // 2^13-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 15 then the word is considered "dirty"
            // and a 32-bit bitmap will be used instead.
            if (end_rle - start_rle < 15) { // 15*2 = 30 and 1 bit for class type = 31 bits
                // Increment position
                i += 15 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 15 < 2*n_samples ? start_rle + 15 : 2*n_samples;
                // Update observed path.
                observed_length += start_rle + 15 < 2*n_samples ? 15 : 2*n_samples - start_rle;
                // Assertion.
                assert(end_rle <= 2*n_samples);
                
                // Construct 32-bit bitmap
                uint32_t bitmap = 0;
                int j = start_rle;
                for (/**/; j < end_rle - 1; ++j) {
                    bitmap |= (prev[j] & 3);
                    bitmap <<= 2;
                }
                for (/**/; j < end_rle; ++j) {
                    bitmap |= (prev[j] & 3);
                }
                bitmap <<= 1;
                
                assert((bitmap & 1) == 0);
                // observed_alts += __builtin_popcount(bitmap);
                
                // Debug:
                // std::cerr << "Use Bitmap=" << (end_rle-start_rle) << "(" << start_rle << "," << end_rle << ") " << std::bitset<32>(bitmap) << std::endl;
                
                memcpy(&dst[dst_len], &bitmap, sizeof(uint32_t));
                dst_len += sizeof(uint32_t);

                start_rle = i;
                // Set new reference.
                ref = prev[i];
                // Set run-length. If this is the final object then set to 0
                // for downstream filter.
                rle_len = (end_rle == 2*n_samples) ? 0 : 1;
                // Update cost.
                rle_cost += sizeof(uint32_t);
                continue;
                
            } else {
                // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
                
                // Model:
                // Lowest bit => 1 iff RLE, 0 otherwise
                // If RLE => bit 2 is the reference
                // Remainder => dst (bitmap / run-length)
                uint16_t rle_pack = (rle_len << 3) | ((ref & 3) << 1) | 1;
                memcpy(&dst[dst_len], &rle_pack, sizeof(uint16_t));
                dst_len += sizeof(uint16_t);
                
                // Update observed alts.
                // observed_alts += ((rle_pack >> 1) & 3) * rle_len; // predicate multiply
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
        // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
        ++n_objects;

        uint16_t rle_pack = (rle_len << 3) | ((ref & 3) << 1) | 1;
        memcpy(&dst[dst_len], &rle_pack, sizeof(uint16_t));
        dst_len += sizeof(uint16_t);

        rle_cost += sizeof(uint16_t);
        observed_length += rle_len;
    }
    // observed_alts += (ref == 1) * rle_len;

    // Ascertain correctness:
    // assert(observed_alts == base_pbwt[target].n_queue[1]);
    assert(observed_length == 2*n_samples);

    // ++gt_width[n_objects];

    return n_objects;
}

int GenotypeCompressorRLEBitmap::EncodeRLEBitmap2N2MM(const int target, const uint8_t* data, const int stride) {
    assert(target == 2 || target == 3);

    // RLE word -> 1 bit typing, 2 bit allele, word*8-3 bits for RLE
    uint8_t* dst = buf_wah[target].data;
    size_t& dst_len = buf_wah[target].len;

    // Setup.
    uint8_t ref = BCF_UNPACK_GENOTYPE(data[0]);
    uint32_t rle_len = 1;
    uint32_t start_rle = 0, end_rle = 1;

    // Debug.
    uint32_t rle_cost = 0;
    uint32_t n_objects = 0;
    // uint32_t observed_alts = 0;
    uint32_t observed_length = 0;

    for (int i = 1; i < 2*n_samples; ++i) {
        if (BCF_UNPACK_GENOTYPE(data[i*stride]) != ref || rle_len == 8191) { // 2^13-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 15 then the word is considered "dirty"
            // and a 32-bit bitmap will be used instead.
            if (end_rle - start_rle < 15) { // 15*2 = 30 and 1 bit for class type = 31 bits
                // Increment position
                i += 15 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 15 < 2*n_samples ? start_rle + 15 : 2*n_samples;
                // Update observed path.
                observed_length += start_rle + 15 < 2*n_samples ? 15 : 2*n_samples - start_rle;
                // Assertion.
                assert(end_rle <= 2*n_samples);
                
                // Construct 32-bit bitmap
                uint32_t bitmap = 0;
                int j = start_rle;
                for (/**/; j < end_rle - 1; ++j) {
                    bitmap |= (BCF_UNPACK_GENOTYPE(data[j*stride]) & 3);
                    bitmap <<= 2;
                }
                for (/**/; j < end_rle; ++j) {
                    bitmap |= (BCF_UNPACK_GENOTYPE(data[j*stride]) & 3);
                }
                bitmap <<= 1;
                
                assert((bitmap & 1) == 0);
                // observed_alts += __builtin_popcount(bitmap);
                
                // Debug:
                // std::cerr << "Use Bitmap=" << (end_rle-start_rle) << "(" << start_rle << "," << end_rle << ") " << std::bitset<32>(bitmap) << std::endl;
                
                memcpy(&dst[dst_len], &bitmap, sizeof(uint32_t));
                dst_len += sizeof(uint32_t);

                start_rle = i;
                // Set new reference.
                ref = BCF_UNPACK_GENOTYPE(data[i*stride]);
                // Set run-length. If this is the final object then set to 0
                // for downstream filter.
                rle_len = (end_rle == 2*n_samples) ? 0 : 1;
                // Update cost.
                rle_cost += sizeof(uint32_t);
                continue;
                
            } else {
                // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
                
                // Model:
                // Lowest bit => 1 iff RLE, 0 otherwise
                // If RLE => bit 2 is the reference
                // Remainder => dst (bitmap / run-length)
                uint16_t rle_pack = (rle_len << 3) | ((ref & 3) << 1) | 1;
                memcpy(&dst[dst_len], &rle_pack, sizeof(uint16_t));
                dst_len += sizeof(uint16_t);
                
                // Update observed alts.
                // observed_alts += ((rle_pack >> 1) & 3) * rle_len; // predicate multiply
                // Update observed path.
                observed_length += rle_len;

                // Reset length;
                rle_len = 0;
                // Set new reference.
                ref = BCF_UNPACK_GENOTYPE(data[i*stride]);
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
        // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
        ++n_objects;

        uint16_t rle_pack = (rle_len << 3) | ((ref & 3) << 1) | 1;
        memcpy(&dst[dst_len], &rle_pack, sizeof(uint16_t));
        dst_len += sizeof(uint16_t);

        rle_cost += sizeof(uint16_t);
        observed_length += rle_len;
    }
    // observed_alts += (ref == 1) * rle_len;

    // Ascertain correctness:
    // assert(observed_alts == base_pbwt[target].n_queue[1]);
    assert(observed_length == 2*n_samples);

    // ++gt_width[n_objects];

    return n_objects;
}

int GenotypeCompressorRLEBitmap::EncodeRLEBitmap2NXM(const int target) {
    assert(target == 0 || target == 1);

    // RLE word -> 1 bit typing, 4 bit allele, word*8-5 bits for RLE
    const uint8_t* prev = complex_pbwt[target].prev;
    uint8_t* dst = buf_wah[4+target].data;
    size_t& dst_len = buf_wah[4+target].len;

    // Setup.
    uint8_t ref = prev[0];
    uint32_t rle_len = 1;
    uint32_t start_rle = 0, end_rle = 1;

    // Debug.
    uint32_t rle_cost = 0;
    uint32_t n_objects = 0;
    uint32_t observed_alts = 0;
    uint32_t observed_length = 0;

    for (int i = 1; i < 2*n_samples; ++i) {
        if (prev[i] != ref || rle_len == 134217727) { // 2^27-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 15 then the word is considered "dirty"
            // and a 64-bit bitmap will be used instead.
            if (end_rle - start_rle < 15) { // 15*4 = 60 and 1 bit for class type = 31 bits
                // Increment position
                i += 15 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 15 < 2*n_samples ? start_rle + 15 : 2*n_samples;
                // Update observed path.
                observed_length += start_rle + 15 < 2*n_samples ? 15 : 2*n_samples - start_rle;
                // Assertion.
                assert(end_rle <= 2*n_samples);
                
                // Construct 32-bit bitmap
                uint64_t bitmap = 0;
                uint32_t j = start_rle;
                for (/**/; j < end_rle - 1; ++j) {
                    bitmap |= (prev[j] & 15);
                    bitmap <<= 4;
                }
                for (/**/; j < end_rle; ++j) {
                    bitmap |= (prev[j] & 15);
                }
                bitmap <<= 1;

                assert((bitmap & 1) == 0);
                // observed_alts += __builtin_popcount(bitmap);
                
                // Debug:
                // std::cerr << "Use Bitmap=" << (end_rle-start_rle) << "(" << start_rle << "," << end_rle << ") " << std::bitset<32>(bitmap) << std::endl;
                
                memcpy(&dst[dst_len], &bitmap, sizeof(uint64_t));
                dst_len += sizeof(uint64_t);

                start_rle = i;
                // Set new reference.
                ref = prev[i];
                // Set run-length. If this is the final object then set to 0
                // for downstream filter.
                rle_len = (end_rle == 2*n_samples) ? 0 : 1;
                // Update cost.
                rle_cost += sizeof(uint64_t);
                continue;
                
            } else {
                // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
                
                // Model:
                // Lowest bit => 1 iff RLE, 0 otherwise
                // If RLE => bit 2 is the reference
                // Remainder => dst (bitmap / run-length)
                uint32_t rle_pack = (rle_len << 5) | ((ref & 15) << 1) | 1;
                memcpy(&dst[dst_len], &rle_pack, sizeof(uint32_t));
                dst_len += sizeof(uint32_t);
                
                // Update observed alts.
                // observed_alts += ((rle_pack >> 1) & 3) * rle_len; // predicate multiply
                
                // Update observed path.
                observed_length += rle_len;

                // Reset length;
                rle_len = 0;
                // Set new reference.
                ref = prev[i];
                // Update cost.
                rle_cost += sizeof(uint32_t);
                // Update start position.
                start_rle = i;
            }
        }
        ++rle_len; 
    }

    // Add only if last element is an unfinished RLE.
    if (rle_len) {
        // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
        ++n_objects;

        uint32_t rle_pack = (rle_len << 5) | ((ref & 15) << 1) | 1;
        memcpy(&dst[dst_len], &rle_pack, sizeof(uint32_t));
        dst_len += sizeof(uint32_t);

        rle_cost += sizeof(uint32_t);
        observed_length += rle_len;
    }
    // observed_alts += (ref == 1) * rle_len;

    // Ascertain correctness:
    // assert(observed_alts == complex_pbwt[target].n_queue[1]);
    assert(observed_length == 2*n_samples);

    // ++gt_width[n_objects];

    return n_objects;
}

int GenotypeCompressorRLEBitmap::EncodeRLEBitmap2NXM(const int target, const uint8_t* data, const int stride) {
    assert(target == 0 || target == 1);

    // RLE word -> 1 bit typing, 4 bit allele, word*8-5 bits for RLE
    uint8_t* dst = buf_wah[4+target].data;
    size_t& dst_len = buf_wah[4+target].len;

    // Setup.
    uint8_t ref = BCF_UNPACK_GENOTYPE_GENERAL(data[0]);
    uint32_t rle_len = 1;
    uint32_t start_rle = 0, end_rle = 1;

    // Debug.
    uint32_t rle_cost = 0;
    uint32_t n_objects = 0;
    uint32_t observed_alts = 0;
    uint32_t observed_length = 0;

    for (int i = 1; i < 2*n_samples; ++i) {
        if (BCF_UNPACK_GENOTYPE_GENERAL(data[i*stride]) != ref || rle_len == 134217727) { // 2^27-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 15 then the word is considered "dirty"
            // and a 64-bit bitmap will be used instead.
            if (end_rle - start_rle < 15) { // 15*4 = 60 and 1 bit for class type = 31 bits
                // Increment position
                i += 15 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 15 < 2*n_samples ? start_rle + 15 : 2*n_samples;
                // Update observed path.
                observed_length += start_rle + 15 < 2*n_samples ? 15 : 2*n_samples - start_rle;
                // Assertion.
                assert(end_rle <= 2*n_samples);
                
                // Construct 32-bit bitmap
                uint64_t bitmap = 0;
                uint32_t j = start_rle;
                for (/**/; j < end_rle - 1; ++j) {
                    bitmap |= (BCF_UNPACK_GENOTYPE_GENERAL(data[j*stride]) & 15);
                    bitmap <<= 4;
                }
                for (/**/; j < end_rle; ++j) {
                    bitmap |= (BCF_UNPACK_GENOTYPE_GENERAL(data[j*stride]) & 15);
                }
                bitmap <<= 1;

                assert((bitmap & 1) == 0);
                // observed_alts += __builtin_popcount(bitmap);
                
                // Debug:
                // std::cerr << "Use Bitmap=" << (end_rle-start_rle) << "(" << start_rle << "," << end_rle << ") " << std::bitset<32>(bitmap) << std::endl;
                
                memcpy(&dst[dst_len], &bitmap, sizeof(uint64_t));
                dst_len += sizeof(uint64_t);

                start_rle = i;
                // Set new reference.
                ref = BCF_UNPACK_GENOTYPE_GENERAL(data[i*stride]);
                // Set run-length. If this is the final object then set to 0
                // for downstream filter.
                rle_len = (end_rle == 2*n_samples) ? 0 : 1;
                // Update cost.
                rle_cost += sizeof(uint64_t);
                continue;
                
            } else {
                // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
                
                // Model:
                // Lowest bit => 1 iff RLE, 0 otherwise
                // If RLE => bit 2 is the reference
                // Remainder => dst (bitmap / run-length)
                uint32_t rle_pack = (rle_len << 5) | ((ref & 15) << 1) | 1;
                memcpy(&dst[dst_len], &rle_pack, sizeof(uint32_t));
                dst_len += sizeof(uint32_t);
                
                // Update observed alts.
                // observed_alts += ((rle_pack >> 1) & 3) * rle_len; // predicate multiply
                
                // Update observed path.
                observed_length += rle_len;

                // Reset length;
                rle_len = 0;
                // Set new reference.
                ref = BCF_UNPACK_GENOTYPE_GENERAL(data[i*stride]);
                // Update cost.
                rle_cost += sizeof(uint32_t);
                // Update start position.
                start_rle = i;
            }
        }
        ++rle_len; 
    }

    // Add only if last element is an unfinished RLE.
    if (rle_len) {
        // std::cerr << "Use RLE=" << rle_len << ":" << (int)ref << "(" << start_rle << "," << end_rle << ")" << std::endl;
        ++n_objects;

        uint32_t rle_pack = (rle_len << 5) | ((ref & 15) << 1) | 1;
        memcpy(&dst[dst_len], &rle_pack, sizeof(uint32_t));
        dst_len += sizeof(uint32_t);

        rle_cost += sizeof(uint32_t);
        observed_length += rle_len;
    }
    // observed_alts += (ref == 1) * rle_len;

    // Ascertain correctness:
    // assert(observed_alts == complex_pbwt[target].n_queue[1]);
    assert(observed_length == 2*n_samples);

    // ++gt_width[n_objects];

    return n_objects;
}

#if DEBUG_WAH
int GenotypeCompressorRLEBitmap::DebugWAH(uint8_t* in, size_t len_in, uint8_t* ref_data, size_t n_cycles, int pbwt_sym, bitmap_debug_decode decode_fn, const uint8_t* lookup_fn) {
    if (len_in) {
        std::shared_ptr<GenotypeDecompressorRLEBitmap> debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(in, len_in, n_cycles, n_samples);
        if (permute_pbwt)
            debug1->InitPbwt(pbwt_sym);

        //TWK_BCF_GT_UNPACK[(A >> 1)]

        uint8_t* out = new uint8_t[n_samples];
        uint32_t n_debug = 0;
        uint32_t debug_offset = 0;
        if (permute_pbwt) {
            for (int i = 0; i < n_cycles; ++i) {
                int32_t alts1 = CALL_MEMBER_FN(*debug1,decode_fn)(out,1);
                debug1->pbwt->ReverseUpdate(out);
                // assert(alts1 == debug1->pbwt->n_queue[1]);
                for (int j = 0; j < n_samples; ++j) {
                    assert(lookup_fn[ref_data[debug_offset + j] >> 1] == debug1->pbwt->prev[j]);
                }

                debug_offset += n_samples;
                ++n_debug;
            }
        } else {
            for (int i = 0; i < n_cycles; ++i) {
                int32_t alts1 = CALL_MEMBER_FN(*debug1,decode_fn)(out,1);
                for (int j = 0; j < n_samples; ++j) {
                    // assert(BCF_UNPACK_GENOTYPE_GENERAL(ref_data[debug_offset + j]) == out[j]);
                    assert(lookup_fn[ref_data[debug_offset + j] >> 1] == out[j]);
                }

                debug_offset += n_samples;
                ++n_debug;
            }
        }
        
        delete[] out;
        return 1;
    }
    return 1;
}
#endif
}