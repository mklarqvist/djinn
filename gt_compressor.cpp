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
        std::cerr << "input is not divisible by 2" << std::endl;
        return 0;
    }

    // Keep track of input bytes.
    bytes_in += fmt->p_len;
    bytes_in_vcf += fmt->p_len * 2 - 1; // (char)(sep)(char)(tab)

    // if (bcf->n_allele < 16 && bcf->n_allele > 2) {
    //     for (int i = 0; i < bcf->n_allele; ++i) {
    //         std::cerr << bcf->d.allele[i] << ",";
    //     }
    //     std::cerr << std::endl;
    // }

    // Todo: extend beyond 2N
    return(Encode2N(fmt->p, fmt->p_len, bcf->n_allele));
}

int GenotypeCompressor::Encode(uint8_t* data, const int32_t n_data, const int32_t ploidy, const int32_t n_alleles) {
    if (data == nullptr) return 0;
    if (n_data == 0) return 0;
    if (n_alleles <= 0) return 0;

    if (ploidy != 2) {
        std::cerr << "only diploid data is currently supported" << std::endl;
        return 0;
    }

    bytes_in += n_data;
    bytes_in_vcf += n_data * 2 - 1;

    return(Encode2N(data, n_data, n_alleles));
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

GenotypeCompressorModelling::GenotypeCompressorModelling(int64_t n_s) : GenotypeCompressor(n_s),
    models(nullptr)
{
    strategy = CompressionStrategy::CONTEXT;
    base_models[0].Construct(n_samples*2, 2);
    base_models[1].Construct(n_samples, 2);
    base_models[2].Construct(n_samples, 4);
    base_models[3].Construct(n_samples, 4);
    base_models_complex[0].Construct(n_samples, 16);
    base_models_complex[1].Construct(n_samples, 16);

    base_models[0].StartEncoding();
    base_models[1].StartEncoding();
    base_models[2].StartEncoding();
    base_models[3].StartEncoding();
    base_models_complex[0].StartEncoding();
    base_models_complex[1].StartEncoding();

    base_model_bitmaps[0].Construct(1, 2);
    base_model_bitmaps[1].Construct(1, 2);
    base_model_bitmaps[0].StartEncoding();
    base_model_bitmaps[1].StartEncoding();

    pack1_context = 0; pack2_context = 0;
    ppm1.encoder.buffer = new uint8_t[10000000];
    ppm1.encoder.dat = ppm1.encoder.buffer;
    ppm2.encoder.buffer = new uint8_t[10000000];
    ppm2.encoder.dat = ppm2.encoder.buffer;

    pack1bin_context = 0; pack2bin_context = 0;
    ppm_bin1.encoder.buffer = new uint8_t[10000000];
    ppm_bin1.encoder.dat = ppm_bin1.encoder.buffer;
    ppm_bin2.encoder.buffer = new uint8_t[10000000];
    ppm_bin2.encoder.dat = ppm_bin2.encoder.buffer;

    bytes_out2 = 0; bytes_out3 = 0;

#if DEBUG_PBWT
    for (int i = 0; i < 2; ++i) {
        debug_bins[i].resize(n_samples*block_size + 65536);
    }
#endif

    mref = std::make_shared<GeneralModel>(2, 1024);
    mlog_rle = std::make_shared<GeneralModel>(32, 32768); // 2^(4+1)
    mrle = std::make_shared<GeneralModel>(256, 32, 24, 1); // 2^5
    mrle2_1 = std::make_shared<GeneralModel>(256, 32);
    mrle2_2 = std::make_shared<GeneralModel>(256, 32);
    mrle4_1 = std::make_shared<GeneralModel>(256, 32);
    mrle4_2 = std::make_shared<GeneralModel>(256, 32);
    mrle4_3 = std::make_shared<GeneralModel>(256, 32);
    mrle4_4 = std::make_shared<GeneralModel>(256, 32);
    dirty_wah = std::make_shared<GeneralModel>(256, 65536, 16, 8);
    mtype = std::make_shared<GeneralModel>(2, 512);

    mref->StartEncoding();
    mlog_rle->StartEncoding();
    mrle->StartEncoding();
    mrle2_1->StartEncoding();
    mrle2_2->StartEncoding();
    mrle4_1->StartEncoding();
    mrle4_2->StartEncoding();
    mrle4_3->StartEncoding();
    mrle4_4->StartEncoding();
    dirty_wah->StartEncoding();
    mtype->StartEncoding();
    bytes_out4 = 0;

    // gtshark
    ctx_model = std::make_shared<djinn_gt_ctx>();
    ctx_model->data = std::make_shared<djinn_ctx_buf_t>(10000000);
    ctx_model->rce = std::make_shared<CRangeEncoder>();
    ctx_model->rce->Start();
    out_gts = 0;
}

GenotypeCompressorModelling::~GenotypeCompressorModelling() { delete[] models;}

bool GenotypeCompressorModelling::CheckLimit() const {
    return (processed_lines_local == block_size || 
            base_models[0].range_coder->OutSize() > 9000000 || 
            base_models[1].range_coder->OutSize() > 9000000 || 
            base_models[2].range_coder->OutSize() > 9000000 || 
            base_models[3].range_coder->OutSize() > 9000000 ||
            base_models_complex[0].range_coder->OutSize() > 9000000 ||
            base_models_complex[1].range_coder->OutSize() > 9000000 ||
            buf_raw.len > 9000000);
}

int GenotypeCompressorModelling::Encode2N(bcf1_t* bcf, const bcf_hdr_t* hdr) {
    if (bcf == NULL) return 0;
    if (hdr == NULL) return 0;
    
    const bcf_fmt_t* fmt = bcf_get_fmt(hdr, bcf, "GT");
    if (fmt == NULL) return 0;
    if (fmt->p_len / n_samples != 2) {
        std::cerr << "input is not divisible by 2" << std::endl;
        return 0;
    }

    if (bcf->n_allele == 2) return(Encode2N2M(fmt->p, fmt->p_len)); // biallelic
    else if (bcf->n_allele < 14) {

        int replaced = RemapGenotypeEOV(fmt->p, fmt->p_len);
        if (replaced) {
            // std::cerr << "replaced: skipping " << replaced << std::endl;
            return Encode2N2MM(fmt->p, fmt->p_len); // 2N2MM supports missing+EOV in RLE-bitmap mode.
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

int GenotypeCompressorModelling::Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles) {
    if (n_alleles == 2) return(Encode2N2M(data, n_data)); // biallelic
    else if (n_alleles < 14) return(Encode2NXM(data, n_data, n_alleles));  // #alleles < 14
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

int GenotypeCompressorModelling::EncodeRLE(uint64_t ref, uint32_t len) {
    mref->EncodeSymbol(ref&1);
    uint32_t log_length = ilog2(len);
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

// Wrapper for 2N2M
int GenotypeCompressorModelling::Encode2N2M(uint8_t* data, const int32_t n_data) {
    // Todo: assert genotypes are set for this variant.
    int replaced = RemapGenotypeEOV(data, n_data);
    if (replaced) {
        std::cerr << "replaced: divert to missing with " << replaced << std::endl;
        return Encode2NXM(data ,n_data, 4); // todo: see if this is true
    }

    // for (int i = 0; i < 256; ++i) {
        // if (alts[i]) std::cerr << i << ":" << alts[i] << ",";
    // }
    // std::cerr << std::endl;

    if (alts[14] == 0 && alts[15] == 0) { // No missing values.
        // std::cerr << "divert 2N2MC" << std::endl;
        return Encode2N2MC(data, n_data);
    } else { // Having missing values.
        // std::cerr << "using extended model" << std::endl;
        // std::cerr << "divert 2N2MM" << std::endl;
        return Encode2N2MM(data, n_data);
    }

    return 1;
}

// 2N2M complete
int GenotypeCompressorModelling::Encode2N2MC(uint8_t* data, const int32_t n_data) {
#if DEBUG_PBWT
    assert(debug_pbwt[0].UpdateDigestStride(&data[0], n_data, 2));
    assert(debug_pbwt[1].UpdateDigestStride(&data[1], n_data, 2));
#endif

    ++base_models[0];
    // ++base_models[1];
    
    if (permute_pbwt) {
        
        base_models[0].pbwt->Update(data, 1);
        // base_models[1].pbwt->Update(&data[1], 2);

// #if 1
        // Approach: split range [0, N-1] into M bins and compute which bins have
        // >0 alts present.
        // int n_steps = std::ceil((float)2*n_samples / 128);
        // const uint32_t step_size = std::ceil((float)2*n_samples / n_steps);
        // uint32_t n_bins = std::ceil((float)n_steps/64);
        // uint64_t* bins1 = new uint64_t[n_bins]; // todo: fix memory
        // memset(bins1, 0, sizeof(uint64_t)*n_bins);

        // for (int i = 0; i < 2*n_samples; ++i) {
        //     if (base_models[0].pbwt->prev[i]) {
        //         bins1[i/step_size/64] |= (1L << ((i/step_size) % 64));
        //         // std::cerr << i << "," << i/step_size/64 << "," << ((i/step_size)%64) << " step=" << step_size << " " << std::bitset<64>((1L << ((i/step_size) % 64))) << std::endl;
        //     }
        // }
        // std::cerr << std::bitset<64>(bins1[0]) << " " << std::bitset<64>(bins1[1]) << std::endl;

        uint64_t n_wah = std::ceil((float)2*n_samples / 64);
        uint64_t* wah_bitmaps = new uint64_t[n_wah];
        memset(wah_bitmaps, 0, n_wah*sizeof(uint64_t));

        for (int i = 0; i < 2*n_samples; ++i) {
            if (base_models[0].pbwt->prev[i]) {
                wah_bitmaps[i / 64] |= 1L << (i % 64);
            }
        }

        // RLE WAHs
        // uint64_t* wah_ref_ptr = &wah_bitmaps[0];
        uint64_t wah_ref = wah_bitmaps[0];
        uint64_t wah_run = 1;
        uint64_t observed_alts = 0;
        for (int i = 1; i < n_wah; ++i) {
            if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint64_t>::max()) || (wah_ref != wah_bitmaps[i])) {
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
                wah_ref = wah_bitmaps[i];
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
        // std::cerr << observed_alts << "/" << base_models[0].pbwt->n_queue[1] << std::endl;
        if (observed_alts != base_models[0].pbwt->n_queue[1]) {
            for (int i = 0; i < n_wah; ++i) {
                std::cerr << std::bitset<64>(wah_bitmaps[i]) << " ";
            }
            std::cerr << std::endl;
            exit(1);
        }
        // assert(observed_alts == base_models[0].pbwt->n_queue[1]);

        //

        // temp
        // uint8_t ref = base_models[0].pbwt->prev[0];
        // uint32_t n_run = 1;
        // for (int i = 1; i < 2*n_samples; ++i) {
        //     if (ref != base_models[0].pbwt->prev[i] || n_run == 4294967296) { // run has to be < 2^32
        //         ctx_model->encode_run_len(ref, n_run);

        //         ref = base_models[0].pbwt->prev[i];
        //         n_run = 0;
        //     }
        //     ++n_run;
        // }

        // if (n_run) {
        //     ctx_model->encode_run_len(ref, n_run);
        // }
        // //
        delete[] wah_bitmaps;
        
        ++processed_lines_local;
        ++processed_lines;
        return 1;

        /*
        uint8_t pack1 = 0; uint8_t n_pack1 = 0;

        // base_models[0].ResetContext();
        // base_model_bitmaps[0].ResetContext();
        uint32_t offset = 0, offset_end = 0;
        uint32_t n_obs = 0;

        for (int i = 0; i < n_steps; ++i) {
            // std::cerr << i << "," << i/64 << "," << (i%64) << std::endl;
            if (bins1[i/64] & (1L << (i%64))) {
                // std::cerr << "heerre" << std::endl;
                // base_model_bitmaps[0].EncodeSymbol(1);
                offset_end = offset + step_size < 2*n_samples ? offset + step_size : 2*n_samples;
                for (int j = offset; j < offset_end; ++j) {
                    // base_models[0].EncodeSymbol(base_models[0].pbwt->prev[j]);

                    if (n_pack1 % 8 == 0 && n_pack1 != 0) {
                        ppm1.Encode(pack1);
                        // pack1_context <<= 8;
                        // pack1_context |= pack1;
                        pack1 = 0; n_pack1 = 0;

                        // if (pack1_context != 0) std::cerr << "(int)" << std::bitset<32>(pack1_context) << std::endl;
                        
                    } 
                    // std::cerr << "adding: " << (int)base_models[0].pbwt->prev[j] << "@" << j << std::endl;
                    n_obs += (base_models[0].pbwt->prev[j] == 1);
                    pack1 |= base_models[0].pbwt->prev[j] << n_pack1;
                    ++n_pack1;
                    
                }
            } 
            // else base_model_bitmaps[0].EncodeSymbol(0);
            offset += step_size;
        }

        // std::cerr << n_obs << "==" << base_models[0].pbwt->n_queue[1] << std::endl;
        // assert(n_obs == base_models[0].pbwt->n_queue[1]);

#if DEBUG_PBWT
// assert(debug_bins[0].UpdateDigest((uint8_t*)bins1, n_bins*sizeof(uint64_t)));

for (int j = 0; j < n_bins; ++j) {
    uint64_t ref = bins1[j];
    for (int i = 0; i < 8; ++i) {
        debug_bins[0].buffer[debug_bins[0].len] = (ref & 255);
        ++debug_bins[0].len;
    //    assert(debug_bins[0].UpdateDigest(ref & 255));
       ref >>= 8;
    }
}
#endif

        for (int j = 0; j < n_bins; ++j) {
            // std::cerr << std::bitset<64>(bins1[j]) << " ";
            uint64_t ref = bins1[j];
            for (int i = 0; i < 8; ++i) {
                ppm_bin1.Encode(ref & 255);
                // pack1bin_context <<= 8;
                // pack1bin_context |= (ref & 255);
                ref >>= 8;
            }
        }
        // std::cerr << std::endl;

    
        if (n_pack1) {
            ppm1.Encode(pack1);
            // pack1_context <<= 8;
            // pack1_context |= pack1;
        }

        delete[] bins1;

        // uint64_t bins2 = 0;
        // for (int i = 0; i < n_samples; ++i) {
        //     if (base_models[1].pbwt->prev[i]) {
        //         bins2 |= (1 << (i/step_size));
        //     }
        // }

        // uint8_t pack2 = 0; uint8_t n_pack2 = 0;

        // // base_models[1].ResetContext();
        // // base_model_bitmaps[1].ResetContext();
        // offset = 0, offset_end = 0;
 
        // for (int i = 0; i < n_steps; ++i) {
        //     if (bins2 & (1 << i)) {
        //         // base_model_bitmaps[1].EncodeSymbol(1);
        //         offset_end = offset + step_size < n_samples ? offset + step_size : n_samples;
        //         for (int j = offset; j < offset_end; ++j) {
        //             // base_models[1].EncodeSymbol(base_models[1].pbwt->prev[j]);
        //             if (n_pack2 % 8 == 0 && n_pack2 != 0) {
        //                 ppm2.Encode(pack2, pack2_context);
        //                 pack2_context <<= 8;
        //                 pack2_context |= pack2;
        //                 pack2 = 0; n_pack2 = 0;
        //                 // if (pack1_context != 0) std::cerr << "(int)" << std::bitset<32>(pack1_context) << std::endl;
                        
        //             } else {
        //                 // std::cerr << "adding" << std::endl;
        //                 pack2 |= base_models[1].pbwt->prev[j] << n_pack2;
        //                 ++n_pack2;
        //             }
        //         }
        //     } 
        //     // else base_model_bitmaps[1].EncodeSymbol(0);
        //     offset += step_size;
        // }

        // if (n_pack2) {
        //     ppm2.Encode(pack2, pack2_context);
        //     pack2_context <<= 8;
        //     pack2_context |= pack2;
        // }

        // for (int i = 0; i < 8; ++i) {
        //     ppm_bin2.Encode(bins2 & 255, pack2bin_context);
        //     pack2bin_context <<= 8;
        //     pack2bin_context |= (bins2 & 255);
        //     bins2 <<= 8;
        // }

#else
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
*/
    } // end permute pbwt
    else {
#if 1
        // Approach: split range [0, N-1] into M bins and compute which bins have
        // >0 alts present.
        int n_steps = 16;
        const uint32_t step_size = std::ceil((float)n_samples / n_steps);
        uint64_t bins1 = 0;
        for (int i = 0, j = 0; i < n_samples; ++i, j += 2) {
            if (BCF_UNPACK_GENOTYPE(data[j])) {
                bins1 |= (1 << (i/step_size));
            }
        }

        base_models[0].ResetContext();
        base_model_bitmaps[0].ResetContext();
        uint32_t offset = 0, offset_end = 0;
        for (int i = 0; i < n_steps; ++i) {
            if (bins1 & (1 << i)) {
                base_model_bitmaps[0].EncodeSymbol(1);
                offset_end = offset + step_size < n_samples ? offset + step_size : n_samples;
                for (int j = offset; j < offset_end; ++j) {
                    base_models[0].EncodeSymbol(BCF_UNPACK_GENOTYPE(data[0+j*2]));
                }
            } else base_model_bitmaps[0].EncodeSymbol(0);
            offset += step_size;
        }

        uint64_t bins2 = 0;
        for (int i = 0, j = 1; i < n_samples; ++i, j += 2) {
            if (BCF_UNPACK_GENOTYPE(data[j])) {
                bins2 |= (1 << (i/step_size));
            }
        }

        base_models[1].ResetContext();
        base_model_bitmaps[1].ResetContext();
        offset = 0, offset_end = 0;
        for (int i = 0; i < n_steps; ++i) {
            if (bins2 & (1 << i)) {
                base_model_bitmaps[1].EncodeSymbol(1);
                offset_end = offset + step_size < n_samples ? offset + step_size : n_samples;
                for (int j = offset; j < offset_end; ++j) {
                    base_models[1].EncodeSymbol(BCF_UNPACK_GENOTYPE(data[1+j*2]));
                }
            } else base_model_bitmaps[1].EncodeSymbol(0);
            offset += step_size;
        }
#else
        base_models[0].ResetContext();
        for (int j = 0; j < n_samples; ++j) {
            base_models[0].EncodeSymbol(BCF_UNPACK_GENOTYPE(data[0+j*2]));
        }

        base_models[1].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            base_models[1].EncodeSymbol(BCF_UNPACK_GENOTYPE(data[1+i*2]));
        }
#endif
    }

    ++processed_lines_local;
    ++processed_lines;
    return 1;
    
}

// 2N2M with missing
int GenotypeCompressorModelling::Encode2N2MM(uint8_t* data, const int32_t n_data) {
#if DEBUG_PBWT
    assert(debug_pbwt[2].UpdateDigestStride(&data[0], n_data, 2));
    assert(debug_pbwt[3].UpdateDigestStride(&data[1], n_data, 2));
#endif

    ++base_models[2];
    ++base_models[3];

    if (permute_pbwt) {
        // std::cerr << "Encode2N2MM" << std::endl;
        
        base_models[2].pbwt->Update(&data[0], 2);
        base_models[3].pbwt->Update(&data[1], 2);

        base_models[2].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            // assert(base_models[2].pbwt->prev[i] < 4);
            base_models[2].EncodeSymbol(base_models[2].pbwt->prev[i]);
        }

        base_models[3].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            // assert(base_models[3].pbwt->prev[i] < 4);
            base_models[3].EncodeSymbol(base_models[3].pbwt->prev[i]);
        }
    } // end if permute pbwt
    else {
        base_models[2].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            // assert(base_models[2].pbwt->prev[i] < 4);
            assert(BCF_UNPACK_GENOTYPE(data[0+i*2]) < 4);
            base_models[2].EncodeSymbol(BCF_UNPACK_GENOTYPE(data[0+i*2]));
        }

        base_models[3].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            // assert(base_models[3].pbwt->prev[i] < 4);
            assert(BCF_UNPACK_GENOTYPE(data[1+i*2]) < 4);
            base_models[3].EncodeSymbol(BCF_UNPACK_GENOTYPE(data[1+i*2]));
        }
    } // end no permute pbwt

    ++processed_lines_local;
    ++processed_lines;
    return 1;
}

// 2N any M (up to 16)
int GenotypeCompressorModelling::Encode2NXM(uint8_t* data, const int32_t n_data, const int32_t n_alleles) {

#if DEBUG_PBWT
    assert(debug_pbwt[4].UpdateDigestStride(&data[0], n_data, 2));
    assert(debug_pbwt[5].UpdateDigestStride(&data[1], n_data, 2));
#endif

    ++base_models_complex[0];
    ++base_models_complex[1];

    // std::cerr << "Ecndoe2nXM: " << n_alleles << std::endl;
    if (permute_pbwt) {
        // Todo: assert genotypes are set for this variant.
        base_models_complex[0].pbwt->UpdateGeneral(&data[0], 2);
        base_models_complex[1].pbwt->UpdateGeneral(&data[1], 2);

        base_models_complex[0].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            // assert(base_models_complex[0].pbwt->prev[i] < 16);
            base_models_complex[0].EncodeSymbol(base_models_complex[0].pbwt->prev[i]);
        }

        base_models_complex[1].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            // assert(base_models_complex[1].pbwt->prev[i] < 16);
            base_models_complex[1].EncodeSymbol(base_models_complex[1].pbwt->prev[i]);
        }
    } // end permute pbwt
    else {
        base_models_complex[0].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            // assert(base_models_complex[0].pbwt->prev[i] < 16);
            base_models_complex[0].EncodeSymbol(BCF_UNPACK_GENOTYPE_GENERAL(data[0+i*2]));
        }

        base_models_complex[1].ResetContext();
        for (int i = 0; i < n_samples; ++i) {
            // assert(base_models_complex[1].pbwt->prev[i] < 16);
            base_models_complex[1].EncodeSymbol(BCF_UNPACK_GENOTYPE_GENERAL(data[1+i*2]));
        }
    }

    ++processed_lines_local;
    ++processed_lines;
    return 1;
}

int GenotypeCompressorModelling::Compress(djinn_block_t*& block) {
    if (block == nullptr) {
        block = new djinn_ctx_block_t();
        block->type = djinn_block_t::BlockType::CONTEXT;
        block->data = new djinn_ctx_t();
    }
    djinn_ctx_t* data_out = (djinn_ctx_t*)block->data;
    data_out->reset();

#if DEBUG_PBWT
    // Finish digests.
    for (int i = 0; i < 6; ++i) {
        if (debug_pbwt[i].len) {
            debug_pbwt[i].FinalizeDigest();
            std::cerr << "[SHA512] Digest-" << i << "=" << std::hex << (int)debug_pbwt[i].digest[0];
            for (int j = 1; j < 64; ++j) std::cerr << std::hex << (int)debug_pbwt[i].digest[j];
            std::cerr << std::dec << std::endl;
        }
    }

    debug_bins[0].FinalizeDigest();
    std::cerr << "[SHA512] Digest-bins-" << 0 << "=" << std::hex << (int)debug_bins[0].digest[0];
    for (int j = 1; j < 64; ++j) std::cerr << std::hex << (int)debug_bins[0].digest[j];
    std::cerr << std::dec << std::endl;
    std::cerr << "Len=" << debug_bins[0].len << "->" << debug_bins[0].len/processed_lines_local << "->" << debug_bins[0].len/processed_lines_local/8 << std::endl;
#endif

    size_t p1   = base_models[0].FinishEncoding();
    size_t p2   = base_models[1].FinishEncoding();
    size_t p1E  = base_models[2].FinishEncoding();
    size_t p2E  = base_models[3].FinishEncoding();
    size_t p2X  = base_models_complex[0].FinishEncoding();
    size_t p2X2 = base_models_complex[1].FinishEncoding();
    size_t extra1 = base_model_bitmaps[0].FinishEncoding();
    size_t extra2 = base_model_bitmaps[1].FinishEncoding();

    ppm1.encoder.Flush();
    ppm2.encoder.Flush();
    ppm_bin1.encoder.Flush();
    ppm_bin2.encoder.Flush();
    size_t s_pp1m = ppm1.encoder.dat-ppm1.encoder.buffer;
    size_t s_pp2m = ppm2.encoder.dat-ppm2.encoder.buffer;
    size_t s_pp1m_bin1 = ppm_bin1.encoder.dat-ppm_bin1.encoder.buffer;
    size_t s_pp2m_bin2 = ppm_bin2.encoder.dat-ppm_bin2.encoder.buffer;
    
    std::cerr << "ppm1=" << (int)s_pp1m << " (" << (float)p1/s_pp1m << "-fold)" << std::endl;
    std::cerr << "ppm2=" << (int)s_pp2m << " (" << (float)p2/s_pp2m << "-fold)" << std::endl;
    std::cerr << "ppm1bin=" << (int)s_pp1m_bin1 << " (" << (float)extra1/s_pp1m_bin1 << "-fold)" << std::endl;
    std::cerr << "ppm2bin=" << (int)s_pp2m_bin2 << " (" << (float)extra2/s_pp2m_bin2 << "-fold)" << std::endl;
    
    size_t smref = mref->FinishEncoding();
    size_t smlrle = mlog_rle->FinishEncoding();
    size_t smrle = mrle->FinishEncoding();
    size_t smrle2_1 = mrle2_1->FinishEncoding();
    size_t smrle2_2 = mrle2_2->FinishEncoding();

    size_t smrle4_1 = mrle4_1->FinishEncoding();
    size_t smrle4_2 = mrle4_2->FinishEncoding();
    size_t smrle4_3 = mrle4_3->FinishEncoding();
    size_t smrle4_4 = mrle4_4->FinishEncoding();
    
    size_t sdirty = dirty_wah->FinishEncoding();
    size_t stype = mtype->FinishEncoding();

    std::cerr << "[TEST] REF=" << smref << " LOG-RLE=" << smlrle << " RLE-1=" << smrle << " RLE-2=" << smrle2_1 << "," << smrle2_2 << " RLE-4=" << smrle4_1 << "," << smrle4_2 << "," << smrle4_3 << "," << smrle4_4 << " dirty=" << sdirty << "==" << smref + smlrle + smrle + smrle2_1 + smrle2_2 + smrle4_1 + smrle4_2 + smrle4_3 + smrle4_4 + sdirty + stype << std::endl;
    std::cerr << "[TEST] TYPE=" << stype << std::endl; 
    bytes_out4 += smref + smlrle + smrle + smrle2_1 + smrle2_2 + smrle4_1 + smrle4_2 + smrle4_3 + smrle4_4 + sdirty + stype;
    mref->StartEncoding();
    mlog_rle->StartEncoding();
    mrle->StartEncoding();
     mrle2_1->StartEncoding();
  mrle2_2->StartEncoding();

   mrle4_1->StartEncoding();
   mrle4_2->StartEncoding();
  mrle4_3->StartEncoding();
    mrle4_4->StartEncoding();

 dirty_wah->StartEncoding();
   mtype->StartEncoding();

    // ctx_model->rce->End();
    // out_gts += ctx_model->data->pos - ctx_model->data->buffer;
    // std::cerr << "[GTSHARK] " << bytes_in << "->" << out_gts << " (" << (double)bytes_in/out_gts << "-fold ubcf, " << (double)bytes_in_vcf/out_gts << "-fold vcf)" << std::endl;
    // ctx_model->data->pos = ctx_model->data->buffer;
    
    // ctx_model = std::make_shared<djinn_gt_ctx>();
    // ctx_model->data = std::make_shared<djinn_ctx_buf_t>(10000000);
    // ctx_model->rce = new CRangeEncoder();
    // ctx_model->rce->Start();

    
#if DEBUG_PBWT
    TPPM test2; 
    test2.decoder.buffer = ppm_bin1.encoder.buffer;
    test2.decoder.dat = ppm_bin1.encoder.buffer;
    test2.decoder.Init();

    uint32_t tester = 0;
    for (int i = 0; i < debug_bins[0].len; ++i) {
        uint8_t decoded = test2.Decode();
        if(debug_bins[0].buffer[i] != decoded) {
            std::cerr << i << ": " << (int)debug_bins[0].buffer[i] << "!=" << (int)decoded << std::endl;
            exit(1);
        }
        // ++tester;
        // std::cerr << decoded << ",";
    }
    std::cerr << "[TEST]" << "PASSED" << std::endl;
    debug_bins[0].len = 0;
#endif

    ppm1.reset();
    ppm2.reset();
    ppm_bin1.reset();
    ppm_bin2.reset();
    // ppm1.encoder.dat = ppm1.encoder.buffer;
    // ppm2.encoder.dat = ppm2.encoder.buffer;
    // ppm_bin1.encoder.dat = ppm_bin1.encoder.buffer;
    // ppm_bin2.encoder.dat = ppm_bin2.encoder.buffer;
    pack1bin_context = 0;
    pack2bin_context = 0;
    pack1_context = 0;
    pack2_context = 0;

    bytes_out2 += s_pp1m + s_pp2m + s_pp1m_bin1 + s_pp2m_bin2;
    std::cerr << "[MODEL2] " << bytes_in << "->" << bytes_out2 << " (" << (double)bytes_in/bytes_out2 << "-fold ubcf, " << (double)bytes_in_vcf/bytes_out2 << "-fold vcf)" << std::endl;
    std::cerr << "[MODEL-RLEC] " << bytes_in << "->" << bytes_out4 << " (" << (double)bytes_in/bytes_out4 << "-fold ubcf, " << (double)bytes_in_vcf/bytes_out4 << "-fold vcf)" << std::endl;

    // size_t praw = ZstdCompress(buf_raw.data, buf_raw.len,
    //                            buf_compress.data, buf_compress.capacity(),
    //                            20);

    if (base_models[0].n_additions) {
        assert(base_models[1].n_additions);
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MC1].SetDataReference(base_models[0].buffer, p1);
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MC1].n = base_models[0].n_additions;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MC1].n_c = p1;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MC1].n_v = base_models[0].n_variants;
        assert(memcmp(base_models[0].buffer, data_out->ctx_models[DJINN_CTX_MODEL_DI2MC1].vptr, p1) == 0);

        data_out->ctx_partitions[DJINN_CTX_MODEL_DI2MC_PART1].SetDataReference(base_model_bitmaps[0].buffer, extra1);
        data_out->ctx_partitions[DJINN_CTX_MODEL_DI2MC_PART1].n = base_model_bitmaps[0].n_additions;
        data_out->ctx_partitions[DJINN_CTX_MODEL_DI2MC_PART1].n_c = extra1;
        assert(memcmp(base_model_bitmaps[0].buffer, data_out->ctx_partitions[DJINN_CTX_MODEL_DI2MC_PART1].vptr, extra1) == 0);

        std::cerr << processed_lines_local << "," << base_models[0].n_additions << "," << p1 << "," << base_model_bitmaps[0].n_additions << "," << extra1 << std::endl;
    }

    if (base_models[1].n_additions) {
        assert(base_models[0].n_additions);

        data_out->ctx_models[DJINN_CTX_MODEL_DI2MC2].SetDataReference(base_models[1].buffer, p2);
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MC2].n = base_models[1].n_additions;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MC2].n_c = p2;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MC2].n_v = base_models[1].n_variants;

        data_out->ctx_partitions[DJINN_CTX_MODEL_DI2MC_PART2].SetDataReference(base_model_bitmaps[1].buffer, extra2);
        data_out->ctx_partitions[DJINN_CTX_MODEL_DI2MC_PART2].n = base_model_bitmaps[1].n_additions;
        data_out->ctx_partitions[DJINN_CTX_MODEL_DI2MC_PART2].n_c = extra2;

        std::cerr << processed_lines_local << "," << base_models[1].n_additions << "," << p2 << "," << base_model_bitmaps[1].n_additions << "," << extra2 << std::endl;
    }

    if (base_models[2].n_additions) {
        assert(base_models[3].n_additions);
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MM1].SetDataReference(base_models[2].buffer, p1E);
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MM1].n = base_models[2].n_additions;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MM1].n_c = p1E;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MM1].n_v = base_models[2].n_variants;
    }

    if (base_models[3].n_additions) {
        assert(base_models[2].n_additions);

        data_out->ctx_models[DJINN_CTX_MODEL_DI2MM2].SetDataReference(base_models[3].buffer, p2E);
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MM2].n = base_models[3].n_additions;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MM2].n_c = p2E;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2MM2].n_v = base_models[3].n_variants;
    }

    if (base_models_complex[0].n_additions) {
        assert(base_models_complex[1].n_additions);
        data_out->ctx_models[DJINN_CTX_MODEL_DI2X1].SetDataReference(base_models_complex[0].buffer, p2X);
        data_out->ctx_models[DJINN_CTX_MODEL_DI2X1].n = base_models_complex[0].n_additions;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2X1].n_c = p2X;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2X1].n_v = base_models_complex[0].n_variants;
    }

    if (base_models_complex[1].n_additions) {
        assert(base_models_complex[0].n_additions);

        data_out->ctx_models[DJINN_CTX_MODEL_DI2X2].SetDataReference(base_models_complex[1].buffer, p2X2);
        data_out->ctx_models[DJINN_CTX_MODEL_DI2X2].n = base_models_complex[1].n_additions;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2X2].n_c = p2X2;
        data_out->ctx_models[DJINN_CTX_MODEL_DI2X2].n_v = base_models_complex[1].n_variants;
    }

    block->n_rcds = processed_lines_local;
    block->ctrl = data_out->GetController();
    block->p_len = data_out->SerializedSize();
    djinn_ctx_block_t* oblock = (djinn_ctx_block_t*)block;
    djinn_ctx_ctrl_t* ctrl = (djinn_ctx_ctrl_t*)&oblock->ctrl;
    ctrl->pbwt = permute_pbwt;
    // oblock->Serialize(std::cout);
    // std::cout.write((char*)base_models[0].buffer, p1);
    // std::cout.flush();
    // std::cerr << "Written size=" << block->p_len << " ctrl=" << std::bitset<16>(block->ctrl) << std::endl;

#if DEBUG_PBWT
    // for (int i = 0; i < 4; ++i) std::cerr << "variants=" << base_models[i].n_variants << std::endl;
    // for (int i = 0; i < 2; ++i) std::cerr << "variants=" << base_models_complex[i].n_variants << std::endl;
    DebugContext(base_models[0].buffer, p1,  base_model_bitmaps[0].buffer, extra1, debug_pbwt[0].buffer, base_models[0].n_variants, 2, &GenotypeDecompressorContext::Decode, TWK_BCF_GT_UNPACK);
    DebugContext(base_models[1].buffer, p2,  base_model_bitmaps[1].buffer, extra2, debug_pbwt[1].buffer, base_models[1].n_variants, 2, &GenotypeDecompressorContext::Decode, TWK_BCF_GT_UNPACK);
    DebugContext(base_models[2].buffer, p1E, debug_pbwt[2].buffer, base_models[2].n_variants, 4, &GenotypeDecompressorContext::Decode2, TWK_BCF_GT_UNPACK);
    DebugContext(base_models[3].buffer, p2E, debug_pbwt[3].buffer, base_models[3].n_variants, 4, &GenotypeDecompressorContext::Decode2, TWK_BCF_GT_UNPACK);
    DebugContext(base_models_complex[0].buffer, p2X,  debug_pbwt[4].buffer, base_models_complex[0].n_variants, 16, &GenotypeDecompressorContext::Decode2, TWK_BCF_GT_UNPACK_GENERAL);
    DebugContext(base_models_complex[1].buffer, p2X2, debug_pbwt[5].buffer, base_models_complex[1].n_variants, 16, &GenotypeDecompressorContext::Decode2, TWK_BCF_GT_UNPACK_GENERAL);
#endif

    // base_models[0].Reset();
    // base_models[1].Reset();
    // base_models[2].Reset();
    // base_models[3].Reset();
    // base_models_complex[0].Reset();
    // base_models_complex[1].Reset();
    // base_model_bitmaps[0].Reset();
    // base_model_bitmaps[1].Reset();

    // base_models[0].StartEncoding();
    // base_models[1].StartEncoding();
    // base_models[2].StartEncoding();
    // base_models[3].StartEncoding();
    // base_models_complex[0].StartEncoding();
    // base_models_complex[1].StartEncoding();
    // base_model_bitmaps[0].StartEncoding();
    // base_model_bitmaps[1].StartEncoding();

    bytes_out += p1 + p2 + p1E + p2E + p2X + p2X2 + extra1 + extra2;
    // bytes_out += s_pp1m + s_pp2m + extra1 + extra2;
#if DEBUG_SIZE
    std::cerr << "[WRITE] Variants=" << processed_lines_local << " 2N2MC=" << p1 << "," << p2 << " 2N2MM=" << p1E << "," << p2E << " 2NXM=" << p2X << "," << p2X2 << " SKIP=" << extra1 << "," << extra2 << std::endl;
    std::cerr << "[PROGRESS] " << bytes_in << "->" << bytes_out << " (" << (double)bytes_in/bytes_out << "-fold ubcf, " << (double)bytes_in_vcf/bytes_out << "-fold vcf)" << std::endl;
#endif
    processed_lines_local = 0;
    buf_raw.reset();
    
#if DEBUG_PBWT
        // Reset digests.
        for (int i = 0; i < 6; ++i) {
            debug_pbwt[i].reset();
        }
        debug_bins[0].reset();
#endif

    return oblock->size();
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
    
    base_pbwt[0].Initiate(n_s, 2);
    base_pbwt[1].Initiate(n_s, 2);
    base_pbwt[2].Initiate(n_s, 4);
    base_pbwt[3].Initiate(n_s, 4);

    complex_pbwt[0].Initiate(n_s, 16);
    complex_pbwt[1].Initiate(n_s, 16);
}

GenotypeCompressorRLEBitmap::~GenotypeCompressorRLEBitmap() { }

int GenotypeCompressorRLEBitmap::Encode2N(uint8_t* data, const int32_t n_data, const int32_t n_alleles) {
    if (data == nullptr) return -1;
    if (n_data == 0) return -2;
    if (n_alleles == 0) return -3;

    // std::cerr << "n_alleles=" << n_alleles << std::endl;

    if (n_alleles == 2) return(Encode2N2M(data, n_data)); // biallelic
    else if (n_alleles < 14) return(Encode2NXM(data, n_data, n_alleles));  // alleles < 14 and 2 reserved for MISSING and EOV
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

int GenotypeCompressorRLEBitmap::Encode2N2M(uint8_t* data, const int32_t n_data) {
    // std::cerr << "2N2M" << std::endl;

    // Todo: assert genotypes are set for this variant.
    int replaced = RemapGenotypeEOV(data, n_data);
    if (replaced) {
        // std::cerr << "replaced: skipping " << replaced << std::endl;
        return Encode2N2MM(data, n_data); // 2N2MM supports missing+EOV in RLE-bitmap mode.
    }

    // for (int i = 0; i < 256; ++i) {
    //     if (alts[i]) std::cerr << i << ":" << alts[i] << ",";
    // }
    // std::cerr << std::endl;

    if (alts[15] == 0 && alts[14] == 0) { // No missing values and no EOV values
        return Encode2N2MC(data, n_data);
    } else { // Having missing values.
        // std::cerr << "using extended model" << std::endl;
        return Encode2N2MM(data, n_data);
    }

    return -1;
}

int GenotypeCompressorRLEBitmap::Encode2N2MC(uint8_t* data, const int32_t n_data) {
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
        base_pbwt[0].Update(&data[0], 2);
        base_pbwt[1].Update(&data[1], 2);

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

    if (permute_pbwt) {
        EncodeRLEBitmap2N2MC(1);
    } else {
        EncodeRLEBitmap2N2MC(1, &data[1], 2);
    }

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

bool GenotypeCompressorRLEBitmap::CheckLimit() const {
    return (processed_lines_local == block_size || 
            buf_wah[0].len > 9000000 || 
            buf_wah[1].len > 9000000 ||
            buf_wah[2].len > 9000000 ||
            buf_raw.len > 9000000 );
}

int GenotypeCompressorRLEBitmap::Encode2N2MM(uint8_t* data, const int32_t n_data) {
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
        base_pbwt[2].Update(&data[0], 2);
        base_pbwt[3].Update(&data[1], 2);

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
    
    if (permute_pbwt) {
        int ret2 = EncodeRLEBitmap2N2MM(3);
    } else {
        int ret2 = EncodeRLEBitmap2N2MM(3, &data[1], 2);
    }

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

int GenotypeCompressorRLEBitmap::Encode2NXM(uint8_t* data, const int32_t n_data, const int32_t n_alleles) { 
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
        complex_pbwt[0].UpdateGeneral(&data[0], 2);
        complex_pbwt[1].UpdateGeneral(&data[1], 2);
        
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

    if (permute_pbwt) {
        EncodeRLEBitmap2NXM(1);
    } else {
        EncodeRLEBitmap2NXM(1, &data[1], 2);
    }

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

int GenotypeCompressorRLEBitmap::Compress(djinn_block_t*& block) {
    if (block == nullptr) {
        block = new djinn_wah_block_t();
        block->type = djinn_block_t::BlockType::WAH;
        block->data = new djinn_wah_t();
    } else {
        // Previous block type was not WAH
        if (block->type != djinn_block_t::BlockType::WAH) {
            delete block;
            block = new djinn_wah_block_t();
            block->type = djinn_block_t::BlockType::WAH;
            block->data = new djinn_wah_t();
        }
    }
    djinn_wah_t* data_out = (djinn_wah_t*)block->data;
    data_out->reset();

#if DEBUG_PBWT
    // Finish digests.
    for (int i = 0; i < 6; ++i) {
        if (debug_pbwt[i].len) {
            debug_pbwt[i].FinalizeDigest();
            std::cerr << "[SHA512] Digest-" << i << "=" << std::hex << (int)debug_pbwt[i].digest[0];
            for (int j = 1; j < 64; ++j) std::cerr << std::hex << (int)debug_pbwt[i].digest[j];
            std::cerr << std::dec << std::endl;
        }
    }
#endif

#if DEBUG_PBWT
    DebugWAH(buf_wah[0].data, buf_wah[0].len, debug_pbwt[0].buffer, n_variants[0], 2,  &GenotypeDecompressorRLEBitmap::Decode2N2MC, TWK_BCF_GT_UNPACK);
    DebugWAH(buf_wah[1].data, buf_wah[1].len, debug_pbwt[1].buffer, n_variants[1], 2,  &GenotypeDecompressorRLEBitmap::Decode2N2MC, TWK_BCF_GT_UNPACK);
    DebugWAH(buf_wah[2].data, buf_wah[2].len, debug_pbwt[2].buffer, n_variants[2], 4,  &GenotypeDecompressorRLEBitmap::Decode2N2MM, TWK_BCF_GT_UNPACK);
    DebugWAH(buf_wah[3].data, buf_wah[3].len, debug_pbwt[3].buffer, n_variants[3], 4,  &GenotypeDecompressorRLEBitmap::Decode2N2MM, TWK_BCF_GT_UNPACK);
    DebugWAH(buf_wah[4].data, buf_wah[4].len, debug_pbwt[4].buffer, n_variants[4], 16, &GenotypeDecompressorRLEBitmap::Decode2NXM,  TWK_BCF_GT_UNPACK_GENERAL);
    DebugWAH(buf_wah[5].data, buf_wah[5].len, debug_pbwt[5].buffer, n_variants[5], 16, &GenotypeDecompressorRLEBitmap::Decode2NXM,  TWK_BCF_GT_UNPACK_GENERAL);
#endif

    // Compress and add to djinn block.
    for (int i = 0; i < 6; ++i) {
        if (buf_wah[i].len) {
            size_t praw = 0;
            if (strategy == CompressionStrategy::ZSTD) {
                praw = ZstdCompress(buf_wah[i].data, buf_wah[i].len,
                                    buf_compress.data, buf_compress.capacity(),
                                    20);

                bytes_out_zstd1 += praw;
                // std::cerr << i << ":Zstd->compressed: " << buf_wah[i].len << "->" << praw << " (" << (float)buf_wah[i].len/praw << ")" << std::endl;
            }

            if (strategy == CompressionStrategy::LZ4) {
                praw = Lz4Compress(buf_wah[i].data, buf_wah[i].len,
                                buf_compress.data, buf_compress.capacity(),
                                9);
                
                // std::cerr << i << ":Lz4->compressed: " << buf_wah[i].len << "->" << praw << " (" << (float)buf_wah[i].len/praw << ")" << std::endl;
                bytes_out_lz4 += praw;       
            }

            data_out->wah_models[i].SetData(buf_compress.data, praw); // copy data because the internal buffer is immediately reused
            data_out->wah_models[i].n = buf_wah[i].len;
            data_out->wah_models[i].n_c = praw;
            data_out->wah_models[i].n_v = n_variants[i];
#if DEBUG_WAH
            assert(memcmp(buf_compress.data, data_out->wah_models[i].vptr, praw) == 0);
#endif
            // Reset
            buf_wah[i].len = 0;
        }
    }

#if DEBUG_PBWT
        // Reset digests.
        for (int i = 0; i < 6; ++i) {
            debug_pbwt[i].reset();
        }
#endif

    base_pbwt[0].reset();
    base_pbwt[1].reset();
    base_pbwt[2].reset();
    base_pbwt[3].reset();
    complex_pbwt[0].reset();
    complex_pbwt[1].reset();

    processed_lines_local = 0;
    buf_raw.reset();
    for (int i = 0; i < 6; ++i) n_variants[i] = 0;

    //
    block->n_rcds = processed_lines_local;
    block->ctrl = data_out->GetController();
    block->p_len = data_out->SerializedSize();
    djinn_wah_block_t* oblock = (djinn_wah_block_t*)block;
    djinn_wah_ctrl_t* ctrl = (djinn_wah_ctrl_t*)&oblock->ctrl;
    ctrl->pbwt = permute_pbwt;
    // std::cerr << "Written size=" << block->p_len << " ctrl=" << std::bitset<16>(block->ctrl) << std::endl;
    //
    
    // Debug compression:
#if DEBUG_SIZE
    std::cerr << "[PROGRESS ZSTD] " << bytes_in << "->" << bytes_out_zstd1 << " (" << (double)bytes_in/bytes_out_zstd1 << "-fold ubcf, " << (double)bytes_in_vcf/bytes_out_zstd1 << "-fold vcf)" << std::endl;
    std::cerr << "[PROGRESS LZ4] "  << bytes_in << "->" << bytes_out_lz4 << " (" << (double)bytes_in/bytes_out_lz4 << "-fold ubcf, " << (double)bytes_in_vcf/bytes_out_lz4 << "-fold vcf)" << std::endl;
#endif
    return oblock->size();
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

    for (int i = 1; i < n_samples; ++i) {
        if (prev[i] != ref || rle_len == 16383) { // 2^14-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 31 then the word is considered "dirty"
            // and a 32-bit bitmap will be used instead.
            if (end_rle - start_rle < 31) {
                // Increment position
                i += 31 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 31 < n_samples ? start_rle + 31 : n_samples;
                // Update observed path.
                observed_length += start_rle + 31 < n_samples ? 31 : n_samples - start_rle;
                // Assertion.
                assert(end_rle <= n_samples);
                
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
                rle_len = (end_rle == n_samples) ? 0 : 1;
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
    assert(observed_length == n_samples);

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

    for (int i = 1; i < n_samples; ++i) {
        if (BCF_UNPACK_GENOTYPE(data[i*stride]) != ref || rle_len == 16383) { // 2^14-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 31 then the word is considered "dirty"
            // and a 32-bit bitmap will be used instead.
            if (end_rle - start_rle < 31) {
                // Increment position
                i += 31 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 31 < n_samples ? start_rle + 31 : n_samples;
                // Update observed path.
                observed_length += start_rle + 31 < n_samples ? 31 : n_samples - start_rle;
                // Assertion.
                assert(end_rle <= n_samples);
                
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
                rle_len = (end_rle == n_samples) ? 0 : 1;
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
    assert(observed_length == n_samples);

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

    for (int i = 1; i < n_samples; ++i) {
        if (prev[i] != ref || rle_len == 8191) { // 2^13-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 15 then the word is considered "dirty"
            // and a 32-bit bitmap will be used instead.
            if (end_rle - start_rle < 15) { // 15*2 = 30 and 1 bit for class type = 31 bits
                // Increment position
                i += 15 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 15 < n_samples ? start_rle + 15 : n_samples;
                // Update observed path.
                observed_length += start_rle + 15 < n_samples ? 15 : n_samples - start_rle;
                // Assertion.
                assert(end_rle <= n_samples);
                
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
                rle_len = (end_rle == n_samples) ? 0 : 1;
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
    assert(observed_length == n_samples);

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

    for (int i = 1; i < n_samples; ++i) {
        if (BCF_UNPACK_GENOTYPE(data[i*stride]) != ref || rle_len == 8191) { // 2^13-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 15 then the word is considered "dirty"
            // and a 32-bit bitmap will be used instead.
            if (end_rle - start_rle < 15) { // 15*2 = 30 and 1 bit for class type = 31 bits
                // Increment position
                i += 15 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 15 < n_samples ? start_rle + 15 : n_samples;
                // Update observed path.
                observed_length += start_rle + 15 < n_samples ? 15 : n_samples - start_rle;
                // Assertion.
                assert(end_rle <= n_samples);
                
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
                rle_len = (end_rle == n_samples) ? 0 : 1;
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
    assert(observed_length == n_samples);

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

    for (int i = 1; i < n_samples; ++i) {
        if (prev[i] != ref || rle_len == 134217727) { // 2^27-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 15 then the word is considered "dirty"
            // and a 64-bit bitmap will be used instead.
            if (end_rle - start_rle < 15) { // 15*4 = 60 and 1 bit for class type = 31 bits
                // Increment position
                i += 15 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 15 < n_samples ? start_rle + 15 : n_samples;
                // Update observed path.
                observed_length += start_rle + 15 < n_samples ? 15 : n_samples - start_rle;
                // Assertion.
                assert(end_rle <= n_samples);
                
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
                rle_len = (end_rle == n_samples) ? 0 : 1;
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
    assert(observed_length == n_samples);

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

    for (int i = 1; i < n_samples; ++i) {
        if (BCF_UNPACK_GENOTYPE_GENERAL(data[i*stride]) != ref || rle_len == 134217727) { // 2^27-1
            end_rle = i;
            ++n_objects;

            // If the run length is < 15 then the word is considered "dirty"
            // and a 64-bit bitmap will be used instead.
            if (end_rle - start_rle < 15) { // 15*4 = 60 and 1 bit for class type = 31 bits
                // Increment position
                i += 15 - (end_rle - start_rle);
                // Make sure the end position is within range.
                end_rle = start_rle + 15 < n_samples ? start_rle + 15 : n_samples;
                // Update observed path.
                observed_length += start_rle + 15 < n_samples ? 15 : n_samples - start_rle;
                // Assertion.
                assert(end_rle <= n_samples);
                
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
                rle_len = (end_rle == n_samples) ? 0 : 1;
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
    assert(observed_length == n_samples);

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