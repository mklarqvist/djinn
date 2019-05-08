#include "gt_compressor.h"

#include <cstdio>//printf: debug
#include "gt_decompressor.h" //debug

/*======   Base model   ======*/

GenotypeCompressor::GenotypeCompressor(int64_t n_s) :
    n_samples(n_s),
    block_size(8192),
    processed_lines(0),
    processed_lines_local(0),
    bytes_in(0), bytes_out(0),
    strategy(CompressionStrategy::LZ4)
{
    buf_compress.resize(10000000);
    buf_raw.resize(10000000);

#if DEBUG_PBWT
    for (int i = 0; i < 6; ++i) {
        debug_pbwt[i].resize(n_s*block_size + 65536);
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

    // if (bcf->n_allele < 16 && bcf->n_allele > 2) {
    //     for (int i = 0; i < bcf->n_allele; ++i) {
    //         std::cerr << bcf->d.allele[i] << ",";
    //     }
    //     std::cerr << std::endl;
    // }

    // Todo: extend beyond 2N
    return(Encode2N(fmt->p, fmt->p_len, bcf->n_allele));
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
    base_models[0].Construct(n_samples, 2);
    base_models[1].Construct(n_samples, 2);
    base_models[2].Construct(n_samples, 3);
    base_models[3].Construct(n_samples, 3);
    base_models_complex[0].Construct(n_samples, 16);
    base_models_complex[1].Construct(n_samples, 16);

    base_models[0].StartEncoding();
    base_models[1].StartEncoding();
    base_models[2].StartEncoding();
    base_models[3].StartEncoding();
    base_models_complex[0].StartEncoding();
    base_models_complex[1].StartEncoding();

    nonsense[0].Construct(1, 2);
    nonsense[1].Construct(1, 2);
    nonsense[0].StartEncoding();
    nonsense[1].StartEncoding();
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

//
// int GenotypeCompressorModelling::Encode(bcf1_t* bcf, const bcf_hdr_t* hdr) {
//     if (bcf == NULL) return 0;
//     if (hdr == NULL) return 0;
    
//     const bcf_fmt_t* fmt = bcf_get_fmt(hdr, bcf, "GT");
//     if (fmt == NULL) return 0;
//     if (fmt->p_len / n_samples != 2) {
//         std::cerr << "input is not divisible by 2" << std::endl;
//         return 0;
//     }

//     bytes_in += bcf->d.fmt->p_len;

//     // Todo: extend beyond 2N
//     return(Encode2N(bcf,hdr));
// }

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
    else if (bcf->n_allele < 4) {

        int replaced = RemapGenotypeEOV(fmt->p, fmt->p_len);
        
        // std::cerr << "second path taken for " << bcf->n_allele << std::endl; 
        return(Encode2NXM(fmt->p, fmt->p_len, bcf->n_allele));  // #alleles < 4
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
    else if (n_alleles < 4) return(Encode2NXM(data, n_data, n_alleles));  // #alleles < 4
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
int GenotypeCompressorModelling::Encode2N2M(uint8_t* data, const int32_t n_data) {
    if (CheckLimit()) {
        Compress();
    }

    // Todo: assert genotypes are set for this variant.
    int replaced = RemapGenotypeEOV(data, n_data);
    if (replaced) {
        std::cerr << "replaced: divert to missing with " << replaced << std::endl;
        return Encode2NXM(data ,n_data, 4); // todo: see if this is true
    }

    if (alts[14] == 0 && alts[15] == 0) { // No missing values.
        return Encode2N2MC(data, n_data);
    } else { // Having missing values.
        // std::cerr << "using extended model" << std::endl;
        return Encode2N2MM(data, n_data);
    }

    return 1;
}

// 2N2M complete
int GenotypeCompressorModelling::Encode2N2MC(uint8_t* data, const int32_t n_data) {
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
int GenotypeCompressorModelling::Encode2N2MM(uint8_t* data, const int32_t n_data) {
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
int GenotypeCompressorModelling::Encode2NXM(uint8_t* data, const int32_t n_data, const int32_t n_alleles) {
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

int GenotypeCompressorModelling::Compress() {
    // flush: temp
    int p1  = base_models[0].FinishEncoding();
    int p2  = base_models[1].FinishEncoding();
    int p1E = base_models[2].FinishEncoding();
    int p2E = base_models[3].FinishEncoding();
    int p2X = base_models_complex[0].FinishEncoding();
    int p2X2 = base_models_complex[1].FinishEncoding();
    int praw = ZstdCompress(buf_raw.data, buf_raw.len,
                            buf_compress.data, buf_compress.capacity(),
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

/*======   RLE-bitmap approach   ======*/

GenotypeCompressorRLEBitmap::GenotypeCompressorRLEBitmap(int64_t n_s) : GenotypeCompressor(n_s),
    bytes_in1(0), bytes_out_zstd1(0), bytes_out_lz4(0)
{
    strategy = CompressionStrategy::LZ4;

    gt_width.resize(n_s, 0);
    buf_wah[0].resize(10000000);
    buf_wah[1].resize(10000000);
    buf_wah[2].resize(10000000);
    // buf_raw.resize(10000000);

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
    else if (n_alleles < 16) return(Encode2NXM(data, n_data, n_alleles));  // alleles < 16
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
    if (CheckLimit()) {
        Compress();
    }

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

    // Update diploid, biallelic, no-missing PBWT for haplotype 1 and 2.
    base_pbwt[0].Update(&data[0], 2);
    base_pbwt[1].Update(&data[1], 2);

#if DEBUG_WAH
    uint32_t buf_wah_before = buf_wah[0].len;
#endif
    EncodeRLEBitmap2N2MC(0);

#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    std::shared_ptr<GenotypeDecompressorRLEBitmap> debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[0].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts1 = debug1->Decode2N2MC(buf_compress.data);
    assert(alts1 == base_pbwt[0].n_queue[1]);
    // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 1.
    for (int i = 0; i < n_samples; ++i) {
        // std::cerr << i << " -> " << (int)base_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
        assert(base_pbwt[0].prev[i] == buf_compress.data[i]);
    }
    buf_wah_before = buf_wah[0].len;
#endif

    EncodeRLEBitmap2N2MC(1);

#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[0].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts2 = debug1->Decode2N2MC(buf_compress.data);
    assert(alts2 == base_pbwt[1].n_queue[1]);
    // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 2.
    for (int i = 0; i < n_samples; ++i) {
        // std::cerr << i << " -> " << (int)base_pbwt[1].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
        assert(base_pbwt[1].prev[i] == buf_compress.data[i]);
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
    // std::cerr << "2N2MM" << std::endl;

    // Update diploid, biallelic, no-missing PBWT for haplotype 1 and 2.
    base_pbwt[2].Update(&data[0], 2);
    base_pbwt[3].Update(&data[1], 2);

#if DEBUG_WAH
    uint32_t buf_wah_before = buf_wah[1].len;
#endif
    int ret1 = EncodeRLEBitmap2N2MM(2);
#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    std::shared_ptr<GenotypeDecompressorRLEBitmap> debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[1].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts1 = debug1->Decode2N2MM(buf_compress.data);
    // assert(alts1 == complex_pbwt[0].n_queue[1]);

    // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 1.
    for (int i = 0; i < n_samples; ++i) {
        // std::cerr << i << " -> " << (int)base_pbwt[2].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
        assert(base_pbwt[2].prev[i] == buf_compress.data[i]);
    }
    buf_wah_before = buf_wah[1].len;
#endif
    
    int ret2 = EncodeRLEBitmap2N2MM(3);
#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[1].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts2 = debug1->Decode2N2MM(buf_compress.data);
    // assert(alts1 == complex_pbwt[0].n_queue[1]);

    // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 1.
    for (int i = 0; i < n_samples; ++i) {
        // std::cerr << i << " -> " << (int)base_pbwt[2].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
        assert(base_pbwt[3].prev[i] == buf_compress.data[i]);
    }
    buf_wah_before = buf_wah[1].len;
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

    // memset(alts, 0, 256*sizeof(uint32_t));
    
    // for (int i = 0; i < 2*n_samples; ++i) {
    //     uint8_t val = BCF_UNPACK_GENOTYPE_GENERAL(data[i]);
    //     ++alts[val];
    // }

    // for (int i = 0; i < 256; ++i) {
    //     if (alts[i]) std::cerr << i << ":" << alts[i] << ",";
    // }
    // std::cerr << std::endl;

    // Todo: assert n_alleles < 16 AND n_alleles + missing + EOV < 16
    complex_pbwt[0].UpdateGeneral(&data[0], 2);
    complex_pbwt[1].UpdateGeneral(&data[1], 2);

#if DEBUG_WAH
    uint32_t buf_wah_before = buf_wah[2].len;
#endif
    EncodeRLEBitmap2NXM(0);
#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    std::shared_ptr<GenotypeDecompressorRLEBitmap> debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[2].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts1 = debug1->Decode2NXM(buf_compress.data);
    // assert(alts1 == complex_pbwt[0].n_queue[1]);

    // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 1.
    for (int i = 0; i < n_samples; ++i) {
        // std::cerr << i << " -> " << (int)complex_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
        assert(complex_pbwt[0].prev[i] == buf_compress.data[i]);
    }
    buf_wah_before = buf_wah[2].len;
#endif
    EncodeRLEBitmap2NXM(1);
#if DEBUG_WAH // Debug RLE-bitmap
    // Spawn copy of decompressor for RLE bitmaps.
    debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(&buf_wah[2].data[buf_wah_before], n_samples, 1, n_samples);
    // Decode RLE bitmaps and ascertain that the number of alts is the same as the input.
    int32_t alts2 = debug1->Decode2NXM(buf_compress.data);
    // assert(alts1 == complex_pbwt[0].n_queue[1]);

    // Iteratively ascertain that the input data is equal to the compressed->decompressed data for haplotype 1.
    for (int i = 0; i < n_samples; ++i) {
        // std::cerr << i << " -> " << (int)complex_pbwt[0].prev[i] << "==" << (int)buf_compress.data[i] << std::endl;
        assert(complex_pbwt[1].prev[i] == buf_compress.data[i]);
    }
#endif

    ++processed_lines_local;
    ++processed_lines;

    return 1;
}

int GenotypeCompressorRLEBitmap::Compress() {
    // RLE-bitmap hybrid.
    if (buf_wah[0].len) {
#if DEBUG_PBWT
        assert(debug_pbwt[0].len/n_samples == debug_pbwt[1].len/n_samples);
        uint32_t n_cycles = debug_pbwt[0].len/n_samples;

        std::shared_ptr<GenotypeDecompressorRLEBitmap> debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(buf_wah[0].data, buf_wah[0].len, 2*n_cycles, n_samples);
        PBWT pbwt1(n_samples, 2);
        PBWT pbwt2(n_samples, 2);

        uint8_t* out = new uint8_t[n_samples];
        uint32_t n_debug = 0;
        uint32_t debug_offset = 0;
        for (int i = 0; i < n_cycles; ++i) {
            int32_t alts1 = debug1->Decode2N2MC(out);
            pbwt1.ReverseUpdate(out);
            assert(alts1 == pbwt1.n_queue[1]);
            for (int j = 0; j < n_samples; ++j) {
                assert(BCF_UNPACK_GENOTYPE(debug_pbwt[0].buffer[debug_offset + j]) == pbwt1.prev[j]);
            }

            int32_t alts2 = debug1->Decode2N2MC(out);
            pbwt2.ReverseUpdate(out);
            assert(alts2 == pbwt2.n_queue[1]);
            for (int j = 0; j < n_samples; ++j) {
                assert(BCF_UNPACK_GENOTYPE(debug_pbwt[1].buffer[debug_offset + j]) == pbwt2.prev[j]);
            }

            debug_offset += n_samples;
            ++n_debug;
        }
        
        delete[] out;
        std::cerr << "debug=" << n_debug << std::endl;
#endif

        if (strategy == CompressionStrategy::ZSTD) {
            int praw = ZstdCompress(buf_wah[0].data, buf_wah[0].len,
                                    buf_compress.data, buf_compress.capacity(),
                                    20);

            bytes_out_zstd1 += praw;
            std::cerr << "Zstd->compressed: " << buf_wah[0].len << "->" << praw << " (" << (float)buf_wah[0].len/praw << ")" << std::endl;
            
            // Temp: emit data to cout.
            std::cout.write((char*)&processed_lines_local, sizeof(uint32_t));
            std::cout.write((char*)&buf_wah[0].len, sizeof(size_t));
            std::cout.write((char*)&praw, sizeof(size_t));
            std::cout.write((char*)buf_compress.data, praw);
            std::cout.flush(); 
        }

        if (strategy == CompressionStrategy::LZ4) {
            size_t plz4 = Lz4Compress(buf_wah[0].data, buf_wah[0].len,
                                      buf_compress.data, buf_compress.capacity(),
                                      9);
            
            // std::cerr << "writing=" << buf_wah[0].len << " and " << plz4 << std::endl;
            std::cerr << "0:Lz4->compressed: " << buf_wah[0].len << "->" << plz4 << " (" << (float)buf_wah[0].len/plz4 << ")" << std::endl;
            bytes_out_lz4 += plz4; 

            // Temp: emit data to cout.
            std::cout.write((char*)&processed_lines_local, sizeof(uint32_t));
            std::cout.write((char*)&buf_wah[0].len, sizeof(size_t));
            std::cout.write((char*)&plz4, sizeof(size_t));
            std::cout.write((char*)buf_compress.data, plz4);
            std::cout.flush();        
        }
        // Reset
        buf_wah[0].len = 0;
    }

    if (buf_wah[1].len) {
        if (strategy == CompressionStrategy::ZSTD) {
            int praw = ZstdCompress(buf_wah[1].data, buf_wah[1].len,
                                    buf_compress.data, buf_compress.capacity(),
                                    20);

            bytes_out_zstd1 += praw;
            std::cerr << "Zstd->compressed: " << buf_wah[1].len << "->" << praw << " (" << (float)buf_wah[1].len/praw << ")" << std::endl;
            
            // Temp: emit data to cout.
            std::cout.write((char*)&processed_lines_local, sizeof(uint32_t));
            std::cout.write((char*)&buf_wah[1].len, sizeof(size_t));
            std::cout.write((char*)&praw, sizeof(size_t));
            std::cout.write((char*)buf_compress.data, praw);
            std::cout.flush(); 
        }

        if (strategy == CompressionStrategy::LZ4) {
            size_t plz4 = Lz4Compress(buf_wah[1].data, buf_wah[1].len,
                                      buf_compress.data, buf_compress.capacity(),
                                      9);
            
            // std::cerr << "writing=" << buf_wah[0].len << " and " << plz4 << std::endl;
            std::cerr << "1:Lz4->compressed: " << buf_wah[1].len << "->" << plz4 << " (" << (float)buf_wah[1].len/plz4 << ")" << std::endl;
            bytes_out_lz4 += plz4; 

            // Temp: emit data to cout.
            std::cout.write((char*)&processed_lines_local, sizeof(uint32_t));
            std::cout.write((char*)&buf_wah[1].len, sizeof(size_t));
            std::cout.write((char*)&plz4, sizeof(size_t));
            std::cout.write((char*)buf_compress.data, plz4);
            std::cout.flush();        
        }
        // Reset
        buf_wah[1].len = 0;
    }

    if (buf_wah[2].len) {
#if DEBUG_PBWT
        assert(debug_pbwt[4].len/n_samples == debug_pbwt[5].len/n_samples);
        uint32_t n_cycles = debug_pbwt[4].len/n_samples;
        std::shared_ptr<GenotypeDecompressorRLEBitmap> debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(buf_wah[2].data, buf_wah[2].len, 2*n_cycles, n_samples);
        PBWT pbwt1(n_samples, 16);
        PBWT pbwt2(n_samples, 16);

        uint8_t* out = new uint8_t[n_samples];
        uint32_t n_debug = 0;
        uint32_t debug_offset = 0;

        for (int i = 0; i < n_cycles; ++i) {
            int32_t alts1 = debug1->Decode2NXM(out);
            pbwt1.ReverseUpdate(out);
            // assert(alts1 == pbwt1.n_queue[1]);
            for (int j = 0; j < n_samples; ++j) {
                // std::cerr << i << "->" << j << "/" << n_samples << ": " << (int)BCF_UNPACK_GENOTYPE_GENERAL(debug_pbwt[4].buffer[debug_offset + j]) << " == " << (int)pbwt1.prev[j] << std::endl;
                assert(BCF_UNPACK_GENOTYPE_GENERAL(debug_pbwt[4].buffer[debug_offset + j]) == pbwt1.prev[j]);
            }

            int32_t alts2 = debug1->Decode2NXM(out);
            pbwt2.ReverseUpdate(out);
            for (int j = 0; j < n_samples; ++j) {
                // std::cerr << i << "->" << j << "/" << n_samples << ": " << (int)BCF_UNPACK_GENOTYPE_GENERAL(debug_pbwt[5].buffer[debug_offset + j]) << " == " << (int)pbwt2.prev[j] << std::endl;
                assert(BCF_UNPACK_GENOTYPE_GENERAL(debug_pbwt[5].buffer[debug_offset + j]) == pbwt2.prev[j]);
            }

            debug_offset += n_samples;
            ++n_debug;
        }
        
        delete[] out;
        std::cerr << "debug=" << n_debug << std::endl;
#endif
        if (strategy == CompressionStrategy::ZSTD) {
            int praw = ZstdCompress(buf_wah[2].data, buf_wah[2].len,
                                    buf_compress.data, buf_compress.capacity(),
                                    20);

            bytes_out_zstd1 += praw;
            std::cerr << "Zstd->compressed: " << buf_wah[2].len << "->" << praw << " (" << (float)buf_wah[2].len/praw << ")" << std::endl;
            
            // Temp: emit data to cout.
            std::cout.write((char*)&processed_lines_local, sizeof(uint32_t));
            std::cout.write((char*)&buf_wah[2].len, sizeof(size_t));
            std::cout.write((char*)&praw, sizeof(size_t));
            std::cout.write((char*)buf_compress.data, praw);
            std::cout.flush(); 
        }

        if (strategy == CompressionStrategy::LZ4) {
            size_t plz4 = Lz4Compress(buf_wah[2].data, buf_wah[2].len,
                                      buf_compress.data, buf_compress.capacity(),
                                      9);
            
            // std::cerr << "writing=" << buf_wah[0].len << " and " << plz4 << std::endl;
            std::cerr << "2:Lz4->compressed: " << buf_wah[2].len << "->" << plz4 << " (" << (float)buf_wah[2].len/plz4 << ")" << std::endl;
            bytes_out_lz4 += plz4; 

            // Temp: emit data to cout.
            std::cout.write((char*)&processed_lines_local, sizeof(uint32_t));
            std::cout.write((char*)&buf_wah[2].len, sizeof(size_t));
            std::cout.write((char*)&plz4, sizeof(size_t));
            std::cout.write((char*)buf_compress.data, plz4);
            std::cout.flush();        
        }
        // Reset
        buf_wah[2].len = 0;
    }
    

#if DEBUG_PBWT
        for (int i = 0; i < 6; ++i) {
            if (debug_pbwt[i].len) {
                debug_pbwt[i].FinalizeDigest();
                std::cerr << "Digest-" << i << "=" << std::hex << (int)debug_pbwt[i].digest[0];
                for (int j = 1; j < 64; ++j) std::cerr << std::hex << (int)debug_pbwt[i].digest[j];
                std::cerr << std::dec << std::endl;
            }
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
    // bytes_out += p1 + p2 + p1E + p2E + p2X + p2X2 + praw + extra1 + extra2;
    // bytes_out += p1 + p2 + p1E + p2E + extra1 + extra2;
    // std::cerr << "[PROGRESS] " << bytes_in << "->" << bytes_out << " (" << (double)bytes_in/bytes_out << "-fold)" << std::endl;

    // Debug compression:
    std::cerr << "[PROGRESS ZSTD] " << bytes_in1 << "->" << bytes_out_zstd1 << " (" << (double)bytes_in1/bytes_out_zstd1 << "-fold)" << std::endl;
    std::cerr << "[PROGRESS LZ4] " << bytes_in1 << "->" << bytes_out_lz4 << " (" << (double)bytes_in1/bytes_out_lz4 << "-fold)" << std::endl;
}

int GenotypeCompressorRLEBitmap::EncodeRLEBitmap2N2MC(const int target) {
    assert(target == 0 || target == 1);

    bytes_in1 += n_samples;

    // RLE word -> 1 bit typing, 1 bit allele, word*8-2 bits for RLE
    const uint8_t* prev = base_pbwt[target].prev;
    uint8_t* dst = buf_wah[0].data;
    size_t& dst_len = buf_wah[0].len;

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

    ++gt_width[n_objects];

    return n_objects;
}

int GenotypeCompressorRLEBitmap::EncodeRLEBitmap2N2MM(const int target) {
    assert(target == 2 || target == 3);

    bytes_in1 += n_samples;

    // RLE word -> 1 bit typing, 2 bit allele, word*8-3 bits for RLE
    const uint8_t* prev = base_pbwt[target].prev;
    uint8_t* dst = buf_wah[1].data;
    size_t& dst_len = buf_wah[1].len;

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

    ++gt_width[n_objects];

    return n_objects;
}

int GenotypeCompressorRLEBitmap::EncodeRLEBitmap2NXM(const int target) {
    assert(target == 0 || target == 1);

    bytes_in1 += n_samples;

    // RLE word -> 1 bit typing, 4 bit allele, word*8-5 bits for RLE
    const uint8_t* prev = complex_pbwt[target].prev;
    uint8_t* dst = buf_wah[2].data;
    size_t& dst_len = buf_wah[2].len;

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

    ++gt_width[n_objects];

    return n_objects;
}