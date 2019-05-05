// #include "r32x16b_avx2.h"
#include "range_coder.h"
#include "frequency_model.h"
#include "vcf_reader.h"

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
                if (print) printf("AC=%d\n",n_alts);
                assert(n_alts <= n_samples);
                n_alts = 0;
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
                if (print) printf("AC=%d\n",n_alts);
                assert(n_alts <= n_samples);
                n_alts = 0;
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

    std::cerr << "Number of variants=" << n_variants << std::endl;

    return n_variants;
}

#include "pbwt.h"

// #include "model_ad.h"

#include <zstd.h>
#include <zstd_errors.h>

// temp
#include <fstream>
#include <iostream>
#include <chrono>

int ZstdCompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
    int ret = ZSTD_compress(out, out_capacity, in, n_in, c_level);
    return(ret);
}

int ZstdDecompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
    int ret = ZSTD_decompress(out, out_capacity, in, n_in);
    return(ret);
}

#include "lz4.h" // lz4
#include "lz4hc.h"

int Lz4Compress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity, const int32_t c_level = 1) {
    int compressed_data_size = LZ4_compress_HC((const char*)in, (char*)out, n_in, out_capacity, c_level);

    if (compressed_data_size < 0) {
        std::cerr << "A negative result from LZ4_compress_default indicates a failure trying to compress the data.  See exit code (echo $?) for value returned." << std::endl;
        exit(1);
    }
        
    if (compressed_data_size == 0) {
        std::cerr << "A result of 0 means compression worked, but was stopped because the destination buffer couldn't hold all the information." << std::endl;
        exit(1);
    }


    return(compressed_data_size);
}

int Lz4Decompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
    int32_t decompressed_size = LZ4_decompress_safe((char*)in, (char*)out, n_in, out_capacity);
    // int decompressed_data_size = LZ4_decompress((const char*)in, (char*)out, n_in, out_capacity);
    if (decompressed_size < 0) {
        std::cerr << "A negative result from LZ4_decompress_safe indicates a failure trying to decompress the data.  See exit code (echo $?) for value returned." << std::endl;
        exit(1);
    }
    if (decompressed_size == 0) {
        std::cerr << "I'm not sure this function can ever return 0.  Documentation in lz4.h doesn't indicate so.";
        exit(1);
    }


    return(decompressed_size);
}


#include "gt_compressor.h"

void ReadVcfGT (const std::string& filename) {
    std::unique_ptr<tachyon::io::VcfReader> reader = tachyon::io::VcfReader::FromFile(filename,8);
    if (reader.get() == nullptr) {
        std::cerr << "failed read" << std::endl;
        exit(1);
    }

    // std::cerr << "samples=" << reader->n_samples_ << std::endl;

    // AD: temp
    // FormatAlelleDepth fmt_ad(10000,reader->n_samples_,true);

    GenotypeCompressorModelling gtperm(reader->n_samples_);

    // temp
    // unsigned char* buf_in1 = new unsigned char[20000000];
    // unsigned char* buf_in2 = new unsigned char[20000000];
    // unsigned char* buf_out = new unsigned char[25000000];
    // unsigned int n_buf1_in = 0, n_buf2_in = 0;
    // unsigned int n_buf_out = 0;

    // uint64_t in1 = 0, in2 = 0;
    
    // While there are bcf records available.
    while (reader->Next()) {
        //const char* chrom = bcf_seqname(hr,line) ;
        //if (!p->chrom) p->chrom = strdup (chrom) ;
        //else if (strcmp (chrom, p->chrom)) break ;
        //int pos = line->pos; // bcf coordinates are 0-based
        //char *ref, *REF;
        //ref = REF = strdup(line->d.allele[0]);
        //while ( (*ref = toupper(*ref)) ) ++ref ;

        //
        // const uint8_t* gts = reader->bcf1_->d.fmt[0].p;
        // uint8_t gt1 = 0, gt2 = 0;
        // for (int i = 0; i + 16 <= 2*reader->n_samples_; i += 16) {
        //     for (int j = 0; j + 2 < 16; j += 2) {
        //         gt1 <<= 1;
        //         gt1 = (gt1 | ((gts[i+j] >> 1) - 1));
        //     }

        //     for (int j = 1; j + 2 < 16; j += 2) {
        //         gt2 <<= 1;
        //         gt2 = (gt2 | ((gts[i+j] >> 1) - 1));
        //     }

        //     buf_in1[n_buf1_in++] = gt1;
        //     buf_in2[n_buf2_in++] = gt2;
        //     // std::cerr << std::bitset<8>(gt1) << " " << std::bitset<8>(gt2) << std::endl;
        //     gt1 = 0, gt2 = 0;
        // }

        // if (n_buf1_in > (20000000 - 10000)) {
        //     unsigned int out_size = 25000000;
        //     std::cerr << "here=" << n_buf1_in << std::endl;
        //     memset(buf_out, 0, 25000000);
        //     uint8_t* ret = rans_compress_O0_32x16(buf_in1, n_buf1_in, buf_out, &out_size);
        //     // int ret = rans_test(buf_in1, n_buf1_in, buf_out, &out_size);
        //     // std::cerr << "ret=" << ret << std::endl;
        //     std::cerr << "cleft= " << n_buf1_in*8 << "->" << out_size << " (" << (float)(n_buf1_in*8)/out_size << "-fold)" << std::endl;
        //     in1 += out_size;
        //     out_size = 25000000;
        //     memset(buf_out, 0, 25000000);
        //     ret = rans_compress_O0_32x16(buf_in2, n_buf2_in, buf_out, &out_size);
        //     std::cerr << "cright= " << n_buf2_in*8 << "->" << out_size << " (" << (float)(n_buf2_in*8)/out_size << "-fold)" << std::endl;
        //     in2 += out_size;

        //     n_buf1_in = 0;
        //     n_buf2_in = 0;
        // }
        //

        // std::cerr << reader->bcf1_->pos+1 << std::endl;
        gtperm.Encode(reader->bcf1_, reader->header_);
    }

    // std::cerr << "total=" << in1 << "," << in2 << std::endl;

    // delete[] buf_in1;
    // delete[] buf_in2;
    // delete[] buf_out;

    gtperm.Compress(); // Final
    // std::cerr << "Final=" << gtperm.bytes_in << "->" << gtperm.bytes_out << " (" << (double)gtperm.bytes_in/gtperm.bytes_out << ")" << std::endl;
    // for (int j = 0; j < gtperm.gt_width.size(); ++j) {
    //     if (gtperm.gt_width[j] != 0) std::cout << j << "\t" << gtperm.gt_width[j] << std::endl;
    // }
}


int main(int argc, char** argv) {
#if 0
    std::ifstream f("/media/mdrk/NVMe/gt_hrc11.bin2.lz4", std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) {
        std::cerr << "Failed to open" << std::endl;
        return 1;
    }
    uint64_t filesize = f.tellg();
    f.seekg(0);

    size_t uncompressed_size = 0, compressed_size = 0;
    size_t dest_capacity = 10000000;
    uint8_t* buffer     = new uint8_t[dest_capacity];
    uint8_t* out_buffer = new uint8_t[dest_capacity];
    int64_t n_variants  = 0;

    uint64_t tot_decomp = 0;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    while (f.good()) {
        f.read((char*)&uncompressed_size, sizeof(size_t));
        f.read((char*)&compressed_size, sizeof(size_t));
        f.read((char*)buffer, compressed_size);
        std::cerr << uncompressed_size << "," << compressed_size << " -> " << f.tellg() << "/" << filesize << std::endl;
        int ret = Lz4Decompress(buffer, compressed_size, out_buffer, dest_capacity);
        assert(ret == uncompressed_size);
        n_variants += DebugRLEBitmap(out_buffer, uncompressed_size, 32470);

        if (f.tellg() == filesize) break;
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[LZ4] Time elapsed " << time_span.count() << " ms " << tot_decomp << " variants=" << n_variants/2 << " (" << n_variants*32470 << " genotypes: " << n_variants*32470/((double)time_span.count()/1000)/1e9 << " billion/s)" << std::endl;

    delete[] buffer;
    delete[] out_buffer;

    return 0;
#endif

    if (argc == 1) return(1);
    else {
        std::cerr << std::string(argv[1]) << std::endl;
        ReadVcfGT(std::string(argv[1]));
        return(0);
    }
}
