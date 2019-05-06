// #include "r32x16b_avx2.h"
#include "range_coder.h"
#include "frequency_model.h"
#include "vcf_reader.h"

#include "pbwt.h"

// #include "model_ad.h"

// temp
#include <fstream>
#include <iostream>
#include <chrono>

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

    // GenotypeCompressorModelling gtperm(reader->n_samples_);
    // GenotypeCompressorRLEBitmap gtperm2(reader->n_samples_);
    GTCompressor gtcomp;
    gtcomp.SetStrategy(GTCompressor::CompressionStrategy::CONTEXT_MODEL, reader->n_samples_);

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
        // gtperm.Encode(reader->bcf1_, reader->header_);
        // gtperm2.Encode(reader->bcf1_, reader->header_);
        gtcomp.Encode(reader->bcf1_, reader->header_);
    }

    // std::cerr << "total=" << in1 << "," << in2 << std::endl;

    // delete[] buf_in1;
    // delete[] buf_in2;
    // delete[] buf_out;

    // gtperm.Compress(); // Final
    // gtperm2.Compress(); // Final
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
