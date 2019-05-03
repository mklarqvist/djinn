// #include "r32x16b_avx2.h"
#include "range_coder.h"
#include "frequency_model.h"
#include "vcf_reader.h"
#include "pbwt.h"

// #include "model_ad.h"

#include <zstd.h>
#include <zstd_errors.h>

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

// int Lz4Decompress(const uint8_t* in, uint32_t n_in, uint8_t* out, uint32_t out_capacity) {
//     int ret = ZSTD_decompress(out, out_capacity, in, n_in);
//     return(ret);
// }


#include "gt_compressor.h"

void ReadVcfGT (const std::string& filename) {
    std::unique_ptr<tachyon::io::VcfReader> reader = tachyon::io::VcfReader::FromFile(filename,8);
    if (reader.get() == nullptr) {
        std::cerr << "failed read" << std::endl;
        exit(1);
    }

    // AD: temp
    // FormatAlelleDepth fmt_ad(10000,reader->n_samples_,true);

    GenotypeCompressor gtperm(reader->n_samples_);

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

         gtperm.Encode(reader->bcf1_, reader->header_);
    }

    // std::cerr << "total=" << in1 << "," << in2 << std::endl;

    // delete[] buf_in1;
    // delete[] buf_in2;
    // delete[] buf_out;

    gtperm.Compress(); // final
    std::cerr << "Final=" << gtperm.bytes_in << "->" << gtperm.bytes_out << " (" << (double)gtperm.bytes_in/gtperm.bytes_out << ")" << std::endl;
    for (int j = 0; j < gtperm.gt_width.size(); ++j) {
        if (gtperm.gt_width[j] != 0) std::cout << j << "\t" << gtperm.gt_width[j] << std::endl;
    }

}

int main(int argc, char** argv) {
    if (argc == 1) return(1);
    else {
        std::cerr << std::string(argv[1]) << std::endl;
        ReadVcfGT(std::string(argv[1]));
        return(0);
    }
}
