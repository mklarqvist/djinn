#include <cstdlib>
#include <cmath>
#include <errno.h>

#include <memory>
#include <iostream>

#include <bitset>
#include <vector>

#include "range_coder.h"
#include "frequency_model.h"
#include "vcf_reader.h"
#include "pbwt.h"

#include "model_ad.h"

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
    
    // While there are bcf records available.
    while (reader->Next()) {
        //const char* chrom = bcf_seqname(hr,line) ;
        //if (!p->chrom) p->chrom = strdup (chrom) ;
        //else if (strcmp (chrom, p->chrom)) break ;
        //int pos = line->pos; // bcf coordinates are 0-based
        //char *ref, *REF;
        //ref = REF = strdup(line->d.allele[0]);
        //while ( (*ref = toupper(*ref)) ) ++ref ;

        // if (reader->bcf1_->n_allele == 2) fmt_ad.AddRecord(fmt);
         gtperm.Encode(reader->bcf1_, reader->header_);
    }

    gtperm.Compress(); // final
    std::cerr << "Final=" << gtperm.bytes_in << "->" << gtperm.bytes_out << " (" << (double)gtperm.bytes_in/gtperm.bytes_out << ")" << std::endl;

}

int main(int argc, char** argv) {
    if (argc == 1) return(1);
    else {
        std::cerr << std::string(argv[1]) << std::endl;
        ReadVcfGT(std::string(argv[1]));
        return(0);
    }
}
