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
#include "gt_decompressor.h"

void ReadVcfGT (const std::string& filename) {
    std::unique_ptr<tachyon::io::VcfReader> reader = tachyon::io::VcfReader::FromFile(filename,8);
    if (reader.get() == nullptr) {
        std::cerr << "failed read" << std::endl;
        exit(1);
    }

    std::cerr << "Samples in VCF=" << reader->n_samples_ << std::endl;

    // AD: temp
    // FormatAlelleDepth fmt_ad(10000,reader->n_samples_,true);

    // GenotypeCompressorModelling gtperm(reader->n_samples_);
    // GenotypeCompressorRLEBitmap gtperm2(reader->n_samples_);
    GTCompressor gtcomp;
    gtcomp.SetStrategy(GTCompressor::CompressionStrategy::CONTEXT_MODEL, reader->n_samples_);
    // gtcomp.SetStrategy(GenotypeCompressor::CompressionStrategy::LZ4);
    
    // While there are bcf records available.
    while (reader->Next()) {
        //const char* chrom = bcf_seqname(hr,line) ;
        //if (!p->chrom) p->chrom = strdup (chrom) ;
        //else if (strcmp (chrom, p->chrom)) break ;
        //int pos = line->pos; // bcf coordinates are 0-based
        //char *ref, *REF;
        //ref = REF = strdup(line->d.allele[0]);
        //while ( (*ref = toupper(*ref)) ) ++ref ;

        // std::cerr << reader->bcf1_->pos+1 << std::endl;
        // gtperm.Encode(reader->bcf1_, reader->header_);
        // gtperm2.Encode(reader->bcf1_, reader->header_);
        gtcomp.Encode(reader->bcf1_, reader->header_);
    }

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
    uint32_t n_variants_block  = 0;

    uint64_t tot_decomp = 0;
    uint8_t* out = new uint8_t[32470];
    uint8_t* vcf_buffer = new uint8_t[32470*8];
    PBWT pbwt1(32470, 2);
    PBWT pbwt2(32470, 2);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    while (f.good()) {
        f.read((char*)&n_variants_block,  sizeof(uint32_t));
        f.read((char*)&uncompressed_size, sizeof(size_t));
        f.read((char*)&compressed_size,   sizeof(size_t));
        f.read((char*)buffer, compressed_size);
        std::cerr << uncompressed_size << "," << compressed_size << " -> " << f.tellg() << "/" << filesize << std::endl;
        int ret = Lz4Decompress(buffer, compressed_size, out_buffer, dest_capacity);
        assert(ret == uncompressed_size);
        // n_variants += DebugRLEBitmap(out_buffer, uncompressed_size, 32470, true);

        //
        std::shared_ptr<GenotypeDecompressorRLEBitmap> debug1 = std::make_shared<GenotypeDecompressorRLEBitmap>(out_buffer, uncompressed_size, 2*n_variants_block, 32470);
        
        pbwt1.reset();
        pbwt2.reset();

        uint32_t n_debug = 0;
        // std::cerr << std::dec;
        while (true) {
            int32_t alts1 = debug1->DecodeRLEBitmap(out);
            if (alts1 < 0) break;
            pbwt1.ReverseUpdate(out);
            assert(alts1 == pbwt1.n_queue[1]);
            int32_t alts2 = debug1->DecodeRLEBitmap(out);
            pbwt2.ReverseUpdate(out);
            assert(alts2 == pbwt2.n_queue[1]);
            ++n_variants; ++n_debug;
            // std::cerr << "AC=" << alts1+alts2 << std::endl;

            vcf_buffer[0] = pbwt1.prev[0] + '0';
            vcf_buffer[1] = '|';
            vcf_buffer[2] = pbwt2.prev[0] + '0';
            
            int j = 3;
            for (int i = 1; i < 32470; ++i, j += 4) {
                vcf_buffer[j+0] = '\t';
                vcf_buffer[j+1] = pbwt1.prev[i] + '0';
                vcf_buffer[j+2] = '|';
                vcf_buffer[j+3] = pbwt2.prev[i] + '0';
            }
            vcf_buffer[j++] = '\n';
            std::cout.write((char*)vcf_buffer, j);

            // int j = 0;
            // for (int i = 0; i < 32470; ++i, j += 2) {
            //     vcf_buffer[j+0] = pbwt1.prev[i];
            //     vcf_buffer[j+1] = pbwt2.prev[i];
            // }
            // std::cout.write((char*)vcf_buffer, j);

        }
        assert(n_variants_block == n_debug);
       
        // std::cerr << "debug=" << n_debug << std::endl;
        //

        if (f.tellg() == filesize) break;
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[LZ4] Time elapsed " << time_span.count() << " ms " << tot_decomp << " variants=" << n_variants/2 << " (" << n_variants*32470 << " genotypes: " << n_variants*32470/((double)time_span.count()/1000)/1e9 << " billion/s)" << std::endl;

    delete[] buffer;
    delete[] out_buffer;
    delete[] out;
    delete[] vcf_buffer;

    return EXIT_SUCCESS;
#endif

    if (argc == 1) return(EXIT_FAILURE);
    else {
        std::cerr << std::string(argv[1]) << std::endl;
        ReadVcfGT(std::string(argv[1]));
        return(EXIT_SUCCESS);
    }
}
