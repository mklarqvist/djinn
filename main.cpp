#include "getopt.h"

#include "range_coder.h"
#include "frequency_model.h"
#include "vcf_reader.h"

#include "pbwt.h"

// #include "model_ad.h"

// temp
#include <fstream>
#include <iostream>
#include <chrono>
// #include <thread>

#include "gt_compressor.h"
#include "gt_decompressor.h"

void ReadVcfGT (const std::string& filename) {
    std::unique_ptr<djinn::VcfReader> reader = djinn::VcfReader::FromFile(filename);
    if (reader.get() == nullptr) {
        std::cerr << "failed read" << std::endl;
        exit(1);
    }

    std::cerr << "Samples in VCF=" << reader->n_samples_ << std::endl;

    // AD: temp
    // FormatAlelleDepth fmt_ad(10000,reader->n_samples_,true);

    // GenotypeCompressorModelling gtperm(reader->n_samples_);
    // GenotypeCompressorRLEBitmap gtperm2(reader->n_samples_);
    
    djinn::GTCompressor gtcomp;
    // gtcomp.SetStrategy(djinn::GTCompressor::CompressionStrategy::CONTEXT_MODEL, reader->n_samples_);
    gtcomp.SetStrategy(djinn::GTCompressor::CompressionStrategy::RLE_BITMAP, reader->n_samples_);
    gtcomp.SetStrategy(djinn::GenotypeCompressor::CompressionStrategy::LZ4);
    // HaplotypeCompressor hcomp(2*reader->n_samples_);
    uint64_t n_lines = 0;
    // gtcomp.permute_pbwt = false;

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
        
        if (n_lines % (8196) == 0 && n_lines != 0) gtcomp.Compress();
        gtcomp.Encode(reader->bcf1_, reader->header_);
        ++n_lines;

        // hcomp.Encode(reader->bcf1_, reader->header_);
    }

    // Compress final.
    gtcomp.Compress();
    // hcomp.PrintSizes();

    // gtperm.Compress(); // Final
    // gtperm2.Compress(); // Final
    // std::cerr << "Final=" << gtperm.bytes_in << "->" << gtperm.bytes_out << " (" << (double)gtperm.bytes_in/gtperm.bytes_out << ")" << std::endl;
    // for (int j = 0; j < gtperm.gt_width.size(); ++j) {
    //     if (gtperm.gt_width[j] != 0) std::cout << j << "\t" << gtperm.gt_width[j] << std::endl;
    // }
}

int DecompressTest(const std::string& file, int type) {
    uint32_t n_samples = 32470;

    std::ifstream f(file, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) {
        std::cerr << "Failed to open" << std::endl;
        return 1;
    }
    uint64_t filesize = f.tellg();
    f.seekg(0);

    std::cerr << "filesize=" << filesize << "@" << f.tellg() << std::endl;

    size_t uncompressed_size = 0, compressed_size = 0, compressed_size2 = 0;
    size_t dest_capacity = 10000000;
    uint8_t* buffer     = new uint8_t[dest_capacity];
    uint8_t* out_buffer = new uint8_t[dest_capacity];
    int64_t n_variants  = 0;
    uint32_t n_variants_block  = 0;

    uint64_t tot_decomp = 0;
    uint8_t* out = new uint8_t[n_samples];
    uint8_t* vcf_buffer = new uint8_t[n_samples*8];
    djinn::PBWT pbwt1(n_samples, 2);
    djinn::PBWT pbwt2(n_samples, 2);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    if (type & 1) {
        uint8_t* buffer2     = new uint8_t[dest_capacity];
        uint8_t* out_buffer2 = new uint8_t[dest_capacity];
        uint8_t* buffer_partition1 = new uint8_t[dest_capacity];
        uint8_t* buffer_partition2 = new uint8_t[dest_capacity];
        size_t bitmap1 = 0;
        size_t compressed_size_bitmap1 = 0;
        size_t uncompressed_size_bitmap1 = 0;
        size_t bitmap2 = 0;
        size_t compressed_size_bitmap2 = 0;
        size_t uncompressed_size_bitmap2 = 0;

        while (f.good()) {
            std::cerr << f.tellg() << "/" << filesize << std::endl;
            f.read((char*)&n_variants_block, sizeof(uint32_t));
            f.read((char*)&uncompressed_size, sizeof(size_t));
            f.read((char*)&compressed_size, sizeof(size_t));
            f.read((char*)buffer, compressed_size);
            // std::cerr << "[STEP 1] " << uncompressed_size << "," << compressed_size << " -> " << f.tellg() << "/" << filesize << std::endl;
            f.read((char*)&uncompressed_size_bitmap1, sizeof(size_t));
            f.read((char*)&compressed_size_bitmap1, sizeof(size_t));
            f.read((char*)buffer_partition1, compressed_size_bitmap1);
            // std::cerr << "[STEP 2] " << uncompressed_size_bitmap1 << "," << compressed_size_bitmap1 << " -> " << f.tellg() << "/" << filesize << std::endl;

            // second
            f.read((char*)&n_variants_block, sizeof(uint32_t));
            f.read((char*)&uncompressed_size, sizeof(size_t));
            f.read((char*)&compressed_size2, sizeof(size_t));
            f.read((char*)buffer2, compressed_size2);
            // std::cerr << "[STEP 3] " << uncompressed_size << "," << compressed_size << " -> " << f.tellg() << "/" << filesize << std::endl;
            f.read((char*)&uncompressed_size_bitmap2, sizeof(size_t));
            f.read((char*)&compressed_size_bitmap2, sizeof(size_t));
            f.read((char*)buffer_partition2, compressed_size_bitmap2);
            // std::cerr << "[STEP 4] " << uncompressed_size_bitmap2 << "," << compressed_size_bitmap2 << " -> " << f.tellg() << "/" << filesize << std::endl;

            n_variants += n_variants_block;

            // Decode
            std::shared_ptr<djinn::GenotypeDecompressorContext> debug1 = 
                std::make_shared<djinn::GenotypeDecompressorContext>(
                    (uint8_t*)buffer, compressed_size, 
                    (uint8_t*)buffer_partition1, compressed_size_bitmap1, 
                    n_variants_block, n_samples, 2);

            std::shared_ptr<djinn::GenotypeDecompressorContext> debug2 = 
                std::make_shared<djinn::GenotypeDecompressorContext>(
                    (uint8_t*)buffer2, compressed_size2, 
                    (uint8_t*)buffer_partition2, compressed_size_bitmap2, 
                    n_variants_block, n_samples, 2);

            for (int i = 0; i < n_variants_block; ++i) {
                debug1->Decode(out);
                pbwt1.ReverseUpdate(out);
                debug2->Decode(out);
                pbwt2.ReverseUpdate(out);
            }

            if (f.tellg() == filesize) break;

            pbwt1.reset();
            pbwt2.reset();
        }

        delete[] buffer2;
        delete[] out_buffer2;
        delete[] buffer_partition1;
        delete[] buffer_partition2;
    }

    else {

        while (f.good()) {
            f.read((char*)&n_variants_block,  sizeof(uint32_t));
            f.read((char*)&uncompressed_size, sizeof(size_t));
            f.read((char*)&compressed_size,   sizeof(size_t));
            f.read((char*)buffer, compressed_size);
            std::cerr << uncompressed_size << "," << compressed_size << " -> " << f.tellg() << "/" << filesize << std::endl;
            int ret;
            if ((type >> 1) & 1)
                ret = djinn::Lz4Decompress(buffer, compressed_size, out_buffer, dest_capacity);
            else if ((type >> 2) & 1)
                ret = djinn::ZstdDecompress(buffer, compressed_size, out_buffer, dest_capacity);

            assert(ret == uncompressed_size);
            // n_variants += DebugRLEBitmap(out_buffer, uncompressed_size, n_samples, true);

            //
            std::shared_ptr<djinn::GenotypeDecompressorRLEBitmap> debug1 = std::make_shared<djinn::GenotypeDecompressorRLEBitmap>(out_buffer, uncompressed_size, 2*n_variants_block, n_samples);
            
            pbwt1.reset();
            pbwt2.reset();

            uint32_t n_debug = 0;
            // std::cerr << std::dec;
            while (true) {
                int32_t alts1 = debug1->Decode2N2MC(out);
                pbwt1.ReverseUpdate(out);
                assert(alts1 == pbwt1.n_queue[1]);
                int32_t alts2 = debug1->Decode2N2MC(out);
                pbwt2.ReverseUpdate(out);
                assert(alts2 == pbwt2.n_queue[1]);
                
                // std::cerr << "AC=" << alts1+alts2 << std::endl;

                vcf_buffer[0] = pbwt1.prev[0] + '0';
                vcf_buffer[1] = '|';
                vcf_buffer[2] = pbwt2.prev[0] + '0';
                
                int j = 3;
                for (int i = 1; i < n_samples; ++i, j += 4) {
                    vcf_buffer[j+0] = '\t';
                    vcf_buffer[j+1] = pbwt1.prev[i] + '0';
                    vcf_buffer[j+2] = '|';
                    vcf_buffer[j+3] = pbwt2.prev[i] + '0';
                }
                vcf_buffer[j++] = '\n';
                std::cout.write((char*)vcf_buffer, j);

                // int j = 0;
                // for (int i = 0; i < n_samples; ++i, j += 2) {
                //     vcf_buffer[j+0] = pbwt1.prev[i];
                //     vcf_buffer[j+1] = pbwt2.prev[i];
                // }
                
                // std::cout.write((char*)vcf_buffer, j);
                // std::cout.flush();

                if (++n_debug == n_variants_block) break;
            }
            assert(n_variants_block == n_debug);
            n_variants += n_variants_block;
        
            // std::cerr << "debug=" << n_debug << std::endl;
            //

            if (f.tellg() == filesize) break;
        }
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cerr << "[DECOMPRESS] Time elapsed " << time_span.count() << " ms " << tot_decomp << " variants=" << n_variants/2 << " (" << n_variants*n_samples << " genotypes: " << n_variants*n_samples/((double)time_span.count()/1000)/1e9 << " billion/s)" << std::endl;

    delete[] buffer;
    delete[] out_buffer;
    delete[] out;
    delete[] vcf_buffer;

    return EXIT_SUCCESS;
}


int main(int argc, char** argv) {
    int option_index = 0;
	static struct option long_options[] = {
		{"input",  required_argument, 0,  'i' },
        {"compress",  required_argument, 0,  'c' },
        {"decompress",  required_argument, 0,  'd' },
		{"zstd",   optional_argument, 0,  'z' },
		{"lz4",    optional_argument, 0,  'l' },
		{"modelling",optional_argument, 0,  'm' },
		{0,0,0,0}
	};

    std::string input;
    bool zstd = false;
    bool lz4 = true;
    bool context = false;
    bool compress = true;
    bool decompress = false;

    int c;
    while ((c = getopt_long(argc, argv, "i:zlcdm?", long_options, &option_index)) != -1){
		switch (c){
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			input = std::string(optarg);
			break;
		
        case 'z': zstd = true; lz4 = false; context = false; break;
        case 'l': zstd = false; lz4 = true; context = false; break;
        case 'm': zstd = false; lz4 = false; context = true; break;
        case 'c': compress = true; decompress = false; break;
        case 'd': compress = false; decompress = true; break;

		default:
			std::cerr << "Unrecognized option: " << (char)c << std::endl;
			return(1);
		}
	}

	if(input.length() == 0){
		std::cerr << "No input value specified..." << std::endl;
		return(1);
	}

    if (compress) {
        ReadVcfGT(input);
        return EXIT_SUCCESS;
    }

    if (decompress) {
        uint8_t type = (zstd << 2) | (lz4 << 1) | (context << 0);
        return (DecompressTest(input, type));
    }

    return EXIT_FAILURE;
}
