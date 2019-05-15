#include "getopt.h"

#include "vcf_reader.h"
#include "gt_compressor.h"
#include "gt_decompressor.h"

// temp
#include <fstream>
#include <iostream>
#include <chrono>

int ReadVcfGT(const std::string& filename, int type, bool permute = true) {
    std::unique_ptr<djinn::VcfReader> reader = djinn::VcfReader::FromFile(filename);
    if (reader.get() == nullptr) {
        std::cerr << "failed read" << std::endl;
        return -1;
    }

    std::cerr << "Samples in VCF=" << reader->n_samples_ << std::endl;

    djinn::djinn_hdr_t hdr;
    hdr.base_ploidy = 2;
    hdr.n_samples = reader->n_samples_;
    hdr.version[0] = 0; hdr.version[1] = 1; hdr.version[2] = 0;
    std::cout.write((char*)hdr.version, sizeof(uint8_t)*3);
    std::cout.write((char*)&hdr.base_ploidy, sizeof(uint8_t));
    std::cout.write((char*)&hdr.n_samples, sizeof(int64_t));
    std::cout.flush();

    djinn::GTCompressor gtcomp;
    if (type & 1) {
        gtcomp.SetStrategy(djinn::GTCompressor::CompressionStrategy::CONTEXT_MODEL, reader->n_samples_);
    }
    else {
        gtcomp.SetStrategy(djinn::GTCompressor::CompressionStrategy::RLE_BITMAP, reader->n_samples_);
        if (type & 2) gtcomp.SetStrategy(djinn::GenotypeCompressor::CompressionStrategy::LZ4);
        else if (type & 4) gtcomp.SetStrategy(djinn::GenotypeCompressor::CompressionStrategy::ZSTD);
        else {
            std::cerr << "error" << std::endl;
            return -2;
        }
    }
    gtcomp.SetPermutePbwt(permute);
    djinn::djinn_block_t* block = nullptr; // output djinn block
    
    uint64_t n_lines = 0;
    uint32_t n_blocks = 8196; // Number of variants per data block.
    uint8_t* output_data = new uint8_t[n_blocks*reader->n_samples_*2];

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
        
        if (n_lines % n_blocks == 0 && n_lines != 0) {
            gtcomp.Compress(block);
            // block->Serialize(std::cout);
            int ret = block->Serialize(output_data);
            std::cerr << "Data=" << ret << std::endl;
            std::cout.write((char*)output_data, ret);
        }
        gtcomp.Encode(reader->bcf1_, reader->header_);
        ++n_lines;
    }

    // Compress final.
    gtcomp.Compress(block);
    // block->Serialize(std::cout);
    int ret = block->Serialize(output_data);
    std::cerr << "Data=" << ret << std::endl;
    std::cout.write((char*)output_data, ret);
    std::cout.flush();
    delete block;

    return n_lines;
}

int DecompressExample(const std::string& file, int type) {
    std::ifstream f(file, std::ios::in | std::ios::binary | std::ios::ate);
    if (f.good() == false) {
        std::cerr << "Failed to open" << std::endl;
        return 1;
    }
    uint64_t filesize = f.tellg();
    f.seekg(0);

    // std::cerr << "filesize=" << filesize << "@" << f.tellg() << std::endl;

    djinn::djinn_hdr_t hdr;
    f.read((char*)hdr.version, sizeof(uint8_t)*3);
    f.read((char*)&hdr.base_ploidy, sizeof(uint8_t));
    f.read((char*)&hdr.n_samples, sizeof(int64_t)); 

    std::cerr << hdr.n_samples << std::endl;

    size_t dest_capacity = 10000000;

    uint8_t** buffers = new uint8_t*[6];
    for (int i = 0; i < 6; ++i) {
        buffers[i] = new uint8_t[dest_capacity];
    }
    uint8_t* line = new uint8_t[hdr.n_samples*2];
    uint8_t* vcf_buffer = new uint8_t[hdr.n_samples*8];

    while (f.good()) {
        djinn::djinn_block_t block;
        // std::cerr << f.tellg() << "/" << filesize << std::endl;
        int ret = block.Deserialize(f);
        if (block.type == djinn::djinn_block_t::BlockType::WAH) {
            djinn::djinn_wah_block_t* bl = (djinn::djinn_wah_block_t*)&block;
            djinn::djinn_wah_t* d = (djinn::djinn_wah_t*)bl->data;
            for (int i = 0; i < 6; ++i) {
                if (d->wah_models[i].vptr_len) {
                    int ret = 0;
                    if ((type >> 1) & 1)
                        ret = djinn::Lz4Decompress(d->wah_models[i].vptr, d->wah_models[i].vptr_len, buffers[i], dest_capacity);
                    else if ((type >> 2) & 1)
                        ret = djinn::ZstdDecompress(d->wah_models[i].vptr, d->wah_models[i].vptr_len, buffers[i], dest_capacity);
                    // std::cerr << "ret=" << ret << std::endl;
                }
            }
            djinn::djinn_wah_ctrl_t* ctrl = (djinn::djinn_wah_ctrl_t*)&bl->ctrl;
            std::cerr << "using pbwt=" << ctrl->pbwt << " " << f.tellg() << "/" << filesize << std::endl;
            
            // Unpack WAH 2N2MC
            std::shared_ptr<djinn::GenotypeDecompressorRLEBitmap> debug1 = 
                std::make_shared<djinn::GenotypeDecompressorRLEBitmap>(
                    buffers[0], d->wah_models[0].n, 
                    d->wah_models[0].n_v, hdr.n_samples);

            std::shared_ptr<djinn::GenotypeDecompressorRLEBitmap> debug2 = 
                std::make_shared<djinn::GenotypeDecompressorRLEBitmap>(
                    buffers[1], d->wah_models[1].n, 
                    d->wah_models[1].n_v, hdr.n_samples);

            if (ctrl->pbwt) { // PBWT
                debug1->InitPbwt(2);
                debug2->InitPbwt(2);

                 assert(d->wah_models[0].n_v == d->wah_models[1].n_v);
                for (int i = 0; i < d->wah_models[0].n_v; ++i) {
                    int ret1 = debug1->Decode2N2MC(line);
                    // debug1->pbwt->ReverseUpdate(line);
                    int ret2 = debug2->Decode2N2MC(line);
                    // debug2->pbwt->ReverseUpdate(line);
                    std::cout << "AC=" << ret1 + ret2 << std::endl;
                }
                
            } 
            else { // no PBWT
                assert(d->wah_models[0].n_v == d->wah_models[1].n_v);
                for (int i = 0; i < d->wah_models[0].n_v; ++i) {
                    debug1->Decode2N2MC(&line[0], 2);
                    debug2->Decode2N2MC(&line[1], 2);
                    // 'line' is now the correct output string
                    // encoded as 0, 1, ..., n_alts with 14 and 15 encoding
                    // for missing and EOV, respectively.
                    //
                    // Emit VCF line. Ignores missing and EOV in this example.
                    // To cover these cases then map 14->'.' and 15->EOV
                    vcf_buffer[0] = line[0] + '0';
                    vcf_buffer[1] = '|';
                    vcf_buffer[2] = line[1] + '0';
                    
                    int j = 3;
                    for (int i = 2; i < 2*hdr.n_samples; i += 2, j += 4) {
                        vcf_buffer[j+0] = '\t';
                        vcf_buffer[j+1] = line[i+0] + '0';
                        vcf_buffer[j+2] = '|';
                        vcf_buffer[j+3] = line[i+1] + '0';
                    }
                    vcf_buffer[j++] = '\n';
                    std::cout.write((char*)vcf_buffer, j);
                }
            } // no PBWT

        } // end type WAH

        std::cerr << "[AFTER] " << f.tellg() << "/" << filesize << std::endl;
        if (ret <= 0) exit(1);
        // exit(1);


        if (f.tellg() == filesize) break;
    }

    for (int i = 0; i < 6; ++i) delete[] buffers[i];
    delete[] buffers;
    delete[] line;
    delete[] vcf_buffer;

    return 1;
}

void usage() {
    printf("djinn\n");
    printf("   -i STRING input file (vcf,vcf.gz,bcf, or ubcf (required))\n");
    printf("   -c BOOL   compress file\n");
    printf("   -d BOOL   decompress file\n");
    printf("   -z BOOL   compress with RLE-hybrid + ZSTD-19\n");
    printf("   -l BOOL   compress with RLE-hybrid + LZ4-HC-9\n");
    printf("   -m BOOL   compress with context modelling\n");
    printf("   -p BOOL   permute data with PBWT\n");
    printf("   -P BOOL   do not permute data with PBWT\n\n");
    printf("Examples:\n");
    printf("  djinn -clpi file.bcf > /dev/null\n");
    printf("  djinn -czPi file.bcf > /dev/null\n");
    printf("  djinn -cmi file.bcf > /dev/null\n\n");
}

int main(int argc, char** argv) {
    if (argc == 1) {
        usage();
        return EXIT_SUCCESS;
    }

    int option_index = 0;
	static struct option long_options[] = {
		{"input",  required_argument, 0,  'i' },
        {"compress",  required_argument, 0,  'c' },
        {"decompress",  required_argument, 0,  'd' },
		{"zstd",   optional_argument, 0,  'z' },
		{"lz4",    optional_argument, 0,  'l' },
		{"modelling",optional_argument, 0,  'm' },
        {"permute",  optional_argument, 0,  'p' },
        {"no-permute",  optional_argument, 0,  'P' },
		{0,0,0,0}
	};

    std::string input;
    bool zstd = false;
    bool lz4 = true;
    bool context = false;
    bool compress = true;
    bool decompress = false;
    bool permute = true;

    int c;
    while ((c = getopt_long(argc, argv, "i:zlcdmpP?", long_options, &option_index)) != -1){
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
        case 'p': permute = true; break;
        case 'P': permute = false; break;

		default:
			std::cerr << "Unrecognized option: " << (char)c << std::endl;
			return 1;
		}
	}

	if(input.length() == 0){
		std::cerr << "No input value specified..." << std::endl;
		return 1;
	}

    assert(zstd || lz4 || context);
    int type = (zstd << 2) | (lz4 << 1) | (context << 0);

    if (compress) {
        return(ReadVcfGT(input, type, permute));
    }

    if (decompress) {
        return(DecompressExample(input, type));
    }

    return EXIT_FAILURE;
}
