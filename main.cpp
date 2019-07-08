#include "getopt.h"

#include "vcf_reader.h"
#include "ctx_model.h"
#include "ewah_model.h"

// temp
#include <algorithm>//sort
#include <fstream>
#include <iostream>
#include <chrono>//time
#include <random>

// Definition for microsecond timer.
typedef std::chrono::high_resolution_clock::time_point clockdef;

std::string MicroPrettyString(uint32_t ms) {
    std::string time;
    if (ms < 1000) return time + std::to_string(ms) + "us";
    else return time + std::to_string(ms/1000) + "ms";
}

int ReadVcfGT(const std::string& filename, int type, bool permute = true) {
    std::unique_ptr<djinn::VcfReader> reader = djinn::VcfReader::FromFile(filename);
    if (reader.get() == nullptr) {
        std::cerr << "failed read" << std::endl;
        return -1;
    }

    std::cerr << "Samples in VCF=" << reader->n_samples_ << std::endl;
    
    uint64_t n_lines = 0;
    uint32_t nv_blocks = 8192; // Number of variants per data block.
    uint32_t n_blocks = 0;

    // djinn::HaplotypeCompressor hc(reader->n_samples_);
    // uint32_t hc_lines = 0;

    // While there are bcf records available.
    clockdef t1 = std::chrono::high_resolution_clock::now();

#if 0
    // setup random
    int64_t n_fake_samples = 50000;
    uint32_t n_fake_sites  = 25000000;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<uint32_t> dis(0, n_fake_samples-1);
    std::uniform_int_distribution<uint32_t> freq_dis(0, 1000);
    uint8_t* rand_vec = new uint8_t[n_fake_samples];
    uint8_t* decode_buf = new uint8_t[10000000];

    // std::string temp_file = "/Users/Mivagallery/Downloads/djn_debug.bin";
    std::string temp_file = "/media/mdrk/08dcb478-5359-41f4-97c8-469190c8a034/djn_debug.bin";

    bool reset_models = false;
    // djinn::djinn_ctx_model djn_ctx;
    djinn::djinn_ctx_model* djn_ctx = new djinn::djinn_ctx_model;
    djn_ctx->StartEncoding(permute, reset_models);

    uint64_t data_in = 0, data_in_vcf = 0;
    uint64_t ewah_out = 0, ctx_out = 0;
    uint64_t ctx_out_progress = 0;

    clockdef line_t1 = std::chrono::high_resolution_clock::now();
    clockdef line_t2 = std::chrono::high_resolution_clock::now();
    uint64_t line_time = 0;
    
    for(int i = 0; i < n_fake_sites; ++i) {
        // std::cerr << "i=" << i << "/" << n_fake_sites << std::endl;
        if (n_lines % nv_blocks == 0 && n_lines != 0) {
            djn_ctx->FinishEncoding();
            int decode_ret2 = djn_ctx->Serialize(std::cout);
            djn_ctx->StartEncoding(permute, reset_models);
            ++n_blocks;
            ctx_out += decode_ret2;

            std::cerr << "[PROGRESS][CTX] In uBCF: " << data_in << "->" << ctx_out 
                << " (" << (double)data_in/ctx_out << "-fold) In VCF: " << data_in_vcf << "->" << ctx_out 
                << " (" << (double)data_in_vcf/ctx_out << "-fold)" << std::endl;

            // auto time_span = std::chrono::duration_cast<std::chrono::microseconds>(line_t2 - line_t1);
            std::cerr << "[PROGRESS] Time per line: " << (double)line_time/n_lines << "us" << std::endl;
            line_time = 0;
        }

        memset(rand_vec, 0, n_fake_samples);
        const int n_alts = freq_dis(gen);
        for (int j = 0; j < n_alts; ++j) {
            rand_vec[dis(gen)] = 1;
        }
        line_t1 = std::chrono::high_resolution_clock::now();
        int ret = djn_ctx->Encode(rand_vec, n_fake_samples, 2, 2);
        line_t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<uint64_t, std::micro> time_span = std::chrono::duration_cast<std::chrono::microseconds>(line_t2 - line_t1);
        line_time += time_span.count();

        ++n_lines;
        data_in += n_fake_samples;
        data_in_vcf += 2*n_fake_samples - 1;

        ctx_out_progress = ctx_out + djn_ctx->GetCurrentSize();
        assert(ctx_out_progress >= 0);
        std::cerr << "line-" << (n_lines%nv_blocks) << " alts=" << n_alts << " time=" << MicroPrettyString(time_span.count()) << ". VCF: " << data_in_vcf << "->" << ctx_out_progress 
                << " (" << (double)data_in_vcf/ctx_out_progress << "-fold)" << std::endl;
    
    }
    delete[] rand_vec;
#else
    std::cerr << "============== ENCODING ===============" << std::endl;

    bool reset_models = true;
    djinn::djinn_model* djn_ctx = new djinn::djinn_ctx_model;
    djn_ctx->StartEncoding(permute, reset_models);

    djinn::djinn_ewah_model djn_ewah;
    djn_ewah.StartEncoding(permute, reset_models);
    // djn_ewah.codec = djinn::CompressionStrategy::LZ4;
    
    std::string temp_file = "/media/mdrk/08dcb478-5359-41f4-97c8-469190c8a034/djn_debug.bin";
    // std::string temp_file = "/Users/Mivagallery/Downloads/djn_debug.bin";
    std::ofstream test_write(temp_file, std::ios::out | std::ios::binary);
    if (test_write.good() == false) {
        std::cerr << "could not open outfile handle" << std::endl;
        return -3;
    }
    uint8_t* decode_buf = new uint8_t[10000000];

    uint64_t data_in = 0, data_in_vcf = 0;
    uint64_t ewah_out = 0, ctx_out = 0;

    while (reader->Next()) {
        if (n_lines % nv_blocks == 0 && n_lines != 0) {
            djn_ctx->FinishEncoding();
            int decode_ret = djn_ctx->Serialize(decode_buf);
            // int decode_ret2 = djn_ctx->Serialize(std::cout);
            std::cerr << "[WRITING CTX] " << decode_ret << "b" << std::endl;
            test_write.write((char*)decode_buf, decode_ret);
            ++n_blocks;
            ctx_out += decode_ret;
            // ctx_out += decode_ret2;

            // EWAH
            int ret_ewah = djn_ewah.FinishEncoding();
            decode_ret = djn_ewah.Serialize(decode_buf);
            std::cerr << "[WRITING EWAH] " << decode_ret << "b" << std::endl;
            // test_write.write((char*)decode_buf, decode_ret);

            ewah_out += decode_ret;
            std::cerr << "[PROGRESS][CTX] In uBCF: " << data_in << "->" << ctx_out 
                << " (" << (double)data_in/ctx_out << "-fold) In VCF: " << data_in_vcf << "->" << ctx_out 
                << " (" << (double)data_in_vcf/ctx_out << "-fold)" << std::endl;
            std::cerr << "[PROGRESS][EWAH] In uBCF: " << data_in << "->" << ewah_out 
                << " (" << (double)data_in/ewah_out << "-fold) In VCF: " << data_in_vcf << "->" << ewah_out 
                << " (" << (double)data_in_vcf/ewah_out << "-fold)" << std::endl;

            djn_ewah.StartEncoding(permute, reset_models);
            djn_ctx->StartEncoding(permute, reset_models);
        }

        if (reader->bcf1_ == NULL)   return -1;
        if (reader->header_ == NULL) return -2;

        const bcf_fmt_t* fmt = bcf_get_fmt(reader->header_, reader->bcf1_, "GT");
        if (fmt == NULL) return 0;
        
        data_in += fmt->p_len;
        data_in_vcf += 2*fmt->p_len - 1;
        int ret = djn_ctx->EncodeBcf(fmt->p, fmt->p_len, fmt->n, reader->bcf1_->n_allele);
        assert(ret>0);
        int ret2 = djn_ewah.EncodeBcf(fmt->p, fmt->p_len, fmt->n, reader->bcf1_->n_allele);
        assert(ret2>0);

        ++n_lines;
    }

    // Compress final.
    djn_ctx->FinishEncoding();
    int decode_ret = djn_ctx->Serialize(decode_buf);
    // int decode_ret2 = djn_ctx->Serialize(std::cout);
    std::cerr << "[WRITING CTX] " << decode_ret << "b" << std::endl;
    assert(decode_ret > 0);
    test_write.write((char*)decode_buf, decode_ret);
    ++n_blocks;
    ctx_out += decode_ret;
    // ctx_out += decode_ret2;

     // EWAH
    int ret_ewah = djn_ewah.FinishEncoding();
    assert(ret_ewah > 0);
    decode_ret = djn_ewah.Serialize(decode_buf);
    std::cerr << "[WRITING EWAH] " << decode_ret << "b" << std::endl;
    // djn_ewah.StartEncoding(permute, reset_models);
    // test_write.write((char*)decode_buf, decode_ret);
    ewah_out += decode_ret;

    std::cerr << "[PROGRESS][CTX] In uBCF: " << data_in << "->" << ctx_out 
        << " (" << (double)data_in/ctx_out << "-fold) In VCF: " << data_in_vcf << "->" << ctx_out 
        << " (" << (double)data_in_vcf/ctx_out << "-fold)" << std::endl;
    std::cerr << "[PROGRESS][EWAH] In uBCF: " << data_in << "->" << ewah_out 
        << " (" << (double)data_in/ewah_out << "-fold) In VCF: " << data_in_vcf << "->" << ewah_out 
        << " (" << (double)data_in_vcf/ewah_out << "-fold)" << std::endl;    

    // Close handle.
    test_write.flush();
    test_write.close();
#endif
    delete djn_ctx;

    // Decode test
    std::cerr << "============== DECODING ===============" << std::endl;

    std::ifstream test_read(temp_file, std::ios::in | std::ios::binary | std::ios::ate);
    if (test_read.good() == false) {
        std::cerr << "could not open infile handle" << std::endl;
        return -3;
    }
    uint64_t filesize = test_read.tellg();
    test_read.seekg(0);

    djinn::djinn_model* djn_decode = new djinn::djinn_ctx_model;
    // djinn::djinn_ctx_model djn_ctx_decode;
    uint8_t* vcf_out_buffer = new uint8_t[4*reader->n_samples_+65536];
    uint32_t len_vcf = 0;
    djinn::djinn_ewah_model djn_ewah_decode;

    djinn::djinn_variant_t* variant = nullptr;

    clockdef t1_decode = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_blocks; ++i) {
        std::cerr << "============== BLOCK-" << i << "===============" << std::endl;
        // Deserialize data and construct the djinn_ctx_model instance. This approach
        // involves copying the data internally and can be slower compared to reading
        // and constructing directly from a file stream.
        // int decode_ctx_ret = djn_ctx_decode.Deserialize(decode_buf);
        int decode_ctx_ret = djn_decode->Deserialize(test_read);

        // Initiate decoding.
        djn_decode->StartDecoding();

        clockdef t1 = std::chrono::high_resolution_clock::now();
        // Cycle over variants in the block.
        for (int i = 0; i < djn_decode->n_variants; ++i) {
            // std::cerr << "Decoding " << i << "/" << djn_decode->n_variants << std::endl;
            
            /*
            int objs = djn_decode->DecodeNextRaw(variant);
            
            uint32_t of = 0;
            for (int j = 0; j < variant->d->n_ewah; ++j) {
                std::cerr << variant->d->ewah[j]->clean << ":" << variant->d->ewah[j]->ref;
                for (int k = 0; k < variant->d->ewah[j]->dirty; ++k, ++of) {
                    std::cerr << "," << std::bitset<32>(*variant->d->dirty[of]);
                }
                std::cerr << std::endl;
            }
            */
            
            // int objs = djn_decode->DecodeNext(variant);
            int objs = djn_decode->DecodeNextRaw(variant);
            assert(objs > 0);

            std::cerr << variant->ploidy << "," << variant->n_allele << "," << variant->data_len << "," << variant->errcode << "," << variant->unpacked << ": data=";
            for (int i = 0; i < variant->d->n_ewah; ++i) {
                std::cerr << variant->d->ewah[i]->clean << ":" << variant->d->ewah[i]->ref << "+" << variant->d->ewah[i]->dirty << ", ";
            }
            std::cerr << std::endl;

            // Write Vcf-encoded data to a local buffer and then write to standard out.
            // for (int j = 0; j < variant->data_len; j += variant->ploidy) {
            //     // Todo: phasing
            //     vcf_out_buffer[len_vcf++] = (char)variant->data[j+0] + '0';
            //     vcf_out_buffer[len_vcf++] = '|';
            //     vcf_out_buffer[len_vcf++] = (char)variant->data[j+1] + '0';
            //     vcf_out_buffer[len_vcf++] = '\t';
            // }
            // vcf_out_buffer[len_vcf++] = '\n';
            // std::cout.write((char*)vcf_out_buffer, len_vcf);
            // len_vcf = 0;
        }
        // Timings per block
        clockdef t2 = std::chrono::high_resolution_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
        std::cerr << "[Decode] Decoded " << djn_decode->n_variants << " records in " 
            << time_span.count() << "us (" << (double)time_span.count()/djn_decode->n_variants << "us/record)" << std::endl;
    }
    // Timings total
    clockdef t2_decode = std::chrono::high_resolution_clock::now();
    auto time_span_decode = std::chrono::duration_cast<std::chrono::milliseconds>(t2_decode - t1_decode);
    std::cerr << "[Decode] Decoded " << n_lines << " records in " << time_span_decode.count() 
        << "ms (" << (double)time_span_decode.count()/n_lines << "ms/record)" << std::endl;

    delete[] decode_buf;
    delete[] vcf_out_buffer;
    delete variant;
    delete djn_decode;

    return n_lines;
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
    printf("   -P BOOL   do NOT permute data with PBWT\n\n");
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
        // return(DecompressExample(input, type));
        return EXIT_FAILURE;
    }

    return EXIT_FAILURE;
}
