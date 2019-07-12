/*
* Copyright (c) 2019 Marcus D. R. Klarqvist
* Author(s): Marcus D. R. Klarqvist
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*   http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing,
* software distributed under the License is distributed on an
* "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
* KIND, either express or implied.  See the License for the
* specific language governing permissions and limitations
* under the License.
*/
#ifndef DJINN_EXAMPLE_ENCODE_H_
#define DJINN_EXAMPLE_ENCODE_H_

#include <chrono>//time
#include <random>//dist*
#include <fstream> // Support for read/write.
#include <djinn.h> // Djinn data models.

/**
 * In this example we will generate random data for import into Djinn. As an example,
 * we will generate random genotypes at 500,000 sites for 100,000 alleles. Using a
 * diploid model this will result in 25,000 diploid genotypes.
 * 
 * @param input_file   Input file string: file path or "-" to read from stdin
 * @param output_file  Output file string: file path or "-" to write to stdout
 * @param type         1: ctx model, 2; LZ4-EWAH, 4: ZSTD-EWAH
 * @param permute      Use PBWT preprocessor
 * @param reset_models Reset models for each block (random access)
 * @return int         Returns the number of imported variants when successful or a negative value otherwise.
 */
int ImportRandom(std::string input_file,   // input file: "-" for stdin
                 std::string output_file,  // output file: "-" for stdout
                 const uint32_t type,      // 1: ctx model, 2; LZ4-EWAH, 4: ZSTD-EWAH
                 const bool permute = true,// PBWT preprocessor
                 const bool reset_models = true) // Reset models for each block (random access)
{
    // setup random
    int64_t  n_fake_samples = 100000;
    uint32_t n_fake_sites   = 500000;
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<uint32_t> dis(0, n_fake_samples-1);
    std::uniform_int_distribution<uint32_t> freq_dis(0, 1000); // ~U(0,1000)
    // Allocate data for random byte array.
    uint8_t* rand_vec   = new uint8_t[n_fake_samples];

    // Setup
    uint64_t n_lines   = 0;    // Keep track of how many variants we've imported
    uint32_t nv_blocks = 8192; // Number of desired variants per data block.
    uint32_t n_blocks  = 0;    // Keep track of how many data blocks we've processed.

    djinn::djinn_model* djn_ctx = nullptr;
    if ((type >> 0) & 1)      djn_ctx = new djinn::djinn_ctx_model();
    else if ((type >> 1) & 1) djn_ctx = new djinn::djinn_ewah_model(djinn::CompressionStrategy::LZ4,  1);
    else if ((type >> 2) & 1) djn_ctx = new djinn::djinn_ewah_model(djinn::CompressionStrategy::ZSTD, 1);
    djn_ctx->StartEncoding(permute, reset_models);
    
    // Open file stream (or file handle) depending on the passed argument.
    bool own_stream = false;
    std::ostream* out_stream = nullptr;
    if (output_file == "-") out_stream = &std::cout; // standard out (pipe)
    else { // file stream (to disk)
        out_stream = new std::ofstream(output_file, std::ios::out | std::ios::binary);
        if (out_stream->good() == false) {
            std::cerr << "Could not open output handle \"" << output_file << "\"!" << std::endl;
            return -3;
        }
        own_stream = true;
    }

    // Cumulators to print our progress.
    uint64_t data_in = 0, data_in_vcf = 0;
    uint64_t ctx_out = 0;
    
    for(int i = 0; i < n_fake_sites; ++i) {
        // Every "nv_block" of variants we finish encoding and serialize
        // the encoded data to the output stream.
        if (n_lines % nv_blocks == 0 && n_lines != 0) {
            djn_ctx->FinishEncoding();
            int decode_ret2 = djn_ctx->Serialize(*out_stream);
            djn_ctx->StartEncoding(permute, reset_models);
            ++n_blocks;
            ctx_out += decode_ret2;

            std::cerr << "[PROGRESS] In uBCF: " << data_in << "->" << ctx_out 
                << " (" << (double)data_in/ctx_out << "-fold) In VCF: " << data_in_vcf << "->" << ctx_out 
                << " (" << (double)data_in_vcf/ctx_out << "-fold)" << std::endl;
        }

        // Generate random data such that each byte is a dictionary-encoded 
        // ref/alt allele as in Bcf. For example, 0 is the reference allele,
        // and 1 is the first alternative allele, and 2 is the second
        // alternative allele, and so on.
        memset(rand_vec, 0, n_fake_samples); // Initiate array to all reference.
        const int n_alts = freq_dis(gen); // Draw a random number of alternative alleles.
        for (int j = 0; j < n_alts; ++j) // Place #n_alts of alternative alleles.
            rand_vec[dis(gen)] = 1; // Place an alt allele at a random position.

        // Encode random data. Encode function is parameterized as
        // (data array, length of data array, ploidy, number of alternative alleles)
        int ret = djn_ctx->Encode(rand_vec, n_fake_samples, 2, 2);
        ++n_lines;
        data_in     += n_fake_samples; // Keep track of total input data.
        data_in_vcf += 2*n_fake_samples - 1; // Track equivalent Vcf size.
    }
    // Cleanup.
    delete[] rand_vec;
    // Return the number of imported variants.
    return n_lines;
}

#endif