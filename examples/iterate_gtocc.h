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
#ifndef DJINN_EXAMPLE_ITERATE_OCC_H_
#define DJINN_EXAMPLE_ITERATE_OCC_H_

#include <fstream> // Support for read/write.
#include <djinn.h> // Djinn data models.
#include <random>


/**
 * In this example we will read data using the provided support class VcfReader
 * that interacts directly with Htslib data structures to consume either 
 * Vcf/Bcf/Vcf.gz/BGZF-based files either from standard in (pipe) or from a 
 * file handle (from disk).
 * 
 * 
 * @param input_file   Input file string: file path or "-" to read from stdin
 * @param output_file  Output file string: file path or "-" to write to stdout
 * @param type         1: ctx model, 2; LZ4-EWAH, 4: ZSTD-EWAH
 * @param permute      Use PBWT preprocessor
 * @param reset_models Reset models for each block (random access)
 * @return int         Returns non-negative value when successful or a negative value otherwise.
 */
int IterateOcc(std::string input_file, int model, int n_sampling) {
    bool own_stream = false;
    uint64_t filesize = 0;
    std::istream* in_stream = nullptr;
    if (input_file == "-") in_stream = &std::cin;
    else {
        in_stream = new std::ifstream(input_file, std::ios::in | std::ios::binary | std::ios::ate);
        if (in_stream->good() == false) {
            std::cerr << "could not open infile handle" << std::endl;
            return -2;
        }
        filesize = in_stream->tellg();
        in_stream->seekg(0);
    }

    djinn::djinn_model* djn_decode = nullptr;
    if (model == 1) djn_decode = new djinn::djinn_ctx_model();
    else if(model == 2 || model == 4) djn_decode = new djinn::djinn_ewah_model();
    else {
        std::cerr << "unknown model: " << model << std::endl;
        return -3;
    }

    uint32_t n_lines   = 0;
    djinn::djinn_variant_t* variant = nullptr;
    djinn::djinn_variant_t* variant2 = nullptr;
    

    // Construct Occ table
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<uint32_t> dis(0, 5096-1);

    djinn::djinn_occ occ(5096);
    std::vector<uint32_t> target_ids;
    for (int i = 0; i < n_sampling; ++i) {
        target_ids.push_back(dis(gen));
    }
    assert(occ.AddGroup(target_ids));
    occ.BuildTable();

    while (true) {
        int decode_ctx_ret = djn_decode->Deserialize(*in_stream);
        if (decode_ctx_ret <= 0) break; // exit condition

        djn_decode->StartDecoding();
        for (int i = 0; i < djn_decode->n_variants; ++i, ++n_lines) {
            int objs = djn_decode->DecodeNextRaw(variant);
            assert(objs > 0);
            int ret = variant->Slice(occ, 0, variant2, false);
        }
    }
    delete variant;
    delete djn_decode;

    return n_lines;
}

#endif