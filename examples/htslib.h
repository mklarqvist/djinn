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
#ifndef DJINN_EXAMPLE_HTSLIB_H_
#define DJINN_EXAMPLE_HTSLIB_H_

#include <fstream> // Support for read/write.
#include <djinn.h> // Djinn data models.
#include <vcf_reader.h> // VcfReader support class for reading Htslib-based files.
                        // Compiling with this header requires the htslib library.

/**
 * In this example we will read data using the provided support class VcfReader
 * that interacts directly with Htslib data structures to consume either 
 * Vcf/Bcf/Vcf.gz/BGZF files either from stdin (pipe) or from a 
 * file handle (disk).
 * 
 * 
 * @param input_file   Input file string: file path or "-" to read from stdin
 * @param output_file  Output file string: file path or "-" to write to stdout
 * @param type         1: ctx model, 2; LZ4-EWAH, 4: ZSTD-EWAH
 * @param permute      Use PBWT preprocessor
 * @param reset_models Reset models for each block (random access)
 * @return int         Returns the number of imported variants when successful or a negative value otherwise.
 */
int ImportHtslib(std::string input_file,   // input file: "-" for stdin
                 std::string output_file,  // output file: "-" for stdout
                 const uint32_t type,      // 1: ctx model, 2; LZ4-EWAH, 4: ZSTD-EWAH
                 const bool permute = true,// PBWT preprocessor
                 const bool reset_models = true) // Reset models for each block (random access)
{
    // VcfReader use a singleton pattern: call the djinn::VcfReader::FromFile
    // function to get the instance.
    std::unique_ptr<djinn::VcfReader> reader = djinn::VcfReader::FromFile(input_file);
    
    // If the file or stream could not be opened we exit here.
    if (reader.get() == nullptr) {
        std::cerr << "Could not open input handle \"" << input_file << "\"!" << std::endl;
        return -1;
    }

    // Print out number of samples listed in the bcf header.
    std::cerr << "Samples in VCF file: " << reader->n_samples_ << std::endl;
    
    // Setup
    uint64_t n_lines   = 0;    // Keep track of how many variants we've imported
    uint32_t nv_blocks = 8192; // Number of desired variants per data block.
    uint32_t n_blocks  = 0;    // Keep track of how many data blocks we've processed.

    djinn::djinn_model* djn_ctx = nullptr;
    if ((type >> 0) & 1)      djn_ctx = new djinn::djinn_ctx_model();
    else if ((type >> 1) & 1) djn_ctx = new djinn::djinn_ewah_model(djinn::CompressionStrategy::LZ4,  9);
    else if ((type >> 2) & 1) djn_ctx = new djinn::djinn_ewah_model(djinn::CompressionStrategy::ZSTD, 21);
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
    uint64_t data_in = 0, data_in_vcf = 0, model_out = 0;

    while (reader->Next()) {
        // When we have decoded nv_blocks of variants we will stop encoding data
        // by calling FinishEncoding and then serialize the final encoded object
        // to the output stream.
        if (n_lines % nv_blocks == 0 && n_lines != 0) {
            // Calling FinisheEncoding is REQUIRED before either Serializing and
            // writing or decompressing.
            djn_ctx->FinishEncoding();
            int serial_size = djn_ctx->Serialize(*out_stream);
            ++n_blocks;
            model_out += serial_size;

            std::cerr << "[PROGRESS] In uBCF: " << data_in << "->" << model_out 
                << " (" << (double)data_in/model_out << "-fold) In VCF: " << data_in_vcf << "->" << model_out 
                << " (" << (double)data_in_vcf/model_out << "-fold)" << std::endl;

            djn_ctx->StartEncoding(permute, reset_models);
        }

        // Error handling: if either bcf1_t or bcf_hdr_t pointers are NULL then
        // a problem has occured.
        if (reader->bcf1_   == NULL) return -2;
        if (reader->header_ == NULL) return -3;

        // Retrieve pointer to FORMAT field that holds GT data.
        const bcf_fmt_t* fmt = bcf_get_fmt(reader->header_, reader->bcf1_, "GT");
        if (fmt == NULL) continue;
        
        // Encode from htslib Bcf encoding by passing the arguments:
        // p: pointer to genotype data array
        // p_len: length of data
        // n: stride size (number of bytes per individual = base ploidy)
        // n_allele: number of alleles
        int ret = djn_ctx->EncodeBcf(fmt->p, fmt->p_len, fmt->n, reader->bcf1_->n_allele);
        assert(ret>0);

        // Update input bytes for uBcf and Vcf
        data_in     += fmt->p_len; // uBcf
        data_in_vcf += 2*fmt->p_len - 1; // Vcf: this is true only for diploid data with #alleles < 10
        ++n_lines; // Number of variants processed
    }

    // Compress final data.
    djn_ctx->FinishEncoding();
    int serial_size = djn_ctx->Serialize(*out_stream);
    assert(serial_size > 0);
    ++n_blocks;
    model_out += serial_size;

    std::cerr << "[PROGRESS] In uBCF: " << data_in << "->" << model_out 
        << " (" << (double)data_in/model_out << "-fold) In VCF: " << data_in_vcf << "->" << model_out 
        << " (" << (double)data_in_vcf/model_out << "-fold)" << std::endl;

    // Close handle and clean up.
    out_stream->flush();
    
    if (own_stream) {
        ((std::ofstream*)out_stream)->close();
        delete out_stream;
    }

    return n_lines;
}

#endif