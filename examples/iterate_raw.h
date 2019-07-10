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
#ifndef DJINN_EXAMPLE_ITERATE_RAW_H_
#define DJINN_EXAMPLE_ITERATE_RAW_H_

#include <fstream> // Support for read/write.
#include <djinn.h> // Djinn data models.

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
int IterateRaw(std::string input_file) {
    
    return 1;
}

#endif