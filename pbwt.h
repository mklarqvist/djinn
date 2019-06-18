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
#ifndef PBWT_H_
#define PBWT_H_

#include <cstdint>//types
#include <cstring>//memcpy
#include <cassert>//assert
#include <cmath>//ceil
#include <memory>//pointers
#include <vector>//vector

#include "djinn.h"
#include "frequency_model.h"

#include <iostream>//temp

namespace djinn {

/*======   PBWT   ======*/

/*
 *--------------------------------------------------------------------------
 * Basic implementation of the positional Burrow-Wheeler transform (PBWT)
 * as described in:
 * 
 * Richard Durbin, Efficient haplotype matching and storage using the positional 
 * Burrows–Wheeler transform (PBWT), Bioinformatics, Volume 30, Issue 9, 
 * 1 May 2014, Pages 1266–1272, https://doi.org/10.1093/bioinformatics/btu014
 *
 * Permutes an input vector of alleles such that haplotypes up to the current
 * position is sorted according to their reverse prefix. It should be noted
 * that the output is the input vector of alleles reverse-prefix sorted up
 * to the previous position.
 * 
 * This implementation handles arbitrarily large alphabets.
 *
 *--------------------------------------------------------------------------
 * Usage instructions
 * 
 * Basic usage for diploid biallelic:
 * PBWT pbwt;
 * pbwt.Initiate(5008, 2); // Initiate class with 2504 samples and alphabet size of 2
 * pbwt.UpdateBcf(data); // data is an array of Bcf-encoded values
 * 
 * pbwt.ppa[0-5007] now stores the permuted order.
 *  
 * Advanced usage for diploid biallelic data storing haplotypes separately:
 * PBWT pbwt1, pbwt2;
 * pbwt1.Initiate(2504, 2); // Initiate class with 2504 samples and alphabet size of 2
 * pbwt2.Initiate(2504, 2);
 * pbwt1.Update(&data[0], 2); // Update PBWT with a stride size of 2 starting at offset 0 
 * pbwt2.Update(&data[1], 2); // Starting at offset 1
 * 
 * pbwt1.ppa[0-2503] now stores the permuted order for haplotype 1.
 * pbwt2.ppa[0-2503] now stores the permuted order for haplotype 2.
 *--------------------------------------------------------------------------
 */
class PBWT {
public:
    PBWT();
    PBWT(int64_t n_samples, int n_symbols);
    ~PBWT();

    // Disallow all forms of copying and moving.
    PBWT(const PBWT& other) = delete;
    PBWT(PBWT&& other) noexcept = delete;
    PBWT& operator=(const PBWT& other) = delete;
    PBWT& operator=(PBWT&& other) noexcept = delete;

    void Initiate(int64_t n_s, int n_sym);
    void Reset();
    
    int Update(const uint8_t* arr, uint32_t stride = 1);
    int UpdateBcf(const uint8_t* arr, uint32_t stride = 1);
    int UpdateBcfGeneral(const uint8_t* arr, uint32_t stride = 1);
    
    // Encode WAH in 63-bits and return data.
    int UpdateBcfWah(const uint8_t* arr, uint8_t* out, uint32_t stride = 1);
    // Encode WAH in 64-bits with archetype data.
    int UpdateBcfWah2(const uint8_t* arr, uint8_t* data, uint32_t* archetype, uint32_t stride = 1);

    int ReverseUpdate(const uint8_t* arr);
    int ReverseUpdateBitmap(const uint8_t* arr);

    std::string ToPrettyString() const;

public:
    int        n_symbols; // universe of symbols (number of unique symbols)
    int64_t    n_samples; // number of samples (free interpretation)
    uint64_t   n_steps; // number of updates made (debugging)
    uint8_t*   prev; // previous output array
    uint32_t*  ppa; // current PPA
    uint32_t*  n_queue; // number of elements in each positional queue
    uint32_t** queue; // the positional queues themselves
    uint64_t*  prev_bitmap; // bitmap version
};

}

#endif