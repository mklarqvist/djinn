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

#include "frequency_model.h"

#include <iostream>//temp

namespace djinn {

// Map missing to 2, 1->0, 2->1, and EOV -> 3.
const uint8_t TWK_BCF_GT_UNPACK[65] = {2,0,1,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3};
const uint8_t TWK_BCF_GT_PACK[3]   = {1, 2, 0};
// const uint8_t TWK_BCF_GT_UNPACK_GENERAL[16] = {15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

// MISSING -> 14, EOV -> 15, other values as normal
const uint8_t TWK_BCF_GT_UNPACK_GENERAL[131] = 
{14,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,15,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129};

const uint8_t TWK_BCF_GT_UNPACK_GENERAL_REV[131] = 
{0,1,2,3,4,5,6,7,8,9,10,11,12,13,0,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,15,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129};


#define BCF_UNPACK_GENOTYPE(A) TWK_BCF_GT_UNPACK[(A) >> 1]
#define BCF_UNPACK_GENOTYPE_GENERAL(A) TWK_BCF_GT_UNPACK_GENERAL[(A) >> 1]

/*======   Canonical representation   ======*/

class PBWT {
public:
    PBWT();
    PBWT(int64_t n_samples, int n_symbols);
    ~PBWT();

    void Initiate(int64_t n_s, int n_sym);
    void reset();
    int Update(const int* arr);
    int Update(const uint8_t* arr, uint32_t stride = 1);
    int UpdateGeneral(const uint8_t* arr, uint32_t stride = 1);
    // Todo:
    // Update the PBWT given the following positions.
    // Only useful for biallelic PBWTs.
    int Update(const uint32_t* alt_pos, const uint32_t len);

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

/*======   Higher order   ======*/

#define MODEL_SIZE 65536

class GeneralPBWTModel {
public:
    GeneralPBWTModel() noexcept;
    GeneralPBWTModel(int64_t n_samples, int n_symbols);
    ~GeneralPBWTModel();

    void Construct(int64_t n_samples, int n_symbols);
    void ResetModels();
    void ResetPBWT();
    void ResetContext();
    void Reset();
    void ResetExceptPBWT();
    int FinishEncoding();
    int FinishDecoding();
    void StartEncoding();
    void StartDecoding(uint8_t* data);
    void EncodeSymbol(const uint16_t symbol);
    uint16_t DecodeSymbol();

public:
    int max_model_symbols;
    int model_context_shift;
    uint32_t model_context;
    std::shared_ptr<PBWT> pbwt;
    std::shared_ptr<RangeCoder> range_coder;
    std::vector < std::shared_ptr<FrequencyModel> > models;
    size_t n_buffer;
    uint8_t* buffer; // fix
    size_t n_additions;
};

}

#endif