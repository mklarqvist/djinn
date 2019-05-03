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

const uint8_t TWK_BCF_GT_UNPACK[3] = {2, 0, 1};
const uint8_t TWK_BCF_GT_UNPACK_GENERAL[16] = {15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
#define BCF_UNPACK_GENOTYPE(A) TWK_BCF_GT_UNPACK[(A >> 1)]
#define BCF_UNPACK_GENOTYPE_GENERAL(A) TWK_BCF_GT_UNPACK_GENERAL[(A >> 1)]

/*======   Canonical representation   ======*/

class PBWT {
public:
    PBWT(int64_t n_samples, int n_symbols);
    ~PBWT();

    void Initiate(int64_t n_samples);
    void reset();
    int Update(const int* arr);
    int Update(const uint8_t* arr, uint32_t stride = 1);
    int UpdateGeneral(const uint8_t* arr, uint32_t stride = 1);

    std::string ToPrettyString() const;

public:
    const int  n_symbols; // universe of symbols (number of unique symbols)
    int64_t    n_samples; // number of samples (free interpretation)
    uint64_t   n_steps; // number of updates made (debugging)
    uint8_t*   prev; // previous output array
    uint32_t*  ppa; // current PPA
    uint32_t*  n_queue; // number of elements in each positional queue
    uint32_t** queue; // the positional queues themselves
};

/*======   Higher order   ======*/

#define MODEL_SIZE 256

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
    int FinishEncoding();
    void StartEncoding();
    void EncodeSymbol(const uint16_t symbol);

public:
    int max_model_symbols;
    int model_context_shift;
    uint32_t model_context;
    std::shared_ptr<PBWT> pbwt;
    std::shared_ptr<pil::RangeCoder> range_coder;
    std::vector < std::shared_ptr<pil::FrequencyModel> > models;
    size_t n_buffer;
    uint8_t* buffer; // fix
};

#endif