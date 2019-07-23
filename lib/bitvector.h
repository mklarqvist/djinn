
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
#ifndef DJINN_TEMP_H_
#define DJINN_TEMP_H_

#include <cmath> //ceil
#include <cassert> //assert
#include <cstring> //memset

namespace djinn {

struct djn_bv_t {
public:
    djn_bv_t() : bv_len(0), s_len(0), p_cap(0), p_free(false), bv(nullptr) {}
    ~djn_bv_t() {
        if (p_free) {
            aligned_free(bv);
            // delete bv;
        }
    }

    int Encode(const uint8_t* data, const uint32_t data_len) {
        if (data == nullptr) return 0;
        if (data_len == 0) return 0;
        
        s_len = data_len;
        bv_len = std::ceil((double)data_len/64);
        
        if (bv == nullptr) {
            // std::cerr << "allocating nullptr" << std::endl;
            p_cap = std::ceil((double)data_len/64)+32;
            // bv = new uint64_t[p_cap];
            bv = (uint64_t*)aligned_malloc(64, p_cap);
            p_free = true;
        }

        if (bv_len > p_cap*sizeof(uint64_t)) {
            // std::cerr << "resizing: " << bv_len << ">" << (p_cap*sizeof(uint64_t)) << std::endl;
            aligned_free(bv);
            p_cap = std::ceil((double)data_len/64)+32;
            // bv = new uint64_t[p_cap];
            bv = (uint64_t*)aligned_malloc(64, p_cap);
            p_free = true;
        }

        memset(bv, 0, bv_len*sizeof(uint64_t));

        for (int i = 0; i < data_len; ++i) {
            assert(data[i] < 2);
            bv[i / 64] |= (uint64_t)data[i] << (i % 64);
        }

        return 1;
    }

    int Add(const djn_bv_t& other) {
        if (other.bv == nullptr) return 0;
        if (other.s_len == 0) return 0;

        if (bv == nullptr) {
            // std::cerr << "allocating nullptr" << std::endl;
            p_cap = std::ceil((double)other.s_len/64)+32;
            // bv = new uint64_t[p_cap];
            bv = (uint64_t*)aligned_malloc(64, p_cap);
            p_free = true;
        }

        if (other.bv_len + bv_len > p_cap*sizeof(uint64_t)) {
            // std::cerr << "resizing: " << bv_len << ">" << (p_cap*sizeof(uint64_t)) << std::endl;
            aligned_free(bv);
            p_cap = bv_len + bv_len + 32;
            // bv = new uint64_t[p_cap];bv_len - (s_len / 64)
            bv = (uint64_t*)aligned_malloc(64, p_cap);
            p_free = true;
        }

        // Todo


        return 1;
    }

    uint8_t operator[](const uint32_t p) const { return (bv[p / 64] & (1ULL << (p % 64))); }

public:
    uint32_t bv_len, s_len; // data length in 64-bits, actual length
    uint32_t p_cap:31, p_free:1; // capacity (memory allocated), flag for data ownership
    uint64_t* bv;
};

}

#endif