#include "pbwt.h"

// temp
#include <iostream>
#include <bitset>
#include <cmath>//ceil

namespace djinn {

PBWT::PBWT() :
    n_symbols(0),
    n_samples(0),
    n_steps(0),
    prev(nullptr),
    ppa(nullptr),
    n_queue(nullptr),
    queue(nullptr),
    prev_bitmap(nullptr)
{

}

PBWT::PBWT(int64_t n_samples, int n_symbols) :
    n_symbols(n_symbols),
    n_samples(n_samples),
    n_steps(0),
    prev(new uint8_t[n_samples]),
    ppa(new uint32_t[n_samples]),
    n_queue(new uint32_t[n_symbols]),
    queue(new uint32_t*[n_symbols]),
    prev_bitmap(nullptr)
{
    assert(n_symbols > 1);

    for (int i = 0; i < n_symbols; ++i)
        queue[i] = new uint32_t[n_samples];

    memset(n_queue, 0, sizeof(uint32_t)*n_symbols);
    Reset();
}

PBWT::~PBWT() {
    delete[] prev;
    delete[] ppa;

    for (int i = 0; i < n_symbols; ++i)
        delete[] queue[i];

    delete[] queue;
    delete[] n_queue;
    delete[] prev_bitmap;
}

void PBWT::Initiate(int64_t n_s, int n_sym) {
    assert(n_sym > 1);

    n_samples = n_s;
    n_steps = 0;
    delete[] prev; delete[] ppa;
    delete[] n_queue;
    if (queue != nullptr) {
        for (int i = 0; i < n_symbols; ++i)
            delete[] queue[i];
    }
    n_symbols = n_sym;

    prev = new uint8_t[n_samples];
    ppa = new uint32_t[n_samples];
    n_queue = new uint32_t[n_symbols];
    queue = new uint32_t*[n_symbols];

    for (int i = 0; i < n_symbols; ++i)
        queue[i] = new uint32_t[n_samples];

    memset(n_queue, 0, sizeof(uint32_t)*n_symbols);
    Reset();
}

void PBWT::Reset() {
    if (ppa != nullptr) {
        for (int i = 0; i < n_samples; ++i)
            ppa[i] = i;
    }
    if (prev != nullptr) {
        memset(prev, 0, sizeof(uint8_t)*n_samples);
    }
    memset(n_queue, 0, sizeof(uint32_t)*n_symbols);
}

int PBWT::UpdateBcf(const uint8_t* arr, uint32_t stride) {
    // Reset queues.
    memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

    for (int i = 0; i < n_samples; ++i) {
        const uint8_t& gt = BCF_UNPACK_GENOTYPE(arr[ppa[i] * stride]);
        // if (gt >= n_symbols) {
        //     std::cerr << "error=" << (int)arr[ppa[i]*stride] << "->" << (int)gt << ">=" << n_symbols << std::endl;
        // }
        assert(gt < n_symbols);
        queue[gt][n_queue[gt]++] = ppa[i];
        prev[i] = gt;
        // std::cerr << " " << (int)gt;
    }
    // std::cerr << std::endl << std::endl;

    uint32_t of = 0;
    for (int j = 0; j < n_symbols; ++j) {
        for (uint32_t i = 0; i < n_queue[j]; ++i, ++of)
            ppa[of] = queue[j][i];
    }
    //for (int i = 0; i < n_q2; ++i, ++of) ppa[of] = q2[i];
    //std::cerr << "of=" << of << "/" << n_samples << std::endl;
    assert(of == n_samples);
    // Debug: data is sorted at this point.
    // std::cerr << ToPrettyString() << std::endl;
    // std::cerr << "Alts=" << n_queue[1] << std::endl;
    // for (int i = 0; i < n_samples; ++i) {
    //   std::cerr << " " << (int)BCF_UNPACK_GENOTYPE(arr[ppa[i]*stride]);
    // }
    // std::cerr << std::endl;
    ++n_steps;

    return(1);
}

int PBWT::UpdateBcfGeneral(const uint8_t* arr, uint32_t stride) {
    // Reset queues.
    memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

    for (int i = 0; i < n_samples; ++i) {
        const uint8_t& gt = BCF_UNPACK_GENOTYPE_GENERAL(arr[ppa[i] * stride]);
        if (gt >= n_symbols) {
            std::cerr << "error=" << (int)gt << "/" << n_symbols << std::endl;
        }
        assert(gt < n_symbols);
        for (int j = 0; j < n_symbols; ++j) {
            if (gt == j)
                queue[j][n_queue[j]++] = ppa[i];
        }
        prev[i] = gt;
    }

    uint32_t of = 0;
    for (int j = 0; j < n_symbols; ++j) {
        for (uint32_t i = 0; i < n_queue[j]; ++i, ++of)
            ppa[of] = queue[j][i];
    }
    assert(of == n_samples);
    ++n_steps;

    return(1);
}

int PBWT::Update(const uint8_t* arr, uint32_t stride) {
    // Reset queues.
    memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

    for (int i = 0; i < n_samples; ++i) {
        const uint8_t& gt = arr[ppa[i] * stride];
        if (gt >= n_symbols) {
            std::cerr << "error=" << (int)gt << "/" << n_symbols << std::endl;
        }
        assert(gt < n_symbols);
        for (int j = 0; j < n_symbols; ++j) {
            if (gt == j)
                queue[j][n_queue[j]++] = ppa[i];
        }
        prev[i] = gt;
    }

    uint32_t of = 0;
    for (int j = 0; j < n_symbols; ++j) {
        for (uint32_t i = 0; i < n_queue[j]; ++i, ++of)
            ppa[of] = queue[j][i];
    }
    assert(of == n_samples);
    ++n_steps;

    return(1);
}

std::string PBWT::ToPrettyString() const {
    std::string ret = "n=" + std::to_string(n_samples) + " {";
    ret += std::to_string(ppa[0]);
    for (int i = 1; i < n_samples; ++i) {
        ret += "," + std::to_string(ppa[i]);
    }
    ret += "}";
    return(ret);
}

int PBWT::ReverseUpdate(const uint8_t* arr) {
    memset(prev, 0, n_samples); // O(n)
    memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

    // Restore + update PPA
    // uint32_t n_skips = 0;
    for (int i = 0; i < n_samples; ++i) { // Worst case O(n), average case O(n) with a smallish constant.
        queue[arr[i]][n_queue[arr[i]]++] = ppa[i];
        if (arr[i] == 0) {
            // ++n_skips;
            continue;
        }
        prev[ppa[i]] = arr[i]; // Unpermute data.
    }

    // Merge PPA queues.
    uint32_t of = 0;
    for (int j = 0; j < n_symbols; ++j) { // O(n)
        memcpy(&ppa[of], queue[j], sizeof(uint32_t)*n_queue[j]);
        // for (int i = 0; i < n_queue[j]; ++i, ++of) {
        //     ppa[of] = queue[j][i];
        // }
        of += n_queue[j];
    }
    assert(of == n_samples);
    ++n_steps;

    // std::cerr << "nskips=" << n_skips << std::endl;

    return 1;
}

int PBWT::ReverseUpdateEWAH(const uint8_t* arr, const uint32_t len) {
    assert(n_samples > 0);
    assert(n_symbols > 0);

    memset(prev, 0, n_samples); // O(n)
    memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

    uint32_t local_offset = 0;
    uint64_t n_s_obs = 0;
    while (local_offset < len) {
        djinn_ewah_t* ewah = (djinn_ewah_t*)&arr[local_offset];
        local_offset += sizeof(djinn_ewah_t);

        // std::cerr << "ewah=" << ewah->ref << "," << ewah->clean << "," << ewah->dirty << std::endl;
        
        // Clean
        uint64_t to = n_s_obs + ewah->clean * 32 > n_samples ? n_samples : n_s_obs + ewah->clean * 32;
        // std::cerr << "clean=" << n_s_obs << "->" << to << std::endl; 
        for (int i = n_s_obs; i < to; ++i) {
            queue[ewah->ref & 1][n_queue[ewah->ref & 1]++] = ppa[i];
        }
        n_s_obs = to;

        // Loop over dirty bitmaps.
        // std::cerr << "dirty=" << ewah->dirty << std::endl;
        for (int i = 0; i < ewah->dirty; ++i) {
            to = n_s_obs + 32 > n_samples ? n_samples - n_s_obs : 32;
            // std::cerr << "dirty steps=" << to << " -> " << n_s_obs << "-" << n_s_obs+to << "/" << n_samples << std::endl;
            assert(n_s_obs < n_samples);
            
            uint32_t dirty = *((uint32_t*)(&arr[local_offset]));
            for (int j = 0; j < to; ++j) {
                queue[dirty & 1][n_queue[dirty & 1]++] = ppa[n_s_obs];
                prev[ppa[n_s_obs]] = (dirty & 1); // update prev when non-zero
                dirty >>= 1;
                ++n_s_obs;
            }
            local_offset += sizeof(uint32_t);
            // n_s_obs += to;
        }
        // local_offset += ewah->dirty * sizeof(uint32_t);
        // std::cerr << "pbwt offset=" << local_offset << "/" << len << " with=" << n_s_obs << "/" << n_samples << std::endl;
        assert(local_offset <= len);
    }
    // std::cerr << "obs" << n_s_obs << "/" << n_samples << std::endl;
    assert(n_s_obs == n_samples);

    // Restore + update PPA
    // uint32_t n_skips = 0;
    // for (int i = 0; i < n_samples; ++i) { // Worst case O(n), average case O(n) with a smallish constant.
    //     queue[arr[i]][n_queue[arr[i]]++] = ppa[i];
    //     if (arr[i] == 0) {
    //         // ++n_skips;
    //         continue;
    //     }
    //     prev[ppa[i]] = arr[i]; // Unpermute data.
    // }

    // Merge PPA queues.
    uint32_t of = 0;
    for (int j = 0; j < n_symbols; ++j) { // O(n)
        // std::cerr << "queue=" << j << ": " << n_queue[j] << std::endl;
        memcpy(&ppa[of], queue[j], sizeof(uint32_t)*n_queue[j]);
        // for (int i = 0; i < n_queue[j]; ++i, ++of) {
        //     ppa[of] = queue[j][i];
        // }
        of += n_queue[j];
    }
    // std::cerr << "of=" << of << "/" << n_samples << std::endl;
    assert(of == n_samples);
    ++n_steps;

    // std::cerr << "nskips=" << n_skips << std::endl;

    return 1;
}

int PBWT::ReverseUpdateEWAH(const uint8_t* arr, const uint32_t len, uint8_t* ret) {
    assert(n_samples > 0);
    assert(n_symbols > 0);

    memset(ret, 0, n_samples); // O(n)
    memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

    uint32_t local_offset = 0;
    uint64_t n_s_obs = 0;
    while (local_offset < len) {
        djinn_ewah_t* ewah = (djinn_ewah_t*)&arr[local_offset];
        local_offset += sizeof(djinn_ewah_t);

        // std::cerr << "ewah=" << ewah->ref << "," << ewah->clean << "," << ewah->dirty << std::endl;
        
        // Clean words.
        uint64_t to = n_s_obs + ewah->clean * 32 > n_samples ? n_samples : n_s_obs + ewah->clean * 32;
        // std::cerr << "clean=" << n_s_obs << "->" << to << std::endl; 
        for (int i = n_s_obs; i < to; ++i) {
            queue[ewah->ref & 1][n_queue[ewah->ref & 1]++] = ppa[i];
            ret[ppa[i]] = (ewah->ref & 1); // update prev when non-zero
        }
        n_s_obs = to;

        // Loop over dirty bitmaps.
        // std::cerr << "dirty=" << ewah->dirty << std::endl;
        for (int i = 0; i < ewah->dirty; ++i) {
            to = n_s_obs + 32 > n_samples ? n_samples : n_s_obs + 32;
            // std::cerr << "dirty steps=" << n_s_obs << "->" << to << ": " << to-n_s_obs << std::endl;
            assert(n_s_obs < n_samples);
            assert(to <= n_samples);
            
            uint32_t dirty = *((uint32_t*)(&arr[local_offset])); // copy
            for (int j = n_s_obs; j < to; ++j) {
                queue[dirty & 1][n_queue[dirty & 1]++] = ppa[j];
                ret[ppa[j]] = (dirty & 1); // update prev when non-zero
                dirty >>= 1;
                // ++n_s_obs;
            }
            local_offset += sizeof(uint32_t);
            n_s_obs = to;
        }
        // local_offset += ewah->dirty * sizeof(uint32_t);
        // std::cerr << "pbwt offset=" << local_offset << "/" << len << " with=" << n_s_obs << "/" << n_samples << std::endl;
        assert(local_offset <= len);
    }
    // std::cerr << "obs" << n_s_obs << "/" << n_samples << std::endl;
    assert(n_s_obs == n_samples);

    // Restore + update PPA
    // uint32_t n_skips = 0;
    // for (int i = 0; i < n_samples; ++i) { // Worst case O(n), average case O(n) with a smallish constant.
    //     queue[arr[i]][n_queue[arr[i]]++] = ppa[i];
    //     if (arr[i] == 0) {
    //         // ++n_skips;
    //         continue;
    //     }
    //     prev[ppa[i]] = arr[i]; // Unpermute data.
    // }

    // Merge PPA queues.
    uint32_t of = 0;
    for (int j = 0; j < n_symbols; ++j) { // O(n)
        // std::cerr << "queue=" << j << ": " << n_queue[j] << std::endl;
        memcpy(&ppa[of], queue[j], sizeof(uint32_t)*n_queue[j]);
        // for (int i = 0; i < n_queue[j]; ++i, ++of) {
        //     ppa[of] = queue[j][i];
        // }
        of += n_queue[j];
    }
    // std::cerr << "of=" << of << "/" << n_samples << std::endl;
    assert(of == n_samples);
    ++n_steps;

    // std::cerr << "nskips=" << n_skips << std::endl;

    return 1;
}

int PBWT::ReverseUpdateEWAHNm(const uint8_t* arr, const uint32_t len, uint8_t* ret) {
    assert(n_samples > 0);
    assert(n_symbols > 0);

    memset(ret, 0, n_samples); // O(n)
    memset(n_queue, 0, sizeof(uint32_t)*n_symbols);

    // std::cerr << ToPrettyString() << std::endl;

    uint32_t local_offset = 0;
    uint64_t n_s_obs = 0;

    while (local_offset < len) {
        djinn_ewah_t* ewah = (djinn_ewah_t*)&arr[local_offset];
        local_offset += sizeof(djinn_ewah_t);

        // std::cerr << "[PBWT::ReverseUpdateEWAHNm] ewah=" << ewah->ref << "," << ewah->clean << "," << ewah->dirty << std::endl;
        
        // Clean words.
        uint64_t to = n_s_obs + ewah->clean * 8 > n_samples ? n_samples : n_s_obs + ewah->clean * 8;
        // std::cerr << "[PBWT::ReverseUpdateEWAHNm] clean=" << n_s_obs << "->" << to << std::endl; 
        for (int i = n_s_obs; i < to; ++i) {
            queue[ewah->ref & 15][n_queue[ewah->ref & 15]++] = ppa[i];
            // if (ppa[i] == 0) {
                // std::cerr << "ppa[i]=0 -> " << (int)ewah->ref << " for i=" << i << std::endl;;
            // }
            ret[ppa[i]] = (ewah->ref & 15); // update prev when non-zero
        }
        n_s_obs = to;

        // Loop over dirty bitmaps.
        // std::cerr << "dirty=" << ewah->dirty << std::endl;
        for (int i = 0; i < ewah->dirty; ++i) {
            to = n_s_obs + 8 > n_samples ? n_samples : n_s_obs + 8;
            // std::cerr << "dirty step-" << i << "/" << ewah->dirty <<  ": " << n_s_obs << "->" << to << ": " << to-n_s_obs << std::endl;
            assert(n_s_obs < n_samples);
            assert(to <= n_samples);
            
            uint32_t dirty = *((uint32_t*)(&arr[local_offset])); // copy
            // std::cerr << "[PBWT::ReverseUpdateEWAHNm] Dirty: " << i << "/" << ewah->dirty << ": " << std::bitset<32>(dirty) << std::endl;

            for (int j = n_s_obs; j < to; ++j) {
                queue[dirty & 15][n_queue[dirty & 15]++] = ppa[j];
                // if (ppa[j] == 0) {
                //     std::cerr << "ppa[j]=0 -> " << (int)(dirty & 15) << " for j=" << j << std::endl;;
                // }
                ret[ppa[j]] = (dirty & 15); // update prev when non-zero
                // std::cerr << "\t" << std::bitset<4>(dirty&15) << std::endl;
                dirty >>= 4;
                // ++n_s_obs;
            }
            local_offset += sizeof(uint32_t);
            n_s_obs = to;
        }
        // local_offset += ewah->dirty * sizeof(uint32_t);
        // std::cerr << "pbwt offset=" << local_offset << "/" << len << " with=" << n_s_obs << "/" << n_samples << std::endl;
        assert(local_offset <= len);
    }
    // std::cerr << "obs" << n_s_obs << "/" << n_samples << std::endl;
    // if (n_s_obs != n_samples) std::cerr << "[PBWT][ERROR] obs=" << n_s_obs << " expected=" << n_samples << std::endl;
    assert(n_s_obs == n_samples);

    // Restore + update PPA
    // uint32_t n_skips = 0;
    // for (int i = 0; i < n_samples; ++i) { // Worst case O(n), average case O(n) with a smallish constant.
    //     queue[arr[i]][n_queue[arr[i]]++] = ppa[i];
    //     if (arr[i] == 0) {
    //         // ++n_skips;
    //         continue;
    //     }
    //     prev[ppa[i]] = arr[i]; // Unpermute data.
    // }

    // Merge PPA queues.
    uint32_t of = 0;
    for (int j = 0; j < n_symbols; ++j) { // O(n)
        // std::cerr << "[PBWT::ReverseUpdateEWAHNm] queue=" << j << ": " << n_queue[j] << std::endl;
        memcpy(&ppa[of], queue[j], sizeof(uint32_t)*n_queue[j]);
        // for (int i = 0; i < n_queue[j]; ++i, ++of) {
        //     ppa[of] = queue[j][i];
        // }
        of += n_queue[j];
    }
    // std::cerr << "of=" << of << "/" << n_samples << std::endl;
    assert(of == n_samples);
    ++n_steps;

    // std::cerr << "nskips=" << n_skips << std::endl;

    return 1;
}

}