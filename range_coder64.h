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
#ifndef RANGE_CODER64_H_
#define RANGE_CODER64_H_

#include <cstdint>
#include <cstddef>

#include <iostream>
#include <cassert>

namespace djinn {

/*	Code for range coding, derived from public domain work by Dmitry Subbotin
    Modified to use 64-bit integer maths, for increased precision

    Note :	Cannot be used at full 'capacity' as the interface still takes DWord parameters (not QWord)
            This is done to maintain uniformity in interface across all entropy coders, feel free to 
            change this.
    
    author : Sachin Garg
*/
class RangeCoder64
{
protected:
    RangeCoder64();
    static constexpr uint64_t Top = (uint64_t)1<<56;
    static constexpr uint64_t Bottom = (uint64_t)1<<48;

    uint64_t Low, Range;

public:
    static constexpr uint64_t MaxRange = (uint64_t)1<<48;
};

class RangeEncoder64 : public RangeCoder64
{
private:
    bool Flushed;

public:
    RangeEncoder64();
    ~RangeEncoder64();

    void EncodeRange(uint8_t*& pos, uint32_t SymbolLow, uint32_t SymbolHigh, uint32_t TotalRange);
    void Flush(uint8_t*& pos);
};

class RangeDecoder64 : public RangeCoder64
{
private:
    uint64_t Code;

public:
    RangeDecoder64();

    uint32_t GetCurrentCount(uint32_t TotalRange);
    void RemoveRange(uint32_t SymbolLow, uint32_t SymbolHigh, uint32_t TotalRange);
};

class SimpleModel64 {
public:
    struct s64_model {
        s64_model() {
            for(int i = 0; i < 4; ++i) Freq[i] = i;
        }

        void reset() {
            for(int i = 0; i < 4; ++i) Freq[i] = i;
        }

        uint32_t Freq[4];
    };

    struct s64_model_rle {
        s64_model_rle() {
            for(int i = 0; i < 257; ++i) Freq[i] = i;
        }

        void reset() {
            for(int i = 0; i < 257; ++i) Freq[i] = i;
        }

        uint32_t Freq[257];
    };

public:
    SimpleModel64() : 
        buf(new uint8_t[10000000]), pos(buf), l_buf(0), 
        context(0), models(new s64_model[128]), 
        rle_ref(-1), n_rle(0), rle_buf(new uint8_t[10000000]), l_rle_buf(0),
        context_rle(0), models_rle(new s64_model_rle[65536]) 
    {
        // for(int i = 0; i < 257; ++i) Freq[i] = i;
    }

    ~SimpleModel64() {
        delete[] buf;
        delete[] models;
        delete[] rle_buf;
    }

    static inline 
    void Rescale(uint32_t* Frequency) {
        for(int i = 1; i < 4; ++i) {
            Frequency[i] /= 2; // Frequency[i] >>= 1
            if (Frequency[i] <= Frequency[i-1])
                Frequency[i] = Frequency[i-1] + 1;
        }
    }

    void Encode(uint8_t symbol) {
        if (rle_ref == -1) rle_ref = symbol;

        if (symbol != rle_ref) {
            if (n_rle >= 16) {
                // std::cerr << "rle@" << n_rle << "|" << (int)rle_ref << std::endl;

                // Encode RLE as multiples of 16
                uint32_t special = n_rle / 16;
                uint16_t residual = n_rle % 16;
                // std::cerr << n_rle << " special=" << (int)special << " and " << residual << std::endl;

                while (special) {
                    Encode2(2); // escape symbol
                    Encode2(rle_ref); // ref
                    uint32_t enc = special <= 15 ? special : 15;
                    EncodeRLE(enc);
                    special -= enc;
                }
                for (int i = 0; i < residual; ++i) {
                    Encode2(rle_ref);
                }


                // if (n_rle < 127) {
                //     Encode2(2); // escape symbol
                //     uint8_t val = (n_rle << 1) | (rle_ref & 1);
                //     rle_buf[l_rle_buf] = val;
                //     l_rle_buf += sizeof(uint8_t);
                // } else {
                
                //     while (n_rle) {
                //         uint16_t n_rle_emit = n_rle < 32767 ? n_rle : 32767;
                //         Encode2(2); // escape symbol
                //         uint16_t val = (n_rle_emit << 1) | (rle_ref & 1);
                //         *reinterpret_cast<uint16_t*>(&rle_buf[l_rle_buf]) = val;
                //         l_rle_buf += sizeof(uint16_t);
                //         n_rle -= n_rle_emit;
                //     }
                // }
                
            } else {
                assert(n_rle != 0);
                // std::cerr << "Encoding=" << n_rle << "|" << (int)rle_ref << std::endl;
                for (int i = 0; i < n_rle; ++i) {
                    Encode2(rle_ref);
                }
            }
            n_rle = 0;
            rle_ref = symbol;
        }

        ++n_rle;
    }

    void Encode2(uint8_t symbol) {
        uint32_t* F = models[context].Freq;

        context <<= 2;
        context |= (symbol & 3);
        context &= (128-1);

        encoder.EncodeRange(pos, F[symbol], F[symbol+1], F[3]);
        for(int j = symbol+1; j < 4; ++j) {
            ++F[j];
        } 
        
        if (F[3] >= RangeEncoder64::MaxRange) {
            std::cerr << "rescale" << std::endl;
            Rescale(F);
        }
    }

    void EncodeRLE(uint8_t symbol) {
        assert(symbol < 256);
        uint32_t* F = models_rle[context].Freq;

        context_rle <<= 8;
        context_rle |= (symbol & 255);
        context_rle &= (65536 - 1);

        encoder_rle.EncodeRange(pos, F[symbol], F[symbol+1], F[256]);
        for(int j = symbol+1; j < 257; ++j) {
            ++F[j];
        } 
        
        if (F[256] >= RangeEncoder64::MaxRange) {
            std::cerr << "rescale" << std::endl;
            Rescale(F);
        }
    }

    void Flush() {
        if (n_rle >= 16) {
            // std::cerr << "rle@" << n_rle << "|" << (int)rle_ref << std::endl;
            
            // Encode RLE as multiples of 16
            uint32_t special = n_rle / 16;
            uint16_t residual = n_rle % 16;
            // std::cerr << n_rle << " special=" << (int)special << " and " << residual << std::endl;

            while (special) {
                Encode2(2); // escape symbol
                Encode2(rle_ref); // ref
                uint32_t enc = special <= 255 ? special : 255;
                EncodeRLE(enc);
                special -= enc;
            }
            for (int i = 0; i < residual; ++i) {
                Encode2(rle_ref);
            }

            // if (n_rle < 127) {
            //     Encode2(2); // escape symbol
            //     uint8_t val = (n_rle << 1) | (rle_ref & 1);
            //     rle_buf[l_rle_buf] = val;
            //     l_rle_buf += sizeof(uint8_t);
            // } else {
            
            //     while (n_rle) {
            //         uint16_t n_rle_emit = n_rle < 32767 ? n_rle : 32767;
            //         Encode2(3); // escape symbol
            //         uint16_t val = (n_rle_emit << 1) | (rle_ref & 1);
            //         *reinterpret_cast<uint16_t*>(&rle_buf[l_rle_buf]) = val;
            //         l_rle_buf += sizeof(uint16_t);
            //         n_rle -= n_rle_emit;
            //     }
            // }
            
        } else {
            if (n_rle) {
                assert(n_rle != 0);
                assert(rle_ref == 0 || rle_ref == 1);
                // std::cerr << "Encoding=" << n_rle << "|" << (int)rle_ref << std::endl;
                for (int i = 0; i < n_rle; ++i) {
                    Encode2(rle_ref);
                }
            }
        }
        
        encoder.Flush(pos);
    }

    void reset() { 
        pos = buf; 
        l_buf = 0; 
        for (int i = 0; i < 128; ++i) models[i].reset();
        l_rle_buf = 0;
        n_rle = 0; rle_ref = -1;
        for (int i = 0; i < 65536; ++i) models_rle[i].reset();
        context_rle = 0;
    }

    size_t size() { return pos - buf; }

public:
    uint8_t* buf;
    uint8_t* pos; // pointer offset to buf
    size_t l_buf;

    int8_t rle_ref;
    uint16_t n_rle;
    uint8_t* rle_buf;
    size_t l_rle_buf; 

    uint32_t context;
    s64_model* models;

    uint32_t context_rle;
    s64_model_rle* models_rle;

    // uint32_t Freq[257];
    RangeEncoder64 encoder;
    RangeEncoder64 encoder_rle;
};

}

#endif