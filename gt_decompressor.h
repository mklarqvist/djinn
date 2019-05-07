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
#ifndef GT_DECOMPRESSOR_H_
#define GT_DECOMPRESSOR_H_

class GenotypeDecompressor {
public:
    GenotypeDecompressor(const uint8_t* in, const int64_t len_in, const int32_t variants, const int64_t n_s) :
        data(in), l_data(len_in),
        n_variants(variants), n_samples(n_s),
        cur_variant(0), cur_offset(0)
    {
        
    }

    virtual ~GenotypeDecompressor() {}

    /**
     * @brief Provide a contiguous section of (uncompressed) data.
     * 
     * @param in 
     * @param len_in 
     * @param variants 
     * @return true 
     * @return false 
     */
    virtual bool SetData(const uint8_t* in, const int64_t len_in, const int32_t variants) {
        if (len_in == 0 || variants == 0) return false;
        data = in;
        l_data = len_in;
        n_variants = variants;
        cur_offset = 0;
        cur_variant = 0;
        return true;
    }

    virtual bool Next() =0;
    virtual bool DecodeNext() =0;
    virtual bool DecodeCurrent() =0;
    virtual bool GetGenotypeArray(uint8_t* data) =0;
    virtual bool GetGenotypeArrayCopy(uint8_t*& data) =0;

public:
    const uint8_t* data; // not owned
    int64_t l_data;
    int32_t n_variants;
    int64_t n_samples;
    int32_t cur_variant;
    int32_t cur_offset;
    // Todo:
    // PBWT models here
};

class GenotypeDecompressorRLEBitmap : public GenotypeDecompressor {
public:
    GenotypeDecompressorRLEBitmap(const uint8_t* in, const int64_t len_in, const int32_t variants, const int64_t n_s) : GenotypeDecompressor(in, len_in, variants, n_s) {}

    bool Next() override { return false; }
    bool DecodeNext() override { return false; }
    bool DecodeCurrent() override { return false; }
    bool GetGenotypeArray(uint8_t* data) override { return false; }
    bool GetGenotypeArrayCopy(uint8_t*& data) override { return false; }

    // Returns the number of alternative alleles if >= 0, otherwise is an error.
    int DecodeRLEBitmap(uint8_t* out) {
        if (out == nullptr) return -1;
        if (cur_variant >= n_variants) return -2;
        if (cur_offset >= l_data) return -3;
        
        uint32_t offset = cur_offset;
        uint32_t n_run = 0;
        uint32_t out_offset = 0;

        uint32_t n_alts = 0;

        while (true) {
            // Looking at most-significant byte for the target compression type.
            const uint8_t type = (data[offset] & 1);
            if (type == 0) {
                assert((*reinterpret_cast<const uint32_t*>(&data[offset]) & 1) == 0);
                n_alts += __builtin_popcount(*reinterpret_cast<const uint32_t*>(&data[offset]) >> 1);

                const uint32_t ulimit = out_offset + 31 >= n_samples ? n_samples - out_offset : 31;
                uint32_t val = *reinterpret_cast<const uint32_t*>(&data[offset]);
                for (int i = 0; i < ulimit; ++i) {
                    val >>= 1;
                    // out[out_offset++] = (val & 1);
                    out[out_offset + (ulimit - i - 1)] = (val & 1); // unpack in reverse order
                }
                out_offset += ulimit;
                offset += sizeof(uint32_t);
                
                n_run  += 31;
                if (n_run >= n_samples) {
                    break;
                }
            }
            else if (type == 1) {
                uint16_t val = *reinterpret_cast<const uint16_t*>(&data[offset]);
                const uint16_t run_length = (val >> 2);
                n_run  += run_length;
                n_alts += ((val >> 1) & 1) * (val >> 2);

                const uint8_t ref = (val >> 1) & 1;
                // std::cerr << "addRLE=" << run_length << "->" << out_offset+run_length << "/" << n_samples << std::endl;
                memset(&out[out_offset], ref, run_length);
                out_offset += run_length;
                assert((*reinterpret_cast<const uint16_t*>(&data[offset]) & 1) == 1);
                offset += sizeof(uint16_t);

                if (n_run >= n_samples) {
                    break;
                }
            }
        }

        cur_offset = offset;
        ++cur_variant;

        // std::cerr << "n_alts=" << n_alts << std::endl;

        return n_alts;
    }

};

#endif
