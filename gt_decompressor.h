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

namespace djinn {

// class GenotypeDecompressor {
// public:
//     GenotypeDecompressor(uint8_t* in, const int64_t len_in, const int32_t variants, const int64_t n_s) :
//         data(in), l_data(len_in),
//         n_variants(variants), n_samples(n_s),
//         cur_variant(0), cur_offset(0)
//     {
        
//     }

//     virtual ~GenotypeDecompressor() {}

//     /**
//      * @brief Provide a contiguous section of (uncompressed) data.
//      * 
//      * @param in 
//      * @param len_in 
//      * @param variants 
//      * @return true 
//      * @return false 
//      */
//     virtual bool SetData(uint8_t* in, const int64_t len_in, const int32_t variants) {
//         if (len_in == 0 || variants == 0) return false;
//         data = in;
//         l_data = len_in;
//         n_variants = variants;
//         cur_offset = 0;
//         cur_variant = 0;
//         return true;
//     }

//     virtual bool Next() =0;
//     virtual bool DecodeNext() =0;
//     virtual bool DecodeCurrent() =0;
//     virtual bool GetGenotypeArray(uint8_t* data) =0;
//     virtual bool GetGenotypeArrayCopy(uint8_t*& data) =0;

// public:
//     uint8_t* data; // not owned
//     int64_t l_data;
//     int32_t n_variants;
//     int64_t n_samples;
//     int32_t cur_variant;
//     int32_t cur_offset;
// };

// class GenotypeDecompressorRLEBitmap : public GenotypeDecompressor {
// public:
//     GenotypeDecompressorRLEBitmap(uint8_t* in, const int64_t len_in, const int32_t variants, const int64_t n_s) : GenotypeDecompressor(in, len_in, variants, n_s) {}
    
//     int InitPbwt(int n_sym) {
//         assert(n_sym > 1);
//         pbwt = std::make_shared<PBWT>(n_samples, n_sym);
//     }

//     bool Next() override { return false; }
//     bool DecodeNext() override { return false; }
//     bool DecodeCurrent() override { return false; }
//     bool GetGenotypeArray(uint8_t* data) override { return false; }
//     bool GetGenotypeArrayCopy(uint8_t*& data) override { return false; }

//     /**
//      * These subroutines return the [0, 1, ..., n] encodings. If bcf-style encoding
//      * is desired then wrap the output through the BcfEncode() subroutine.
//      * 
//      * Todo: Decode2N* for 
//      */

//     // Returns the number of alternative alleles if >= 0, otherwise is an error.
//     int Decode2N2MC(uint8_t* out, const uint32_t stride = 1) {
//         if (out == nullptr) return -1;
//         if (cur_variant >= n_variants) return -2;
//         if (cur_offset >= l_data) return -3;
        
//         uint32_t offset = cur_offset;
//         uint32_t n_run = 0;
//         uint32_t out_offset = 0;
//         // Number of alternative alleles observed.
//         uint32_t n_alts = 0;

//         while (true) {
//             // Looking at most-significant byte for the target compression type.
//             const uint8_t type = (data[offset] & 1);
//             if (type == 0) {
//                 assert((*reinterpret_cast<const uint32_t*>(&data[offset]) & 1) == 0);
//                 n_alts += __builtin_popcount(*reinterpret_cast<const uint32_t*>(&data[offset]) >> 1);

//                 const uint32_t ulimit = out_offset + 31 >= n_samples ? n_samples - out_offset : 31;
//                 uint32_t val = *reinterpret_cast<const uint32_t*>(&data[offset]);
//                 for (int i = 0; i < ulimit; ++i) {
//                     val >>= 1;
//                     out[(out_offset + (ulimit - i - 1))*stride] = (val & 1); // unpack in reverse order
//                 }
//                 out_offset += ulimit;
//                 offset += sizeof(uint32_t);
                
//                 n_run  += 31;
//                 if (n_run >= n_samples) {
//                     break;
//                 }
//             }
//             else if (type == 1) {
//                 uint16_t val = *reinterpret_cast<const uint16_t*>(&data[offset]);
//                 const uint16_t run_length = (val >> 2);
//                 n_run  += run_length;
//                 n_alts += ((val >> 1) & 1) * (val >> 2);

//                 const uint8_t ref = (val >> 1) & 1;
//                 // std::cerr << "addRLE=" << run_length << "->" << out_offset+run_length << "/" << n_samples << std::endl;
//                 memset(&out[out_offset*stride], ref, run_length);
//                 out_offset += run_length;
//                 assert((*reinterpret_cast<const uint16_t*>(&data[offset]) & 1) == 1);
//                 offset += sizeof(uint16_t);

//                 if (n_run >= n_samples) {
//                     break;
//                 }
//             }
//         }

//         cur_offset = offset;
//         ++cur_variant;

//         return n_alts;
//     }

//     int Decode2N2MM(uint8_t* out, const uint32_t stride = 1) {
//         if (out == nullptr) return -1;
//         if (cur_variant >= n_variants) return -2;
//         if (cur_offset >= l_data) return -3;
        
//         uint32_t offset = cur_offset;
//         uint32_t n_run = 0;
//         uint32_t out_offset = 0;
//         // Number of alternative alleles observed.
//         uint32_t n_alts = 0;

//         while (true) {
//             // Looking at most-significant byte for the target compression type.
//             const uint8_t type = (data[offset] & 1);
//             if (type == 0) {
//                 uint32_t val = *reinterpret_cast<const uint32_t*>(&data[offset]);
//                 assert((val & 1) == 0);
//                 // std::cerr << "bitmap " << std::bitset<32>(val) << std::endl;
//                 // n_alts += __builtin_popcount(*reinterpret_cast<const uint64_t*>(&data[offset]) >> 1);

//                 const uint32_t ulimit = out_offset + 15 >= n_samples ? n_samples - out_offset : 15;
//                 val >>= 1;
//                 for (int i = 0; i < ulimit; ++i) {
//                     out[(out_offset + (ulimit - i - 1))*stride] = (val & 3); // unpack in reverse order
//                     val >>= 2;
//                 }
//                 out_offset += ulimit;
//                 offset += sizeof(uint32_t);
                
//                 n_run  += 15;
//                 if (n_run >= n_samples) {
//                     break;
//                 }
//             }
//             else if (type == 1) {
//                 uint16_t val = *reinterpret_cast<const uint16_t*>(&data[offset]);
//                 assert((val & 1) == 1);
//                 const uint16_t run_length = (val >> 3);
//                 // std::cerr << "rle " << run_length << ":" << ((val >> 1) & 3) << std::endl;
//                 n_run  += run_length;
//                 // n_alts += ((val >> 1) & 1) * (val >> 2);

//                 const uint8_t ref = (val >> 1) & 3;
//                 // std::cerr << "addRLE=" << run_length << "->" << out_offset+run_length << "/" << n_samples << std::endl;
//                 memset(&out[out_offset*stride], ref, run_length);
//                 out_offset += run_length;
//                 offset += sizeof(uint16_t);

//                 if (n_run >= n_samples) {
//                     break;
//                 }
//             }
//         }

//         cur_offset = offset;
//         ++cur_variant;

//         return 1;
//     }

//     int Decode2NXM(uint8_t* out, const uint32_t stride = 1) {
//         if (out == nullptr) return -1;
//         if (cur_variant >= n_variants) return -2;
//         if (cur_offset >= l_data) return -3;
        
//         uint32_t offset = cur_offset;
//         uint32_t n_run = 0;
//         uint32_t out_offset = 0;
//         // Number of alternative alleles observed.
//         uint32_t n_alts = 0;

//         while (true) {
//             // Looking at most-significant byte for the target compression type.
//             const uint8_t type = (data[offset] & 1);
//             if (type == 0) {
//                 uint64_t val = *reinterpret_cast<const uint64_t*>(&data[offset]);
//                 assert((val & 1) == 0);
//                 // std::cerr << "bitmap " << std::bitset<64>(val) << std::endl;
//                 // n_alts += __builtin_popcount(*reinterpret_cast<const uint64_t*>(&data[offset]) >> 1);

//                 const uint64_t ulimit = out_offset + 15 >= n_samples ? n_samples - out_offset : 15;
//                 val >>= 1;
//                 for (int i = 0; i < ulimit; ++i) {
//                     out[(out_offset + (ulimit - i - 1))*stride] = (val & 15); // unpack in reverse order
//                     val >>= 4;
//                 }
//                 out_offset += ulimit;
//                 offset += sizeof(uint64_t);
                
//                 n_run  += 15;
//                 if (n_run >= n_samples) {
//                     break;
//                 }
//             }
//             else if (type == 1) {
//                 uint32_t val = *reinterpret_cast<const uint32_t*>(&data[offset]);
//                 assert((val & 1) == 1);
//                 const uint32_t run_length = (val >> 5);;
//                 // std::cerr << "rle " << run_length << ":" << ((val >> 1) & 15) << std::endl;
//                 n_run  += run_length;
//                 // n_alts += ((val >> 1) & 1) * (val >> 2);

//                 const uint8_t ref = (val >> 1) & 15;
//                 // std::cerr << "addRLE=" << run_length << "->" << out_offset+run_length << "/" << n_samples << std::endl;
//                 memset(&out[out_offset*stride], ref, run_length);
//                 out_offset += run_length;
//                 offset += sizeof(uint32_t);

//                 if (n_run >= n_samples) {
//                     break;
//                 }
//             }
//         }

//         cur_offset = offset;
//         ++cur_variant;

//         return 1;
//     }

//     std::shared_ptr<PBWT> pbwt;
// };

// class GenotypeDecompressorContext : public GenotypeDecompressor {
// public:
//     GenotypeDecompressorContext(uint8_t* in, const int64_t len_in, 
//         const int32_t variants, const int64_t n_s, const uint32_t n_sym) : 
//             GenotypeDecompressor(in, len_in, variants, n_s),
//             in_partition(nullptr), len_in_partition(0), partition_size(16) 
//     {
//         model = std::make_shared<GeneralPBWTModel>();
//         if (n_sym) model->Construct(n_samples, n_sym);
//         model->StartDecoding(in);
//     }

//     GenotypeDecompressorContext(uint8_t* in, const int64_t len_in, 
//         uint8_t* in_partition, const int64_t len_in_partition, 
//         const int32_t variants, const int64_t n_s, const uint32_t n_sym) : 
//             GenotypeDecompressor(in, len_in, variants, n_s),
//             in_partition(in_partition), len_in_partition(len_in_partition), 
//             partition_size(16)
//     {
//         model = std::make_shared<GeneralPBWTModel>();
//         if (n_sym) model->Construct(n_samples, n_sym);
//         model->StartDecoding(in);
//         partition = std::make_shared<GeneralPBWTModel>();
//         partition->Construct(1,2);
//         partition->StartDecoding(in_partition);
//     }
    
//     ~GenotypeDecompressorContext() {}

//     int InitPbwt(int n_sym) {
//         assert(n_sym > 1);
//         if (model.get() != nullptr) model->Construct(n_samples, n_sym);
//     }

//     bool Next() override { return false; }
//     bool DecodeNext() override { return false; }
//     bool DecodeCurrent() override { return false; }
//     bool GetGenotypeArray(uint8_t* data) override { return false; }
//     bool GetGenotypeArrayCopy(uint8_t*& data) override { return false; }

//     int Decode(uint8_t* out) {
//         if (out == nullptr) return -1;

//         partition->ResetContext();
//         model->ResetContext();
//         memset(out, 0, n_samples);

//         const uint32_t step_size = std::ceil((float)n_samples / partition_size);
//         uint32_t l_offset = 0;

//         for (int j = 0; j < partition_size; ++j) {
//             uint16_t ret = partition->DecodeSymbol();
//             uint32_t lim = l_offset + step_size > n_samples ? n_samples - l_offset : step_size;
//             if (ret) {
//                 uint32_t nonzero = 0;
//                 for (int k = 0; k < lim; ++k) {
//                     int retgt = model->DecodeSymbol();
//                     nonzero += retgt;
//                     out[l_offset + k] = retgt;
//                 }
//                 assert(nonzero > 0);
//             }
//             l_offset += lim;
//         }
//         assert(l_offset == n_samples);

//         return 1;
//     }

//     int Decode2(uint8_t* out) {
//         if (out == nullptr) return -1;

//         model->ResetContext();
//         memset(out, 0, n_samples);

//         for (int j = 0; j < n_samples; ++j) {
//             out[j] = model->DecodeSymbol();
//         }

//         return 1;
//     }

// public:
//     uint8_t* in_partition;
//     int64_t len_in_partition;
//     uint32_t partition_size;
//     std::shared_ptr<GeneralPBWTModel> model;
//     std::shared_ptr<GeneralPBWTModel> partition;
// };

}

#endif
