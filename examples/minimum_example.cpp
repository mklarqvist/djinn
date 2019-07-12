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
#include <fstream> // Support for read/write.
#include <bitset> // for std::bitset
#include <djinn.h> // Djinn data models.

/**
 * In this example we will generate random data and import it into Djinn. As an example,
 * we will generate random genotypes for 100,000 samples at 500,000 sites . Using a
 * diploid model this will result in 50,000 diploid genotypes.
 */
int main() {
    const uint32_t n_samples = 250000;
    uint8_t* rand_vec = new uint8_t[n_samples]; // 250k haplotypes
    djinn::djinn_model* enc_model = new djinn::djinn_ewah_model(djinn::CompressionStrategy::LZ4, 1);
    enc_model->StartEncoding(true, true); // use PBWT, reset models.
    
    memset(rand_vec, 0, n_samples); // Set all alleles to reference (0)
    rand_vec[21] = 1;
    rand_vec[10048] = 1;
    // Encode this data array.
    int ret = enc_model->Encode(rand_vec, n_samples, 2, 2);

    memset(rand_vec, 0, n_samples); // Set all alleles to reference (0)
    rand_vec[2112]  = 1;
    rand_vec[7182]  = 1;
    rand_vec[92817] = 1;
    rand_vec[92818] = 1;

    // Encode this data array.
    ret = enc_model->Encode(rand_vec, n_samples, 2, 2);

    // Encode a site with three (3) alleles.
    memset(rand_vec, 0, n_samples); // Set all alleles to reference (0).
    rand_vec[1239]  = 1; // First alternative allele.
    rand_vec[1922]  = 1;
    rand_vec[92817] = 2; // Second alternative allele.
    rand_vec[92818] = 1;

    // Encode this data array.
    ret = enc_model->Encode(rand_vec, n_samples, 2, 3);

    // Finish encoding.
    enc_model->FinishEncoding();

    // Allocate memory to store serialized object.
    uint8_t* buffer = new uint8_t[10000];
    int serial_size = enc_model->Serialize(buffer);

    // Print out compression ratio compared to uBCF.
    std::cerr << "[COMPRESS] In uBCF: " << 3*n_samples << "->" << serial_size 
              << " (" << (double)(3*n_samples)/serial_size << "-fold)" << std::endl;

    delete enc_model;

    djinn::djinn_model* dec_model = new djinn::djinn_ewah_model();
    dec_model->Deserialize(buffer); // Read data back into model.
    dec_model->StartDecoding(); // Start decoding the data. Calling this function
          // is required prior to calling any downstream decoding functions.

    // This data structure is the core return component in Djinn. Depending
    // on the decompression context (DecodeNext or DecodeNextRaw) this
    // component will be populated with either a Bcf-styled byte array or
    // a byte array of EWAH objects (see below). Data allocation is managed
    // internally during decompression and automatically reused if possible.
    djinn::djinn_variant_t* variant = nullptr;

    // Simple use case: iterate over variants in a data block and return
    // the bit-exact data. This approach is generally much slower compared
    // to operating directly in EWAH-space but is considerably simpler.
    std::cerr << "There are " << dec_model->n_variants << " variants in the data block" << std::endl;
    for (int i = 0; i < dec_model->n_variants; ++i) {
        // The DecodeNext function decompress and decode the next variant site.
        // This procedure will unpermute the PBWT-based preprocessor, if necessary.
        // Therefore, the returned data is bit-exact to the input data.
        int dec_ret = dec_model->DecodeNext(variant);
        // The returned data array is always equal to the length of the input
        // data (number of samples). Here we assert that this is true.
        assert(variant->data_len == n_samples);
        // Print out the sample numbers that do not encode the reference allele.
        // These positions are identifical to the encoded positions above.
        //
        // Accessing data in this format is trivial:
        // variant->data[j] is the j-th allele in original order.
        for (int j = 0; j < n_samples; ++j) {
            if (variant->data[j] != 0) 
                std::cerr << j << ":" << (char)(variant->data[j] + '0') << ",";
        }
        std::cerr << std::endl;
    }

    delete dec_model;

    // Start over for advanced use case.
    dec_model = new djinn::djinn_ewah_model();
    dec_model->Deserialize(buffer); // Read data back into model.
    dec_model->StartDecoding(); // Start decoding the data. Calling this function
          // is required prior to calling any downstream decoding functions.

    // Advanced use case by operating directly on EWAH-encoded data. Operating
    // in this space is considerably faster but require additional effort.
    std::cerr << "There are " << dec_model->n_variants << " variants in the data block" << std::endl;
    for (int i = 0; i < dec_model->n_variants; ++i) {
        int dec_ret = dec_model->DecodeNextRaw(variant);
        std::cerr << "Variant-" << i+1 << "/" << dec_model->n_variants << " has " 
           << variant->ploidy << " base ploidy and " << variant->n_allele << " different alleles" << std::endl;
        // Raw decoding updates the djinn_variant_t structure and properly
        // populates the unpacked support structure djn_variant_dec_t
        // accessibly at variant->d. This additional information includes
        // the number of EWAH objects returned and pointers to them and the
        // (potentially) matched dirty words. The variable variant->unpacked
        // will be set to DJN_UN_EWAH when returning EWAH-compressed data.
        //
        // It is **very** important to understand that _if_ the input data
        // was encoded using the PBWT preprocessor then the returned data
        // is in the **permuted** order. This means that the sample order
        // will be (potentially) different at each variant site. In many
        // applications this does not matter and operating directly 
        // at this compressed level result in considerable performance 
        // gains across the board.
        if (variant->unpacked == DJN_UN_EWAH) {
            std::cerr << "EWAH objects: " << variant->d->n_ewah << std::endl;
            for (int j = 0; j < variant->d->n_ewah; ++j) {
                std::cerr << variant->d->ewah[j]->clean << " words with reference " << variant->d->ewah[j]->ref
                << " followed by " << variant->d->ewah[j]->dirty << " words." << std::endl;
                std::cerr << "Dirty: ";
                for (int k = 0; k < variant->d->ewah[j]->dirty; ++k) {
                    // Pointers are stored to the first dirty bitmap of each
                    // EWAH object such that dirty[j] is the first pointer
                    // and dirty[j]+1 is the second dirty word for EWAH[j],
                    // etc.
                    std::cerr << std::bitset<32>(*(variant->d->dirty[j]+k));
                }
                std::cerr << std::endl;
            }
        }
        // A **very** importent concept to understand is that there are 
        // situations where there will be *more* alleles in the enocodings 
        // than there are samples. This happens because the smallest word
        // size we can use is 32-bits wide. if your sample number is NOT
        // divisible by 32 or 8 then there will be ceil(#samples/32) 
        // encoded samples.
        //
        // In this example there will be 250016 encoded alleles for the 
        // 250000 input samples. It is up to the user to make sure that
        // decoding the EWAH data terminates when reaching the expected
        // output values.
        //
        // Another **very** important concept is that EWAH encodings may
        // be returned as either 1-bit or 4-bit encodings depending on how
        // many alternative alleles are present at the given locus. Which
        // bit-width-encoding is used can be checked at variant->d->dirty_type
        // that is set to one of the DJN_DIRTY_* values. The default machine
        // word is 32-bits and therefore contain either 32/1 = 32 alleles
        // or 32/4 = 8 alleles per word.
        //
        // For the last variant in this example, there will not be any excess
        // alleles stores as 250000 mod 8 is 0.
        uint32_t n_obs = 0;
        const uint32_t factor = (variant->d->dirty_type == DJN_DIRTY_2MC ? 32 : 8);
        for (int j = 0; j < variant->d->n_ewah; ++j) {
            n_obs += (variant->d->ewah[j]->clean + variant->d->ewah[j]->dirty) * factor;
        }
        std::cerr << "Observed alleles: " << n_obs << "/" << n_samples << std::endl;
    }

    // Cleanup.
    delete[] buffer;
    delete[] rand_vec;
    delete dec_model;

    return EXIT_SUCCESS;
}