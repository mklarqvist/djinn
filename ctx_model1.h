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

#include <limits>

#include "djinn.h"
#include "pbwt.h"

namespace djinn {

template<typename T>
uint32_t round_log2(T x)
{
	uint32_t r = 0;

	for (; x; ++r)
		x >>= 1;

	return r;
}

/*
EncodeBcf(); // input data formatted according to bcf spec
Encode() // any input data properly formatted
Decode() // get a proper variant back
DecodeRaw() // get uncompressed but encoded and unpermuted data back

StartEncoding()
StopEncoding()
StartDecoding()
StopDecoding() // does nothing
*/

struct djinn_ctx_model_t {
    djinn_ctx_model_t() :
        use_pbwt(true), init(false), unused(0),
        p(nullptr),
        p_len(0), p_cap(0), p_free(false),
        n_variants(0)
    {
       
    }

    ~djinn_ctx_model_t() {
        if (p_free) delete[] p;
    }

    void Initiate2mc() {
        range_coder = std::make_shared<RangeCoder>();
        mref = std::make_shared<GeneralModel>(2, 4096, range_coder);
        mlog_rle = std::make_shared<GeneralModel>(32, 4096, range_coder);
        mrle = std::make_shared<GeneralModel>(256, 1024, 24, 1, range_coder);
        mrle2_1 = std::make_shared<GeneralModel>(256, 1024, range_coder);
        mrle2_2 = std::make_shared<GeneralModel>(256, 1024, range_coder);
        mrle4_1 = std::make_shared<GeneralModel>(256, 1024, range_coder);
        mrle4_2 = std::make_shared<GeneralModel>(256, 1024, range_coder);
        mrle4_3 = std::make_shared<GeneralModel>(256, 1024, range_coder);
        mrle4_4 = std::make_shared<GeneralModel>(256, 1024, range_coder);
        dirty_wah = std::make_shared<GeneralModel>(256, 65536, 24, 1, range_coder);
        mtype = std::make_shared<GeneralModel>(2, 1024, range_coder);
    }

    void InitiateNm() {
        range_coder = std::make_shared<RangeCoder>();
        mref = std::make_shared<GeneralModel>(16, 4096, range_coder); // order-3
        mlog_rle = std::make_shared<GeneralModel>(32, 4096, range_coder);
        mrle = std::make_shared<GeneralModel>(256, 256, 24, 1, range_coder);
        mrle2_1 = std::make_shared<GeneralModel>(256, 256, range_coder);
        mrle2_2 = std::make_shared<GeneralModel>(256, 256, range_coder);
        mrle4_1 = std::make_shared<GeneralModel>(256, 256, range_coder);
        mrle4_2 = std::make_shared<GeneralModel>(256, 256, range_coder);
        mrle4_3 = std::make_shared<GeneralModel>(256, 256, range_coder);
        mrle4_4 = std::make_shared<GeneralModel>(256, 256, range_coder);
        dirty_wah = std::make_shared<GeneralModel>(256, 65536, 24, 1, range_coder);
        mtype = std::make_shared<GeneralModel>(2, 1024, range_coder);
    }

    void reset() {
        if (mref.get() != nullptr) mref->Reset();
        if (mlog_rle.get() != nullptr) mlog_rle->Reset();
        if (mrle.get() != nullptr) mrle->Reset();
        if (mrle2_1.get() != nullptr) mrle2_1->Reset(); 
        if (mrle2_2.get() != nullptr) mrle2_2->Reset();
        if (mrle4_1.get() != nullptr) mrle4_1->Reset(); 
        if (mrle4_2.get() != nullptr) mrle4_2->Reset(); 
        if (mrle4_3.get() != nullptr) mrle4_3->Reset(); 
        if (mrle4_4.get() != nullptr) mrle4_4->Reset();
        if (dirty_wah.get() != nullptr) dirty_wah->Reset();
        if (mtype.get() != nullptr) mtype->Reset();
        if (mref.get() != nullptr) mref->Reset();
        pbwt.reset();
        n_variants = 0; p_len = 0;
    }

    int StartEncoding(bool use_pbwt, bool reset = false) {
        if (range_coder.get() == nullptr) return -1;
        if (p_cap == 0) { // initiate a buffer if there is none
            delete[] p;
            p = new uint8_t[10000000];
            p_len = 0;
            p_cap = 10000000;
            p_free = true;
        }
        if (reset) this->reset();
        this->use_pbwt = use_pbwt;
        range_coder->SetOutput(p);
        range_coder->StartEncode();
        return 1;
    }

    size_t FinishEncoding() {
        if (range_coder.get() == nullptr) return -1;
        range_coder->FinishEncode();
        return range_coder->OutSize();
    }

    int StartDecoding(uint8_t* data, bool reset = false) {
        // If resetting the model_2mc.
        if (range_coder.get() == nullptr) return -1;
        if (reset) this->reset();

        range_coder->SetInput(data);
        range_coder->StartDecode();
        return 1;
    }

    uint8_t use_pbwt: 1, init: 1, unused: 6;
    PBWT pbwt;
    std::shared_ptr<RangeCoder> range_coder; // shared range coder
    std::shared_ptr<GeneralModel> mref, mlog_rle, mrle, mrle2_1, mrle2_2, mrle4_1, mrle4_2, mrle4_3, mrle4_4; // rle models
    std::shared_ptr<GeneralModel> dirty_wah; // dirty word model
    std::shared_ptr<GeneralModel> mtype; // mtype model
    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1;
    uint32_t n_variants;
};

class djinn_ctx_model {
public:
    djinn_ctx_model() : 
        n_samples(0), 
        n_samples_wah(0), n_wah(0), wah_bitmaps(nullptr), 
        p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
        range_coder(std::make_shared<RangeCoder>()), 
        marchetype(std::make_shared<GeneralModel>(2, 1024, range_coder))
    {
        model_2mc.Initiate2mc();
        model_nm.InitiateNm();
    }

    djinn_ctx_model(uint64_t n_s) : 
        n_samples(n_s), 
        n_wah(std::ceil((float)2*n_samples / 64)), 
        n_samples_wah(n_wah*64), 
        wah_bitmaps(new uint64_t[n_wah]), 
        p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
        range_coder(std::make_shared<RangeCoder>()), 
        marchetype(std::make_shared<GeneralModel>(2, 1024, range_coder))
    {
        model_2mc.Initiate2mc();
        model_nm.InitiateNm();
    }
    
    ~djinn_ctx_model() { 
        delete[] wah_bitmaps;
        if (p_free) delete[] p;
    }

    void SetSamples(int64_t n_s) { 
        n_samples = n_s;
        n_wah = std::ceil((float)n_samples / 64);
        n_samples_wah = n_wah*64;
    } 

    int EncodeBcf(uint8_t* data, uint8_t alt_alleles) {
        if (n_samples == 0) return -1;
        if (data == nullptr) return -2;

        // uint64_t alts = 0;
        uint32_t hist_alts[256] = {0};
        for (int i = 0; i < n_samples; ++i) {
            ++hist_alts[BCF_UNPACK_GENOTYPE_GENERAL(data[i])];
        }

        if (alt_alleles <= 2 && hist_alts[14] == 0 && hist_alts[15] == 0) { // does not check for missingness
            

            if (hist_alts[1] < 10) { // dont update if < 10 alts
                for (int i = 0; i < n_samples; ++i) {
                    model_2mc.pbwt.prev[i] = BCF_UNPACK_GENOTYPE(data[model_2mc.pbwt.ppa[i]]);
                }
            } else {
                model_2mc.pbwt.Update(data, 1);
            }

            marchetype->EncodeSymbol(0); // add archtype as 2mc
            return Encode2mc(model_2mc.pbwt.prev, n_samples);
        } else {
            model_nm.pbwt.UpdateGeneral(data, 1); // otherwise
            marchetype->EncodeSymbol(1);
            return EncodeNm(model_nm.pbwt.prev, n_samples);
        }
    }

    int Encode2mc(uint8_t* data, uint32_t len) {
        if (data == nullptr) return -1;
        if (n_samples == 0) return -2;
        
        if (wah_bitmaps == nullptr) {
            n_wah = std::ceil((float)n_samples*4 / 64); // up to 4 bits per allele
            wah_bitmaps = new uint64_t[n_wah];
        }
        
        memset(wah_bitmaps, 0, n_wah*sizeof(uint64_t));

        for (int i = 0; i < n_samples; ++i) {
            if (data[i]) {
                wah_bitmaps[i / 64] |= 1L << (i % 64);
            }
        }

        return EncodeWah(wah_bitmaps, n_wah/4);
    }

    int EncodeNm(uint8_t* data, uint32_t len) {
        if (data == nullptr) return -1;
        if (n_samples == 0) return -2;
        
        if (wah_bitmaps == nullptr) {
            n_wah = std::ceil((float)n_samples*4 / 64);
            wah_bitmaps = new uint64_t[n_wah];
        }
        
        memset(wah_bitmaps, 0, n_wah*sizeof(uint64_t));

        for (int i = 0; i < n_samples; ++i) {
            // std::cerr << (int)data[i] << " ";
            // std::cerr <<(4*(i % 16)) << std::endl;
            wah_bitmaps[i / 16] |= (uint64_t)data[i] << (4*(i % 16));
        }

        // for (int i = 0; i < n_wah; ++i) std::cerr << std::bitset<64>(wah_bitmaps[i]) << " ";
        // std::cerr << std::endl;

        return EncodeWahNm(wah_bitmaps, n_wah);
        // return -1;
    }

    // Return raw, potentially permuted, WAH-like encoding + archetype dict
    int DecodeRaw(uint8_t* data) {
        if (data == nullptr) return -1;

        int64_t n_samples_obs = 0;
        while(true) {
            uint8_t type = model_2mc.mtype->DecodeSymbol();
            // std::cerr << "Type=" << (int)type << std::endl;
            if (type == 0) { // bitmaps
                uint64_t wah = 0;
                // std::cerr << "Bitmap=";
                for (int i = 0; i < 8; ++i) {
                    wah |= model_2mc.dirty_wah->DecodeSymbol();
                    // std::cerr << std::bitset<8>(wah&255);
                    wah <<= 8;
                }
                // std::cerr << std::endl;
                n_samples_obs += 64;

            } else { // is RLE
                // Decode an RLE
                uint64_t ref = 1; uint32_t len = 1;
                DecodeWahRLE(ref, len);
                // std::cerr << "decoded RLE=" << ref << " len=" << len << std::endl;
                n_samples_obs += len*64;
            }

            if (n_samples_obs == n_samples_wah) {
                std::cerr << "obs=" << n_samples_obs << "/" << n_samples_wah << std::endl;
                break;
            }
            if (n_samples_obs > n_samples_wah) {
                std::cerr << "Decompression corruption!" << std::endl;
                exit(1);
            }
        }

        return 1;
    }

    int EncodeWah(uint64_t* wah, uint32_t len) { // input WAH-encoded data
        if (wah == nullptr) return -1;

        ++model_2mc.n_variants;

        uint64_t wah_ref = wah[0];
        uint64_t wah_run = 1;

        for (int i = 1; i < len; ++i) {
            if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint64_t>::max()) || (wah_ref != wah[i])) {
                if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint64_t>::max()) || wah_run == 1) {
                    model_2mc.mtype->EncodeSymbol(0);
                    
                    // std::cerr << "Bitmap=" << std::bitset<64>(wah_ref) << std::endl;
                    for (int i = 0; i < 8; ++i) {
                        model_2mc.dirty_wah->EncodeSymbol(wah_ref & 255);
                        wah_ref >>= 8;
                    }
                } else {
                    model_2mc.mtype->EncodeSymbol(1);
                    EncodeWahRLE(wah_ref, wah_run, &model_2mc);
                }
                wah_run = 0;
                wah_ref = wah[i];
            }
            ++wah_run;
        }

        if (wah_run) {
            if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint64_t>::max()) || wah_run == 1) {
                model_2mc.mtype->EncodeSymbol(0);
                
                // std::cerr << "Bitmap=" << std::bitset<64>(wah_ref) << std::endl;
                for (int i = 0; i < 8; ++i) {
                    model_2mc.dirty_wah->EncodeSymbol(wah_ref & 255);
                    wah_ref >>= 8;
                }
                
            } else {
                model_2mc.mtype->EncodeSymbol(1);
                EncodeWahRLE(wah_ref, wah_run, &model_2mc);
            }
        }

        return 1;
    }

    int EncodeWahNm(uint64_t* wah, uint32_t len) { // input WAH-encoded data
        if (wah == nullptr) return -1;

        ++model_nm.n_variants;

        uint64_t wah_ref = wah[0];
        uint64_t wah_run = 1;

        for (int i = 1; i < len; ++i) {
            if ((wah_ref != 1229782938247303441ULL) || (wah_ref != wah[i])) {// 000100010001
                if ((wah_ref != wah_bitmaps[i]) || wah_run == 1) {
                    model_nm.mtype->EncodeSymbol(0);
                    
                    // std::cerr << "Bitmap=" << std::bitset<64>(wah_ref) << std::endl;
                    for (int i = 0; i < 8; ++i) {
                        model_nm.dirty_wah->EncodeSymbol(wah_ref & 255);
                        wah_ref >>= 8;
                    }
                } else {
                    model_nm.mtype->EncodeSymbol(1);
                    EncodeWahRLE_nm(wah_ref, wah_run, &model_nm);
                    // std::cerr << "RLE=" << wah_ref << "," << wah_run << std::endl;
                }
                wah_run = 0;
                wah_ref = wah[i];
            }
            ++wah_run;
        }

        if (wah_run) {
            if ((wah_ref != 1229782938247303441ULL) || wah_run == 1) {
                model_nm.mtype->EncodeSymbol(0);
                
                // std::cerr << "Bitmap=" << std::bitset<64>(wah_ref) << std::endl;
                for (int i = 0; i < 8; ++i) {
                    model_nm.dirty_wah->EncodeSymbol(wah_ref & 255);
                    wah_ref >>= 8;
                }
                
            } else {
                model_nm.mtype->EncodeSymbol(1);
                EncodeWahRLE_nm(wah_ref, wah_run, &model_nm);
                // std::cerr << "RLE=" << wah_ref << "," << wah_run << std::endl;
            }
        }

        return 1;
    }

    int EncodeWahRLE(uint64_t ref, uint32_t len, djinn_ctx_model_t* model) {
        model->mref->EncodeSymbol(ref&1);
        uint32_t log_length = round_log2(len);
        model->mlog_rle->EncodeSymbol(log_length);

        // std::cerr << "RLE=" << (ref&1) << ",len=" << len << ",log=" << log_length << std::endl;

        if (log_length < 2) {
            std::cerr << "single=" << len << "," << (ref&1) << std::endl;
        }
        else if (log_length <= 8) {
            assert(len < 256);

            model->mrle->model_context <<= 1;
            model->mrle->model_context |= (ref & 1);
            model->mrle->model_context <<= 4;
            model->mrle->model_context |= log_length;
            model->mrle->model_context &= model->mrle->model_ctx_mask;
            // std::cerr << "inserting=" << (len&255) << std::endl;
            model->mrle->EncodeSymbolNoUpdate(len & 255);

        } else if (log_length <= 16) {
            assert(len < 65536);

            model->mrle2_1->model_context <<= 1;
            model->mrle2_1->model_context |= (ref & 1);
            model->mrle2_1->model_context <<= 4;
            model->mrle2_1->model_context |= log_length;
            model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
            model->mrle2_1->EncodeSymbolNoUpdate(len & 255);
            // std::cerr << "inserting1=" << (len&255) << std::endl;
            len >>= 8;
            model->mrle2_2->model_context <<= 1;
            model->mrle2_2->model_context |= (ref & 1);
            model->mrle2_2->model_context <<= 4;
            model->mrle2_2->model_context |= log_length;
            model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
            model->mrle2_2->EncodeSymbolNoUpdate(len & 255);
            // std::cerr << "inserting2=" << (len&255) << std::endl;
        } else {
            model->mrle4_1->model_context <<= 1;
            model->mrle4_1->model_context |= (ref & 1);
            model->mrle4_1->model_context <<= 4;
            model->mrle4_1->model_context |= log_length;
            model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
            model->mrle4_1->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            model->mrle4_2->model_context <<= 1;
            model->mrle4_2->model_context |= (ref & 1);
            model->mrle4_2->model_context <<= 4;
            model->mrle4_2->model_context |= log_length;
            model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
            model->mrle4_2->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            model->mrle4_3->model_context <<= 1;
            model->mrle4_3->model_context |= (ref & 1);
            model->mrle4_3->model_context <<= 4;
            model->mrle4_3->model_context |= log_length;
            model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
            model->mrle4_3->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            model->mrle4_4->model_context <<= 1;
            model->mrle4_4->model_context |= (ref & 1);
            model->mrle4_4->model_context <<= 4;
            model->mrle4_4->model_context |= log_length;
            model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
            model->mrle4_4->EncodeSymbolNoUpdate(len & 255);
        }

        return 1;
    }

    int EncodeWahRLE_nm(uint64_t ref, uint32_t len, djinn_ctx_model_t* model) {
        model->mref->EncodeSymbol(ref&15);
        uint32_t log_length = round_log2(len);
        model->mlog_rle->EncodeSymbol(log_length);

        // std::cerr << "RLE=" << (ref&1) << ",len=" << len << ",log=" << log_length << std::endl;

        if (log_length < 2) {
            std::cerr << "single=" << len << "," << (ref&1) << std::endl;
        }
        else if (log_length <= 8) {
            assert(len < 256);

            model->mrle->model_context = (ref & 15);
            model->mrle->model_context <<= 4;
            model->mrle->model_context |= log_length;
            model->mrle->model_context &= model->mrle->model_ctx_mask;
            // std::cerr << "inserting=" << (len&255) << std::endl;
            model->mrle->EncodeSymbolNoUpdate(len & 255);

        } else if (log_length <= 16) {
            assert(len < 65536);

            model->mrle2_1->model_context = (ref & 15);
            model->mrle2_1->model_context <<= 4;
            model->mrle2_1->model_context |= log_length;
            model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
            model->mrle2_1->EncodeSymbolNoUpdate(len & 255);
            // std::cerr << "inserting1=" << (len&255) << std::endl;
            len >>= 8;
            model->mrle2_2->model_context = (ref & 15);
            model->mrle2_2->model_context <<= 4;
            model->mrle2_2->model_context |= log_length;
            model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
            model->mrle2_2->EncodeSymbolNoUpdate(len & 255);
            // std::cerr << "inserting2=" << (len&255) << std::endl;
        } else {
            model->mrle4_1->model_context = (ref & 15);
            model->mrle4_1->model_context <<= 4;
            model->mrle4_1->model_context |= log_length;
            model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
            model->mrle4_1->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            model->mrle4_2->model_context = (ref & 15);
            model->mrle4_2->model_context <<= 4;
            model->mrle4_2->model_context |= log_length;
            model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
            model->mrle4_2->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            model->mrle4_3->model_context = (ref & 15);
            model->mrle4_3->model_context <<= 4;
            model->mrle4_3->model_context |= log_length;
            model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
            model->mrle4_3->EncodeSymbolNoUpdate(len & 255);
            len >>= 8;
            model->mrle4_4->model_context = (ref & 15);
            model->mrle4_4->model_context <<= 4;
            model->mrle4_4->model_context |= log_length;
            model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
            model->mrle4_4->EncodeSymbolNoUpdate(len & 255);
        }

        return 1;
    }

    int DecodeWahRLE(uint64_t& ref, uint32_t& len) {
        ref = model_2mc.mref->DecodeSymbol();
        uint32_t log_length = model_2mc.mlog_rle->DecodeSymbol();
        // std::cerr << "ref=" << ref << " log=" << log_length << std::endl;

        if (log_length < 2) {
            std::cerr << "single=" << len << "," << (ref&1) << std::endl;
        }
        else if (log_length <= 8) {
            model_2mc.mrle->model_context <<= 1;
            model_2mc.mrle->model_context |= (ref & 1);
            model_2mc.mrle->model_context <<= 4;
            model_2mc.mrle->model_context |= log_length;
            model_2mc.mrle->model_context &= model_2mc.mrle->model_ctx_mask;
            len = model_2mc.mrle->DecodeSymbolNoUpdate();
        } else if (log_length <= 16) {
            model_2mc.mrle2_1->model_context <<= 1;
            model_2mc.mrle2_1->model_context |= (ref & 1);
            model_2mc.mrle2_1->model_context <<= 4;
            model_2mc.mrle2_1->model_context |= log_length;
            model_2mc.mrle2_1->model_context &= model_2mc.mrle2_1->model_ctx_mask;
            len = model_2mc.mrle2_1->DecodeSymbolNoUpdate();
            model_2mc.mrle2_2->model_context <<= 1;
            model_2mc.mrle2_2->model_context |= (ref & 1);
            model_2mc.mrle2_2->model_context <<= 4;
            model_2mc.mrle2_2->model_context |= log_length;
            model_2mc.mrle2_2->model_context &= model_2mc.mrle2_2->model_ctx_mask;
            len |= (uint32_t)model_2mc.mrle2_2->DecodeSymbolNoUpdate() << 8;
        } else {
            model_2mc.mrle4_1->model_context <<= 1;
            model_2mc.mrle4_1->model_context |= (ref & 1);
            model_2mc.mrle4_1->model_context <<= 4;
            model_2mc.mrle4_1->model_context |= log_length;
            model_2mc.mrle4_1->model_context &= model_2mc.mrle4_1->model_ctx_mask;
            len = model_2mc.mrle4_1->DecodeSymbolNoUpdate();
            model_2mc.mrle4_2->model_context <<= 1;
            model_2mc.mrle4_2->model_context |= (ref & 1);
            model_2mc.mrle4_2->model_context <<= 4;
            model_2mc.mrle4_2->model_context |= log_length;
            model_2mc.mrle4_2->model_context &= model_2mc.mrle4_2->model_ctx_mask;
            len |= (uint32_t)model_2mc.mrle4_2->DecodeSymbolNoUpdate() << 8;
            model_2mc.mrle4_3->model_context <<= 1;
            model_2mc.mrle4_3->model_context |= (ref & 1);
            model_2mc.mrle4_3->model_context <<= 4;
            model_2mc.mrle4_3->model_context |= log_length;
            model_2mc.mrle4_3->model_context &= model_2mc.mrle4_3->model_ctx_mask;
            len |= (uint32_t)model_2mc.mrle4_3->DecodeSymbolNoUpdate() << 16;
            model_2mc.mrle4_4->model_context <<= 1;
            model_2mc.mrle4_4->model_context |= (ref & 1);
            model_2mc.mrle4_4->model_context <<= 4;
            model_2mc.mrle4_4->model_context |= log_length;
            model_2mc.mrle4_4->model_context &= model_2mc.mrle4_4->model_ctx_mask;
            len |= (uint32_t)model_2mc.mrle4_4->DecodeSymbolNoUpdate() << 24;
        }

        return 1;
    }

    int Compress(uint8_t* out); // compress ctx_block to buffer
    int Compress(djinn_ctx_block_t*& block); // compress ctx_block to given block

    void StartEncoding(bool use_pbwt, bool reset = false) {
        if (use_pbwt) {
            if (model_2mc.pbwt.n_symbols == 0) 
                model_2mc.pbwt.Initiate(n_samples, 2);

            if (model_nm.pbwt.n_symbols == 0) 
                model_nm.pbwt.Initiate(n_samples, 16);
        }

        // If resetting the model_2mc.
        if (reset) {
            model_2mc.reset();
            model_nm.reset();
        }

        // Local range coder
        range_coder->SetOutput(p);
        range_coder->StartEncode();

        model_2mc.StartEncoding(use_pbwt, reset);
        model_nm.StartEncoding(use_pbwt, reset);
    }

    size_t FinishEncoding() {
        if (range_coder.get() == nullptr) return -1;
        range_coder->FinishEncode();
        model_2mc.FinishEncoding();
        model_nm.FinishEncoding();
        std::cerr << range_coder->OutSize() << " and " << model_2mc.range_coder->OutSize() << " and " << model_nm.range_coder->OutSize() << std::endl;
        return range_coder->OutSize() + model_2mc.range_coder->OutSize() + model_nm.range_coder->OutSize();
    }

    void StartDecoding(uint8_t* data, bool reset = false) {
        if (p_free) delete[] p;
        p = data; p_free = false;
        range_coder->SetInput(p);
        range_coder->StartDecode();
        model_2mc.StartDecoding(p, reset);
        model_nm.StartDecoding(p, reset);
    }

public:
    int64_t n_samples;
    
    int64_t n_samples_wah;
    uint64_t n_wah;
    uint64_t* wah_bitmaps;

    uint8_t *p;     // data
    uint32_t p_len; // data length
    uint32_t p_cap:31, p_free:1;
    uint32_t n_variants;
    
    std::shared_ptr<RangeCoder> range_coder;
    std::shared_ptr<GeneralModel> marchetype; // 0 for 2MC, 1 for 2M, 2 else
    djinn_ctx_model_t model_2mc;
    djinn_ctx_model_t model_2m;
    djinn_ctx_model_t model_nm;
};

}