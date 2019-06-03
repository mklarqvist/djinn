#include "ctx_model.h"

namespace djinn {

djinn_ctx_model_t::djinn_ctx_model_t() :
    use_pbwt(true), init(false), unused(0),
    p(nullptr),
    p_len(0), p_cap(0), p_free(false),
    n_variants(0)
{
    
}

djinn_ctx_model_t::~djinn_ctx_model_t() {
    if (p_free) delete[] p;
}

void djinn_ctx_model_t::Initiate2mc() {
    range_coder = std::make_shared<RangeCoder>();
    mref = std::make_shared<GeneralModel>(2, 1024, 24, 32, range_coder);
    mlog_rle = std::make_shared<GeneralModel>(32, 4096, 24, 16, range_coder);
    mrle = std::make_shared<GeneralModel>(256, 256, 24, 32, range_coder);
    mrle2_1 = std::make_shared<GeneralModel>(256, 256, 24, 8, range_coder);
    mrle2_2 = std::make_shared<GeneralModel>(256, 256, 24, 8, range_coder);
    mrle4_1 = std::make_shared<GeneralModel>(256, 256, 24, 8, range_coder);
    mrle4_2 = std::make_shared<GeneralModel>(256, 256, 24, 8, range_coder);
    mrle4_3 = std::make_shared<GeneralModel>(256, 256, 24, 8, range_coder);
    mrle4_4 = std::make_shared<GeneralModel>(256, 256, 24, 8, range_coder);
    dirty_wah = std::make_shared<GeneralModel>(256, 256, 20, 32, range_coder);
    mtype = std::make_shared<GeneralModel>(2, 1024, 24, 1, range_coder);

    // for (int i = 0; i < 1500; ++i) {
    //     // for (int j = 0; j < 256; ++j) {
    //         dirty_wah->models[i]->total_frequency += 256;
    //         dirty_wah->models[i]->F[256].Freq += 256;
    //     // }
    // }
    // for (int i = 62000; i < 65536; ++i) {
    //     dirty_wah->models[i]->total_frequency += 32;
    //     dirty_wah->models[i]->F[255].Freq += 32;
    // }

    for (int i = 0; i < 256; ++i) {
            dirty_wah->models[i]->total_frequency += 32;
            dirty_wah->models[i]->F[1].Freq += 32;
            
        // dirty_wah->models[i]->total_frequency += 10000;
        // dirty_wah->models[i]->F[1].Freq += 10000;
    }
    
    // range_coder = std::make_shared<RangeCoder>();
    // mref = std::make_shared<GeneralModel>(2, 4096, range_coder);
    // mlog_rle = std::make_shared<GeneralModel>(32, 4096, range_coder);
    // mrle = std::make_shared<GeneralModel>(256, 1024, range_coder);
    // mrle2_1 = std::make_shared<GeneralModel>(256, 1024, range_coder);
    // mrle2_2 = std::make_shared<GeneralModel>(256, 1024, range_coder);
    // mrle4_1 = std::make_shared<GeneralModel>(256, 1024, range_coder);
    // mrle4_2 = std::make_shared<GeneralModel>(256, 1024, range_coder);
    // mrle4_3 = std::make_shared<GeneralModel>(256, 1024, range_coder);
    // mrle4_4 = std::make_shared<GeneralModel>(256, 1024, range_coder);
    // dirty_wah = std::make_shared<GeneralModel>(256, 65536, range_coder);
    // mtype = std::make_shared<GeneralModel>(2, 1024, range_coder);
}

void djinn_ctx_model_t::InitiateNm() {
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

void djinn_ctx_model_t::reset() {
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
    pbwt.Reset();
    n_variants = 0; p_len = 0;
}

int djinn_ctx_model_t::StartEncoding(bool use_pbwt, bool reset) {
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

size_t djinn_ctx_model_t::FinishEncoding() {
    if (range_coder.get() == nullptr) return -1;
    range_coder->FinishEncode();
    return range_coder->OutSize();
}

int djinn_ctx_model_t::StartDecoding(uint8_t* data, bool reset) {
    // If resetting the model_2mc.
    if (range_coder.get() == nullptr) return -1;
    if (reset) this->reset();
    if (data == nullptr) return 0; // or result in corruption as range coder immediately loads data

    if (p_free) delete[] p;
    p = data; p_free = false;

    range_coder->SetInput(p);
    range_coder->StartDecode();
    return 1;
}

///////////////////////////

djinn_ctx_model::djinn_ctx_model() : 
    n_samples(0), 
    n_wah(0), n_samples_wah(0), wah_bitmaps(nullptr), 
    p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
    range_coder(std::make_shared<RangeCoder>()), 
    marchetype(std::make_shared<GeneralModel>(2, 1024, range_coder))
{
    model_2mc.Initiate2mc();
    model_nm.InitiateNm();
}

djinn_ctx_model::djinn_ctx_model(uint64_t n_s) : 
    n_samples(n_s), 
    n_wah(std::ceil((float)n_samples*4 / 64)), 
    n_samples_wah((n_wah*64)/4), 
    wah_bitmaps(new uint64_t[n_wah]), 
    p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
    range_coder(std::make_shared<RangeCoder>()), 
    marchetype(std::make_shared<GeneralModel>(2, 1024, range_coder))
{
    model_2mc.Initiate2mc();
    model_nm.InitiateNm();
}

djinn_ctx_model::~djinn_ctx_model() { 
    delete[] wah_bitmaps;
    if (p_free) delete[] p;
}

void djinn_ctx_model::SetSamples(int64_t n_s) { 
    n_samples = n_s;
    n_wah = std::ceil((float)n_samples / 64) * 4;
    n_samples_wah = (n_wah * 64) / 4;
} 

int djinn_ctx_model::EncodeBcf(uint8_t* data, uint8_t alt_alleles) {
    if (n_samples == 0) return -1;
    if (data == nullptr) return -2;

    // Todo: check that input data is equal to length
    // otherwise if length/2 (haploid)
    // otherwise store but no PBWT regardless of desired or not

    memset(hist_alts, 0, 256*sizeof(uint32_t));
    for (int i = 0; i < n_samples; ++i) {
        ++hist_alts[BCF_UNPACK_GENOTYPE_GENERAL(data[i])];
    }

    // for (int i = 0; i < 256; ++i) {
    //     if (hist_alts[i]) std::cerr << i << ":" << hist_alts[i] << ",";
    // }
    // std::cerr << " data2m=" << model_2mc.range_coder->OutSize() << " dataNm=" << model_nm.range_coder->OutSize() << std::endl;

    // Biallelic, no missing, and no special EOV symbols.
    if (alt_alleles <= 2 && hist_alts[14] == 0 && hist_alts[15] == 0) {
        marchetype->EncodeSymbol(0); // add archtype as 2mc

        if (hist_alts[2] < 10) { // dont update if < 10 alts
            for (int i = 0; i < n_samples; ++i) {
                model_2mc.pbwt.prev[i] = BCF_UNPACK_GENOTYPE(data[model_2mc.pbwt.ppa[i]]);
                // model_2mc.pbwt.prev[i] = BCF_UNPACK_GENOTYPE(data[i]);
            }
        } else {
            model_2mc.pbwt.UpdateBcf(data, 1);
        }

        return Encode2mc(model_2mc.pbwt.prev, n_samples);
    } else {
        marchetype->EncodeSymbol(1);
        model_nm.pbwt.UpdateBcfGeneral(data, 1); // otherwise
        return EncodeNm(model_nm.pbwt.prev, n_samples);
    }
}

int djinn_ctx_model::Encode2mc(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;
    
    if (wah_bitmaps == nullptr) {
        n_wah = std::ceil((float)n_samples / 64) * 4;
        n_samples_wah = (n_wah * 64) / 4;
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

int djinn_ctx_model::EncodeNm(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;
    
    if (wah_bitmaps == nullptr) {
        n_wah = std::ceil((float)n_samples / 64) * 4;
        n_samples_wah = (n_wah * 64) / 4;
        wah_bitmaps = new uint64_t[n_wah];
    }
    
    memset(wah_bitmaps, 0, n_wah*sizeof(uint64_t));

    for (int i = 0; i < n_samples; ++i) {
        wah_bitmaps[i / 16] |= (uint64_t)data[i] << (4*(i % 16));
    }

    // for (int i = 0; i < n_wah; ++i) std::cerr << std::bitset<64>(wah_bitmaps[i]) << " ";
    // std::cerr << std::endl;

    return EncodeWahNm(wah_bitmaps, n_wah);
    // return -1;
}

int djinn_ctx_model::DecodeNext(uint8_t* data) {
    if (data == nullptr) return -1;
    uint8_t type = marchetype->DecodeSymbol();
    
    switch(type) {
    case 0: return DecodeRaw(data);
    case 1: return DecodeRaw_nm(data);
    default: std::cerr << "decoding error: " << type << std::endl; return -1;
    }

    return -1;
}

// Return raw, potentially permuted, WAH-like encoding + archetype dict
int djinn_ctx_model::DecodeRaw(uint8_t* data) {
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
            DecodeWahRLE(ref, len, &model_2mc);
            // std::cerr << "decoded RLE=" << ref << " len=" << len << std::endl;
            n_samples_obs += len*64;
        }

        if (n_samples_obs == n_samples_wah) {
            // std::cerr << "obs=" << n_samples_obs << "/" << n_samples_wah << std::endl;
            break;
        }
        if (n_samples_obs > n_samples_wah) {
            std::cerr << "Decompression corruption: " << n_samples_obs << "/" << n_samples_wah  << std::endl;
            exit(1);
        }
    }

    return 1;
}

int djinn_ctx_model::DecodeRaw_nm(uint8_t* data) {
    if (data == nullptr) return -1;

    int64_t n_samples_obs = 0;
    while(true) {
        uint8_t type = model_nm.mtype->DecodeSymbol();
        // std::cerr << "Type=" << (int)type << std::endl;
        if (type == 0) { // bitmaps
            uint64_t wah = 0;
            // std::cerr << "Bitmap=";
            for (int i = 0; i < 8; ++i) {
                wah |= model_nm.dirty_wah->DecodeSymbol();
                // std::cerr << std::bitset<8>(wah&255);
                wah <<= 8;
            }
            // std::cerr << std::endl;
            n_samples_obs += 16;

        } else { // is RLE
            // Decode an RLE
            uint64_t ref = 1; uint32_t len = 1;
            DecodeWahRLE_nm(ref, len, &model_nm);
            // std::cerr << "decoded RLE=" << ref << " len=" << len << std::endl;
            n_samples_obs += len*16;
        }

        if (n_samples_obs == n_samples_wah) {
            // std::cerr << "obs=" << n_samples_obs << "/" << n_samples_wah << std::endl;
            break;
        }
        if (n_samples_obs > n_samples_wah) {
            std::cerr << "Decompression corruption: " << n_samples_obs << "/" << n_samples_wah  << std::endl;
            exit(1);
        }
    }

    return 1;
}

int djinn_ctx_model::EncodeWah(uint64_t* wah, uint32_t len) { // input WAH-encoded data
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

    ++n_variants;
    return 1;
}

int djinn_ctx_model::EncodeWahNm(uint64_t* wah, uint32_t len) { // input WAH-encoded data
    if (wah == nullptr) return -1;

    ++model_nm.n_variants;

    // std::cerr << "refnm=" << std::bitset<64>(wah[0]) << std::endl;

    uint64_t wah_ref = wah[0];
    uint64_t wah_run = 1;

    for (int i = 1; i < len; ++i) {
        if ((wah_ref != ref_bits[1] && wah_ref != ref_bits[2] && wah_ref != ref_bits[0]) || (wah_ref != wah[i])) {
            if ((wah_ref != wah_bitmaps[i]) || wah_run == 1) {
                // std::cerr << "Bitmap=" << std::bitset<64>(wah_ref) << " = " << wah_ref << std::endl;
                model_nm.mtype->EncodeSymbol(0);
                for (int i = 0; i < 8; ++i) {
                    model_nm.dirty_wah->EncodeSymbol(wah_ref & 255);
                    wah_ref >>= 8;
                }
            } else {
                // std::cerr << "RLE=" << wah_ref << "(" << (wah_ref&15) << "):" << wah_run << std::endl;
                model_nm.mtype->EncodeSymbol(1);
                EncodeWahRLE_nm(wah_ref, wah_run, &model_nm);
            }
            wah_run = 0;
            wah_ref = wah[i];
        }
        ++wah_run;
    }

    if (wah_run) {
        if ((wah_ref != ref_bits[1] && wah_ref != ref_bits[2] && wah_ref != ref_bits[0]) || wah_run == 1) {
            // std::cerr << "Bitmap=" << std::bitset<64>(wah_ref) << " = " << wah_ref << std::endl;
            model_nm.mtype->EncodeSymbol(0);
            
            for (int i = 0; i < 8; ++i) {
                model_nm.dirty_wah->EncodeSymbol(wah_ref & 255);
                wah_ref >>= 8;
            }
            
        } else {
            // std::cerr << "RLE=" << wah_ref << "(" << (wah_ref&15) << "):" << wah_run << std::endl;
            model_nm.mtype->EncodeSymbol(1);
            EncodeWahRLE_nm(wah_ref, wah_run, &model_nm);
        }
    }

    ++n_variants;
    return 1;
}

int djinn_ctx_model::EncodeWahRLE(uint64_t ref, uint32_t len, djinn_ctx_model_t* model) {
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
        model->mrle->EncodeSymbolNoUpdate(len & 255);

    } else if (log_length <= 16) {
        assert(len < 65536);
        model->mrle2_1->model_context <<= 1;
        model->mrle2_1->model_context |= (ref & 1);
        model->mrle2_1->model_context <<= 4;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        model->mrle2_1->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle2_2->model_context <<= 1;
        model->mrle2_2->model_context |= (ref & 1);
        model->mrle2_2->model_context <<= 4;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        model->mrle2_2->EncodeSymbolNoUpdate(len & 255);
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

int djinn_ctx_model::EncodeWahRLE_nm(uint64_t ref, uint32_t len, djinn_ctx_model_t* model) {
    model->mref->EncodeSymbol(ref & 15);
    uint32_t log_length = round_log2(len);
    model->mlog_rle->EncodeSymbol(log_length);

    // std::cerr << "nmRLE=" << (ref&15) << ",len=" << len << ",log=" << log_length << std::endl;

    if (log_length < 2) {
        std::cerr << "single=" << len << "," << (ref & 15) << std::endl;
    }
    else if (log_length <= 8) {
        assert(len < 256);
        model->mrle->model_context = (ref & 15);
        model->mrle->model_context <<= 4;
        model->mrle->model_context |= log_length;
        model->mrle->model_context &= model->mrle->model_ctx_mask;
        model->mrle->EncodeSymbolNoUpdate(len & 255);

    } else if (log_length <= 16) {
        assert(len < 65536);
        model->mrle2_1->model_context = (ref & 15);
        model->mrle2_1->model_context <<= 4;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        model->mrle2_1->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle2_2->model_context = (ref & 15);
        model->mrle2_2->model_context <<= 4;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        model->mrle2_2->EncodeSymbolNoUpdate(len & 255);
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

int djinn_ctx_model::DecodeWahRLE(uint64_t& ref, uint32_t& len, djinn_ctx_model_t* model) {
    ref = model->mref->DecodeSymbol();
    uint32_t log_length = model->mlog_rle->DecodeSymbol();
    // std::cerr << "ref=" << ref << " log=" << log_length << std::endl;

    if (log_length < 2) {
        std::cerr << "single=" << len << "," << (ref&1) << std::endl;
    }
    else if (log_length <= 8) {
        model->mrle->model_context <<= 1;
        model->mrle->model_context |= (ref & 1);
        model->mrle->model_context <<= 4;
        model->mrle->model_context |= log_length;
        model->mrle->model_context &= model->mrle->model_ctx_mask;
        len = model->mrle->DecodeSymbolNoUpdate();
    } else if (log_length <= 16) {
        model->mrle2_1->model_context <<= 1;
        model->mrle2_1->model_context |= (ref & 1);
        model->mrle2_1->model_context <<= 4;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        len = model->mrle2_1->DecodeSymbolNoUpdate();
        model->mrle2_2->model_context <<= 1;
        model->mrle2_2->model_context |= (ref & 1);
        model->mrle2_2->model_context <<= 4;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        len |= (uint32_t)model->mrle2_2->DecodeSymbolNoUpdate() << 8;
    } else {
        model->mrle4_1->model_context <<= 1;
        model->mrle4_1->model_context |= (ref & 1);
        model->mrle4_1->model_context <<= 4;
        model->mrle4_1->model_context |= log_length;
        model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
        len = model->mrle4_1->DecodeSymbolNoUpdate();
        model->mrle4_2->model_context <<= 1;
        model->mrle4_2->model_context |= (ref & 1);
        model->mrle4_2->model_context <<= 4;
        model->mrle4_2->model_context |= log_length;
        model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
        len |= (uint32_t)model->mrle4_2->DecodeSymbolNoUpdate() << 8;
        model->mrle4_3->model_context <<= 1;
        model->mrle4_3->model_context |= (ref & 1);
        model->mrle4_3->model_context <<= 4;
        model->mrle4_3->model_context |= log_length;
        model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
        len |= (uint32_t)model->mrle4_3->DecodeSymbolNoUpdate() << 16;
        model->mrle4_4->model_context <<= 1;
        model->mrle4_4->model_context |= (ref & 1);
        model->mrle4_4->model_context <<= 4;
        model->mrle4_4->model_context |= log_length;
        model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
        len |= (uint32_t)model->mrle4_4->DecodeSymbolNoUpdate() << 24;
    }

    return 1;
}

int djinn_ctx_model::DecodeWahRLE_nm(uint64_t& ref, uint32_t& len, djinn_ctx_model_t* model) {
    ref = model->mref->DecodeSymbol();
    uint32_t log_length = model->mlog_rle->DecodeSymbol();
    // std::cerr << "ref=" << ref << " log=" << log_length << std::endl;

    if (log_length < 2) {
        std::cerr << "single=" << len << "," << (ref&15) << std::endl;
    }
    else if (log_length <= 8) {
        model->mrle->model_context = (ref & 15);
        model->mrle->model_context <<= 4;
        model->mrle->model_context |= log_length;
        model->mrle->model_context &= model->mrle->model_ctx_mask;
        len = model->mrle->DecodeSymbolNoUpdate();
    } else if (log_length <= 16) {
        model->mrle2_1->model_context = (ref & 15);
        model->mrle2_1->model_context <<= 4;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        len = model->mrle2_1->DecodeSymbolNoUpdate();
        model->mrle2_2->model_context = (ref & 15);
        model->mrle2_2->model_context <<= 4;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        len |= (uint32_t)model->mrle2_2->DecodeSymbolNoUpdate() << 8;
    } else {
        model->mrle4_1->model_context = (ref & 15);
        model->mrle4_1->model_context <<= 4;
        model->mrle4_1->model_context |= log_length;
        model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
        len = model->mrle4_1->DecodeSymbolNoUpdate();
        model->mrle4_2->model_context = (ref & 15);
        model->mrle4_2->model_context <<= 4;
        model->mrle4_2->model_context |= log_length;
        model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
        len |= (uint32_t)model->mrle4_2->DecodeSymbolNoUpdate() << 8;
        model->mrle4_3->model_context = (ref & 15);
        model->mrle4_3->model_context <<= 4;
        model->mrle4_3->model_context |= log_length;
        model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
        len |= (uint32_t)model->mrle4_3->DecodeSymbolNoUpdate() << 16;
        model->mrle4_4->model_context = (ref & 15);
        model->mrle4_4->model_context <<= 4;
        model->mrle4_4->model_context |= log_length;
        model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
        len |= (uint32_t)model->mrle4_4->DecodeSymbolNoUpdate() << 24;
    }

    return 1;
}

void djinn_ctx_model::StartEncoding(bool use_pbwt, bool reset) {
    if (use_pbwt) {
        if (model_2mc.pbwt.n_symbols == 0) 
            model_2mc.pbwt.Initiate(n_samples, 2);

        if (model_nm.pbwt.n_symbols == 0) 
            model_nm.pbwt.Initiate(n_samples, 16);
    }

    // If resetting the model_2mc.
    if (reset) {
        std::cerr << "resetting" << std::endl;
        model_2mc.reset();
        model_nm.reset();
    } else {
        model_2mc.n_variants = 0;
        model_nm.n_variants = 0;
    }

    n_variants = 0;

    // Local range coder
    range_coder->SetOutput(p);
    range_coder->StartEncode();

    model_2mc.StartEncoding(use_pbwt, reset);
    model_nm.StartEncoding(use_pbwt, reset);
}

size_t djinn_ctx_model::FinishEncoding() {
    if (range_coder.get() == nullptr) return -1;
    range_coder->FinishEncode();
    model_2mc.FinishEncoding();
    model_nm.FinishEncoding();

    size_t s_rc  = range_coder->OutSize();
    size_t s_2mc = model_2mc.range_coder->OutSize();
    size_t s_nm  = model_nm.range_coder->OutSize();
    std::cerr << s_rc << " and " << s_2mc << " and " << s_nm << std::endl;
    
    return s_rc + s_2mc + s_nm;
}

int djinn_ctx_model::StartDecoding(djinn_block_t* block, bool reset) {
    if (block == nullptr) return  -1;
    if (block->type != djinn_block_t::BlockType::CONTEXT) return -2;
    
    djinn_ctx_t* data_out = (djinn_ctx_t*)block->data;

    if (p_free) delete[] p;
    p = data_out->ctx_models[0].vptr; p_free = false;
    p_len = data_out->ctx_models[0].vptr_len;

    n_variants = block->n_rcds;
    range_coder->SetInput(p);
    range_coder->StartDecode();
    model_2mc.StartDecoding(data_out->ctx_models[1].vptr, reset);
    model_2mc.n_variants = data_out->ctx_models[1].n_v;
    model_nm.StartDecoding(data_out->ctx_models[3].vptr, reset);
    model_nm.n_variants = data_out->ctx_models[3].n_v;
    return 1;
}

int djinn_ctx_model::GetBlockReference(djinn_block_t*& block) {
    if (block != nullptr) {
        if (block->type != djinn_block_t::BlockType::CONTEXT) {
            delete block;
            block = nullptr;
        } else { // block type is correct context: if data is set then reset it otherwise initiate
            if (block->data != nullptr) {
                block->data->reset();
            } else {
                block->data = new djinn_ctx_t();
            }
        }
    } else { // data is nullptr
        block = new djinn_ctx_block_t;
        block->type = djinn_block_t::BlockType::CONTEXT;
        block->data = new djinn_ctx_t();
    }
    // Internal pointer.
    djinn_ctx_t* data_out = (djinn_ctx_t*)block->data;

    // Archetype
    data_out->ctx_models[0].SetDataReference(p, range_coder->OutSize());
    data_out->ctx_models[0].n = 0;
    data_out->ctx_models[0].n_c = range_coder->OutSize();
    data_out->ctx_models[0].n_v = n_variants;

    djinn_ctx_model_t* target_model[3] = {&model_2mc, &model_2m, &model_nm};
    for (int i = 0; i < 3; ++i) {
        // std::cerr << target_model[i]->n_variants << std::endl;
        if (target_model[i]->n_variants) {
            data_out->ctx_models[i+1].SetDataReference(target_model[i]->p, target_model[i]->range_coder->OutSize());
            data_out->ctx_models[i+1].n = 0;
            data_out->ctx_models[i+1].n_c = target_model[i]->range_coder->OutSize();
            data_out->ctx_models[i+1].n_v = target_model[i]->n_variants;
        }
    }

    block->n_rcds = n_variants;
    block->p_len = data_out->SerializedSize();
    djinn_ctx_ctrl_t* ctrl = (djinn_ctx_ctrl_t*)&block->ctrl;
    ctrl->pbwt = model_2mc.use_pbwt;
    ctrl->dimc = model_2mc.n_variants > 0;
    ctrl->dim  = model_2m.n_variants > 0;
    ctrl->nm   = model_nm.n_variants > 0;

    return 1;
}

}