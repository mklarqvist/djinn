#include <cstring> //memcpy
#include "ctx_model.h"

namespace djinn {

/*======   Canonical representation   ======*/

GeneralModel::GeneralModel() noexcept :
    max_model_symbols(0),
    model_context_shift(0),
    model_context(0), model_ctx_mask(0),
    n_buffer(10000000), buffer(new uint8_t[n_buffer]),
    n_additions(0)
{}

GeneralModel::GeneralModel(int n_symbols, int model_size) :
    max_model_symbols(n_symbols),
    model_context_shift(ceil(log2(n_symbols))),
    model_context(0), model_ctx_mask(model_size - 1),
    range_coder(std::make_shared<RangeCoder>()),
    n_buffer(10000000),
    buffer(new uint8_t[n_buffer]),
    n_additions(0)
{
    assert(n_symbols > 1);
    models.resize(model_size);
    for (int i = 0; i < model_size; ++i) models[i] = std::make_shared<FrequencyModel>();

    Reset();
}

int GeneralModel::Initiate(int n_symbols, int model_size) {
    max_model_symbols = n_symbols;
    model_context_shift = ceil(log2(n_symbols));
    model_context = 0; 
    model_ctx_mask = model_size - 1;
    range_coder = std::make_shared<RangeCoder>();
    delete[] buffer;
    n_buffer = 10000000;
    buffer = new uint8_t[n_buffer];
    n_additions = 0;
    models.clear();

    assert(n_symbols > 1);
    models.resize(model_size);
    for (int i = 0; i < model_size; ++i) models[i] = std::make_shared<FrequencyModel>();

    Reset();
    return 1;
}


GeneralModel::GeneralModel(int n_symbols, int model_size, std::shared_ptr<RangeCoder> rc) :
    max_model_symbols(n_symbols),
    model_context_shift(ceil(log2(n_symbols))),
    model_context(0), model_ctx_mask(model_size - 1),
    range_coder(rc),
    n_buffer(0),
    buffer(nullptr),
    n_additions(0)
{
    assert(n_symbols > 1);
    models.resize(model_size);
    for (int i = 0; i < model_size; ++i) models[i] = std::make_shared<FrequencyModel>();

    Reset();
}

int GeneralModel::Initiate(int n_symbols, int model_size, std::shared_ptr<RangeCoder> rc) {
    max_model_symbols = n_symbols;
    model_context_shift = ceil(log2(n_symbols));
    model_context = 0; 
    model_ctx_mask = model_size - 1;
    range_coder = rc;
    delete[] buffer;
    n_buffer = 0;
    buffer = nullptr;
    n_additions = 0;
    models.clear();

    assert(n_symbols > 1);
    models.resize(model_size);
    for (int i = 0; i < model_size; ++i) models[i] = std::make_shared<FrequencyModel>();

    Reset();
    return 1;
}

GeneralModel::GeneralModel(int n_symbols, int model_size, int shift, int step) :
    max_model_symbols(n_symbols),
    model_context_shift(ceil(log2(n_symbols))),
    model_context(0), model_ctx_mask(model_size - 1),
    range_coder(std::make_shared<RangeCoder>()),
    n_buffer(10000000),
    buffer(new uint8_t[n_buffer]),
    n_additions(0)
{
    assert(n_symbols > 1);
    models.resize(model_size);
    for (int i = 0; i < model_size; ++i) {
        models[i] = std::make_shared<FrequencyModel>();
        models[i]->max_total = (1 << shift) - step;
        models[i]->step_size = step;
    }

    Reset();
}

int GeneralModel::Initiate(int n_symbols, int model_size, int shift, int step) {
    max_model_symbols = n_symbols;
    model_context_shift = ceil(log2(n_symbols));
    model_context = 0; 
    model_ctx_mask = model_size - 1;
    range_coder = std::make_shared<RangeCoder>();
    delete[] buffer;
    n_buffer = 10000000;
    buffer = new uint8_t[n_buffer];
    n_additions = 0;
    models.clear();

    assert(n_symbols > 1);
    models.resize(model_size);
    for (int i = 0; i < model_size; ++i) {
        models[i] = std::make_shared<FrequencyModel>();
        models[i]->max_total = (1 << shift) - step;
        models[i]->step_size = step;
    }

    Reset();
    return 1;
}

GeneralModel::GeneralModel(int n_symbols, int model_size, int shift, int step, std::shared_ptr<RangeCoder> rc) :
    max_model_symbols(n_symbols),
    model_context_shift(ceil(log2(n_symbols))),
    model_context(0), model_ctx_mask(model_size - 1),
    range_coder(rc),
    n_buffer(0),
    buffer(nullptr),
    n_additions(0)
{
    assert(n_symbols > 1);
    models.resize(model_size);
    for (int i = 0; i < model_size; ++i) {
        models[i] = std::make_shared<FrequencyModel>();
        models[i]->max_total = (1 << shift) - step;
        models[i]->step_size = step;
    }

    Reset();
}

int GeneralModel::Initiate(int n_symbols, int model_size, int shift, int step, std::shared_ptr<RangeCoder> rc) {
    max_model_symbols = n_symbols;
    model_context_shift = ceil(log2(n_symbols));
    model_context = 0; 
    model_ctx_mask = model_size - 1;
    range_coder = rc;
    delete[] buffer;
    n_buffer = 0;
    buffer = nullptr;
    n_additions = 0;
    models.clear();

    assert(n_symbols > 1);
    models.resize(model_size);
    for (int i = 0; i < model_size; ++i) {
        models[i] = std::make_shared<FrequencyModel>();
        models[i]->max_total = (1 << shift) - step;
        models[i]->step_size = step;
    }

    Reset();
    return 1;
}

GeneralModel::~GeneralModel() {
    delete[] buffer;
}

void GeneralModel::Reset() {
    n_additions = 0;
    ResetModels();
    ResetContext();
}

void GeneralModel::ResetModels() {
    for (int i = 0; i < models.size(); ++i)
        models[i]->Initiate(max_model_symbols, max_model_symbols);
}

void GeneralModel::ResetContext() { model_context = 0; }

int GeneralModel::FinishEncoding() {
    range_coder->FinishEncode();
    int ret = range_coder->OutSize();
    return(ret);
}

int GeneralModel::FinishDecoding() {
    range_coder->FinishDecode();
    return 1;
}

void GeneralModel::StartEncoding() {
    range_coder->SetOutput(buffer);
    range_coder->StartEncode();
}

void GeneralModel::StartDecoding(uint8_t* data) {
    range_coder->SetInput(data);
    range_coder->StartDecode();
}

void GeneralModel::EncodeSymbol(const uint16_t symbol) {
    models[model_context]->EncodeSymbol(range_coder.get(), symbol);
    model_context <<= model_context_shift;
    model_context |= symbol;
    model_context &= model_ctx_mask;
    _mm_prefetch((const char *)(models[model_context].get()), _MM_HINT_T0);
    ++n_additions;
}

void GeneralModel::EncodeSymbolNoUpdate(const uint16_t symbol) {
    models[model_context]->EncodeSymbol(range_coder.get(), symbol);

    ++n_additions;
}

uint16_t GeneralModel::DecodeSymbol() {
    uint16_t symbol = models[model_context]->DecodeSymbol(range_coder.get());
    model_context <<= model_context_shift;
    model_context |= symbol;
    model_context &= model_ctx_mask;
    assert(model_context < models.size());
    _mm_prefetch((const char *)(models[model_context].get()), _MM_HINT_T0);
    return symbol;
}

uint16_t GeneralModel::DecodeSymbolNoUpdate() {
    uint16_t symbol = models[model_context]->DecodeSymbol(range_coder.get());
    return symbol;
}

/*======   Context container   ======*/

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
    mref = std::make_shared<GeneralModel>(2, 512, 18, 32, range_coder);
    mlog_rle = std::make_shared<GeneralModel>(16, 32, 18, 16, range_coder);  // 2^32 max -> log2() = 5
    mrle = std::make_shared<GeneralModel>(256, 64, 18, 32, range_coder); // 1 bit ref, 5 bit alt
    mrle2_1 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    mrle2_2 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    mrle4_1 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    mrle4_2 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    mrle4_3 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    mrle4_4 = std::make_shared<GeneralModel>(256, 64, 18, 8, range_coder);
    dirty_wah = std::make_shared<GeneralModel>(256, 256, 18, 16, range_coder);
    mtype = std::make_shared<GeneralModel>(2, 512, 18, 1, range_coder);
}

void djinn_ctx_model_t::InitiateNm() {
    range_coder = std::make_shared<RangeCoder>();
    mref = std::make_shared<GeneralModel>(2, 512, 18, 32, range_coder);
    mlog_rle = std::make_shared<GeneralModel>(64, 32, 18, 16, range_coder);  // 4 bit ref, log2(32) = 5 bit -> 2^9 = 512
    mrle = std::make_shared<GeneralModel>(256, 512, 24, 1, range_coder); // 4 bits ref + 5 bits log2(run) -> 2^9
    mrle2_1 = std::make_shared<GeneralModel>(256, 512, range_coder);
    mrle2_2 = std::make_shared<GeneralModel>(256, 512, range_coder);
    mrle4_1 = std::make_shared<GeneralModel>(256, 512, range_coder);
    mrle4_2 = std::make_shared<GeneralModel>(256, 512, range_coder);
    mrle4_3 = std::make_shared<GeneralModel>(256, 512, range_coder);
    mrle4_4 = std::make_shared<GeneralModel>(256, 512, range_coder);
    dirty_wah = std::make_shared<GeneralModel>(256, 256, 18, 32, range_coder);
    mtype = std::make_shared<GeneralModel>(2, 512, 18, 1, range_coder);
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

/*======   Varaint context model   ======*/

djinn_ctx_model::djinn_ctx_model() : 
    n_wah(0), n_samples_wah(0), wah_bitmaps(nullptr), 
    p(new uint8_t[1000000]), p_len(0), p_cap(1000000), p_free(true),
    range_coder(std::make_shared<RangeCoder>()), 
    marchetype(std::make_shared<GeneralModel>(2, 1024, range_coder))
{
    model_2mc.Initiate2mc();
    model_nm.InitiateNm();
}

djinn_ctx_model::djinn_ctx_model(uint64_t n_s) : 
    djinn_model(n_s),
    n_wah(std::ceil((float)n_samples*4 / 32)),
    n_samples_wah((n_wah*32)/4), 
    wah_bitmaps(new uint32_t[n_wah]), 
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
    n_wah = std::ceil((float)n_samples / 32) * 4;
    n_samples_wah = (n_wah * 32) / 4;
} 

int djinn_ctx_model::EncodeBcf(uint8_t* data, uint8_t alt_alleles, const bool permute) {
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
    // std::cerr << " data2m=" << model_2mc.range_coder->OutSize() << " dataNm=" << model_nm.range_coder->OutSize() << " permute=" << permute << std::endl;

    // Biallelic, no missing, and no special EOV symbols.
    if (alt_alleles <= 2 && hist_alts[14] == 0 && hist_alts[15] == 0) {
        marchetype->EncodeSymbol(0); // add archtype as 2mc

        // TODO: we need to remember what the minumum alt cutoff was for the unpermuting
        // process!
        if (hist_alts[2] < 10 || permute == false) { // dont update if < 10 alts
            for (int i = 0; i < n_samples; ++i) {
                model_2mc.pbwt.prev[i] = BCF_UNPACK_GENOTYPE(data[model_2mc.pbwt.ppa[i]]);
            }
        } else {
            model_2mc.pbwt.UpdateBcf(data, 1);
        }

        return Encode2mc(model_2mc.pbwt.prev, n_samples);
    } else {
        marchetype->EncodeSymbol(1);
        
        if (permute)
            model_nm.pbwt.UpdateBcfGeneral(data, 1); // otherwise
        else {
            for (int i = 0; i < n_samples; ++i) {
                model_nm.pbwt.prev[i] = BCF_UNPACK_GENOTYPE_GENERAL(data[model_nm.pbwt.ppa[i]]);
            }
        }
        return EncodeNm(model_nm.pbwt.prev, n_samples);
    }
}

int djinn_ctx_model::Encode2mc(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;
    
    if (wah_bitmaps == nullptr) {
        n_wah = std::ceil((float)n_samples / 32) * 4;
        n_samples_wah = (n_wah * 32) / 4;
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    // uint32_t alts = 0;
    for (int i = 0; i < n_samples; ++i) {
        if (data[i]) {
            wah_bitmaps[i / 32] |= 1L << (i % 32);
            // ++alts;
        }
    }
    // std::cerr << "encoding alts=" << alts << std::endl;

    return EncodeWah(wah_bitmaps, n_wah/4); // 4 bytes in a 32-bit word
}

int djinn_ctx_model::EncodeNm(uint8_t* data, uint32_t len) {
    if (data == nullptr) return -1;
    if (n_samples == 0) return -2;
    
    if (wah_bitmaps == nullptr) {
        n_wah = std::ceil((float)n_samples / 32) * 4;
        n_samples_wah = (n_wah * 32) / 4;
        wah_bitmaps = new uint32_t[n_wah];
    }
    
    memset(wah_bitmaps, 0, n_wah*sizeof(uint32_t));

    for (int i = 0; i < n_samples; ++i) {
        wah_bitmaps[i / 8] |= (uint64_t)data[i] << (4*(i % 8));
    }

    // for (int i = 0; i < n_wah; ++i) std::cerr << std::bitset<64>(wah_bitmaps[i]) << " ";
    // std::cerr << std::endl;

    return EncodeWahNm(wah_bitmaps, n_wah);
    // return -1;
}

int djinn_ctx_model::DecodeNext(uint8_t* data, size_t& len) {
    if (data == nullptr) return -1;
    uint8_t type = marchetype->DecodeSymbol();
    
    // Todo:
    // Step 1: DecodeRaw
    // Step 2: if PBWT -> unpermute
    // Step 3: inflate (either raw, bcf-style, or bitmap encodings)

    int objs = 0;
    switch(type) {
    case 0: objs = DecodeRaw(data, len); break;
    case 1: objs = DecodeRaw_nm(data, len); break;
    default: std::cerr << "decoding error: " << type << std::endl; return -1;
    }

    // Compute alts


    if (objs <= 0) return -1;
    switch(type) {
        case 0: model_2mc.pbwt.ReverseUpdateEWAH(data, len); break;
        case 1: model_nm.pbwt.ReverseUpdateEWAH(data, len); break;
    }

    return objs;
}

int djinn_ctx_model::DecodeNext(uint8_t* ewah_data, size_t& ret_ewah, uint8_t* ret_buffer, size_t& ret_len) {
    if (ewah_data == nullptr) return -1;
    if (ret_buffer == nullptr) return -1;

    // Decode stream archetype.
    uint8_t type = marchetype->DecodeSymbol();

    int objs = 0;
    switch(type) {
    case 0: objs = DecodeRaw(ewah_data, ret_ewah); break;
    case 1: objs = DecodeRaw_nm(ewah_data, ret_ewah); break;
    default: std::cerr << "decoding error: " << type << std::endl; return -1;
    }

    if (objs <= 0) return -1;

    std::cerr << "alts=" << hist_alts[0] << "," << hist_alts[1] << std::endl;

    if (hist_alts[1] >= 10) {
        switch(type) {
            case 0: model_2mc.pbwt.ReverseUpdateEWAH(ewah_data, ret_ewah, ret_buffer); break;
            case 1: model_nm.pbwt.ReverseUpdateEWAH(ewah_data, ret_ewah, ret_buffer); break;
        }
    } else {
        // Todo: finish me
        // Unpack EWAH into literals
        uint32_t local_offset = 0;
        uint32_t ret_pos = 0;
        for (int j = 0; j < objs; ++j) {
            djinn_ewah_t* ewah = (djinn_ewah_t*)&ewah_data[local_offset];
            local_offset += sizeof(djinn_ewah_t);
            // std::cerr << "ewah=" << ewah->ref << "," << ewah->clean << "," << ewah->dirty << std::endl;
            uint32_t to = ret_pos + ewah->clean*32 > n_samples ? n_samples : ret_pos + ewah->clean*32;
            for (int i = ret_pos; i < to; ++i) {
                ret_buffer[i] = (ewah->ref & 1);
            }
            ret_pos = to;

            for (int i = 0; i < ewah->dirty; ++i) {
                to = ret_pos + 32 > n_samples ? n_samples : ret_pos + 32;
                
                uint32_t dirty = *((uint32_t*)(&ewah_data[local_offset])); // copy
                for (int j = ret_pos; j < to; ++j) {
                    ret_buffer[j] = (dirty & 1);
                    dirty >>= 1;
                }
                ++local_offset;
                ret_pos = to;
            }

            // local_offset += ewah->dirty * sizeof(uint32_t);
            // vals += ewah->dirty + ewah->clean;

            // std::cerr << "local=" << local_offset << "/" << offset << std::endl;
            // assert(local_offset <= offset);
        }
        assert(ret_pos == n_samples);
        
    }

    return objs;
}

int djinn_ctx_model::DecodeNextRaw(uint8_t* data, size_t& len) {
    if (data == nullptr) return -1;
    uint8_t type = marchetype->DecodeSymbol();
    
    switch(type) {
    case 0: return DecodeRaw(data, len);
    case 1: return DecodeRaw_nm(data, len);
    default: std::cerr << "decoding error: " << type << std::endl; return -1;
    }

    // Never reached.
    return -1;
}

// Return raw, potentially permuted, EWAH encoding
int djinn_ctx_model::DecodeRaw(uint8_t* data, size_t& len) {
    if (data == nullptr) return -1;

    int64_t n_samples_obs = 0;
    int objects = 0;

    // Emit empty EWAH marker.
    djinn_ewah_t* ewah = (djinn_ewah_t*)&data[len]; 
    ewah->reset();
    len += sizeof(djinn_ewah_t);

    // uint32_t ref_alt[2] = {0};

    // Compute als
    memset(hist_alts, 0, 256*sizeof(uint32_t));

    while(true) {
        uint8_t type = model_2mc.mtype->DecodeSymbol();

        // std::cerr << "Type=" << (int)type << std::endl;
        if (type == 0) { // bitmaps
            ++ewah->dirty;

            // std::cerr << "Bitmap=";
            // std::cerr << "bitmap" << std::endl;
            for (int i = 0; i < 4; ++i) {
                data[len] = model_2mc.dirty_wah->DecodeSymbol();
                hist_alts[1] += __builtin_popcount(data[len]);
                ++len;
                // wah |= model_2mc.dirty_wah->DecodeSymbol();
                // std::cerr << std::bitset<8>(wah&255);
                // wah <<= 8;
            }
            // std::cerr << std::endl;
            n_samples_obs += 32;

        } else { // is RLE
            // Emit new EWAH marker
            if (ewah->clean > 0 || ewah->dirty > 0) {
                ewah = (djinn_ewah_t*)&data[len];
                ewah->reset();
                len += sizeof(djinn_ewah_t);
                ++objects;
            }

            // Decode an RLE
            uint32_t ref = 1; uint32_t len = 1;
            DecodeWahRLE(ref, len, &model_2mc);
            // std::cerr << "decoded RLE=" << ref << " len=" << len << std::endl;
            ewah->ref = ref & 1;
            ewah->clean = len;

            // ref_alt[ewah->ref&1] += ewah->clean*32;
            hist_alts[ewah->ref] += ewah->clean * 32;

            n_samples_obs += len*32;
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

    if (ewah->clean > 0 || ewah->dirty > 0) {
        ++objects;
        // std::cerr << "final=" << std::bitset<64>(*ewah) << " ref=" << ewah_ref << ",run=" << ewah_run << ",dirty=" << ewah_dirty << std::endl;
        // ref_alt[ewah->ref&1] += ewah->clean*32;
        hist_alts[ewah->ref] += ewah->clean * 32;
        // ++len;
        // Todo: dirty
    }

    // std::cerr << "ewah decode ref_alt=" << ref_alt[0] << "," << ref_alt[1] << std::endl;

    return objects;
}

int djinn_ctx_model::DecodeRaw_nm(uint8_t* data, size_t& len) {
    if (data == nullptr) return -1;

    int64_t n_samples_obs = 0;
    int objects = 0;

    // Emit empty EWAH marker.
    djinn_ewah_t* ewah = (djinn_ewah_t*)&data[len]; 
    ewah->reset();
    len += sizeof(djinn_ewah_t);

    while(true) {
        uint8_t type = model_nm.mtype->DecodeSymbol();
        // data[len++] = type; // store archetype

        // std::cerr << "Type=" << (int)type << std::endl;
        if (type == 0) { // bitmaps
            ++ewah->dirty;
            // uint64_t wah = 0;
            // std::cerr << "Bitmap=";
            for (int i = 0; i < 4; ++i) {
                data[len++] = model_nm.dirty_wah->DecodeSymbol();
                // wah |= model_nm.dirty_wah->DecodeSymbol();
                // std::cerr << std::bitset<8>(wah&255);
                // wah <<= 8;
            }
            // std::cerr << std::endl;
            n_samples_obs += 8;

        } else { // is RLE
            // Emit new EWAH marker
            if (ewah->clean > 0 || ewah->dirty > 0) {
                ewah = (djinn_ewah_t*)&data[len];
                ewah->reset();
                len += sizeof(djinn_ewah_t);
                ++objects;
            }

            // Decode an RLE
            uint32_t ref = 1; uint32_t len = 1;
            DecodeWahRLE_nm(ref, len, &model_nm);
            // std::cerr << "decoded RLE=" << ref << " len=" << len << std::endl;
            ewah->ref = ref & 1;
            ewah->clean = len;

            n_samples_obs += len*8;
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

    if (ewah->clean > 0 || ewah->dirty > 0) {
        ++objects;
        // std::cerr << "final=" << std::bitset<64>(*ewah) << " ref=" << ewah_ref << ",run=" << ewah_run << ",dirty=" << ewah_dirty << std::endl;
    }

    return 1;
}

int djinn_ctx_model::EncodeWah(uint32_t* wah, uint32_t len) { // input WAH-encoded data
    if (wah == nullptr) return -1;
    ++model_2mc.n_variants;

    uint32_t wah_ref = wah[0];
    uint32_t wah_run = 1;

    // uint32_t ref_alt[2] = {0};
    // uint32_t wah_ref_pos = 0;

    for (int i = 1; i < len; ++i) {
        if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint32_t>::max()) || (wah_ref != wah[i])) {
            if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint32_t>::max()) || wah_run == 1) {
                model_2mc.mtype->EncodeSymbol(0);
                
                // std::cerr << "Bitmap=" << std::bitset<32>(wah_ref) << ": " << __builtin_popcount(wah_ref) << " at pos=" << wah_ref_pos << std::endl;
                for (int i = 0; i < 4; ++i) {
                    model_2mc.dirty_wah->EncodeSymbol(wah_ref & 255);
                    // std::cerr << "bitmap-" << i << " " << std::bitset<32>(wah_ref * 255)
                    // ref_alt[1] += __builtin_popcount(wah_ref & 255);
                    wah_ref >>= 8;
                }
            } else {
                model_2mc.mtype->EncodeSymbol(1);
                // ref_alt[wah_ref & 1] += wah_run * 32;
                EncodeWahRLE(wah_ref, wah_run, &model_2mc);
            }
            wah_run = 0;
            wah_ref = wah[i];
            // wah_ref_pos = i;
        }
        ++wah_run;
    }

    if (wah_run) {
        // std::cerr << "wah run left=" << wah_run << std::endl;
        if ((wah_ref != 0 && wah_ref != std::numeric_limits<uint32_t>::max()) || wah_run == 1) {
            model_2mc.mtype->EncodeSymbol(0);
            
            // std::cerr << "Bitmap=" << std::bitset<32>(wah_ref) << ": " << __builtin_popcount(wah_ref) << " at pos=" << wah_ref_pos << std::endl;
            for (int i = 0; i < 4; ++i) {
                model_2mc.dirty_wah->EncodeSymbol(wah_ref & 255);
                // ref_alt[1] += __builtin_popcount(wah_ref & 255);
                wah_ref >>= 8;
            }
            
        } else {
            model_2mc.mtype->EncodeSymbol(1);
            // ref_alt[wah_ref & 1] += wah_run * 32;
            EncodeWahRLE(wah_ref, wah_run, &model_2mc);
        }
    }

    // std::cerr << "encoding wah ref_alt=" << ref_alt[0] << "," << ref_alt[1] << std::endl;

    ++n_variants;
    return 1;
}

int djinn_ctx_model::EncodeWahNm(uint32_t* wah, uint32_t len) { // input WAH-encoded data
    if (wah == nullptr) return -1;

    ++model_nm.n_variants;

    // std::cerr << "refnm=" << std::bitset<32>(wah[0]) << std::endl;

    uint32_t wah_ref = wah[0];
    uint64_t wah_run = 1;

    for (int i = 1; i < len; ++i) {
        if ((wah_ref != ref_bits[1] && wah_ref != ref_bits[2] && wah_ref != ref_bits[0]) || (wah_ref != wah[i])) {
            if ((wah_ref != wah_bitmaps[i]) || wah_run == 1) {
                // std::cerr << "Bitmap=" << std::bitset<64>(wah_ref) << " = " << wah_ref << std::endl;
                model_nm.mtype->EncodeSymbol(0);
                for (int i = 0; i < 4; ++i) {
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
            
            for (int i = 0; i < 4; ++i) {
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

int djinn_ctx_model::EncodeWahRLE(uint32_t ref, uint32_t len, djinn_ctx_model_t* model) {
    // Length cannot be 0 so remove 1
    // len -= 1;

    model->mref->EncodeSymbol(ref&1);
    uint32_t log_length = round_log2(len);
    model->mlog_rle->EncodeSymbol(log_length);

    // std::cerr << "RLE=" << (ref&1) << ",len=" << len << ",log=" << log_length << std::endl;

    if (log_length < 2) {
        // std::cerr << "single=" << len << "," << (ref&1) << std::endl;
    }
    else if (log_length <= 8) {
        assert(len < 256);
        model->mrle->model_context <<= 1;
        model->mrle->model_context |= (ref & 1);
        model->mrle->model_context <<= 5;
        model->mrle->model_context |= log_length;
        model->mrle->model_context &= model->mrle->model_ctx_mask;
        model->mrle->EncodeSymbolNoUpdate(len & 255);

    } else if (log_length <= 16) {
        assert(len < 65536);
        model->mrle2_1->model_context <<= 1;
        model->mrle2_1->model_context |= (ref & 1);
        model->mrle2_1->model_context <<= 5;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        model->mrle2_1->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle2_2->model_context <<= 1;
        model->mrle2_2->model_context |= (ref & 1);
        model->mrle2_2->model_context <<= 5;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        model->mrle2_2->EncodeSymbolNoUpdate(len & 255);
    } else {
        model->mrle4_1->model_context <<= 1;
        model->mrle4_1->model_context |= (ref & 1);
        model->mrle4_1->model_context <<= 5;
        model->mrle4_1->model_context |= log_length;
        model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
        model->mrle4_1->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_2->model_context <<= 1;
        model->mrle4_2->model_context |= (ref & 1);
        model->mrle4_2->model_context <<= 5;
        model->mrle4_2->model_context |= log_length;
        model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
        model->mrle4_2->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_3->model_context <<= 1;
        model->mrle4_3->model_context |= (ref & 1);
        model->mrle4_3->model_context <<= 5;
        model->mrle4_3->model_context |= log_length;
        model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
        model->mrle4_3->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_4->model_context <<= 1;
        model->mrle4_4->model_context |= (ref & 1);
        model->mrle4_4->model_context <<= 5;
        model->mrle4_4->model_context |= log_length;
        model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
        model->mrle4_4->EncodeSymbolNoUpdate(len & 255);
    }

    return 1;
}

int djinn_ctx_model::EncodeWahRLE_nm(uint32_t ref, uint32_t len, djinn_ctx_model_t* model) {
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
        model->mrle->model_context <<= 5;
        model->mrle->model_context |= log_length;
        model->mrle->model_context &= model->mrle->model_ctx_mask;
        model->mrle->EncodeSymbolNoUpdate(len & 255);

    } else if (log_length <= 16) {
        assert(len < 65536);
        model->mrle2_1->model_context = (ref & 15);
        model->mrle2_1->model_context <<= 5;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        model->mrle2_1->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle2_2->model_context = (ref & 15);
        model->mrle2_2->model_context <<= 5;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        model->mrle2_2->EncodeSymbolNoUpdate(len & 255);
    } else {
        model->mrle4_1->model_context = (ref & 15);
        model->mrle4_1->model_context <<= 5;
        model->mrle4_1->model_context |= log_length;
        model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
        model->mrle4_1->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_2->model_context = (ref & 15);
        model->mrle4_2->model_context <<= 5;
        model->mrle4_2->model_context |= log_length;
        model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
        model->mrle4_2->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_3->model_context = (ref & 15);
        model->mrle4_3->model_context <<= 5;
        model->mrle4_3->model_context |= log_length;
        model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
        model->mrle4_3->EncodeSymbolNoUpdate(len & 255);
        len >>= 8;
        model->mrle4_4->model_context = (ref & 15);
        model->mrle4_4->model_context <<= 5;
        model->mrle4_4->model_context |= log_length;
        model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
        model->mrle4_4->EncodeSymbolNoUpdate(len & 255);
    }

    return 1;
}

int djinn_ctx_model::DecodeWahRLE(uint32_t& ref, uint32_t& len, djinn_ctx_model_t* model) {
    ref = model->mref->DecodeSymbol();
    uint32_t log_length = model->mlog_rle->DecodeSymbol();
    // std::cerr << "ref=" << ref << " log=" << log_length << std::endl;

    if (log_length < 2) {
        std::cerr << "single=" << len << "," << (ref&1) << std::endl;
    }
    else if (log_length <= 8) {
        model->mrle->model_context <<= 1;
        model->mrle->model_context |= (ref & 1);
        model->mrle->model_context <<= 5;
        model->mrle->model_context |= log_length;
        model->mrle->model_context &= model->mrle->model_ctx_mask;
        len = model->mrle->DecodeSymbolNoUpdate();
    } else if (log_length <= 16) {
        model->mrle2_1->model_context <<= 1;
        model->mrle2_1->model_context |= (ref & 1);
        model->mrle2_1->model_context <<= 5;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        len = model->mrle2_1->DecodeSymbolNoUpdate();
        model->mrle2_2->model_context <<= 1;
        model->mrle2_2->model_context |= (ref & 1);
        model->mrle2_2->model_context <<= 5;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        len |= (uint32_t)model->mrle2_2->DecodeSymbolNoUpdate() << 8;
    } else {
        model->mrle4_1->model_context <<= 1;
        model->mrle4_1->model_context |= (ref & 1);
        model->mrle4_1->model_context <<= 5;
        model->mrle4_1->model_context |= log_length;
        model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
        len = model->mrle4_1->DecodeSymbolNoUpdate();
        model->mrle4_2->model_context <<= 1;
        model->mrle4_2->model_context |= (ref & 1);
        model->mrle4_2->model_context <<= 5;
        model->mrle4_2->model_context |= log_length;
        model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
        len |= (uint32_t)model->mrle4_2->DecodeSymbolNoUpdate() << 8;
        model->mrle4_3->model_context <<= 1;
        model->mrle4_3->model_context |= (ref & 1);
        model->mrle4_3->model_context <<= 5;
        model->mrle4_3->model_context |= log_length;
        model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
        len |= (uint32_t)model->mrle4_3->DecodeSymbolNoUpdate() << 16;
        model->mrle4_4->model_context <<= 1;
        model->mrle4_4->model_context |= (ref & 1);
        model->mrle4_4->model_context <<= 5;
        model->mrle4_4->model_context |= log_length;
        model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
        len |= (uint32_t)model->mrle4_4->DecodeSymbolNoUpdate() << 24;
    }

    return 1;
}

int djinn_ctx_model::DecodeWahRLE_nm(uint32_t& ref, uint32_t& len, djinn_ctx_model_t* model) {
    ref = model->mref->DecodeSymbol();
    uint32_t log_length = model->mlog_rle->DecodeSymbol();
    // std::cerr << "ref=" << ref << " log=" << log_length << std::endl;

    if (log_length < 2) {
        std::cerr << "single=" << len << "," << (ref&15) << std::endl;
    }
    else if (log_length <= 8) {
        model->mrle->model_context = (ref & 15);
        model->mrle->model_context <<= 5;
        model->mrle->model_context |= log_length;
        model->mrle->model_context &= model->mrle->model_ctx_mask;
        len = model->mrle->DecodeSymbolNoUpdate();
    } else if (log_length <= 16) {
        model->mrle2_1->model_context = (ref & 15);
        model->mrle2_1->model_context <<= 5;
        model->mrle2_1->model_context |= log_length;
        model->mrle2_1->model_context &= model->mrle2_1->model_ctx_mask;
        len = model->mrle2_1->DecodeSymbolNoUpdate();
        model->mrle2_2->model_context = (ref & 15);
        model->mrle2_2->model_context <<= 5;
        model->mrle2_2->model_context |= log_length;
        model->mrle2_2->model_context &= model->mrle2_2->model_ctx_mask;
        len |= (uint32_t)model->mrle2_2->DecodeSymbolNoUpdate() << 8;
    } else {
        model->mrle4_1->model_context = (ref & 15);
        model->mrle4_1->model_context <<= 5;
        model->mrle4_1->model_context |= log_length;
        model->mrle4_1->model_context &= model->mrle4_1->model_ctx_mask;
        len = model->mrle4_1->DecodeSymbolNoUpdate();
        model->mrle4_2->model_context = (ref & 15);
        model->mrle4_2->model_context <<= 5;
        model->mrle4_2->model_context |= log_length;
        model->mrle4_2->model_context &= model->mrle4_2->model_ctx_mask;
        len |= (uint32_t)model->mrle4_2->DecodeSymbolNoUpdate() << 8;
        model->mrle4_3->model_context = (ref & 15);
        model->mrle4_3->model_context <<= 5;
        model->mrle4_3->model_context |= log_length;
        model->mrle4_3->model_context &= model->mrle4_3->model_ctx_mask;
        len |= (uint32_t)model->mrle4_3->DecodeSymbolNoUpdate() << 16;
        model->mrle4_4->model_context = (ref & 15);
        model->mrle4_4->model_context <<= 5;
        model->mrle4_4->model_context |= log_length;
        model->mrle4_4->model_context &= model->mrle4_4->model_ctx_mask;
        len |= (uint32_t)model->mrle4_4->DecodeSymbolNoUpdate() << 24;
    }

    return 1;
}

void djinn_ctx_model::StartEncoding(bool use_pbwt, bool reset) {
    if (use_pbwt) {
        if (model_2mc.pbwt.n_symbols == 0) {
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_2mc.pbwt.Initiate(n_samples, 2);
        }

        if (model_nm.pbwt.n_symbols == 0) {
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_nm.pbwt.Initiate(n_samples, 16);
        }
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
    // std::cerr << s_rc << " and " << s_2mc << " and " << s_nm << std::endl;
    
    return s_rc + s_2mc + s_nm;
}

int djinn_ctx_model::StartDecoding(djinn_block_t* block, bool reset) {
    if (block == nullptr) return  -1;
    if (block->type != djinn_block_t::BlockType::CONTEXT) return -2;
    
    djinn_ctx_t* data_out = (djinn_ctx_t*)block->data;

    if (p_free) delete[] p; // if there is data and we own it then delete
    p = data_out->ctx_models[0].vptr; p_free = false; // set data but no ownership
    p_len = data_out->ctx_models[0].vptr_len; // set data length

    n_variants = block->n_rcds;
    range_coder->SetInput(p);
    range_coder->StartDecode();
    model_2mc.StartDecoding(data_out->ctx_models[1].vptr, reset);
    model_2mc.n_variants = data_out->ctx_models[1].n_v;
    model_nm.StartDecoding(data_out->ctx_models[3].vptr, reset);
    model_nm.n_variants = data_out->ctx_models[3].n_v;

    bool use_pbwt = true; // fix this
    if (use_pbwt) {
        if (model_2mc.pbwt.n_symbols == 0) {
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_2mc.pbwt.Initiate(n_samples, 2);
        }

        if (model_nm.pbwt.n_symbols == 0) {
            if (n_samples == 0) {
                std::cerr << "illegal: no sample number set!" << std::endl;
            }
            model_nm.pbwt.Initiate(n_samples, 16);
        }
    }

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