#include <cassert>
#include <iostream>//debug

#include "frequency_model.h"

namespace djinn {

FrequencyModel::FrequencyModel() :
    n_symbols(0), step_size(1), max_total((1 << 16) - step_size),
    total_frequency(0),
    F(nullptr)
{

}

FrequencyModel::FrequencyModel(int nsym, int max_sym) :
    n_symbols(nsym), step_size(1), max_total((1 << 16) - step_size),
    total_frequency(0),
    F(new SymFreqs[n_symbols+1])
{
    assert(n_symbols >= 1);
    Initiate(nsym, max_sym);
}

FrequencyModel::~FrequencyModel() {
    delete[] F;
}

void FrequencyModel::Initiate(int nsym, int max_sym) {
    assert(nsym >= 1);
    n_symbols = nsym;
    delete[] F;
    F = new SymFreqs[n_symbols+1]; // Add 1 for loop terminator in normalize().
    Initiate(max_sym);
}

void FrequencyModel::Initiate(int max_sym) {
    assert(n_symbols >= 1);
    int i = 0;

    for (/**/; i < max_sym; i++) {
        F[i].Symbol = i;
        F[i].Freq   = 1;
    }

    for (/**/; i < n_symbols; i++) {
        F[i].Symbol = i;
        F[i].Freq   = 0;
    }

    total_frequency = max_sym;
    sentinel.Symbol = 0;
    sentinel.Freq   = max_total; // Always first; simplifies sorting.

    F[n_symbols].Freq = 0; // terminates normalize() loop. See below.
}

void FrequencyModel::Normalize() {
    SymFreqs* s;

    /* Faster than F[i].Freq for 0 <= i < n_symbols */
    total_frequency = 0;
    for (s = F; s->Freq; s++) {
        s->Freq -= s->Freq >> 1;
        total_frequency += s->Freq;
    }
}

void FrequencyModel::EncodeSymbol(RangeCoder* rc, uint16_t sym) {
    SymFreqs* s = F;
    uint32_t AccFreq = 0;

    while (s->Symbol != sym) {
        AccFreq += s++->Freq;
        _mm_prefetch((uint8_t*)(s+1), _MM_HINT_T0);
    }

    // std::cerr << "encoding symbol: " << sym << " with params " << n_symbols << "," << max_total << "," << step_size << " -> current state (" << (void*)this << ") " << AccFreq << "/" << s->Freq << "/" << total_frequency << " P=" << ((double)s->Freq/total_frequency) << std::endl;
    rc->Encode(AccFreq, s->Freq, total_frequency);

    s->Freq += step_size;
    total_frequency += step_size;

    if (total_frequency > max_total)
        Normalize();

    // Keep approx sorted
    if(s != F) { // Prevent s[-1] to segfault when s == F
        if (s[0].Freq > s[-1].Freq) {
            SymFreqs t = s[0];
            s[0]  = s[-1];
            s[-1] = t;
        }
    }
}

void FrequencyModel::EncodeSymbol(uint16_t sym) {
    SymFreqs* s = F;
    uint32_t AccFreq = 0;

    while (s->Symbol != sym) {
        AccFreq += s++->Freq;
        _mm_prefetch((uint8_t*)(s+1), _MM_HINT_T0);
    }

    s->Freq += step_size;
    total_frequency += step_size;

    if (total_frequency > max_total)
        Normalize();

    // Keep approx sorted
    if(s != F) { // Prevent s[-1] to segfault when s == F
        if (s[0].Freq > s[-1].Freq) {
            SymFreqs t = s[0];
            s[0]  = s[-1];
            s[-1] = t;
        }
    }
}

double FrequencyModel::GetP(uint16_t sym) const {
    SymFreqs* s = F;
    uint32_t AccFreq = 0;

    while (s->Symbol != sym) {
        AccFreq += s++->Freq;
        _mm_prefetch((uint8_t*)(s+1), _MM_HINT_T0);
    }
    return s->Freq / (double)total_frequency;
}

uint16_t FrequencyModel::DecodeSymbol(RangeCoder *rc) {
    SymFreqs* s = F;
    uint32_t freq = rc->GetFreq(total_frequency);
    uint32_t AccFreq;

    for (AccFreq = 0; (AccFreq += s->Freq) <= freq; s++)
        _mm_prefetch((uint8_t*)s, _MM_HINT_T0);

    AccFreq -= s->Freq;

    // std::cerr << "decoding symbol with " << n_symbols << "," << SHIFT << "," << STEP << " -> " << AccFreq << "/" << s->Freq << "/" << total_frequency << std::endl;
    rc->Decode(AccFreq, s->Freq, total_frequency);
    s->Freq += step_size;
    total_frequency += step_size;

    if (total_frequency > max_total)
        Normalize();

    // Keep approx sorted
    if(s != F) { // Prevent s[-1] to segfault when s == F
        if (s[0].Freq > s[-1].Freq) {
            SymFreqs t = s[0];
            s[0] = s[-1];
            s[-1] = t;
            return t.Symbol;
        }
    }

    return s->Symbol;
}


/*======   Canonical representation   ======*/

GeneralModel::GeneralModel() noexcept :
    max_model_symbols(0),
    model_context_shift(0),
    model_context(0), model_ctx_mask(0),
    // n_buffer(10000000), buffer(new uint8_t[n_buffer]),
    n_additions(0)
{}

GeneralModel::GeneralModel(int n_symbols, int model_size) :
    max_model_symbols(n_symbols),
    model_context_shift(ceil(log2(n_symbols))),
    model_context(0), model_ctx_mask(model_size - 1),
    range_coder(std::make_shared<RangeCoder>()),
    // n_buffer(10000000),
    // buffer(new uint8_t[n_buffer]),
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
    // delete[] buffer;
    // n_buffer = 10000000;
    // buffer = new uint8_t[n_buffer];
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
    // n_buffer(0),
    // buffer(nullptr),
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
    // delete[] buffer;
    // n_buffer = 0;
    // buffer = nullptr;
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
    // n_buffer(10000000),
    // buffer(new uint8_t[n_buffer]),
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
    // range_coder = std::make_shared<RangeCoder>();
    if (range_coder.get() == nullptr) range_coder = std::make_shared<RangeCoder>();
    // delete[] buffer;
    // n_buffer = 10000000;
    // buffer = new uint8_t[n_buffer];
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
    // n_buffer(0),
    // buffer(nullptr),
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
    // delete[] buffer;
    // n_buffer = 0;
    // buffer = nullptr;
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
    // delete[] buffer;
}

void GeneralModel::Reset() {
    n_additions = 0;
    std::cerr << "[GeneralModel::Reset] Resetting" << std::endl;
    ResetModels();
    std::cerr << "post models" << std::endl;
    ResetContext();
    std::cerr << "post context" << std::endl;
}

void GeneralModel::ResetModels() {
    std::cerr << "[GeneralModel::ResetModels] #models=" << models.size() << std::endl;
    for (int i = 0; i < models.size(); ++i) {
        assert(models[i].get() != nullptr);
        std::cerr << "reset model#" << i << "/" << models.size() << std::endl;
        models[i]->Initiate(max_model_symbols, max_model_symbols);
    }
}

void GeneralModel::ResetContext() { model_context = 0; }

int GeneralModel::FinishEncoding() {
    std::cerr << "should not be used" << std::endl;
    exit(1);
    // range_coder->FinishEncode();
    // int ret = range_coder->OutSize();
    // return(ret);
    return -1;
}

int GeneralModel::FinishDecoding() {
    range_coder->FinishDecode();
    return 1;
}

void GeneralModel::StartEncoding() {
    std::cerr << "should not be used" << std::endl;
    exit(1);
    // range_coder->SetOutput(buffer);
    // range_coder->StartEncode();
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

}
