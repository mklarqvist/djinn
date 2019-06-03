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

}
