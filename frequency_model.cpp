#include <cassert>
#include <iostream>//debug

#include "frequency_model.h"

namespace djinn {

FrequencyModel::FrequencyModel() :
    NSYM(0), STEP(2), SHIFT(16),
    TotFreq(0),
    F(nullptr)
{
    
}

FrequencyModel::FrequencyModel(int nsym, int max_sym) :
    NSYM(nsym), STEP(2), SHIFT(13),
    TotFreq(0), F(new SymFreqs[NSYM+1])
{
    assert(NSYM >= 1);
    Initiate(max_sym);
}

FrequencyModel::~FrequencyModel() {
    delete[] F;
}

void FrequencyModel::Initiate(int nsym, int max_sym) {
    assert(nsym >= 1);
    NSYM = nsym;
    delete[] F;
    F = new SymFreqs[NSYM+1]; // Add 1 for loop terminator in normalize().
    Initiate(max_sym);
}

void FrequencyModel::Initiate(int max_sym) {
    assert(NSYM >= 1);
    int i = 0;

    for (/**/; i < max_sym; i++) {
        F[i].Symbol = i;
        F[i].Freq   = 1;
    }

    for (/**/; i < NSYM; i++) {
        F[i].Symbol = i;
        F[i].Freq   = 0;
    }

    TotFreq         = max_sym;
    sentinel.Symbol = 0;
    sentinel.Freq   = ((1 << SHIFT) - STEP); // Always first; simplifies sorting.

    F[NSYM].Freq    = 0; // terminates normalize() loop. See below.
}

void FrequencyModel::Normalize() {
    SymFreqs* s;

    /* Faster than F[i].Freq for 0 <= i < NSYM */
    TotFreq = 0;
    for (s = F; s->Freq; s++) {
        s->Freq -= s->Freq >> 1;
        TotFreq += s->Freq;
    }
}

void FrequencyModel::EncodeSymbol(RangeCoder* rc, uint16_t sym) {
    // std::cerr << "encoding symbol: " << sym << " with " << NSYM << "," << SHIFT << "," << STEP << std::endl;
    SymFreqs* s = F;
    uint32_t AccFreq = 0;

    while (s->Symbol != sym) {
        AccFreq += s++->Freq;
        _mm_prefetch((uint8_t*)(s+1), _MM_HINT_T0);
    }

    rc->Encode(AccFreq, s->Freq, TotFreq);

    s->Freq += STEP;
    TotFreq += STEP;

    if (TotFreq > ((1 << SHIFT) - STEP)) {
        //std::cerr << "normalzing" << std::endl;
        Normalize();
    }

    /* Keep approx sorted */
    if(s != F) { // Prevent s[-1] to segfault when s == F
        if (s[0].Freq > s[-1].Freq) {
            SymFreqs t = s[0];
            s[0]  = s[-1];
            s[-1] = t;
        }
    }
}

uint16_t FrequencyModel::DecodeSymbol(RangeCoder *rc) {
    SymFreqs* s = F;
    uint32_t freq = rc->GetFreq(TotFreq);
    uint32_t AccFreq;

    for (AccFreq = 0; (AccFreq += s->Freq) <= freq; s++)
        _mm_prefetch((uint8_t*)s, _MM_HINT_T0);

    AccFreq -= s->Freq;

    rc->Decode(AccFreq, s->Freq, TotFreq);
    s->Freq += STEP;
    TotFreq += STEP;

    if (TotFreq > ((1 << SHIFT) - STEP))
        Normalize();

    /* Keep approx sorted */
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
