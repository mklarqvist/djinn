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
/*
 * Copyright (c) 2013-2019 Genome Research Ltd.
 * Author(s): James Bonfield
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *    Institute nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific
 *    prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Code rewritten to templated C++11 and commented by Marcus D. R. Klarqvist
 */

#ifndef FREQUENCY_MODEL_H_
#define FREQUENCY_MODEL_H_

#ifdef __SSE__
#include <xmmintrin.h>
#else
#define _mm_prefetch(a,b)
#endif

#include "range_coder.h"

namespace pil {

/*
 *--------------------------------------------------------------------------
 * A simple frequency model.
 *
 * Define NSYM to be an integer value before including this file.
 * It will then generate types and functions specific to that
 * maximum number of symbols.
 *
 * This keeps a list of symbols and their frequencies, approximately
 * sorted by symbol frequency. We allow for a single symbol to periodically
 * move up the list when emitted, effectively doing a single step of
 * bubble sort periodically. This means it's largely the same complexity
 * irrespective of alphabet size.
 * It's more efficient on strongly biased distributions than random data.
 *
 * There is no escape symbol, so the model is tailored to relatively
 * stationary samples (although we do have occasional normalisation to
 * avoid frequency counters getting too high).
 *--------------------------------------------------------------------------
 */

template <int NSYM, int STEP = 64, int SHIFT = 11>
class FrequencyModel {
private:
    static constexpr int PIL_FREQUENCY_MODEL_STEP = STEP;
    static constexpr int PIL_FREQUENCY_MODEL_MAX = (1 << SHIFT) - STEP;

    typedef struct {
        uint16_t Freq;
        uint16_t Symbol;
    } SymFreqs;

public:
    FrequencyModel();
    FrequencyModel(int max_sym);
    void Initiate(int max_sym);
    void Normalize();
    void EncodeSymbol(RangeCoder* rc, uint16_t sym);
    void EncodeSymbol(uint16_t sym, uint32_t weight = 1);
    uint16_t DecodeSymbol(RangeCoder* rc);

private:
    uint32_t TotFreq;  // Total frequency

    // Array of Symbols approximately sorted by Freq.
    SymFreqs sentinel, F[NSYM+1];
};

// impl

template <int NSYM, int STEP, int SHIFT>
FrequencyModel<NSYM, STEP, SHIFT>::FrequencyModel() :
    TotFreq(0)
{
    static_assert(NSYM >= 1, "NSYM must be >= 1");
    Initiate(NSYM);
}

template <int NSYM, int STEP, int SHIFT>
FrequencyModel<NSYM, STEP, SHIFT>::FrequencyModel(int max_sym) { Initiate(max_sym); }

template <int NSYM, int STEP, int SHIFT>
void FrequencyModel<NSYM, STEP, SHIFT>::Initiate(int max_sym) {
    int i;

    for (i = 0; i < max_sym; i++) {
        F[i].Symbol = i;
        F[i].Freq   = 1;
    }

    for (; i < NSYM; i++) {
        F[i].Symbol = i;
        F[i].Freq   = 0;
    }

    TotFreq         = max_sym;
    sentinel.Symbol = 0;
    sentinel.Freq   = PIL_FREQUENCY_MODEL_MAX; // Always first; simplifies sorting.

    F[NSYM].Freq    = 0; // terminates normalize() loop. See below.
}

template <int NSYM, int STEP, int SHIFT>
void FrequencyModel<NSYM, STEP, SHIFT>::Normalize() {
    SymFreqs *s;

    /* Faster than F[i].Freq for 0 <= i < NSYM */
    TotFreq = 0;
    for (s = F; s->Freq; s++) {
        s->Freq -= s->Freq >> 1;
        TotFreq += s->Freq;
    }
}

template <int NSYM, int STEP, int SHIFT>
void FrequencyModel<NSYM, STEP, SHIFT>::EncodeSymbol(RangeCoder* rc, uint16_t sym) {
    SymFreqs *s = F;
    uint32_t AccFreq  = 0;

    while (s->Symbol != sym) {
        AccFreq += s++->Freq;
        _mm_prefetch((uint8_t*)(s+1), _MM_HINT_T0);
    }

    rc->Encode(AccFreq, s->Freq, TotFreq);
    s->Freq += PIL_FREQUENCY_MODEL_STEP;
    TotFreq += PIL_FREQUENCY_MODEL_STEP;

    if (TotFreq > PIL_FREQUENCY_MODEL_MAX)
        Normalize();

    /* Keep approx sorted */
    if (s[0].Freq > s[-1].Freq) {
        SymFreqs t = s[0];
        s[0] = s[-1];
        s[-1] = t;
    }
}

template <int NSYM, int STEP, int SHIFT>
void FrequencyModel<NSYM, STEP, SHIFT>::EncodeSymbol(uint16_t sym, uint32_t weight) {
    SymFreqs *s = F;
    uint32_t AccFreq  = 0;

    while (s->Symbol != sym) {
        AccFreq += s++->Freq;
        _mm_prefetch((uint8_t*)(s+1), _MM_HINT_T0);
    }

    s->Freq += PIL_FREQUENCY_MODEL_STEP;
    TotFreq += PIL_FREQUENCY_MODEL_STEP;

    if (TotFreq > PIL_FREQUENCY_MODEL_MAX)
        Normalize();

    /* Keep approx sorted */
    if (s[0].Freq > s[-1].Freq) {
        SymFreqs t = s[0];
        s[0]  = s[-1];
        s[-1] = t;
    }
}

template <int NSYM, int STEP, int SHIFT>
uint16_t FrequencyModel<NSYM, STEP, SHIFT>::DecodeSymbol(RangeCoder *rc) {
    SymFreqs* s = F;
    uint32_t freq = rc->GetFreq(TotFreq);
    uint32_t AccFreq;

    for (AccFreq = 0; (AccFreq += s->Freq) <= freq; s++)
        _mm_prefetch((uint8_t*)s, _MM_HINT_T0);

    AccFreq -= s->Freq;

    rc->Decode(AccFreq, s->Freq, TotFreq);
    s->Freq += PIL_FREQUENCY_MODEL_STEP;
    TotFreq += PIL_FREQUENCY_MODEL_STEP;

    if (TotFreq > PIL_FREQUENCY_MODEL_MAX)
        Normalize();

    /* Keep approx sorted */
    if (s[0].Freq > s[-1].Freq) {
        SymFreqs t = s[0];
        s[0] = s[-1];
        s[-1] = t;
        return t.Symbol;
    }

    return s->Symbol;
}


}

#endif /* FREQUENCY_MODEL_H_ */
