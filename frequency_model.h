/*
* Copyright (c) 2019
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
 * Code rewritten to C++11 and commented by Marcus D. R. Klarqvist
 */

#ifndef FREQUENCY_MODEL_H_
#define FREQUENCY_MODEL_H_

#ifdef __SSE__
#include <xmmintrin.h>
#else
#define _mm_prefetch(a,b)
#endif

#include <vector> // vector
#include <memory> // shared_ptr
#include <cmath> // log2

namespace djinn {

// Based on Subbotin Range Coder.
class RangeCoder {
private:
    static constexpr uint32_t TopValue = 1 << 24;
	static constexpr uint32_t Mask32 = (uint32_t)-1;

public:
    RangeCoder() : low(0), buffer(0), range(Mask32){}
	virtual ~RangeCoder() {}

public:
    void SetInput(uint8_t* in) { out_buf = in_buf = in; }
    void SetOutput(uint8_t* out) { in_buf = out_buf = out; }
    char* GetInput() { return (char *)in_buf; }
    char* GetOutput() { return (char *)out_buf; }
    size_t OutSize() { return out_buf - in_buf; }
    size_t InSize() { return in_buf - out_buf; }

	void StartEncode() {
		low = 0;
		range = Mask32;
	}

	void Encode(uint32_t cumFreq, uint32_t symFreq, uint32_t totalFreqSum) {
        assert(range > totalFreqSum);
		range /= totalFreqSum;
        low += range * cumFreq;
        range *= symFreq;

        while(range < TopValue) {
            assert(range != 0);
            // range = 0x00ffffff..
            // low/high may be matching
            //       eg 88332211/88342211 (range 00010000)
            // or differing
            //       eg 88ff2211/89002211 (range 00010000)
            //
            // If the latter, we need to reduce range down
            // such that high=88ffffff.
            // Eg. top-1      == 00ffffff
            //     low|top-1  == 88ffffff
            //     ...-low    == 0000ddee
            if ( (uint8_t)((low ^ (low + range)) >> 56) )
                range = (((uint32_t)(low) | (TopValue - 1)) - (uint32_t)(low));
            *out_buf++ = low >> 56, range <<= 8, low <<= 8;
        }
	}

	void FinishEncode() {
		for (int i = 0; i < 8; i++) {
            *out_buf++ = (uint8_t)(low >> 56);
			low <<= 8;
		}
	}

	void StartDecode() {
		buffer = 0;
		for (uint32_t i = 1; i <= 8; ++i) {
			buffer |= (uint64_t)*in_buf << (64 - i * 8);
            ++in_buf;
		}

		low = 0;
		range = Mask32;
	}

	inline uint32_t GetFreq(uint32_t totalFreq) {
		assert(totalFreq != 0);
		return (uint32_t) (buffer / (range /= totalFreq));
	}

	void Decode(uint32_t lowEnd, uint32_t symFreq, uint32_t /*totalFreq_*/) {
		uint32_t r = lowEnd * range;
		buffer -= r;
		low += r;
		range *= symFreq;

		while (range < TopValue) {
			if ( (uint8_t)((low ^ (low + range)) >> 56) )
				range = (((uint32_t)(low) | (TopValue - 1)) - (uint32_t)(low));

			buffer = (buffer << 8) + *in_buf++;
			low <<= 8, range <<= 8;
		}
	}

	void FinishDecode() {}

public:
    uint64_t low, buffer;
	uint32_t range;
    uint8_t* in_buf;
    uint8_t* out_buf;
};

/*
 *--------------------------------------------------------------------------
 * A simple frequency model.
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

class FrequencyModel {
private:
    struct SymFreqs {
        bool operator<(const SymFreqs& other) const { return(Symbol < other.Symbol); }
        
        uint32_t Freq;
        uint16_t Symbol;
    };

public:
    FrequencyModel();
    FrequencyModel(int nsym, int max_sym);
    FrequencyModel(int nsym, int max_sym, int step, int shift);
    ~FrequencyModel();

    void Initiate(int max_sym);
    void Initiate(int nsym, int max_sym);
    void Normalize();
    void EncodeSymbol(RangeCoder* rc, uint16_t sym);
    uint16_t DecodeSymbol(RangeCoder* rc);
    void EncodeSymbol(uint16_t sym);
    double GetP(uint16_t sym) const;

public:
    uint32_t n_symbols;
	uint32_t step_size;
	uint32_t max_total;
	uint32_t total_frequency;

    // Array of Symbols approximately sorted by Freq.
    SymFreqs sentinel;
    SymFreqs* F;
};

/*======   Context model container   ======*/

class GeneralModel {
public:
    GeneralModel() noexcept;
    GeneralModel(int n_symbols, int model_size);
    GeneralModel(int n_symbols, int model_size, std::shared_ptr<RangeCoder> rc);
    GeneralModel(int n_symbols, int model_size, int shift, int step);
    GeneralModel(int n_symbols, int model_size, int shift, int step, std::shared_ptr<RangeCoder> rc);
    ~GeneralModel();

    int Initiate(int n_symbols, int model_size);
    int Initiate(int n_symbols, int model_size, std::shared_ptr<RangeCoder> rc);
    int Initiate(int n_symbols, int model_size, int shift, int step);
    int Initiate(int n_symbols, int model_size, int shift, int step, std::shared_ptr<RangeCoder> rc);

    int FinishEncoding();
    int FinishDecoding();
    void StartEncoding();
    void StartDecoding(uint8_t* data);

    void EncodeSymbol(const uint16_t symbol);
    void EncodeSymbolNoUpdate(const uint16_t symbol);
    
    uint16_t DecodeSymbol();
    uint16_t DecodeSymbolNoUpdate();

    void ResetModels();
    void ResetContext();
    void Reset();

public:
    int max_model_symbols;
    int model_context_shift;
    uint32_t model_context, model_ctx_mask;
    std::shared_ptr<RangeCoder> range_coder;
    std::vector < std::shared_ptr<FrequencyModel> > models;
    size_t n_additions; // number of updates performed
    size_t n_buffer; // buffer size
    uint8_t* buffer; // buffer. todo: fixme
};

}

#endif /* FREQUENCY_MODEL_H_ */
