#include "range_coder.h"

#include <stdio.h> //frpint
#include <stdlib.h>
#include <cstdlib> // abort

#ifdef __SSE__
#include <xmmintrin.h>
#else
#define _mm_prefetch(a,b)
#endif

#define  DO(n) int _;for (_=0; _<n; _++)
#define  TOP   (1<<24)

namespace pil {

void RangeCoder::StartEncode()
{
    low = 0;
    range = (uint32_t)-1;
}

void RangeCoder::StartDecode()
{
    code = low = 0;
    range = (uint32_t)-1;
    DO(8) code = (code << 8) | *in_buf++;
}

void RangeCoder::FinishEncode()
{
    DO(8) (*out_buf++ = low >> 56), low <<= 8;
}

void RangeCoder::FinishDecode() {}

#if DIV_RCP
static uint64_t rcp[65536];
void RangeCoder::build_rcp_freq(void) {
    int i;
    for (i = 1; i < 65536; i++)
    rcp[i] = (((uint64_t)1<<32)) / i;
}
#else
void RangeCoder::build_rcp_freq(void) {}
#endif

void RangeCoder::Encode(uint32_t cumFreq, uint32_t freq, uint32_t totFreq)
{
    //fprintf(stderr, "                       RC %d+%d of %d\n", cumFreq, freq, totFreq);

    // division-less doesn't help much in this case as only a single division anyway
    //uint32_t d = (rcp[totFreq] * range)>>32;
    //d += (range - totFreq * d >= totFreq);
    //low  += cumFreq * (range = d);

    low   += cumFreq * (range /= totFreq);
    range *= freq;

    if (cumFreq + freq > totFreq) {
        fprintf(stderr, "cumFreq %d + freq %d > tot %d\n", cumFreq, freq, totFreq);
        std::abort();
    }

    while( range<TOP ) {
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
            range = (((uint32_t)(low) | (TOP - 1)) - (uint32_t)(low));
        *out_buf++ = low >> 56, range <<= 8, low <<= 8;
    }
}

uint32_t RangeCoder::GetFreq(uint32_t totFreq) {
#if DIV_RCP
    // Replacing two divisions by one is beneficial as they get stuck waiting.
    // 2.53s
    uint32_t d = (rcp[totFreq] * range)>>32;
    d += (range - totFreq * d >= totFreq);
    return code/(range=d);
#else
    // 2.67s
    return code/(range/=totFreq);
#endif
}

void RangeCoder::Decode(uint32_t cumFreq, uint32_t freq, uint32_t totFreq)
{
    //fprintf(stderr, "                       RC %d+%d of %d\n", cumFreq, freq, totFreq);

    uint32_t temp = cumFreq * range;
    low   += temp;
    code  -= temp;
    range *= freq;

    while( range<TOP ) {
        if ( (uint8_t)((low ^ (low + range)) >> 56) )
            range = (((uint32_t)(low) | (TOP - 1)) - (uint32_t)(low));
        code = (code << 8) | *in_buf++, range <<= 8, low <<= 8;
    }
}

}
