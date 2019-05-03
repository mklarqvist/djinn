/*
 * AVX2 edition using 32 RANS states.  This uses a shared pointer for the
 * compressed buffer.
 *
 * TODO: implement SIMD version of the order-1 encoder (decoder is done).
 */

#include "r32x16b_avx2.h"
#include <x86intrin.h>
#include <cstring>

#ifndef TF_SHIFT
#  define TF_SHIFT 12
#endif
#define TOTFREQ (1<<TF_SHIFT)

// 9 is considerably faster on some data sets due to reduced table size.
#ifndef TF_SHIFT_O1
#  define TF_SHIFT_O1 9
#endif
#define TOTFREQ_O1 (1<<TF_SHIFT_O1)

/*
 * These tables are 8k.  On older systems with small L1 cache, this may be
 * a problem.
 *
 * #define PM(a,b,c,d,e,f,g,h) ((a<<0)|(b<<4)|(c<<8)|(d<<12)|(e<<16)|(f<<20)|(g<<24)|(h<<28))
 *
 * Instead of permute via
 *    __m256i idx1 = _mm256_loadu_si256((const __m256i*)permute[imask1]);
 *
 * we can pack the indices and shift them back again
 *   __m256i idx1 = _mm256_srlv_epi32(_mm256_set1_epi32(permute2[imask1]),
 *                                    _mm256_set_epi32(28,24,20,16,12,8,4,0));
 *
 * However on my Haswell system this slows down r32x16b_avx2 from 1440 to
 * 1200 MB/s decode speeds.
 * It's much closer for order-1 decoder, but still doesn't help.
 *
 * The encoder side seems to make no difference either way or be very marginal.
 */

#define _ 9
static uint32_t permute[256][8] __attribute__((aligned(32))) = { // reverse binary bit order
  { _,_,_,_,_,_,_,_,},
  { 0,_,_,_,_,_,_,_,},
  { _,0,_,_,_,_,_,_,},
  { 0,1,_,_,_,_,_,_,},
  { _,_,0,_,_,_,_,_,},
  { 0,_,1,_,_,_,_,_,},
  { _,0,1,_,_,_,_,_,},
  { 0,1,2,_,_,_,_,_,},
  { _,_,_,0,_,_,_,_,},
  { 0,_,_,1,_,_,_,_,},
  { _,0,_,1,_,_,_,_,},
  { 0,1,_,2,_,_,_,_,},
  { _,_,0,1,_,_,_,_,},
  { 0,_,1,2,_,_,_,_,},
  { _,0,1,2,_,_,_,_,},
  { 0,1,2,3,_,_,_,_,},
  { _,_,_,_,0,_,_,_,},
  { 0,_,_,_,1,_,_,_,},
  { _,0,_,_,1,_,_,_,},
  { 0,1,_,_,2,_,_,_,},
  { _,_,0,_,1,_,_,_,},
  { 0,_,1,_,2,_,_,_,},
  { _,0,1,_,2,_,_,_,},
  { 0,1,2,_,3,_,_,_,},
  { _,_,_,0,1,_,_,_,},
  { 0,_,_,1,2,_,_,_,},
  { _,0,_,1,2,_,_,_,},
  { 0,1,_,2,3,_,_,_,},
  { _,_,0,1,2,_,_,_,},
  { 0,_,1,2,3,_,_,_,},
  { _,0,1,2,3,_,_,_,},
  { 0,1,2,3,4,_,_,_,},
  { _,_,_,_,_,0,_,_,},
  { 0,_,_,_,_,1,_,_,},
  { _,0,_,_,_,1,_,_,},
  { 0,1,_,_,_,2,_,_,},
  { _,_,0,_,_,1,_,_,},
  { 0,_,1,_,_,2,_,_,},
  { _,0,1,_,_,2,_,_,},
  { 0,1,2,_,_,3,_,_,},
  { _,_,_,0,_,1,_,_,},
  { 0,_,_,1,_,2,_,_,},
  { _,0,_,1,_,2,_,_,},
  { 0,1,_,2,_,3,_,_,},
  { _,_,0,1,_,2,_,_,},
  { 0,_,1,2,_,3,_,_,},
  { _,0,1,2,_,3,_,_,},
  { 0,1,2,3,_,4,_,_,},
  { _,_,_,_,0,1,_,_,},
  { 0,_,_,_,1,2,_,_,},
  { _,0,_,_,1,2,_,_,},
  { 0,1,_,_,2,3,_,_,},
  { _,_,0,_,1,2,_,_,},
  { 0,_,1,_,2,3,_,_,},
  { _,0,1,_,2,3,_,_,},
  { 0,1,2,_,3,4,_,_,},
  { _,_,_,0,1,2,_,_,},
  { 0,_,_,1,2,3,_,_,},
  { _,0,_,1,2,3,_,_,},
  { 0,1,_,2,3,4,_,_,},
  { _,_,0,1,2,3,_,_,},
  { 0,_,1,2,3,4,_,_,},
  { _,0,1,2,3,4,_,_,},
  { 0,1,2,3,4,5,_,_,},
  { _,_,_,_,_,_,0,_,},
  { 0,_,_,_,_,_,1,_,},
  { _,0,_,_,_,_,1,_,},
  { 0,1,_,_,_,_,2,_,},
  { _,_,0,_,_,_,1,_,},
  { 0,_,1,_,_,_,2,_,},
  { _,0,1,_,_,_,2,_,},
  { 0,1,2,_,_,_,3,_,},
  { _,_,_,0,_,_,1,_,},
  { 0,_,_,1,_,_,2,_,},
  { _,0,_,1,_,_,2,_,},
  { 0,1,_,2,_,_,3,_,},
  { _,_,0,1,_,_,2,_,},
  { 0,_,1,2,_,_,3,_,},
  { _,0,1,2,_,_,3,_,},
  { 0,1,2,3,_,_,4,_,},
  { _,_,_,_,0,_,1,_,},
  { 0,_,_,_,1,_,2,_,},
  { _,0,_,_,1,_,2,_,},
  { 0,1,_,_,2,_,3,_,},
  { _,_,0,_,1,_,2,_,},
  { 0,_,1,_,2,_,3,_,},
  { _,0,1,_,2,_,3,_,},
  { 0,1,2,_,3,_,4,_,},
  { _,_,_,0,1,_,2,_,},
  { 0,_,_,1,2,_,3,_,},
  { _,0,_,1,2,_,3,_,},
  { 0,1,_,2,3,_,4,_,},
  { _,_,0,1,2,_,3,_,},
  { 0,_,1,2,3,_,4,_,},
  { _,0,1,2,3,_,4,_,},
  { 0,1,2,3,4,_,5,_,},
  { _,_,_,_,_,0,1,_,},
  { 0,_,_,_,_,1,2,_,},
  { _,0,_,_,_,1,2,_,},
  { 0,1,_,_,_,2,3,_,},
  { _,_,0,_,_,1,2,_,},
  { 0,_,1,_,_,2,3,_,},
  { _,0,1,_,_,2,3,_,},
  { 0,1,2,_,_,3,4,_,},
  { _,_,_,0,_,1,2,_,},
  { 0,_,_,1,_,2,3,_,},
  { _,0,_,1,_,2,3,_,},
  { 0,1,_,2,_,3,4,_,},
  { _,_,0,1,_,2,3,_,},
  { 0,_,1,2,_,3,4,_,},
  { _,0,1,2,_,3,4,_,},
  { 0,1,2,3,_,4,5,_,},
  { _,_,_,_,0,1,2,_,},
  { 0,_,_,_,1,2,3,_,},
  { _,0,_,_,1,2,3,_,},
  { 0,1,_,_,2,3,4,_,},
  { _,_,0,_,1,2,3,_,},
  { 0,_,1,_,2,3,4,_,},
  { _,0,1,_,2,3,4,_,},
  { 0,1,2,_,3,4,5,_,},
  { _,_,_,0,1,2,3,_,},
  { 0,_,_,1,2,3,4,_,},
  { _,0,_,1,2,3,4,_,},
  { 0,1,_,2,3,4,5,_,},
  { _,_,0,1,2,3,4,_,},
  { 0,_,1,2,3,4,5,_,},
  { _,0,1,2,3,4,5,_,},
  { 0,1,2,3,4,5,6,_,},
  { _,_,_,_,_,_,_,0,},
  { 0,_,_,_,_,_,_,1,},
  { _,0,_,_,_,_,_,1,},
  { 0,1,_,_,_,_,_,2,},
  { _,_,0,_,_,_,_,1,},
  { 0,_,1,_,_,_,_,2,},
  { _,0,1,_,_,_,_,2,},
  { 0,1,2,_,_,_,_,3,},
  { _,_,_,0,_,_,_,1,},
  { 0,_,_,1,_,_,_,2,},
  { _,0,_,1,_,_,_,2,},
  { 0,1,_,2,_,_,_,3,},
  { _,_,0,1,_,_,_,2,},
  { 0,_,1,2,_,_,_,3,},
  { _,0,1,2,_,_,_,3,},
  { 0,1,2,3,_,_,_,4,},
  { _,_,_,_,0,_,_,1,},
  { 0,_,_,_,1,_,_,2,},
  { _,0,_,_,1,_,_,2,},
  { 0,1,_,_,2,_,_,3,},
  { _,_,0,_,1,_,_,2,},
  { 0,_,1,_,2,_,_,3,},
  { _,0,1,_,2,_,_,3,},
  { 0,1,2,_,3,_,_,4,},
  { _,_,_,0,1,_,_,2,},
  { 0,_,_,1,2,_,_,3,},
  { _,0,_,1,2,_,_,3,},
  { 0,1,_,2,3,_,_,4,},
  { _,_,0,1,2,_,_,3,},
  { 0,_,1,2,3,_,_,4,},
  { _,0,1,2,3,_,_,4,},
  { 0,1,2,3,4,_,_,5,},
  { _,_,_,_,_,0,_,1,},
  { 0,_,_,_,_,1,_,2,},
  { _,0,_,_,_,1,_,2,},
  { 0,1,_,_,_,2,_,3,},
  { _,_,0,_,_,1,_,2,},
  { 0,_,1,_,_,2,_,3,},
  { _,0,1,_,_,2,_,3,},
  { 0,1,2,_,_,3,_,4,},
  { _,_,_,0,_,1,_,2,},
  { 0,_,_,1,_,2,_,3,},
  { _,0,_,1,_,2,_,3,},
  { 0,1,_,2,_,3,_,4,},
  { _,_,0,1,_,2,_,3,},
  { 0,_,1,2,_,3,_,4,},
  { _,0,1,2,_,3,_,4,},
  { 0,1,2,3,_,4,_,5,},
  { _,_,_,_,0,1,_,2,},
  { 0,_,_,_,1,2,_,3,},
  { _,0,_,_,1,2,_,3,},
  { 0,1,_,_,2,3,_,4,},
  { _,_,0,_,1,2,_,3,},
  { 0,_,1,_,2,3,_,4,},
  { _,0,1,_,2,3,_,4,},
  { 0,1,2,_,3,4,_,5,},
  { _,_,_,0,1,2,_,3,},
  { 0,_,_,1,2,3,_,4,},
  { _,0,_,1,2,3,_,4,},
  { 0,1,_,2,3,4,_,5,},
  { _,_,0,1,2,3,_,4,},
  { 0,_,1,2,3,4,_,5,},
  { _,0,1,2,3,4,_,5,},
  { 0,1,2,3,4,5,_,6,},
  { _,_,_,_,_,_,0,1,},
  { 0,_,_,_,_,_,1,2,},
  { _,0,_,_,_,_,1,2,},
  { 0,1,_,_,_,_,2,3,},
  { _,_,0,_,_,_,1,2,},
  { 0,_,1,_,_,_,2,3,},
  { _,0,1,_,_,_,2,3,},
  { 0,1,2,_,_,_,3,4,},
  { _,_,_,0,_,_,1,2,},
  { 0,_,_,1,_,_,2,3,},
  { _,0,_,1,_,_,2,3,},
  { 0,1,_,2,_,_,3,4,},
  { _,_,0,1,_,_,2,3,},
  { 0,_,1,2,_,_,3,4,},
  { _,0,1,2,_,_,3,4,},
  { 0,1,2,3,_,_,4,5,},
  { _,_,_,_,0,_,1,2,},
  { 0,_,_,_,1,_,2,3,},
  { _,0,_,_,1,_,2,3,},
  { 0,1,_,_,2,_,3,4,},
  { _,_,0,_,1,_,2,3,},
  { 0,_,1,_,2,_,3,4,},
  { _,0,1,_,2,_,3,4,},
  { 0,1,2,_,3,_,4,5,},
  { _,_,_,0,1,_,2,3,},
  { 0,_,_,1,2,_,3,4,},
  { _,0,_,1,2,_,3,4,},
  { 0,1,_,2,3,_,4,5,},
  { _,_,0,1,2,_,3,4,},
  { 0,_,1,2,3,_,4,5,},
  { _,0,1,2,3,_,4,5,},
  { 0,1,2,3,4,_,5,6,},
  { _,_,_,_,_,0,1,2,},
  { 0,_,_,_,_,1,2,3,},
  { _,0,_,_,_,1,2,3,},
  { 0,1,_,_,_,2,3,4,},
  { _,_,0,_,_,1,2,3,},
  { 0,_,1,_,_,2,3,4,},
  { _,0,1,_,_,2,3,4,},
  { 0,1,2,_,_,3,4,5,},
  { _,_,_,0,_,1,2,3,},
  { 0,_,_,1,_,2,3,4,},
  { _,0,_,1,_,2,3,4,},
  { 0,1,_,2,_,3,4,5,},
  { _,_,0,1,_,2,3,4,},
  { 0,_,1,2,_,3,4,5,},
  { _,0,1,2,_,3,4,5,},
  { 0,1,2,3,_,4,5,6,},
  { _,_,_,_,0,1,2,3,},
  { 0,_,_,_,1,2,3,4,},
  { _,0,_,_,1,2,3,4,},
  { 0,1,_,_,2,3,4,5,},
  { _,_,0,_,1,2,3,4,},
  { 0,_,1,_,2,3,4,5,},
  { _,0,1,_,2,3,4,5,},
  { 0,1,2,_,3,4,5,6,},
  { _,_,_,0,1,2,3,4,},
  { 0,_,_,1,2,3,4,5,},
  { _,0,_,1,2,3,4,5,},
  { 0,1,_,2,3,4,5,6,},
  { _,_,0,1,2,3,4,5,},
  { 0,_,1,2,3,4,5,6,},
  { _,0,1,2,3,4,5,6,},
  { 0,1,2,3,4,5,6,7,},
};

static uint32_t permutec[256][8] __attribute__((aligned(32))) = { // reverse binary bit order
  { _,_,_,_,_,_,_,_,},
  { _,_,_,_,_,_,_,0,},
  { _,_,_,_,_,_,_,1,},
  { _,_,_,_,_,_,0,1,},
  { _,_,_,_,_,_,_,2,},
  { _,_,_,_,_,_,0,2,},
  { _,_,_,_,_,_,1,2,},
  { _,_,_,_,_,0,1,2,},
  { _,_,_,_,_,_,_,3,},
  { _,_,_,_,_,_,0,3,},
  { _,_,_,_,_,_,1,3,},
  { _,_,_,_,_,0,1,3,},
  { _,_,_,_,_,_,2,3,},
  { _,_,_,_,_,0,2,3,},
  { _,_,_,_,_,1,2,3,},
  { _,_,_,_,0,1,2,3,},
  { _,_,_,_,_,_,_,4,},
  { _,_,_,_,_,_,0,4,},
  { _,_,_,_,_,_,1,4,},
  { _,_,_,_,_,0,1,4,},
  { _,_,_,_,_,_,2,4,},
  { _,_,_,_,_,0,2,4,},
  { _,_,_,_,_,1,2,4,},
  { _,_,_,_,0,1,2,4,},
  { _,_,_,_,_,_,3,4,},
  { _,_,_,_,_,0,3,4,},
  { _,_,_,_,_,1,3,4,},
  { _,_,_,_,0,1,3,4,},
  { _,_,_,_,_,2,3,4,},
  { _,_,_,_,0,2,3,4,},
  { _,_,_,_,1,2,3,4,},
  { _,_,_,0,1,2,3,4,},
  { _,_,_,_,_,_,_,5,},
  { _,_,_,_,_,_,0,5,},
  { _,_,_,_,_,_,1,5,},
  { _,_,_,_,_,0,1,5,},
  { _,_,_,_,_,_,2,5,},
  { _,_,_,_,_,0,2,5,},
  { _,_,_,_,_,1,2,5,},
  { _,_,_,_,0,1,2,5,},
  { _,_,_,_,_,_,3,5,},
  { _,_,_,_,_,0,3,5,},
  { _,_,_,_,_,1,3,5,},
  { _,_,_,_,0,1,3,5,},
  { _,_,_,_,_,2,3,5,},
  { _,_,_,_,0,2,3,5,},
  { _,_,_,_,1,2,3,5,},
  { _,_,_,0,1,2,3,5,},
  { _,_,_,_,_,_,4,5,},
  { _,_,_,_,_,0,4,5,},
  { _,_,_,_,_,1,4,5,},
  { _,_,_,_,0,1,4,5,},
  { _,_,_,_,_,2,4,5,},
  { _,_,_,_,0,2,4,5,},
  { _,_,_,_,1,2,4,5,},
  { _,_,_,0,1,2,4,5,},
  { _,_,_,_,_,3,4,5,},
  { _,_,_,_,0,3,4,5,},
  { _,_,_,_,1,3,4,5,},
  { _,_,_,0,1,3,4,5,},
  { _,_,_,_,2,3,4,5,},
  { _,_,_,0,2,3,4,5,},
  { _,_,_,1,2,3,4,5,},
  { _,_,0,1,2,3,4,5,},
  { _,_,_,_,_,_,_,6,},
  { _,_,_,_,_,_,0,6,},
  { _,_,_,_,_,_,1,6,},
  { _,_,_,_,_,0,1,6,},
  { _,_,_,_,_,_,2,6,},
  { _,_,_,_,_,0,2,6,},
  { _,_,_,_,_,1,2,6,},
  { _,_,_,_,0,1,2,6,},
  { _,_,_,_,_,_,3,6,},
  { _,_,_,_,_,0,3,6,},
  { _,_,_,_,_,1,3,6,},
  { _,_,_,_,0,1,3,6,},
  { _,_,_,_,_,2,3,6,},
  { _,_,_,_,0,2,3,6,},
  { _,_,_,_,1,2,3,6,},
  { _,_,_,0,1,2,3,6,},
  { _,_,_,_,_,_,4,6,},
  { _,_,_,_,_,0,4,6,},
  { _,_,_,_,_,1,4,6,},
  { _,_,_,_,0,1,4,6,},
  { _,_,_,_,_,2,4,6,},
  { _,_,_,_,0,2,4,6,},
  { _,_,_,_,1,2,4,6,},
  { _,_,_,0,1,2,4,6,},
  { _,_,_,_,_,3,4,6,},
  { _,_,_,_,0,3,4,6,},
  { _,_,_,_,1,3,4,6,},
  { _,_,_,0,1,3,4,6,},
  { _,_,_,_,2,3,4,6,},
  { _,_,_,0,2,3,4,6,},
  { _,_,_,1,2,3,4,6,},
  { _,_,0,1,2,3,4,6,},
  { _,_,_,_,_,_,5,6,},
  { _,_,_,_,_,0,5,6,},
  { _,_,_,_,_,1,5,6,},
  { _,_,_,_,0,1,5,6,},
  { _,_,_,_,_,2,5,6,},
  { _,_,_,_,0,2,5,6,},
  { _,_,_,_,1,2,5,6,},
  { _,_,_,0,1,2,5,6,},
  { _,_,_,_,_,3,5,6,},
  { _,_,_,_,0,3,5,6,},
  { _,_,_,_,1,3,5,6,},
  { _,_,_,0,1,3,5,6,},
  { _,_,_,_,2,3,5,6,},
  { _,_,_,0,2,3,5,6,},
  { _,_,_,1,2,3,5,6,},
  { _,_,0,1,2,3,5,6,},
  { _,_,_,_,_,4,5,6,},
  { _,_,_,_,0,4,5,6,},
  { _,_,_,_,1,4,5,6,},
  { _,_,_,0,1,4,5,6,},
  { _,_,_,_,2,4,5,6,},
  { _,_,_,0,2,4,5,6,},
  { _,_,_,1,2,4,5,6,},
  { _,_,0,1,2,4,5,6,},
  { _,_,_,_,3,4,5,6,},
  { _,_,_,0,3,4,5,6,},
  { _,_,_,1,3,4,5,6,},
  { _,_,0,1,3,4,5,6,},
  { _,_,_,2,3,4,5,6,},
  { _,_,0,2,3,4,5,6,},
  { _,_,1,2,3,4,5,6,},
  { _,0,1,2,3,4,5,6,},
  { _,_,_,_,_,_,_,7,},
  { _,_,_,_,_,_,0,7,},
  { _,_,_,_,_,_,1,7,},
  { _,_,_,_,_,0,1,7,},
  { _,_,_,_,_,_,2,7,},
  { _,_,_,_,_,0,2,7,},
  { _,_,_,_,_,1,2,7,},
  { _,_,_,_,0,1,2,7,},
  { _,_,_,_,_,_,3,7,},
  { _,_,_,_,_,0,3,7,},
  { _,_,_,_,_,1,3,7,},
  { _,_,_,_,0,1,3,7,},
  { _,_,_,_,_,2,3,7,},
  { _,_,_,_,0,2,3,7,},
  { _,_,_,_,1,2,3,7,},
  { _,_,_,0,1,2,3,7,},
  { _,_,_,_,_,_,4,7,},
  { _,_,_,_,_,0,4,7,},
  { _,_,_,_,_,1,4,7,},
  { _,_,_,_,0,1,4,7,},
  { _,_,_,_,_,2,4,7,},
  { _,_,_,_,0,2,4,7,},
  { _,_,_,_,1,2,4,7,},
  { _,_,_,0,1,2,4,7,},
  { _,_,_,_,_,3,4,7,},
  { _,_,_,_,0,3,4,7,},
  { _,_,_,_,1,3,4,7,},
  { _,_,_,0,1,3,4,7,},
  { _,_,_,_,2,3,4,7,},
  { _,_,_,0,2,3,4,7,},
  { _,_,_,1,2,3,4,7,},
  { _,_,0,1,2,3,4,7,},
  { _,_,_,_,_,_,5,7,},
  { _,_,_,_,_,0,5,7,},
  { _,_,_,_,_,1,5,7,},
  { _,_,_,_,0,1,5,7,},
  { _,_,_,_,_,2,5,7,},
  { _,_,_,_,0,2,5,7,},
  { _,_,_,_,1,2,5,7,},
  { _,_,_,0,1,2,5,7,},
  { _,_,_,_,_,3,5,7,},
  { _,_,_,_,0,3,5,7,},
  { _,_,_,_,1,3,5,7,},
  { _,_,_,0,1,3,5,7,},
  { _,_,_,_,2,3,5,7,},
  { _,_,_,0,2,3,5,7,},
  { _,_,_,1,2,3,5,7,},
  { _,_,0,1,2,3,5,7,},
  { _,_,_,_,_,4,5,7,},
  { _,_,_,_,0,4,5,7,},
  { _,_,_,_,1,4,5,7,},
  { _,_,_,0,1,4,5,7,},
  { _,_,_,_,2,4,5,7,},
  { _,_,_,0,2,4,5,7,},
  { _,_,_,1,2,4,5,7,},
  { _,_,0,1,2,4,5,7,},
  { _,_,_,_,3,4,5,7,},
  { _,_,_,0,3,4,5,7,},
  { _,_,_,1,3,4,5,7,},
  { _,_,0,1,3,4,5,7,},
  { _,_,_,2,3,4,5,7,},
  { _,_,0,2,3,4,5,7,},
  { _,_,1,2,3,4,5,7,},
  { _,0,1,2,3,4,5,7,},
  { _,_,_,_,_,_,6,7,},
  { _,_,_,_,_,0,6,7,},
  { _,_,_,_,_,1,6,7,},
  { _,_,_,_,0,1,6,7,},
  { _,_,_,_,_,2,6,7,},
  { _,_,_,_,0,2,6,7,},
  { _,_,_,_,1,2,6,7,},
  { _,_,_,0,1,2,6,7,},
  { _,_,_,_,_,3,6,7,},
  { _,_,_,_,0,3,6,7,},
  { _,_,_,_,1,3,6,7,},
  { _,_,_,0,1,3,6,7,},
  { _,_,_,_,2,3,6,7,},
  { _,_,_,0,2,3,6,7,},
  { _,_,_,1,2,3,6,7,},
  { _,_,0,1,2,3,6,7,},
  { _,_,_,_,_,4,6,7,},
  { _,_,_,_,0,4,6,7,},
  { _,_,_,_,1,4,6,7,},
  { _,_,_,0,1,4,6,7,},
  { _,_,_,_,2,4,6,7,},
  { _,_,_,0,2,4,6,7,},
  { _,_,_,1,2,4,6,7,},
  { _,_,0,1,2,4,6,7,},
  { _,_,_,_,3,4,6,7,},
  { _,_,_,0,3,4,6,7,},
  { _,_,_,1,3,4,6,7,},
  { _,_,0,1,3,4,6,7,},
  { _,_,_,2,3,4,6,7,},
  { _,_,0,2,3,4,6,7,},
  { _,_,1,2,3,4,6,7,},
  { _,0,1,2,3,4,6,7,},
  { _,_,_,_,_,5,6,7,},
  { _,_,_,_,0,5,6,7,},
  { _,_,_,_,1,5,6,7,},
  { _,_,_,0,1,5,6,7,},
  { _,_,_,_,2,5,6,7,},
  { _,_,_,0,2,5,6,7,},
  { _,_,_,1,2,5,6,7,},
  { _,_,0,1,2,5,6,7,},
  { _,_,_,_,3,5,6,7,},
  { _,_,_,0,3,5,6,7,},
  { _,_,_,1,3,5,6,7,},
  { _,_,0,1,3,5,6,7,},
  { _,_,_,2,3,5,6,7,},
  { _,_,0,2,3,5,6,7,},
  { _,_,1,2,3,5,6,7,},
  { _,0,1,2,3,5,6,7,},
  { _,_,_,_,4,5,6,7,},
  { _,_,_,0,4,5,6,7,},
  { _,_,_,1,4,5,6,7,},
  { _,_,0,1,4,5,6,7,},
  { _,_,_,2,4,5,6,7,},
  { _,_,0,2,4,5,6,7,},
  { _,_,1,2,4,5,6,7,},
  { _,0,1,2,4,5,6,7,},
  { _,_,_,3,4,5,6,7,},
  { _,_,0,3,4,5,6,7,},
  { _,_,1,3,4,5,6,7,},
  { _,0,1,3,4,5,6,7,},
  { _,_,2,3,4,5,6,7,},
  { _,0,2,3,4,5,6,7,},
  { _,1,2,3,4,5,6,7,},
  { 0,1,2,3,4,5,6,7,},
};
#undef _


/*-----------------------------------------------------------------------------
 * Memory to memory compression functions.
 *
 * These are original versions without any manual loop unrolling. They
 * are easier to understand, but can be up to 2x slower.
 */

#define MAGIC_RANS 8

static void hist8(unsigned char *in, unsigned int in_size, int F0[256]) {
    int F1[256+MAGIC_RANS] = {0}, F2[256+MAGIC_RANS] = {0}, F3[256+MAGIC_RANS] = {0};
    int F4[256+MAGIC_RANS] = {0}, F5[256+MAGIC_RANS] = {0}, F6[256+MAGIC_RANS] = {0}, F7[256+MAGIC_RANS] = {0};
    int i, i4 = ((in_size-4) & ~7)/4; // permits vnext
    uint32_t *in4 = (uint32_t *)in;
    uint32_t vnext = i4 ? in4[0] : 0;
    for (i = 0; i < i4; i+=2) {
	uint32_t v = vnext; vnext = in4[i+1];
	F0[(unsigned char)(v>> 0)]++;
	F1[(unsigned char)(v>> 8)]++;
	F2[(unsigned char)(v>>16)]++;
	F3[(unsigned char)(v>>24)]++;
	v = vnext; vnext = in4[i+2];
	F4[(unsigned char)(v>> 0)]++;
	F5[(unsigned char)(v>> 8)]++;
	F6[(unsigned char)(v>>16)]++;
	F7[(unsigned char)(v>>24)]++;
    }

    i *= 4;
    while (i < in_size)
	F0[in[i++]]++;

    for (i = 0; i < 256; i++)
	F0[i] += F1[i] + F2[i] + F3[i] + F4[i] + F5[i] + F6[i] + F7[i];
}

unsigned int rans_compress_bound_32x16(unsigned int size, int order, int *tab) {
  int tabsz = order == 0
    ?     257*3 + 4 + NX*4 + NX*4
    : 257*257*3 + 4 + NX*4 + NX*4;
  if (tab) *tab = tabsz;
  return 1.05*size + NX*4 + tabsz;
}


static __m256i _mm256_mulhi_epu32(__m256i a, __m256i b) {
    // Multiply bottom 4 items and top 4 items together.
    __m256i ab_hm = _mm256_mul_epu32(_mm256_srli_epi64(a, 32), _mm256_srli_epi64(b, 32));
    __m256i ab_lm = _mm256_srli_epi64(_mm256_mul_epu32(a, b), 32);

    // Shift to get hi 32-bit of each 64-bit product
    ab_hm = _mm256_and_si256(ab_hm, _mm256_set1_epi64x((uint64_t)0xffffffff00000000));

    return _mm256_or_si256(ab_lm, ab_hm);
}

__m256i _mm256_div_epi32(__m256i a, __m256i b) {
    __m256d a1  = _mm256_cvtepi32_pd(_mm256_extracti128_si256(a, 0));
    __m256d b1  = _mm256_cvtepi32_pd(_mm256_extracti128_si256(b, 0));
    __m128i ab1 = _mm256_cvttpd_epi32(_mm256_div_pd(a1, b1));

    __m256d a2  = _mm256_cvtepi32_pd(_mm256_extracti128_si256(a, 1));
    __m256d b2  = _mm256_cvtepi32_pd(_mm256_extracti128_si256(b, 1));
    __m128i ab2 = _mm256_cvttpd_epi32(_mm256_div_pd(a2, b2));

    return _mm256_inserti128_si256(_mm256_castsi128_si256(ab1), ab2, 1);
}

// Simulated gather.  This is sometimes faster as it can run on other ports.
static inline __m256i _mm256_i32gather_epi32x(int *b, __m256i idx, int size) {
    int c[8] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i *)c, idx);
    return _mm256_set_epi32(b[c[7]], b[c[6]], b[c[5]], b[c[4]], b[c[3]], b[c[2]], b[c[1]], b[c[0]]);
}

unsigned char *rans_compress_O0_32x16(unsigned char *in, unsigned int in_size, unsigned char *out, unsigned int *out_size) {
    unsigned char *cp, *out_end;
    RansEncSymbol syms[256];
    RansState ransN[NX] __attribute__((aligned(32)));
    uint8_t* ptr;
    int F[256+MAGIC_RANS] = {0}, i, j, tab_size, rle, x, fsum = 0;
    int m = 0, M = 0;
    int tabsz;
    int bound = rans_compress_bound_32x16(in_size,0,&tabsz), z;

    //uint64_t sym_mf[256], sym_bfs[256];
    //uint32_t SD[256], SC[256], SB[256], SA[256];
    uint32_t SB[256], SA[256], SD[256], SC[256];

    if (!out) {
	*out_size = bound;
	out = (unsigned char*)malloc(*out_size);
    }
    if (!out || bound > *out_size) {
	return NULL;
    }

    ptr = out_end = out + bound;

    // Compute statistics
    hist8(in, in_size, F);
    double p = (double)TOTFREQ/(double)in_size;

    // Normalise so T[i] == TOTFREQ
    for (m = M = j = 0; j < 256; j++) {
	if (!F[j])
	    continue;

	if (m < F[j])
	    m = F[j], M = j;

	if ((F[j] = F[j]*p+0.499) == 0)
	    F[j] = 1;
	fsum += F[j];
    }

    fsum++; // not needed, but can't remove without removing assert x<TOTFREQ (in old code)
    int adjust = TOTFREQ - fsum;
    if (adjust > 0) {
	F[M] += adjust;
    } else if (adjust < 0) {
	if (F[M] > -adjust) {
	    F[M] += adjust;
	} else {
	    adjust += F[M]-1;
	    F[M] = 1;
	    for (j = 0; adjust && j < 256; j++) {
		if (F[j] < 2) continue;

		int d = F[j] > -adjust;
		int m = d ? adjust : 1-F[j];
		F[j]   += m;
		adjust -= m;
	    }
	}
    }

    //printf("F[%d]=%d\n", M, F[M]);
    assert(F[M]>0);

    // Encode statistics.
    cp = out+4;

    for (x = rle = j = 0; j < 256; j++) {
	if (F[j]) {
	    // j
	    if (rle) {
		rle--;
	    } else {
		*cp++ = j;
		if (!rle && j && F[j-1])  {
		    for(rle=j+1; rle<256 && F[rle]; rle++)
			;
		    rle -= j+1;
		    *cp++ = rle;
		}
		//fprintf(stderr, "%d: %d %d\n", j, rle, N[j]);
	    }
	    
	    // F[j]
	    if (F[j]<128) {
		*cp++ = F[j];
	    } else {
		*cp++ = 128 | (F[j]>>8);
		*cp++ = F[j]&0xff;
	    }
	    RansEncSymbolInit(&syms[j], x, F[j], TF_SHIFT);
	    //SB[j] = syms[j].x_max;
	    //SA[i] = syms[j].rcp_freq;
	    //SD[i] = (syms[j].cmpl_freq<<0) | (syms[j].rcp_shift<<16);
	    //SC[i] = syms[j].bias;
	    x += F[j];
	}
    }
    *cp++ = 0;

    //write(2, out+4, cp-(out+4));
    tab_size = cp-out;

    for (z = 0; z < NX; z++)
      RansEncInit(&ransN[z]);

    z = i = in_size&(NX-1);
    while (z-- > 0)
      RansEncPutSymbol(&ransN[z], &ptr, &syms[in[in_size-(i-z)]]);

    for (i = 0; i < 256; i++) {
	//sym_mf[i] = syms[i].rcp_freq | (((uint64_t)syms[i].x_max)<<32);
	SB[i] = syms[i].x_max;
	SA[i] = syms[i].rcp_freq;
	SD[i] = (syms[i].cmpl_freq<<0) | (syms[i].rcp_shift<<16);
	SC[i] = syms[i].bias;
    }

    uint16_t *ptr16 = (uint16_t *)ptr;

#define LOAD1(a,b) __m256i a##1 = _mm256_loadu_si256((__m256i *)&b[0]);
#define LOAD2(a,b) __m256i a##2 = _mm256_loadu_si256((__m256i *)&b[8]);
#define LOAD3(a,b) __m256i a##3 = _mm256_loadu_si256((__m256i *)&b[16]);
#define LOAD4(a,b) __m256i a##4 = _mm256_loadu_si256((__m256i *)&b[24]);
#define LOAD(a,b) LOAD1(a,b);LOAD2(a,b);LOAD3(a,b);LOAD4(a,b)

#define STORE1(a,b) _mm256_store_si256((__m256i *)&b[0],  a##1);
#define STORE2(a,b) _mm256_store_si256((__m256i *)&b[8],  a##2);
#define STORE3(a,b) _mm256_store_si256((__m256i *)&b[16], a##3);
#define STORE4(a,b) _mm256_store_si256((__m256i *)&b[24], a##4);
#define STORE(a,b) STORE1(a,b);STORE2(a,b);STORE3(a,b);STORE4(a,b)

    LOAD(Rv, ransN);

    for (i=(in_size &~(NX-1)); i>0; i-=NX) {
	uint8_t *c = &in[i-32];

// Set vs gather methods of loading data.
// Gather is faster, but can only schedule a few to run in parallel.
#define SET1(a,b) __m256i a##1 = _mm256_set_epi32(b[c[ 7]], b[c[ 6]], b[c[ 5]], b[c[ 4]], b[c[ 3]], b[c[ 2]], b[c[ 1]], b[c[ 0]])
#define SET2(a,b) __m256i a##2 = _mm256_set_epi32(b[c[15]], b[c[14]], b[c[13]], b[c[12]], b[c[11]], b[c[10]], b[c[ 9]], b[c[ 8]])
#define SET3(a,b) __m256i a##3 = _mm256_set_epi32(b[c[23]], b[c[22]], b[c[21]], b[c[20]], b[c[19]], b[c[18]], b[c[17]], b[c[16]])
#define SET4(a,b) __m256i a##4 = _mm256_set_epi32(b[c[31]], b[c[30]], b[c[29]], b[c[28]], b[c[27]], b[c[26]], b[c[25]], b[c[24]])
#define SET(a,b) SET1(a,b);SET2(a,b);SET3(a,b);SET4(a,b)

	// Renorm:
	// if (x > x_max) {*--ptr16 = x & 0xffff; x >>= 16;}
	SET(xmax, SB);
	__m256i cv1 = _mm256_cmpgt_epi32(Rv1, xmax1);
	__m256i cv2 = _mm256_cmpgt_epi32(Rv2, xmax2);
	__m256i cv3 = _mm256_cmpgt_epi32(Rv3, xmax3);
	__m256i cv4 = _mm256_cmpgt_epi32(Rv4, xmax4);

	// Store bottom 16-bits at ptr16
	unsigned int imask1 = _mm256_movemask_ps((__m256)cv1);
	unsigned int imask2 = _mm256_movemask_ps((__m256)cv2);
	unsigned int imask3 = _mm256_movemask_ps((__m256)cv3);
	unsigned int imask4 = _mm256_movemask_ps((__m256)cv4);

	__m256i idx1 = _mm256_loadu_si256((const __m256i*)permutec[imask1]);
	__m256i idx2 = _mm256_loadu_si256((const __m256i*)permutec[imask2]);
	__m256i idx3 = _mm256_loadu_si256((const __m256i*)permutec[imask3]);
	__m256i idx4 = _mm256_loadu_si256((const __m256i*)permutec[imask4]);

	// Permute; to gather together the rans states that need flushing
	__m256i V1 = _mm256_permutevar8x32_epi32(_mm256_and_si256(Rv1, cv1), idx1);
	__m256i V2 = _mm256_permutevar8x32_epi32(_mm256_and_si256(Rv2, cv2), idx2);
	__m256i V3 = _mm256_permutevar8x32_epi32(_mm256_and_si256(Rv3, cv3), idx3);
	__m256i V4 = _mm256_permutevar8x32_epi32(_mm256_and_si256(Rv4, cv4), idx4);
	
	// We only flush bottom 16 bits, to squash 32-bit states into 16 bit.
	V1 = _mm256_and_si256(V1, _mm256_set1_epi32(0xffff));
	V2 = _mm256_and_si256(V2, _mm256_set1_epi32(0xffff));
	V3 = _mm256_and_si256(V3, _mm256_set1_epi32(0xffff));
	V4 = _mm256_and_si256(V4, _mm256_set1_epi32(0xffff));
	__m256i V12 = _mm256_packus_epi32(V1, V2);
	__m256i V34 = _mm256_packus_epi32(V3, V4);

	// It's BAba order, want BbAa so shuffle.
	V12 = _mm256_permute4x64_epi64(V12, 0xd8);
	V34 = _mm256_permute4x64_epi64(V34, 0xd8);

	// Now we have bottom N 16-bit values in each V12/V34 to flush
	__m128i f =  _mm256_extractf128_si256(V34, 1);
	_mm_storeu_si128((__m128i *)(ptr16-8), f);
	ptr16 -= _mm_popcnt_u32(imask4);

	f =  _mm256_extractf128_si256(V34, 0);
	_mm_storeu_si128((__m128i *)(ptr16-8), f);
	ptr16 -= _mm_popcnt_u32(imask3);

	f =  _mm256_extractf128_si256(V12, 1);
	_mm_storeu_si128((__m128i *)(ptr16-8), f);
	ptr16 -= _mm_popcnt_u32(imask2);

	f =  _mm256_extractf128_si256(V12, 0);
	_mm_storeu_si128((__m128i *)(ptr16-8), f);
	ptr16 -= _mm_popcnt_u32(imask1);

	__m256i Rs;
	Rs = _mm256_srli_epi32(Rv1, 16); Rv1 = _mm256_blendv_epi8(Rv1, Rs, cv1);
	Rs = _mm256_srli_epi32(Rv2, 16); Rv2 = _mm256_blendv_epi8(Rv2, Rs, cv2);
	Rs = _mm256_srli_epi32(Rv3, 16); Rv3 = _mm256_blendv_epi8(Rv3, Rs, cv3);
	Rs = _mm256_srli_epi32(Rv4, 16); Rv4 = _mm256_blendv_epi8(Rv4, Rs, cv4);

	// Cannot trivially replace the multiply as mulhi_epu32 doesn't exist (only mullo).
	// However we can use _mm256_mul_epu32 twice to get 64bit results (half our lanes)
	// and shift/or to get the answer.
	//
	// (AVX512 allows us to hold it all in 64-bit lanes and use mullo_epi64
	// plus a shift.  KNC has mulhi_epi32, but not sure if this is available.)

	SET(rfv,   SA);

	rfv1 = _mm256_mulhi_epu32(Rv1, rfv1);
	rfv2 = _mm256_mulhi_epu32(Rv2, rfv2);
	rfv3 = _mm256_mulhi_epu32(Rv3, rfv3);
	rfv4 = _mm256_mulhi_epu32(Rv4, rfv4);

	SET(SDv,   SD);

	__m256i shiftv1 = _mm256_srli_epi32(SDv1, 16);
	__m256i shiftv2 = _mm256_srli_epi32(SDv2, 16);
	__m256i shiftv3 = _mm256_srli_epi32(SDv3, 16);
	__m256i shiftv4 = _mm256_srli_epi32(SDv4, 16);

	shiftv1 = _mm256_sub_epi32(shiftv1, _mm256_set1_epi32(32));
	shiftv2 = _mm256_sub_epi32(shiftv2, _mm256_set1_epi32(32));
	shiftv3 = _mm256_sub_epi32(shiftv3, _mm256_set1_epi32(32));
	shiftv4 = _mm256_sub_epi32(shiftv4, _mm256_set1_epi32(32));

	__m256i qv1 = _mm256_srlv_epi32(rfv1, shiftv1);
	__m256i qv2 = _mm256_srlv_epi32(rfv2, shiftv2);

	__m256i freqv1 = _mm256_and_si256(SDv1, _mm256_set1_epi32(0xffff));
	__m256i freqv2 = _mm256_and_si256(SDv2, _mm256_set1_epi32(0xffff));
	qv1 = _mm256_mullo_epi32(qv1, freqv1);
	qv2 = _mm256_mullo_epi32(qv2, freqv2);

	__m256i qv3 = _mm256_srlv_epi32(rfv3, shiftv3);
	__m256i qv4 = _mm256_srlv_epi32(rfv4, shiftv4);

	__m256i freqv3 = _mm256_and_si256(SDv3, _mm256_set1_epi32(0xffff));
	__m256i freqv4 = _mm256_and_si256(SDv4, _mm256_set1_epi32(0xffff));
	qv3 = _mm256_mullo_epi32(qv3, freqv3);
	qv4 = _mm256_mullo_epi32(qv4, freqv4);

	SET(biasv, SC);

	qv1 = _mm256_add_epi32(qv1, biasv1);
	qv2 = _mm256_add_epi32(qv2, biasv2);
	qv3 = _mm256_add_epi32(qv3, biasv3);
	qv4 = _mm256_add_epi32(qv4, biasv4);

	Rv1 = _mm256_add_epi32(Rv1, qv1);
	Rv2 = _mm256_add_epi32(Rv2, qv2);
	Rv3 = _mm256_add_epi32(Rv3, qv3);
	Rv4 = _mm256_add_epi32(Rv4, qv4);
    }

    STORE(Rv, ransN);

    ptr = (uint8_t *)ptr16;

    for (z = NX-1; z >= 0; z--)
      RansEncFlush(&ransN[z], &ptr);

    // Finalise block size and return it
    *out_size = (out_end - ptr) + tab_size;

    cp = out;
    *cp++ = (in_size>> 0) & 0xff;
    *cp++ = (in_size>> 8) & 0xff;
    *cp++ = (in_size>>16) & 0xff;
    *cp++ = (in_size>>24) & 0xff;

    memmove(out + tab_size, ptr, out_end-ptr);

    return out;
}

typedef struct {
    unsigned char R[TOTFREQ];
} ari_decoder;

unsigned char *rans_uncompress_O0_32x16(unsigned char *in, unsigned int in_size,
				       unsigned char *out, unsigned int *out_size) {
    /* Load in the static tables */
    unsigned char *cp = in + 4;
    int i, j, x, y, out_sz, rle;
    //uint16_t sfreq[TOTFREQ+32];
    //uint16_t sbase[TOTFREQ+32]; // faster to use 32-bit on clang
    uint8_t  ssym [TOTFREQ+64]; // faster to use 16-bit on clang

    uint32_t s3[TOTFREQ] __attribute__((aligned(32))); // For TF_SHIFT <= 12

    out_sz = ((in[0])<<0) | ((in[1])<<8) | ((in[2])<<16) | ((in[3])<<24);
    if (!out) {
	out = (unsigned char*)malloc(out_sz);
	*out_size = out_sz;
    }
    if (!out || out_sz > *out_size)
	return NULL;

    // Precompute reverse lookup of frequency.
    rle = x = y = 0;
    j = *cp++;
    do {
	int F, C;
	if ((F = *cp++) >= 128) {
	    F &= ~128;
	    F = ((F & 127) << 8) | *cp++;
	}
	C = x;

        for (y = 0; y < F; y++) {
	  ssym [y + C] = j;
	  //sfreq[y + C] = F;
	  //sbase[y + C] = y;
	  s3[y+C] = (((uint32_t)F)<<(TF_SHIFT+8))|(y<<8)|j;
        }
	x += F;

	if (!rle && j+1 == *cp) {
	    j = *cp++;
	    rle = *cp++;
	} else if (rle) {
	    rle--;
	    j++;
	} else {
	    j = *cp++;
	}
    } while(j);

    assert(x < TOTFREQ);

    int z;

    RansState R[NX] __attribute__((aligned(32)));
    for (z = 0; z < NX; z++)
      RansDecInit(&R[z], &cp);

    uint16_t *sp = (uint16_t *)cp;

    int out_end = (out_sz&~(NX-1));
    const uint32_t mask = (1u << TF_SHIFT)-1;

    __m256i maskv  = _mm256_set1_epi32(mask); // set mask in all lanes
    LOAD(Rv, R);

    for (i=0; i < out_end; i+=NX) {
	//for (z = 0; z < NX; z++)
	//  m[z] = R[z] & mask;
	__m256i masked1 = _mm256_and_si256(Rv1, maskv);
	__m256i masked2 = _mm256_and_si256(Rv2, maskv);

	//  S[z] = s3[m[z]];
	__m256i Sv1 = _mm256_i32gather_epi32x((int *)s3, masked1, sizeof(*s3));
	__m256i Sv2 = _mm256_i32gather_epi32x((int *)s3, masked2, sizeof(*s3));

	//  f[z] = S[z]>>(TF_SHIFT+8);
	__m256i fv1 = _mm256_srli_epi32(Sv1, TF_SHIFT+8);
	__m256i fv2 = _mm256_srli_epi32(Sv2, TF_SHIFT+8);

	//  b[z] = (S[z]>>8) & mask;
	__m256i bv1 = _mm256_and_si256(_mm256_srli_epi32(Sv1, 8), maskv);
	__m256i bv2 = _mm256_and_si256(_mm256_srli_epi32(Sv2, 8), maskv);

	//  s[z] = S[z] & 0xff;
	__m256i sv1 = _mm256_and_si256(Sv1, _mm256_set1_epi32(0xff));
	__m256i sv2 = _mm256_and_si256(Sv2, _mm256_set1_epi32(0xff));

	//  R[z] = f[z] * (R[z] >> TF_SHIFT) + b[z];
	Rv1 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(Rv1,TF_SHIFT),fv1),bv1);
	Rv2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(Rv2,TF_SHIFT),fv2),bv2);

	// Tricky one:  out[i+z] = s[z];
	//             ---h---g ---f---e  ---d---c ---b---a
	//             ---p---o ---n---m  ---l---k ---j---i
	// packs_epi32 -p-o-n-m -h-g-f-e  -l-k-j-i -d-c-b-a
	// permute4x64 -p-o-n-m -l-k-j-i  -h-g-f-e -d-c-b-a
	// packs_epi16 ponmlkji ponmlkji  hgfedcba hgfedcba
	sv1 = _mm256_packus_epi32(sv1, sv2);
	sv1 = _mm256_permute4x64_epi64(sv1, 0xd8);
	__m256i Vv1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
	sv1 = _mm256_packus_epi16(sv1, sv1);

	// c =  R[z] < RANS_BYTE_L;
	__m256i renorm_mask1 = _mm256_xor_si256(Rv1, _mm256_set1_epi32(0x80000000));
	__m256i renorm_mask2 = _mm256_xor_si256(Rv2, _mm256_set1_epi32(0x80000000));
	renorm_mask1 = _mm256_cmpgt_epi32(_mm256_set1_epi32(RANS_BYTE_L-0x80000000), renorm_mask1);
	renorm_mask2 = _mm256_cmpgt_epi32(_mm256_set1_epi32(RANS_BYTE_L-0x80000000), renorm_mask2);
	
	// y = (R[z] << 16) | V[z];
	unsigned int imask1 = _mm256_movemask_ps((__m256)renorm_mask1);
	__m256i idx1 = _mm256_loadu_si256((const __m256i*)permute[imask1]);
	__m256i Yv1 = _mm256_slli_epi32(Rv1, 16);
	Vv1 = _mm256_permutevar8x32_epi32(Vv1, idx1);
	__m256i Yv2 = _mm256_slli_epi32(Rv2, 16);

	// Shuffle the renorm values to correct lanes and incr sp pointer
	unsigned int imask2 = _mm256_movemask_ps((__m256)renorm_mask2);
	sp += _mm_popcnt_u32(imask1);

	__m256i idx2 = _mm256_loadu_si256((const __m256i*)permute[imask2]);
	__m256i Vv2 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
	sp += _mm_popcnt_u32(imask2);

	Yv1 = _mm256_or_si256(Yv1, Vv1);
	Vv2 = _mm256_permutevar8x32_epi32(Vv2, idx2);
	Yv2 = _mm256_or_si256(Yv2, Vv2);

	// R[z] = c ? Y[z] : R[z];
	Rv1 = _mm256_blendv_epi8(Rv1, Yv1, renorm_mask1);
	Rv2 = _mm256_blendv_epi8(Rv2, Yv2, renorm_mask2);

	// ------------------------------------------------------------

	//  m[z] = R[z] & mask;
	//  S[z] = s3[m[z]];
	__m256i masked3 = _mm256_and_si256(Rv3, maskv);
	__m256i Sv3 = _mm256_i32gather_epi32x((int *)s3, masked3, sizeof(*s3));

	*(uint64_t *)&out[i+0] = _mm256_extract_epi64(sv1, 0);
	*(uint64_t *)&out[i+8] = _mm256_extract_epi64(sv1, 2);

	__m256i masked4 = _mm256_and_si256(Rv4, maskv);
	__m256i Sv4 = _mm256_i32gather_epi32x((int *)s3, masked4, sizeof(*s3));

	//  f[z] = S[z]>>(TF_SHIFT+8);
	__m256i fv3 = _mm256_srli_epi32(Sv3, TF_SHIFT+8);
	__m256i fv4 = _mm256_srli_epi32(Sv4, TF_SHIFT+8);

	//  b[z] = (S[z]>>8) & mask;
	__m256i bv3 = _mm256_and_si256(_mm256_srli_epi32(Sv3, 8), maskv);
	__m256i bv4 = _mm256_and_si256(_mm256_srli_epi32(Sv4, 8), maskv);

	//  s[z] = S[z] & 0xff;
	__m256i sv3 = _mm256_and_si256(Sv3, _mm256_set1_epi32(0xff));
	__m256i sv4 = _mm256_and_si256(Sv4, _mm256_set1_epi32(0xff));

	//  R[z] = f[z] * (R[z] >> TF_SHIFT) + b[z];
	Rv3 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(Rv3,TF_SHIFT),fv3),bv3);
	Rv4 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(Rv4,TF_SHIFT),fv4),bv4);

	// Tricky one:  out[i+z] = s[z];
	//             ---h---g ---f---e  ---d---c ---b---a
	//             ---p---o ---n---m  ---l---k ---j---i
	// packs_epi32 -p-o-n-m -h-g-f-e  -l-k-j-i -d-c-b-a
	// permute4x64 -p-o-n-m -l-k-j-i  -h-g-f-e -d-c-b-a
	// packs_epi16 ponmlkji ponmlkji  hgfedcba hgfedcba
	sv3 = _mm256_packus_epi32(sv3, sv4);
	sv3 = _mm256_permute4x64_epi64(sv3, 0xd8);
	__m256i renorm_mask3 = _mm256_xor_si256(Rv3, _mm256_set1_epi32(0x80000000));
	__m256i renorm_mask4 = _mm256_xor_si256(Rv4, _mm256_set1_epi32(0x80000000));
	sv3 = _mm256_packus_epi16(sv3, sv3);
	// c =  R[z] < RANS_BYTE_L;

	renorm_mask3 = _mm256_cmpgt_epi32(_mm256_set1_epi32(RANS_BYTE_L-0x80000000), renorm_mask3);
	renorm_mask4 = _mm256_cmpgt_epi32(_mm256_set1_epi32(RANS_BYTE_L-0x80000000), renorm_mask4);
	
	*(uint64_t *)&out[i+16] = _mm256_extract_epi64(sv3, 0);
	*(uint64_t *)&out[i+24] = _mm256_extract_epi64(sv3, 2);

	// y = (R[z] << 16) | V[z];
	__m256i Vv3 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
	__m256i Yv3 = _mm256_slli_epi32(Rv3, 16);
	unsigned int imask3 = _mm256_movemask_ps((__m256)renorm_mask3);
	__m256i idx3 = _mm256_loadu_si256((const __m256i*)permute[imask3]);

	// Shuffle the renorm values to correct lanes and incr sp pointer
	Vv3 = _mm256_permutevar8x32_epi32(Vv3, idx3);
	__m256i Yv4 = _mm256_slli_epi32(Rv4, 16);
	unsigned int imask4 = _mm256_movemask_ps((__m256)renorm_mask4);
	sp += _mm_popcnt_u32(imask3);

	__m256i idx4 = _mm256_loadu_si256((const __m256i*)permute[imask4]);
	__m256i Vv4 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));

	//Vv = _mm256_and_si256(Vv, renorm_mask);  (blend does the AND anyway)
	Yv3 = _mm256_or_si256(Yv3, Vv3);
	Vv4 = _mm256_permutevar8x32_epi32(Vv4, idx4);
	Yv4 = _mm256_or_si256(Yv4, Vv4);

	sp += _mm_popcnt_u32(imask4);

	// R[z] = c ? Y[z] : R[z];
	Rv3 = _mm256_blendv_epi8(Rv3, Yv3, renorm_mask3);
	Rv4 = _mm256_blendv_epi8(Rv4, Yv4, renorm_mask4);
    }

    STORE(Rv, R);
    //_mm256_store_si256((__m256i *)&R[0], Rv1);
    //_mm256_store_si256((__m256i *)&R[8], Rv2);
    //_mm256_store_si256((__m256i *)&R[16], Rv3);
    //_mm256_store_si256((__m256i *)&R[24], Rv4);

//#pragma omp simd
//	for (z = 0; z < NX; z++) {
//	  uint32_t m = R[z] & mask;
//	  R[z] = sfreq[m] * (R[z] >> TF_SHIFT) + sbase[m];
//	  out[i+z] = ssym[m];
//	  uint32_t c = R[z] < RANS_BYTE_L;  // NX16=>166MB/s
//	  uint32_t y = (R[z] << 16) | *spN[z];
//	  spN[z] += c ? 1 : 0;
//	  R[z]    = c ? y : R[z];
//
//	}
//    }

    for (z = out_sz & (NX-1); z-- > 0; )
      out[out_end + z] = ssym[R[z] & mask];

    *out_size = out_sz;
    return out;
}

static void hist1_4(unsigned char *in, unsigned int in_size,
		    int F0[256][256], int *T0) {
    int T1[256+MAGIC_RANS] = {0}, T2[256+MAGIC_RANS] = {0}, T3[256+MAGIC_RANS] = {0};
    unsigned int idiv4 = in_size/4;
    int i;
    unsigned char c0, c1, c2, c3;

    unsigned char *in0 = in + 0;
    unsigned char *in1 = in + idiv4;
    unsigned char *in2 = in + idiv4*2;
    unsigned char *in3 = in + idiv4*3;

    unsigned char last_0 = 0, last_1 = in1[-1], last_2 = in2[-1], last_3 = in3[-1];
    //unsigned char last_0 = 0, last_1 = 0, last_2 = 0, last_3 = 0;

    unsigned char *in0_end = in1;

    while (in0 < in0_end) {
	F0[last_0][c0 = *in0++]++;
	T0[last_0]++;
	last_0 = c0;

	F0[last_1][c1 = *in1++]++;
	T1[last_1]++;
	last_1 = c1;

	F0[last_2][c2 = *in2++]++;
	T2[last_2]++;
	last_2 = c2;

	F0[last_3][c3 = *in3++]++;
	T3[last_3]++;
	last_3 = c3;
    }

    while (in3 < in + in_size) {
	F0[last_3][c3 = *in3++]++;
	T3[last_3]++;
	last_3 = c3;
    }

    for (i = 0; i < 256; i++) {
	T0[i]+=T1[i]+T2[i]+T3[i];
    }
}

unsigned char *rans_compress_O1_32x16(unsigned char *in, unsigned int in_size,
				     unsigned char *out, unsigned int *out_size) {
    unsigned char *cp, *out_end;
    unsigned int tab_size, rle_i, rle_j;
    RansEncSymbol syms[256][256];
    int bound = rans_compress_bound_32x16(in_size,1, NULL), z;
    RansState ransN[NX] __attribute__((aligned(32)));

    if (!out) {
	*out_size = bound;
	out = (unsigned char*)malloc(*out_size);
    }
    if (!out || bound > *out_size)
	return NULL;

    out_end = out + bound;
    cp = out+4;

    int F[256][256] = {{0}}, T[256+MAGIC_RANS] = {0}, i, j;

    //memset(F, 0, 256*256*sizeof(int));
    //memset(T, 0, 256*sizeof(int));

    hist1_4(in, in_size, F, T);

    for (z = 1; z < NX; z++)
	F[0][in[z*(in_size/NX)]]++;
    T[0]+=NX-1;

    // Normalise so T[i] == TOTFREQ
    for (rle_i = i = 0; i < 256; i++) {
	int t2, m, M;
	unsigned int x;

	if (T[i] == 0)
	    continue;

	//uint64_t p = (TOTFREQ * TOTFREQ) / t;
	double p = ((double)TOTFREQ_O1)/T[i];

	for (t2 = m = M = j = 0; j < 256; j++) {
	    if (!F[i][j])
		continue;

	    if (m < F[i][j])
		m = F[i][j], M = j;

	    if ((F[i][j] *= p) <= 0)
	        F[i][j] = 1;
	    t2 += F[i][j];
	}

	//t2++;

	int adjust = TOTFREQ_O1-t2;
	if (adjust > 0) {
	    // Boost most common
	    F[i][M] += adjust;
	} else if (adjust < 0) {
	    // Reduce highest and distribute remainder
	    if (F[i][M] > -adjust) {
		F[i][M] += adjust;
	    } else {
		adjust += F[i][M]-1;
		F[i][M] = 1;

		for (j = 0; adjust && j < 256; j++) {
		    if (F[i][j] < 2) continue;

		    int d = F[i][j] > -adjust;
		    int m = d ? adjust : 1-F[i][j];
		    F[i][j]   += m;
		    adjust -= m;
		}
	    }
	}

	// Store frequency table
	// i
	if (rle_i) {
	    rle_i--;
	} else {
	    *cp++ = i;
	    // FIXME: could use order-0 statistics to observe which alphabet
	    // symbols are present and base RLE on that ordering instead.
	    if (i && T[i-1]) {
		for(rle_i=i+1; rle_i<256 && T[rle_i]; rle_i++)
		    ;
		rle_i -= i+1;
		*cp++ = rle_i;
	    }
	}

	int *F_i_ = F[i];
	x = 0;
	rle_j = 0;
	for (j = 0; j < 256; j++) {
	    if (F_i_[j]) {
		//fprintf(stderr, "F[%d][%d]=%d, x=%d\n", i, j, F_i_[j], x);

		// j
		if (rle_j) {
		    rle_j--;
		} else {
		    *cp++ = j;
		    if (!rle_j && j && F_i_[j-1]) {
			for(rle_j=j+1; rle_j<256 && F_i_[rle_j]; rle_j++)
			    ;
			rle_j -= j+1;
			*cp++ = rle_j;
		    }
		}

		// F_i_[j]
		if (F_i_[j]<128) {
 		    *cp++ = F_i_[j];
		} else {
		    *cp++ = 128 | (F_i_[j]>>8);
		    *cp++ = F_i_[j]&0xff;
		}

		RansEncSymbolInit(&syms[i][j], x, F_i_[j], TF_SHIFT_O1);
		x += F_i_[j];
	    }
	}
	*cp++ = 0;
    }
    *cp++ = 0;

    //write(2, out+4, cp-(out+4));
    tab_size = cp - out;
    assert(tab_size < 257*257*3);
    
    for (z = 0; z < NX; z++)
      RansEncInit(&ransN[z]);

    uint8_t* ptr = out_end;

    int isz4 = in_size/NX;
    int iN[NX];
    for (z = 0; z < NX; z++)
	iN[z] = (z+1)*isz4-2;

    unsigned char lN[NX];
    for (z = 0; z < NX; z++)
	lN[z] = in[iN[z]+1];

    // Deal with the remainder
    z = NX-1;
    lN[z] = in[in_size-1];
    for (iN[z] = in_size-2; iN[z] > NX*isz4-2; iN[z]--) {
	unsigned char c = in[iN[z]];
	RansEncPutSymbol(&ransN[z], &ptr, &syms[c][lN[z]]);
	lN[z] = c;
    }

    uint16_t *ptr16 = (uint16_t *)ptr;

    LOAD(Rv, ransN);

    for (; iN[0] >= 0; ) {
	uint32_t c[NX];

	// Gather all the symbol values together in adjacent arrays.
	// Better to just use raw set?
	RansEncSymbol *sN[NX];
	for (z = 0; z < NX; z++)
	    sN[z] = &syms[c[z] = in[iN[z]]][lN[z]];

#define SET1x(a,b,x) __m256i a##1 = _mm256_set_epi32(b[ 7]->x, b[ 6]->x, b[ 5]->x, b[ 4]->x, b[ 3]->x, b[ 2]->x, b[ 1]->x, b[ 0]->x)
#define SET2x(a,b,x) __m256i a##2 = _mm256_set_epi32(b[15]->x, b[14]->x, b[13]->x, b[12]->x, b[11]->x, b[10]->x, b[ 9]->x, b[ 8]->x)
#define SET3x(a,b,x) __m256i a##3 = _mm256_set_epi32(b[23]->x, b[22]->x, b[21]->x, b[20]->x, b[19]->x, b[18]->x, b[17]->x, b[16]->x)
#define SET4x(a,b,x) __m256i a##4 = _mm256_set_epi32(b[31]->x, b[30]->x, b[29]->x, b[28]->x, b[27]->x, b[26]->x, b[25]->x, b[24]->x)
#define SETx(a,b,x) SET1x(a,b,x);SET2x(a,b,x);SET3x(a,b,x);SET4x(a,b,x)

	// ------------------------------------------------------------
	//	for (z = NX-1; z >= 0; z--) {
	//	    if (ransN[z] >= x_max[z]) {
	//		*--ptr16 = ransN[z] & 0xffff;
	//		ransN[z] >>= 16;
	//	    }
	//	}
	//LOAD(xmax,x_max);
	SETx(xmax, sN, x_max);
        __m256i cv1 = _mm256_cmpgt_epi32(Rv1, xmax1);
        __m256i cv2 = _mm256_cmpgt_epi32(Rv2, xmax2);
        __m256i cv3 = _mm256_cmpgt_epi32(Rv3, xmax3);
        __m256i cv4 = _mm256_cmpgt_epi32(Rv4, xmax4);

        // Store bottom 16-bits at ptr16                                                     
        //                                                                                   
        // for (z = NX-1; z >= 0; z--) {                                                     
        //     if (cond[z]) *--ptr16 = (uint16_t )(ransN[z] & 0xffff);                       
        // }                                                                                 
        unsigned int imask1 = _mm256_movemask_ps((__m256)cv1);
        unsigned int imask2 = _mm256_movemask_ps((__m256)cv2);
        unsigned int imask3 = _mm256_movemask_ps((__m256)cv3);
        unsigned int imask4 = _mm256_movemask_ps((__m256)cv4);

        __m256i idx1 = _mm256_loadu_si256((const __m256i*)permutec[imask1]);
        __m256i idx2 = _mm256_loadu_si256((const __m256i*)permutec[imask2]);
        __m256i idx3 = _mm256_loadu_si256((const __m256i*)permutec[imask3]);
        __m256i idx4 = _mm256_loadu_si256((const __m256i*)permutec[imask4]);

        // Permute; to gather together the rans states that need flushing                    
        __m256i V1 = _mm256_permutevar8x32_epi32(_mm256_and_si256(Rv1, cv1), idx1);
        __m256i V2 = _mm256_permutevar8x32_epi32(_mm256_and_si256(Rv2, cv2), idx2);
        __m256i V3 = _mm256_permutevar8x32_epi32(_mm256_and_si256(Rv3, cv3), idx3);
        __m256i V4 = _mm256_permutevar8x32_epi32(_mm256_and_si256(Rv4, cv4), idx4);

        // We only flush bottom 16 bits, to squash 32-bit states into 16 bit.                
        V1 = _mm256_and_si256(V1, _mm256_set1_epi32(0xffff));
        V2 = _mm256_and_si256(V2, _mm256_set1_epi32(0xffff));
        V3 = _mm256_and_si256(V3, _mm256_set1_epi32(0xffff));
        V4 = _mm256_and_si256(V4, _mm256_set1_epi32(0xffff));
        __m256i V12 = _mm256_packus_epi32(V1, V2);
        __m256i V34 = _mm256_packus_epi32(V3, V4);

        // It's BAba order, want BbAa so shuffle.                                            
        V12 = _mm256_permute4x64_epi64(V12, 0xd8);
        V34 = _mm256_permute4x64_epi64(V34, 0xd8);
        // Now we have bottom N 16-bit values in each V12/V34 to flush                       
        __m128i f =  _mm256_extractf128_si256(V34, 1);
        _mm_storeu_si128((__m128i *)(ptr16-8), f);
        ptr16 -= _mm_popcnt_u32(imask4);

        f =  _mm256_extractf128_si256(V34, 0);
        _mm_storeu_si128((__m128i *)(ptr16-8), f);
        ptr16 -= _mm_popcnt_u32(imask3);

        f =  _mm256_extractf128_si256(V12, 1);
        _mm_storeu_si128((__m128i *)(ptr16-8), f);
        ptr16 -= _mm_popcnt_u32(imask2);

        f =  _mm256_extractf128_si256(V12, 0);
        _mm_storeu_si128((__m128i *)(ptr16-8), f);
        ptr16 -= _mm_popcnt_u32(imask1);

        __m256i Rs;
        Rs = _mm256_srli_epi32(Rv1, 16); Rv1 = _mm256_blendv_epi8(Rv1, Rs, cv1);
        Rs = _mm256_srli_epi32(Rv2, 16); Rv2 = _mm256_blendv_epi8(Rv2, Rs, cv2);
        Rs = _mm256_srli_epi32(Rv3, 16); Rv3 = _mm256_blendv_epi8(Rv3, Rs, cv3);
        Rs = _mm256_srli_epi32(Rv4, 16); Rv4 = _mm256_blendv_epi8(Rv4, Rs, cv4);

	// ------------------------------------------------------------
	// uint32_t q = (uint32_t) (((uint64_t)ransN[z] * rcp_freq[z]) >> rcp_shift[z]);
	// ransN[z] = ransN[z] + bias[z] + q * cmpl_freq[z];

        // Cannot trivially replace the multiply as mulhi_epu32 doesn't exist (only mullo).  
        // However we can use _mm256_mul_epu32 twice to get 64bit results (half our lanes)   
        // and shift/or to get the answer.                                                   
        //                                                                                   
        // (AVX512 allows us to hold it all in 64-bit lanes and use mullo_epi64              
        // plus a shift.  KNC has mulhi_epi32, but not sure if this is available.)           

	SETx(rfv,   sN, rcp_freq);

        rfv1 = _mm256_mulhi_epu32(Rv1, rfv1);
        rfv2 = _mm256_mulhi_epu32(Rv2, rfv2);
        rfv3 = _mm256_mulhi_epu32(Rv3, rfv3);
        rfv4 = _mm256_mulhi_epu32(Rv4, rfv4);

	SETx(SDv,   sN, SD);

        __m256i shiftv1 = _mm256_srli_epi32(SDv1, 16);
        __m256i shiftv2 = _mm256_srli_epi32(SDv2, 16);
        __m256i shiftv3 = _mm256_srli_epi32(SDv3, 16);
        __m256i shiftv4 = _mm256_srli_epi32(SDv4, 16);

        shiftv1 = _mm256_sub_epi32(shiftv1, _mm256_set1_epi32(32));
        shiftv2 = _mm256_sub_epi32(shiftv2, _mm256_set1_epi32(32));
        shiftv3 = _mm256_sub_epi32(shiftv3, _mm256_set1_epi32(32));
        shiftv4 = _mm256_sub_epi32(shiftv4, _mm256_set1_epi32(32));

        __m256i qv1 = _mm256_srlv_epi32(rfv1, shiftv1);
        __m256i qv2 = _mm256_srlv_epi32(rfv2, shiftv2);

        __m256i freqv1 = _mm256_and_si256(SDv1, _mm256_set1_epi32(0xffff));
        __m256i freqv2 = _mm256_and_si256(SDv2, _mm256_set1_epi32(0xffff));
        qv1 = _mm256_mullo_epi32(qv1, freqv1);
        qv2 = _mm256_mullo_epi32(qv2, freqv2);

        __m256i qv3 = _mm256_srlv_epi32(rfv3, shiftv3);
        __m256i qv4 = _mm256_srlv_epi32(rfv4, shiftv4);

        __m256i freqv3 = _mm256_and_si256(SDv3, _mm256_set1_epi32(0xffff));
        __m256i freqv4 = _mm256_and_si256(SDv4, _mm256_set1_epi32(0xffff));
        qv3 = _mm256_mullo_epi32(qv3, freqv3);
        qv4 = _mm256_mullo_epi32(qv4, freqv4);

	SETx(biasv, sN, bias);

        qv1 = _mm256_add_epi32(qv1, biasv1);
	qv2 = _mm256_add_epi32(qv2, biasv2);
        qv3 = _mm256_add_epi32(qv3, biasv3);
        qv4 = _mm256_add_epi32(qv4, biasv4);

	Rv1 = _mm256_add_epi32(Rv1, qv1);
        Rv2 = _mm256_add_epi32(Rv2, qv2);
        Rv3 = _mm256_add_epi32(Rv3, qv3);
        Rv4 = _mm256_add_epi32(Rv4, qv4);

	for (z = 0; z < NX; z++) {
	    //uint32_t q = (uint32_t) (((uint64_t)ransN[z] * rcp_freq[z]) >> rcp_shift[z]);
	    //ransN[z] = ransN[z] + bias[z] + q * cmpl_freq[z];
	    
	    lN[z] = c[z];
	    iN[z]--;
	}
    }

    STORE(Rv, ransN);

    ptr = (uint8_t *)ptr16;

    for (z = NX-1; z>=0; z--)
	RansEncPutSymbol(&ransN[z], &ptr, &syms[0][lN[z]]);

    for (z = NX-1; z>=0; z--)
	RansEncFlush(&ransN[z], &ptr);

    *out_size = (out_end - ptr) + tab_size;

    cp = out;
    *cp++ = (in_size>> 0) & 0xff;
    *cp++ = (in_size>> 8) & 0xff;
    *cp++ = (in_size>>16) & 0xff;
    *cp++ = (in_size>>24) & 0xff;

    memmove(out + tab_size, ptr, out_end-ptr);

    return out;
}

unsigned char *rans_uncompress_O1_32x16(unsigned char *in, unsigned int in_size,
				       unsigned char *out, unsigned int *out_size) {
    /* Load in the static tables */
    unsigned char *cp = in + 4;
    int i, j = -999, x, y, out_sz, rle_i, rle_j;
    uint32_t s3[256][TOTFREQ_O1] __attribute__((aligned(32)));
    
    //memset(D, 0, 256*sizeof(*D));

    out_sz = ((in[0])<<0) | ((in[1])<<8) | ((in[2])<<16) | ((in[3])<<24);
    if (!out) {
	out = (unsigned char*)malloc(out_sz);
	*out_size = out_sz;
    }
    if (!out || out_sz > *out_size)
	return NULL;

    //fprintf(stderr, "out_sz=%d\n", out_sz);

    //i = *cp++;
    rle_i = 0;
    i = *cp++;
    do {
	rle_j = x = y = 0;
	j = *cp++;
	do {
	    int F, C;
	    if ((F = *cp++) >= 128) {
		F &= ~128;
		F = ((F & 127) << 8) | *cp++;
	    }
	    C = x;

	    //fprintf(stderr, "i=%d j=%d F=%d C=%d\n", i, j, D[i].F[j], D[i].C[j]);

	    if (!F)
		F = TOTFREQ_O1;

	    for (y = 0; y < F; y++) {
		s3[i][y+C] = (((uint32_t)F)<<(TF_SHIFT_O1+8))|(y<<8)|j;
	    }

	    x += F;
	    assert(x <= TOTFREQ_O1);

	    if (!rle_j && j+1 == *cp) {
		j = *cp++;
		rle_j = *cp++;
	    } else if (rle_j) {
		rle_j--;
		j++;
	    } else {
		j = *cp++;
	    }
	} while(j);

	if (!rle_i && i+1 == *cp) {
	    i = *cp++;
	    rle_i = *cp++;
	} else if (rle_i) {
	    rle_i--;
	    i++;
	} else {
	    i = *cp++;
	}
    } while (i);

    // Precompute reverse lookup of frequency.

    RansState ransN[NX] __attribute__((aligned(32)));
    int z;
    uint8_t *ptr = cp;
    for (z = 0; z < NX; z++)
	RansDecInit(&ransN[z], &ptr);

    int isz4 = out_sz/NX;
    int lN[NX] = {0}, iN[NX];
    for (z = 0; z < NX; z++)
	iN[z] = z*isz4;

    uint16_t *sp = (uint16_t *)ptr;
    const uint32_t mask = (1u << TF_SHIFT_O1)-1;
    __m256i maskv  = _mm256_set1_epi32(mask);
    LOAD(Rv, ransN);
    LOAD(Lv, lN);

    union {
	unsigned char tbuf[32][32];
	uint64_t tbuf64[32][4];
    } u;
    int tidx = 0;

    for (; iN[0] < isz4; ) {
	// m[z] = ransN[z] & mask;
	__m256i masked1 = _mm256_and_si256(Rv1, maskv);
	__m256i masked2 = _mm256_and_si256(Rv2, maskv);

	//  S[z] = s3[lN[z]][m[z]];
	Lv1 = _mm256_slli_epi32(Lv1, TF_SHIFT_O1);
	masked1 = _mm256_add_epi32(masked1, Lv1);

	Lv2 = _mm256_slli_epi32(Lv2, TF_SHIFT_O1);
	masked2 = _mm256_add_epi32(masked2, Lv2);

	__m256i masked3 = _mm256_and_si256(Rv3, maskv);
	__m256i masked4 = _mm256_and_si256(Rv4, maskv);

	Lv3 = _mm256_slli_epi32(Lv3, TF_SHIFT_O1);
	masked3 = _mm256_add_epi32(masked3, Lv3);

	Lv4 = _mm256_slli_epi32(Lv4, TF_SHIFT_O1);
	masked4 = _mm256_add_epi32(masked4, Lv4);

	__m256i Sv1 = _mm256_i32gather_epi32x((int *)&s3[0][0], masked1, sizeof(s3[0][0]));
	__m256i Sv2 = _mm256_i32gather_epi32x((int *)&s3[0][0], masked2, sizeof(s3[0][0]));

	//  f[z] = S[z]>>(TF_SHIFT_O1+8);
	__m256i fv1 = _mm256_srli_epi32(Sv1, TF_SHIFT_O1+8);
	__m256i fv2 = _mm256_srli_epi32(Sv2, TF_SHIFT_O1+8);

	__m256i Sv3 = _mm256_i32gather_epi32x((int *)&s3[0][0], masked3, sizeof(s3[0][0]));
	__m256i Sv4 = _mm256_i32gather_epi32x((int *)&s3[0][0], masked4, sizeof(s3[0][0]));

	__m256i fv3 = _mm256_srli_epi32(Sv3, TF_SHIFT_O1+8);
	__m256i fv4 = _mm256_srli_epi32(Sv4, TF_SHIFT_O1+8);

	//  b[z] = (S[z]>>8) & mask;
	__m256i bv1 = _mm256_and_si256(_mm256_srli_epi32(Sv1, 8), maskv);
	__m256i bv2 = _mm256_and_si256(_mm256_srli_epi32(Sv2, 8), maskv);
	__m256i bv3 = _mm256_and_si256(_mm256_srli_epi32(Sv3, 8), maskv);
	__m256i bv4 = _mm256_and_si256(_mm256_srli_epi32(Sv4, 8), maskv);

	//  s[z] = S[z] & 0xff;
	__m256i sv1 = _mm256_and_si256(Sv1, _mm256_set1_epi32(0xff));
	__m256i sv2 = _mm256_and_si256(Sv2, _mm256_set1_epi32(0xff));
	__m256i sv3 = _mm256_and_si256(Sv3, _mm256_set1_epi32(0xff));
	__m256i sv4 = _mm256_and_si256(Sv4, _mm256_set1_epi32(0xff));

	//  R[z] = f[z] * (R[z] >> TF_SHIFT_O1) + b[z];
	Rv1 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(Rv1,TF_SHIFT_O1),fv1),bv1);
	Rv2 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(Rv2,TF_SHIFT_O1),fv2),bv2);


	//for (z = 0; z < NX; z++) lN[z] = c[z];
	Lv1 = sv1;
	Lv2 = sv2;

	sv1 = _mm256_packus_epi32(sv1, sv2);
	sv1 = _mm256_permute4x64_epi64(sv1, 0xd8);

	// Start loading next batch of normalised states
	__m256i Vv1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));

	sv1 = _mm256_packus_epi16(sv1, sv1);

	// out[iN[z]] = c[z];  // simulate scatter
	// RansDecRenorm(&ransN[z], &ptr);	
	__m256i renorm_mask1 = _mm256_xor_si256(Rv1, _mm256_set1_epi32(0x80000000));
	__m256i renorm_mask2 = _mm256_xor_si256(Rv2, _mm256_set1_epi32(0x80000000));

	renorm_mask1 = _mm256_cmpgt_epi32(_mm256_set1_epi32(RANS_BYTE_L-0x80000000), renorm_mask1);
	renorm_mask2 = _mm256_cmpgt_epi32(_mm256_set1_epi32(RANS_BYTE_L-0x80000000), renorm_mask2);

	unsigned int imask1 = _mm256_movemask_ps((__m256)renorm_mask1);
	__m256i idx1 = _mm256_loadu_si256((const __m256i*)permute[imask1]);
	__m256i Yv1 = _mm256_slli_epi32(Rv1, 16);
	__m256i Yv2 = _mm256_slli_epi32(Rv2, 16);

	unsigned int imask2 = _mm256_movemask_ps((__m256)renorm_mask2);
	Vv1 = _mm256_permutevar8x32_epi32(Vv1, idx1);
	sp += _mm_popcnt_u32(imask1);

	u.tbuf64[tidx][0] = _mm256_extract_epi64(sv1, 0);
	u.tbuf64[tidx][1] = _mm256_extract_epi64(sv1, 2);

	__m256i idx2 = _mm256_loadu_si256((const __m256i*)permute[imask2]);
	__m256i Vv2 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
	sp += _mm_popcnt_u32(imask2);
	Vv2 = _mm256_permutevar8x32_epi32(Vv2, idx2);

	//Vv = _mm256_and_si256(Vv, renorm_mask);  (blend does the AND anyway)
	Yv1 = _mm256_or_si256(Yv1, Vv1);
	Yv2 = _mm256_or_si256(Yv2, Vv2);

	Rv1 = _mm256_blendv_epi8(Rv1, Yv1, renorm_mask1);
	Rv2 = _mm256_blendv_epi8(Rv2, Yv2, renorm_mask2);

	/////////////////////////////////////////////////////////////////////

	// Start loading next batch of normalised states
	__m256i Vv3 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));

	//  R[z] = f[z] * (R[z] >> TF_SHIFT_O1) + b[z];
	Rv3 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(Rv3,TF_SHIFT_O1),fv3),bv3);
	Rv4 = _mm256_add_epi32(_mm256_mullo_epi32(_mm256_srli_epi32(Rv4,TF_SHIFT_O1),fv4),bv4);

	//for (z = 0; z < NX; z++) lN[z] = c[z];
	Lv3 = sv3;
	Lv4 = sv4;

	// out[iN[z]] = c[z];  // simulate scatter
	// RansDecRenorm(&ransN[z], &ptr);	
	__m256i renorm_mask3 = _mm256_xor_si256(Rv3, _mm256_set1_epi32(0x80000000));
	__m256i renorm_mask4 = _mm256_xor_si256(Rv4, _mm256_set1_epi32(0x80000000));

	renorm_mask3 = _mm256_cmpgt_epi32(_mm256_set1_epi32(RANS_BYTE_L-0x80000000), renorm_mask3);
	renorm_mask4 = _mm256_cmpgt_epi32(_mm256_set1_epi32(RANS_BYTE_L-0x80000000), renorm_mask4);

	__m256i Yv3 = _mm256_slli_epi32(Rv3, 16);
	__m256i Yv4 = _mm256_slli_epi32(Rv4, 16);

	unsigned int imask3 = _mm256_movemask_ps((__m256)renorm_mask3);
	unsigned int imask4 = _mm256_movemask_ps((__m256)renorm_mask4);
	__m256i idx3 = _mm256_loadu_si256((const __m256i*)permute[imask3]);
	sp += _mm_popcnt_u32(imask3);
	Vv3 = _mm256_permutevar8x32_epi32(Vv3, idx3);

	sv3 = _mm256_packus_epi32(sv3, sv4);
	sv3 = _mm256_permute4x64_epi64(sv3, 0xd8);
	sv3 = _mm256_packus_epi16(sv3, sv3);

	u.tbuf64[tidx][2] = _mm256_extract_epi64(sv3, 0);
	u.tbuf64[tidx][3] = _mm256_extract_epi64(sv3, 2);

	iN[0]++;
	if (++tidx == 32) {
	    iN[0]-=32;

	    for (z = 0; z < NX; z++) {
		// replace by gathers?
		*(uint64_t *)&out[iN[z]] =
		    ((uint64_t)(u.tbuf[0][z])<< 0) + ((uint64_t)(u.tbuf[1][z])<< 8) +
		    ((uint64_t)(u.tbuf[2][z])<<16) + ((uint64_t)(u.tbuf[3][z])<<24) +
		    ((uint64_t)(u.tbuf[4][z])<<32) + ((uint64_t)(u.tbuf[5][z])<<40) +
		    ((uint64_t)(u.tbuf[6][z])<<48) + ((uint64_t)(u.tbuf[7][z])<<56);
		*(uint64_t *)&out[iN[z]+8] =
		    ((uint64_t)(u.tbuf[8+0][z])<< 0) + ((uint64_t)(u.tbuf[8+1][z])<< 8) +
		    ((uint64_t)(u.tbuf[8+2][z])<<16) + ((uint64_t)(u.tbuf[8+3][z])<<24) +
		    ((uint64_t)(u.tbuf[8+4][z])<<32) + ((uint64_t)(u.tbuf[8+5][z])<<40) +
		    ((uint64_t)(u.tbuf[8+6][z])<<48) + ((uint64_t)(u.tbuf[8+7][z])<<56);
		*(uint64_t *)&out[iN[z]+16] =
		    ((uint64_t)(u.tbuf[16+0][z])<< 0) + ((uint64_t)(u.tbuf[16+1][z])<< 8) +
		    ((uint64_t)(u.tbuf[16+2][z])<<16) + ((uint64_t)(u.tbuf[16+3][z])<<24) +
		    ((uint64_t)(u.tbuf[16+4][z])<<32) + ((uint64_t)(u.tbuf[16+5][z])<<40) +
		    ((uint64_t)(u.tbuf[16+6][z])<<48) + ((uint64_t)(u.tbuf[16+7][z])<<56);
		*(uint64_t *)&out[iN[z]+24] =
		    ((uint64_t)(u.tbuf[24+0][z])<< 0) + ((uint64_t)(u.tbuf[24+1][z])<< 8) +
		    ((uint64_t)(u.tbuf[24+2][z])<<16) + ((uint64_t)(u.tbuf[24+3][z])<<24) +
		    ((uint64_t)(u.tbuf[24+4][z])<<32) + ((uint64_t)(u.tbuf[24+5][z])<<40) +
		    ((uint64_t)(u.tbuf[24+6][z])<<48) + ((uint64_t)(u.tbuf[24+7][z])<<56);
		iN[z] += 32;
	    }

	    tidx = 0;
	}

	__m256i idx4 = _mm256_loadu_si256((const __m256i*)permute[imask4]);
	__m256i Vv4 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));

	//Vv = _mm256_and_si256(Vv, renorm_mask);  (blend does the AND anyway)
	Yv3 = _mm256_or_si256(Yv3, Vv3);
	Vv4 = _mm256_permutevar8x32_epi32(Vv4, idx4);
	Yv4 = _mm256_or_si256(Yv4, Vv4);

	sp += _mm_popcnt_u32(imask4);

	Rv3 = _mm256_blendv_epi8(Rv3, Yv3, renorm_mask3);
	Rv4 = _mm256_blendv_epi8(Rv4, Yv4, renorm_mask4);
    }

    STORE(Rv, ransN);
    //_mm256_store_si256((__m256i *)&ransN[ 0], Rv1);
    //_mm256_store_si256((__m256i *)&ransN[ 8], Rv2);
    //_mm256_store_si256((__m256i *)&ransN[16], Rv3);
    //_mm256_store_si256((__m256i *)&ransN[24], Rv4);
    STORE(Lv, lN);
    //_mm256_store_si256((__m256i *)&lN[ 0], Lv1);
    //_mm256_store_si256((__m256i *)&lN[ 8], Lv2);
    //_mm256_store_si256((__m256i *)&lN[16], Lv3);
    //_mm256_store_si256((__m256i *)&lN[24], Lv4);
    ptr = (uint8_t *)sp;

    if (1) {
	iN[0]-=tidx;
	int T;
	for (z = 0; z < NX; z++)
	    for (T = 0; T < tidx; T++)
		out[iN[z]++] = u.tbuf[T][z];
    }

    // Remainder
    z = NX-1;
    for (; iN[z] < out_sz; ) {
	uint32_t m = ransN[z] & ((1u<<TF_SHIFT_O1)-1);
	uint32_t S = s3[lN[z]][m];
	unsigned char c = S & 0xff;
	out[iN[z]++] = c;
	ransN[z] = (S>>(TF_SHIFT_O1+8)) * (ransN[z]>>TF_SHIFT_O1) +
	    ((S>>8) & ((1u<<TF_SHIFT_O1)-1));
	RansDecRenorm(&ransN[z], &ptr);
	lN[z] = c;
    }
    
    *out_size = out_sz;

    return out;
}

int rans_test(unsigned char* __restrict__ in, unsigned int in_size, unsigned char* __restrict__ out, unsigned int* out_size) {
    return 45;
}