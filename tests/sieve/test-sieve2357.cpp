#include "cado.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>
#ifdef HAVE_SSSE3
#include "tmmintrin.h"
#endif
#ifdef HAVE_AVX2
#include "immintrin.h"
#endif
//#define DO_TIMING 1
#ifdef DO_TIMING
// The Jevents library is part of the PMU tools
// https://github.com/andikleen/pmu-tools
#ifdef HAVE_JEVENTS
#define USE_JEVENTS 1
#include "rdtsc.h"
#endif
#endif
#include "utils.h"
#include "tests_common.h"
#include "las-sieve2357.hpp"

template<typename T>
class gettypename {
public:
    static constexpr const char * name = "unknown";
};

template<>
class gettypename<unsigned long> {
public:
    static constexpr const char * name = "unsigned long";
};

template<>
class gettypename<unsigned char> {
public:
    static constexpr const char * name = "unsigned char";
};

#ifdef HAVE_SSSE3
template<>
class gettypename<__m128i> {
public:
    static constexpr const char * name = "__m128i";
};
#endif

#ifdef HAVE_AVX2
template<>
class gettypename<__m256i> {
public:
    static constexpr const char * name = "__m256i";
};
#endif


template <typename SIMDTYPE, typename ELEMTYPE>
bool test_bcaststride(const unsigned int offset, const unsigned int stride) {
  const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
  bool ok = true;
  SIMDTYPE w;
  
  w = bcaststride<SIMDTYPE, ELEMTYPE>(1, offset, stride);

  for (size_t i = 0; i < N; i++) {
    if ((i % stride == offset) != ((ELEMTYPE *)&w)[i]) {
      fprintf(stderr,
              "test_bcaststride<%s, %s>(v=1, offset=%u, stride=%u) wrong: byte[%zu] == %hhu\n",
              gettypename<SIMDTYPE>::name, gettypename<ELEMTYPE>::name,
              offset, stride, i, ((ELEMTYPE *)&w)[i]);
      ok = false;
    }
  }
  return ok;
}

template <typename SIMDTYPE, typename ELEMTYPE>
bool test_bcaststride_many()
{
  const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
  const unsigned int strides[] = {2, 3, 4, 5, 7, 8, 16, 32};
  const size_t nr_strides = sizeof(strides) / sizeof(strides[0]);
  bool ok = true;
  
  for (size_t i = 0; i < nr_strides && strides[i] < N; i++) {
    const unsigned int stride = strides[i];
    for (unsigned int offset = 0; offset < stride; offset++) {
      ok &= test_bcaststride<SIMDTYPE, ELEMTYPE>(offset, stride);
    }
  }
  return ok;
}

template <typename SIMDTYPE, typename ELEMTYPE>
bool test(const unsigned long iter)
{
  const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
  bool ok = true;
  ok &= test_bcaststride_many<SIMDTYPE, ELEMTYPE>();

  const size_t arraysize = false ? N*3*5*7 : 65536;
  if(arraysize % N != 0) {abort();}
#define ARRAY_ON_HEAP 1
#if ARRAY_ON_HEAP
  SIMDTYPE *sievearray = (SIMDTYPE *) malloc_aligned(arraysize, sizeof(SIMDTYPE));
  ELEMTYPE *sievearray2 = (ELEMTYPE *) malloc_aligned(arraysize, sizeof(ELEMTYPE));
#else
  alignas(sizeof(SIMDTYPE)) SIMDTYPE sievearray[arraysize / N];
  ELEMTYPE sievearray2[arraysize];
#endif
  sieve2357_prime_t all_primes[] = {
    /* p q idx logp */
    {3, 1, 0, 2},
    {2, 2, 1, 1},
    {2, 4, 3, 2},
    {2, 8, 0, 6},
    {2, 16, 14, 12},
    {2, 32, 21, 24},
    {3, 3, 2, 3},
    {3, 3, 2, 3},
    {5, 5, 4, 4},
    {7, 7, 5, 5},
  };
  const size_t nr_all_primes = sizeof(all_primes) / sizeof(all_primes[1]);
  sieve2357_prime_t use_primes[nr_all_primes + 1]; /* +1 for end marker */

  /* Skip those prime (powers) this SIMD type can't sieve */
  size_t j = 0;
  for (size_t i = 0; i < nr_all_primes; i++) {
    if (sieve2357_can_sieve<SIMDTYPE, ELEMTYPE>(all_primes[i].p, all_primes[i].q)) {
      use_primes[j++] = all_primes[i];
    }
  }
  use_primes[j] = sieve2357_prime_t{0,0,0,0};

#ifdef DO_TIMING
  sieve2357<SIMDTYPE, ELEMTYPE>(sievearray, arraysize, use_primes);
  sieve2357<SIMDTYPE, ELEMTYPE>(sievearray, arraysize, use_primes);
  start_timing();
#endif
  for (unsigned long i = 0; i < iter; i++) {
    sieve2357<SIMDTYPE, ELEMTYPE>(sievearray, arraysize, use_primes);
  }
#ifdef DO_TIMING
  end_timing();
  printf("%lu calls of sieve2357<%s, %s>(%zu) took %lu cycles, %lu per call\n",
         iter, gettypename<SIMDTYPE>::name, gettypename<ELEMTYPE>::name, arraysize,
         (unsigned long) get_diff_timing(), (unsigned long) get_diff_timing() / iter);
#endif

  /* Do the same thing again, using simple sieve code */
  memset(sievearray2, 0, arraysize);
  for (sieve2357_prime_t *p = use_primes; p->p != 0; p++) {
    for (size_t i = p->idx; i < arraysize; i += p->q) {
      sievearray2[i] += p->logp;
    }
  }

  /* Compare the two arrays */
  unsigned long nr_errors = 0;
  for (size_t i = 0; i < arraysize && nr_errors < 10; i++) {
    if (((ELEMTYPE *)&sievearray[0])[i] != (sievearray2)[i]) {
      printf("Mismatch, sievearray<%s,%s>[%zu] = %hhu, sievearray2[%zu] = %hhu\n",
             gettypename<SIMDTYPE>::name, gettypename<ELEMTYPE>::name, i,
             ((ELEMTYPE *)&sievearray[0])[i], i, sievearray2[i]);
      nr_errors++;
    }
  }
#if ARRAY_ON_HEAP
  free(sievearray);
  free(sievearray2);
#endif
  return nr_errors == 0;
}

int main(int argc, const char **argv)
{
  unsigned long iter = 1;
  bool ok = true;
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);

#ifdef DO_TIMING
  init_timing();
#endif
  ok &= test<unsigned long, unsigned char>(iter);
#ifdef HAVE_SSSE3
  ok &= test<__m128i, unsigned char>(iter);
#endif
#ifdef HAVE_AVX2
  ok &= test<__m256i, unsigned char>(iter);
#endif
#ifdef DO_TIMING
  clear_timing();
#endif
  tests_common_clear();
  exit (ok ? EXIT_SUCCESS : EXIT_FAILURE);
}
