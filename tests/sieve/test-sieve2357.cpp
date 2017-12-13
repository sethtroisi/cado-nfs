#include "cado.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>
#ifdef HAVE_SSE2
#include "emmintrin.h"
#endif
#ifdef HAVE_AVX2
#include "immintrin.h"
#endif
//#define USE_JEVENTS 1
#include "rdtsc.h"
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

#ifdef HAVE_SSE2
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
              "test_bcaststride<%s, %s>(v=1, offset=%u, stride=%u) wrong: byte[%zu] == %hu\n",
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
void test(const unsigned long iter)
{
  const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
  bool ok = true;
  ok &= test_bcaststride_many<SIMDTYPE, ELEMTYPE>();

  const size_t arraysize = false ? N*3*5*7 : 65536;
  if(arraysize % N != 0) {exit(EXIT_FAILURE);}
  SIMDTYPE sievearray[arraysize / N];
  ELEMTYPE sievearray2[arraysize];
  sieve2357_prime_t primes[10] = {
    /* p q idx logp */
    {2, 2, 1, 1},
    {2, 4, 3, 2},
    {2, 8, 0, 6},
    {3, 3, 2, 3},
    {3, 3, 2, 3},
    {5, 5, 4, 4},
    {7, 7, 5, 5},
  };
  
  start_timing();
  for (unsigned long i = 0; i < iter; i++) {
    sieve2357<SIMDTYPE, ELEMTYPE>(sievearray, arraysize, primes);
  }
  end_timing();
  if (1)
    printf("%lu calls of sieve2357<unsigned char, %zu>(%zu) took %lu cycles, %lu per call\n",
           iter, N, arraysize, (unsigned long) get_diff_timing(), (unsigned long) get_diff_timing() / iter);

  /* Do the same thing again, using simple sieve code */
  memset(sievearray2, 0, arraysize);
  for (sieve2357_prime_t *p = primes; p->p != 0; p++) {
    for (size_t i = p->idx; i < arraysize; i += p->q) {
      sievearray2[i] += p->logp;
    }
  }

  /* Compare the two arrays */
  for (size_t i = 0; i < arraysize; i++) {
    if (((ELEMTYPE *)&sievearray)[i] != (sievearray2)[i]) {
      printf("Mismatch, sievearray[%zu] = %hu, sievearray2[%zu] = %hu\n",
             i, ((ELEMTYPE *)&sievearray)[i], i, sievearray2[i]);
    }
  }
}

int main(int argc, const char **argv)
{
  unsigned long iter = 1;
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);

  init_timing();
  test<unsigned long, unsigned char>(iter);
#ifdef HAVE_SSE2
  test<__m128i, unsigned char>(iter);
#endif
#ifdef HAVE_AVX2
  test<__m256i, unsigned char>(iter);
#endif
  clear_timing();
  tests_common_clear();
  exit (EXIT_SUCCESS);
}
