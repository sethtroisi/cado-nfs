#ifndef SMALLSET_H_
#define SMALLSET_H_

#if defined(HAVE_SSE2)

#include <immintrin.h>
#include <cassert>
#include <vector>
#include <stdint.h>
#include "macros.h"

template <typename RT, typename ST>
static inline
RT mm_set1_m(const ST & d);

template <typename T>
static inline
unsigned mm_movemask_epi8(const T d);

template <>
inline
unsigned mm_movemask_epi8<__m128i>(const __m128i d)
{
  return _mm_movemask_epi8(d);
}

template <int S, typename T>
T mm_cmpeq_epi(const T a, const T b);

template<>
__m128i mm_cmpeq_epi<1, __m128i>(const __m128i a, const __m128i b)
{
  return _mm_cmpeq_epi8(a, b);
}

template<>
__m128i mm_cmpeq_epi<2, __m128i>(const __m128i a, const __m128i b)
{
  return _mm_cmpeq_epi16(a, b);
}

template<>
__m128i mm_cmpeq_epi<4, __m128i>(const __m128i a, const __m128i b)
{
  return _mm_cmpeq_epi32(a, b);
}

#if defined(HAVE_AVX2)

template <>
inline
unsigned mm_movemask_epi8<__m256i>(const __m256i d)
{
  return _mm256_movemask_epi8(d);
}

template<>
__m256i mm_cmpeq_epi<1, __m256i>(const __m256i a, const __m256i b)
{
  return _mm256_cmpeq_epi8(a, b);
}

template<>
__m256i mm_cmpeq_epi<2, __m256i>(const __m256i a, const __m256i b)
{
  return _mm256_cmpeq_epi16(a, b);
}

template<>
__m256i mm_cmpeq_epi<4, __m256i>(const __m256i a, const __m256i b)
{
  return _mm256_cmpeq_epi32(a, b);
}

/* Inline assembly functions for the VPBROADCAST{B,W,D,Q} instructions with
   a memory operand as source. The GCC intrinsics appear to be implemented
   in a way that causes GCC to use only VPBROADCAST? instructions with an
   SSE2 register as operand, so it generates additional code to fetch the
   source data from memory, writes this data to a second memory location,
   loads the SSE2 register from that second memory location (with store
   forwarding stall due to width increase), and then does the broadcast.
   Naturally, this ruins performance. */

template <>
inline
__m128i mm_set1_m<__m128i, uint8_t>(const uint8_t & d)
{
  __m128i r;
  __asm__(
    "vpbroadcastb %1, %0"
    : "=x" (r)
    : "m" (d)
  );
  return r;
}

template <>
inline
__m128i mm_set1_m<__m128i, uint16_t>(const uint16_t & d)
{
  __m128i r;
  __asm__(
    "vpbroadcastw %1, %0"
    : "=x" (r)
    : "m" (d)
  );
  return r;
}

template <>
inline
__m128i mm_set1_m<__m128i, uint32_t>(const uint32_t & d)
{
  __m128i r;
  __asm__(
    "vpbroadcastd %1, %0"
    : "=x" (r)
    : "m" (d)
  );
  return r;
}

template <>
inline
__m128i mm_set1_m<__m128i, uint64_t>(const uint64_t & d)
{
  __m128i r;
  __asm__(
    "vpbroadcastq %1, %0"
    : "=x" (r)
    : "m" (d)
  );
  return r;
}

template <>
inline
__m256i mm_set1_m<__m256i, uint8_t>(const uint8_t & d)
{
  __m256i r;
  __asm__(
    "vpbroadcastb %1, %0"
    : "=x" (r)
    : "m" (d)
  );
  return r;
}

template <>
inline
__m256i mm_set1_m<__m256i, uint16_t>(const uint16_t & d)
{
  __m256i r;
  __asm__(
    "vpbroadcastw %1, %0"
    : "=x" (r)
    : "m" (d)
  );
  return r;
}

template <>
inline
__m256i mm_set1_m<__m256i, uint32_t>(const uint32_t & d)
{
  __m256i r;
  __asm__(
    "vpbroadcastd %1, %0"
    : "=x" (r)
    : "m" (d)
  );
  return r;
}

template <>
inline
__m256i mm_set1_m<__m256i, uint64_t>(const uint64_t & d)
{
  __m256i r;
  __asm__(
    "vpbroadcastq %1, %0"
    : "=x" (r)
    : "m" (d)
  );
  return r;
}

#else

/* Without AVX2, there's no VPBROADCAST instruction, and we use the SSE2
   intrinsics as a (slow) fallback, so the code still runs on SSE2-only
   machines. */

template <>
inline
__m128i mm_set1_m<__m128i, uint8_t>(const uint8_t & d)
{
  return _mm_set1_epi8(d);
}

template <>
inline
__m128i mm_set1_m<__m128i, uint16_t>(const uint16_t & d)
{
  return _mm_set1_epi16(d);
}

template <>
inline
__m128i mm_set1_m<__m128i, uint32_t>(const uint32_t & d)
{
  return _mm_set1_epi32(d);
}

#endif /* if defined(HAVE_AVX2) else */


template <int SIZE, typename ELEMENTTYPE>
class smallset {
public:
#if defined(HAVE_AVX2)
  typedef __m256i storagetype;
#else
  typedef __m128i storagetype;
#endif

private:
  storagetype items[SIZE];

public:
  /* The maximum number of items, each of type ELEMENTTYPE, which we can
     store in "items" */
  static const size_t nr_items = SIZE * sizeof(storagetype) / sizeof(ELEMENTTYPE);

  /* "data" contains "len" entries, each of type ELEMENTTYPE */
  smallset(const ELEMENTTYPE *data, const size_t len) {

    /* Temp array so we can apply padding if necessary. We could also pad with
       clever SSE mixing operations, but don't. Doing so would probably be
       cleaner. */
    ELEMENTTYPE tmp_data[nr_items];

    /* We cannot handle the empty set if SIZE > 0, as we use padding to fill
       unused slots in items, and in an empty set there is nothing we could
       use for the padding. */
    ASSERT_ALWAYS((SIZE == 0 || 0 < len) && len <= nr_items);

    /* Copy the members of the set */
    size_t i;
    for (i = 0; i < len; i++)
      tmp_data[i] = data[i];

    /* If data does not fill entries completely, then the last item in
       data gets duplicated for padding */
    for ( ; i < nr_items; i++)
      tmp_data[i] = tmp_data[i - 1];

    /* Copy everything to items */
    for (i = 0; i < SIZE; i++)
      items[i] = ((storagetype *)&tmp_data)[i];
  }

  smallset(const std::vector<ELEMENTTYPE> &data) {
    ELEMENTTYPE tmp_data[nr_items];
    /* We want the ASSERT_ALWAYS below to instruct gcc that we can't have
     * an empty for loop below (the first one). gcc cannot infer that if
     * we reason on data.size() only, because it does not propagate
     * values the same way, it seems. So we store the length in a
     * variable first.
     */
    const size_t len = data.size();
    ASSERT_ALWAYS((SIZE == 0 || 0 < len) && len <= nr_items);

    size_t i;
    for (i = 0; i < len; i++)
      tmp_data[i] = data[i];

    for ( ; i < nr_items; i++)
      tmp_data[i] = tmp_data[i - 1];

    for (i = 0; i < SIZE; i++)
      items[i] = ((storagetype *)&tmp_data)[i];
  }

  /* Returns true if "item" is in the set, and false otherwise */
  inline bool contains(const ELEMENTTYPE &item) const {
    if (SIZE == 0)
      return false;

    return contains(mm_set1_m<storagetype, ELEMENTTYPE>(item));
  }

  /* Returns true if "item" is in the set, and false otherwise */
  __attribute__((optimize("unroll-all-loops")))
  inline bool contains(const storagetype pattern) const {
    if (SIZE == 0) {
      return false;
    } else if (SIZE <= 8) {
      unsigned u = mm_movemask_epi8(mm_cmpeq_epi<sizeof(ELEMENTTYPE)>(items[0], pattern));
      if (SIZE > 1) {u |= mm_movemask_epi8(mm_cmpeq_epi<sizeof(ELEMENTTYPE)>(items[1], pattern));}
      if (SIZE > 2) {u |= mm_movemask_epi8(mm_cmpeq_epi<sizeof(ELEMENTTYPE)>(items[2], pattern));}
      if (SIZE > 3) {u |= mm_movemask_epi8(mm_cmpeq_epi<sizeof(ELEMENTTYPE)>(items[3], pattern));}
      if (SIZE > 4) {u |= mm_movemask_epi8(mm_cmpeq_epi<sizeof(ELEMENTTYPE)>(items[4], pattern));}
      if (SIZE > 5) {u |= mm_movemask_epi8(mm_cmpeq_epi<sizeof(ELEMENTTYPE)>(items[5], pattern));}
      if (SIZE > 6) {u |= mm_movemask_epi8(mm_cmpeq_epi<sizeof(ELEMENTTYPE)>(items[6], pattern));}
      if (SIZE > 7) {u |= mm_movemask_epi8(mm_cmpeq_epi<sizeof(ELEMENTTYPE)>(items[7], pattern));}
      return u != 0;
    } else /* SIZE > 8 */ {
      unsigned u = mm_movemask_epi8(mm_cmpeq_epi<sizeof(ELEMENTTYPE)>(items[0], pattern));
      for (size_t i = 1; i < SIZE; i++)
        u |= mm_movemask_epi8(mm_cmpeq_epi<sizeof(ELEMENTTYPE)>(items[i], pattern));
      return u != 0;
    }
  }
};

#endif /* if defined(HAVE_SSE2) */
#endif /* ifndef SMALLSET_H_ */
