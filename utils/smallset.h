#ifndef SMALLSET_H_
#define SMALLSET_H_

#ifdef HAVE_SSE2

#include <xmmintrin.h>
#include <cassert>
#include <vector>
#include <stdint.h>
#include "macros.h"

/* Static methods that don't depend on smallset template parameters */
class smallset_tools {
public:
  template <unsigned size_x, unsigned idx_x>
  static inline __m128i
  extract_pattern(const __m128i updates)
  {
    ASSERT_ALWAYS((idx_x + 1) + size_x <= 16);

    if (size_x == 2) {
      if (idx_x < 4) {
        /* Fill lower 4 words with x */
        __m128i pattern = _mm_shufflelo_epi16(updates, idx_x * 0x55);
        /* Interleave lower 4 words with themselves to fill all 8 words */
        return _mm_unpacklo_epi16 (pattern, pattern);
      } else {
        /* Fill upper 4 words with x */
        __m128i pattern = _mm_shufflehi_epi16(updates, (idx_x - 4) * 0x55);
        /* Interleave upper 4 words with themselves to fill all 8 words */
        return _mm_unpackhi_epi16 (pattern, pattern);
      }
    } else if (size_x == 4) {
        return _mm_shuffle_epi32(updates, idx_x * 0x55);
    } else {
      /* Not implemented for bytes atm. */
      abort();
    }
  }
};

template <int SIZE, typename ELEMENTTYPE>
class smallset : public smallset_tools {
  /* Data type could be made a template parameter if desired */
  __m128i items[SIZE];

public:
  /* The maximum number of items, each of type ELEMENTTYPE, which we can
     store in "items" */
  static const size_t nr_items = SIZE * sizeof(__m128i) / sizeof(ELEMENTTYPE);

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
      items[i] = ((__m128i *)&tmp_data)[i];
  }

  smallset(const std::vector<ELEMENTTYPE> &data) {
    ELEMENTTYPE tmp_data[nr_items];
    /* We want the ASSERT_ALWAYS below to instruct gcc that we can't have
     * an empty for loop below (the first one). gcc cannot infer that if
     * we reason on data.size() only, because it does not propagate
     * values the same way, it seems. So we store the length in a
     * variable first.
     */
    size_t len = data.size();
    ASSERT_ALWAYS((SIZE == 0 || 0 < len) && len <= nr_items);

    size_t i;
    for (i = 0; i < len; i++)
      tmp_data[i] = data[i];

    for ( ; i < nr_items; i++)
      tmp_data[i] = tmp_data[i - 1];

    for (i = 0; i < SIZE; i++)
      items[i] = ((__m128i *)&tmp_data)[i];
  }

  /* Returns true if "item" is in the set, and false otherwise */
  inline bool contains(const ELEMENTTYPE item) const {
    if (SIZE == 0)
      return false;

    if (sizeof(ELEMENTTYPE) == 1) {
      const __m128i c = _mm_set1_epi8(item);
      unsigned u = _mm_movemask_epi8(_mm_cmpeq_epi8(items[0], c));
      /* TODO: Unroll this */
      for (size_t i = 1; i < SIZE; i++)
        u |= _mm_movemask_epi8(_mm_cmpeq_epi8(items[i], c));
      return u != 0;
    } else if (sizeof(ELEMENTTYPE) == 2) {
      const __m128i c = _mm_set1_epi16(item);
      unsigned u = _mm_movemask_epi8(_mm_cmpeq_epi16(items[0], c));
      /* TODO: Unroll this */
      for (size_t i = 1; i < SIZE; i++)
        u |= _mm_movemask_epi8(_mm_cmpeq_epi16(items[i], c));
      return u != 0;
    } else if (sizeof(ELEMENTTYPE) == 4) {
      const __m128i c = _mm_set1_epi32(item);
      unsigned u = _mm_movemask_epi8(_mm_cmpeq_epi32(items[0], c));
      /* TODO: Unroll this */
      for (size_t i = 1; i < SIZE; i++)
        u |= _mm_movemask_epi8(_mm_cmpeq_epi32(items[i], c));
      return u != 0;
    } else {
      abort();
    }
  }
};

#endif
#endif
