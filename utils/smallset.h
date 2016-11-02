#ifndef SMALLSET_H_
#define SMALLSET_H_

#ifdef HAVE_SSE2

#include <xmmintrin.h>
#include <cassert>
#include <vector>

template <int SIZE, typename ELEMENTTYPE>
class smallset {
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
    assert((SIZE == 0 || 0 < len) && len <= nr_items);

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
    assert((SIZE == 0 || 0 < data.size()) && data.size() <= nr_items);

    size_t i;
    for (i = 0; i < data.size(); i++)
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
