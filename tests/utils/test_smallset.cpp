#include "cado_config.h"
#include <cstdio>
#include <string.h>
#include "smallset.h"

#ifdef HAVE_SSE2


template <int SIZE, typename ELEMENTTYPE, int IDX>
void
check_pattern(smallset<SIZE, ELEMENTTYPE> &set MAYBE_UNUSED, size_t n, __m128i data, size_t i MAYBE_UNUSED)
{
   if (n > IDX) {
      __m128i pattern MAYBE_UNUSED = smallset_tools::extract_pattern<sizeof(ELEMENTTYPE),IDX>(data);
      assert(set.contains(pattern) == (i == IDX));
    }
}

template <int SIZE, typename ELEMENTTYPE>
void
test_smallset()
{
  const size_t nr_items = smallset<SIZE, ELEMENTTYPE>::nr_items;
  ELEMENTTYPE items[nr_items];
  for (size_t i = 0; i < nr_items; i++)
    items[i] = i;

  if (SIZE == 0) {
    smallset<SIZE, ELEMENTTYPE> set(items, 0);
    assert(!set.contains(0));
  }

    for (size_t i = 1; i <= nr_items; i++) {
    smallset<SIZE, ELEMENTTYPE> set(items, i);
    assert(set.contains(0));
    assert(set.contains(i-1));
    assert(!set.contains(i));
    __m128i pattern MAYBE_UNUSED;
    if (sizeof(ELEMENTTYPE) == 1)
      pattern = _mm_set1_epi8(2);
    else if(sizeof(ELEMENTTYPE) == 2)
      pattern = _mm_set1_epi16(2);
    else if(sizeof(ELEMENTTYPE) == 4)
      pattern = _mm_set1_epi32(2);
    assert(set.contains(pattern) != (i <= 2));
  }

  if (SIZE > 0 && sizeof(ELEMENTTYPE) > 1) {
    ELEMENTTYPE arr[sizeof(__m128i) / sizeof(ELEMENTTYPE)];
    const size_t n = sizeof(__m128i) / sizeof(ELEMENTTYPE);
    for (size_t i = 0; i < n; i++)
      arr[i] = (ELEMENTTYPE) i;
    smallset<SIZE, ELEMENTTYPE> set(arr, n);

    for (size_t j = 0; j < n - 1; j++) {
      for (size_t i = 0; i < n; i++) {
        memset(arr, 255, sizeof(__m128i));
        arr[i] = (ELEMENTTYPE) j;
        __m128i data = *(__m128i *)arr;

        check_pattern<SIZE, ELEMENTTYPE, 0>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 1>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 2>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 3>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 4>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 5>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 6>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 7>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 8>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 9>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 9>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 10>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 11>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 12>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 13>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 14>(set, n, data, i);
        check_pattern<SIZE, ELEMENTTYPE, 15>(set, n, data, i);
      }
    }
  }
}

template <typename ELEMENTTYPE>
void
test_one_type()
{
  test_smallset<0, ELEMENTTYPE>();
  test_smallset<1, ELEMENTTYPE>();
  test_smallset<2, ELEMENTTYPE>();
  test_smallset<3, ELEMENTTYPE>();
  test_smallset<4, ELEMENTTYPE>();
  test_smallset<5, ELEMENTTYPE>();
}

int main()
{
  test_one_type<unsigned char>();
  test_one_type<unsigned short>();
  test_one_type<unsigned int>();
  return(0);
}

#else

int main()
{
  return(0);
}

#endif
