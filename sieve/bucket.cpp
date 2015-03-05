#include "cado.h"
#include <inttypes.h>
#include <stdlib.h>   // for malloc and friends
#include <string.h>   // for memcpy
#include <gmp.h>
#if defined(HAVE_SSE2)
#include <emmintrin.h>
#endif
#include "bucket.h"
#include "portability.h"
#include "memory.h"
#include "las-config.h"
#include "iqsort.h"
#include "verbose.h"
#include "ularith.h"

/* sz is the size of a bucket for an array of buckets. In bytes, a bucket
   size is sz * sr, with sr = sizeof of one element of the bucket (a record).
   If sz * sr is a multiple of the page, and if the buckets array is filled
   *UNIFORMLY & REGULARLY*, all the buckets in the array have the possibility
   to reach their end of the memory pages at the same time.
   In this case, the system must manage multiple TLB misses & must provide
   multiple physical pages of memory in same time = really a bad idea!
   So, in this case, I increase lightly sz.
   NB: I use in fact sz*sr*4 and not sz*sr, because an half or a quart
   of the buckets which reach at the same time the end of their page is
   sufficient to slow down the code.
*/
size_t
bucket_misalignment(const size_t sz, const size_t sr) {
  size_t size; 
  if ((sz * sr * 4) & (pagesize() - 1))
    size = sz;
  else 
    size = sz + 8 / sr + ((8 % sr) ? 1 : 0);
#ifdef HAVE_SSE2
  /* Round up to a multiple of CACHELINESIZE to make SSE2 aligned accesses
     work */
  size = iceildiv(size, 16U) * 16U;
#endif
  return size;
}

/* Not really interesting to do this in SSE2: this function is called max 6 times... */
void
bucket_start_init(void **ps, void **eps, size_t init, size_t add) {
  for ( ; ps + 8 <= eps; ps += 8) {
    ps[0] = (void *) init; init += add;
    ps[1] = (void *) init; init += add;
    ps[2] = (void *) init; init += add;
    ps[3] = (void *) init; init += add;
    ps[4] = (void *) init; init += add;
    ps[5] = (void *) init; init += add;
    ps[6] = (void *) init; init += add;
    ps[7] = (void *) init; init += add;
  }
  for (; ps < eps; *ps++ = (void *) init, init += add);
}

/* Set the write pointers of the normal/kilo/mega buckets, and the read
   pointers of the normal buckets, back to the respective bucket start,
   and set nr_slices back to 0. */
static void
re_init_bucket_array(bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA) {
  if (mBA->bucket_start) aligned_medium_memcpy (mBA->bucket_write, mBA->bucket_start, mBA->size_b_align);
  if (kBA->bucket_start) aligned_medium_memcpy (kBA->bucket_write, kBA->bucket_start, kBA->size_b_align);
  aligned_medium_memcpy (BA->bucket_write, BA->bucket_start, BA->size_b_align);
  aligned_medium_memcpy (BA->bucket_read,  BA->bucket_start, BA->size_b_align);
  BA->nr_slices = 0;
}

static void
init_bucket_array_common(const uint32_t n_bucket, const uint64_t size_bucket, const slice_index_t prealloc_slices, bucket_array_t *BA)
{
  BA->n_bucket = n_bucket;
  BA->bucket_size = size_bucket;
  BA->size_b_align = ((sizeof(void *) * BA->n_bucket + 0x3F) & ~((size_t) 0x3F));
  BA->nr_slices = 0;
  BA->alloc_slices = prealloc_slices;
  BA->bucket_write = (bucket_update_t **) malloc_pagealigned (BA->size_b_align);
  BA->bucket_start = (bucket_update_t **) malloc_aligned (BA->size_b_align, 0x40);
  BA->bucket_read = (bucket_update_t **) malloc_aligned (BA->size_b_align, 0x40);
  BA->slice_index = (slice_index_t  *) malloc_check (BA->alloc_slices * sizeof(slice_index_t));
  BA->slice_start = (bucket_update_t **) malloc_aligned (BA->size_b_align * BA->alloc_slices, 0x40);
}

/* This function is called only in the one pass sort in big buckets sieve.
   The parameter diff_logp should be size_arr_logp, the number of different logp
   in the corresponding fb_iterators.
 */
static void
init_bucket_array(const uint32_t n_bucket, const uint64_t size_bucket, const slice_index_t prealloc_slices, bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA)
{
  init_bucket_array_common(n_bucket, size_bucket, prealloc_slices, BA);
  BA->big_size = BA->n_bucket * BA->bucket_size * sizeof(bucket_update_t);

  verbose_output_print(0, 3, "# Allocating %zu bytes for %" PRIu32 " buckets of %" PRIu64 " update entries of %zu bytes each\n",
                       BA->big_size, BA->n_bucket, BA->bucket_size, sizeof(bucket_update_t));
  uint8_t *big_data = (uint8_t *) physical_malloc (BA->big_size, 1);

  bucket_start_init((void **) BA->bucket_start, (void **) (BA->bucket_start + BA->n_bucket),
		    (size_t) big_data, BA->bucket_size * sizeof(bucket_update_t));

  memset (kBA, 0, sizeof(*kBA));

  memset (mBA, 0, sizeof(*mBA));

  re_init_bucket_array(BA, kBA, mBA);
}

static void 
init_k_bucket_array_common(bucket_array_t *BA, k_bucket_array_t *kBA)
{
  kBA->n_bucket = (BA->n_bucket >> 8) + ((unsigned char) BA->n_bucket != 0 ? 1 : 0); 
  kBA->bucket_size = bucket_misalignment (BA->bucket_size << 8, sizeof(k_bucket_update_t));
  kBA->size_b_align = (sizeof(void *) * kBA->n_bucket + 0x3F) & ~((size_t) 0x3F);
  kBA->bucket_write = (k_bucket_update_t **) malloc_pagealigned (kBA->size_b_align);
  kBA->bucket_start = (k_bucket_update_t **) malloc_aligned (kBA->size_b_align, 0x40);
  kBA->logp_idx = (k_bucket_update_t **) malloc_aligned (kBA->size_b_align * BA->alloc_slices, 0x40);
}

/* This function is called only in the two passes sort in big buckets sieve.
 */
static void 
init_k_bucket_array(const uint32_t n_bucket, const uint64_t size_bucket, const slice_index_t nr_slices, bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA)
{
  init_bucket_array_common(n_bucket, size_bucket, nr_slices, BA);
  init_k_bucket_array_common(BA, kBA);
  
  BA->big_size = kBA->bucket_size * (kBA->n_bucket * sizeof(k_bucket_update_t) + sizeof(bucket_update_t));

#ifdef PRINT_ALLOC
  printf("# Allocating %zu bytes for %" PRIu32 " kilo-buckets of %" PRIu64 " k-update entries of %zu bytes each and one bucket of %" PRIu64 " update entries of %zu bytes each\n",
         BA->big_size, kBA->n_bucket, kBA->bucket_size, sizeof(k_bucket_update_t), kBA->bucket_size, sizeof(bucket_update_t));
#endif
  uint8_t *big_data = (uint8_t *) physical_malloc (BA->big_size, 1);

  bucket_start_init((void **) BA->bucket_start, (void **) (BA->bucket_start + BA->n_bucket),
		    (size_t) big_data, BA->bucket_size * sizeof(bucket_update_t));

  bucket_start_init((void **) kBA->bucket_start, (void **) (kBA->bucket_start + kBA->n_bucket),
		    (size_t) big_data + kBA->bucket_size * sizeof(bucket_update_t),
		    kBA->bucket_size * sizeof(k_bucket_update_t));

  memset (mBA, 0, sizeof(*mBA));

  re_init_bucket_array(BA, kBA, mBA);
}

/* This function is called only in the three passes sort in big buckets sieve.
 */
static void
init_m_bucket_array(const uint32_t n_bucket, const uint64_t size_bucket, const slice_index_t nr_slices, bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA)
{
  init_bucket_array_common(n_bucket, size_bucket, nr_slices, BA);
  init_k_bucket_array_common(BA, kBA);
  
  mBA->n_bucket = (kBA->n_bucket >> 8) + ((unsigned char) kBA->n_bucket != 0 ? 1 : 0); 
  mBA->bucket_size = bucket_misalignment (kBA->bucket_size << 8, sizeof(m_bucket_update_t));

  if (mBA->n_bucket == 1) /* Only for tests */
    mBA->bucket_size = kBA->n_bucket * (kBA->bucket_size + (kBA->bucket_size >> 2));

  mBA->size_b_align = (sizeof(void *) * mBA->n_bucket + 0x3F) & ~((size_t) 0x3F);
  mBA->bucket_write = (m_bucket_update_t **) malloc_pagealigned (mBA->size_b_align);
  mBA->bucket_start = (m_bucket_update_t **) malloc_aligned (mBA->size_b_align, 0x40);
  mBA->logp_idx = (m_bucket_update_t **) malloc_aligned (mBA->size_b_align * BA->alloc_slices, 0x40);

  BA->big_size = mBA->bucket_size * (mBA->n_bucket * sizeof(m_bucket_update_t) + sizeof(k_bucket_update_t)) + kBA->bucket_size * sizeof(bucket_update_t);

#ifdef PRINT_ALLOC
  printf("# Allocating %zu bytes for %" PRIu32 " mega-buckets of %" PRIu64 " update entries each\n",
         BA->big_size, mBA->n_bucket, mBA->bucket_size);

  printf("# Allocating %zu bytes for %" PRIu32 " kilo-buckets of %" PRIu64 " k-update entries of %zu bytes each and one bucket of %" PRIu64 " update entries of %zu bytes each\n",
         BA->big_size, kBA->n_bucket, kBA->bucket_size, sizeof(k_bucket_update_t), kBA->bucket_size, sizeof(bucket_update_t));
#endif
  uint8_t *big_data = (uint8_t *) physical_malloc (BA->big_size, 1);

  bucket_start_init((void **) BA->bucket_start, (void **) (BA->bucket_start + BA->n_bucket),
		    (size_t) big_data, BA->bucket_size * sizeof(bucket_update_t));
  
  bucket_start_init((void **) kBA->bucket_start, (void **) (kBA->bucket_start + kBA->n_bucket),
		    (size_t) big_data + kBA->bucket_size * sizeof(bucket_update_t),
		    kBA->bucket_size * sizeof(k_bucket_update_t));

  bucket_start_init((void **) mBA->bucket_start, (void **) (mBA->bucket_start + mBA->n_bucket),
		    (size_t) big_data + kBA->bucket_size * sizeof(bucket_update_t) +
		    mBA->bucket_size * sizeof(k_bucket_update_t),
		    mBA->bucket_size * sizeof(m_bucket_update_t));

  re_init_bucket_array(BA, kBA, mBA);
}

void
init_buckets(bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA,
             const double max_bucket_fill_ratio, const uint32_t nb_buckets)
{
  const size_t bucket_region = (size_t) 1 << LOG_BUCKET_REGION;
  const size_t bucket_size = bucket_misalignment((size_t) (max_bucket_fill_ratio * bucket_region), sizeof(bucket_update_t));
  /* The previous buckets are identical ? */
  if (BA->n_bucket == nb_buckets && BA->bucket_size == bucket_size) {
    /* Yes; so (bucket_write & bucket_read) = bucket_start; nr_slices = 0 */
    re_init_bucket_array(BA, kBA, mBA);
    /* Buckets are ready to be filled */
  } else {
    /* No. We free the buckets, if we have already malloc them. */
    if (BA->n_bucket)
      clear_buckets(BA, kBA, mBA);
    /* We (re)create the buckets */
    if (nb_buckets < THRESHOLD_K_BUCKETS) {
      init_bucket_array   (nb_buckets, bucket_size, BA->initial_slice_alloc, BA, kBA, mBA);
    } else if (nb_buckets < THRESHOLD_M_BUCKETS) {
      init_k_bucket_array (nb_buckets, bucket_size, BA->initial_slice_alloc, BA, kBA, mBA);
    } else {
      init_m_bucket_array (nb_buckets, bucket_size, BA->initial_slice_alloc, BA, kBA, mBA);
    }
  }
}


void
clear_buckets(bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA)
{
  /* Never free mBA->bucket_start[0]. This is done by free BA->bucket_start[0] */
  if (mBA->logp_idx) free_aligned(mBA->logp_idx);
  if (mBA->bucket_write) free_pagealigned(mBA->bucket_write);
  if (mBA->bucket_start) free_aligned(mBA->bucket_start);
  memset(mBA, 0, sizeof(*mBA));

  /* Never free kBA->bucket_start[0]. This is done by free BA->bucket_start[0] */
  if (kBA->logp_idx) free_aligned(kBA->logp_idx);
  if (kBA->bucket_write) free_pagealigned(kBA->bucket_write);
  if (kBA->bucket_start) free_aligned(kBA->bucket_start);
  memset(kBA, 0, sizeof(*kBA));

  physical_free (BA->bucket_start[0], BA->big_size); /* Always = big_data in all BA, kBA, mBA init functions */
  free (BA->slice_index);
  free_aligned(BA->slice_start);
  free_aligned(BA->bucket_read);
  free_aligned(BA->bucket_start);
  free_pagealigned(BA->bucket_write);
  memset(BA, 0, sizeof(*BA));
}


void bucket_array_t::realloc_slice_start(const size_t extra_space)
{
  const size_t old_nr_entries = alloc_slices;
  alloc_slices += extra_space;
  verbose_output_print(0, 3, "# Reallocating BA->slice_start from %zu entries to %zu entries\n",
                       nr_slices, old_nr_entries);

  const size_t old_size = size_b_align * old_nr_entries;
  const size_t new_size = size_b_align * alloc_slices;
  slice_start = (bucket_update_t **) realloc_aligned(slice_start, old_size, new_size, 0x40);
  ASSERT_ALWAYS(slice_start != NULL);
  slice_index = (slice_index_t *) realloc(slice_index, alloc_slices * sizeof(slice_index_t));
  ASSERT_ALWAYS(slice_index != NULL);
}


/* Returns how full the fullest bucket is */
double
buckets_max_full (const bucket_array_t BA)
{
  unsigned int max = 0;
  for (unsigned int i = 0; i < BA.n_bucket; ++i)
    {
      unsigned int j = nb_of_updates (BA, i);
      if (max < j) max = j;
    }
  return (double) max / (double) BA.bucket_size;
}

/* A compare function suitable for sorting updates in order of ascending x
   with qsort() */
int
bucket_cmp_x (const bucket_complete_update_t *a, const bucket_complete_update_t *b)
{
  if (a->x < b->x)
    return -1;
  if (a->x == b->x)
    return 0;
  return 1;
}

template <class UPDATE_TYPE>
void
bucket_single<UPDATE_TYPE>::sort()
{
//  qsort (start, write - start, sizeof (bucket_complete_update_t),
//	 (int(*)(const void *, const void *)) &bucket_cmp_x);
#define islt(a,b) ((a)->x < (b)->x)
  QSORT(UPDATE_TYPE, start, write - start, islt);
#undef islt  
}

void
bucket_primes_t::purge (const bucket_array_t BA, 
              const int i, const fb_part *fb, const unsigned char *S)
{
  ASSERT_ALWAYS(BA.nr_slices == 0 || BA.begin(i, 0) == BA.bucket_start[i]);

  for (slice_index_t i_slice = 0; i_slice < BA.nr_slices; i_slice++) {
    const slice_index_t slice_index = BA.slice_index[i_slice];
    const bucket_update_t *it = BA.begin(i, i_slice);
    const bucket_update_t * const end_it = BA.end(i, i_slice);

    for ( ; it != end_it ; it++) {
      if (UNLIKELY(S[it->x] != 255)) {
        const fb_slice_interface *slice = fb->get_slice(slice_index);
        ASSERT_ALWAYS(slice != NULL);
        fbprime_t p = slice->get_prime(it->hint);
        bucket_prime_t buc = {it->x, p};
        push_update(buc);
      }
    }
  }
}

void
bucket_array_complete::purge (const bucket_array_t BA, 
              const int i, const unsigned char *S)
{
  ASSERT_ALWAYS(BA.nr_slices == 0 || BA.begin(i, 0) == BA.bucket_start[i]);

  for (slice_index_t i_slice = 0; i_slice < BA.nr_slices; i_slice++) {
    const slice_index_t slice_index = BA.slice_index[i_slice];
    const bucket_update_t *it = BA.begin(i, i_slice);
    const bucket_update_t * const end_it = BA.end(i, i_slice);

    for ( ; it != end_it ; it++) {
      if (UNLIKELY(S[it->x] != 255)) {
        bucket_complete_update_t buc = {it->x, slice_index, it->hint};
        push_update(buc);
      }
    }
  }
}

template class bucket_single<bucket_prime_t>;
template class bucket_single<bucket_complete_update_t>;


void
sieve_checksum::update(const unsigned int other)
{
    unsigned long r;
    ularith_addmod_ul_ul(&r, checksum, other, checksum_prime);
    checksum = r;
}

void
sieve_checksum::update(const unsigned char *data, const size_t len)
{
    mpz_t mb;
    unsigned int new_checksum;

    mpz_init(mb);
    mpz_import(mb, len, -1, sizeof(unsigned char), -1, 0, data);
    new_checksum = mpz_tdiv_ui(mb, checksum_prime);
    mpz_clear(mb);
    this->update(new_checksum);
}
