#include "cado.h"
#include "bucket.h"
#include "portability.h"
#include "memory.h"

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
uint64_t
bucket_misalignment(const uint64_t sz, const size_t sr) {
  if (LIKELY((sz * (sr << 2)) & (pagesize() - 1)))
    return sz;
  else 
    return sz + 8 / sr + ((8 % sr) ? 1 : 0);
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

void
re_init_bucket_array(bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA) {
  if (mBA->bucket_start) aligned_medium_memcpy (mBA->bucket_write, mBA->bucket_start, mBA->size_b_align);
  if (kBA->bucket_start) aligned_medium_memcpy (kBA->bucket_write, kBA->bucket_start, kBA->size_b_align);
  aligned_medium_memcpy (BA->bucket_write, BA->bucket_start, BA->size_b_align);
  aligned_medium_memcpy (BA->bucket_read,  BA->bucket_start, BA->size_b_align);
  BA->nr_logp = 0;
}

/* This function is called only in the one pass sort in big buckets sieve.
   The parameter diff_logp should be size_arr_logp, the number of different logp
   in the corresponding fb_iterators.
 */
void
init_bucket_array(const uint32_t n_bucket, const uint64_t size_bucket, const unsigned char diff_logp, bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA)
{
  long pagesz = pagesize();
  BA->n_bucket = n_bucket;
  BA->bucket_size = size_bucket;
  BA->size_b_align = ((sizeof(void *) * BA->n_bucket + 0x3F) & ~((size_t) 0x3F));
  BA->nr_logp = 0;
  BA->size_arr_logp = diff_logp;
  BA->bucket_write = (bucket_update_t **) malloc_aligned (BA->size_b_align, pagesz);
  BA->bucket_start = (bucket_update_t **) malloc_aligned (BA->size_b_align, 0x40);
  BA->bucket_read = (bucket_update_t **) malloc_aligned (BA->size_b_align, 0x40);
  BA->logp_val = (unsigned char *) malloc_check (BA->size_arr_logp);
  BA->logp_idx = (bucket_update_t **) malloc_aligned (BA->size_b_align * BA->size_arr_logp, 0x40);

  uint8_t *big_data = physical_malloc (BA->n_bucket * BA->bucket_size * sizeof(bucket_update_t), 1);

  bucket_start_init((void **) BA->bucket_start, (void **) (BA->bucket_start + BA->n_bucket),
		    (size_t) big_data, BA->bucket_size * sizeof(bucket_update_t));

  memset (kBA, 0, sizeof(*kBA));

  memset (mBA, 0, sizeof(*mBA));

  re_init_bucket_array(BA, kBA, mBA);
}

/* This function is called only in the two passes sort in big buckets sieve.
 */
void 
init_k_bucket_array(const uint32_t n_bucket, const uint64_t size_bucket, const unsigned char diff_logp, bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA)
{
  BA->n_bucket = n_bucket;
  BA->bucket_size = size_bucket;
  BA->size_b_align = ((sizeof(void *) * BA->n_bucket + 0x3F) & ~((size_t) 0x3F));
  BA->nr_logp = 0;
  BA->size_arr_logp = diff_logp;
  BA->bucket_write = (bucket_update_t **) malloc_pagealigned (BA->size_b_align);
  BA->bucket_start = (bucket_update_t **) malloc_aligned (BA->size_b_align, 0x40);
  BA->bucket_read = (bucket_update_t **) malloc_aligned (BA->size_b_align, 0x40);
  BA->logp_val = (unsigned char *) malloc_check (BA->size_arr_logp);
  BA->logp_idx = (bucket_update_t **) malloc_aligned (BA->size_b_align * BA->size_arr_logp, 0x40);
  
  kBA->n_bucket = (BA->n_bucket >> 8) + ((unsigned char) BA->n_bucket != 0 ? 1 : 0); 
  kBA->bucket_size = bucket_misalignment (BA->bucket_size << 8, sizeof(k_bucket_update_t));
  kBA->size_b_align = (sizeof(void *) * kBA->n_bucket + 0x3F) & ~((size_t) 0x3F);
  kBA->bucket_write = (k_bucket_update_t **) malloc_pagealigned (kBA->size_b_align);
  kBA->bucket_start = (k_bucket_update_t **) malloc_aligned (kBA->size_b_align, 0x40);
  kBA->logp_idx = (k_bucket_update_t **) malloc_aligned (kBA->size_b_align * BA->size_arr_logp, 0x40);

  uint8_t *big_data = physical_malloc (kBA->bucket_size * (kBA->n_bucket * sizeof(k_bucket_update_t) + sizeof(bucket_update_t)), 1);

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
void
init_m_bucket_array(const uint32_t n_bucket, const uint64_t size_bucket, const unsigned char diff_logp, bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA)
{
  BA->n_bucket = n_bucket;
  BA->bucket_size = size_bucket;
  BA->size_b_align = ((sizeof(void *) * BA->n_bucket + 0x3F) & ~((size_t) 0x3F));
  BA->nr_logp = 0;
  BA->size_arr_logp = diff_logp;
  BA->bucket_write = (bucket_update_t **) malloc_pagealigned (BA->size_b_align);
  BA->bucket_start = (bucket_update_t **) malloc_aligned (BA->size_b_align, 0x40);
  BA->bucket_read = (bucket_update_t **) malloc_aligned (BA->size_b_align, 0x40);
  BA->logp_val = (unsigned char *) malloc_check (BA->size_arr_logp);
  BA->logp_idx = (bucket_update_t **) malloc_aligned (BA->size_b_align * BA->size_arr_logp, 0x40);
  
  kBA->n_bucket = (BA->n_bucket >> 8) + ((unsigned char) BA->n_bucket != 0 ? 1 : 0); 
  kBA->bucket_size = bucket_misalignment (BA->bucket_size << 8, sizeof(k_bucket_update_t));
  kBA->size_b_align = (sizeof(void *) * kBA->n_bucket + 0x3F) & ~((size_t) 0x3F);
  kBA->bucket_write = (k_bucket_update_t **) malloc_pagealigned (kBA->size_b_align);
  kBA->bucket_start = (k_bucket_update_t **) malloc_aligned (kBA->size_b_align, 0x40);
  kBA->logp_idx = (k_bucket_update_t **) malloc_aligned (kBA->size_b_align * BA->size_arr_logp, 0x40);

  mBA->n_bucket = (kBA->n_bucket >> 8) + ((unsigned char) kBA->n_bucket != 0 ? 1 : 0); 
  mBA->bucket_size = bucket_misalignment (kBA->bucket_size << 8, sizeof(m_bucket_update_t));

  if (mBA->n_bucket == 1) /* Only for tests */
    mBA->bucket_size = kBA->n_bucket * (kBA->bucket_size + (kBA->bucket_size >> 2));

  mBA->size_b_align = (sizeof(void *) * mBA->n_bucket + 0x3F) & ~((size_t) 0x3F);
  mBA->bucket_write = (m_bucket_update_t **) malloc_pagealigned (mBA->size_b_align);
  mBA->bucket_start = (m_bucket_update_t **) malloc_aligned (mBA->size_b_align, 0x40);
  mBA->logp_idx = (m_bucket_update_t **) malloc_aligned (mBA->size_b_align * BA->size_arr_logp, 0x40);

  uint8_t *big_data = physical_malloc (mBA->bucket_size * (mBA->n_bucket * sizeof(m_bucket_update_t) + sizeof(k_bucket_update_t)) + kBA->bucket_size * sizeof(bucket_update_t), 1);

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
clear_bucket_array(bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA)
{
  /* Never free mBA->bucket_start[0]. This is done by free BA->bucket_start[0] */
  if (mBA->logp_idx) free_aligned(mBA->logp_idx, 0x40);
  if (mBA->bucket_write) free_pagealigned(mBA->bucket_write);
  if (mBA->bucket_start) free_aligned(mBA->bucket_start, 0x40);
  memset(mBA, 0, sizeof(*mBA));

  /* Never free kBA->bucket_start[0]. This is done by free BA->bucket_start[0] */
  if (kBA->logp_idx) free_aligned(kBA->logp_idx, 0x40);
  if (kBA->bucket_write) free_pagealigned(kBA->bucket_write);
  if (kBA->bucket_start) free_aligned(kBA->bucket_start, 0x40);
  memset(kBA, 0, sizeof(*kBA));

  free (BA->bucket_start[0]); /* Always = big_data in all BA, kBA, mBA init functions */
  free (BA->logp_val);
  free_aligned(BA->logp_idx, 0x40);
  free_aligned(BA->bucket_read, 0x40);
  free_aligned(BA->bucket_start, 0x40);
  free_pagealigned(BA->bucket_write);
  memset(BA, 0, sizeof(*BA));
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

bucket_primes_t
init_bucket_primes (const int size)
{
  bucket_primes_t BP;
  BP.size = size;
  BP.start = (bucket_prime_t *) malloc_check (size * sizeof(bucket_prime_t));
  BP.read = BP.start;
  BP.write = BP.start;
  return BP;
}

void
clear_bucket_primes (bucket_primes_t *BP)
{
  free (BP->start);
  BP->start = NULL;
  BP->read = NULL;
  BP->write = NULL;
  BP->size = 0;
}

/* A compare function suitable for sorting updates in order of ascending x
   with qsort() */
int
bucket_cmp_x (const bucket_prime_t *a, const bucket_prime_t *b)
{
  if (a->x < b->x)
    return -1;
  if (a->x == b->x)
    return 0;
  return 1;
}

void
bucket_sortbucket (bucket_primes_t *BP)
{
  qsort (BP->start, BP->write - BP->start, sizeof (bucket_prime_t), 
	 (int(*)(const void *, const void *)) &bucket_cmp_x);
}


/* Copy only those bucket entries where x yields a sieve report.
 * These entries get sorted, to speed up trial division. 
 * Due to the purging and sorting, it will not be possible to
 * reconstruct the correct p from its low 16 bits, so the
 * reconstruction is done here and the full p is stored in the output.
 */

#ifdef BUCKET_CAREFUL_DECODE
#ifdef BUCKET_ENCODE3
#define PURGE_BUCKET_HEART(A) do {					\
    prime_hint_t up = (u + (A))->p;					\
    if (UNLIKELY(up < last_p)) phigh += BUCKET_P_WRAP;			\
    uint32_t decoded = phigh + bucket_decode_prime(up); last_p = up;	\
    if (UNLIKELY(decoded * 0xCCCCCCCDU <= 0x33333333U)) { /* Divisible by 5? */ \
      decoded += BUCKET_P_WRAP; phigh += BUCKET_P_WRAP; }		\
    uint16_t ux = (u + (A))->x;						\
    if (UNLIKELY(S[ux] != 255)) {					\
      bucket_prime_t bp;						\
      bp.x = ux; bp.p = decoded; push_bucket_prime (BP, bp); }		\
  } while (0)
#else
#define PURGE_BUCKET_HEART(A) do {					\
    prime_hint_t up = (u + (A))->p;					\
    if (UNLIKELY(up < last_p)) phigh += BUCKET_P_WRAP;			\
    uint32_t decoded = phigh + bucket_decode_prime(up); last_p = up;	\
    if (UNLIKELY(decoded * 0xAAAAAAABU <= 0x55555555U)) { /* Divisible by 3? */ \
      decoded += BUCKET_P_WRAP; phigh += BUCKET_P_WRAP; }		\
    uint16_t ux = (u + (A))->x;						\
    if (UNLIKELY(S[ux] != 255)) {					\
      bucket_prime_t bp;						\
      bp.x = ux; bp.p = decoded; push_bucket_prime (BP, bp); }		\
  } while (0)
#endif
#else
#define PURGE_BUCKET_HEART(A) do {					\
    prime_hint_t up = (u + (A))->p;					\
    if (UNLIKELY(up < last_p)) phigh += BUCKET_P_WRAP;			\
    last_p = up;							\
    uint16_t ux = (u + (A))->x;						\
    if (UNLIKELY(S[ux] != 255)) {					\
      bucket_prime_t bp;						\
      bp.x = ux; bp.p = phigh + bucket_decode_prime(up);		\
      push_bucket_prime (BP, bp); }					\
  } while (0)
#endif

void
purge_bucket (bucket_primes_t *BP, const bucket_array_t BA, 
              const int i, const unsigned char *S)
{
  bucket_update_t *u = BA.bucket_start[i], *end_u = BA.bucket_write[i];
  prime_hint_t last_p = 0;
  uint32_t phigh = 0;

  for (; u + 16 <= end_u; u += 16) {
#ifdef HAVE_SSE2
    _mm_prefetch((uint8_t *) u + 0x100, _MM_HINT_T0);
#endif
    PURGE_BUCKET_HEART( 0); PURGE_BUCKET_HEART( 1); PURGE_BUCKET_HEART( 2); PURGE_BUCKET_HEART( 3);
    PURGE_BUCKET_HEART( 4); PURGE_BUCKET_HEART( 5); PURGE_BUCKET_HEART( 6); PURGE_BUCKET_HEART( 7);
    PURGE_BUCKET_HEART( 8); PURGE_BUCKET_HEART( 9); PURGE_BUCKET_HEART(10); PURGE_BUCKET_HEART(11);
    PURGE_BUCKET_HEART(12); PURGE_BUCKET_HEART(13); PURGE_BUCKET_HEART(14); PURGE_BUCKET_HEART(15);
  }
  for (; u < end_u; ++u) PURGE_BUCKET_HEART(0);
}
