#include "cado.h"
#include "bucket.h"
#include "portability.h"

bucket_array_t
init_bucket_array(const int n_bucket, const int bucket_size)
{
    bucket_array_t BA;
    int i;
    BA.bucket_size = bucket_size;
    BA.n_bucket = n_bucket;

    BA.bucket_start = (bucket_update_t **)
      malloc_pagealigned (n_bucket * sizeof(bucket_update_t *));
    if (BA.bucket_start == NULL)
      {
        fprintf (stderr, "Error, cannot allocate memory\n");
        exit (EXIT_FAILURE);
      }
    BA.bucket_write = (bucket_update_t **)
      malloc_check (n_bucket * sizeof(bucket_update_t *));
    BA.bucket_read  = (bucket_update_t **)
      malloc_check (n_bucket * sizeof(bucket_update_t *));

#ifdef  ONE_BIG_MALLOC
    {
      size_t alloc = (size_t)n_bucket * (size_t)bucket_size
          * sizeof(bucket_update_t);
      BA.bucket_start[0] = (bucket_update_t *) malloc_check (alloc);
    }
#endif

    for (i = 0; i < n_bucket; ++i) {
        // TODO: shall we ensure here that those pointer do not differ by
        // a large power of 2, to take into account the associativity of
        // L1 cache ?
#ifdef  ONE_BIG_MALLOC
        BA.bucket_start[i] = BA.bucket_start[0] + i * bucket_size;
#else
        BA.bucket_start[i] = (bucket_update_t *)
          malloc_check (bucket_size * sizeof(bucket_update_t));
#endif
        BA.bucket_write[i] = BA.bucket_start[i];
        BA.bucket_read[i] = BA.bucket_start[i];
    }
    BA.logp_val = NULL;
    BA.logp_idx = NULL;
    BA.nr_logp = 0;
    BA.last_logp = 0;
    return BA;
}

void
clear_bucket_array(bucket_array_t BA)
{
#ifdef  ONE_BIG_MALLOC
    free (BA.bucket_start[0]);
#else
    int i;
    for (i = 0; i < BA.n_bucket; ++i)
      free (BA.bucket_start[i]);
#endif
    free_pagealigned(BA.bucket_start, BA.n_bucket*sizeof(bucket_update_t *));
    free(BA.bucket_write);
    BA.bucket_write = NULL;
    free(BA.bucket_read);
    BA.bucket_read = NULL;
    free(BA.logp_val);
    BA.logp_val = NULL;
    free(BA.logp_idx);
    BA.logp_idx = NULL;
}

/* Returns how full the fullest bucket is */
double
buckets_max_full (const bucket_array_t BA)
{
  int i, max = 0;
  for (i = 0; i < BA.n_bucket; ++i)
    {
      int j = nb_of_updates (BA, i);
      if (max < j)
        max = j;
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

void
purge_bucket (bucket_primes_t *BP, bucket_array_t BA, 
              const int i, const unsigned char *S)
{
  bucket_update_t *u;
  uint16_t last_p = 0;
  uint32_t phigh = 0;
  bucket_prime_t bp;

  for (u = BA.bucket_start[i] ; u < BA.bucket_write[i]; u++)
    {
      uint32_t decoded;
      if (u->p < last_p)
	phigh += BUCKET_P_WRAP;
      last_p = u->p;
      decoded = phigh + bucket_decode_prime(u->p);
#ifdef BUCKET_CAREFUL_DECODE
      if (
#ifndef BUCKET_ENCODE3
          decoded * 0xAAAAAAABU <= 0x55555555U /* Divisible by 3? */
#else
          decoded * 0xCCCCCCCDU <= 0x33333333U /* Divisible by 5? */
#endif
        ) {
        decoded += BUCKET_P_WRAP;
        phigh += BUCKET_P_WRAP;
      }
#endif
      if (S[u->x] != 255)
        {
	  bp.p = decoded;
          bp.x = u->x;
          push_bucket_prime (BP, bp);
	}
    }
}
