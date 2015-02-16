#ifndef BUCKET_H_
#define BUCKET_H_

/*
 * Bucket sieving: radix-sort sieve updates as they are created.
 */

#include <stdint.h>
#include "cado-endian.h"
// #define SAFE_BUCKETS
#ifdef SAFE_BUCKETS
#include <stdio.h>
#include "portability.h"
#endif

/* If BUCKET_ENCODE3 is defined, we encode primes p as
   (floor(p/6) * 2 + p%3-1) % 2^16,
   i.e., we assume that p is odd and not divisible by 3,
   and store the residue class modulo 3 in the LSB. */
// #define BUCKET_ENCODE3

/* We store ((p-1)/2) % 2^16 or (floor(p/6) * 2 + p%3-1) % 2^16 in the bucket 
   updates and reconstruct p from that. In rare cases a wrap-around is 
   missed and p does not divide the norm; in such a case we try 
   p + k*P_WRAP for successive values of k until we find one that divides */
#ifdef BUCKET_ENCODE3
#define BUCKET_P_WRAP (3U*(1U<<16))
#else
#define BUCKET_P_WRAP (1U<<17)
#endif

/* If BUCKET_CAREFUL_DECODE is defined, purge_bucket() tests whether 
   reconstructed primes are divisible by 3 (or 5 if BUCKET_ENCODE3 is 
   defined) which indicates that a wrap-around was missed, and if so adds 
   BUCKET_P_WRAP. This reduces the number of incorrectly reconstructed 
   primes but is a bit slower. */
// #define BUCKET_CAREFUL_DECODE

/*
 * This bucket module provides a way to store elements (that are called
 * updates), while partially sorting them, according to some criterion (to
 * be defined externally): the incoming data is stored into several
 * buckets. The user says for each data to which bucket it belongs. This
 * module is supposed to perform this storage in a cache-friendly way and
 * so on...
 */

/*
 * Main available commands are 
 *   push_bucket_update(i, x)  :   add x to bucket number i
 *   get_next_bucket_update(i) :   iterator to read contents of bucket number i
 *
 * See the MAIN FUNCTIONS section below for the complete interface, with
 * prototypes of exported functions.
 */

/********  Data structure for the contents of buckets **************/

/* In principle, the typedef for the bucket_update can be changed without
 * affecting the rest of the code. 
 */

/*
 * For the moment, we store the bucket updates and a 16-bit field
 * that can contain, for instance, the low bits of p.
 */

/* define prime_hint_t to uint32_t to use 32-bit prime hints */
typedef uint32_t prime_hint_t;

/* THE ORDER IS MANDATORY! */
typedef struct {
#ifdef CADO_LITTLE_ENDIAN
    prime_hint_t p;
    uint16_t x;
#else
    uint16_t x;
    prime_hint_t p;
#endif
} bucket_update_t;

/* Same than previous, but 1 byte more for indirection.
   Careful: this byte is the first in little endian, and
   the last byte in big endian.
   k = kilo bucket */
typedef uint8_t k_bucket_update_t[1 + sizeof(bucket_update_t)];

/* Same than previous, but with 2 bytes (double 1 byte indirections).
   m = mega bucket */
typedef uint8_t m_bucket_update_t[2 + sizeof(bucket_update_t)];

/* After a bucket is purged, it is usually not possible to reconstruct
 * the primes from their lower 16 bits any more, so purging does the 
 * reconstruction and stores the remaining updates with the complete primes 
 */

typedef struct {
    uint16_t x;
    uint32_t p;
} bucket_prime_t;


/*
 * will be used as a sentinel
 */
static const bucket_update_t LAST_UPDATE = {0,0};

/******** Bucket array typedef **************/
typedef struct {
  bucket_update_t ** bucket_write;    // Contains pointers to first empty
                                      // location in each bucket
  bucket_update_t ** bucket_start;    // Contains pointers to beginning of
                                      // buckets
  bucket_update_t ** bucket_read;     // Contains pointers to first unread
                                      // location in each bucket
  unsigned char    * logp_val;        // Each time logp changes, the
                                      // new value is added here (including
                                      // the initial logp value used) 
  bucket_update_t ** logp_idx;        // For each logp value there are
                                      // n_bucket pointers, each pointer
                                      // tells where in the corresponding 
                                      // bucket that logp values starts 
                                      // being used
  uint32_t           n_bucket;        // Number of buckets
  uint64_t           bucket_size;     // The allocated size of one bucket.
  size_t             size_b_align;    // (sizeof(void *) * n_bucket + 63) & ~63
                                      // to align bucket_* on cache line
  size_t big_size;                    // size of bucket update memory
  unsigned int       nr_logp;         // Number of different logp
  unsigned int       size_arr_logp;   // size array logp_val & idx. 256 by default
                                      // or less if the number of different logp
                                      // is known
} bucket_array_t;

/* Same than previous, for kilo & mega buckets in fill_in_k/m_buckets.
   These buckets are used only for 2 and 3 sort passes buckets */ 
typedef struct {
  k_bucket_update_t ** bucket_write;
  k_bucket_update_t ** bucket_start;
  k_bucket_update_t ** logp_idx;
  uint32_t             n_bucket;
  uint64_t             bucket_size;
  size_t               size_b_align;
} k_bucket_array_t;

typedef struct {
  m_bucket_update_t ** bucket_write;
  m_bucket_update_t ** bucket_start;
  m_bucket_update_t ** logp_idx;
  uint32_t             n_bucket;
  uint64_t             bucket_size;
  size_t               size_b_align;
} m_bucket_array_t;

/* What does the logp entry do?
   In order to avoid having to store the logp value inside the 
   bucket_update_t, we use the fact that logp changes very rarely,
   and only store pointers to the locations inside the buckets where
   logp changes.
   Assuming that very small primes are not bucket sieved, and that
   factor base entries are sorted in order of increasing logp, 
   \footnote{with maybe a few exceptions where maintainging strict order 
     is too awkward, but note that 
     log(largest-bucket-sieved-p) - log(smallest-bucket-sieved-p)
     is really a rather small value, and the number of exceptions 
     should not be very much larger}
   we store a pointer for *every* bucket when logp changes, even though 
   not necessarily every bucket receives an update with that logp. We do 
   this because we'd like to have a single last_logp value to compare to, 
   rather than n_bucket different ones.
   
   Thus there is an array logp_val with the logp values, starting with 
   the first logp value to use, plus one entry whenever the logp value 
   changes. For each entry in logp_val, logp_idx has a pointer to an 
   array of n_bucket pointers, each pointer giving the address inside
   the corresponding bucket from where on this new logp values should be 
   used.
   I.e., logp_val[j] should be used for the updates in bucket[i] whose 
   addresses are >= logp_idx[j*n_bucket+i], but less than 
     logp_idx[(j+1)*n_bucket+i] if j+1 < nr_logp, bucket_write otherwise.
*/


/* Similar, but stores info containing the complete prime instead of
   only the low 16 bits, and has only one bucket */

typedef struct {
  bucket_prime_t *start;
  bucket_prime_t *read;
  bucket_prime_t *write;
  int size;
} bucket_primes_t;

/******** MAIN FUNCTIONS **************/

#ifdef __cplusplus
extern "C" {
#endif

/* Set an allocated array of <nb_buckets> buckets each having size
 * max_bucket_fill_ratio * BUCKET_REGION.
 * Must be freed with clear_buckets().
 */
extern void init_buckets(bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA,
                         double max_bucket_fill_ratio, uint32_t nb_buckets);
/* Only one (clever) function to clear the buckets. */
extern void clear_buckets(bucket_array_t *BA, k_bucket_array_t *kBA, m_bucket_array_t *mBA);

/* Main writing function: appends update to bucket number i.
 * If SAFE_BUCKETS is not #defined, then there is no checking that there is
 * enough room for the update. This could lead to a segfault, with the
 * current implementation!
 */
static inline void
push_bucket_update(bucket_array_t BA, const int i, 
                   const bucket_update_t update);

/* Main reading iterator: returns the first unread update in bucket i.
 * If SAFE_BUCKETS is not #defined, there is no checking that you are reading
 * at most as much as what you wrote. If keeping the count while pushing
 * is impossible, the following functions can help:
 * - get the count after the last push using nb_of_updates()
 * - put a sentinel after your bucket using push_sentinel() and check if
 * the value returned by get_next_bucket_update() is LAST_UPDATE.
 * - call the is_end() functions that tells you if you can apply
 *   get_next_bucket_update().
 */
static inline bucket_update_t
get_next_bucket_update(bucket_array_t BA, const int i);
static inline int
nb_of_updates(const bucket_array_t BA, const int i);
static inline void
push_sentinel(bucket_array_t BA, const int i);
static inline int
is_end(const bucket_array_t BA, const int i);

/* If you need to read twice a bucket, you can rewind it: */
static inline void
rewind_bucket(bucket_array_t BA, const int i);

/* If you want to read the most recently read update again, rewind by 1: */
static inline void
rewind_bucket_by_1 (bucket_array_t BA, const int i);

/* If you want to access updates in a non-sequential way: */
static inline bucket_update_t
get_kth_bucket_update(const bucket_array_t BA, const int i, const int k);

/* Functions for handling entries with x and complete prime p */

extern bucket_primes_t init_bucket_primes (const int size);

extern void clear_bucket_primes (bucket_primes_t *BP);

static inline void
push_bucket_prime (bucket_primes_t *BP, const bucket_prime_t prime);

static inline bucket_prime_t
get_next_bucket_prime (bucket_primes_t *BP);

extern void purge_bucket (bucket_primes_t *BP, bucket_array_t BA, const int i, const unsigned char *S);

extern size_t bucket_misalignment(const size_t sz, const size_t sr);

/* We also forward-define some auxiliary functions which are defined in
 * bucket.c (alongside with the non-inlined functions already listed).
 */
extern double buckets_max_full (const bucket_array_t BA);
extern void bucket_sortbucket (bucket_primes_t *BP);

#ifdef __cplusplus
}
#endif

/******** Bucket array implementation **************/

#include "utils/misc.h"

static inline void
push_bucket_update(bucket_array_t BA, const int i, 
                   const bucket_update_t update)
{
    *(BA.bucket_write[i])++ = update; /* Pretty! */
#ifdef SAFE_BUCKETS
    if (BA.bucket_start[i] + BA.bucket_size <= BA.bucket_write[i]) {
        fprintf(stderr, "# Warning: hit end of bucket nb %d\n", i);
        BA.bucket_write[i]--;
    }
#endif
}

static inline void
rewind_bucket(bucket_array_t BA, const int i)
{
    BA.bucket_read[i] = BA.bucket_start[i];
}

static inline void
rewind_bucket_by_1 (bucket_array_t BA, const int i)
{
  if (BA.bucket_read[i] > BA.bucket_start[i])
    BA.bucket_read[i]--;
}

static inline bucket_update_t
get_next_bucket_update(bucket_array_t BA, const int i)
{
    bucket_update_t rep = *(BA.bucket_read[i])++;
#ifdef SAFE_BUCKETS
    if (BA.bucket_read[i] > BA.bucket_write[i]) {
        fprintf(stderr, "# Warning: reading too many updates in bucket nb %d\n", i);
        BA.bucket_read[i]--;
        return LAST_UPDATE;
    }
#endif
    return rep;
}

static inline bucket_update_t
get_kth_bucket_update(const bucket_array_t BA, const int i, const int k)
{
    bucket_update_t rep = (BA.bucket_start[i])[k];
#ifdef SAFE_BUCKETS
    if (BA.bucket_start[i] + k >= BA.bucket_write[i]) {
        fprintf(stderr, "# Warning: reading outside valid updates in bucket nb %d\n", i);
        return LAST_UPDATE;
    }
#endif
    return rep;
}

static inline int
nb_of_updates(const bucket_array_t BA, const int i)
{
    return (BA.bucket_write[i] - BA.bucket_start[i]);
}

static inline void
push_sentinel(bucket_array_t BA, const int i)
{
    push_bucket_update(BA, i, LAST_UPDATE);
}

static inline int
is_end(const bucket_array_t BA, const int i)
{
    return (BA.bucket_read[i] == BA.bucket_write[i]);
}


static inline bucket_prime_t
get_next_bucket_prime (bucket_primes_t *BP)
{
  return *BP->read++;
}

static inline void
push_bucket_prime (bucket_primes_t *BP, const bucket_prime_t p)
{
  *BP->write++ = p;
}
                   

static inline int
bucket_primes_is_end(const bucket_primes_t *BP)
{
  return (BP->read == BP->write);
}

static inline void
rewind_primes_by_1 (bucket_primes_t *BP)
{
  if (BP->read > BP->start)
    BP->read--;
}


/* Remove some redundancy form the stored primes, e.g., remove the low
   bit which is always 1, or if BUCKET_ENCODE3 is set, store floor(p/6)*2 and
   the LSB telling whether it was 1 or 5 (mod 6). */
static inline prime_hint_t
bucket_encode_prime (uint32_t p)
{
  if (sizeof(prime_hint_t) == sizeof(uint32_t))
    return (prime_hint_t) p;
#ifdef BUCKET_ENCODE3
  return (prime_hint_t)(p/3U); /* This happens to work: if p=6k+1, we store
                                2k; if p=6k+5, we store 2k+1 */
#else
  return (prime_hint_t)(p/2U);
#endif
}

static inline uint32_t
bucket_decode_prime (prime_hint_t h)
{
  if (sizeof(prime_hint_t) == sizeof(uint32_t))
    return (uint32_t) h;
#ifdef BUCKET_ENCODE3
  return 3U * (uint32_t) h + 1U + ((uint32_t) h & 1U);
  /* if p was 6k+1, we stored h=2k, thus we want 3h+1;
     if p was 6k+5, we stored h=2k+1, thus we want 3h+2 */
#else
  return (((uint32_t) h) << 1) + 1U;
#endif
}
#endif	/* BUCKET_H_ */
