#ifndef __BUCKETS_H__
#define __BUCKETS_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>



/* Bucket updates.
 *
 * Each bucket update contains:
 * - pos:  the position (i, j) of the update, represented using an unsigned
 *         integer of BUCKET_UPDATE_POS_BITS bits.
 * - hint: the most significant coefficients of the prime p corresponding to
 *         the update, represented using an unsigned integer of
 *         BUCKET_UPDATE_HINT_BITS bits.
 *
 * The actual size in memory of a bucket update is BUCKET_UPDATE_POS_BITS +
 * BUCKET_UPDATE_HINT_BITS, rounded up to a multiple of 8 bits.
 *****************************************************************************/

// Size of the two fields and of a bucket update, in bits.
#define   BUCKET_UPDATE_POS_BITS  17
#define   BUCKET_UPDATE_HINT_BITS  7
#define __BUCKET_UPDATE_BITS     (BUCKET_UPDATE_POS_BITS + \
                                  BUCKET_UPDATE_HINT_BITS)

// In memory, a bucket update will just be an array of bytes.
typedef struct {
  uint8_t __d[(__BUCKET_UPDATE_BITS+7) / 8];
} bucket_update_t;

// Smallest unsigned integer in which a bucket update will fit.
#if   __BUCKET_UPDATE_BITS <= 16
typedef uint16_t __bucket_update_ui_t;
#elif __BUCKET_UPDATE_BITS <= 32
typedef uint32_t __bucket_update_ui_t;
#else
typedef uint64_t __bucket_update_ui_t;
#endif

// Internal type for accessing the fields of a bucket update.
typedef union {
    bucket_update_t    data;
  __bucket_update_ui_t ui;
} __bucket_update_t;

// Types for positions and hints as exported by the functions.
typedef unsigned pos_t;
typedef unsigned hint_t;

// Construct a bucket update from a position and a prime hint.
static inline
bucket_update_t bucket_update_set(pos_t pos, hint_t hint)
{ __bucket_update_t t = { .ui = pos | (hint << BUCKET_UPDATE_POS_BITS) };
  return t.data; }

// Retrieve the position from a bucket update.
static inline
pos_t  bucket_update_get_pos (bucket_update_t update)
{ __bucket_update_t t = { .data = update };
  return t.ui & (((pos_t)1<<BUCKET_UPDATE_POS_BITS)-1); }

// Retrieve the prime hint from a bucket update.
static inline
hint_t bucket_update_get_hint(bucket_update_t update)
{ __bucket_update_t t = { .data = update };
  return (t.ui >> BUCKET_UPDATE_POS_BITS) &
         (((hint_t)1<<BUCKET_UPDATE_HINT_BITS)-1); }



/* Array of buckets.
 *****************************************************************************/

typedef struct {
  unsigned          n;
  unsigned          size;
  bucket_update_t **begin;
  bucket_update_t **iter;
  bucket_update_t **end;
} __bucket_array_struct;

typedef       __bucket_array_struct  bucket_array_t[1];
typedef       __bucket_array_struct *bucket_array_ptr;
typedef const __bucket_array_struct *bucket_array_srcptr;

static inline
void bucket_array_init(bucket_array_ptr buckets, unsigned n, unsigned size)
{
  buckets->n     = n;
  buckets->size  = size;
  buckets->begin = (bucket_update_t **)malloc(n * sizeof(bucket_update_t *));
  buckets->iter  = (bucket_update_t **)malloc(n * sizeof(bucket_update_t *));
  buckets->end   = (bucket_update_t **)malloc(n * sizeof(bucket_update_t *));
  ASSERT_ALWAYS(buckets->begin != NULL);
  ASSERT_ALWAYS(buckets->iter  != NULL);
  ASSERT_ALWAYS(buckets->end   != NULL);
  for (unsigned i = 0; i < n; ++i) {
    buckets->begin[i] = buckets->iter[i] = buckets->end[i] = 
      (bucket_update_t *)malloc(size * sizeof(bucket_update_t));
    ASSERT_ALWAYS(buckets->begin[i] != NULL);
  }
}

static inline
void bucket_array_clear(bucket_array_ptr buckets)
{
  for (unsigned i = 0; i < buckets->n; ++i)
    free(buckets->begin[i]);
  free(buckets->begin);
  free(buckets->iter);
  free(buckets->end);
}

static inline
unsigned bucket_size(bucket_array_srcptr buckets, unsigned i)
{ return buckets->end[i] - buckets->begin[i]; }

static inline
int bucket_is_done(bucket_array_srcptr buckets, unsigned i)
{ return buckets->iter[i] == buckets->end[i]; }

static inline
void bucket_push(bucket_array_ptr buckets, unsigned i,
                 const bucket_update_t update)
{ *buckets->end[i]++ = update; }

static inline
bucket_update_t bucket_next(bucket_array_ptr buckets, unsigned i)
{ return *buckets->iter[i]++; }

#endif  /* __BUCKETS_H__ */
