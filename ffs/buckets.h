#ifndef __BUCKETS_H__
#define __BUCKETS_H__

#include <stdint.h>

#include "types.h"
#include "fb.h"



/* Array of buckets.
 *****************************************************************************/

// Forward declaration of bucket updates.
typedef struct __update_packed_struct update_packed_t;

typedef struct {
  unsigned          n;
  unsigned          max_size;
  unsigned          min_degp, max_degp;
  update_packed_t **start;
  update_packed_t **degp_end;
} __buckets_struct;

typedef       __buckets_struct  buckets_t[1];
typedef       __buckets_struct *buckets_ptr;
typedef const __buckets_struct *buckets_srcptr;

// Initialize structure and allocate buckets.
void buckets_init(buckets_ptr buckets, unsigned I, unsigned J,
                  unsigned max_size, unsigned min_degp, unsigned max_degp);

// Clean up memory.
void buckets_clear(buckets_ptr buckets);

// Return the size of a bucket region.
unsigned bucket_region_size();

// Print information about the buckets.
void print_bucket_info(buckets_srcptr buckets);

// Fill the buckets with updates corresponding to divisibility by elements of
// the factor base.
void buckets_fill(buckets_ptr buckets, factor_base_srcptr FB,
                  sublat_srcptr sublat, unsigned I, unsigned J);

// Apply all the updates from a given bucket to the sieve region S.
void bucket_apply(uint8_t *S, buckets_srcptr buckets, unsigned k);

#endif  /* __BUCKETS_H__ */
