/*
 * Some bucket stuff.
 */

#include <stdlib.h>   // for malloc and friends
#include <unistd.h>   // for getpagesize
#include <stdint.h>

//#define SAFE_BUCKETS

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

/* In principal, the typedef for the bucket_update can be changed without
 * affecting the rest of the code. 
 */

/*
 * For the moment, we will keep the bucket updates aligned by adding an
 * 8-bit field that can contain, for instance, the low bits of p.
 *
 * TODO:
 * If the memory pressure becomes too high with this, we can remove this
 * p_low field and pack the updates as follows:
 *    [ x0 ] [ x1 ] [ x2 ] [ x3 ] [logp0] [logp1] [logp2] [logp3]
 * in the bucket.
 *
 * This will be slightly more tricky to store/load bucket updates, but
 * the additional cost should be negligible.
 */

typedef struct {
    uint16_t x;
    char logp;
    uint8_t p_low;  // to keep alignment and help trial div.
} bucket_update_t;

/*
 * will be used as a sentinel
 */
static const bucket_update_t LAST_UPDATE = {0,0,0};

/******** Bucket array typedef **************/

/*
 * A bucket is just an array of bucket_update_t.
 * We are going to manipulate several of them simultaneously when
 * sieving. In order to reduce the cache pressure, we are going to
 * introduce another structure that handles an array of bucket
 * efficiently.
 *
 * When doing a (long) sequence of pushes, only the bucket_start array is
 * used: pointer to it should be in a register, and the set of pointer
 * will probably fit in one page-size (please align!), so that the pressure
 * on the TLB is minimal. The other parts should not really be used
 * during this memory intensive part.
 * (except if one want to check for overflow, which can probably be
 * avoided by having large enough bucket_size).
 */

typedef struct {
    bucket_update_t ** bucket_start;    // Contains pointers to begining of
                                        // buckets.
    bucket_update_t ** bucket_write;    // Contains pointers to first empty
                                        // location in each bucket.
    bucket_update_t ** bucket_read;     // Contains pointers to first unread
                                        // location in each bucket.
    int bucket_size;                    // The allocated size of one bucket.
    int n_bucket;                       // Number of buckets.
} bucket_array_t;

/*
 * Notes for future improvements:
 *
 * 1) Double buckets.
 * If we can not afford enough buckets (limited by TLB) and we don't want
 * to make the bucket_region larger than L1, we can have a double bucket
 * system: A first level of buckets, that correspond to very large
 * bucket_regions, and a second level of buckets fitting in L1.
 * We collect all updates in the buckets of the first level, then we
 * process one first level bucket at a time, applying its update to the
 * corresponding second level buckets, and then applying to the L1
 * bucket_regions.
 *
 * 2) Buffered buckets.
 * One can have a buffer that contains one cache-line for each bucket.
 * Once a cache-line is full of updates, it can be moved in a
 * write-combined (non-temporal) way to the memory location of the
 * bucket. Details are left to the reader ;-)
 *
 */


/******** MAIN FUNCTIONS **************/


/* Returns an allocated array of <n_bucket> buckets each having size
 * <bucket_size>. Must be freed with clear_bucket_array().
 * It also put pointers in position ready for read/writes.
 */
static inline bucket_array_t
init_bucket_array(const int n_bucket, const int bucket_size);

static inline void
clear_bucket_array(bucket_array_t BA);

/* Main writing function: appends update to bucket number i.
 * If SAFE_BUCKETS is not #defined, then there is no checking that there is
 * enough room for the update. This could lead to a segfault, with the
 * current implementation!
 */
static inline void
push_bucket_update(bucket_array_t BA, int i, bucket_update_t update);

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
get_next_bucket_update(bucket_array_t BA, int i);
static inline int
nb_of_updates(bucket_array_t BA, int i);
static inline void
push_sentinel(bucket_array_t BA, int i);
static inline int
is_end(bucket_array_t BA, int i);

/* If you need to read twice a bucket, you can rewind it: */
static inline void
rewind_bucket(bucket_array_t BA, int i);

/* If you want to access updates in a non-sequential way: */
static inline bucket_update_t
get_kth_bucket_update(bucket_array_t BA, int i, int k);



/******** Bucket array implementation **************/

static inline void * 
malloc_check(size_t x) {
    void *p;
    p = malloc(x);
    ASSERT_ALWAYS(p != NULL);
    return p;
}

static inline void *
malloc_pagealigned(size_t x) {
    int sz = getpagesize();
    void *p;
    int ret = posix_memalign(&p, sz, x);
    ASSERT_ALWAYS(ret == 0);
    return p;
}

static inline bucket_array_t
init_bucket_array(const int n_bucket, const int bucket_size)
{
    bucket_array_t BA;
    int i;
    BA.bucket_size = bucket_size;
    BA.n_bucket = n_bucket;

    BA.bucket_start = (bucket_update_t **)malloc_pagealigned(n_bucket*sizeof(bucket_update_t *));
    BA.bucket_write = (bucket_update_t **)malloc_check(n_bucket*sizeof(bucket_update_t *));
    BA.bucket_read  = (bucket_update_t **)malloc_check(n_bucket*sizeof(bucket_update_t *));

    for (i = 0; i < n_bucket; ++i) {
        // TODO: shall we ensure here that those pointer do not differ by
        // a large power of 2, to take into account the associativity of
        // L1 cache ?
        // TODO: would it be better to have a single big malloc for all
        // the bucket_start[i] ?
        BA.bucket_start[i] = (bucket_update_t *)malloc_check(bucket_size*sizeof(bucket_update_t));
        BA.bucket_write[i] = BA.bucket_start[i];
        BA.bucket_read[i] = BA.bucket_start[i];
    }
    return BA;
}

static inline void
clear_bucket_array(bucket_array_t BA)
{
    int i;
    for (i = 0; i < BA.n_bucket; ++i)
        free(BA.bucket_start[i]);
    free(BA.bucket_start);
    free(BA.bucket_write);
    free(BA.bucket_read);
}

static inline void
push_bucket_update(bucket_array_t BA, int i, bucket_update_t update)
{
    *(BA.bucket_write[i])++ = update;
#ifdef SAFE_BUCKETS
    if (BA.bucket_start[i] + BA.bucket_size <= BA.bucket_write[i]) {
        fprintf(stderr, "# Warning: hit end of bucket nb %d\n", i);
        BA.bucket_write[i]--;
    }
#endif
}

static inline void
rewind_bucket(bucket_array_t BA, int i)
{
    BA.bucket_read[i] = BA.bucket_start[i];
}

static inline bucket_update_t
get_next_bucket_update(bucket_array_t BA, int i)
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
get_kth_bucket_update(bucket_array_t BA, int i, int k)
{
    bucket_update_t rep = (BA.bucket_start[i])[k];
#ifdef SAFE_BUCKETS
    if (k >= BA.bucket_write[i]) {
        fprintf(stderr, "# Warning: reading outside valid updates in bucket nb %d\n", i);
        return LAST_UPDATE;
    }
#endif
    return rep;
}

static inline int
nb_of_updates(bucket_array_t BA, int i)
{
    return (BA.bucket_write[i] - BA.bucket_start[i]);
}

static inline void
push_sentinel(bucket_array_t BA, int i)
{
    push_bucket_update(BA, i, LAST_UPDATE);
}

static inline int
is_end(bucket_array_t BA, int i)
{
    return (BA.bucket_read[i] == BA.bucket_write[i]);
}

