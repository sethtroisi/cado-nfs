#include "cado.h"
#include <stdarg.h>
#include <gmp.h>
#include "las-types.hpp"
#include "fb.hpp"

/*  reordering of the small factor base
 *
 * We split the small factor base in several non-overlapping, contiguous
 * zones:
 *
 *      - powers of 2 (up until the pattern sieve limit)
 *      - powers of 3 (up until the pattern sieve limit)
 *      - trialdiv primes (not powers)
 *      - resieved primes
 *      (- powers of trialdiv primes)
 *      - rest.
 *
 * Problem: bad primes may in fact be pattern sieved, and we might want
 * to pattern-sieve more than just the ``do it always'' cases where p is
 * below the pattern sieve limit.
 *
 * The answer to this is that such primes are expected to be very very
 * rare, so we don't really bother. If we were to do something, we could
 * imagine setting up a schedule list for projective primes -- e.g. a
 * priority queue. But it feels way overkill.
 *
 * Note that the pre-treatment (splitting the factor base in chunks) can
 * be done once and for all.
 */

#if 0
static size_t count_roots(const std::vector<fb_entry_general> &v)
{
    size_t count = 0;
    for(std::vector<fb_entry_general>::const_iterator it = v.begin();
        it != v.end(); it++) {
        count += it->nr_roots;
    }
    return count;
}
#endif

