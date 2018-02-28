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

/* {{{ sieve_info::{init,clear}_fb_smallsieved */
void sieve_info::init_fb_smallsieved(int side)
{
    /* TODO: This is just a copy. We should kill that, and use the field
     * from the slicing directly. */
    sieve_info & si(*this);
    if (!si.sides[side].fb) return;
    const fb_factorbase::slicing * fbs = si.sides[side].fbs;
    /* put the resieved primes first. */
    auto s = std::make_shared<std::vector<fb_entry_general>>();
    std::vector<fb_entry_general> const & RS(fbs->small_sieve_entries.resieved);
    s->insert(s->end(), RS.begin(), RS.end());
    si.sides[side].resieve_start_offset = 0;
    si.sides[side].resieve_end_offset = s->size();
    std::vector<fb_entry_general> const & R(fbs->small_sieve_entries.rest);
    s->insert(s->end(), R.begin(), R.end());
    si.sides[side].fb_smallsieved = s;
}

/* }}} */

