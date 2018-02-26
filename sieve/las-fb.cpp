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
    sieve_info & si(*this);

    if (!si.sides[side].fb) return;

    /* We go through all the primes in FB part 0 and sort them into one of
       5 vectors, and some we discard */

    /* FIXME: that static separation between FB parts is artificial, and
     * actually causes problems. We want to get rid of it, because our
     * small sieve limit depends on logI, and logI depends on the special q.
     */

    /* alloc the 6 vectors */
    enum {POW2, POW3, TD, RS, REST, SKIPPED};
    std::vector<std::vector<fb_entry_general>> pieces(6);

    fb_part *small_part = si.sides[side].fb->get_part(0);
    ASSERT_ALWAYS(small_part->is_only_general());
    const fb_vector<fb_entry_general> *small_entries = small_part->get_general_vector();

    fbprime_t plim = si.conf.bucket_thresh;
    fbprime_t costlim = si.conf.td_thresh;

    const size_t pattern2_size = sizeof(unsigned long) * 2;
    for (fb_vector<fb_entry_general>::const_iterator it = small_entries->begin();
         it != small_entries->end();
         it++)
    {
        /* The extra conditions on powers of 2 and 3 are related to how
         * pattern-sieving is done.
         *
         * FIXME: there is some duplicated logic between here and
         * small_sieve_init.
         */
        if (it->q < si.conf.skipped) {
            pieces[SKIPPED].push_back(*it);
        } else if (it->p == 2 && it->q <= pattern2_size) {
            pieces[POW2].push_back(*it);
        } else if (it->q == 3) { /* Currently only q=3 is pattern sieved */
            pieces[POW3].push_back(*it);
        } else if (it->k != 1) {
            /* prime powers always go into "rest" */
            pieces[REST].push_back(*it);
        } else if (it->q <= plim && it->q <= costlim * it->nr_roots) {
            pieces[TD].push_back(*it);
        } else if (it->q <= plim) {
            pieces[RS].push_back(*it);
        } else {
            abort();
        }
    }

    /* put the resieved primes first. */
    auto s = std::make_shared<std::vector<fb_entry_general>>(pieces[RS]);
    si.sides[side].resieve_start_offset = 0;
    si.sides[side].resieve_end_offset = s->size();

    /* now the rest. The small_sieve ctor will filter out the exceptional
     * cases (proj, powers, pattern sieve, etc), and sort the remaining
     * primes (which may or may not be already sorted, depending on the
     * shape of the polynomial -- we're getting them as per the slice
     * ordering, which may be a bit messy).
     */
    for (int i = 0; i < 6; i++) {
        if ((i == RS) || (i == SKIPPED))
            continue;
        s->insert(s->end(), pieces[i].begin(), pieces[i].end());
    }
    si.sides[side].fb_smallsieved = s;
}

/* }}} */

/* {{{ Print some statistics about the factor bases
 * This also fills the field si.sides[*]->max_bucket_fill_ratio, which
 * is used to verify that per-thread allocation for buckets is
 * sufficient.
 * */
void sieve_info::print_fb_statistics(int side)
{
    side_info & s(sides[side]);
    if (!s.fb) return;
    for (int i_part = 0; i_part < FB_MAX_PARTS; i_part++)
    {
        size_t nr_primes, nr_roots;
        double weight;
        s.fb->get_part(i_part)->count_entries(&nr_primes, &nr_roots, &weight);
        if (nr_primes != 0 || weight != 0.) {
            verbose_output_print(0, 1, "# Number of primes in side-%d factor base part %d = %zu\n",
                    side, i_part, nr_primes);
            verbose_output_print(0, 1, "# Number of prime ideals in side-%d factor base part %d = %zu\n",
                    side, i_part, nr_roots);
            verbose_output_print(0, 1, "# Weight of primes in side-%d factor base part %d = %0.5g\n",
                    side, i_part, weight);
        }
        s.max_bucket_fill_ratio[i_part] = weight;
    }
}
/*}}}*/



