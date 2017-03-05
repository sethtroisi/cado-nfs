#include "cado.h"
#include "las-types.hpp"

void sieve_info::init_trialdiv(int side)
{
    sieve_info & si(*this);
    /* Our trial division needs odd divisors, 2 is handled by mpz_even_p().
       If the FB primes to trial divide contain 2, we skip over it.
       We assume that if 2 is in the list, it is the first list entry,
       and that it appears at most once. */

    /* XXX This function consider the full contents of the list
     * si.sides[side].fb to be trialdiv primes, therefore that
     * sieve_info_split_bucket_fb_for_threads must have already been
     * called, so that the only primes which are still in s->fb are the
     * trial-divided ones */
    sieve_info::side_info & s(si.sides[side]);
    unsigned long pmax = MIN((unsigned long) si.conf.bucket_thresh,
                             trialdiv_get_max_p());
    std::vector<unsigned long> trialdiv_primes;
    s.fb->extract_bycost(trialdiv_primes, pmax, si.conf.td_thresh);
    std::sort(trialdiv_primes.begin(), trialdiv_primes.end());
    size_t n = trialdiv_primes.size();
    int skip2 = n > 0 && trialdiv_primes[0] == 2;
    s.trialdiv_data = std::shared_ptr<trialdiv_divisor_t>(trialdiv_init(&trialdiv_primes.front() + skip2, n - skip2), trialdiv_clear);
}
