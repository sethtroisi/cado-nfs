#include "cado.h"
#include <iterator>
#include "las-types.hpp"

template<typename T>
unsigned long
append_prime_list (T inserter, prime_info pi, unsigned long pmax, cxx_mpz_poly const & f, int minroots = 1)
{
    unsigned long p;
    ASSERT_ALWAYS(minroots <= f->deg);
    if (f->deg == 1) {
        for (; (p = getprime_mt (pi)) < pmax; )
            *inserter++ = p;
    } else {
        for (; (p = getprime_mt (pi)) < pmax; )
            if (mpz_divisible_ui_p (f->coeff[f->deg], p) ||
                    mpz_poly_roots_ulong (NULL, f, p) >= minroots)
                *inserter++ = p;
    }
    return p;
}

void sieve_info::init_trialdiv(int side)
{
    sieve_info & si(*this);

    /* Our trial division needs odd divisors, 2 is handled by mpz_even_p().
       If the FB primes to trial divide contain 2, we skip over it.
       We assume that if 2 is in the list, it is the first list entry,
       and that it appears at most once. */

    sieve_info::side_info & s(si.sides[side]);
    unsigned long pmax = MIN((unsigned long) s.fbK.thresholds[0],
                             trialdiv_get_max_p());
    std::vector<unsigned long> trialdiv_primes = s.fbs->small_sieve_entries.skipped;

    /* Maybe we can use the factor base. If we have one, of course ! */
    unsigned long pmax_sofar = 0;
    if (s.fbs) {
        for(auto const & pp : s.fbs->small_sieve_entries.rest) {
            if (pp.k > 1) continue;
            trialdiv_primes.push_back(pp.p);
        }
        cxx_mpz zz(trialdiv_primes.back());
        mpz_nextprime(zz, zz);
        pmax_sofar = MIN(pmax, mpz_get_ui(zz));
    }
    if (pmax_sofar < pmax) {
        /* we need some more. */
        prime_info pi;
        prime_info_init(pi);
        unsigned long p;
        /* first seek to the end of the fb. */
        for ( ; (p = getprime_mt (pi)) < pmax_sofar ; );
        cxx_mpz_poly const & f(si.cpoly()->pols[side]);
        for(int minroots = 1 ; minroots <= f->deg ; minroots++) {
            p = append_prime_list(std::back_inserter(trialdiv_primes),
                    pi, MIN(pmax, minroots * s.fbK.td_thresh), f, minroots);
        }
        prime_info_clear (pi);
    }

    ASSERT(std::is_sorted(trialdiv_primes.begin(), trialdiv_primes.end()));
    // std::sort(trialdiv_primes.begin(), trialdiv_primes.end());
    
    size_t n = trialdiv_primes.size();
    int skip2 = n > 0 && trialdiv_primes[0] == 2;
    s.trialdiv_data = std::shared_ptr<trialdiv_divisor_t>(trialdiv_init(&trialdiv_primes.front() + skip2, n - skip2), trialdiv_clear);
}
