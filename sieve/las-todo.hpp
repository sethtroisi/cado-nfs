#ifndef LAS_TODO_HPP_
#define LAS_TODO_HPP_

#include <gmp.h>
#include "las-base.hpp"
#include "cxx_mpz.hpp"

struct las_todo_entry {
  cxx_mpz p; /* this is the 'special-q', despite the 'p' name... */
    /* even for a rational side, the field below is used, since
     * it is needed for the initialization of the q-lattice. All callers
     * of las_todo_push must therefore make sure that a proper argument
     * is provided.
     */
    cxx_mpz r;
    int side;
    bool prime_sq;
    std::vector<uint64_t> prime_factors;

    /* some fields which are specific to the descent */
    int depth;
    int iteration;      /* number of times we failed on this prime */

    /* Default ctor is just to enable sieve_info default ctor */
    las_todo_entry() : side(0), depth(0), iteration(0) {}
    /* Empty p, r is used for "closing brace" */
    las_todo_entry(const int _side, const int _depth) : side(_side), depth(_depth), iteration(0) { }
    las_todo_entry(const mpz_t _p, const mpz_t _r, const int _side, const int _depth = 0, const int _iteration = 0) {
      set(_p, _r, _side, _depth, _iteration);
    }

    // Given a *prime* ell, check whether ell is coprime to current
    // entry.
    bool is_coprime_to(unsigned long ell) const {
        if (prime_sq)
            return (mpz_cmp_ui(p, ell) != 0);
        else {
            for (unsigned int i = 0; i < prime_factors.size(); ++i) {
                if (prime_factors[i] == ell)
                    return false;
            }
            return true;
        }
    }
private:
  // TODO: this function becomes a bit long for a header file.
  // But for the moment we don't have a las-todo.cpp, so...
  void set(const mpz_t _p, const mpz_t _r, const int _side, const int _depth, const int _iteration) {
      mpz_set(p, _p);
      mpz_set(r, _r);
      side = _side;
      depth = _depth;
      iteration = _iteration;
      prime_sq = mpz_probab_prime_p(p, 2);
      if (prime_sq) return;

      // Need to pre-compute the prime factors of the special-q in order to
      // skip them while sieving.
      prime_info pi;
      prime_info_init (pi);
      unsigned long f = 2;
      mpz_t B;
      mpz_init(B);
      mpz_sqrt(B, p);
      unsigned long bound = mpz_get_ui(B);
      mpz_clear(B);
      mpz_t pp;
      mpz_init_set(pp, p);
      while ((f <= bound) && (mpz_cmp_ui(pp, 1) > 0)) {
          if (mpz_divisible_ui_p(pp, f)) {
              mpz_divexact_ui(pp, pp, f);
              // Powers are not allowed in special-q
              ASSERT_ALWAYS(!mpz_divisible_ui_p(pp, f));
              prime_factors.push_back(f);
              if (mpz_probab_prime_p(pp, 2)) {
                  ASSERT_ALWAYS(mpz_fits_ulong_p(pp));
                  prime_factors.push_back(mpz_get_ui(pp));
                  mpz_set_ui(pp, 1);
              }
          }
          f = getprime_mt(pi);
      }
      prime_info_clear (pi);
      mpz_clear(pp);
  }
};

#endif
