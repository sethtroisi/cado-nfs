#ifndef LAS_TODO_H_
#define LAS_TODO_H_

#include <gmp.h>
#include <stack>
#include "las-base.hpp"

struct las_todo_entry {
    mpz_t p;
    /* even for a rational side, the field below is used, since
     * it is needed for the initialization of the q-lattice. All callers
     * of las_todo_push must therefore make sure that a proper argument
     * is provided.
     */
    mpz_t r;
    int side;
    bool prime_sq;
    std::vector<uint64_t> prime_factors;

    /* some fields which are specific to the descent */
    int depth;
    int iteration;      /* number of times we failed on this prime */

    /* FIXME: Default constructor should be removed! Currently in just to make
       sieve_info_s without a constructor work */
    las_todo_entry() : side(0), depth(0), iteration(0) {
      alloc();
    }
    /* Empty p, r is used for "closing brace" */
    las_todo_entry(const int _side, const int _depth) : side(_side), depth(_depth), iteration(0) {
      alloc();
    }
    las_todo_entry(const mpz_t _p, const mpz_t _r, const int _side, const int _depth = 0, const int _iteration = 0) {
      alloc();
      set(_p, _r, _side, _depth, _iteration);
    }
    las_todo_entry(const las_todo_entry &other) {
      /* Delegating constructors is C++11 :( */
      alloc();
      set(other.p, other.r, other.side, other.depth, other.iteration);
    }
    ~las_todo_entry() {
      mpz_clear(p);
      mpz_clear(r);
    }
    las_todo_entry& operator=(const las_todo_entry &other) {
      if (this != &other)
        set(other.p, other.r, other.side, other.depth, other.iteration);
      return *this;
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
  void alloc() {
      mpz_init(p);
      mpz_init(r);
  }
  // TODO: this function becomes a bit long for a header file.
  // But for the moment we don't have a las-todo.cpp, so...
  void set(const mpz_t _p, const mpz_t _r, const int _side, const int _depth, const int _iteration) {
      mpz_set(p, _p);
      mpz_set(r, _r);
      side = _side;
      depth = _depth;
      iteration = _iteration;
      prime_sq = bool(mpz_probab_prime_p(p, 2));
      if (!prime_sq) {
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
  }
};

typedef std::stack<las_todo_entry *> las_todo_stack;

#endif
