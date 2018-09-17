#ifndef COFAC_LIST_H
#define COFAC_LIST_H

#include <stdlib.h>
#include <inttypes.h>
#include <gmp.h>
#include <list>
#include <array>
#include "cxx_mpz.hpp"
#include "facul.hpp"
#include "facul_doit.hpp"
#include "utils.h"
#include "relation.hpp"

/* structure to compute on-line a product tree, avoiding to first compute a
   list of mpz_t (which might take too much memory) */
typedef struct {
  mpz_t *l;     /* the value stored is l[0] * l[1] * ... * l[size-1],
                   where l[0] is the product of n[0] elements, l[1] is
                   the product of n[1] elements, ..., with n[0]=0 or 1,
                   n[1]=0 or 2, ..., n[k]=0 or 2^k */
  unsigned long *n;
  size_t size;
} mpz_product_tree_t;
typedef mpz_product_tree_t mpz_product_tree[1];

struct cofac_candidate {
  int64_t a;
  uint64_t b;
  std::array<cxx_mpz, 2> cofactor;
  cxx_mpz sq;   /* special-q */
  int sqside;   /* most often n times the same thing.
                   Have to do that since otherwise the las todo mode breaks.
                   There ought to be a better way, but this one's a no-brain.
                   Well, we have mpz's around anyway... */
  cofac_candidate() = default;
  cofac_candidate(int64_t a, uint64_t b, std::array<cxx_mpz,2> & cofactor, cxx_mpz const & sq, int sqside)
      : a(a), b(b), cofactor(std::move(cofactor)), sq(sq), sqside(sqside)
      {}
};

typedef std::list<cofac_candidate> cofac_list;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * These three functions add to to the double& extra_time argument the
 * cpu time (RUSAGE_THREAD, seconds_thread()) spent in openmp helper
 * threads, NOT counting the time spent in the main thread.
 */
size_t find_smooth (
        cofac_list & l,
        std::array<cxx_mpz, 2> & batchP,
        int batchlpb[2], int lpb[2], int batchmfb[2],
        FILE *out,
        int nthreads MAYBE_UNUSED, double &);

std::list<relation> factor (cofac_list const &, cado_poly_srcptr, int[2], int[2], FILE*, int, double&);
void create_batch_file (const char*, mpz_t, unsigned long, unsigned long,
                        mpz_poly, FILE*, int, double &);

#ifdef __cplusplus
}
#endif

#endif /* COFAC_LIST_H */
