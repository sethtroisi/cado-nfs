#ifndef ANT_HPP_
#define ANT_HPP_

#include "mpz_poly.h"
#include "mpz_mat.h"
#include "cxx_mpz.hpp"

#include <vector>
#include <utility>

/* Return a basis for a p-maximal order of the number field defined by
 * the polynomial f.
 *
 * The basis of the order is written in lower triangular positive row hnf
 * (magma returns lower triangular centered row hnf, with sometimes negative
 * coefficients).
 */

cxx_mpq_mat p_maximal_order(cxx_mpz_poly const& f, unsigned long p);

std::vector<std::pair<cxx_mpz_mat, int>> factorization_of_prime(cxx_mpq_mat & B, cxx_mpz_poly const& g, unsigned long p, gmp_randstate_t state);

struct ideal_comparator {
    typedef std::pair<cxx_mpz_mat,int> Im_t;
    bool operator()(Im_t const& a, Im_t const& b) const {
        int r = mpz_mat_cmp(a.first, b.first);
        if (r) return r < 0;
        return a.second < b.second;     /* should never happen */
    }
};

#endif	/* ANT_HPP_ */
