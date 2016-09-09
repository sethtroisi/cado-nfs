#ifndef ANT_HPP_
#define ANT_HPP_

#include "mpz_poly.h"
#include "mpz_mat.h"
#include "cxx_mpz.hpp"

/* Return a basis for a p-maximal order of the number field defined by
 * the polynomial f.
 *
 * TODO: explain the normalization choices for the matrix.
 */

cxx_mpq_mat p_maximal_order(cxx_mpz_poly const& f, unsigned long p);

#endif	/* ANT_HPP_ */
