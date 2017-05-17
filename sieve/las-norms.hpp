#ifndef LAS_NORMS_HPP_
#define LAS_NORMS_HPP_

#include <stdint.h>
#include "las-types.hpp"
#include "double_poly.h"


/* Initialize lognorms for the bucket_region number J. It's a wrapper.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime.
 * The sieve area S must be preallocated with at least (BUCKET_REGION +
 * MEMSET_MIN) space. Indeed, the algorithm that is used might write a
 * bit beyond the meaningful area.
 */
void init_norms_bucket_region (unsigned char *S, uint32_t J, sieve_info& si, unsigned int side, unsigned int smart);

double get_maxnorm_alg (double_poly_srcptr src_poly, const double X, const double Y);

void sieve_info_init_norm_data_sq (sieve_info& si, unsigned long q);

  
/* These functions are internals. Don't use them. Use the wrapper above.
   It's need to declare them here for units & coverage tests.
 */
void init_degree_one_norms_bucket_region_internal     (unsigned char *S, uint32_t J, uint32_t I, double scale, cxx_double_poly const &, double *cexp2);
void init_exact_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, cxx_double_poly const & fijd);
void init_smart_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, cxx_double_poly const & fijd, std::vector<smart_norm_root> const & roots);
void init_norms_roots_internal (cxx_double_poly const &, double max_abs_root, double precision, std::vector<smart_norm_root> & roots);
void init_norms_roots (sieve_info & si, unsigned int side);

#endif	/* LAS_NORMS_HPP_ */
