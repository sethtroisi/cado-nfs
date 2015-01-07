#ifndef POLYSELECT_AUXILIARY_H_
#define POLYSELECT_AUXILIARY_H_

/* header file for auxiliary routines for polyselect

Copyright 2008, 2009, 2010, 2013 Emmanuel Thome, Paul Zimmermann

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <float.h>  /* for DBL_MAX */
#include <string.h>

#define DEFAULT_INCR 60 /* we want a positive integer with many divisors,
                           other values are 210, 2310, 30030, 510510, 9699690,
                           223092870 */

#define SKEWNESS_DEFAULT_PREC 10

/* prime bounds for the computation of alpha */
#define ALPHA_BOUND_SMALL  100
#define ALPHA_BOUND       2000

#define MURPHY_K 1000

extern double bound_f, bound_g, area;

#define mpz_add_si(a,b,c)                       \
  if (c >= 0) mpz_add_ui (a, b, c);             \
  else mpz_sub_ui (a, b, -(c))

#define mpz_submul_si(a,b,c)                    \
  if (c >= 0) mpz_submul_ui (a, b, c);          \
  else mpz_addmul_ui (a, b, -(c))
  
#include "cado_poly.h"
#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

double L2_lognorm (mpz_poly_ptr, double);
double L2_skewness (mpz_poly_ptr, int);
double L2_skew_lognorm (mpz_poly_ptr, int);

/* alpha */
double special_val0 (mpz_poly_ptr, unsigned long);
double special_valuation (mpz_poly_ptr f, unsigned long p, mpz_t disc);
double special_valuation_affine (mpz_poly_ptr f, unsigned long p, mpz_t disc);
double get_alpha (mpz_poly_ptr, unsigned long);
double get_biased_alpha_projective (mpz_poly_ptr f, unsigned long B);

/* poly info, being called in order */
void print_cadopoly_fg (FILE*, mpz_t*, int, mpz_t*, int, mpz_t);
double print_cadopoly (FILE*, cado_poly);
void print_cadopoly_extra (FILE*, cado_poly, int, char**, double);
double print_poly_fg (mpz_poly_ptr, mpz_t*, mpz_t, int);
void discriminant (mpz_t, mpz_t*, const int);
long rotate_aux (mpz_t *f, mpz_t b, mpz_t m, long k0, long k, unsigned int t);
void rotate_auxg_z (mpz_t*, const mpz_t, const mpz_t, const mpz_t, unsigned int);
void do_translate_z (mpz_poly_ptr f, mpz_t *g, const mpz_t k);


void cado_poly_fprintf_with_info (FILE *, cado_poly_ptr, const char *);
void cado_poly_fprintf_with_info_and_MurphyE (FILE *fp, cado_poly_ptr, double,
                                              double, double, double,
                                              const char *);

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_AUXILIARY_H_ */

