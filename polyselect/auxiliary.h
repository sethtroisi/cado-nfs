/* header file for auxiliary routines for polyselect

Copyright 2008, 2009, 2010 Emmanuel Thome, Paul Zimmermann

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

#define MAX_DEGREE 6

#if 1 /* use L2-norm with integral over whole sieving region */
#define LOGNORM  L2_lognorm
#define SKEWNESS L2_skewness
#else /* use 1-norm relative to coefficients */
#define LOGNORM  l1_lognorm
#define SKEWNESS l1_skewness
#endif

#define SKEWNESS_DEFAULT_PREC 10

/* prime bounds for the computation of alpha */
#define ALPHA_BOUND_SMALL  100
#define ALPHA_BOUND       2000

#define mpz_add_si(a,b,c)                       \
  if (c >= 0) mpz_add_ui (a, b, c);             \
  else mpz_sub_ui (a, b, -(c))

#define mpz_submul_si(a,b,c)                    \
  if (c >= 0) mpz_submul_ui (a, b, c);          \
  else mpz_addmul_ui (a, b, -(c))
  
void mpz_ndiv_qr (mpz_t, mpz_t, mpz_t, mpz_t);
void generate_base_mb (cado_poly, mpz_t, mpz_t);
double l1_skewness (mpz_t*, int, int);
double l1_lognorm (mpz_t*, unsigned long, double);
double L2_lognorm (mpz_t*, unsigned long, double);
double L2_skewness (mpz_t*, int, int);
/* rotation */
double get_alpha (mpz_t*, const int, unsigned long);
void discriminant (mpz_t, mpz_t*, const int);
long rotate_aux (mpz_t *f, mpz_t b, mpz_t m, long k0, long k);
long rotate_aux1 (mpz_t *f, mpz_t b, mpz_t m, long j0, long j);
double rotate (mpz_t*, int, unsigned long, mpz_t, mpz_t, long*, long*, int,int);
void print_poly (FILE*, cado_poly, int, char**, double, int);
long translate (mpz_t*, int, mpz_t*, mpz_t, mpz_t, int);
void optimize (mpz_t*, int, mpz_t*, int);
void rotate_bounds (mpz_t *f, int d, mpz_t b, mpz_t m, long *K0, long *K1, long *J0, long *J1, int verbose);

/********************* data structures for first phase ***********************/

typedef struct {
  /* the linear polynomial is b*x-m, with lognorm logmu */
  mpz_t b;
  mpz_t m;
  double logmu;
} m_logmu_t;

m_logmu_t* m_logmu_init (unsigned long);
void m_logmu_clear (m_logmu_t*, unsigned long);
int m_logmu_insert (m_logmu_t*, unsigned long, unsigned long*, mpz_t, mpz_t,
                    double, char*);


