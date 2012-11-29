#ifndef __FQPOL_H__
#define __FQPOL_H__

#include "fq.h"

// Polynomials over Fq[x], represented as polynomials in x whose
// coefficients lie in Fq (using the fq_t type).
// By convention, deg(0) = -1.
typedef struct {
  int       deg;
  unsigned  alloc;
  fq_t  *coeffs;
} __fqpol_struct;

typedef       __fqpol_struct  fqpol_t[1];
typedef       __fqpol_struct *fqpol_ptr;
typedef const __fqpol_struct *fqpol_srcptr;


void fqpol_init(fqpol_ptr r);
void fqpol_clear(fqpol_ptr r);
int fqpol_deg(fqpol_srcptr r, MAYBE_UNUSED fq_info_srcptr Fq);
int fqpol_is_zero(fqpol_srcptr r, MAYBE_UNUSED fq_info_srcptr Fq);
void fqpol_set_zero(fqpol_ptr r, MAYBE_UNUSED fq_info_srcptr Fq);
void fqpol_set_ffspol(fqpol_ptr r, ffspol_srcptr f, fq_info_srcptr Fq);
void fqpol_set(fqpol_ptr r, fqpol_srcptr p, fq_info_srcptr Fq);
void fqpol_set_ti(fqpol_ptr r, unsigned int n, fq_info_srcptr Fq);
void fqpol_add(fqpol_ptr r, fqpol_srcptr p, fqpol_srcptr q, fq_info_srcptr Fq);
void fqpol_sub(fqpol_ptr r, fqpol_srcptr p, fqpol_srcptr q, fq_info_srcptr Fq);
void fqpol_mul(fqpol_ptr r, fqpol_srcptr p, fqpol_srcptr q, fq_info_srcptr Fq);
void fqpol_sqr(fqpol_ptr r, fqpol_srcptr p, fq_info_srcptr Fq);
void fqpol_smul(fqpol_ptr r, fqpol_srcptr p, fq_srcptr c, fq_info_srcptr Fq);
void fqpol_sdiv(fqpol_ptr r, fqpol_srcptr p, fq_srcptr c, fq_info_srcptr Fq);
void fqpol_mul_ti(fqpol_ptr r, fqpol_srcptr p, unsigned i, fq_info_srcptr Fq);
void fqpol_get_coeff(fq_ptr r, fqpol_srcptr p, unsigned i, fq_info_srcptr Fq);
void fqpol_set_coeff(fqpol_ptr r, fq_srcptr x, unsigned i, fq_info_srcptr Fq);
int fqpol_divrem(fqpol_ptr q, fqpol_ptr r,
                fqpol_srcptr a, fqpol_srcptr b, fq_info_srcptr Fq);
int fqpol_div(fqpol_ptr q, fqpol_srcptr a, fqpol_srcptr b, fq_info_srcptr Fq);
int fqpol_rem(fqpol_ptr r, fqpol_srcptr a, fqpol_srcptr b, fq_info_srcptr Fq);
void fqpol_gcd(fqpol_ptr r, fqpol_srcptr p, fqpol_srcptr q, fq_info_srcptr Fq);

// roots must be allocated by caller.
// return the number of (distinct) roots. 
int fqpol_roots(fq_t * roots, fqpol_srcptr f, fq_info_srcptr Fq);
int fqpol_is_split(fq_t * roots, fqpol_srcptr f, fq_info_srcptr Fq);
int fqpol_is_irreducible(fqpol_srcptr f, fq_info_srcptr Fq);
int fqpol_root_multiplicity(fqpol_srcptr f, fq_srcptr root, fq_info_srcptr Fq);

void fqpol_out(FILE *file, fqpol_srcptr p, fq_info_srcptr Fq);


#endif   /* __FQPOL_H__ */
