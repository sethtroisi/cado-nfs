#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "int64_poly.h"
#include "sieving_interval.h"
#include "macros.h"
#include "utils_int64.h"

/* #include "utils/gcd.h" */
/* #include <gmp.h> */
/* #include "mpz_poly.h" */

/*
  --- Private functions. ---
*/

/* static inline int64_t int64_poly_lc(int64_poly_ptr f) */
/* { */
/*     assert(f->deg >= 0); */
/*     return f->coeff[f->deg]; */
/* } */

/* static inline int64_t int64_poly_lc_const(int64_poly_srcptr f) */
/* { */
/*     assert(f->deg >= 0); */
/*     return f->coeff[f->deg]; */
/* } */

/*
  --- Public functions. ---
*/

void int64_poly_init(int64_poly_ptr f, int d)
{
  f->deg = -1;
  if (d < 0) {
    f->alloc = 0;
    f->coeff = (int64_t *) NULL;
  } else {
    f->alloc = d + 1;
    f->coeff = malloc(sizeof(int64_t) * (d + 1));
    FATAL_ERROR_CHECK (f->coeff == NULL, "not enough memory");
    /* //WARNING: usage of memset. */
    /* memset(f->coeff, 0, sizeof(int)); */
    for (int i = 0; i <= d; ++i) {
      f->coeff[i] = 0;
    }
  }
}

void int64_poly_realloc(int64_poly_ptr f, int nc)
{
  if (f->alloc < nc)
  {
    f->coeff = (int64_t*) realloc (f->coeff, nc * sizeof (int64_t));
    FATAL_ERROR_CHECK (f->coeff == NULL, "not enough memory");
    /* //WARNING: usage of memset. */
    /* memset(f->coeff + f->alloc, 0, sizeof(int)); */
    for (int i = f->alloc; i < nc; i++) {
        f->coeff[i] = 0;
    }
    f->alloc = nc;
  }
}

void int64_poly_copy(int64_poly_ptr g, int64_poly_srcptr f)
{
  g->deg = f->deg;
  for (int i = f->deg; i >= 0; --i) {
    int64_poly_setcoeff(g, i, f->coeff[i]);
  }
}

void int64_poly_swap (int64_poly_ptr f, int64_poly_ptr g)
{
  int i;
  int64_t * t;

  i = f->alloc;
  f->alloc = g->alloc;
  g->alloc = i;
  i = f->deg;
  f->deg = g->deg;
  g->deg = i;
  t = f->coeff;
  f->coeff = g->coeff;
  g->coeff = t;
}

void int64_poly_clear(int64_poly_ptr f)
{
  free(f->coeff);
  memset(f, 0, sizeof(int64_poly_t));
  f->deg = -1;
}

void int64_poly_cleandeg(int64_poly_ptr f, int deg)
{
  ASSERT(deg >= -1);
  while ((deg >= 0) && (f->coeff[deg] == 0))
    deg--;
  f->deg = deg;
}

void int64_poly_set_zero(int64_poly_ptr f)
{
  f->deg = -1;
}

void int64_poly_setcoeff(int64_poly_ptr f, int i, int64_t coeff)
{
  int64_poly_realloc (f, i + 1);
  f->coeff[i] = coeff;
  if (i >= f->deg) {
    int64_poly_cleandeg (f, i);
  }
}

void int64_poly_get_coeff(int64_t * res, int i, int64_poly_srcptr f)
{
  ASSERT_ALWAYS( f->deg == -1 ||  f->deg>=i );
  if (i > f->deg) {
    * res = 0;
  } else {
    * res = f->coeff[i];
  }
}

void int64_poly_set_xi(int64_poly_ptr f, int i)
{
  int64_poly_realloc (f, i + 1);
  for(int j = 0 ; j < i + 1 ; j++) {
    f->coeff[j] = (j == i);
  }
  f->deg = i;
}

void int64_poly_set_bxi(int64_poly_ptr f, int i, int64_t b)
{
  int j = 0;
  int64_poly_realloc (f, i + 1);
  for( ; j < i ; j++) {
    f->coeff[j] = 0;
  }
  f->coeff[j + 1] = b;
  f->deg = i;
}

void int64_poly_printf(int64_poly_srcptr f)
{
  int i;
  if (f->deg == -1) {
    printf("0\n");
    return;
  }
  else if (f->deg == 0) {
    printf("%ld\n", f->coeff[0]);
    return;
  }
  printf("%ld", f->coeff[0]);
  for (i = 1; i <= f->deg; i++) {
    if (f->coeff[i] > 0) {
      printf("+%ld*x^%d", f->coeff[i], i);
    } else {
      printf("%ld*x^%d", f->coeff[i], i);
    }
  }
  printf("\n");
}

void int64_poly_fprintf(FILE * file, int64_poly_srcptr f)
{
  int i;
  if (f->deg == -1) {
    fprintf(file, "0\n");
    return;
  }
  else if (f->deg == 0) {
    printf("%ld\n", f->coeff[0]);
    return;
  }
  fprintf(file, "%ld", f->coeff[0]);
  for (i = 1; i <= f->deg; i++) {
    if (f->coeff[i] > 0) {
      fprintf(file, "+%ld*x^%d", f->coeff[i], i);
    } else {
      fprintf(file, "%ld*x^%d", f->coeff[i], i);
    }
  }
  fprintf(file, "\n");
}

void int64_poly_add_one_sieving_region(int64_poly_ptr f,
                                       sieving_interval_srcptr H)
{
  unsigned int k = 0;
  while(f->coeff[k] == H->h[k]) {
    f->coeff[k] = -H->h[k];
    k++;
  }
  if (k > H->t) {
    int64_poly_set_zero(f);
  } else {
    f->coeff[k] = f->coeff[k] + 1;
  }
}

void int64_poly_max(int64_t * max, int64_poly_srcptr f)
{
  int i = 1;
  * max = f->coeff[0];
  for ( ; i < f->deg + 1; i++) {
    * max = MAX(* max, f->coeff[i]);
  }
}

void int64_poly_infinite_norm(uint64_t * in, int64_poly_srcptr f)
{
  int i = 1;
  uint64_t tmp;
  int64_abs(in, f->coeff[0]);
  for ( ; i < f->deg + 1; i++) {
    int64_abs(&tmp, f->coeff[i]);
    * in = MAX(* in, tmp);
  }
}

/* int int64_poly_is_zero(int64_poly_t f) */
/* { */
/*   int i = 0; */
/*   for ( ; i <= f.deg; i ++) { */
/*     if (f.coeff[i] != 0) { */
/*       return 0; */
/*     } */
/*   } */
/*   return 1; */
/* } */

/* void int64_poly_delete_unsignificant(int64_poly_t * f) */
/* { */
/*   int k = f->deg; */
/*   int i = 0; */
/*   int64_poly_t tmp; */
/*   while (f->coeff[k] == 0) { */
/*     k--; */
/*   } */
/*   int64_poly_init(&tmp, k); */
/*   for ( ; i < k + 1; i++) { */
/*     int64_poly_setcoeff(i, f->coeff[i], &tmp); */
/*   } */
/*   int64_poly_copy(tmp, f); */
/* } */

/* int64_t int64_poly_resultant(int64_poly_t A, int64_poly_t B) */
/* { */
/*   int64_t a; */
/*   int64_t b; */
/*   int64_t g; */
/*   int64_t h; */
/*   int64_t s; */
/*   int64_t t; */
/*   int d; */
/*   int64_poly_t r; */
/*   int64_t tmp; */

/*   if (int64_poly_is_zero(A) || int64_poly_is_zero(B)) { */
/*     return (int64_t) 0; */
/*   } */

/*   /\* printf("A: "); *\/ */
/*   /\* int64_poly_printf(A); *\/ */
/*   /\* printf("B: "); *\/ */
/*   /\* int64_poly_printf(B); *\/ */

/*   a = int64_poly_content(A); */
/*   b = int64_poly_content(B); */

/*   int64_poly_divide_constant(a, &A); */
/*   int64_poly_divide_constant(b, &B); */

/*   g = 1; */
/*   h = 1; */
/*   s = 1; */
/*   t = pow_int64_t(a, (int64_t) B.deg) * pow_int64_t(b, (int64_t) A.deg); */

/*   /\* printf("g: %ld, h: %ld, s: %ld, t: %ld, a: %ld, b: %ld\n", g, h, s, t, a, b); *\/ */
/*   /\* printf("A: "); *\/ */
/*   /\* int64_poly_printf(A); *\/ */
/*   /\* printf("B: "); *\/ */
/*   /\* int64_poly_printf(B); *\/ */

/*   if (A.deg < B.deg) { */
/*     int64_poly_swap(&A, &B); */

/*     /\* printf("A: "); *\/ */
/*     /\* int64_poly_printf(A); *\/ */
/*     /\* printf("B: "); *\/ */
/*     /\* int64_poly_printf(B); *\/ */

/*     if ((A.deg % 2) == 1 && (B.deg % 2) == 1) { */
/*       /\* printf("s\n"); *\/ */
/*       s = -1; */
/*     } */
/*   } */

/*   while (B.deg > 0) { */
/*     d = A.deg - B.deg; */

/*     /\* printf("A: "); *\/ */
/*     /\* int64_poly_printf(A); *\/ */
/*     /\* printf("B: "); *\/ */
/*     /\* int64_poly_printf(B); *\/ */
/*     /\* printf("d: %d\n", d); *\/ */

/*     if ((A.deg % 2) == 1 && (B.deg % 2) == 1) { */
/*       /\* printf("s\n"); *\/ */
/*       s = -s; */
/*     } */
/*     int64_poly_pseudo_remainder(A, B, &r); */

/*     /\* printf("R: "); *\/ */
/*     /\* int64_poly_printf(r); *\/ */

/*     int64_poly_copy(B, &A); */

/*     /\* printf("A: "); *\/ */
/*     /\* int64_poly_printf(A); *\/ */

/*     tmp = g * pow_int64_t(h, d); */

/*     /\* printf("tmp: %ld\n", tmp); *\/ */

/*     int64_poly_divide_constant(tmp, &r); */
/*     int64_poly_copy(r, &B); */

/*     /\* printf("B: "); *\/ */
/*     /\* int64_poly_printf(B); *\/ */

/*     g = int64_poly_leading_coefficient(A); */
/*     h = pow_int64_t(g, d) / pow_int64_t(h, d-1); */

/*     /\* printf("h: %ld\n", h); *\/ */

/*     int64_poly_clear(&r); */
/*   } */

/*   /\* printf("B: "); *\/ */
/*   /\* int64_poly_printf(B); *\/ */
/*   h = pow_int64_t(int64_poly_leading_coefficient(B), A.deg) / pow_int64_t(h, */
/*   A.deg - 1); */
/*   /\* printf("h: %ld\n", h); *\/ */
/*   return s * t * h; */
/* } */

/* int64_t int64_poly_leading_coefficient(int64_poly_t f) */
/* { */
/*   int64_poly_delete_unsignificant(&f); */
/*   //TODO: This patch is maybe strange. */
/*   if (f.deg == -1) { */
/*       return 0; */
/*   } */
/*   return f.coeff[f.deg]; */
/* } */

/* void int64_poly_multiply_constant(int64_t a, int64_poly_t * f) */
/* { */
/*   int i = 0; */
/*   for ( ; i < f->deg + 1; i++) { */
/*     f->coeff[i] = a * f->coeff[i]; */
/*   } */
/*   int64_poly_delete_unsignificant(f); */
/* } */

/* void int64_poly_divide_constant(int64_t a, int64_poly_t * f) */
/* { */
/*   int i = 0; */
/*   for ( ; i < f->deg + 1; i++) { */
/*     f->coeff[i] = f->coeff[i] / a; */
/*   } */
/*   int64_poly_delete_unsignificant(f); */
/* } */

/* void int64_poly_swap(int64_poly_t * f, int64_poly_t * g) */
/* { */
/*   int64_poly_t tmp; */
/*   int64_poly_copy(* f, &tmp); */
/*   int64_poly_copy(* g, f); */
/*   int64_poly_copy(tmp, g); */
/*   int64_poly_clear(&tmp); */
/* } */

/* void int64_poly_add_poly(int64_poly_t a, int64_poly_t b, int64_poly_t * c) */
/* { */
/*   int d = MAX(a.deg, b.deg); */
/*   int i = 0; */
/*   if (d != a.deg) { */
/*     int64_poly_swap(&a, &b); */
/*   } */
/*   int64_poly_copy(a, c); */
/*   for ( ; i < b.deg + 1; i++) { */
/*     int64_poly_setcoeff(i, c->coeff[i] + b.coeff[i], c); */
/*   } */
/*   int64_poly_delete_unsignificant(c); */
/* } */

/* void int64_poly_mult_xi(int64_poly_t * a, int i) */
/* { */
/*   int64_poly_t tmp; */
/*   int k = 0; */
/*   int64_poly_copy(* a, &tmp); */
/*   int64_poly_init(a, tmp.deg + i); */
/*   for ( ; k < tmp.deg + 1 ; k ++) { */
/*     int64_poly_setcoeff(k + i, int64_poly_getcoeff(k, tmp), a); */
/*   } */
/*   int64_poly_clear(&tmp); */
/*   int64_poly_delete_unsignificant(a); */
/* } */

/* void int64_poly_mult_bxi(int64_poly_t * a, int64_t b, int i) */
/* { */
/*   int64_poly_mult_xi(a, i); */
/*   int64_poly_multiply_constant(b, a); */
/*   int64_poly_delete_unsignificant(a); */
/* } */

/* void int64_poly_mult_poly(int64_poly_t a, int64_poly_t b, int64_poly_t * */
/* 				c) */
/* { */
/*   int i = 0; */
/*   int64_poly_t d; */
/*   int64_poly_delete_unsignificant(&a); */
/*   int64_poly_delete_unsignificant(&b); */
/*   int64_poly_init(c, 1); */
/*   for ( ; i < a.deg + 1; i++) { */
/*     int64_poly_copy(b, &d); */
/*     int64_poly_mult_bxi(&d, a.coeff[i], i); */
/*     int64_poly_add_poly(* c, d, c); */
/*   } */
/*   int64_poly_clear(&d); */
/* } */

/* void int64_poly_division_with_remainder(int64_poly_t a, int64_poly_t b, */
/* int64_poly_t * quo, int64_poly_t * rem) */
/* { */
/*   int i = a.deg - b.deg; */
/*   int64_t value; */
/*   int64_poly_t tmp; */
/*   int64_poly_init(quo, i); */
/*   int64_poly_copy(a, rem); */
/*   for ( ; i >=0; i--) { */
/*     if (rem->deg == b.deg + i) { */
/*       value = int64_poly_leading_coefficient(* rem); */
/*       int64_poly_setcoeff(i, value, quo); */
/*       int64_poly_copy(b, &tmp); */
/*       int64_poly_mult_bxi(&tmp, value, i); */
/*       int64_poly_multiply_constant((int64_t) -1, &tmp); */
/*       int64_poly_add_poly(* rem, tmp, rem); */
/*     } else { */
/*       int64_poly_setcoeff(i, 0, quo); */
/*     } */
/*   } */
/* } */

/* int64_t int64_poly_content(int64_poly_t f) */
/* { */
/*   int64_poly_delete_unsignificant(&f); */
/*   int64_t tmp = f.coeff[0]; */
/*   int k = 1; */
/*   for ( ; k < f.deg + 1; k++) { */
/*     tmp = gcd_int64_t(tmp, f.coeff[k]); */
/*   } */
/*   return tmp; */
/* } */

/* void int64_poly_primitive_part(int64_poly_t * f) */
/* { */
/*   int64_t content = int64_poly_content(* f); */
/*   int64_poly_divide_constant(content, f); */
/* } */

/* /\* */
/*   This algorithm is also described in Cohen, page 111, pseudo-division. */
/*   See also vZZG, page 191, primitive euclidean algorithm. */
/* *\/ */
/* void int64_poly_pseudo_division(int64_poly_t a, int64_poly_t b, int64_poly_t * */
/* 				quo, int64_poly_t * rem) */
/* { */
/*   int64_t d = int64_poly_leading_coefficient(b); */
/*   int m = a.deg; */
/*   int n = b.deg; */
/*   int e = m - n + 1; */
/*   int64_t q; */
/*   int64_poly_t s; */
/*   int64_poly_copy(a, rem); */
/*   int64_poly_init(quo, a.deg - b.deg); */
/*   while(rem->deg >= b.deg) { */
/*     int64_poly_init(&s, rem->deg - n); */
/*     int64_poly_setcoeff(rem->deg - n, int64_poly_leading_coefficient(* rem), */
/* 			&s); */
/*     int64_poly_multiply_constant(d, quo); */
/*     int64_poly_add_poly(* quo, s, quo); */
/*     int64_poly_multiply_constant(d, rem); */
/*     int64_poly_multiply_constant((int64_t) -1, &s); */
/*     int64_poly_mult_poly(b, s, &s); */
/*     int64_poly_add_poly(* rem, s, rem); */
/*     e--; */
/*     int64_poly_clear(&s); */
/*   } */
/*   q = pow_int64_t(d, (int64_t) e); */
/*   int64_poly_multiply_constant(q, quo); */
/*   int64_poly_multiply_constant(q, rem); */
/* } */

/* void int64_poly_pseudo_remainder(int64_poly_t a, int64_poly_t b, int64_poly_t * */
/* 				rem) */
/* { */
/*   int64_t d = int64_poly_leading_coefficient(b); */
/*   int m = a.deg; */
/*   int n = b.deg; */
/*   int e = m - n + 1; */
/*   int64_t q; */
/*   int64_poly_t s; */
/*   int64_poly_copy(a, rem); */
/*   while(rem->deg >= b.deg) { */
/*     int64_poly_init(&s, rem->deg - n); */
/*     int64_poly_setcoeff(rem->deg - n, int64_poly_leading_coefficient(* rem), */
/* 			&s); */
/*     int64_poly_multiply_constant(d, rem); */
/*     int64_poly_multiply_constant((int64_t) -1, &s); */
/*     int64_poly_mult_poly(b, s, &s); */
/*     int64_poly_add_poly(* rem, s, rem); */
/*     e--; */
/*     int64_poly_clear(&s); */
/*   } */
/*   q = pow_int64_t(d, (int64_t) e); */
/*   int64_poly_multiply_constant(q, rem); */
/* } */

/* //TODO: link with mpz_poly.c */
/* /\* void int64_poly_to_mpz_poly(int64_poly_t a, mpz_poly_t * b) *\/ */
/* /\* { *\/ */
/* /\*   int i = 0; *\/ */
/* /\*   int64_poly_delete_unsignificant(&a); *\/ */
/* /\*   mpz_poly_init(b, a.deg); *\/ */
/* /\*   for ( ; i < a.deg + 1; i++) { *\/ */
/* /\*     mpz_poly_setcoeff(b, i, (mpz_t) a.coeff[i]); *\/ */
/* /\*   } *\/ */
/* /\* } *\/ */

void int64_poly_to_mpz_poly(mpz_poly_ptr a, int64_poly_srcptr b)
{
  mpz_t coeff;
  mpz_init(coeff);
  if (b->deg != -1) {
    for (int i = 0; i <= b->deg; i++) {
      mpz_poly_setcoeff_int64(a, i, b->coeff[i]);
    }
  }
  mpz_clear(coeff);
}

