#ifndef IDEAL_H
#define IDEAL_H

#include <stdint.h>
#include "cado.h"
#include "utils.h"

/*
  Represent an ideal.
*/
typedef struct {
  uint64_t r;
  mpz_poly_t h;
} s_ideal_t;

typedef s_ideal_t ideal_t[1];
typedef s_ideal_t * ideal_ptr;
typedef const s_ideal_t * ideal_srcptr;

typedef struct {
  ideal_t ideal;
  mpz_t * Tr; // Keep mpz_t because of the explosion of g0^i.
  unsigned char log;
} s_ideal_1_t;

//1 for ideal (r, h) with deg(h) = 1.
typedef s_ideal_1_t ideal_1_t[1];
typedef s_ideal_1_t * ideal_1_ptr;
typedef const s_ideal_1_t * ideal_1_srcptr;

//u for upper deg (h) > 1.
typedef struct {
  ideal_t ideal;
  mpz_t ** Tr; // Keep mpz_t because of big coefficients.
  unsigned char log;
} s_ideal_u_t;

typedef s_ideal_u_t ideal_u_t[1];
typedef s_ideal_u_t * ideal_u_ptr;
typedef const s_ideal_u_t * ideal_u_srcptr;

//pr for projective root.
typedef struct {
  ideal_t ideal;
  mpz_t * Tr; // Keep mpz_t because of big coefficients.
  unsigned char log;
} s_ideal_pr_t;

typedef s_ideal_pr_t ideal_pr_t[1];
typedef s_ideal_pr_t * ideal_pr_ptr;
typedef const s_ideal_pr_t * ideal_pr_srcptr;

/* Ideal: (r, h) */

/*
  Initialize an ideal.

  ideal: an ideal.
*/
void ideal_init(ideal_ptr ideal);

/*
  Set an ideal.

  ideal: ideal_ptr, the ideal we want to set.
  r: r is equal to the r of the ideal.
  h: h is equal to the h of the ideal.
*/
void ideal_set_part(ideal_ptr ideal, uint64_t r, mpz_poly_srcptr h);

/*
  Delete an ideal.

  ideal: the ideal we want to delete.
*/
void ideal_clear(ideal_ptr ideal);

/*
  Write an ideal in a file.

  file: the file.
  ideal: the ideal we want to write.
*/
void ideal_fprintf(FILE * file, ideal_srcptr ideal);

/* Ideal with h of degree 1 and Tr */

/*
  Initialize an ideal_1.

  ideal: an ideal.
*/
void ideal_1_init(ideal_1_ptr ideal);

/*
  Set an ideal_1, compute Tr and the logarithm.

  ideal: the ideal we want to set.
  r: r is equal to the r of the ideal.
  h: h is equal to the h of the ideal.
  t: t is the dimension of the lattice we explore.
*/
void ideal_1_set_part(ideal_1_ptr ideal, uint64_t r, mpz_poly_srcptr h,
                 unsigned int t);

/*
  Delete an ideal_1.

  ideal: the ideal we want to delete.
  t: t is the dimension of the lattice we explore.
*/
void ideal_1_clear(ideal_1_ptr ideal, unsigned int t);

/*
  Write an ideal_1.

  ideal: the ideal we want to delete.
  t: t is the dimension of the lattice we explore.
*/
void ideal_1_fprintf(FILE * file, ideal_1_srcptr ideal, unsigned int t);

/* Ideal with h of degree u > 1 and Tr */

/*
  Initialize an ideal_u.

  ideal: an ideal.
*/
void ideal_u_init(ideal_u_ptr ideal);

/*
  Set an ideal_u, compute Tr and the logarithm.

  ideal: the ideal we want to set.
  r: r is equal to the r of the ideal.
  h: h is equal to the h of the ideal.
  t: t is the dimension of the lattice we explore.
*/
void ideal_u_set_part(ideal_u_ptr ideal, uint64_t r, mpz_poly_srcptr h,
                      unsigned int t);

/*
  Delete an ideal_u.

  ideal: ideal_u_srcptr, the ideal we want to delete.
  t: t is the dimension of the lattice we explore.
*/
void ideal_u_clear(ideal_u_ptr ideal, unsigned int t);

/*
  Write an ideal_u.

  file: the file in which we write.
  ideal: ideal_u_srcptr, the ideal we want to write.
  t: t is the dimension of the lattice we explore.
*/
void ideal_u_fprintf(FILE * file, ideal_u_srcptr ideal, unsigned int t);

/* Ideal projective root and Tr */

void ideal_pr_init(ideal_pr_ptr ideal);

/*
  Set an ideal_pr, compute Tr and the logarithm. The polynomial is x.

  ideal: the ideal we want to set.
  r: r is equal to the r of the ideal.
  t: t is the dimension of the lattice we explore.
*/
void ideal_pr_set_part(ideal_pr_ptr ideal, uint64_t r, unsigned int t);

/*
  Delete an ideal_pr.

  ideal: ideal_pr_srcptr, the ideal we want to delete.
  t: t is the dimension of the lattice we explore.
*/
void ideal_pr_clear(ideal_pr_ptr ideal, unsigned int t);

/*
  Write an ideal_pr.

  file: the file in which we write.
  ideal: ideal_pr_srcptr, the ideal we want to write.
  t: t is the dimension of the lattice we explore.
*/
void ideal_pr_fprintf(FILE * file, ideal_pr_srcptr ideal, unsigned int t);

#endif /* IDEAL_H */
