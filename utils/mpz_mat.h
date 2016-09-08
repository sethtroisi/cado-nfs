#ifndef MPZ_MAT_H_
#define MPZ_MAT_H_

#include <stdio.h>
#include <stdint.h>
#include <gmp.h>
#include "macros.h"
#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* typedefs */

// mpz_t : multiple precision Z-elements ; same for mpq_t
// a matric on mpz_t type
struct mpz_mat_s {
    mpz_t *x;
    unsigned int m,n;
};
typedef struct mpz_mat_s mpz_mat[1];
typedef struct mpz_mat_s * mpz_mat_ptr;
typedef const struct mpz_mat_s * mpz_mat_srcptr;

// same on mpq_t type
struct mpq_mat_s {
    mpq_t *x;
    unsigned int m,n;
};
typedef struct mpq_mat_s mpq_mat[1];
typedef struct mpq_mat_s * mpq_mat_ptr;
typedef const struct mpq_mat_s * mpq_mat_srcptr;

/* entry access*/
mpz_ptr mpz_mat_entry(mpz_mat_ptr M, unsigned int i, unsigned int j);
mpz_srcptr mpz_mat_entry_const(mpz_mat_srcptr M, unsigned int i, unsigned int j);
mpq_ptr mpq_mat_entry(mpq_mat_ptr M, unsigned int i, unsigned int j);
mpq_srcptr mpq_mat_entry_const(mpq_mat_srcptr M, unsigned int i, unsigned int j);

/* init/clear/realloc*/
void mpz_mat_init(mpz_mat_ptr M, unsigned int m, unsigned int n);
void mpz_mat_clear(mpz_mat_ptr M);
void mpz_mat_realloc(mpz_mat_ptr M, unsigned int m, unsigned int n);
void mpq_mat_init(mpq_mat_ptr M, unsigned int m, unsigned int n);
void mpq_mat_clear(mpq_mat_ptr M);
void mpq_mat_realloc(mpq_mat_ptr M, unsigned int m, unsigned int n);

/* operations on submatrices, and swaps*/
void mpz_mat_submat_swap(
        mpz_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpz_mat_ptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn);
void mpz_mat_swap(mpz_mat_ptr A, mpz_mat_ptr B);
void mpq_mat_submat_swap(
        mpq_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpq_mat_ptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn);
void mpq_mat_swap(mpq_mat_ptr A, mpq_mat_ptr B);

/* set/set_ui*/
void mpz_mat_set(mpz_mat_ptr dst, mpz_mat_srcptr src);
void mpz_mat_set_ui(mpz_mat_ptr M, unsigned int a);
void mpq_mat_set_ui(mpq_mat_ptr M, unsigned int a);
void mpq_mat_set(mpq_mat_ptr dst, mpq_mat_srcptr src);
void mpz_mat_urandomm(mpz_mat_ptr M, gmp_randstate_t state, mpz_srcptr p);
void mpq_mat_urandomm(mpq_mat_ptr M, gmp_randstate_t state, mpz_srcptr p);

/* Joining matrices */
// I didn't write all possible functions there, only those I needed, but it could be completed
void mpz_mat_vertical_join(mpz_mat_ptr N, mpz_mat_srcptr M1, mpz_mat_srcptr M2);
void mpq_mat_vertical_join(mpq_mat_ptr N, mpq_mat_srcptr M1, mpq_mat_srcptr M2);

/* determinant and trace */
// We assume that M is square
void mpz_mat_trace(mpz_ptr t, mpz_mat_srcptr M);
// We assume that M is triangular (and square)
void mpz_mat_determinant_triangular(mpz_ptr d, mpz_mat_srcptr M);
// We assume that M is square
void mpq_mat_trace(mpq_ptr t, mpq_mat_srcptr M);
// We assume that M is triangular (and square)
void mpq_mat_determinant_triangular(mpq_ptr d, mpq_mat_srcptr M);
void mpz_mat_transpose(mpz_mat_ptr D, mpz_mat_srcptr M);

/* miscellaneous */
int mpq_mat_numden(mpz_mat_ptr num, mpz_ptr den, mpq_mat_srcptr M);
void mpz_mat_to_mpq_mat(mpq_mat_ptr N, mpz_mat_srcptr M);
void mpz_mat_mod_ui(mpz_mat_ptr dst, mpz_mat_srcptr src, unsigned int p);

/*  row-level operations */
// Return 1 if the k-th line of M is null, 0 otherwise
int mpz_mat_isnull_row(mpz_mat_srcptr M, unsigned int k);
void mpz_mat_swaprows(mpz_mat_ptr M, unsigned int i0, unsigned int i1);
void mpq_mat_swaprows(mpq_mat_ptr M, unsigned int i0, unsigned int i1);
/* apply a circular shift on rows [i0...i0+k[ */
void mpz_mat_rotatedownrows(mpz_mat_ptr M, unsigned int i0, unsigned int k);
void mpq_mat_rotatedownrows(mpq_mat_ptr M, unsigned int i0, unsigned int k);
/* put row perm[k] in row k */
void mpz_mat_permuterows(mpz_mat_ptr M, unsigned int * perm);
void mpq_mat_permuterows(mpq_mat_ptr M, unsigned int * perm);
/* add lambda times row i1 to row i0 */
void mpz_mat_addmulrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda);
void mpz_mat_addmulrow_mod(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda, mpz_srcptr p);
/* add lambda times row i1 to row i0 */
void mpq_mat_addmulrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1, mpq_srcptr lambda);
/* subtract lambda times row i1 to row i0 */
void mpz_mat_submulrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda);
void mpz_mat_submulrow_mod(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda, mpz_srcptr p);
/* subtract lambda times row i1 to row i0 */
void mpq_mat_submulrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1, mpq_srcptr lambda);
/* add row i1 to row i0 */
void mpz_mat_addrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1);
/* add row i1 to row i0 */
void mpq_mat_addrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1);
/* subtract row i1 to row i0 */
void mpz_mat_subrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1);
/* subtract row i1 to row i0 */
void mpq_mat_subrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1);
/* multiply row i0 by lambda */
void mpz_mat_mulrow(mpz_mat_ptr M, unsigned int i0, mpz_srcptr lambda);
void mpz_mat_mulrow_mod(mpz_mat_ptr M, unsigned int i0, mpz_srcptr lambda, mpz_srcptr p);
/* multiply row i0 by lambda */
void mpq_mat_mulrow(mpq_mat_ptr M, unsigned int i0, mpq_srcptr lambda);

/* this computes an additive combination of n rows into row [didx] of the
 * initial matrix. We require that this destination row be cleared
 * initially.
 */
void mpz_mat_combinerows(mpz_mat_ptr M, unsigned int didx, unsigned int sidx,
        mpz_srcptr * lambda, unsigned int n);
/*  */
/* I/O*/
void mpz_mat_fprint(FILE * stream, mpz_mat_srcptr M);
void mpq_mat_fprint(FILE * stream, mpq_mat_srcptr M);

/* multiplication */
void mpz_mat_multiply(mpz_mat_ptr D, mpz_mat_srcptr A, mpz_mat_srcptr B);
void mpz_mat_multiply_mod_ui(mpz_mat_ptr D, mpz_mat_srcptr A, mpz_mat_srcptr B, unsigned int p);
void mpz_mat_multiply_by_mpz(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_ptr k);
void mpz_mat_multiply_by_si(mpz_mat_ptr C, mpz_mat_srcptr A, int k);
void mpz_mat_multiply_by_ui(mpz_mat_ptr C, mpz_mat_srcptr A, unsigned int k);
// Returns A^n for n >= 2 ; assume A is a square matrix ; it's possible to use the same variable for A and B, but you lose the contents of A
void mpz_mat_power_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned int n);
// Returns A^n mod p for n >= 2 ; assume A is a square matrix ; it's possible to use the same variable for A and B, but you lose the contents of A
void mpz_mat_power_ui_mod_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned int n, unsigned int p);
void mpq_mat_multiply(mpq_mat_ptr D, mpq_mat_srcptr A, mpq_mat_srcptr B);
void mpq_mat_multiply_by_mpq(mpq_mat_ptr C, mpq_mat_srcptr A, mpq_ptr k);
void mpq_mat_multiply_by_si(mpq_mat_ptr C, mpq_mat_srcptr A, int k);
void mpq_mat_multiply_by_ui(mpq_mat_ptr C, mpq_mat_srcptr A, unsigned int k);
void mpz_mat_add(mpz_mat_ptr D, mpz_mat_srcptr A, mpz_mat_srcptr B);
/* XXX what does that do ??? */
void mpz_poly_eval_mpz_mat(mpz_mat_ptr D, mpz_mat_srcptr M, mpz_poly_srcptr f);
void mpz_poly_eval_mpz_mat_mod_ui(mpz_mat_ptr D, mpz_mat_srcptr M, mpz_poly_srcptr f, unsigned int p);

/*  gaussian reduction over the rationals
 * this is a backend for row gaussian reduction. T receives the
 * transformation matrix. M is modified in place.
 * this includes reordering of the rows.
 *
 * this works reasonably well over Fp, but for rationals we suffer from
 * coefficient blowup.
 *
 * T may be NULL in case we don't care.
 */
void mpq_gauss_backend(mpq_mat_ptr M, mpq_mat_ptr T);
void mpz_gauss_backend_mod_ui(mpz_mat_ptr M, mpz_mat_ptr T, unsigned long p);


/* Gaussian reduction over Z/pZ
 *
 * We sort of return something even in the case where p is not prime.
 *
 * M is the reduced matrix (modification in place of the input), T the
 * transformation matrix (so that T * input-M = output-M)
 */
void mpz_gauss_backend_mod(mpz_mat_ptr M, mpz_mat_ptr T, mpz_srcptr p);

/* this computes the row hnf on a column C.
 * The result is always of the form C'=(gcd(C),0,0,...,0). The transform
 * matrix T such that C'=T*C is most interesting.  a is modified in
 * place.
 */
void mpz_gcd_many(mpz_mat_ptr dT, mpz_ptr * a, unsigned int n);

/* return +1 or -1, which is the determinant of the transformation matrix T */
int mpz_hnf_backend(mpz_mat_ptr M, mpz_mat_ptr T);

/* kernel*/
// This is supposed to compute the Kernel of M mod p and to store it in the matrix K. If r is the rank of M, and M is a square matrix n*n, K is a n*(n-r) matrix
void mpz_mat_kernel(mpz_mat_ptr K, mpz_mat_srcptr M, unsigned int p);

void mpq_mat_invert(mpq_mat_ptr dst, mpq_mat_srcptr src);

int mpq_mat_eq(mpq_mat_srcptr A, mpq_mat_srcptr B);
int mpq_mat_cmp(mpq_mat_srcptr M, mpq_mat_srcptr N);
void mpq_mat_fprint_as_mpz(FILE* f, mpq_mat_srcptr M);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
struct cxx_mpz_mat {
    mpz_mat x;
    cxx_mpz_mat() { mpz_mat_init(x, 0, 0); }
    ~cxx_mpz_mat() { mpz_mat_clear(x); }
    cxx_mpz_mat(cxx_mpz_mat const & o) {
        mpz_mat_init(x, 0, 0);
        mpz_mat_set(x, o.x);
    }
    cxx_mpz_mat & operator=(cxx_mpz_mat const & o) {
        mpz_mat_set(x, o.x);
        return *this;
    }
#if __cplusplus >= 201103L
    cxx_mpz_mat(cxx_mpz_mat && o) {
        mpz_mat_init(x, 0, 0);
        mpz_mat_swap(x, o.x);
    }
    cxx_mpz_mat& operator=(cxx_mpz_mat && o) {
        mpz_mat_swap(x, o.x);
        return *this;
    }
#endif
    operator mpz_mat_ptr() { return x; }
    operator mpz_mat_srcptr() const { return x; }
    mpz_mat_ptr operator->() { return x; }
    mpz_mat_srcptr operator->() const { return x; }
};
struct cxx_mpq_mat {
    mpq_mat x;
    cxx_mpq_mat() { mpq_mat_init(x, 0, 0); }
    ~cxx_mpq_mat() { mpq_mat_clear(x); }
    cxx_mpq_mat(cxx_mpq_mat const & o) {
        mpq_mat_init(x, 0, 0);
        mpq_mat_set(x, o.x);
    }
    cxx_mpq_mat & operator=(cxx_mpq_mat const & o) {
        mpq_mat_set(x, o.x);
        return *this;
    }
#if __cplusplus >= 201103L
    cxx_mpq_mat(cxx_mpq_mat && o) {
        mpq_mat_init(x, 0, 0);
        mpq_mat_swap(x, o.x);
    }
    cxx_mpq_mat& operator=(cxx_mpq_mat && o) {
        mpq_mat_swap(x, o.x);
        return *this;
    }
#endif
    /* As a bonus, allow conversion and assignment from mpz matrices,
     * too. Move ctor and move assignment do not make sense in that case,
     * of course. */
    cxx_mpq_mat(cxx_mpz_mat const & o) {
        mpq_mat_init(x, 0, 0);
        mpz_mat_to_mpq_mat(x, o.x);
    }
    cxx_mpq_mat & operator=(cxx_mpz_mat const & o) {
        mpz_mat_to_mpq_mat(x, o.x);
        return *this;
    }
    operator mpq_mat_ptr() { return x; }
    operator mpq_mat_srcptr() const { return x; }
    mpq_mat_ptr operator->() { return x; }
    mpq_mat_srcptr operator->() const { return x; }
};
#endif

#endif	/* MPZ_MAT_H_ */
