#include <iostream>
#include <math.h>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <string>
#include <errno.h>
#include <gmp.h>
#include <stdint.h>
#include <sys/resource.h>	/* for getrusage */
#include <queue>
#include <list>
#include <time.h>
#include <algorithm>
#include "macros.h"
#include "mpz_poly.h"
#include "mpz_mat.h"
#include "utils/cxx_mpz.hpp"
#include "utils/gmp_aux.h"
#include "ant.hpp"

using namespace std;

/*{{{ conversion of rows and columns to polynomials*/
void mpz_mat_row_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int i)
{
    mpz_poly_realloc(f,M->n);
    unsigned int j;
    for (j = 0 ; j < M->n; j++){
        mpz_poly_setcoeff(f,j,mpz_mat_entry_const(M,i,j));
    }
    mpz_poly_cleandeg(f, M->n - 1);
}

void mpz_mat_row_to_poly_rev(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int i)
{
    mpz_poly_realloc(f,M->n);
    unsigned int j;
    for (j = 0 ; j < M->n; j++){
        mpz_poly_setcoeff(f,j,mpz_mat_entry_const(M,i,M->n-1-j));
    }
    mpz_poly_cleandeg(f, M->n - 1);
}

void mpz_mat_column_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int j)
{
    mpz_poly_realloc(f,M->m);
    unsigned int i;
    for (i = 0 ; i < M->m; i++){
        mpz_poly_setcoeff(f,i,mpz_mat_entry_const(M,i,j));
    }
    mpz_poly_cleandeg(f, M->m - 1);
}

void mpq_mat_row_to_poly(mpz_poly_ptr f, mpz_ptr lcm, mpq_mat_srcptr M, const unsigned int i)
{
    /* read element w[i] as a polynomial. beware, the element might
     * have rational coordinates, that's life ! */
    unsigned int n = M->n;
    ASSERT_ALWAYS(i < M->m);
    mpz_set_ui(lcm, 1);
    for (unsigned int j = 0; j < n; j++) {
        mpz_lcm(lcm, lcm, mpq_denref(mpq_mat_entry_const(M, i, j)));
    }
    mpz_poly_realloc(f, n);
    for (unsigned int j = 0; j < n; j++) {
        mpq_srcptr mij = mpq_mat_entry_const(M, i, j);
        mpz_divexact(f->coeff[j], lcm, mpq_denref(mij));
        mpz_mul(f->coeff[j], f->coeff[j], mpq_numref(mij));
    }
    mpz_poly_cleandeg(f, n-1);
}

void mpq_poly_to_mat_row(mpq_mat_ptr M, const unsigned int i, mpz_poly_srcptr f, mpz_srcptr denom)
{
    ASSERT_ALWAYS(f->deg < (int) M->n);
    for (unsigned int j = 0 ; j < M->n; j++){
        mpq_ptr mij = mpq_mat_entry(M,i,j);
        mpz_poly_getcoeff(mpq_numref(mij), j, f);
        mpz_set(mpq_denref(mij), denom);
        mpq_canonicalize(mij);
    }
}

void mpq_mat_column_to_poly(mpz_poly_ptr f, mpz_ptr lcm, mpq_mat_srcptr M, const unsigned int j)
{
    mpz_poly_realloc(f,M->m);
    mpz_set_si(lcm,1);
    for (unsigned int i = 0 ; i < M->m ; i++) {
        mpz_lcm(lcm,lcm,mpq_denref(mpq_mat_entry_const(M,i,j)));
    }
    for (unsigned int i = 0 ; i < M->m; i++){
        mpq_srcptr mij = mpq_mat_entry_const(M, i, j);
        mpz_divexact(f->coeff[i], lcm, mpq_denref(mij));
        mpz_mul(f->coeff[i], f->coeff[i], mpq_numref(mij));
    }
    mpz_poly_cleandeg(f, M->m - 1);
}
/*}}}*/

/*{{{ commonly used wrappers around HNF functions */
int mpz_mat_hnf_backend_rev(mpz_mat_ptr M, mpz_mat_ptr T) // {{{
{
    /* This is almost like hnf_backend, except that we do it in a
     * different order, which is more suitable for displaying number
     * field elements in a way which ends up being similar to magma's
     *
     * T receives the transformation matrix.
     * M is put into HNF form.
     */
    mpz_mat_reverse_rows(M, M);
    mpz_mat_reverse_columns(M, M);
    int s = mpz_mat_hnf_backend(M, T);
    mpz_mat_reverse_rows(M, M);
    mpz_mat_reverse_columns(M, M);
    if (T) mpz_mat_reverse_rows(T, T);
    if (T) mpz_mat_reverse_columns(T, T);
    if (M->m > M->n) {
        /* we need some swaps... */
        mpz_mat sM;
        mpz_mat_init(sM, M->m, M->n);
        mpz_mat_submat_swap(sM, 0, 0,    M, M->m-M->n, 0, M->n, M->n);
        mpz_mat_submat_swap(sM, M->n, 0, M, 0, 0,         M->m-M->n, M->n);
        mpz_mat_swap(sM, M);
        mpz_mat_clear(sM);
        if (T) {
            mpz_mat sT;
            mpz_mat_init(sT, T->m, T->n);
            mpz_mat_submat_swap(sT, 0, 0,    T, T->m-T->n, 0, T->n, T->n);
            mpz_mat_submat_swap(sT, T->n, 0, T, 0, 0,         T->m-T->n, T->n);
            mpz_mat_swap(sT, T);
            mpz_mat_clear(sT);
        }
        /* While the transformations above had no effect on s (because
         * they compensate), this one has.
         * we have n circular shifts on length m, plus a reversal on m-n.
         * a circular shift on length k is exactly k-1 inversions, so
         * that sums up to n*(m-1) inversions. Then we add
         * (m-n)*(m-n-1)/2 inversions. This is, in total,
         * (m(m-1)+n(n-1))/2 inversions.
         * m*(m-1) is congruent to 2 mod 4 when m is 2 or 3 mod 4
         */
        int ninvs = ((M->m&2)+(M->n&2))/2;
        if (ninvs) s=-s;
    }
    return s;
}//}}}

cxx_mpz_mat join_HNF(cxx_mpz_mat const& K, cxx_mpz const& p)//{{{
{
    cxx_mpz_mat J(K->n, K->n, p);
    cxx_mpz_mat T0;
    cxx_mpz_mat I;

    mpz_mat_vertical_join(I, J, K);
    mpz_mat_hnf_backend(I, T0);
    mpz_mat_submat_swap(I, 0, 0, J, 0, 0, K->n, K->n);
    return J;
}
//}}}

cxx_mpz_mat join_HNF_rev(cxx_mpz_mat const& K, cxx_mpz const& p)//{{{
{
    // Builds the block matrix containing p*identity in the top, and K in
    // the bottom Then computes its HNF and stores it in I
    cxx_mpz_mat J(K->n, K->n, p);
    cxx_mpz_mat T0;
    cxx_mpz_mat I;

    mpz_mat_vertical_join(I, J, K);
    mpz_mat_hnf_backend_rev(I, T0);
    mpz_mat_submat_swap(I, 0, 0, J, 0, 0, K->n, K->n);
    return J;
}//}}}
/*}}}*/

/* {{{ multiplication_table_of_order */
cxx_mpz_mat multiplication_table_of_order(cxx_mpq_mat const& O,
                                  cxx_mpz_poly const& g)
{
    /* Let O be an order, with basis written with respect to the
     * polynomial basis defined by g (of degree denoted below by n). O is
     * thus an n*n matrix.
     *
     * This function computes the integer matrix of size n*n^2 such that
     * the n coordinates at position (i,j*n) to (i,j*n+n-1) are the
     * coordinates of the i-th times the j-th generator of O, expressed
     * as combinations of the generators of O.
     */
    unsigned int n = g->deg;
    ASSERT_ALWAYS(O->m == n);
    ASSERT_ALWAYS(O->n == n);
    cxx_mpq_mat R;
    mpq_mat_inv(R, O);
    cxx_mpz_mat M(n, n * n);

    for(unsigned int i = 0; i < n; i++) {
        cxx_mpz_poly w;
        cxx_mpz dw;
        mpq_mat_row_to_poly(w, dw, O, i);
        cxx_mpq_mat T(n, n);
        for (unsigned int j = 0; j < n; j++) {
            cxx_mpz_poly c;
            cxx_mpz dc;
            mpq_mat_row_to_poly(c, dc, O, j);
            mpz_poly_mul_mod_f(c, c, w, g);
            mpz_mul(dc, dc, dw);
            mpq_poly_to_mat_row(T, j, c, dc);
        }
        mpq_mat_mul(T, T, R);
        cxx_mpz_mat Tz;
        int rc = mpq_mat_numden(Tz, NULL, T);
        ASSERT_ALWAYS(rc);
        for (unsigned int j = 0; j < n; j++) {
            mpz_mat_submat_swap(M, i, j*n, Tz, j, 0, 1, n);
        }
    }
    return M;
}/*}}}*/

/* {{{ multiplication_table_of_ideal*/
cxx_mpz_mat multiplication_table_of_ideal(cxx_mpz_mat const& M,
				  cxx_mpz_mat const& I)
{
    /* Let O be an order of a degree n number field. Let M be the n*n^2
     * multiplication matrix of O.
     *
     * Let I be an ideal of O, given as an n*n matrix with entries
     * expressed as coordinate vectors with respect to the basis O.
     *
     * This function computes the integer matrix of size n*n^2 such that
     * the n coordinates at position (i,j*n) to (i,j*n+n-1) are the
     * coordinates of the i-th generator of O times the j-th generator of
     * I, expressed as combinations of the generators of I.
     */
    unsigned int n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(I->m == n);
    ASSERT_ALWAYS(I->n == n);

    cxx_mpq_mat R;
    mpq_mat_inv(R, cxx_mpq_mat(I));
    cxx_mpz_mat MI(n, n * n);
    for(unsigned int i = 0 ; i < n ; i++) {
        cxx_mpz_mat Tz(n, n);
        for (unsigned int j = 0; j < n; j++) {
            mpz_mat_submat_set(Tz, j, 0, M, i, j*n, 1, n);
        }
        /* We have the matrix of multiplication by w (some generator of
         * O). Row i is w*wi.
         *
         * We need to compute I*T*I^-1 in order to have that represent
         * how multiplication by w affects the generators of I.
         */
        mpz_mat_mul(Tz, I, Tz);
        cxx_mpq_mat T(Tz);
        mpq_mat_mul(T, T, R);
        int rc = mpq_mat_numden(Tz, NULL, T);
        ASSERT_ALWAYS(rc);
        for (unsigned int j = 0; j < n; j++) {
            mpz_mat_submat_swap(MI, i, j*n, Tz, j, 0, 1, n);
        }
    }
    return MI;
}/*}}}*/

/*{{{ multiply_elements_in_order */
cxx_mpz_mat multiply_elements_in_order(cxx_mpz_mat const& M, cxx_mpz_mat const& E, cxx_mpz_mat const& F)
{
    /* Let O be an order in a degree n number field.
     *
     * Given the n*n^2 matrix M of the multiplication within that order,
     * given matrices E and F of size k*n with coordinates of k elements
     * of O, return the matrix of size k*n with i-th row giving
     * coordinates of e_i times f_i in O.
     */
    unsigned int n = M->m;
    unsigned int k = E->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(E->n == n);
    ASSERT_ALWAYS(F->n == n);
    ASSERT_ALWAYS(F->m == k);

    cxx_mpz_mat EM;
    mpz_mat_mul(EM, E, M);
    cxx_mpz_mat R(k, n);
    for(unsigned int ell = 0 ; ell < k ; ell++) {
        for (unsigned int j = 0; j < n; j++) {
            mpz_ptr Rlj = mpz_mat_entry(R, ell, j);
            for(unsigned int i = 0 ; i < n ; i++) {
                mpz_addmul(Rlj,
                        mpz_mat_entry_const(F, ell, i),
                        mpz_mat_entry_const(EM, ell, i*n+j));
            }
        }
    }
    return R;
}/*}}}*/

//{{{ frobenius_matrix
cxx_mpz_mat frobenius_matrix(cxx_mpz_mat const& M, cxx_mpz const& p)
{
    // frobenius_matrix ; utility function for computing the p-radical
    // Takes a matrix B containing the generators (w_0, ... w_{n-1}) of an
    // order, expressed as polynomials in the root of the polynomial g,
    // return the matrix U containing ((w_0)^p, ..., (w_{n-1})^p), and
    // expressed as a linear transformation within the order B.

    /* This version uses the multiplication table of the order (and does all
     * arithmetic mod p), and binary powering. */
    unsigned int n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    cxx_mpz_mat Mp = M;
    mpz_mat_mod_mpz(Mp, Mp, p);

    cxx_mpz_mat E(n, n, 1), F(E);

    int k = mpz_sizeinbase(p, 2) - 1;
    for( ; k-- ; ) {
        F = multiply_elements_in_order(M, F, F);
        if (mpz_tstbit(p, k))
            F = multiply_elements_in_order(M, E, F);
        mpz_mat_mod_mpz(F, F, p);
    }
    return F;
}//}}}

// {{{ cxx_mpz_mat p_radical_of_order
// Stores in I the p-radical of the order whose multiplication matrix is
// given by M. I is expressed with respect to the basis of the order.
cxx_mpz_mat p_radical_of_order(cxx_mpz_mat const& M, cxx_mpz const& p)
{
    unsigned int n = M->m;
    ASSERT_ALWAYS(M->n == n * n);

    cxx_mpz_mat K;

    // Now building the matrix U, containing all generators to the power
    // of p Generators are polynomials, stored in the matrix B
    cxx_mpz_mat T = frobenius_matrix(M, p);

    /* Which power of p ? */
    int k = 1;
    for (cxx_mpz pk = p; mpz_cmp_ui(pk, n) < 0; k++)
        mpz_mul(pk, pk, p);
    mpz_mat_pow_ui_mod_mpz(T, T, k, p);

    // Storing in K a basis of Ker((z -> (z^%p mod g mod p))^k)
    mpz_mat_kernel_mod_mpz(K, T, p);

    // Getting generators of the p radical from Ker(X) by computing HNF
    // of the vertical block matrix (p*Id, K);
    return join_HNF(K, p);
}// }}}

// {{{ cxx_mpq_mat p_maximal_order(cxx_mpz_poly const& f, cxx_mpz const& p)
// Builds one p-maximal-order starting from the matrix B, containing the
// generators of one order (in the basis of alpha^) In the number field
// of the polynomial f, of degree n Stores the generators of this
// p-maximal-order (in the basis of alpha^) in D
cxx_mpq_mat p_maximal_order(cxx_mpz_poly const& f, cxx_mpz const& p)
{
    unsigned int n = f->deg;

    cxx_mpz_poly g;
    mpz_poly_to_monic(g, f);
    cxx_mpq_mat B(n, n, 1);

    cxx_mpq_mat new_D = B;
    cxx_mpq_mat D = B;

    do {
        D = new_D;

        cxx_mpz_mat M = multiplication_table_of_order(D, g);

        // Getting the p-radical of the order generated by D in I.
        cxx_mpz_mat I = p_radical_of_order(M, p);

        // Building the (n,n^2) matrix containing the integers mod p
        // associated to all generators of O
        cxx_mpz_mat M2 = multiplication_table_of_ideal(M, I);
        mpz_mat_mod_mpz(M2, M2, p);
        //printf("M is :\n"); mpz_mat_fprint(stdout, M); printf("\n");
        
        // Computing Ker(M)
        cxx_mpz_mat K;	        // Ker(M)
        mpz_mat_kernel_mod_mpz(K, M2, p);
        // printf("Ker(M) is :\n"); mpz_mat_fprint(stdout, K_M); printf("\n");

        // Getting generators of p*O' by computing HNF of the vertical block matrix (p*Id, K_M);
        cxx_mpz_mat J = join_HNF(K, p);

        // Converting in the basis which is used to express elements of D
        mpq_mat_mul(new_D, cxx_mpq_mat(J), D);

        // Dividing by p
        mpq_mat_div_mpz(new_D, new_D, p);

    } while (D != new_D);

    /* express D back in the basis of f */
    cxx_mpz x;
    mpz_set_ui(x, 1);
    for(unsigned int j = 0 ; j < n ; j++) {
        for(unsigned int i = 0; i < n; i++) {
            mpq_ptr dij = mpq_mat_entry(D, i, j);
            mpz_mul(mpq_numref(dij), mpq_numref(dij), x);
            mpq_canonicalize(dij);
        }
        mpz_mul(x, x, f->coeff[f->deg]);
    }
    // Put D into HNF.
    cxx_mpz_mat Dz;
    cxx_mpz den;
    mpq_mat_numden(Dz, den, D);
    mpz_mat_hnf_backend_rev(Dz, NULL);
    mpq_mat_set_mpz_mat_denom(D, Dz, den);
    return D;
}
//}}}

//{{{ mpz_mat_minpoly_mod_ui
cxx_mpz_poly mpz_mat_minpoly_mod_mpz(cxx_mpz_mat M, cxx_mpz const& p)
{
    ASSERT_ALWAYS(M->m == M->n);
    unsigned long n = M->n;
    cxx_mpz_mat B(n+1, n*n);
    cxx_mpz_mat T(n,n,1);
    for(unsigned int k = 0 ; k < n + 1 ; k++) {
        for(unsigned int i = 0 ; i < n ; i++) {
            mpz_mat_submat_set(B, n - k, i * n, T, i, 0, 1, n);
        }
        mpz_mat_mul_mod_mpz(T, T, M, p);
    }
    cxx_mpz_mat K;
    mpz_mat_kernel_mod_mpz(K, B, p);
    mpz_mat_gauss_backend_mod_mpz(K, NULL, p);

    /* TODO: write down exactly what we need in terms of ordering from
     * mpz_mat_gauss_backend_mod_ui. I take it that we're happy if it
     * computes a RREF, so that the last row is the one with the largest
     * number of column coordinates killed, corresponding to highest
     * degree coefficients absent ?
     */
    cxx_mpz_poly f;
    mpz_mat_row_to_poly_rev(f,K,K->m-1);

    return f;
}
/*}}}*/

/*{{{ matrix_of_multmap */
cxx_mpz_mat matrix_of_multmap(
        cxx_mpz_mat const& M,
        cxx_mpz_mat const& J, 
        cxx_mpz_mat const& c,
        cxx_mpz const& p)
{
    /* Let O be an order, represented here simply by its multiplication
     * matrix. Let J be an m-dimensional subalgebra of O/pO. Let c be an
     * element of J (a linear combination of the rows of J).
     * This returns the m times m matrix which expresses multiplication by c
     * within J.
     */
    /* We write the coordinates in O/pO of all products of the for c*j_k,
     * j_k being one of the generators of J, and denote this matrix by
     * CJ.  Now we look for a matrix M such that M*J = CJ. To do this, we
     * need to do gaussian reduction on transpose(J): Find a matrix K
     * such that K*transpose(J) is verticaljoin(identity_matrix(m),
     * zero_matrix(n-m,m)). Then, taking K' as the first m rows of K, we
     * will have M = CJ * transpose(K').
     */
    unsigned int m = J->m;
    unsigned int n = M->m;
    ASSERT_ALWAYS(M->n == n*n);
    ASSERT_ALWAYS(J->n == n);
    ASSERT_ALWAYS(c->n == n);
    ASSERT_ALWAYS(c->m == 1);
    cxx_mpz_mat Mc;
    cxx_mpz_mat Jt;
    cxx_mpz_mat K;
    mpz_mat_transpose(Jt, J);
    cxx_mpz_mat T;
    mpz_mat_gauss_backend_mod_mpz(Jt, T, p);
    cxx_mpz_mat Kt(m, n);
    mpz_mat_submat_swap(Kt, 0, 0, T, 0, 0, m, n);
    mpz_mat_transpose(Kt, Kt);  /* Kt is now n times m */
    cxx_mpz_mat CJ(m, n);
    for(unsigned int k = 0 ; k < m ; k++) {
        cxx_mpz_mat jk(1,n);
        mpz_mat_submat_set(jk,0,0,J,k,0,1,n);
        cxx_mpz_mat w = multiply_elements_in_order(M, c, jk);
        mpz_mat_submat_swap(CJ,k,0,w,0,0,1,n);
    }
    mpz_mat_mul_mod_mpz(Mc, CJ, Kt, p);
    return Mc;
}
/*}}}*/

/*{{{ factorization_of_polynomial_mod_mpz */
vector<pair<cxx_mpz_poly, int>> factorization_of_polynomial_mod_mpz(cxx_mpz_poly const& f, cxx_mpz const& p, gmp_randstate_t state)
{
    mpz_poly_factor_list lf;
    mpz_poly_factor_list_init(lf);
    mpz_poly_factor(lf,f,p,state);
    vector<pair<cxx_mpz_poly, int>> res(lf->size);
    for(int i = 0 ; i < lf->size ; i++) {
        mpz_poly_swap(res[i].first, lf->factors[i]->f);
        res[i].second = lf->factors[i]->m;
    }
    mpz_poly_factor_list_clear(lf);
    return res;
}
/*}}}*/

/*{{{ template <typename T> void append_move(vector<T> &a, vector<T> &b) */
template <typename T> void append_move(vector<T> &a, vector<T> &b)
{
    a.reserve(a.size() + b.size());
    size_t na = a.size();
    size_t nb = b.size();
    if (&a == &b) {
        for(size_t i = 0 ; i < nb ; i++) {
            a.push_back(b[i]);
        }
    } else {
        a.insert(a.end(), nb, T());
        for(size_t i = 0 ; i < nb ; i++) {
            std::swap(a[na + i], b[i]);
        }
    }
}
/*}}}*/

// {{{ factorization_of_prime
vector<pair<cxx_mpz_mat, int>> factorization_of_prime_inner(
        cxx_mpq_mat const & B,
        cxx_mpz_mat const & M,
        cxx_mpz const& p,
        cxx_mpz_mat const& Ip,
        cxx_mpz_mat const& I,
        cxx_mpz_mat const& J,
        gmp_randstate_t state)
{
    unsigned int m = J->m;
    unsigned int n = J->n;
    ASSERT_ALWAYS(I->m == n);
    ASSERT_ALWAYS(I->n == n);
    /* J represents an m-dimensional subalgebra of O/pO */

    /* Pick a random element of J */
    cxx_mpz_mat c(1, m);
    for(unsigned int i = 0 ; i < m ; i++) {
        mpz_urandomm(mpz_mat_entry(c,0,i), state, p);
    }
    mpz_mat_mul_mod_mpz(c, c, J, p);

    cxx_mpz_mat Mc = matrix_of_multmap(M, J, c, p);

    /* Now we would like to find the minimal polynomial of Mc, which is
     * simply the matrix of multiplication by c. For this, we write an
     * (m+1) times m^2 matrix, and compute a left kernel.
     */
    cxx_mpz_poly Pc = mpz_mat_minpoly_mod_mpz(Mc, p);

    vector<pair<cxx_mpz_poly, int>> facP = factorization_of_polynomial_mod_mpz(Pc, p, state);

    vector<pair<cxx_mpz_mat, int>> ideals;

    vector<cxx_mpz_mat> characteristic_subspaces;

    for(unsigned int i = 0 ; i < facP.size() ; i++) {
        cxx_mpz_poly const& f(facP[i].first);
        /* We need the basis of the kernel of f(Mc) */
        cxx_mpz_mat E;
        mpz_poly_eval_mpz_mat_mod_mpz(E, Mc, f, p);
        mpz_mat_pow_ui_mod_mpz(E, E, facP[i].second, p);
        mpz_mat_kernel_mod_mpz(E, E, p);
        /* This line is just to be exactly in line with what magma says
         */
        mpz_mat_gauss_backend_mod_mpz(E, NULL, p);
        mpz_mat_mul_mod_mpz(E, E, J, p);
        characteristic_subspaces.push_back(E);
    }

    for(unsigned int i = 0 ; i < facP.size() ; i++) {
        unsigned int e = facP[i].second;
        cxx_mpz_poly const& f(facP[i].first);
        cxx_mpz_mat const& Ci(characteristic_subspaces[i]);
        /* We need to find an ideal which is smaller than Ip, and whose
         * generators are p*the generators of Ci, as well as the
         * unmodified generators of I and of the other characteristic
         * subspaces.
         */
        cxx_mpz_mat Ix(n + J->m - Ci->m, n);
        mpz_mat_submat_set(Ix,0,0,I,0,0,n,n);
        for(unsigned int r = n, k = 0; k < facP.size() ; k++) {
            if (k == i) continue;
            cxx_mpz_mat const& Ck(characteristic_subspaces[k]);
            unsigned int mk = Ck->m;
            mpz_mat_submat_set(Ix, r, 0, Ck, 0, 0, mk, n);
            r += mk;
        }
        cxx_mpz_mat Ihead(n, n);
        if (Ci->m == e * f->deg) {
            mpz_mat_vertical_join(Ix, Ix, Ip);
            mpz_mat_hnf_backend_rev(Ix, NULL);
            mpz_mat_submat_swap(Ihead,0,0,Ix,0,0,n,n);
            ideals.push_back(make_pair(Ihead, e));
        } else {
            mpz_mat_hnf_backend_rev(Ix, NULL);
            mpz_mat_submat_swap(Ihead,0,0,Ix,0,0,n,n);
            vector<pair<cxx_mpz_mat, int>> more_ideals;
            more_ideals = factorization_of_prime_inner(B,M,p,Ip,Ihead,Ci,state);
            append_move(ideals, more_ideals);
        }
    }
    sort(ideals.begin(), ideals.end(), ideal_comparator());
    return ideals;
}

vector<pair<cxx_mpz_mat, int>> factorization_of_prime(
        cxx_mpq_mat & B, cxx_mpz_poly const& g,
        cxx_mpz const& p,
        gmp_randstate_t state)
{
    int n = g->deg;
    cxx_mpz_mat M = multiplication_table_of_order(B, g);
    cxx_mpz_mat Ip = p_radical_of_order(M, p);
    return factorization_of_prime_inner(B, M, p, Ip,
            cxx_mpz_mat(n, n, p),
            cxx_mpz_mat(n, n, 1), state);
}
//}}}

// {{{ valuation_helper_for_ideal
//
// compute an uniformizing element for the prime ideal I of the order O,
// where I is above the rational prime p, and O is p-maximal and given by
// means of its multiplication matrix M.
//
// The uniformizing element is such that (a/p)*I is in O, yet a is not in
// p*O.
cxx_mpz_mat valuation_helper_for_ideal(cxx_mpz_mat const& M, cxx_mpz_mat const& I, cxx_mpz const& p)
{
    unsigned int n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(I->m == n);
    ASSERT_ALWAYS(I->n == n);

    /* We begin by something which is very similar to
     * multiplication_table_of_ideal. It's a bit like computing I*M,
     * except that we want it in a different order.
     */
    cxx_mpz_mat MI(n, n * n);
    for(unsigned int i = 0 ; i < n ; i++) {
        cxx_mpz_mat Tz(n, n);
        for (unsigned int j = 0; j < n; j++) {
            mpz_mat_submat_set(Tz, j, 0, M, i, j*n, 1, n);
        }
        mpz_mat_mul_mod_mpz(Tz, I, Tz, p);
        for (unsigned int j = 0; j < n; j++) {
            mpz_mat_submat_swap(MI, i, j*n, Tz, j, 0, 1, n);
        }
    }

    cxx_mpz_mat ker;
    mpz_mat_kernel_mod_mpz(ker, MI, p);
    mpz_mat_hnf_backend(ker, NULL);

    cxx_mpz_mat res(1, n);
    mpz_mat_submat_swap(res, 0, 0, ker, 0, 0, 1, n);
    return res;
}// }}}

// {{{ generate_ideal -- create an ideal from a set of generators
// generators (gens) are here given as elements of the field, not
// elements of the order. In case the ideal is fractional, its
// denominator is also returned. For an integral ideal, the denominator
// is always 1.
pair<cxx_mpz_mat, cxx_mpz> generate_ideal(cxx_mpq_mat const& O, cxx_mpz_mat const& M, cxx_mpq_mat const& gens)
{
    unsigned int n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(O->m == n);
    ASSERT_ALWAYS(O->n == n);
    ASSERT_ALWAYS(gens->n == n);
    cxx_mpq_mat R;
    mpq_mat_inv(R, O);
    cxx_mpq_mat I0q;
    mpq_mat_mul(I0q, gens, R);
    cxx_mpz_mat I0;
    cxx_mpz denom;
    mpq_mat_numden(I0, denom, I0q);
    /* Now create the full list of elements of O which are generated by
     * the products I_i * O_j */
    cxx_mpz_mat products(gens->m * O->m, n);
    for(unsigned int i = 0 ;  i < I0->m ; i++) {
        cxx_mpz_mat a(1,n);
        mpz_mat_submat_set(a,0,0,I0,i,0,1,n);
        mpz_mat_mul(a,a,M);
        for(unsigned int j = 0 ; j < n ; j++){
            mpz_mat_submat_swap(a,0,j*n,products,i*n+j,0,1,n);
        }
    }
    /* And put this in HNF */
    mpz_mat_hnf_backend_rev(products, NULL);
    cxx_mpz_mat I(n,n);
    mpz_mat_submat_swap(I,0,0,products,0,0,n,n);
    return make_pair(I, denom);
}//}}}

int prime_ideal_inertia_degree(cxx_mpz_mat const& I)/*{{{*/
{
    unsigned int n = I->m;
    ASSERT_ALWAYS(I->n == n);
    int f = 0;
    for(unsigned int i = 0 ; i < n ; i++) {
        f += mpz_cmp_ui(mpz_mat_entry_const(I, i, i), 1) != 0;
    }
    return f;
}
/*}}}*/
int valuation_of_ideal_at_prime_ideal(cxx_mpz_mat const& M, cxx_mpz_mat const& I, cxx_mpz_mat const& a, cxx_mpz const& p)/*{{{*/
{
    unsigned int n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(I->m == n);
    ASSERT_ALWAYS(I->n == n);
    ASSERT_ALWAYS(a->m == 1);
    ASSERT_ALWAYS(a->n == n);
    
    cxx_mpz_mat Ia(I);
    int val = 0;
    for(;;val++) {
        for(unsigned int i = 0 ; i < n ; i++) {
            cxx_mpz_mat b(1,n);
            mpz_mat_submat_set(b,0,0,Ia,i,0,1,n);
            b = multiply_elements_in_order(M, a, b);
            mpz_mat_submat_swap(b,0,0,Ia,i,0,1,n);
        }
        if (mpz_mat_p_valuation(Ia, p) < 1)
            return val;
        mpz_mat_divexact_mpz(Ia, Ia, p);
    }
}
/*}}}*/

struct hypercube_walk {/*{{{*/
    struct iterator {
        int B;
        vector<int> v;
        vector<int> speed;
        pair<int, int> last;
        iterator& operator++() {
            for(unsigned int j = 0 ; j < v.size() ; j++) {
                int s = v[j] + speed[j]; 
                if (s >= 0 && s <= B) {
                    v[j] = s; 
                    last = make_pair(j, speed[j]);
                    return *this;
                }
                speed[j]=-speed[j];
            }
            last = make_pair(-1, -1);
            return *this;
        }
        inline bool operator!=(iterator const& x) const { return last != x.last; }
        vector<int> & operator*() { return v; }
    };
    int n,B;
    hypercube_walk(int n, int B) : n(n), B(B) {}
    iterator begin() const {
        iterator z;
        z.B = B; 
        z.v.assign(n, 0);
        z.speed.assign(n, 1);
        return z;
    }
    iterator middle(vector<int> const& x) const {
        ASSERT_ALWAYS(x.size() == (unsigned int) n);
        iterator z;
        z.B = B; 
        z.v = x;
        z.speed.assign(n, 1);
        return z;
    }
    iterator end() const {
        iterator z;
        z.B = B; 
        z.v.assign(n,0);
        z.speed.assign(n, 1);
        z.last = make_pair(-1,-1);
        return z;
    }
};/*}}}*/
 
/* {{{ prime_ideal_two_element */
/* I must be in HNF */
pair<cxx_mpz, cxx_mpz_mat> prime_ideal_two_element(cxx_mpq_mat const& O, cxx_mpz_poly const& f, cxx_mpz_mat const& M, cxx_mpz_mat const& I)
{
    cxx_mpz p;
    unsigned int n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(I->m == n);
    ASSERT_ALWAYS(I->n == n);
    int inertia = 0;
    for(unsigned int i = 0 ; i < n ; i++) {
        mpz_srcptr di = mpz_mat_entry_const(I, i, i);
        if (mpz_cmp_ui(di, 1) != 0) {
            mpz_set(p, di);
            inertia++;
        }
    }
    ASSERT_ALWAYS(mpz_size(p) != 0);

    /* We use a homebrew algorithm, loosely inspired on Cohen's. We
     * prefer an enumeration which favors small generators. Based on
     * [Cohen93, lemma 4.7.9], we're simply going to search for an
     * element whose norm has p-valuation exactly f.
     *
     * Note though that we're only achieving small element size in terms
     * of coordinates on the order.
     */
    /* There's a question about which generators we should combine in
     * order to find our winning generator. Cohen takes "the generators
     * of p*O, and the generators of I". In fact, we can simply take the
     * generators of I as a Z-basis: that generates the right set. Then,
     * we'd better work on some ordering for these. Elements in the HNF
     * form of I with a p on the diagonal have in fact all other elements
     * on that row equal to zero. That is because since I contains p*O,
     * the subtraction of this generating row and the p-th multiple of
     * the corresponding generator has to be zero, or I would not be in
     * HNF. Therefore, these elements are of secondary importance in
     * generating I (e.g.: alone, they can't).
     */
    unsigned int m = n; /* number of generators we take */
    cxx_mpz_mat OI(m,n);
    for(unsigned int i = 0, r = 0, s = 0 ; i < n ; i++) {
        mpz_srcptr di = mpz_mat_entry_const(I, i, i);
        if (mpz_cmp_ui(di, 1) != 0) {
            mpz_mat_submat_set(OI,n - inertia + s++,0,I,i,0,1,n);
        } else {
            mpz_mat_submat_set(OI,r++,0,I,i,0,1,n);
        }
    }

    cxx_mpq_mat OIOq;
    mpq_mat_mul(OIOq, cxx_mpq_mat(OI), O); 
    for(int B = 1 ; ; B++) {
        hypercube_walk H(m, B);
        vector<int> H0(m,0);
        cxx_mpq_mat gen(1, n);
        // H0[n]=1;
        // mpq_mat_submat_set(gen, 0, 0, OIOq, n, 0, 1, n);
        cxx_mpq_mat temp(1, n);
        for(hypercube_walk::iterator it = H.middle(H0) ; it != H.end() ; ++it) {
            int j = it.last.first;
            int s = it.last.second;
            /* add (s=1) or subtract (s=-1) the j-th generator to gen */
            mpq_mat_submat_swap(temp,0,0,OIOq,j,0,1,n);
            if (s==1)
                mpq_mat_add(gen, gen, temp);
            else if (s==-1)
                mpq_mat_sub(gen, gen, temp);
            mpq_mat_submat_swap(temp,0,0,OIOq,j,0,1,n);

            cxx_mpz_poly pgen;
            cxx_mpz dgen;
            cxx_mpz res;
            mpq_mat_row_to_poly(pgen, dgen, gen, 0);
            if (pgen->deg < 0) continue;
            mpz_poly_resultant(res, pgen, f);
            /* We want the absolute norm to have p-valuation equal to the
             * inertia degree.
             * absolute norm is galois norm. galois norm is product of
             * all conjugate of pgen(alpha). Resultant(f,pgen) is
             * lc(f)^deg(pgen) times the galois norm.
             */
            int v = mpz_p_valuation(res, p) - pgen->deg * mpz_p_valuation(f->coeff[f->deg], p) - f->deg * mpz_p_valuation(dgen, p);
            ASSERT_ALWAYS(v >= inertia);
            if (v == inertia) {
                cxx_mpz_mat lambda(1, m);
                vector<int> & v(*it);
                for(unsigned int i = 0 ; i < m ; i++) {
                    mpz_set_si(mpz_mat_entry(lambda,0,i),v[i]);
                }
                mpz_mat_mul(lambda, lambda, OI);
                return make_pair(p, lambda);
            }
        }
    }
}
// }}}

#if 0
// Roots contains elements of P^1 (Z/pZ)
// This function leaves each element only once (e.g. if roots contains 1/1 and 4/4, then 4/4 will be deleted from the vector)
void filter_roots(vector<pair<cxx_mpz, cxx_mpz>>& roots, unsigned long p)
{
    unsigned int i = 0;

    // Filtering the root
    while(i < roots.size()-1){
        // We look for roots[j] equal to roots[i], for j > i
        unsigned int j = i+1;
        while (j < roots.size()){
            cxx_mpz a = roots[j].first;
            cxx_mpz b = roots[j].second;
            
            // At this point, we will multiply a and b with 2, then a and
            // b with 3, etc...  if we find one k such that a*k =
            // roots[i].first and b*k = roots[i].second, then it's the
            // same quotient Thus we can erase a/b
            cxx_mpz a1, b1;
            unsigned int k = 2;
            bool eq = false;
            while( (k < p) && (!eq) ){
                a1 = a;
                b1 = b;
                
                mpz_mul_ui(a1, a1, k);
                mpz_mul_ui(b1, b1, k);
                mpz_mod_ui(a1, a1, p);
                mpz_mod_ui(b1, b1, p);
                if( (mpz_cmp(a1,roots[i].first) == 0) && (mpz_cmp(b1,roots[i].second) == 0) ){
                    eq = true;
                }
                k++;
            }
            
            if(eq){
                roots.erase(roots.begin()+j);
            }
            else{
                j++;
            }
        }
        
        i++;
    }
    
    // Putting each root on the form (_/1)
    for (i = 0 ; i < roots.size() ; i++){
        if(mpz_cmp_ui(roots[i].second,1) != 0){
            cxx_mpz a1, b1;
            unsigned int k = 2;
            bool eq = false;
            while( (k < p) && (!eq) ){
                a1 = roots[i].first;
                b1 = roots[i].second;

                mpz_mul_ui(a1, a1, k);
                mpz_mul_ui(b1, b1, k);
                mpz_mod_ui(a1, a1, p);
                mpz_mod_ui(b1, b1, p);
                if(mpz_cmp_ui(b1,1) == 0 ){
                    eq = true;
                }
                k++;
            }
            mpz_set(roots[i].first,a1);
            mpz_set(roots[i].second,b1);
        }
    }
}

void print_comments_for_badideals_above_p(string& SBAD,
        const unsigned int side, 
        mpq_mat_srcptr order, mpz_poly_srcptr f, 
        vector<pair<cxx_mpq_mat, int>> ideals, const unsigned long p)
{
    ASSERT_ALWAYS(order->m == order->n);
    int n = order->m;
    
    // Computing the inertia degree of each ideal
    vector<int> inertia;
    for(unsigned int i = 0 ; i < ideals.size() ; i++){
        inertia.push_back(prime_ideal_inertia_degree(order, ideals[i].first));
        printf("Inertia degree of ideal %d is %d\n", i, inertia[i]);
    }
    
    
    // Building the ideal jj (J^-1 where J is <1, alpha>^-1)
    cxx_mpz_mat jj_gen;
    cxx_mpq_mat jj;
    
    mpz_mat_realloc(jj_gen,1,n);
    mpz_set_ui(mpz_mat_entry(jj_gen,0,0),1);
    mpz_set_ui(mpz_mat_entry(jj_gen,0,1),1);
    make_ideal(jj,order,f,jj_gen);
       
    
    // Now listing all roots of the homogenous polynomial build from f
    vector<pair<cxx_mpz, cxx_mpz>> rootsp;
    for (unsigned int i = 0 ; i < p ; i++){
        for (unsigned int j = 1 ; j < p ; j++){
            cxx_mpz res, a, b;
            mpz_set_ui(a,i);
            mpz_set_ui(b,j);
            
            mpz_poly_homogeneous_eval_siui(res,f,i,j);
            mpz_mod_ui(res,res,p);
            if(mpz_congruent_ui_p(res,0,p)){
                pair<cxx_mpz,cxx_mpz> new_elem;
                new_elem = make_pair(a,b);
                rootsp.push_back(new_elem);
            }

            
        }
    }
    
    // Since each roots is listed several times under different forms (ex
    // : if 1/1 is root, 2/2 is too, etc...)
    // we filter roots here to keep each root only once
    filter_roots(rootsp, p); 
    
    // Testing if (1 : 0) is a root
    cxx_mpz fd;
    mpz_poly_getcoeff(fd,f->deg,f);
    if(mpz_congruent_ui_p(fd,0,p)){
        cxx_mpz a,b;
        mpz_set_ui(a,1);
        mpz_set_ui(b,0);
        pair<cxx_mpz,cxx_mpz> new_elem;
        new_elem = make_pair(a,b);
        rootsp.push_back(new_elem);
    }
    
    printf("Roots :\n");
    for(unsigned int i = 0 ; i < rootsp.size() ; i++){
        gmp_printf("(%Zd : %Zd)\n", rootsp[i].first, rootsp[i].second);
    }
    
    // Now going to compute the valuation of each ideal <p, a-b*alpha> on
    // each ideals above p, where a/b belong to the set of roots
    for(unsigned int i = 0 ; i < rootsp.size() ; i++){
        cxx_mpz u = rootsp[i].first;
        cxx_mpz v = rootsp[i].second;
        
        // Computing the ideal ii, e.g. <p,(v*alpha-u)>*J
        cxx_mpz_mat gens_ii;
        cxx_mpq_mat ii;
        mpz_mat_realloc(gens_ii,3,n);
        mpz_set_ui(mpz_mat_entry(gens_ii,0,0),p);
        mpz_set_ui(mpz_mat_entry(gens_ii,1,1),p);
        mpz_set(mpz_mat_entry(gens_ii,2,0),u); mpz_mul_si(mpz_mat_entry(gens_ii,2,0),mpz_mat_entry(gens_ii,2,0),-1);
        mpz_set(mpz_mat_entry(gens_ii,2,1),v);
        make_ideal(ii, order, f, gens_ii);

        //printf("ii is\n"); mpq_mat_fprint(stdout, ii); printf("\n");
        
        
        
        
        // Computing the valuation of ii on each prime ideal above p
        vector<int> vals;
        for(unsigned int j = 0 ; j < ideals.size(); j++){
            vals.push_back(valuation_of_ideal_at_prime_ideal(order, f, ii, ideals[j].first, ideals[j].second, p));
            gmp_printf("valuation of (%Zd : %Zd) on ideal %d is : %d\n", u, v, j, vals.back());
        }
        printf("\n");
        
        // Detecting all indices for which vals[i] != 0
        vector<int> indices;
        for(unsigned int j = 0 ; j < vals.size(); j++){
            if(vals[j] != 0){indices.push_back(j);}
        }
        
        // If only one ideal holds the valuation, there's nothing to do
        if(indices.size() == 1){
            continue;
        }
        
        // Normalisation of the root
        int a;
        pair<cxx_mpz, cxx_mpz> Q;
        Q = rootsp[i];
        if(mpz_cmp_ui(Q.second, 0) != 0){
            a = mpz_get_ui(Q.first);
        }
        else{
            a = p;
        }
        
        std::stringstream stream;
        stream << std::hex << p << "," << a;
        std::string result( stream.str() );

        SBAD += (result+":"+to_string(side)+", "+to_string(indices.size()));
        cout << SBAD << endl;
        
    }
}
#endif
