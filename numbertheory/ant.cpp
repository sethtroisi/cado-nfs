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
#include "macros.h"
#include "mpz_poly.h"
#include "mpz_mat.h"
#include "utils/cxx_mpz.hpp"
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
    return s;
}//}}}

cxx_mpz_mat join_HNF(cxx_mpz_mat const& K, const unsigned long p)//{{{
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

cxx_mpz_mat join_HNF_rev(cxx_mpz_mat const& K, const unsigned long p)//{{{
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
cxx_mpz_mat frobenius_matrix(cxx_mpz_mat const& M, unsigned long p)
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
    mpz_mat_mod_ui(Mp, Mp, p);

    cxx_mpz_mat E(n, n, 1), F(E);

    unsigned long k = ((~0UL)>>1) + 1;
    for( ; k > p ; k >>= 1);
    for( ; k >>= 1 ; ) {
        F = multiply_elements_in_order(M, F, F);
        if (p & k)
            F = multiply_elements_in_order(M, E, F);
        mpz_mat_mod_ui(F, F, p);
    }
    return F;
}//}}}

// {{{ cxx_mpz_mat p_radical_of_order
// Stores in I the p-radical of the order whose multiplication matrix is
// given by M. I is expressed with respect to the basis of the order.
cxx_mpz_mat p_radical_of_order(cxx_mpz_mat const& M, unsigned long p)
{
    unsigned int n = M->m;
    ASSERT_ALWAYS(M->n == n * n);

    cxx_mpz_mat K;

    // Now building the matrix U, containing all generators to the power
    // of p Generators are polynomials, stored in the matrix B
    cxx_mpz_mat T = frobenius_matrix(M, p);

    /* Which power of p ? */
    int k = 1;
    for (unsigned int pk = p; pk < n; k++, pk *= p);
    mpz_mat_pow_ui_mod_ui(T, T, k, p);

    // Storing in K a basis of Ker((z -> (z^%p mod g mod p))^k)
    mpz_mat_kernel_mod_ui(K, T, p);

    // Getting generators of the p radical from Ker(X) by computing HNF
    // of the vertical block matrix (p*Id, K);
    return join_HNF(K, p);
}// }}}

// {{{ cxx_mpq_mat p_maximal_order(cxx_mpz_poly const& f, unsigned long p)
// Builds one p-maximal-order starting from the matrix B, containing the
// generators of one order (in the basis of alpha^) In the number field
// of the polynomial f, of degree n Stores the generators of this
// p-maximal-order (in the basis of alpha^) in D
cxx_mpq_mat p_maximal_order(cxx_mpz_poly const& f, unsigned long p)
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
        mpz_mat_mod_ui(M2, M2, p);
        //printf("M is :\n"); mpz_mat_fprint(stdout, M); printf("\n");
        
        // Computing Ker(M)
        cxx_mpz_mat K;	        // Ker(M)
        mpz_mat_kernel_mod_ui(K, M2, p);
        // printf("Ker(M) is :\n"); mpz_mat_fprint(stdout, K_M); printf("\n");

        // Getting generators of p*O' by computing HNF of the vertical block matrix (p*Id, K_M);
        cxx_mpz_mat J = join_HNF(K, p);

        // Converting in the basis which is used to express elements of D
        mpq_mat_mul(new_D, cxx_mpq_mat(J), D);

        // Dividing by p
        mpq_mat_div_ui(new_D, new_D, p);

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
cxx_mpz_poly mpz_mat_minpoly_mod_ui(cxx_mpz_mat M, unsigned long p)
{
    ASSERT_ALWAYS(M->m == M->n);
    unsigned long n = M->n;
    cxx_mpz_mat B(n+1, n*n);
    cxx_mpz_mat T(n,n,1);
    for(unsigned int k = 0 ; k < n + 1 ; k++) {
        for(unsigned int i = 0 ; i < n ; i++) {
            mpz_mat_submat_set(B, n - k, i * n, T, i, 0, 1, n);
        }
        mpz_mat_mul_mod_ui(T, T, M, p);
    }
    cxx_mpz_mat K;
    mpz_mat_kernel_mod_ui(K, B, p);
    mpz_mat_gauss_backend_mod_ui(K, NULL, p);

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
        unsigned long p)
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
    mpz_mat_gauss_backend_mod_ui(Jt, T, p);
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
    mpz_mat_mul_mod_ui(Mc, CJ, Kt, p);
    return Mc;
}
/*}}}*/

/*{{{ factorization_of_polynomial_mod_ui */
vector<pair<cxx_mpz_poly, int>> factorization_of_polynomial_mod_ui(cxx_mpz_poly const& f, unsigned long p, gmp_randstate_t state)
{
    cxx_mpz pz;
    mpz_set_ui(pz,p);
    mpz_poly_factor_list lf;
    mpz_poly_factor_list_init(lf);
    mpz_poly_factor(lf,f,pz,state);
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
        unsigned long p,
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
        mpz_set_ui(mpz_mat_entry(c,0,i), gmp_urandomm_ui(state, p));
    }
    mpz_mat_mul_mod_ui(c, c, J, p);

    cxx_mpz_mat Mc = matrix_of_multmap(M, J, c, p);

    /* Now we would like to find the minimal polynomial of Mc, which is
     * simply the matrix of multiplication by c. For this, we write an
     * (m+1) times m^2 matrix, and compute a left kernel.
     */
    cxx_mpz_poly Pc = mpz_mat_minpoly_mod_ui(Mc, p);

    vector<pair<cxx_mpz_poly, int>> facP = factorization_of_polynomial_mod_ui(Pc, p, state);

    vector<pair<cxx_mpz_mat, int>> ideals;

    vector<cxx_mpz_mat> characteristic_subspaces;

    for(unsigned int i = 0 ; i < facP.size() ; i++) {
        cxx_mpz_poly const& f(facP[i].first);
        /* We need the basis of the kernel of f(Mc) */
        cxx_mpz_mat E;
        mpz_poly_eval_mpz_mat_mod_ui(E, Mc, f, p);
        mpz_mat_pow_ui_mod_ui(E, E, facP[i].second, p);
        mpz_mat_kernel_mod_ui(E, E, p);
        /* This line is just to be exactly in line with what magma says
         */
        mpz_mat_gauss_backend_mod_ui(E, NULL, p);
        mpz_mat_mul_mod_ui(E, E, J, p);
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
            mpz_mat_submat_swap(Ihead,0,0,Ix,Ix->m-n,0,n,n);
            ideals.push_back(make_pair(Ihead, e));
        } else {
            mpz_mat_hnf_backend_rev(Ix, NULL);
            mpz_mat_submat_swap(Ihead,0,0,Ix,Ix->m-n,0,n,n);
            vector<pair<cxx_mpz_mat, int>> more_ideals;
            more_ideals = factorization_of_prime_inner(B,M,p,Ip,Ihead,Ci,state);
            append_move(ideals, more_ideals);
        }
    }
    return ideals;
}

vector<pair<cxx_mpz_mat, int>> factorization_of_prime(
        cxx_mpq_mat & B, cxx_mpz_poly const& g,
        unsigned long p,
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

// Let O be the order generated by elements of G.
// This stores in I the ideal of O generated by gen (supposed to be in
// a-b*alpha form)
//
void make_principal_ideal(mpq_mat_ptr I, mpq_mat_srcptr G, mpz_poly_srcptr g, mpz_mat_srcptr gen)
{
    ASSERT_ALWAYS((G->m == G->n) && (gen->m == 1) && (gen->n > 0) && (gen->n <= G->n) && (g->deg == (int) G->n));
    unsigned int n = G->m;
    
    cxx_mpz_poly aux1;
    mpz_poly_realloc(aux1,n);
    mpq_mat_realloc(I,n,n);
    mpz_mat_row_to_poly(aux1,gen,0);
    
    for(unsigned int j = 0 ; j < n ; j++){
        cxx_mpz_poly aux2, aux3;
        mpz_t denom;
        mpz_init(denom);
        mpz_poly_realloc(aux3, n);
        mpz_poly_realloc(aux2, n);
            
        mpq_mat_row_to_poly(aux2, denom, G, j);
            
        mpz_poly_mul_mod_f(aux3, aux1, aux2, g);
            
        cxx_mpq_mat aux_mat;
        mpq_mat_realloc(aux_mat,1,n);
        mpq_poly_to_mat_row(aux_mat, 0, aux3, denom);
        
        mpq_mat_submat_swap(I,j,0,aux_mat,0,0,1,n);
        
        mpz_clear(denom);
    }

    // hnf_magma_style(I,I);
    // mpz_mat_hnf_backend_rev(I, NULL);
}

// Let O be the order generated by elements of G.
// This stores in I the ideal of O generated by rows of gen
void make_ideal(mpq_mat_ptr I, mpq_mat_srcptr G, mpz_poly_srcptr g, mpz_mat_srcptr gen)
{
    ASSERT_ALWAYS((G->m == G->n) && (gen->m > 0) && (gen->n > 0) && (gen->n <= G->n) && (g->deg == (int) G->n));
    unsigned int n = G->m;
    
    mpq_mat_realloc(I,n*gen->m,n);
    
    for(unsigned int i = 0 ; i < gen->m ; i++){
        cxx_mpz_poly aux1;
        mpz_poly_realloc(aux1,n);
        mpz_mat_row_to_poly(aux1,gen,i);
        for(unsigned int j = 0 ; j < n ; j++){
            cxx_mpz_poly aux2, aux3;
            mpz_t denom;
            mpz_init(denom);
            mpz_poly_realloc(aux3, n);
            mpz_poly_realloc(aux2, n);
            
            mpq_mat_row_to_poly(aux2, denom, G, j);
            mpz_poly_mul_mod_f(aux3, aux1, aux2, g);
            
            cxx_mpq_mat aux_mat;
            mpq_mat_realloc(aux_mat,1,n);
            mpq_poly_to_mat_row(aux_mat, 0, aux3, denom);
        
            mpq_mat_submat_swap(I,j+n*i,0,aux_mat,0,0,1,n);
        
            mpz_clear(denom);
        }
    }
    
    // hnf_magma_style(I,I);
    
    cxx_mpq_mat J;
    mpq_mat_realloc(J,n,n);
    mpq_mat_submat_swap(J,0,0,I,I->m-n,0,n,n);
    mpq_mat_swap(I,J);
}

// G is one p-maximal order (on the basis of alpha^) on the number field
// of g, I one ideal (on alpha^)
// Finds a such that (a/p)*I is in the order of G (a is one element of
// the order) ; stores a/p in h
void valuation_helper_for_ideal(mpq_mat_ptr h, mpq_mat_srcptr G,
        mpz_poly_srcptr g, mpq_mat_srcptr I, const unsigned long p)
{
    ASSERT_ALWAYS((G->m == G->n) && (I->m == I->n) && (G->m == I->m));
    unsigned int n = G->m;
    
    cxx_mpz_mat C, ker;
    mpz_mat_realloc(C,n,n*n);
    
    // Filling C with the products G[i]*I[j] (being polynomials on alpha^)
    for(unsigned int i = 0 ; i < n ; i++){
        for(unsigned int j = 0 ; j < n ; j++){
            cxx_mpz_poly aux1, aux2, aux3;
            mpz_t denom1, denom2, denom3;
            mpz_init(denom1);
            mpz_init(denom2);
            mpz_init(denom3);
            mpz_poly_realloc(aux1,n);
            mpz_poly_realloc(aux2,n);
            mpz_poly_realloc(aux3,n);
            
            mpq_mat_row_to_poly(aux1, denom1, G, i);
            mpq_mat_row_to_poly(aux2, denom2, I, j);
            mpz_poly_cleandeg(aux1, n-1);
            mpz_poly_cleandeg(aux2, n-1);

            mpz_poly_mul_mod_f(aux3, aux1, aux2, g);
            mpz_mul(denom3,denom2,denom1);
            mpz_poly_cleandeg(aux3, n-1);

            cxx_mpq_mat aux_mat;
            mpq_mat_realloc(aux_mat, 1, n);
            mpq_poly_to_mat_row(aux_mat, 0, aux3, denom3);

            cxx_mpq_mat G_inv;
            mpq_mat_realloc(G_inv,n,n);
            mpq_mat_inv(G_inv,G);
            mpq_mat_mul(aux_mat, aux_mat, G_inv);
            
            cxx_mpz_mat aux_mat_int;
            mpz_mat_realloc(aux_mat_int,1,n);
            mpq_mat_numden(aux_mat_int,NULL,aux_mat);
            
            mpz_mat_submat_swap(aux_mat_int, 0, 0, C, i, n*j, 1, n);
            
            mpz_clear(denom1);
            mpz_clear(denom2);
            mpz_clear(denom3);
        }
    }
    
    // Reducing mod p and computing the kernel
    mpz_mat_mod_ui(C, C, p);
    mpz_mat_kernel_mod_ui(ker, C, p);
    mpz_mat_hnf_backend(ker,C); // We don't need C anymore, so we use it to contain the transformation matrix
    
    cxx_mpq_mat ker_rat;
    mpq_t p_inv;
    mpq_init(p_inv);
    mpq_set_ui(p_inv,1,p);
    mpq_mat_realloc(ker_rat, 1, n);
    
    // Turning the vectors of the kernel in the basis of alpha^, and dividing by p
    mpq_mat_set_mpz_mat(ker_rat, ker);
    mpq_mat_mul(ker_rat, ker_rat, G);
    mpq_mat_mul_mpq(ker_rat, ker_rat, p_inv);
        
    mpq_mat_realloc(h, 1, n);
    mpq_mat_submat_swap(h, 0, 0, ker_rat, 0, 0, 1, n);
    
    mpq_clear(p_inv);
}

// g is the polynomial generating the number field
// G contains generators of order O in basis of alpha^
// I is a prime ideal of O above p, on the basis of alpha^
// e is its ramification index
// J is an ideal on the basis of alpha^
// WARNING : I only tested it when J is a principal ideal
int valuation_of_ideal_at_prime_ideal(mpq_mat_srcptr G, mpz_poly_srcptr g, mpq_mat_srcptr J, mpq_mat_srcptr I, const int e, const unsigned long p)
{
    ASSERT_ALWAYS((G->m == G->n) && (G->n == J->m) && (J->m == J->n) && (J->n == I->m) && (I->m == I->n));
    unsigned int n = G->m;
    
    cxx_mpq_mat G_inv, MJ, h;
    cxx_mpz_mat MJ_int;
    mpz_t dJ, denom;
    mpz_init(dJ); // Initial denominator of MJ
    mpz_init(denom); // Denominator of MJ for the loop
    mpq_mat_realloc(G_inv,n,n);
    mpq_mat_inv(G_inv,G);
    mpq_mat_realloc(MJ,n,n); // MJ is the MJ existing in magma
    mpz_mat_realloc(MJ_int,n,n); // The MJ_int here is only MJ*dJ
    
    mpq_mat_mul(MJ, J, G_inv);
    mpq_mat_numden(MJ_int, dJ, MJ); // Now MJ_int = dJ*MJ, like in magma
    mpq_mat_set_mpz_mat(MJ, MJ_int); // We need a rational matrix because we have to compute its denominator

    
    // Now dJ contains the denominator MJ_rat had at the beginning (most
    // of the time, 1) And MJ contains integers, (former values of MJ *
    // dJ)
    
    valuation_helper_for_ideal(h, G, g, I, p);
    //printf("h = \n"); mpq_mat_fprint(stdout, h); printf("\n");
    int v = -1;
    //printf("MJ :\n"); mpq_mat_fprint(stdout, MJ); printf("\n");
    
    mpz_set(denom, dJ);
    do{
        //mpq_mat_fprint(stdout, MJ); printf("\n");
        //printf("v = %d\n", v);
        //gmp_printf("denom is %Zd\n", denom);
        cxx_mpq_mat aux, h_times_MJ_G;
        mpq_mat_realloc(h_times_MJ_G, n, n);
        mpq_mat_realloc(aux,n,n);
        mpq_mat_mul(aux, MJ, G); // This way, you obtain all rows of the form (Vector(K,j),Vector(G))
        // Now, these rows (in aux) are polynomials, and have to be
        // multiplied with h ; then, again in a matrix, and times G^-1

        for(unsigned int j = 0 ; j < n ; j++){
            cxx_mpz_poly aux1, aux2, aux3;
            mpz_t denom1, denom2, denom3;
            mpz_init(denom1);
            mpz_init(denom2);
            mpz_init(denom3);
            mpz_poly_realloc(aux1, n);
            mpz_poly_realloc(aux2, n);
            
            mpq_mat_row_to_poly(aux1, denom1, aux, j);
            mpq_mat_row_to_poly(aux2, denom2, h, 0);
            
            mpz_poly_mul_mod_f(aux3, aux1, aux2, g);
            mpz_mul(denom3, denom1, denom2);
            
            cxx_mpq_mat aux_mat;
            mpq_mat_realloc(aux_mat,1,n);
            mpq_poly_to_mat_row(aux_mat, 0, aux3, denom3);
            
            mpq_mat_submat_swap(h_times_MJ_G,j,0,aux_mat,0,0,1,n);
        
            mpz_clear(denom1);
            mpz_clear(denom2);
            mpz_clear(denom3);
        }
        // h_times_MJ_G contains the polynomials multiplied by h ; now,
        // MJ must be h_times_MJ_G * G^-1 Thus MJ contains only h * MJ
        //printf("h_times_MJ_G :\n"); mpq_mat_fprint(stdout,
        //h_times_MJ_G); printf("\n");
        mpq_mat_mul(MJ, h_times_MJ_G, G_inv);
        //printf("new MJ :\n"); mpq_mat_fprint(stdout, MJ); printf("\n");
        //int k;
        //scanf("%d",&k);
        
        cxx_mpz_mat trash;
        mpz_mat_realloc(trash,n,n);
        mpq_mat_numden(trash, denom, MJ);
        
        v++;
    }while(mpz_cmp_ui(denom,1) == 0); // while denom == 1
    
    int dJ_int = mpz_get_si(dJ);
    int val_dJ_on_p = 0;
    while(dJ_int%p == 0){
        dJ_int = dJ_int/p;
        val_dJ_on_p += 1;
    }
    
    mpz_clear(dJ);
    mpz_clear(denom);
    
    return v-val_dJ_on_p*e;
    
}

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

// Takes G, basis of one order, in alpha^ basis, and I, basis of an
// ideal, in alpha^ basis
// Returns inertia degree of I
unsigned int prime_ideal_inertia_degree(mpq_mat_srcptr G, mpq_mat_srcptr I)
{
    
    ASSERT_ALWAYS((G->m == G->n) && (I->n == I->m) && (I->n == G->m));
    int n = G->m;
    
    cxx_mpq_mat tmp, G_inv;
    cxx_mpz_mat coords;
    mpq_mat_realloc(G_inv, n, n);
    
    mpq_mat_inv(G_inv, G);
    mpq_mat_mul(tmp,I, G_inv);
    
    mpq_mat_numden(coords, NULL, tmp);
    
    mpz_t d;
    mpz_init(d);
    mpz_mat_determinant_triangular(d,coords);
    int a = mpz_get_ui(d);
    int p = mpz_get_ui(mpq_numref(mpq_mat_entry_const(I,0,0)));
    mpz_clear(d);
    
    
    return (int) (log((double) a)/log((double) p));
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

/*
void factorization_of_arbitrary_ideal(vector<pair<cxx_mpq_mat, int>>& ideals, mpz_poly_srcptr g, 
    mpz_mat_srcptr gens, gmp_randstate_t state)
{
    ASSERT_ALWAYS((gens->m == 1) && (gens->n == (unsigned int) g->deg));
    unsigned long int a = mpz_get_ui(mpz_mat_entry_const(gens,0,0));
    unsigned long int b = mpz_get_ui(mpz_mat_entry_const(gens,0,1));
    mpz_t res;
    mpz_init(res);
    
    mpz_poly_homogeneous_eval_siui(res,g,a,b); // Res(g,a-bX)
    
    mpz_clear(res);
}
*/

void sort_matrices(vector<pair<cxx_mpq_mat, int>>& ideals)
{
    unsigned int n = ideals.size();
    for(unsigned int i = n-2 ; (int) i >= 0 ; i--){
        for(unsigned int j = 0 ; j <= i ; j++){
            
            if(mpq_mat_cmp(ideals[j].first,ideals[j+1].first) > 0){
                pair<cxx_mpq_mat, int> k;
                k = ideals[j];
                ideals[j] = ideals[j+1];
                ideals[j+1] = k; 
            }
        }
    }
}

