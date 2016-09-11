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
	mpz_poly_clear(f);
	mpz_poly_init(f,M->n-1);
	unsigned int j;
	for (j = 0 ; j < M->n; j++){
		mpz_poly_setcoeff(f,j,mpz_mat_entry_const(M,i,j));
	}
}

void mpz_mat_row_to_poly_rev(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int i)
{
	mpz_poly_clear(f);
	mpz_poly_init(f,M->n-1);
	unsigned int j;
	for (j = 0 ; j < M->n; j++){
		mpz_poly_setcoeff(f,j,mpz_mat_entry_const(M,i,M->n-1-j));
	}
}

void mpz_mat_column_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int j)
{
	mpz_poly_clear(f);
	mpz_poly_init(f,M->m-1);
	unsigned int i;
	for (i = 0 ; i < M->m; i++){
		mpz_poly_setcoeff(f,i,mpz_mat_entry_const(M,i,j));
	}
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

void mpq_mat_column_to_poly(mpz_poly_ptr f, mpz_ptr denom, mpq_mat_srcptr M, const unsigned int j)
{
	mpz_poly_clear(f);
	mpz_poly_init(f,M->m-1);
	mpz_set_si(denom,1);
	unsigned int i;
    for (i = 0 ; i < M->m ; i++) {
        mpz_lcm(denom,denom,mpq_denref(mpq_mat_entry_const(M,i,j)));
    }
	for (i = 0 ; i < M->m; i++){
		mpq_t aux;
		mpz_t num;
		mpq_init(aux);
		mpz_init(num);
		mpq_set_z(aux,denom);
		mpq_mul(aux,aux,mpq_mat_entry_const(M,i,j));
		mpz_set(num,mpq_numref(aux));
		mpz_poly_setcoeff(f,i,num);
		mpq_clear(aux);
		mpz_init(num);
	}
}
/*}}}*/

/*{{{ commonly used wrappers around HNF functions */
int mpz_hnf_backend_rev(mpz_mat_ptr M, mpz_mat_ptr T) // {{{
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
    int s = mpz_hnf_backend(M, T);
    mpz_mat_reverse_rows(M, M);
    mpz_mat_reverse_columns(M, M);
    mpz_mat_reverse_rows(T, T);
    mpz_mat_reverse_columns(T, T);
    return s;
}//}}}

cxx_mpz_mat join_HNF(cxx_mpz_mat const& K, const unsigned long p)//{{{
{
    cxx_mpz_mat J(K->n, K->n, p);
    cxx_mpz_mat T0;
    cxx_mpz_mat I;

    mpz_mat_vertical_join(I, J, K);
    mpz_hnf_backend(I, T0);
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
    mpz_hnf_backend_rev(I, T0);
    mpz_mat_submat_swap(I, 0, 0, J, 0, 0, K->n, K->n);
    return J;
}//}}}
/*}}}*/

cxx_mpz_mat multiplication_table_of_order(cxx_mpq_mat const& O,/*{{{*/
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

cxx_mpz_mat multiplication_table_of_ideal(cxx_mpz_mat const& M,/*{{{*/
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

cxx_mpz_mat multiply_elements_in_order(cxx_mpz_mat const& M, cxx_mpz_mat const& E, cxx_mpz_mat const& F)/*{{{*/
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

cxx_mpz_mat frobenius_matrix(cxx_mpz_mat const& M, unsigned long p) //{{{
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

// {{{ cxx_mpz_mat p_radical_of_order(cxx_mpq_mat const& B, cxx_mpz_poly const& g, unsigned long p)
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
    mpz_hnf_backend_rev(Dz, cxx_mpz_mat());
    mpq_mat_set_mpz_mat_denom(D, Dz, den);
    return D;
}
//}}}

/* theta is expected to be an element of the order whose basis is given
 * by W. Thus we read theta as a vector with integer coordinates (matrix
 * of size 1*n)
 */
void matrix_of_multiplication_theta_local(mpz_mat_ptr M, mpq_mat_srcptr W, mpz_mat_srcptr theta, mpz_poly_srcptr g, const unsigned long p)
{
    unsigned int n = W->m;
    ASSERT_ALWAYS(W->m == W->n);
    ASSERT_ALWAYS(theta->m == 1 && theta->n == n);
    
    // First let's build the matrix of multiplication by theta, in the
    // basis of (W[0], ... W[n-1])
    mpz_mat times_theta;
    mpq_mat times_theta_rat;
    mpz_mat_init(times_theta,n,n);
    mpq_mat_init(times_theta_rat,n,n);
    
    mpq_mat theta_rat; // Will contain theta in the basis of alpha^
    mpq_mat_init(theta_rat,1,n);
    mpq_mat_set_mpz_mat(theta_rat,theta);
    mpq_mat_mul(theta_rat,theta_rat,W); // Now contains theta in the basis of alpha^ (rational coefficients)

    // Converting theta into one polynomial, with a denominator
    mpz_poly theta_poly;
    mpz_t theta_denom;
    mpz_poly_init(theta_poly,n-1);
    mpz_init(theta_denom);
    mpq_mat_row_to_poly(theta_poly,theta_denom,theta_rat,0);
    
    // Computing theta*w[i] (in the basis of alpha^ ) for each i
    for (unsigned i = 0 ; i < n ; i++) {
        // Converting w[i] (already in the basis of alpha^) into one
        // polynom and one common denominator
        mpz_poly w_poly;
        mpz_poly_init(w_poly,n-1);

        mpz_t w_denom;
        mpz_init(w_denom);

        mpq_mat_row_to_poly(w_poly,w_denom,W,i);

        mpz_t denom;
        mpz_init(denom);
        mpz_mul(denom,theta_denom,w_denom);
        
        mpz_poly res;
        mpz_poly_init(res,n-1);        
        mpz_poly_mul_mod_f(res,theta_poly,w_poly,g);
        mpz_poly_cleandeg(res,n-1);
        
        mpq_poly_to_mat_row(times_theta_rat, i, res, denom);

        mpz_poly_clear(res);
        mpz_clear(w_denom);
        mpz_poly_clear(w_poly);
        mpz_clear(denom);
    }
    
    
    // Now we have to multiply by W^-1, in order to get into the basis of
    // (w[0], ... w[n-1]) again
    mpq_mat W_inv;
    mpq_mat_init(W_inv,n,n);
    mpq_mat_inv(W_inv,W);
    mpq_mat_mul(times_theta_rat,times_theta_rat,W_inv);
    
    // Now we have to convert it into a matrix of integers
    int rc = mpq_mat_numden(times_theta, NULL, times_theta_rat);
    ASSERT_ALWAYS(rc == 1);
    
    // And mod p
    mpz_mat_mod_ui(times_theta,times_theta,p);

    mpz_mat_swap(M, times_theta);

    mpq_mat_clear(W_inv);
    mpz_clear(theta_denom);
    mpz_poly_clear(theta_poly);
    mpq_mat_clear(theta_rat);
    mpq_mat_clear(times_theta_rat);
    mpz_mat_clear(times_theta);
}


// W is the matrix containing the generators of one p-maximal order Ok,
// in the basis of alpha^ (root of g) theta is a row matrix representing
// an element of Ok/pOk, thus a vector on the basis (W[0], ... W[n-1])
// (with coeffs mod p)
// g is the monic polynom defining the number field in which we are
void minimal_poly_of_mul_theta(mpz_poly_ptr f, mpq_mat_srcptr W, mpz_mat_srcptr theta, mpz_poly_srcptr g, const unsigned long p)
{
    unsigned int n = W->m;
    
    mpz_mat times_theta;
    mpz_mat_init(times_theta,n,n);

    matrix_of_multiplication_theta_local(times_theta, W, theta, g, p);
    
    // Now starting to compute the (n+1,n^2) matrix whose kernel will be computed
    mpz_mat M,current;
    mpz_mat_init(M,n+1,n*n);
    mpz_mat_init(current,n,n);
    mpz_mat_set_ui(current,1);
    for (unsigned int k = 0 ; k < n+1 ; k++){
        for (unsigned int i = 0 ; i < n ; i++){
            for (unsigned int j = 0 ; j < n ; j++){
                mpz_set(mpz_mat_entry(M,n-k,j+n*i),mpz_mat_entry(current,i,j));
            }
        }
        mpz_mat_mod_ui(M,M,p);
        mpz_mat_mul(current,current,times_theta);
        mpz_mat_mod_ui(current,current,p);
    }


    // Now computing its kernel
    mpz_mat K;
    mpz_mat_init(K,0,0);
    mpz_mat_kernel_mod_ui(K,M,p);

    mpz_gauss_backend_mod_ui(K,NULL,p);
    
    // Getting the minimal polynomial and verifying that f(M) = 0
    if(K->m == 0){
        mpz_poly_realloc(f,0);
    } else {
        mpz_poly_set_zero(f);
        mpz_mat_row_to_poly_rev(f,K,K->m-1);
    }
    
    mpz_mat_clear(K);
    mpz_mat_clear(current);
    mpz_mat_clear(M);
    mpz_mat_clear(times_theta);
}

// Represents the type of elements found in pick_from, in badideals.mag
struct subspace_ideal {
    cxx_mpz_mat E; // the basis of the subspace E
    cxx_mpz_mat I; // the basis of one ideal
    cxx_mpz_mat V; // same as E, in fact
};

// U and W are matrices containing generators (in rows) of vector subspaces of one big vector space
// This computes the basis of the intersection of those subspaces
void intersection_of_subspaces_mod_ui(mpz_mat_ptr M, mpz_mat_srcptr U, mpz_mat_srcptr W, const unsigned long p)
{
    ASSERT_ALWAYS(U->n == W->n);
    mpz_mat U_t, W_t; // U and W, transposed
    mpz_mat_init(U_t, U->n, U->m);
    mpz_mat_init(W_t, W->n, W->m);
    
    mpz_mat_transpose(U_t,U);
    mpz_mat_transpose(W_t,W);
    
    mpz_mat X, Y; // Left kernel of U_t and W_t
    mpz_mat_init(X, 0, U->n);
    mpz_mat_init(Y, 0, W->n);
    
    mpz_mat_kernel_mod_ui(X,U_t,p); // X = ker(x -> x*(transpose of U)
    mpz_mat_kernel_mod_ui(Y,W_t,p); // Y = ker(y -> y*(transpose of W)
    
    mpz_mat Z, Z_t; // Block matrix, X on the top, Y in the bottom ; and its transposed
    mpz_mat_init(Z,X->m + Y->m, U->n);
    mpz_mat_init(Z_t, U->n, X->m + Y->m);
    mpz_mat_vertical_join(Z,X,Y);
    
    mpz_mat_transpose(Z_t, Z);
    mpz_mat_kernel_mod_ui(M, Z_t, p); // M = ker(z -> z*(transpose of Z);
    
    mpz_mat_clear(Z_t);
    mpz_mat_clear(Z);
    mpz_mat_clear(X);
    mpz_mat_clear(Y);
    mpz_mat_clear(U_t);
    mpz_mat_clear(W_t);
}


// I wasn't sure if hnf_backend now follows the convention of magma
// This takes M on (under magma conventions) and returns its hnf under
// magma conventions
void hnf_magma_style(mpq_mat_ptr D, mpq_mat_srcptr M)
{
    // Will contain the result before going in D
    mpq_mat tmp;
    mpq_mat_init(tmp,M->m,M->n);
    mpq_mat_set(tmp,M);
    
    // Thats not really a transposition. In fact, coeff (i,j) goes to
    // (n-1-i,n-1-j);
    for (unsigned int i = 0 ;  i < M->m/2 ; i++){
        mpq_mat_submat_swap(tmp,i,0,tmp,M->m-1-i,0,1,M->n);
    }
    for (unsigned int j = 0 ;  j < M->n/2 ; j++){
        mpq_mat_submat_swap(tmp,0,j,tmp,0,M->n-1-j,M->m,1);
    }
    
    mpz_t denom;
    mpq_t denom_inv, denom_rat;
    mpz_mat tmp_int, T;
    mpz_init(denom);
    mpq_init(denom_inv);
    mpq_init(denom_rat);
    mpz_mat_init(tmp_int,M->m,M->n);
    mpz_mat_init(T,0,0);
    
    // Doing HNF on integer-matrix
    mpq_mat_numden(tmp_int,denom,tmp);
    mpq_set_ui(denom_inv,1,1);
    mpq_set_z(denom_rat,denom);
    mpq_div(denom_inv,denom_inv,denom_rat);
    mpz_hnf_backend(tmp_int,T);
    mpq_mat_set_mpz_mat(tmp,tmp_int);
    mpq_mat_mul_mpq(tmp,tmp,denom_inv);


    // Inverting the transformation made before
    for (unsigned int j = 0 ;  j < M->n/2 ; j++){
        mpq_mat_submat_swap(tmp,0,j,tmp,0,M->n-1-j,M->m,1);
    }
    for (unsigned int i = 0 ;  i < M->m/2 ; i++){
        mpq_mat_submat_swap(tmp,i,0,tmp,M->m-1-i,0,1,M->n);
    }
    
    mpz_mat_clear(T);
    mpz_mat_clear(tmp_int);
    mpq_clear(denom_rat);
    mpq_clear(denom_inv);
    mpz_clear(denom);
    mpq_mat_set(D,tmp);
    mpq_mat_clear(tmp);
}

void mpq_hnf_backend_rev(mpq_mat_ptr M0, mpz_mat_ptr T0)
{
    // Will contain the result before going in D
    mpq_mat tmp;
    mpq_mat_init(tmp,M0->m,M0->n);
    mpq_mat_set(tmp,M0);
    
    mpz_t denom;
    mpq_t denom_inv, denom_rat;
    mpz_mat tmp_int;
    mpz_init(denom);
    mpq_init(denom_inv);
    mpq_init(denom_rat);
    mpz_mat_init(tmp_int,M0->m,M0->n);
    
    // Doing HNF on integer-matrix
    mpq_mat_numden(tmp_int,denom,tmp);
    mpq_set_ui(denom_inv,1,1);
    mpq_set_z(denom_rat,denom);
    mpq_div(denom_inv,denom_inv,denom_rat);
    mpz_hnf_backend_rev(tmp_int,T0);
    mpq_mat_set_mpz_mat(tmp,tmp_int);
    mpq_mat_mul_mpq(tmp,tmp,denom_inv);

    mpz_mat_clear(tmp_int);
    mpq_clear(denom_rat);
    mpq_clear(denom_inv);
    mpz_clear(denom);
    mpq_mat_set(M0,tmp);
    mpq_mat_clear(tmp);
}

#if 0
void factorization_of_prime(vector<pair<cxx_mpq_mat, int>>& ideals, mpz_poly_srcptr g, const unsigned long p, gmp_randstate_t state)
{
    int n = g->deg;
    
    mpq_mat G, G_inv; // G = Basis of p-maximal-order Ok ; Id_q = identity (in rational domain)
    mpz_mat Ip; // p radical of Ok, in its basis
    cxx_mpz_mat Id_z; // Identity (in integers domain)
    cxx_mpq_mat Ip_rat; // p radical, in the basis of alpha^
    mpq_mat_init(G,n,n);
    mpq_mat_init(G_inv,n,n);
    mpz_mat_init(Ip,n,n);
    mpz_mat_realloc(Id_z,n,n);
    mpq_mat_realloc(Ip_rat,n,n);
    mpz_mat_set_ui(Id_z,1);
    
    // Computing the generators of p-maximal order, storing them in G
    mpq_mat_set(G, p_maximal_order(g,p));
    // Computing the p-radical of the order of G, storing it in Ip
    p_radical_of_order(Ip,G,g,p);
    mpq_mat_inv(G_inv,G);
    
    // Computing the p-radical in the basis of alpha^
    mpq_mat_set_mpz_mat(Ip_rat,Ip);
    mpq_mat_mul(Ip_rat,Ip_rat,G);
    hnf_magma_style(Ip_rat,Ip_rat);
    
    
    
    // The initial value in pick_from
    subspace_ideal initial;
    initial.E = Id_z; // It is the basis of Ok/p*Ok, but written on the basis of Ok/p*Ok ; thus, it's the identity
    initial.I = Id_z; //p*identity ; the ideal we're trying to separate, in the basis of Ok/p*Ok
    initial.V = Id_z; // identity ; current characteristic subspace, in the basis of Ok/p*Ok
    mpz_mat_mul_ui(initial.I,initial.I,p);
    
    // The vector on which the recursion will happen.
    queue<subspace_ideal> pick_from;
    pick_from.push(initial);
    
    while (!pick_from.empty()) {
        mpz_poly f;
        mpz_poly_init(f,n);
        
        // Getting the head element
        subspace_ideal current = pick_from.front();
        pick_from.pop();
        
        // Picking one random element of subspace E
        cxx_mpz_mat c_E, c_O; // c_E : coeffs of linear combination of dim(E) elements (which are given on the basis of O)
                              // c_O : same element, given directly on the basis of O
        mpz_mat_realloc(c_E,1,current.E->m);
        mpz_mat_realloc(c_O,1,n);
        for (unsigned int i = 0 ; i < current.E->m ; i++) {
            mpz_set_ui(mpz_mat_entry(c_E,0,i),gmp_urandomm_ui(state,p));
        }
        mpz_mat_mul_mod_ui(c_O,c_E,current.E,p);
        
        // Finding its minimal polynomial
        cxx_mpz_mat Mc;
        mpz_mat_realloc(Mc,n,n);
        matrix_of_multiplication_theta_local(Mc,G,c_O,g,p);
        minimal_poly_of_mul_theta(f,G,c_O,g,p);
        mpz_poly_cleandeg(f,n);
        
        // Factorization of the minimal polynomial
        mpz_poly_factor_list lf;
        mpz_poly_factor_list_init(lf);
        mpz_t p_0;
        mpz_init(p_0);
        mpz_set_ui(p_0,p);
        mpz_poly_factor(lf,f,p_0,state);
        
        // Building the list of characteristic subspaces
        vector<cxx_mpz_mat> char_subspaces;
        for (int i = 0 ; i < lf->size ; i++){
            // For each factor f (with multiplicity m), we compute (f(Mc))^m mod p
            cxx_mpz_mat res, ker, char_sub;
            mpz_mat_realloc(res,n,n);
            mpz_poly_eval_mpz_mat_mod_ui(res,Mc,lf->factors[i]->f,p);
            mpz_mat_pow_ui_mod_ui(res,res,lf->factors[i]->m,p);
            
            // Computing the kernel
            mpz_mat_kernel_mod_ui(ker,res,p);
            
            
            //Now we have to compute the intersection of Vect(ker) and V
            intersection_of_subspaces_mod_ui(char_sub, current.V, ker, p);
            mpz_gauss_backend_mod_ui(char_sub,NULL,p);
            char_subspaces.push_back(char_sub);
        }
        
        // Purging it from null characteristic subspaces
        mpz_poly_factor_list fac_Pc; // reduced list of factors
        mpz_poly_factor_list_init(fac_Pc);
        vector<cxx_mpz_mat> r_char_subspaces;
        for (int i = 0 ; i < lf->size ; i++){
            if(char_subspaces[i]->m > 0){
                mpz_poly_factor_list_push(fac_Pc, lf->factors[i]->f, lf->factors[i]->m);
                cxx_mpz_mat aux;
                aux = char_subspaces[i];
                r_char_subspaces.push_back(aux);
            }
        }        
        
        // Now finishing the run on the tree
        for (int i = 0 ; i < fac_Pc->size ; i++){
            
            int e = fac_Pc->factors[i]->m;
                
            /* Consider the elements of O which map to zero in O/I^e, and
             * non-zero elsewhere. These generate I !  Those are going in
             * gens*/
            cxx_mpq_mat gens;
            mpq_mat_realloc(gens,n,n);
            
            int current_line = 0;
            for (int j = 0 ; j < fac_Pc->size ; j++){
                unsigned long p1;
                if(i == j)
                    p1 = p;
                else
                    p1 = 1;
                for (unsigned int k = 0 ; k < r_char_subspaces[j]->m ; k++){
                    // Extracting the k-th vector of the basis of j-th
                    // characteristic subspace
                    cxx_mpz_mat aux;
                    cxx_mpz_mat v;
                    cxx_mpq_mat v_rat;
                    mpz_mat_realloc(aux,1,n);
                    mpz_mat_realloc(v,1,n);
                    
                    mpz_mat_submat_swap(aux,0,0,r_char_subspaces[j],k,0,1,n);
                    mpz_mat_set(v,aux);
                    mpz_mat_submat_swap(aux,0,0,r_char_subspaces[j],k,0,1,n);
                    // v now contains this vector, on the basis of O
                    
                    
                    // Multiplying v by G transfers v in the number field
                    // K (rational coefficients, thus v goes into v_rat)
                    mpq_mat_set_mpz_mat(v_rat,v);
                    mpq_mat_mul(v_rat,v_rat,G);
                    mpq_mat_mul_ui(v_rat,v_rat,p1);

                    mpq_mat_submat_swap(gens,current_line,0,v_rat,0,0,1,n);
                    current_line++;
                }
            }
             
             // Transfering current ideal (the one we try to separate)
             // subspace basis of alpha^
             cxx_mpq_mat gcd_with;
             mpq_mat_set_mpz_mat(gcd_with,current.I);
             mpq_mat_mul(gcd_with,gcd_with,G);
             
             //Now computing 3n*n matrix with gens, Ip_rat and gcd_with
             
             cxx_mpq_mat big_matrix, bigger_matrix, Ix, new_ideal_rat;
             cxx_mpz_mat new_ideal;
             mpq_mat_realloc(Ix,n,n);
             mpq_mat_realloc(new_ideal_rat,n,n);
             mpz_mat_realloc(new_ideal,n,n);
             
             mpq_mat_vertical_join(big_matrix,gens,Ip_rat);
             mpq_mat_vertical_join(bigger_matrix,big_matrix,gcd_with);
             hnf_magma_style(bigger_matrix, bigger_matrix);
             mpq_mat_submat_swap(Ix,0,0,bigger_matrix,2*n,0,n,n);
             
             mpq_mat_mul(new_ideal_rat,Ix,G_inv);
             mpq_mat_numden(new_ideal,NULL,new_ideal_rat);
             
             if( (int) r_char_subspaces[i]->m == e*(fac_Pc->factors[i]->f->deg)){
                 pair<cxx_mpq_mat, int> new_pair = make_pair(Ix,e);
                 ideals.push_back(new_pair);
             }
             else{
                 subspace_ideal new_subspace_ideal;
                 new_subspace_ideal.E = r_char_subspaces[i];
                 new_subspace_ideal.I = new_ideal;
                 new_subspace_ideal.V = r_char_subspaces[i];
                 pick_from.push(new_subspace_ideal);
             }
             
            
        }
        
        mpz_clear(p_0);
        mpz_poly_factor_list_clear(fac_Pc);
        mpz_poly_factor_list_clear(lf);
        mpz_poly_clear(f);
    }
    
    mpz_mat_clear(Ip);
    mpq_mat_clear(G_inv);
    mpq_mat_clear(G);
}
#endif

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

    hnf_magma_style(I,I);
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
    
    hnf_magma_style(I,I);
    
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
    mpz_hnf_backend(ker,C); // We don't need C anymore, so we use it to contain the transformation matrix
    
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

