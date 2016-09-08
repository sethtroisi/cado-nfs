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

using namespace std;

/*{{{ timer*/
double
seconds (void)
{
    struct rusage res[1];
    getrusage(RUSAGE_SELF, res);
    uint64_t r;
    r = (uint64_t) res->ru_utime.tv_sec;
    r *= (uint64_t) 1000000UL;
    r += (uint64_t) res->ru_utime.tv_usec;
    return r * 1.0e-6;
}
/*}}}*/
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

void mpq_mat_row_to_poly(mpz_poly_ptr f, mpz_ptr denom, mpq_mat_srcptr M, const unsigned int i)
{
    mpz_poly_realloc(f,M->n-1);
    mpz_set_si(denom,1);
    unsigned int j;
    for (j = 0 ; j < M->n ; j++) {
        mpz_lcm(denom,denom,mpq_denref(mpq_mat_entry_const(M,i,j)));
    }
    for (j = 0 ; j < M->n; j++){
        mpq_t aux;
        mpz_t num;
        mpq_init(aux);
        mpz_init(num);
        mpq_set_z(aux,denom);
        mpq_mul(aux,aux,mpq_mat_entry_const(M,i,j));
        mpz_set(num,mpq_numref(aux));
        mpz_poly_setcoeff(f,j,num);
        mpq_clear(aux);
        mpz_clear(num);
    }
}

void mpq_poly_to_mat_row(mpq_mat_ptr M, const unsigned int i, mpz_poly_srcptr f, mpz_srcptr denom)
{
    ASSERT_ALWAYS(f->deg < (int) M->n);
    mpz_t coeff;
    mpz_init(coeff);
    for (int j = 0 ; j <= f->deg; j++){
        mpz_poly_getcoeff(coeff,j,f);
        mpq_set_num(mpq_mat_entry(M,i,j),coeff);
        mpq_set_den(mpq_mat_entry(M,i,j),denom);
        mpq_canonicalize(mpq_mat_entry(M,i,j));
    }
    mpz_clear(coeff);
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


/* Structure defining a rational polynomial
 * For now it's just here to make it easier to store rational polynomials
 * in a vector Some parts of ant.cpp could be changed with this, but
 * that's no priority for now */
struct cxx_mpq_poly
{
    cxx_mpz_poly num;
    cxx_mpz den;
};

// Prints the polynomial
void print_polynomial(mpz_t * f, int degree)
{
    int i;
    for (i = degree; i >= 0; i--) {
	if (mpz_get_ui(f[i]) != 0) {
	    gmp_printf("%Zd", f[i]);
	    if (i > 0) {
		printf("*x^%d + ", i);
	    }
	}
    }
    printf("\n");
}

// Takes a matrix B containing the generators (w_0, ... w_{n-1}) of an order,
// and returns the matrix U containing ((w_0)^p, ..., (w_{n-1})^p),
// reduced mod g and mod p.
void generators_to_power_p(mpq_mat_ptr U, mpq_mat_srcptr B,
			   mpz_poly_srcptr g, const unsigned int p)
{
    ASSERT_ALWAYS(B->m == B->n);

    mpq_mat_realloc(U, B->m, B->n);

    mpz_t lcm, tmp;
    mpq_t aux1, K_rat;
    mpz_init(lcm);
    mpz_init(tmp);
    mpq_init(aux1);
    mpq_init(K_rat);

    for (unsigned int i = 0; i < U->m; i++) {
	mpz_poly f;
	mpz_poly_init(f, B->n - 1);

	// Putting the LCM of all denominators of coefficients of w[i] in lcm
	mpz_set_si(lcm, 1);
	for (unsigned int j = 0; j < B->n; j++) {
	    mpz_lcm(lcm, lcm, mpq_denref(mpq_mat_entry_const(B, i, j)));
	}

	// Generating the polynomial
	for (unsigned int j = 0; j < B->n; j++) {
	    mpq_set_z(K_rat, lcm);
	    //gmp_printf("%Zd/%Zd * %Zd/%Zd\n", mpq_numref(mpq_mat_entry_const(B,i,j)), mpq_denref(mpq_mat_entry_const(B,i,j)), mpq_numref(K_rat), mpq_denref(K_rat));
	    mpq_mul(aux1, mpq_mat_entry_const(B, i, j), K_rat);
	    ASSERT_ALWAYS(!mpz_cmp_ui(mpq_denref(aux1), 0) == 0);
	    mpz_poly_setcoeff(f, j, mpq_numref(aux1));
	}

	// Computing f^p mod g (it returns (lcm * the corresponding generator)^p
	mpz_poly_power_mod_f(f, f, g, p);

	mpz_pow_ui(lcm, lcm, p);

	// Storing w[i] in i-th row of U
	for (int j = 0; j <= f->deg; j++) {
	    mpz_poly_getcoeff(tmp, j, f);
	    mpq_set_num(mpq_mat_entry(U, i, j), tmp);
	    mpq_set_den(mpq_mat_entry(U, i, j), lcm);
	    mpq_canonicalize(mpq_mat_entry(U, i, j));
	}
	mpz_poly_clear(f);
    }
    mpq_clear(K_rat);
    mpq_clear(aux1);
    mpz_clear(tmp);
    mpz_clear(lcm);
}

/* This is almost like hnf_backend, except that we do it in a different
 * order. We consider columns in reverse order. And we also reverse the
 * order of the first non-zero rows we obtain.
 */
int mpz_hnf_backend_rev(mpz_mat_ptr M0, mpz_mat_ptr T0)
{
    mpz_mat M, T;
    mpz_mat_init(M, M0->m, M0->n);
    mpz_mat_realloc(T0, M0->m, M0->m);
    mpz_mat_init(T, T0->m, T0->n);
    for(unsigned int i = 0; i < M0->m; i++) {
        for(unsigned int j = 0; j < M0->n; j++) {
            mpz_swap(mpz_mat_entry(M,i,j),mpz_mat_entry(M0,i,M0->n-1-j));
        }
    }
    int s = mpz_hnf_backend(M, T);
    unsigned int rank = 0;
    for(unsigned int i = 0; i < M->m; i++) {
        for(unsigned int j = 0; j < M->n; j++) {
            if (mpz_cmp_ui(mpz_mat_entry_const(M,i,j),0) != 0) {
                rank=i+1;
                break;
            }
        }
    }
    ASSERT_ALWAYS(M->m == T->m);
    ASSERT_ALWAYS(M->m == T->n);
    ASSERT_ALWAYS(M->m == T0->m);
    ASSERT_ALWAYS(M->m == T0->n);
    for(unsigned int i = 0; i < rank; i++) {
        for(unsigned int j = 0; j < M->n; j++) {
            mpz_swap(mpz_mat_entry(M0,i,j),mpz_mat_entry(M,rank-1-i,M0->n-1-j));
        }
        for(unsigned int j = 0; j < T->n; j++) {
            mpz_swap(mpz_mat_entry(T0,i,j),mpz_mat_entry(T,rank-1-i,j));
        }
    }

    for(unsigned int i = rank; i < M->m; i++) {
        for(unsigned int j = 0; j < M->n; j++) {
            mpz_swap(mpz_mat_entry(M0,i,j),mpz_mat_entry(M,i,M0->n-1-j));
        }
        for(unsigned int j = 0; j < T->n; j++) {
            mpz_swap(mpz_mat_entry(T0,i,j),mpz_mat_entry(T,i,j));
        }
    }
    if ((rank/2)&1) s = -s;
    mpz_mat_clear(M);
    mpz_mat_clear(T);
    return s;
}



// Builds the block matrix containing p*identity in the top, and K in the bottom
// Then computes its HNF and stores it in I
void join_HNF(mpz_mat_ptr I, mpz_mat_srcptr K, const unsigned int p)
{
    mpz_mat J, T0;
    mpz_mat_init(J, K->n, K->n);
    mpz_mat_init(T0, K->n, K->n);
    mpz_mat_realloc(I, K->n, K->n);

    mpz_mat_set_ui(J, 1);
    mpz_mat_multiply_by_ui(J, J, p);
    mpz_mat_vertical_join(I, J, K);
    mpz_hnf_backend_rev(I, T0);
    mpz_mat_submat_swap(I, 0, 0, J, 0, 0, K->n, K->n);
    mpz_mat_set(I, J);

    mpz_mat_clear(T0);
    mpz_mat_clear(J);
}

// Builds the (n,n^2) matrix containing the integers mod p
// They're obtained with the generators of O and the generators of I_p
// Generators of O are in B (they're given in the basis of alpha^)
// Generators of I_p are in I (in the basis of O)
// The products are computed mod g
void generators_to_integers_mod_p(mpz_mat_ptr M, mpq_mat_srcptr B,
				  mpz_mat_srcptr I_p, mpz_poly_srcptr g,
				  const unsigned int p)
{
    ASSERT_ALWAYS((B->m == B->n) && (I_p->m == I_p->n) && (B->m == I_p->m));
    unsigned int n = B->m;
    mpz_mat I;
    mpq_mat I_rat, I_inv, B_inv;
    mpz_mat_realloc(M, n, n * n);

    mpz_mat_init(I, n, n);
    mpq_mat_init(I_rat, n, n);
    mpq_mat_init(I_inv, n, n);
    mpq_mat_init(B_inv, n, n);

    mpz_mat_set(I, I_p);
    mpz_mat_to_mpq_mat(I_rat, I);
    mpq_mat_invert(I_inv, I_rat);
    mpq_mat_invert(B_inv, B);

    //printf("Ip is :\n"); mpz_mat_fprint(stdout, I); printf("\n");
    //printf("Ip^-1 is :\n"); mpq_mat_fprint(stdout, I_inv); printf("\n");
        
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {

            // Gammma : polynomial of generator of Ip (only the
            // numerators, multiplied by LCM of denominators) c : same
            // for generators of O
            mpz_poly gamma, c, aux;
            // LCM of denominators of generators of Ip, of O, and one
            // auxiliary coefficient
            mpz_t denom_g, denom_c, coeff;
            // Line matrices containing the product gamma*c mod g
            mpz_mat row, aux_row;
            // Auxiliary row matrices, necessary because generators of Ip
            // need to be converted in the basis of alpha^
            mpq_mat row_q, c_mat, res;

            // Initialisation
            mpz_mat_init(row, 1, n);
            mpz_mat_init(aux_row, 1, n);
            mpq_mat_init(row_q, 1, n);
            mpq_mat_init(c_mat, 1, n);
            mpq_mat_init(res, n, n);
            mpz_poly_init(gamma, n);
            mpz_poly_init(c, n);
            mpz_poly_init(aux, n);
            mpz_init(denom_g);
            mpz_init(denom_c);
            mpz_init(coeff);

            // Storing one generator of B in the polynomial gamma and in
            // denom_g
            mpq_mat_row_to_poly(gamma, denom_g, B, i);
            //Storing one generator of I_p in the polyomial c and in
            //denom_c
            mpz_mat_submat_swap(aux_row, 0, 0, I, j, 0, 1, n);
            mpz_mat_set(row, aux_row);
            mpz_mat_submat_swap(aux_row, 0, 0, I, j, 0, 1, n);
            mpz_mat_to_mpq_mat(row_q, row);
            mpq_mat_multiply(c_mat, row_q, B);
            mpq_mat_row_to_poly(c, denom_c, c_mat, 0);

            // Computing gamma*c mod g
            mpz_poly_mul_mod_f(aux, gamma, c, g);
            
            /*
            if((i == 1) && (j == n-1)){
                printf("row_q is : "); mpq_mat_fprint(stdout,row_q); printf("\n");
            }
            */
            
            // Storing the result in row_q
            for (int k = 0; k <= aux->deg; k++) {
                mpz_poly_getcoeff(coeff, k, aux);
                mpq_set_num(mpq_mat_entry(row_q, 0, k), coeff);
                mpz_mul(coeff, denom_c, denom_g);
                mpq_set_den(mpq_mat_entry(row_q, 0, k), coeff);
                mpq_canonicalize(mpq_mat_entry(row_q, 0, k));
            }
            // Since we re-use row_q from before, we have to delete all
            // previous coefficients
            for (int k = aux->deg+1 ; k < (int) n ; k++){
                mpz_set_ui(mpq_numref(mpq_mat_entry(row_q, 0, k)), 0);
                mpz_set_ui(mpq_denref(mpq_mat_entry(row_q, 0, k)), 1);
            }
            

            // Converting row_q (gamma*c mod g) in the basis of I_p (it
            // is supposed to contain only integers)
            mpq_mat_multiply(res, row_q, B_inv);
            mpq_mat_multiply(res, res, I_inv);

            // Extracting the numerators (remember, integers only) into a
            // mpz_mat
            mpq_mat_numden(row, coeff, res);
            
            /*
            if((i == 1) && (j == n-1)){
                printf("gamma is : "); mpz_poly_fprintf(stdout,gamma);
                printf("c is : "); mpz_poly_fprintf(stdout,c);
                printf("aux is : "); mpz_poly_fprintf(stdout,aux);
                printf("row_q is : "); mpq_mat_fprint(stdout,row_q); printf("\n");
                printf("row is : "); mpz_mat_fprint(stdout,row); printf("\n");
            }
            */
            
            // Computing the same matrix, modulo p (it is supposed to be
            // a vector of n integers, associated to gamma)
            mpz_mat_mod_ui(row, row, p);
            mpz_mat_submat_swap(row, 0, 0, M, i, n * j, 1, n);

            mpq_mat_clear(res);
            mpq_mat_clear(c_mat);
            mpq_mat_clear(row_q);
            mpz_mat_clear(row);
            mpz_mat_clear(aux_row);
            mpz_clear(coeff);
            mpz_clear(denom_g);
            mpz_clear(denom_c);
            mpz_poly_clear(gamma);
            mpz_poly_clear(c);
            mpz_poly_clear(aux);
        }
    }
    mpz_mat_clear(I);
    mpq_mat_clear(I_rat);
    mpq_mat_clear(I_inv);
    mpq_mat_clear(B_inv);
}

void read_data(unsigned int *deg, mpz_poly_ptr f, /*mpz_mat_ptr gen,*/
	       FILE * problemfile)
{
    int rc = fscanf(problemfile, "%u", deg);
    ASSERT_ALWAYS(rc == 1);
    
    mpz_poly_realloc(f, *deg + 1);
    mpz_t c;
    mpz_init(c);
    for (unsigned int i = 0; i <= *deg; i++) {
        gmp_fscanf(problemfile, "%Zd", &c);
        mpz_poly_setcoeff (f, i, c);
    }
    mpz_clear(c);

    // Creating one ideal from two coefficients a and b, representing the generator a-b*alpha
    //mpz_mat_realloc(gen, 2, *deg);
    /*
    for(unsigned int i = 0 ; i < gen->m ; i++){
        for(unsigned int j = 0 ; j < gen->n ; j++){
            int a;
            fscanf(problemfile, "%d", &a);
            mpz_set_si(mpz_mat_entry(gen,i,j),a);
        }
    }
    */
    /*
    mpq_mat_realloc(J,*deg,*deg);
    for (unsigned int i = 0; i < *deg; i++) {
        long int denom;
        denom = 1;
        fscanf(problemfile, "%ld", &denom);
        for (unsigned int j = 0; j < *deg; j++) {
            long int c;
            fscanf(problemfile, "%ld", &c);
            mpq_set_si(mpq_mat_entry(J, i, j), c, denom);
            mpq_canonicalize(mpq_mat_entry(J,i,j));
        }
    }*/
}

// Stores in I the p-radical of the order generated by B (in the basis of
// these generators)
void p_radical_of_order(mpz_mat_ptr I, mpq_mat_srcptr B, mpz_poly_srcptr g, const unsigned int p)
{
    ASSERT_ALWAYS(B->m == B->n);
    unsigned int n = B->m;
    
    mpq_mat B_inv, U, T;
    mpz_mat X, K;
    mpz_t den;
    mpq_mat_init(B_inv,n,n);
    mpq_mat_init(U,n,n);
    mpq_mat_init(T,n,n);
    mpz_mat_init(X,n,n);
    mpz_mat_init(K,n,n);
    mpz_init(den);
    
    
    // Inverting B
    mpq_mat_invert(B_inv, B);

    // Now building the matrix U, containing all generators to the power
    // of p Generators are polynomials, stored in the matrix B
    generators_to_power_p(U, B, g, p);

    // Now do U * B^-1, storing in X the application  F : z -> (z^p mod
    // g mod p)
    mpq_mat_multiply(T, U, B_inv);
    mpq_mat_numden(X, den, T);
    mpz_mat_mod_ui(X, X, p);

    /* Which power of p ? */
    int k = 1;
    for (unsigned int pk = p; pk < n; k++, pk *= p);
    mpz_mat_power_ui_mod_ui(X, X, k, p);

    // Storing in K a basis of Ker((z -> (z^%p mod g mod p))^k)
    mpz_mat_kernel(K, X, p);

    // Getting generators of the p radical from Ker(X) by computing HNF
    // of the vertical block matrix (p*Id, K);
    join_HNF(I, K, p);
    
    mpz_clear(den);
    mpz_mat_clear(K);
    mpz_mat_clear(X);
    mpq_mat_clear(T);
    mpq_mat_clear(U);
    mpq_mat_clear(B_inv);
}

// Builds one p-maximal-order starting from the matrix B, containing the
// generators of one order (in the basis of alpha^) In the number field
// of the polynomial f, of degree n Stores the generators of this
// p-maximal-order (in the basis of alpha^) in D
void p_maximal_order(mpq_mat_ptr D, mpz_poly_srcptr f,
		    const unsigned int p)
{

    mpq_mat new_D, B;
    mpz_poly g;
    mpq_t p_inv;
    unsigned int n = f->deg;

    mpq_mat_init(B,n,n);
    mpq_mat_init(new_D, n, n);
    mpz_poly_init(g, n);	// Monic polynomial corresponding to g
    mpq_init(p_inv);		// 1/p



    for (unsigned int i = 0; i < n; i++) {
        int denom;
        denom = 1;
        for (unsigned int j = 0; j < n; j++) {
            int c;
            if(i == j){c = 1;}
            else{c = 0;}
            mpq_set_si(mpq_mat_entry(B, i, j), c, denom);
            mpq_canonicalize(mpq_mat_entry(B,i,j));
        }
    }
    
    
    // Storing in g the monic polynomial such as (fd*alpha) is a root of
    // if and only if alpha is a root of f
    mpz_poly_to_monic(g, f);
    // Initializing D to B
    mpq_mat_set(D, B);
    mpq_mat_set(new_D, B);
    mpq_set_ui(p_inv, 1, p);

    do {
        mpq_mat_set(D, new_D);

        mpq_mat T, U, J2;
        mpz_mat X, K, I, T0, M, K_M, J;
        mpz_t den;

        mpz_mat_init(K_M, n, n);	// Ker(M)
        mpz_mat_init(X, n, n);	// Matrix of the Froebenius morphism
        mpz_mat_init(T0, n, n);
        mpz_mat_init(K, n, n);	// Ker(X)
        mpz_mat_init(M, n, n * n);	// The (n,n^2) Matrix containing all products of generators of O and of its p-radical
        mpz_mat_init(I, n, n);	// Matrix of the p-radical
        mpq_mat_init(T, n, n);
        mpq_mat_init(U, n, n);	// Contains all the generators of D, to the power of p
        mpq_mat_set_ui(T, 1);
        mpz_mat_init(J, n, n);	// Generators of p * new order in the basis of the previous order (mpz_mat version)
        mpq_mat_init(J2, n, n);	// Generators of p * new order in the basis of the previous order (mpq_mat version)
        mpz_init(den);

        // Getting the p-radical of the order generated by D in I.
        p_radical_of_order(I,D,g,p);

        // Building the (n,n^2) matrix containing the integers mod p associated to all genereators of O
        generators_to_integers_mod_p(M, D, I, g, p);
        //printf("M is :\n"); mpz_mat_fprint(stdout, M); printf("\n");
        
        // Computing Ker(M)
        mpz_mat_kernel(K_M, M, p);
        //printf("Ker(M) is :\n"); mpz_mat_fprint(stdout, K_M); printf("\n");

        // Getting generators of p*O' by computing HNF of the vertical block matrix (p*Id, K_M);
        join_HNF(J, K_M, p);
        // Converting into a mpq_mat
        mpz_mat_to_mpq_mat(J2, J);
        // Converting in the basis of alpha^
        mpq_mat_multiply(new_D, J2, D);
        // Dividing by p
        mpq_mat_multiply_by_mpq(new_D, new_D, p_inv);


        mpq_mat_clear(J2);
        mpz_mat_clear(J);
        mpz_mat_clear(K_M);
        mpz_mat_clear(M);
        mpz_mat_clear(T0);
        mpz_mat_clear(K);
        mpz_mat_clear(X);
        mpz_mat_clear(I);
        mpq_mat_clear(T);
        mpq_mat_clear(U);
        mpz_clear(den);

    } while (!mpq_mat_eq(new_D, D));


    // Now we have to make the HNF of D
    // First take the lcm of all denominators in D
    mpz_t denom; // lcm of denominator
    mpq_t denom_q, denom_inv; // denom_q = denom ; denom_inv = 1/denom
    mpz_init(denom);
    mpq_init(denom_q);
    mpq_init(denom_inv);
    mpz_set_ui(denom,1);
    mpq_set_ui(denom_inv,1,1);
    for(unsigned int i = 0 ; i < n ; i++){
        for(unsigned int j = 0 ; j < n ; j++){
            mpz_lcm(denom,denom,mpq_denref(mpq_mat_entry(D,i,j)));
        }
    }
    mpq_set_z(denom_q,denom);
    mpq_div(denom_inv,denom_inv,denom_q);
    
    // Multiplying D by the denom
    mpq_mat_multiply_by_mpq(D,D,denom_q);
    
    // Preparing the HNF in D_z
    mpz_t denom_eq_1;
    mpz_mat D_z,T;
    mpz_mat_init(D_z,n,n);
    mpz_mat_init(T,n,n);
    mpz_init(denom_eq_1);
    
    mpq_mat_numden(D_z,denom_eq_1,D);
    mpz_hnf_backend_rev(D_z,T);
    mpz_mat_to_mpq_mat(D,D_z);
    mpq_mat_multiply_by_mpq(D,D,denom_inv);
    
    mpq_clear(denom_inv);
    mpz_clear(denom_eq_1);
    mpz_mat_clear(D_z);
    mpz_mat_clear(T);
    mpq_clear(denom_q);
    mpz_clear(denom);
    mpq_mat_clear(new_D);
    mpq_mat_clear(B);
    mpz_poly_clear(g);
    mpq_clear(p_inv);
}

/* theta is expected to be an element of the order whose basis is given
 * by W. Thus we read theta as a vector with integer coordinates (matrix
 * of size 1*n)
 */
void matrix_of_multiplication_by_theta_local(mpz_mat_ptr M, mpq_mat_srcptr W, mpz_mat_srcptr theta, mpz_poly_srcptr g, const unsigned int p)
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
    mpz_mat_to_mpq_mat(theta_rat,theta);
    mpq_mat_multiply(theta_rat,theta_rat,W); // Now contains theta in the basis of alpha^ (rational coefficients)

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
    mpq_mat_invert(W_inv,W);
    mpq_mat_multiply(times_theta_rat,times_theta_rat,W_inv);
    
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
void minimal_poly_of_mul_by_theta(mpz_poly_ptr f, mpq_mat_srcptr W, mpz_mat_srcptr theta, mpz_poly_srcptr g, const unsigned int p)
{
    unsigned int n = W->m;
    
    mpz_mat times_theta;
    mpz_mat_init(times_theta,n,n);

    matrix_of_multiplication_by_theta_local(times_theta, W, theta, g, p);
    
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
        mpz_mat_multiply(current,current,times_theta);
        mpz_mat_mod_ui(current,current,p);
    }


    // Now computing its kernel
    mpz_mat K;
    mpz_mat_init(K,0,0);
    mpz_mat_kernel(K,M,p);

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
void intersection_of_subspaces_mod_ui(mpz_mat_ptr M, mpz_mat_srcptr U, mpz_mat_srcptr W, const unsigned int p)
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
    
    mpz_mat_kernel(X,U_t,p); // X = ker(x -> x*(transpose of U)
    mpz_mat_kernel(Y,W_t,p); // Y = ker(y -> y*(transpose of W)
    
    mpz_mat Z, Z_t; // Block matrix, X on the top, Y in the bottom ; and its transposed
    mpz_mat_init(Z,X->m + Y->m, U->n);
    mpz_mat_init(Z_t, U->n, X->m + Y->m);
    mpz_mat_vertical_join(Z,X,Y);
    
    mpz_mat_transpose(Z_t, Z);
    mpz_mat_kernel(M, Z_t, p); // M = ker(z -> z*(transpose of Z);
    
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
    mpz_mat_to_mpq_mat(tmp,tmp_int);
    mpq_mat_multiply_by_mpq(tmp,tmp,denom_inv);


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
    mpz_mat_to_mpq_mat(tmp,tmp_int);
    mpq_mat_multiply_by_mpq(tmp,tmp,denom_inv);

    mpz_mat_clear(tmp_int);
    mpq_clear(denom_rat);
    mpq_clear(denom_inv);
    mpz_clear(denom);
    mpq_mat_set(M0,tmp);
    mpq_mat_clear(tmp);
}

void factorization_of_prime(vector<pair<cxx_mpq_mat, int>>& ideals, mpz_poly_srcptr g, const unsigned int p, gmp_randstate_t state)
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
    p_maximal_order(G,g,p);
    // Computing the p-radical of the order of G, storing it in Ip
    p_radical_of_order(Ip,G,g,p);
    mpq_mat_invert(G_inv,G);
    
    // Computing the p-radical in the basis of alpha^
    mpz_mat_to_mpq_mat(Ip_rat,Ip);
    mpq_mat_multiply(Ip_rat,Ip_rat,G);
    hnf_magma_style(Ip_rat,Ip_rat);
    
    
    
    // The initial value in pick_from
    subspace_ideal initial;
    initial.E = Id_z; // It is the basis of Ok/p*Ok, but written on the basis of Ok/p*Ok ; thus, it's the identity
    initial.I = Id_z; //p*identity ; the ideal we're trying to separate, in the basis of Ok/p*Ok
    initial.V = Id_z; // identity ; current characteristic subspace, in the basis of Ok/p*Ok
    mpz_mat_multiply_by_ui(initial.I,initial.I,p);
    
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
        mpz_mat_multiply_mod_ui(c_O,c_E,current.E,p);
        
        // Finding its minimal polynomial
        cxx_mpz_mat Mc;
        mpz_mat_realloc(Mc,n,n);
        matrix_of_multiplication_by_theta_local(Mc,G,c_O,g,p);
        minimal_poly_of_mul_by_theta(f,G,c_O,g,p);
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
            mpz_mat_power_ui_mod_ui(res,res,lf->factors[i]->m,p);
            
            // Computing the kernel
            mpz_mat_kernel(ker,res,p);
            
            
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
                unsigned int p1;
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
                    mpz_mat_to_mpq_mat(v_rat,v);
                    mpq_mat_multiply(v_rat,v_rat,G);
                    mpq_mat_multiply_by_ui(v_rat,v_rat,p1);

                    mpq_mat_submat_swap(gens,current_line,0,v_rat,0,0,1,n);
                    current_line++;
                }
            }
             
             // Transfering current ideal (the one we try to separate)
             // subspace basis of alpha^
             cxx_mpq_mat gcd_with;
             mpz_mat_to_mpq_mat(gcd_with,current.I);
             mpq_mat_multiply(gcd_with,gcd_with,G);
             
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
             
             mpq_mat_multiply(new_ideal_rat,Ix,G_inv);
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
        mpz_poly_srcptr g, mpq_mat_srcptr I, const unsigned int p)
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
            mpq_mat_invert(G_inv,G);
            mpq_mat_multiply(aux_mat, aux_mat, G_inv);
            
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
    mpz_mat_kernel(ker, C, p);
    mpz_hnf_backend(ker,C); // We don't need C anymore, so we use it to contain the transformation matrix
    
    cxx_mpq_mat ker_rat;
    mpq_t p_inv;
    mpq_init(p_inv);
    mpq_set_ui(p_inv,1,p);
    mpq_mat_realloc(ker_rat, 1, n);
    
    // Turning the vectors of the kernel in the basis of alpha^, and dividing by p
    mpz_mat_to_mpq_mat(ker_rat, ker);
    mpq_mat_multiply(ker_rat, ker_rat, G);
    mpq_mat_multiply_by_mpq(ker_rat, ker_rat, p_inv);
        
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
int valuation_of_ideal_at_prime_ideal(mpq_mat_srcptr G, mpz_poly_srcptr g, mpq_mat_srcptr J, mpq_mat_srcptr I, const int e, const unsigned int p)
{
    ASSERT_ALWAYS((G->m == G->n) && (G->n == J->m) && (J->m == J->n) && (J->n == I->m) && (I->m == I->n));
    unsigned int n = G->m;
    
    cxx_mpq_mat G_inv, MJ, h;
    cxx_mpz_mat MJ_int;
    mpz_t dJ, denom;
    mpz_init(dJ); // Initial denominator of MJ
    mpz_init(denom); // Denominator of MJ for the loop
    mpq_mat_realloc(G_inv,n,n);
    mpq_mat_invert(G_inv,G);
    mpq_mat_realloc(MJ,n,n); // MJ is the MJ existing in magma
    mpz_mat_realloc(MJ_int,n,n); // The MJ_int here is only MJ*dJ
    
    mpq_mat_multiply(MJ, J, G_inv);
    mpq_mat_numden(MJ_int, dJ, MJ); // Now MJ_int = dJ*MJ, like in magma
    mpz_mat_to_mpq_mat(MJ, MJ_int); // We need a rational matrix because we have to compute its denominator

    
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
        mpq_mat_multiply(aux, MJ, G); // This way, you obtain all rows of the form (Vector(K,j),Vector(G))
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
        mpq_mat_multiply(MJ, h_times_MJ_G, G_inv);
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
void filter_roots(vector<pair<cxx_mpz_t, cxx_mpz_t>>& roots, unsigned int p)
{
    unsigned int i = 0;

    // Filtering the root
    while(i < roots.size()-1){
        // We look for roots[j] equal to roots[i], for j > i
        unsigned int j = i+1;
        while (j < roots.size()){
            cxx_mpz_t a = roots[j].first;
            cxx_mpz_t b = roots[j].second;
            
            // At this point, we will multiply a and b with 2, then a and
            // b with 3, etc...  if we find one k such that a*k =
            // roots[i].first and b*k = roots[i].second, then it's the
            // same quotient Thus we can erase a/b
            cxx_mpz_t a1, b1;
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
            cxx_mpz_t a1, b1;
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
    
    mpq_mat_invert(G_inv, G);
    mpq_mat_multiply(tmp,I, G_inv);
    
    mpq_mat_numden(coords, NULL, tmp);
    
    mpz_t d;
    mpz_init(d);
    mpz_mat_determinant_triangular(d,coords);
    int a = mpz_get_ui(d);
    int p = mpz_get_ui(mpq_numref(mpq_mat_entry_const(I,0,0)));
    mpz_clear(d);
    
    
    return (int) (log((double) a)/log((double) p));
}

void print_comments_for_badideals_above_p(string& SBAD, string& SBADINFO, const unsigned int side, 
                mpq_mat_srcptr order, mpz_poly_srcptr f, 
                vector<pair<cxx_mpq_mat, int>> ideals, const unsigned int p)
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
    vector<pair<cxx_mpz_t, cxx_mpz_t>> rootsp;
    for (unsigned int i = 0 ; i < p ; i++){
        for (unsigned int j = 1 ; j < p ; j++){
            cxx_mpz_t res, a, b;
            mpz_set_ui(a,i);
            mpz_set_ui(b,j);
            
            mpz_poly_homogeneous_eval_siui(res,f,i,j);
            mpz_mod_ui(res,res,p);
            if(mpz_congruent_ui_p(res,0,p)){
                pair<cxx_mpz_t,cxx_mpz_t> new_elem;
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
    cxx_mpz_t fd;
    mpz_poly_getcoeff(fd,f->deg,f);
    if(mpz_congruent_ui_p(fd,0,p)){
        cxx_mpz_t a,b;
        mpz_set_ui(a,1);
        mpz_set_ui(b,0);
        pair<cxx_mpz_t,cxx_mpz_t> new_elem;
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
        cxx_mpz_t u = rootsp[i].first;
        cxx_mpz_t v = rootsp[i].second;
        
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
        pair<cxx_mpz_t, cxx_mpz_t> Q;
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

int main(int argc, char *argv[])
{				/*{{{ */
    /*
       unsigned int m = 8;
       unsigned int n = 5;
       if (argc == 3) {
       m = strtoul(argv[1], NULL, 0);
       n = strtoul(argv[2], NULL, 0);
       }
       gmp_randstate_t state;
       gmp_randinit_default(state);

       mpq_mat M;
       mpq_mat T;
       mpz_mat Mz;
       mpz_mat Tz;
       mpz_t p;

       mpq_mat_init(M, m, n);
       mpq_mat_init(T, m, m);

       mpz_mat_init(Mz, m, n);
       mpz_mat_init(Tz, m, m);

       mpz_init(p);

       mpz_set_ui(p, 19);

       if (0) {
       printf("\n\nCas 0.1\n\n");
       mpq_mat_urandomm(M, state, p);
       mpq_mat_fprint(stdout, M);
       printf("\n");
       mpq_gauss_backend(M, T);
       mpq_mat_fprint(stdout, M);
       printf("\n");
       mpq_mat_fprint(stdout, T);
       printf("\n");
       }

       if (0) {
       printf("\n\nCas 0.2\n\n");
       mpz_mat_urandomm(Mz, state, p);
       mpz_mat_fprint(stdout, Mz);
       printf("\n");
       mpz_gauss_backend_mod(Mz, Tz, p);
       mpz_mat_fprint(stdout, Mz);
       printf("\n");
       mpz_mat_fprint(stdout, Tz);
       printf("\n");
       }

       if (1) {
       printf("\n\nCas 1\n\n");
       mpz_mat_realloc(Mz, m, n);
       mpz_mat_urandomm(Mz, state, p);
       mpz_mat_fprint(stdout, Mz); printf("\n");
       double t = seconds();
       mpz_hnf_backend(Mz, Tz);
       t = seconds()-t;
       mpz_mat_fprint(stdout, Mz); printf("\n");
       mpz_mat_fprint(stdout, Tz); printf("\n");

       printf("%1.4f\n", t);
       }

       mpz_clear(p);
       mpq_mat_clear(M);
       mpq_mat_clear(T);
       mpz_mat_clear(Mz);
       mpz_mat_clear(Tz);
       gmp_randclear(state);
     */


    // Here starts my personal work
    // We assume that the polynomial was given in a command line, by giving its coefficient, in reversed order (first the constant coefficient, etc. And the head coefficient in the end);

    /*
       if(argc > 1){

       // Initialisation
       unsigned int degree = argc-2;
       mpz_t f[degree+1];
       mpz_mat mul_alpha, M, N, D, D2;
       mpz_t p;
       mpz_t minus;

       // Storing the coefficients obtained in the command line in the polynomial
       unsigned int i,j, k;
       for(i = 0 ; i <= degree ; i++){
       mpz_init(f[i]);
       mpz_set_str(f[i],argv[i+1],10);
       }
       printf("\nYour polynomial is :\n");
       print_polynomial(f,degree);

       //Initialising each matrix
       mpz_mat_init(mul_alpha,degree,degree);
       mpz_mat_init(M,degree,degree);
       mpz_mat_init(N,degree,degree);
       mpz_mat_init(D,degree,degree);
       mpz_mat_init(D2,2*degree,degree);
       mpz_init(p);
       mpz_init(minus);
       mpz_mat_set_ui(M,1); // M starts as the identity Matrix
       mpz_set_si(minus,-1);


       // Filling the coefficients in mul_alpha, the companion matrix of f
       for(i = 0 ; i < degree ; i++){
       // Setting the left part of the companion matrix
       for(j = 0 ; j < degree-1 ; j++){
       if(i == j+1){ mpz_set_ui(mpz_mat_entry(mul_alpha,i,j),1); }
       else{ mpz_set_ui(mpz_mat_entry(mul_alpha,i,j),0); }
       }

       // Computing the coefficients for the column on the right
       mpz_set(p,f[i]);
       mpz_mul(p,p,minus);
       for(k = 1 ; k <= degree-1-i ; k++){
       mpz_mul(p,p,f[degree]);
       }
       mpz_set(mpz_mat_entry(mul_alpha,i,j),p);
       }



       // Filling the coefficients in D, whose determinant must be computed to get the discriminant.
       for(i = 0 ; i < degree ; i++){
       for(j = 0 ; j <= i ; j++){
       mpz_mat_trace(mpz_mat_entry(D, i-j, j), M);
       }
       mpz_mat_multiply(N,M,mul_alpha);
       mpz_mat_swap(M,N);
       }
       for(j = 1 ; j < degree ; j++){
       for(i = degree-1 ; j <= i ; i--){
       mpz_mat_trace(mpz_mat_entry(D,i,j+(degree-1)-i),M);
       }
       mpz_mat_multiply(N,M,mul_alpha);
       mpz_mat_swap(M,N);
       }


       // Preparing the HNF
       mpz_mat_realloc(M,degree,degree);
       int sign = mpz_hnf_backend(D, M);

       mpz_mat_determinant_triangular(p, D);
       mpz_mul_si(p,p,sign);
       gmp_printf("\nThe discriminant of Z[f_d * alpha] is %Zd\n\n",p);

       mpz_clear(p);
       mpz_clear(minus);
       mpz_mat_clear(mul_alpha);
       mpz_mat_clear(M);
       mpz_mat_clear(N);
       mpz_mat_clear(D);
       mpz_mat_clear(D2);
       for(i = 0 ; i <= degree ; i++){
       mpz_clear(f[i]);
       }

       } */

    // The inputs to this problem are f, one polynomial of degree n, and
    // B, the matrix containing the genereators of one order of the
    // number field obtained with f, as well as p, a prime number

    unsigned long seed = clock();

    for( ; argc > 3 ; ) {
        if (strcmp(argv[1], "-s") == 0 || strcmp(argv[1], "--seed") == 0) {
            seed = atoi(argv[2]);
            argc--,argv++;
            argc--,argv++;
            continue;
        }
        fprintf(stderr, "Usage: ./a.out [options] [filename] [p]\n");
        fprintf(stderr, "Unexpected arg: %s\n", argv[1]);
	exit(EXIT_FAILURE);
    }


    if (argc != 3) {
	fprintf(stderr, "Usage: ./a.out [options] [filename] [p]\n");
	exit(EXIT_FAILURE);
    }

    unsigned int p = strtoul(argv[2], NULL, 0);	//19; //atoi(argv[1]);
    FILE *problemfile = fopen(argv[1], "r");

    if (!problemfile) {
        fprintf(stderr, "%s: %s\n", argv[1], strerror(errno));
        exit(EXIT_FAILURE);
    }

    mpq_mat D;
    mpz_poly f;

    printf("Format: [degree] [coeffs] [coeffs of order basis]\n");

    unsigned int n = 0;
    mpz_poly_init(f, n);
    read_data(&n, f, /*gen,*/ problemfile);	
    fclose(problemfile);
    mpq_mat_init(D, n, n);

    mpz_poly g;
    mpz_poly_init(g, n);
    mpz_poly_to_monic(g, f);
    printf("f  is : ");
    mpz_poly_fprintf(stdout, f);
    printf("\n");
    printf("f^ is : ");
    mpz_poly_fprintf(stdout, g);
    printf("\n");

    
    p_maximal_order(D, f, p);
    
    printf("\n");
    printf("the %d-maximal order is \n", p);
    mpq_mat_fprint(stdout, D);
    printf("\n");

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);
        
    vector<pair<cxx_mpq_mat, int>> ideals;
    factorization_of_prime(ideals, g, p, state);
    for(unsigned int i = 0 ; i < ideals.size() ; i++){
        printf("Ideal :\n");
        mpq_mat_fprint(stdout, ideals[i].first);
        printf("with a multiplicity of %d\n\n",ideals[i].second);
    }

    string SBAD = "";
    string SBADINFO = "";
    print_comments_for_badideals_above_p(SBAD, SBADINFO, 0, D,g,ideals,p);
    
    gmp_randclear(state);
    
    mpz_poly_clear(g);
    mpz_poly_clear(f);
    mpq_mat_clear(D);



}

/*}}}*/
