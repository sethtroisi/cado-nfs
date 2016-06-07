#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <stdint.h>
#include <sys/resource.h>	/* for getrusage */
#include <queue>
#include <time.h>
#include "macros.h"
#include "mpz_poly.h"
#include "mpz_mat.h"

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
void mpz_mat_row_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, unsigned int i)
{
	mpz_poly_clear(f);
	mpz_poly_init(f,M->n-1);
	unsigned int j;
	for (j = 0 ; j < M->n; j++){
		mpz_poly_setcoeff(f,j,mpz_mat_entry_const(M,i,j));
	}
}

void mpz_mat_column_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, unsigned int j)
{
	mpz_poly_clear(f);
	mpz_poly_init(f,M->m-1);
	unsigned int i;
	for (i = 0 ; i < M->m; i++){
		mpz_poly_setcoeff(f,i,mpz_mat_entry_const(M,i,j));
	}
}

void mpq_mat_row_to_poly(mpz_poly_ptr f, mpz_ptr denom, mpq_mat_srcptr M, unsigned int i)
{
	mpz_poly_clear(f);
	mpz_poly_init(f,M->n-1);
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

void mpq_mat_column_to_poly(mpz_poly_ptr f, mpz_ptr denom, mpq_mat_srcptr M, unsigned int j)
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
			   mpz_poly_srcptr g, unsigned int p)
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

// Builds the block matrix containing p*identity in the top, and K in the bottom
// Then computes its HNF and stores it in I
void join_HNF(mpz_mat_ptr I, mpz_mat_srcptr K, unsigned int p)
{
    mpz_mat J, T0;
    mpz_mat_init(J, K->n, K->n);
    mpz_mat_init(T0, K->n, K->n);
    mpz_mat_realloc(I, K->n, K->n);

    mpz_mat_set_ui(J, 1);
    mpz_mat_multiply_by_ui(J, J, p);
    mpz_mat_vertical_join(I, J, K);
    mpz_hnf_backend(I, T0);
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
				  unsigned int p)
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

    for (unsigned int i = 0; i < n; i++) {
	for (unsigned int j = 0; j < n; j++) {

	    // Gammma : polynomial of generator of Ip (only the numerators, multiplied by LCM of denominators)
	    // c : same for generators of O
	    mpz_poly gamma, c, aux;
	    // LCM of denominators of generators of Ip, of O, and one auxiliary coefficient
	    mpz_t denom_g, denom_c, coeff;
	    // Line matrices containing the product gamma*c mod g
	    mpz_mat row, aux_row;
	    // Auxiliary row matrices, necessary because generators of Ip need to be converted in the basis of alpha^
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

	    // Storing one generator of B in the polynomial gamma and in denom_g
	    mpq_mat_row_to_poly(gamma, denom_g, B, i);
	    //Storing one generator of I_p in the polyomial c and in denom_c
	    mpz_mat_submat_swap(aux_row, 0, 0, I, j, 0, 1, n);
	    mpz_mat_set(row, aux_row);
	    mpz_mat_submat_swap(aux_row, 0, 0, I, j, 0, 1, n);
	    mpz_mat_to_mpq_mat(row_q, row);
	    mpq_mat_multiply(c_mat, row_q, B);
	    mpq_mat_row_to_poly(c, denom_c, c_mat, 0);

	    // Computing gamma*c mod g
	    mpz_poly_mul_mod_f(aux, gamma, c, g);
	    // Storing the result in row_q
	    for (int k = 0; k <= aux->deg; k++) {
		mpz_poly_getcoeff(coeff, k, aux);
		mpq_set_num(mpq_mat_entry(row_q, 0, k), coeff);
		mpz_mul(coeff, denom_c, denom_g);
		mpq_set_den(mpq_mat_entry(row_q, 0, k), coeff);
		mpq_canonicalize(mpq_mat_entry(row_q, 0, k));
	    }

	    // Converting row_q (gamma*c mod g) in the basis of I_p (it is supposed to contain only integers)
	    mpq_mat_multiply(res, row_q, B_inv);
	    mpq_mat_multiply(res, res, I_inv);


	    // Extracting the numerators (remember, integers only) into a mpz_mat
	    mpq_mat_numden(row, coeff, res);
	    // Computing the same matrix, modulo p (it is supposed to be a vector of n integers, associated to gamma)
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

void read_data(unsigned int *deg, mpz_poly_ptr f, mpq_mat_ptr B,
	       FILE * problemfile)
{
    fscanf(problemfile, "%u", deg);
    
    mpz_poly_realloc(f, *deg + 1);
    mpz_t c;
    mpz_init(c);
    for (unsigned int i = 0; i <= *deg; i++) {
        gmp_fscanf(problemfile, "%Zd", &c);
        mpz_poly_setcoeff (f, i, c);
    }
    mpz_clear(c);

    mpq_mat_realloc(B, *deg, *deg);
    for (unsigned int i = 0; i < *deg; i++) {
        int denom;
        denom = 1;
        //fscanf(problemfile, "%d", &denom);
        for (unsigned int j = 0; j < *deg; j++) {
            int c;
            if(i == j){c = 1;}
            else{c = 0;}
            //fscanf(problemfile, "%d", &c);
            mpq_set_si(mpq_mat_entry(B, i, j), c, denom);
        }
    }
}

// Stores in I the p-radical of the order generated by B (in the basis of these generators)
void p_radical_of_order(mpz_mat_ptr I, mpq_mat_ptr B, mpz_poly_srcptr g, unsigned int p)
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

    // Now building the matrix U, containing all generators to the power of p
    // Generators are polynomials, stored in the matrix B
    generators_to_power_p(U, B, g, p);

     // Now do U * B^-1, storing in X the application  F : z -> (z^p mod g mod p)
    mpq_mat_multiply(T, U, B_inv);
    mpq_mat_numden(X, den, T);
    mpz_mat_mod_ui(X, X, p);

    /* Which power of p ? */
    int k = 1;
    for (unsigned int pk = p; pk < n; k++, pk *= p);
    mpz_mat_power_ui_mod_ui(X, X, k, p);

    // Storing in K a basis of Ker((z -> (z^%p mod g mod p))^k)
    mpz_mat_kernel(K, X, p);

    // Getting generators of the p radical from Ker(X) by computing HNF of the vertical block matrix (p*Id, K);
    join_HNF(I, K, p);
    
    mpz_clear(den);
    mpz_mat_clear(K);
    mpz_mat_clear(X);
    mpq_mat_clear(T);
    mpq_mat_clear(U);
    mpq_mat_clear(B_inv);
}

// Builds one p-maximal-order starting from the matrix B, containing the generators of one order (in the basis of alpha^)
// In the number field of the polynomial f, of degree n
// Stores the generators of this p-maximal-order (in the basis of alpha^) in D
void p_maximal_order(mpq_mat_ptr D, mpq_mat_srcptr B, mpz_poly_srcptr f,
		     unsigned int p)
{
    ASSERT_ALWAYS(B->m == B->n);

    mpq_mat new_D;
    mpz_poly g;
    mpq_t p_inv;
    unsigned int n = B->n;


    mpq_mat_init(new_D, n, n);
    mpz_poly_init(g, n);	// Monic polynomial corresponding to g
    mpq_init(p_inv);		// 1/p

    // Storing in g the monic polynomial such as (fd*alpha) is a root of if and only if alpha is a root of f
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

        // Computing Ker(M)
        mpz_mat_kernel(K_M, M, p);

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
    mpz_hnf_backend(D_z,T);
    mpz_mat_to_mpq_mat(D,D_z);
    mpq_mat_multiply_by_mpq(D,D,denom_inv);
    
    mpq_clear(denom_inv);
    mpz_clear(denom_eq_1);
    mpz_mat_clear(D_z);
    mpz_mat_clear(T);
    mpq_clear(denom_q);
    mpz_clear(denom);
    mpq_mat_clear(new_D);
    mpz_poly_clear(g);
    mpq_clear(p_inv);
}

// W is the matrix containing the generators of one p-maximal order Ok, in the basis of alpha^ (root of g)
// theta is a row matrix representing an element of Ok/pOk, thus a vector on the basis (W[0], ... W[n-1]) (with coeffs mod p)
// g is the monic polynom defining the number field in which we are
void minimal_poly_of_mul_by_theta(mpz_poly_ptr f, mpq_mat_srcptr W, mpz_mat_srcptr theta, mpz_poly_srcptr g, unsigned int p)
{
    ASSERT_ALWAYS(W->m == W->n);
    unsigned int n = W->m;
    ASSERT_ALWAYS((theta->m == 1) && (theta->n == n));
    
    // First let's build the matrix of multiplication by theta, in the basis of (W[0], ... W[n-1])
    mpz_mat times_theta;
    mpq_mat times_theta_rat;
    mpz_mat_init(times_theta,n,n);
    mpq_mat_init(times_theta_rat,n,n);
    
    mpq_mat theta_rat; // Will contain theta in the basis of alpha^
    mpq_mat_init(theta_rat,1,n);
    mpz_mat_to_mpq_mat(theta_rat,theta);
    mpq_mat_multiply(theta_rat,theta_rat,W); // Now contains theta in the basis of alpha^ (rationnal coefficients)

    // Converting theta into one polynom and one common denominator
    mpz_poly theta_poly;
    mpz_t theta_denom;
    mpz_poly_init(theta_poly,n-1);
    mpz_init(theta_denom);
    mpq_mat_row_to_poly(theta_poly,theta_denom,theta_rat,0);

    //printf("theta : "); mpz_mat_fprint(stdout,theta);
    //printf("theta rat : "); mpq_mat_fprint(stdout,theta_rat);
    printf("theta poly : "); mpz_poly_fprintf(stdout,theta_poly);
    gmp_printf("theta denom = %Zd\n",theta_denom);
    
    // Computing theta*w[i] (in the basis of alpha^ ) for each i
    for (unsigned i = 0 ; i < n ; i++) {
        // Converting w[i] (already in the basis of alpha^) into one polynom and one common denominator
        mpz_poly w_poly;
        mpz_t w_denom;
        mpz_poly_init(w_poly,n-1);
        mpz_init(w_denom);
        mpq_mat_row_to_poly(w_poly,w_denom,W,i);
        
        //printf("w[%d] poly : ",i); mpz_poly_fprintf(stdout,w_poly);
    
        // Computing theta (in the basis of alpha)^ * w[i] mod g, without the common denominators.
        // We divide at the end
        mpz_poly res;
        mpz_poly_init(res,n-1);        
        mpz_poly_mul_mod_f(res,theta_poly,w_poly,g);
        mpz_poly_cleandeg(res,n-1); // In case the degree of the product isn't n-1
        
        //printf("res : "); mpz_poly_fprintf(stdout,res);
        
        for (int j = 0 ; j <= res->deg ; j++) {
            mpz_t coeff, denom;
            mpz_init(coeff);
            mpz_init(denom);
            
            // Storing the product into one row of times_theta_rat, then dividing by the common denominators
            mpz_poly_getcoeff(coeff,j,res);
            mpq_set_num(mpq_mat_entry(times_theta_rat,i,j),coeff);
            mpz_mul(denom,theta_denom,w_denom);
            mpq_set_den(mpq_mat_entry(times_theta_rat,i,j),denom);
            mpq_canonicalize(mpq_mat_entry(times_theta_rat,i,j));
            
            
            mpz_clear(denom);
            mpz_clear(coeff);
        }
        
        mpz_poly_clear(res);
        mpz_clear(w_denom);
        mpz_poly_clear(w_poly);
    }
    
    // Now we have to multiply by W^-1, in order to get into the basis of (w[0], ... w[n-1]) again
    mpq_mat W_inv;
    mpq_mat_init(W_inv,n,n);
    mpq_mat_invert(W_inv,W);
    mpq_mat_multiply(times_theta_rat,times_theta_rat,W_inv);
    
    // Now we have to convert it into a matrix of integers
    mpz_t denom_eq_1;
    mpz_init(denom_eq_1);
    mpq_mat_numden(times_theta,denom_eq_1,times_theta_rat);
    
    // And mod p
    mpz_mat_mod_ui(times_theta,times_theta,p);

    printf("matrix of multiplication by theta, mod %d:\n", p);
    mpz_mat_fprint(stdout,times_theta); printf("\n");
    
    // Now starting to compute the (n,n^2) matrix whose kernel will be computed
    mpz_mat M,current;
    mpz_mat_init(M,n+1,n*n);
    mpz_mat_init(current,n,n);
    mpz_mat_set_ui(current,1);
    for (unsigned int k = 0 ; k < n+1 ; k++){
        for (unsigned int i = 0 ; i < n ; i++){
            for (unsigned int j = 0 ; j < n ; j++){
                mpz_set(mpz_mat_entry(M,k,j+n*i),mpz_mat_entry(current,i,j));
            }
        }
        mpz_mat_mod_ui(M,M,p);
        mpz_mat_multiply(current,current,times_theta);
        mpz_mat_mod_ui(current,current,p);
    }
    
    printf("big matrix, mod %d:\n", p);
    //mpz_mat_fprint(stdout,M); printf("\n");
    
    // Now computing its kernel
    mpz_mat K;
    mpz_mat_init(K,0,0);
    mpz_mat_kernel(K,M,p);

    printf("kernel of (n+1,n^2) matrix, mod %d:\n", p);
    mpz_mat_fprint(stdout,K); printf("\n");

    // Getting the minimal polynomial and verifying that f(M) = 0
    if(K->m == 0){
        mpz_poly_realloc(f,0);
    }
    else{
        mpz_mat_row_to_poly(f,K,0);
    }
    
    // Testing
    mpz_mat test_mat;
    mpz_mat_init(test_mat,n,n);
    mpz_mat_in_poly_mod_ui(test_mat,times_theta,f,p);
    //printf("f(times_theta) :\n");
    //mpz_mat_fprint(stdout,test_mat); printf("\n");
    mpz_mat_clear(test_mat);
    
    
    printf("minimal polynomial of times_theta mod %d:\n", p);
    mpz_poly_fprintf(stdout,f); printf("\n");
    
    mpz_mat_clear(K);
    mpz_mat_clear(current);
    mpz_mat_clear(M);
    mpz_clear(denom_eq_1);
    mpq_mat_clear(W_inv);
    mpz_clear(theta_denom);
    mpz_poly_clear(theta_poly);
    mpq_mat_clear(theta_rat);
    mpq_mat_clear(times_theta_rat);
    mpz_mat_clear(times_theta);
    
}

// Represents the type of elements found in pick_from, in badideals.mag
struct subspace_ideal {
    cxx_mpz_mat E; // the basis of the subspace E
    cxx_mpz_mat I; // the basis of one ideal
    cxx_mpz_mat V;
//    unsigned int dim; // The dimension of E, instead of having the full space like in magma
};

void factorization_of_prime(/*vector<pair<cxx_mpz_mat, int>>& res,*/ mpz_poly_srcptr g, unsigned int p)
{
    int n = g->deg;
    
    mpq_mat G, G_inv, Id_q; // G = Basis of p-maximal-order Ok ; Id_q = identity (in rational domain)
    mpz_mat Ip; // p radical of Ok
    cxx_mpz_mat Id_z; // Identity (in integers domain)
    mpq_mat_init(G,n,n);
    mpq_mat_init(G_inv,n,n);
    mpq_mat_init(Id_q,n,n);
    mpz_mat_init(Ip,n,n);
    mpz_mat_realloc(Id_z,n,n);
    mpq_mat_set_ui(Id_q,1);
    mpz_mat_set_ui(Id_z,1);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,clock());
        
    // Computing the generators of p-maximal order, storing them in G
    p_maximal_order(G,Id_q,g,p);
    // Computing the p-radical of the order of G, storing it in Ip
    p_radical_of_order(Ip,G,g,p);
    mpq_mat_invert(G_inv,G);
    
    
    // The initial value in pick_from
    subspace_ideal initial;
    initial.E = Id_z; // It is the basis of Ok/p*Ok, but written on the basis of Ok/p*Ok ; thus, it's the identity
    initial.V = Id_z; // identity
    initial.I = Id_z; //p*identity
    mpz_mat_multiply_by_ui(initial.I,initial.I,p);
    
    // The vector on which the recursion will happen.
    queue<subspace_ideal> pick_from;
    pick_from.push(initial);
    
    
    printf("W = \n"); mpq_mat_fprint(stdout,G); printf("\n");
    printf("W_inv = \n"); mpq_mat_fprint(stdout,G_inv); printf("\n");
    
    while (pick_from.size() > 0)  {
        mpz_poly f;
        mpz_poly_init(f,n);
        
        // Getting the head element
        subspace_ideal current = pick_from.front();
        pick_from.pop();
        
        // Picking one random element of subspace E
        cxx_mpz_mat c;
        mpz_mat_realloc(c,1,n);
        if (n == 8) {
            mpz_set_ui(mpz_mat_entry(c,0,0),0);//0);
            mpz_set_ui(mpz_mat_entry(c,0,1),1);//0);
            mpz_set_ui(mpz_mat_entry(c,0,2),0);//0);
            mpz_set_ui(mpz_mat_entry(c,0,3),2);//1);
            mpz_set_ui(mpz_mat_entry(c,0,4),1);//2);
            mpz_set_ui(mpz_mat_entry(c,0,5),0);//2);
            mpz_set_ui(mpz_mat_entry(c,0,6),0);//0);
            mpz_set_ui(mpz_mat_entry(c,0,7),0);//1);
        } else {
            for (int i = 0 ; i < n ; i++) {
                mpz_set_ui(mpz_mat_entry(c,0,i),gmp_urandomm_ui(state,p));
            }
        }
        minimal_poly_of_mul_by_theta(f,G,c,g,p);
        mpz_poly_cleandeg(f,n);
        
        printf("Element c is : "); mpz_mat_fprint(stdout,c); printf("\n");
        printf("Minimal polynomial is : "); mpz_poly_fprintf(stdout,f); printf("\n");
        mpz_poly_clear(f);
    }
    
    gmp_randclear(state);
    mpz_mat_clear(Ip);
    mpq_mat_clear(G_inv);
    mpq_mat_clear(G);
    mpq_mat_clear(Id_q);
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

    // The inputs to this problem are f, one polynomial of degree n, and B, the matrix containing the genereators of one order of the number field obtained with f, as well as p, a prime number

    if (argc != 3) {
	fprintf(stderr, "Usage: ./a.out [filename] [p]\n");
	exit(EXIT_FAILURE);
    }

    unsigned int p = strtoul(argv[2], NULL, 0);	//19; //atoi(argv[1]);
    FILE *problemfile = fopen(argv[1], "r");

    mpq_mat B, D;
    mpz_poly f;

    printf("Format: [degree] [coeffs] [coeffs of order basis]\n");

    unsigned int n = 0;
    mpz_poly_init(f, n);
    mpq_mat_init(B, 0, 0);	// The matrix of generators
    read_data(&n, f, B, problemfile);	/* Read the data in the file, and reallocate enough space for f and B
					   And sets the value of n, and fills in f and B */
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

    
    p_maximal_order(D, B, f, p);
    /*
    printf("Starting from\n");
    mpq_mat_fprint(stdout, B);
    printf("\n");
    printf("the %d-maximal order is \n", p);
    mpq_mat_fprint(stdout, D);
    printf("\n");
    */
    
    factorization_of_prime(g,p);
    
    
    mpz_poly_clear(g);
    mpz_poly_clear(f);
    mpq_mat_clear(B);
    mpq_mat_clear(D);



}

/*}}}*/
