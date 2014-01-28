#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fast_rootsieve.h"
#include "auxiliary.h"
#include "portability.h"
#define PLB 2 // Prime lower bound
#define PUB 2 // Prime upper bound
#define DEGREE 4 
#define ALIM 2
#define VERBOSE 1
#define VERBOSE_RS 1
#define MULTI 10

/* Here we compute alpha in different ways in order
 * to be sure that we are computing the same thing.
 *
 * The first method computes alpha for a given
 * polynomial using the methods directly inspired
 * from Emmanuel's CADO-wide paper (proposition 1).
 * We compute them in the most naive way : we 
 * take the current polynomial and apply the formula
 * of proposition 1. We don't exploit the fact
 * that we are computing consecutive polynomials.
 *
 * The second method is rootsieve : 
 * We compute the root property for many polynomials
 * at once. We consider the generic polynomial 
 * h(x) = f(x) + (j*x+k)*g(x), with g(x) = b*x-m,
 * and we compute, for every prime p less than some bound, 
 * the average p-valuation of the values taken by h(x), 
 * for (j,k) belonging to some rectangle of lower 
 * left corner (J0,K0) and upper right corner (J1,K1).
 * Then we compute the difference between the obtained
 * valuation and the average valuation for a random
 * integer. We add up all this contributions to get alpha,
 * the root property (see Emmanuel's CADO paper for 
 * definitions).
 * 
 * */

double bound_f = BOUND_F, bound_g = BOUND_G, area = AREA;

void trivial_method (mpz_t *f, int d, mpz_t b, mpz_t m, unsigned long alim) {
  FILE * true_alpha = fopen("true_alpha.txt","w");
  long J0, K0, J1, K1, k0, j0, k, j; 
  mpz_poly_t F;

  F->coeff = f;
  F->deg =  d;
// We must obtain the bounds, and declare the array.
  rotate_bounds (F, b, m, &K0, &K1, &J0, &J1, VERBOSE, DEFAULT_L2_METHOD);
  j0 = k0 = 0; 

  for(j=J0; j<=J1; j++) {
    j0 = rotate_aux (f, b, m, j0, j, 1);
    for(k=K0; k<=K1; k++) {
      k0 = rotate_aux (f, b, m, k0, k, 0);
			fprintf(true_alpha,"j=%+5ld k=%+5ld alpha=%+2.6f ", j, k, get_alpha (F, alim));
			//print_coeffs(f,d);
			fprintf(true_alpha,"\n");
		}
	}
  
  rotate_aux (f, b, m, j0, 0, 1);
  rotate_aux (f, b, m, k0, 0, 0);
  fclose(true_alpha);
}

void rootsieve_method (mpz_t *f, int d, mpz_t b, mpz_t m, unsigned long alim, int verbose) {
  FILE * rs_alpha = fopen("rs_alpha.txt","w");
	unsigned long p;
  long J0, K0, J1, K1, k, j; 
  mpz_poly_t F;

  F->coeff = f;
  F->deg = d;
// We must obtain the bounds, and declare the array.
  rotate_bounds (F, b, m, &K0, &K1, &J0, &J1, VERBOSE, DEFAULT_L2_METHOD);
  double * A = malloc((K1-K0+1)*sizeof(double));	
  double * P = malloc((K1-K0+1)*sizeof(double));

	double offset;

	// We obtain the list of lattices from rootsieve
  lattice_list aux = rootsieve(f,d,b,m,J0,J1,K0,K1,alim);
  rootsieve_lattice * lattices = aux.lattices; 
	long nb_lattices = aux.length;
	long n,q;
	//fprintf(stderr,"nb_lattices=%ld \n",nb_lattices);
	rootsieve_lattice current_lattice;

  for (j = J0; j<=J1; ++j) {        

    for (k = K0; k <= K1; k++)      
      A[k - K0] = 0.0; // A[k - K0] will store the value alpha(f + k*g) 

    // COMPUTATION OF ALPHA : START
    // The list of lattices fournished by rootsieve have an associated prime
    // which is increasing, and here we use this fact.
    
    n=0; // We start from the zero-th lattice.

    for (p = PLB; p <= alim; p += 1 + (p & 1))
      if (isprime (p)) {
	//fprintf(stderr,"p=%lu lattices[n].p=%lu \n",p,lattices[n].p);
	offset = (log(p))/((double)(p-1));
	for(k = K0; k <= K1; k++)
	  P[k-K0] = offset;

	while(lattices[n].p == p && n < nb_lattices) {
	  //fprintf(stderr,"n=%ld p=%lu \n",n,p);
	  current_lattice = lattices[n];
	  ASSERT_ALWAYS(current_lattice.u0>=J0 && current_lattice.v0>=K0 
			&& current_lattice.vs>=0 && current_lattice.us>=0 && 
			current_lattice.skew>=0);  	  
	  
	  if(current_lattice.u0 == j) {
	    if (verbose>29) {
	      fprintf(stderr,"[u0=%ld v0=%ld us=%ld vs=%ld skew=%ld contrib=%2.7f\n",
		     current_lattice.u0,current_lattice.v0,current_lattice.us,
		     current_lattice.vs,current_lattice.skew,current_lattice.contribution);  
	    }
	    
	    while(current_lattice.v0 <= K1) {
	      if (verbose>19) {
		fprintf(stderr,"j=u0=%ld,k=v0=%ld,contribution=%2.7f.\n",j,current_lattice.v0,current_lattice.contribution);
	      }
	      P[current_lattice.v0 - K0] -= current_lattice.contribution;
	      current_lattice.v0 += current_lattice.vs;
	    }
	    
	    current_lattice.v0 += current_lattice.skew;
	    current_lattice.u0 += current_lattice.us;
	    q = (long)floor(((double)(current_lattice.v0-K0))/((double)current_lattice.vs));
	    current_lattice.v0 -= q*current_lattice.vs;
	
	    if (verbose>29) {
	      fprintf(stderr,"u0=%ld v0=%ld us=%ld vs=%ld skew=%ld contrib=%2.7f]\n",
		     current_lattice.u0,current_lattice.v0,current_lattice.us,
		     current_lattice.vs,current_lattice.skew,current_lattice.contribution);  
	    }
	    ASSERT_ALWAYS(current_lattice.v0>=K0);
	  }
	  
	  n++;
	}
	
	for(k = K0; k <= K1; k++) 
	  A[k-K0] += P[k-K0];

      }

		for(k=K0; k<=K1; k++) 
      fprintf(rs_alpha,"j=%+5ld k=%+5ld alpha=%+2.6f \n", j, k, A[k-K0]);				
      
	}

  fclose(rs_alpha);
	free (A);
	free (P);
	free (lattices);

}



int main() {  
  
  mpz_t m,b;
  int d;
  mpz_t *f;
  long *jmin,*kmin;
  int multi;

  d = DEGREE;  
  multi = MULTI;

  jmin = malloc(multi*sizeof(long));
  kmin = malloc(multi*sizeof(long)); 

  f = malloc((d+1)*sizeof(mpz_t));
  
  mpz_init_set_str(f[4],"162660",10);
  mpz_init_set_str(f[3],"-66166894",10);
  mpz_init_set_str(f[2],"1494525285",10);
  mpz_init_set_str(f[1],"25965416770272",10);
  mpz_init_set_str(f[0],"-387855863207419",10);  

  mpz_init_set_str(b,"29993273",10);
  mpz_init_set_str(m,"25764542337264",10);

  // We add a g(x) term to create the challenge.
  //mpz_add(f[2],f[2],b);
  //mpz_sub(f[1],f[1],m);
  mpz_add(f[1],f[1],b);
  mpz_sub(f[0],f[0],m);

  trivial_method (f, d, b, m, ALIM);
  //rotate_method  (f, d, ALIM, m,  b, jmin, kmin, MULTI, VERBOSE);
	rootsieve_method (f, d, b, m, ALIM, VERBOSE);

  free(jmin);
  free(kmin);  
  free(f);

	return 0;
}

