#include "rspoly.h"
#include<stdlib.h>

int main() {

}

// ************************
// BASIC METHODS (ON PLACE)
// ************************
void rspoly_init(rspoly p,int deg) {
	int i;	

	p->deg = deg;
	p->a = malloc((deg+1)*sizeof(mpz_t));

	for(i=0; i<=deg; i++) 
		mpz_init(p->a[i]);
}

void rspoly_clear(rspoly p) {
	int i;

	for( i=0 ; i<=p->deg ; i++)
		mpz_clear(p->a[i]);

	free(p->a);
}

void rspoly_set(rspoly p, mpz_t * coeffs) {

	// We suppose clist has length at least p->deg+1
	int i;

	for( i=0 ; i<=p->deg ; i++)
     		mpz_set(p->a[i],coeffs[i]);
}

void rspoly_init_set(rspoly p, int deg, mpz_t * coeffs) {
	int i;
   	p->deg = deg;
	p->a = malloc((deg+1)*sizeof(mpz_t));

	for( i=0 ; i<=p->deg ; i++) 
		mpz_init_set(p->a[i],coeffs[i]);	
}

void rspoly_set_to_identity(rspoly p) {
	rspoly_clear(p);
	rspoly_init(p,1);
	// The coefficients of p are zero at this point.
	mpz_set_ui(p->a[1],1);
}


void rspoly_copy(rspoly target, rspoly source) {
	
	if(target->deg != source->deg) {
		rspoly_clear(target);
		rspoly_init(target, source->deg);
	}	
	
	int i;
	for( i=0 ; i<=source->deg; i++) 
                mpz_set(target->a[i],source->a[i]);        

}

void rspoly_trim(rspoly p) {
	int i;
	for (i = p->deg; i >= 1; i--)
		if(mpz_sgn(p->a[i])!=0)
			break;
		else
			mpz_clear(p->a[i]);
	p->deg = i;   
}

int  rspoly_is_constant(rspoly p);
void rspoly_constant_coeff(mpz_t res, rspoly p);
int  rspoly_equal(rspoly p, rspoly q);
int  rspoly_equivalent(rspoly p, rspoly q);

// Algebraic operations
void rspoly_diff(rspoly res, rspoly op);
void rspoly_mod_coeffs(rspoly res, rspoly op, unsigned long m);
void rspoly_div_coeffs(rspoly res, rspoly op, unsigned long m);
void rspoly_coeff_product_si(rspoly res, rspoly op, long s);
void rspoly_sum(rspoly res, rspoly op1, rspoly op2);
void rspoly_prod(rspoly res, rspoly op1, rspoly op2);

