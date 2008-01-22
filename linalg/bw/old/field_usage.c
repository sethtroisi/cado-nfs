#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "lingen_params.h"
#include "params.h"
#include <assert.h>
#include "types.h"
#include "macros.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "variables.h"
#include "structure.h"
#include "gmp-hacks.h"
#include "modulus_hacks.h"
#include "field_def.h"
#include "field_prime.h"
#include "field_quad.h"
#include "field_complex.h"
#include "field_usage.h"

struct field * field_k,* field_l;

/* They are reduced anyway */
void restrict_scalars(mp_limb_t * dest, mp_limb_t * src)
{
	if (!k_is_zero(src+k_size)) {
		printf("weird (scalar restriction)\n");
	}
	k_set(dest,src);
}

void extend_scalars(mp_limb_t * dest, mp_limb_t * src)
{
	k_set(dest,src);
	k_set_zero(dest+k_size);
}

#define QNR_LIMIT	128
#define ROOT_LIMIT	128

void prepare_fields_for_fft(int root_order UNUSED_VARIABLE, int enable_cplx)
{
	mpz_t a;
	mp_limb_t * d;
	int i;

	field_k=new_prime_field(modulus_plain,bw_allocsize);
	printf("Looking for a quadratic non-residue\n");
#ifdef	HAS_NATIVE_FFT
	if (enable_cplx==-1) {
		field_l = field_k;
		printf("Assuming roots of 1 are present in the base field\n");
		return;
	}
#endif	/* HAS_NATIVE_FFT */

	mpz_init(a);
	if (enable_cplx==1) {
		mpz_set_si(a,-1);
		if (mpz_jacobi(a,modulus)==-1) {
			printf("-1 is not a root => "
				"we can use a complex field L=K(sqrt(-1)).\n");
			mpz_clear(a);
			field_l=new_complex_field(field_k);
			return;
		}
	}
	/*
	 * The residue is sought from positive numbers since these
	 * are the ones with the smallest representation.
	 */
	for(i=2;i<QNR_LIMIT;i++) {
		mpz_set_ui(a,i);
		if (mpz_jacobi(a,modulus)==-1) break;
	}
	if (i==QNR_LIMIT) {
		fprintf(stderr,"found no quadratic non-residue before limit\n");
		exit(100);
	}
	printf("Found L=K(sqrt(%d))\n",i);
	mpz_clear(a);
	d=malloc(k_size*sizeof(mp_limb_t));
	k_set_int(d,i);
	field_l=new_quad_field(field_k,d);
	free(d);
}

mp_limb_t * fetch_primitive_root(int root_order)
{
	mp_limb_t *t, *r, *a;
	mpz_t	mgorder;
	int	max_order;
	int	i,j;

	printf("Looking for a primitive 2^%d-th root of 1 in L\n",root_order);

	mpz_init(mgorder);
	mpz_set(mgorder,modulus);
	mpz_mul(mgorder,mgorder,mgorder);
	mpz_sub_ui(mgorder,mgorder,1);
	max_order=mpz_scan1(mgorder,0);
	mpz_tdiv_q_2exp(mgorder,mgorder,max_order);

	assert(mgorder[0]._mp_size>0);

	printf("The biggest possible order (among powers of 2) is 2^%d\n",
		       max_order);

	if (root_order > max_order) {
		fprintf(stderr,"Ahem, it won't make it !\n");
		exit(100);
	}

	
	t=malloc(l_size*sizeof(mp_limb_t));
	r=malloc(l_size*sizeof(mp_limb_t));
	a=malloc(l_size*sizeof(mp_limb_t));

	for(i=0;i<ROOT_LIMIT;i++) {
		mpn_random(t,l_size);
		l_reduce(t);
		l_exp(r,t,mgorder[0]._mp_d,mgorder[0]._mp_size);
		l_set(a,r);
		for(j=0;j<root_order-1;j++) {
			l_mul(t,r,r);
			l_set(r,t);
		}
		if (!l_is_one(r)) break;
		printf("Try #%d : small 2-order\n",i);
	}

	/* We'll get it. */
	for(;j<=max_order && ! l_is_one(r);j++) {
		l_mul(t,r,r);
		l_set(r,t);
	}
	assert(j<=max_order);

	printf("Found a root on the %d-th try (its actual 2-order was %d)\n",
			i,j);

	for(i=root_order;i<j;i++) {
		l_mul(t,a,a);
		l_set(a,t);
	}

	free(r);
	free(t);

	printf("Root is: ");
	l_print(a);
	printf("\n");

	return a;
}
