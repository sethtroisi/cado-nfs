#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <gmp.h>
#include "master_params.h"
#include "params.h"
#include <assert.h>
#include "types.h"
#include "macros.h"
#include "bw_scalar.h"
#include "variables.h"
#include "gmp-hacks.h"
#include "modulus_hacks.h"
#include "auxfuncs.h"


#define OPT_POWEROF2
/* #define CHECK_MULT_BUGS */

#ifdef ACTIVATE_HACKS
static unsigned int order_of_2;
static mp_limb_t high_precomp_mask;

void do_modulus_precomps(void)
{
	int i;
	mp_limb_t * p;

	p = malloc((bw_allocsize+1)*sizeof(mp_limb_t));
	memset(p,0,(bw_allocsize+1)*sizeof(mp_limb_t));
	p[0]=2;
	for(i=1;i<(1<<20)&&(p[0]!=1||mpn_cmp(p+1,zero,bw_allocsize-1)!=0);i++) {
		mpn_lshift(p,p,bw_allocsize+1,1);
		if (p[bw_allocsize] || mpn_cmp(p,modulus_plain,bw_allocsize)>0)
			mpn_sub_n(p,p,modulus_plain,bw_allocsize);
		
	}
	if (i==(1<<20)) {
		order_of_2=0;
		free(p);
		return;
	}

	printf("The multiplicative order of 2 modulo our modulus is %d\n",i);
#ifdef OPT_POWEROF2
	if (i / mp_bits_per_limb != bw_allocsize - 1 && i % mp_bits_per_limb) {
		printf("We cannot take advantage of this fact, though\n");
		order_of_2=0;
	} else {
		order_of_2=i;
		high_precomp_mask=
			(((mp_limb_t)1)<<(order_of_2 % mp_bits_per_limb)) - 1;
	}
#else	/* OPT_POWEROF2 */
	order_of_2=0;
#endif	/* OPT_POWEROF2 */
	free(p);
}

#endif /* ACTIVATE_HACKS */

#ifndef ACTIVATE_HACKS
/* Crap is seldomly efficient. */
void addmul(mp_limb_t * dest, mp_limb_t * s1, mp_limb_t * s2)
{
	mp_limb_t * tmp;

	tmp = FAST_ALLOC((2*bw_allocsize+1)*sizeof(mp_limb_t));

	mpn_mul_n(tmp,s1,s2,bw_allocsize);
	mpn_add(tmp,tmp,2*bw_allocsize,dest,bw_allocsize);

	mpn_tdiv_qr(	tmp+bw_allocsize, dest, 0,
			tmp, 2*bw_allocsize,
			modulus_hbs, bw_allocsize);
	FAST_FREE(tmp);

}
#else
/* Enable specialized code to replace the tdiv_qr. */
void addmul(mp_limb_t * dest, mp_limb_t * s1, mp_limb_t * s2)
{
	mp_limb_t * tmp;
	mp_limb_t  carry0, carry1;

	tmp = FAST_ALLOC((2*bw_allocsize+1)*sizeof(mp_limb_t));

	mpn_mul_n(tmp,s1,s2,bw_allocsize);
	mpn_add(tmp,tmp,2*bw_allocsize,dest,bw_allocsize);

	if (order_of_2 && (mp_bits_per_limb % order_of_2 == 0)) {
		/* Just an addition... */
		if (mpn_add_n(dest,tmp,tmp+bw_allocsize,bw_allocsize)) {
			dest[0]++;
		}
	} else if (order_of_2) {
		carry0=tmp[bw_allocsize-1] & high_precomp_mask;
		carry1=mpn_lshift(	tmp+bw_allocsize-1,
					tmp+bw_allocsize-1,
					bw_allocsize+1,
			mp_bits_per_limb-(order_of_2 % mp_bits_per_limb));
		tmp[bw_allocsize-1] = carry0;
		assert(carry1==0);
		carry1=mpn_add_n(tmp,tmp,tmp+bw_allocsize,bw_allocsize);
		assert(carry1==0);
		/* overflow is at most 1bit . Since we know we have at
		 * least 1bit free space in the limb, it's ok */
				
#ifdef CHECK_MULT_BUGS
		if ((tmp[0] >> (order_of_2 % mp_bits_per_limb))>1) {
			fprintf(stderr,"Multiplication error\n");
			eternal_sleep();
		}
#endif	/* CHECK_MULT_BUGS */
		/* The result of the shift is guaranteed to fit in
		 * bw_allocsize because the product was at most
		 * 2*order_of_2 bits long, and we are facing this product
		 * >> order_of_2 less one limb. */
		if (mpn_cmp(tmp,modulus_plain,bw_allocsize)>=0) {
				mpn_sub_n(dest,tmp,modulus_plain,bw_allocsize);
		} else {
				memcpy(dest,tmp,bw_allocsize*sizeof(mp_limb_t));
		}
	} else {
		mpn_tdiv_qr(	tmp+bw_allocsize, dest, 0,
				tmp, 2*bw_allocsize,
				modulus_hbs, bw_allocsize);
	}
	FAST_FREE(tmp);
}
#endif

#if 0
OUT void addmul(mp_limb_t * dest, mp_limb_t * s1, mp_limb_t * s2)
OUT {
OUT 	mp_limb_t * tmp;
OUT 	mp_limb_t  carry0, carry1;
OUT #ifdef CHECK_MULT_BUGS
OUT 	mp_limb_t * buf2;
OUT #endif	/* CHECK_MULT_BUGS */
OUT 
OUT #if defined(CHECK_MULT_BUGS) || !defined(ACTIVATE_HACKS)
OUT #ifdef TDIV_QR_OVERLAP_OK
OUT 	tmp = alloca((2*bw_allocsize+1)*sizeof(mp_limb_t));
OUT #else	/* TDIV_QR_OVERLAP_OK */
OUT 	mp_limb_t * trash;
OUT 	trash=alloca((bw_allocsize+1)*sizeof(mp_limb_t));
OUT 	tmp = alloca(2*bw_allocsize*sizeof(mp_limb_t));
OUT #endif	/* TDIV_QR_OVERLAP_OK */
OUT #else	/* CHECK_MULT_BUGS || !ACTIVATE_HACKS */
OUT 	tmp = alloca(2*bw_allocsize*sizeof(mp_limb_t));
OUT #endif
OUT 
OUT 	/* One thing to note here: addmul is not the true multiplication
OUT 	 * function anymore. This rather seems to boil down to mul, then
OUT 	 * add...
OUT 	 */
OUT 	/* useless ! memset(tmp,0,2*bw_allocsize*sizeof(mp_limb_t)); */
OUT 	mpn_mul_n(tmp,s1,s2,bw_allocsize);
OUT 	mpn_add(tmp,tmp,2*bw_allocsize,dest,bw_allocsize);
OUT 	/* XXX hack.
OUT 	 * carry==1 iff the most significant bw_alllocsize limbs are
OUT 	 * 0xffffffff ; we know this won't happen because (2^n-1)^2 div
OUT 	 * 2^n=2^n-2, in other words since 2^n is not an eligible modulus
OUT 	 * (not a field...), we are restricted to 0..2^n-1, and even if
OUT 	 * we add 2^n-1 to this, the result won't exceed 2^2n.
OUT 	 */
OUT 
OUT #ifdef CHECK_MULT_BUGS
OUT 	buf2=alloca(bw_allocsize*sizeof(mp_limb_t));
OUT #ifdef TDIV_QR_OVERLAP_OK
OUT 	mpn_tdiv_qr(	tmp+bw_allocsize, buf2, 0,
OUT 			tmp, 2*bw_allocsize,
OUT 			modulus_plain, bw_allocsize);
OUT #else	/* TDIV_QR_OVERLAP_OK */
OUT 	mpn_tdiv_qr(	trash, buf2, 0,
OUT 			tmp, 2*bw_allocsize,
OUT 			modulus_plain, bw_allocsize);
OUT #endif	/* TDIV_QR_OVERLAP_OK */
OUT #endif	/* CHECK_MULT_BUGS */
OUT 	
OUT #ifdef ACTIVATE_HACKS
OUT 	if (order_of_2 && (mp_bits_per_limb % order_of_2 == 0)) {
OUT 		if (mpn_add_n(dest,tmp,tmp+bw_allocsize,bw_allocsize)) {
OUT 			dest[0]++;
OUT 		}
OUT 	} else if (order_of_2) {
OUT 		carry0=tmp[bw_allocsize-1] & high_precomp_mask;
OUT 		carry1=mpn_lshift(	tmp+bw_allocsize-1,
OUT 					tmp+bw_allocsize-1,
OUT 					bw_allocsize+1,
OUT 			mp_bits_per_limb-(order_of_2 % mp_bits_per_limb));
OUT 		tmp[bw_allocsize-1] = carry0;
OUT 		assert(carry1==0);
OUT 		carry1=mpn_add_n(tmp,tmp,tmp+bw_allocsize,bw_allocsize);
OUT 		assert(carry1==0);
OUT 		/* overflow is at most 1bit . Since we know we have at
OUT 		 * least 1bit free space in the limb, it's ok */
OUT 				
OUT #ifdef CHECK_MULT_BUGS
OUT 		if ((tmp[0] >> (order_of_2 % mp_bits_per_limb))>1) {
OUT 			fprintf(stderr,"Multiplication error\n");
OUT 			eternal_sleep();
OUT 		}
OUT #endif	/* CHECK_MULT_BUGS */
OUT 		/* The result of the shift is guaranteed to fit in
OUT 		 * bw_allocsize because the product was at most
OUT 		 * 2*order_of_2 bits long, and we are facing this product
OUT 		 * >> order_of_2 less one limb. */
OUT 		if (mpn_cmp(tmp,modulus_plain,bw_allocsize)>=0) {
OUT 				mpn_sub_n(dest,tmp,modulus_plain,bw_allocsize);
OUT 		} else {
OUT 				memcpy(dest,tmp,bw_allocsize*sizeof(mp_limb_t));
OUT 		}
OUT 	} else
OUT #endif /* ACTIVATE_HACKS */
OUT 		
OUT #ifdef TDIV_QR_OVERLAP_OK
OUT 	mpn_tdiv_qr(	tmp+bw_allocsize, dest, 0,
OUT 			tmp, 2*bw_allocsize,
OUT 			modulus_plain, bw_allocsize);
OUT #else	/* TDIV_QR_OVERLAP_OK */
OUT 	mpn_tdiv_qr(	trash, dest, 0,
OUT 			tmp, 2*bw_allocsize,
OUT 			modulus_plain, bw_allocsize);
OUT #endif	/* TDIV_QR_OVERLAP_OK */
OUT 
OUT #ifdef CHECK_MULT_BUGS
OUT 	if (mpn_cmp(dest,buf2,bw_allocsize)!=0) {
OUT 		fprintf(stderr,"Multiplication error\n");
OUT 		eternal_sleep();
OUT 	}
OUT #endif	/* CHECK_MULT_BUGS */
OUT }
#endif

void
inverse(mp_limb_t * dest, mp_limb_t * src)
{
	mp_limb_t * tmp1,* tmp2;
	mp_limb_t * tmp3,* tmp4;
	mp_size_t s3,s4;

	tmp1=FAST_ALLOC((bw_allocsize+2)*sizeof(mp_limb_t));
	tmp2=FAST_ALLOC((bw_allocsize+2)*sizeof(mp_limb_t));
	tmp3=FAST_ALLOC((bw_allocsize+2)*sizeof(mp_limb_t));
	tmp4=FAST_ALLOC((bw_allocsize+2)*sizeof(mp_limb_t));
	memcpy(tmp1,src,bw_allocsize);
	memcpy(tmp2,modulus_plain,bw_allocsize);
	/* (cosmetic)
	 * memset(tmp1+bw_allocsize,0,2);
	 * memset(tmp2+bw_allocsize,0,2);
	*/
	tmp1[bw_allocsize]=mpn_add_n(tmp1,tmp1,tmp2,bw_allocsize);

	s3=mpn_gcdext(tmp3,tmp4,&s4,tmp1,bw_allocsize,tmp2,bw_allocsize);
	assert(s3==1 && tmp3[0]==1);
	if (s4<0) {
		s4=-s4;
		mpn_tdiv_qr(	tmp1+bw_allocsize,tmp1,0,
				tmp4,s4,
				modulus_plain,bw_allocsize);
		/* It is guaranteed that the remainder will not be zero
		 */
		mpn_sub_n(dest,modulus_plain,tmp1,bw_allocsize);
	} else {
		mpn_tdiv_qr(	tmp1+bw_allocsize,dest,0,
				tmp4,s4,
				modulus_plain,bw_allocsize);
	}
	FAST_FREE(tmp1);
	FAST_FREE(tmp2);
	FAST_FREE(tmp3);
	FAST_FREE(tmp4);
}

