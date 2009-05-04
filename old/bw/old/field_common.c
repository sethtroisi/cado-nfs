#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "lingen_params.h"
#include "params.h"
#include "types.h"
#include "macros.h"
#include "auxfuncs.h"
#include "gmp-hacks.h"
#include "modulus_hacks.h"
#include "field_def.h"
#include "field_common.h"

void generic_field_set(mp_limb_t * a, mp_limb_t * b, struct field * f)
{
	memcpy(a,b,f->size*sizeof(mp_limb_t));
}

int generic_field_eq(mp_limb_t * a, mp_limb_t * b, struct field * f)
{
	f->reduce(a,f);
	f->reduce(b,f);
	return mpn_cmp(a,b,f->size)==0;
}

void generic_field_destroy(struct field * f)
{
	if (f->modulus)	free(f->modulus);
	if (f->minpoly)	free(f->minpoly);
	free(f);
}

void generic_field_exp(mp_limb_t * a, mp_limb_t * r,
		mp_limb_t * x, mp_size_t n,
		struct field * f)
{
	mp_limb_t * u;
	mp_size_t i,j,lead;	/* it is a signed type */
	mp_limb_t mask;

	for(i=n-1;i>=0 && x[i]==0;i--);
	if (i<0) return;

	j=mp_bits_per_limb-1;
	mask=(1UL<<j);

	for(;(x[i]&mask)==0;j--,mask>>=1);

	lead=i*mp_bits_per_limb+j;	/* Ensured. */

	f->set(a,r,f);

	u=malloc(f->size*sizeof(mp_limb_t));

	for(;lead>0;lead--) {
		if (j--==0) {
			i--;
			j=mp_bits_per_limb-1;
			mask=(1UL<<j);
		} else {
			mask>>=1;
		}
	       	if (x[i]&mask) {
			f->mul(u,a,a,f);
			f->mul(a,u,r,f);
		} else {
			f->mul(u,a,a,f);
			f->set(a,u,f);
		}
	}
	free(u);
}
