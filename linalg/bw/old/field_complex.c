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
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "variables.h"
#include "gmp-hacks.h"
#include "modulus_hacks.h"
#include "field_def.h"
#include "field_common.h"

/* almost completely identical to field_quad. There's jsut a little
 * speedup to be hoped for the multiplication, if only this were
 * programmed with more care... */
/* First forward some declarations. */
#define ELT mp_limb_t *
void complex_field_add(ELT, ELT, ELT, struct field *);
void complex_field_neg(ELT, ELT, struct field *);
void complex_field_sub(ELT, ELT, ELT, struct field *);
void complex_field_addto(ELT, ELT, struct field *);
void complex_field_mul(ELT, ELT, ELT, struct field *);
void complex_field_addmul(ELT, ELT, ELT, struct field *);
void complex_field_inv(ELT, ELT, struct field *);
void complex_field_set_int(ELT, int, struct field *);
void complex_field_set_mpn(ELT, ELT, mp_size_t, struct field *);
void complex_field_set_mpz(ELT, mpz_t p, struct field *);
int complex_field_is_int(ELT, int, struct field *);
void complex_field_reduce(ELT, struct field *);
void complex_field_conjugate(ELT, ELT, struct field *);
#undef ELT

/* Fortunately, 2^(kw) is never a prime number for 32|w and k>0.
 * Therefore we can assume that the same storage space is enough to store
 * elements of k as well as the modulus.
 */
void complex_field_add(mp_limb_t * a,
		mp_limb_t * b,
		mp_limb_t * c,
		struct field * f)
{
	/* These first routine work for any vs over the base field */
	int i;
	for(i=0;i<f->size;i+=f->base_field->size) {
		f->base_field->add(a+i,b+i,c+i,f->base_field);
	}
}
void complex_field_neg(mp_limb_t * a,
		mp_limb_t * b,
		struct field * f)
{
	int i;
	for(i=0;i<f->size;i+=f->base_field->size) {
		f->base_field->neg(a+i,b+i,f->base_field);
	}
}

void complex_field_sub(mp_limb_t * a,
		mp_limb_t * b,
		mp_limb_t * c,
		struct field * f)
{
	int i;
	for(i=0;i<f->size;i+=f->base_field->size) {
		f->base_field->sub(a+i,b+i,c+i,f->base_field);
	}
}

void complex_field_addto(mp_limb_t * a,
		mp_limb_t * b,
		struct field * f)
{
	complex_field_add(a,a,b,f);
}

/* overlap allowed */
void complex_field_mul(mp_limb_t * a,
		mp_limb_t * b,
		mp_limb_t * c,
		struct field * f)
{
	mp_limb_t *t;
	mp_size_t pl;
	pl=f->base_field->size;
	t=FAST_ALLOC(f->size*sizeof(mp_limb_t));
	/* Four multiplications... */
	f->base_field->mul(t,b,c,f->base_field);
	f->base_field->mul(t+pl,b+pl,c,f->base_field);
	f->base_field->addmul(t+pl,b,c+pl,f->base_field);
	/* funny trick here to allow overlaps */
	f->base_field->set(a,t,f->base_field);
	f->base_field->mul(t,b+pl,c+pl,f->base_field);
	f->base_field->set(a+pl,t+pl,f->base_field);
	f->base_field->neg(t+pl,t,f->base_field);
	f->base_field->addto(a,t+pl,f->base_field);
	FAST_FREE(t);
}

void complex_field_addmul(mp_limb_t * a,
		mp_limb_t * b,
		mp_limb_t * c,
		struct field * f)
{
	mp_limb_t * t;
	t=FAST_ALLOC(f->size*sizeof(mp_limb_t));
	complex_field_mul(t,b,c,f);
	complex_field_add(a,a,t,f);
	FAST_FREE(t);
}

void complex_field_inv(mp_limb_t * a,
		mp_limb_t * b,
		struct field * f)
{
	mp_limb_t * t, *u;
	mp_size_t pl;
	pl=f->base_field->size;
	t=FAST_ALLOC(f->size*sizeof(mp_limb_t));
	u=FAST_ALLOC(f->size*sizeof(mp_limb_t));
	f->base_field->mul(t,b,b,f->base_field);
	f->base_field->mul(t+pl,b+pl,b+pl,f->base_field);
	f->base_field->addto(t,t+pl,f->base_field);
	complex_field_conjugate(u,b,f);
	f->base_field->inv(t+pl,t,f->base_field);
	f->base_field->mul(a,u,t+pl,f->base_field);
	f->base_field->mul(a+pl,u+pl,t+pl,f->base_field);
	FAST_FREE(t);
	FAST_FREE(u);
}

void complex_field_set_int(mp_limb_t * a,
		int x,
		struct field * f)
{
	f->base_field->set_int(a,x,f->base_field);
	f->base_field->set_int(a+f->base_field->size,0,f->base_field);
}

void complex_field_set_mpn(mp_limb_t * a,
		mp_limb_t * p,
		mp_size_t n,
		struct field * f)
{
	f->base_field->set_mpn(a,p,n,f->base_field);
	f->base_field->set_int(a+f->base_field->size,0,f->base_field);
}

void complex_field_set_mpz(mp_limb_t * a,
		mpz_t p,
		struct field * f)
{
	f->base_field->set_mpz(a,p,f->base_field);
	f->base_field->set_int(a+f->base_field->size,0,f->base_field);
}

int complex_field_is_int(mp_limb_t * a,
		int x,
		struct field * f)
{
	return	f->base_field->is_int(a+f->base_field->size,0,f->base_field) && 
		f->base_field->is_int(a,x,f->base_field);
}

void complex_field_reduce(mp_limb_t * a, struct field * f)
{
	f->base_field->reduce(a,f->base_field);
	f->base_field->reduce(a+f->base_field->size,f->base_field);
}

void complex_field_conjugate(mp_limb_t * a, mp_limb_t * b, struct field * f)
{
	memcpy(a,b,f->size*sizeof(mp_limb_t));
	f->base_field->neg(a+f->base_field->size,a+f->base_field->size,f->base_field);

}

void complex_field_print(mp_limb_t * p, struct field * f)
{
	int z=1;
	printf("(");
	if (!f->base_field->is_int(p,0,f->base_field)) {
		z=0;
		f->base_field->print(p,f->base_field);
	}
	if (!f->base_field->is_int(p+f->base_field->size,0,f->base_field)) {
		if (!z) printf("+");
		z=0;
		f->base_field->print(p+f->base_field->size,f->base_field);
		printf("*i");
	}
	if (z) {
		printf("0");
	}
	printf(")");
}

struct field * new_complex_field(struct field * base)
{
	struct field * res;

	res=malloc(sizeof(struct field));
	res->size=2*base->size;
	res->degree=2;
	res->minpoly	= NULL;
	res->base_field	= base;
	res->add	= complex_field_add;
	res->sub	= complex_field_sub;
	res->neg	= complex_field_neg;
	res->addto	= complex_field_addto;
	res->mul	= complex_field_mul;
	res->addmul	= complex_field_addmul;
	res->inv	= complex_field_inv;
	res->set_int	= complex_field_set_int;
	res->set_mpn	= complex_field_set_mpn;
	res->set_mpz	= complex_field_set_mpz;
	res->is_int	= complex_field_is_int;
	res->exp	= generic_field_exp;
	res->eq		= generic_field_eq;
	res->set	= generic_field_set;
	res->reduce	= complex_field_reduce;
	res->conjugate	= complex_field_conjugate;
	res->print	= complex_field_print;
	res->destroy	= generic_field_destroy;
	res->modulus	= NULL;		/* Does not apply */
	return res;
}
