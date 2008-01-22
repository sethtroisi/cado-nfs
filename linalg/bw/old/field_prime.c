#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <gmp.h>
#include "lingen_params.h"
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


/* Fortunately, 2^(kw) is never a prime number for 32|w and k>0.
 * Therefore we can assume that the same storage space is enough to store
 * elements of k as well as the modulus.
 */

/* First forward some declarations. */
#define ELT mp_limb_t *
void prime_field_add(ELT, ELT, ELT, struct field *);
void prime_field_neg(ELT, ELT, struct field *);
void prime_field_sub(ELT, ELT, ELT, struct field *);
void prime_field_addto(ELT, ELT, struct field *);
void prime_field_mul(ELT, ELT, ELT, struct field *);
void prime_field_addmul(ELT, ELT, ELT, struct field *);
void prime_field_inv(ELT, ELT, struct field *);
void prime_field_set_int(ELT, int, struct field *);
void prime_field_set_mpn(ELT, ELT, mp_size_t, struct field *);
void prime_field_set_mpz(ELT, mpz_t p, struct field *);
int prime_field_is_int(ELT, int, struct field *);
void prime_field_reduce(ELT, struct field *);
void prime_field_conjugate(ELT, ELT, struct field *);
#undef ELT

void prime_field_add(mp_limb_t * a,
		mp_limb_t * b,
		mp_limb_t * c,
		struct field * f)
{
	mp_limb_t * t;
	t=FAST_ALLOC((2+f->size)*sizeof(mp_limb_t));
	t[f->size]=mpn_add_n(t,b,c,f->size);
	prime_field_set_mpn(a,t,1+f->size,f);
	FAST_FREE(t);
}
void prime_field_neg(mp_limb_t * a,
		mp_limb_t * b,
		struct field * f)
{
	prime_field_reduce(b,f);
	mpn_sub_n(a,f->modulus,b,f->size);
}

void prime_field_sub(mp_limb_t * a,
		mp_limb_t * b,
		mp_limb_t * c,
		struct field * f)
{
	mp_limb_t * t;
	t=FAST_ALLOC(f->size*sizeof(mp_limb_t));
	prime_field_neg(t,c,f);
	prime_field_add(a,b,t,f);
	FAST_FREE(t);
}

void prime_field_addto(mp_limb_t * a,
		mp_limb_t * b,
		struct field * f)
{
	prime_field_add(a,a,b,f);
}

void prime_field_mul(mp_limb_t * a,
		mp_limb_t * b,
		mp_limb_t * c,
		struct field * f)
{
	mp_limb_t * t;
	t=FAST_ALLOC((1+(f->size<<1))*sizeof(mp_limb_t));
	mpn_mul_n(t,b,c,f->size);
	prime_field_set_mpn(a,t,f->size<<1,f);
	FAST_FREE(t);
}

void prime_field_addmul(mp_limb_t * a,
		mp_limb_t * b,
		mp_limb_t * c,
		struct field * f)
{
	mp_limb_t * t;
	t=FAST_ALLOC((1+(f->size<<1))*sizeof(mp_limb_t));
	mpn_mul_n(t,b,c,f->size);
	prime_field_set_mpn(t,t,f->size<<1,f);	/* This is allowed */
	prime_field_add(a,a,t,f);
	FAST_FREE(t);
}

void prime_field_inv(mp_limb_t * a,
		mp_limb_t * b,
		struct field * f)
{
	mpz_t az,bz,cz;

	mpz_init(az);
	MPZ_INIT_SET_MPN(bz,b,f->size);
	MPZ_INIT_SET_MPN(cz,f->modulus,f->size);

	if (mpz_invert(az,bz,cz) == 0) {
		mpz_gcd(az,bz,modulus);
		mpz_out_str(stdout,10,az);
		printf(" is a factor of the modulus\n");
		exit(191);
	}
	mpz_clear(bz);
	mpz_clear(cz);
	prime_field_set_mpz(a,az,f);
	mpz_clear(az);
}

void prime_field_set_int(mp_limb_t * a,
		int x,
		struct field * f)
{
	if (x>=0) {
		memset(a+1,0,(f->size-1)*sizeof(mp_limb_t));
		a[0]=x;
#if 0
		/* Hmmm, who cares ? Didn't I say _lazy_ reduction ? */
		if (f->size==1 && x >= f->modulus[0]) {
			a[0]=x % f->modulus[0];
		}
#endif
	} else {	/* then x<=0 */
		x=-x;
		if (f->size==1 && x >= f->modulus[0]) {
			x=x % f->modulus[0];
		}
		a[0]=x;
		prime_field_sub(a,f->modulus,a,f);
	}
}

void prime_field_set_mpn(mp_limb_t * a,
		mp_limb_t * p,
		mp_size_t n,
		struct field * f)
{
	mp_limb_t * t;
	if (n < f->size) {
		memcpy(a,p,n*sizeof(mp_limb_t));
		memset(a+n,0,(f->size-n)*sizeof(mp_limb_t));
	} else {
		/* We use a temporary in order to avoid a problem.
		 * Otherwise, p would have to have enough size to store
		 * __n+1__ limbs... This is a requirement due to tdiv_qr
		 * and the fact that the modulus needs not have highest
		 * bit set.
		 */
		t=FAST_ALLOC((n+1)*sizeof(mp_limb_t));
		mpn_tdiv_qr(t+f->size, t,0, p,n, f->modulus,f->size);
		/* memset(p+f->size,0,(n+1-f->size)*sizeof(mp_limb_t));
		 */
		memcpy(a,t,f->size*sizeof(mp_limb_t));
		FAST_FREE(t);
	}
}

void prime_field_set_mpz(mp_limb_t * a,
		mpz_t p,
		struct field * f)
{
	mp_limb_t * t;
	if (p[0]._mp_size==0)
		prime_field_set_int(a,0,f);
	else if (p[0]._mp_size>0) {
		t=FAST_ALLOC((1+p[0]._mp_size)*sizeof(mp_limb_t));
		memcpy(t,p[0]._mp_d,p[0]._mp_size*sizeof(mp_limb_t));
		prime_field_set_mpn(a,t,p[0]._mp_size,f);
		FAST_FREE(t);
	} else {
		t=FAST_ALLOC((1-p[0]._mp_size)*sizeof(mp_limb_t));
		memcpy(t,p[0]._mp_d,-p[0]._mp_size*sizeof(mp_limb_t));
		prime_field_set_mpn(a,t,-p[0]._mp_size,f);
		mpn_sub_n(a,f->modulus,a,f->size);
		FAST_FREE(t);
	}
}

int prime_field_is_int(mp_limb_t * a,
		int x,
		struct field * f)
{
	int i;
	f->reduce(a,f->base_field);
	for(i=f->size;i<f->size;i++)
		if (a[i]!=0UL) return 0;
	if (x>=0) {
		for(i=1;i<f->size;i++)
			if (a[i]!=0UL) return 0;
		if (f->size==1 && x >= f->modulus[0]) {
			x=x % f->modulus[0];
		}
		return (a[0]==x);
	} else {
		mp_limb_t * t;
		int res;
		x=-x;
		if (f->size==1 && x >= f->modulus[0]) {
			x=x % f->modulus[0];
		}
		t=FAST_ALLOC(f->size * sizeof(mp_limb_t));
		t[0]=x;
		mpn_add_n(t,t,a,f->size);
		res=(memcmp(t,f->modulus,f->size)==0);
		FAST_FREE(t);
		return res;
	}
}

void prime_field_reduce(mp_limb_t * a, struct field * f)
{
	mp_limb_t quot;
	mpn_tdiv_qr(&quot,a,0,a,f->size,f->modulus,f->size);
}

void prime_field_conjugate(mp_limb_t * a, mp_limb_t * b, struct field * f)
{
	generic_field_set(a,b,f);
}

void prime_field_print(mp_limb_t * p, struct field * f)
{
	mpz_t a;
	mpz_init(a);
	MPZ_SET_MPN(a,p,f->size);
	gmp_printf("%2Zd", a);
	/* mpz_out_str(stdout,10,a); */
	mpz_clear(a);
}

struct field * new_prime_field(mp_limb_t * p, mp_size_t n)
{
	struct field * res;

	res=malloc(sizeof(struct field));
	res->size=n;
	res->modulus=malloc(n*sizeof(mp_limb_t));
	memcpy(res->modulus,p,n*sizeof(mp_limb_t));
	res->degree=1;
	res->base_field	= res;
	res->add	= prime_field_add;
	res->sub	= prime_field_sub;
	res->neg	= prime_field_neg;
	res->addto	= prime_field_addto;
	res->mul	= prime_field_mul;
	res->addmul	= prime_field_addmul;
	res->inv	= prime_field_inv;
	res->set_int	= prime_field_set_int;
	res->set_mpn	= prime_field_set_mpn;
	res->set_mpz	= prime_field_set_mpz;
	res->is_int	= prime_field_is_int;
	res->exp	= generic_field_exp;
	res->eq		= generic_field_eq;
	res->set	= generic_field_set;
	res->reduce	= prime_field_reduce;
	res->conjugate	= prime_field_conjugate;
	res->print	= prime_field_print;
	res->destroy	= generic_field_destroy;
	res->minpoly	= NULL;		/* Does not apply */

	return res;
}

