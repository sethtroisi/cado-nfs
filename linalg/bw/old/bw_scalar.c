#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"
#include <assert.h>
#include <alloca.h>
#define _PROTO(x) x
#include <gmp.h>
#include "macros.h"
#include "types.h"
#include "auxfuncs.h"
#include "bw_scalar.h"


/* XXX Watch out here ; we must *NOT* include any header likely to define
 * HARDCODE_PARAMS, because all parameters here are referred to as
 * plain-named (not computed_blah) */

static int init_done;
mpz_t modulus;
bw_scalar xxmodulus_hbs;
mp_limb_t * two_n_plus_2_mod_p;

mp_size_t bw_allocsize;
mp_size_t bw_longsize;
mp_size_t bw_filesize;
mp_size_t shift_amount;

bw_scalar zero;

int modulus_nfactors=0;
mpz_t * modulus_fact=NULL;
/* relies on GMP */

static int register_factor(mpz_t f)
{
	modulus_nfactors++;
	modulus_fact=realloc(modulus_fact,modulus_nfactors*sizeof(mpz_t));
	mpz_init(modulus_fact[modulus_nfactors-1]);
	mpz_set(modulus_fact[modulus_nfactors-1],f);

	return modulus_nfactors;
}

static void set_allocsize(void)
{
	for(bw_allocsize=abs(modulus->_mp_size);
			bw_allocsize>1 && modulus->_mp_d[bw_allocsize-1]==0;
			bw_allocsize--);

	bw_longsize=bw_allocsize*2+4;

	bw_filesize=bw_allocsize;
#if (mp_bits_per_limb<RAW_DATA_ALIGN)
	bw_filesize+=(-bw_filesize)&((RAW_DATA_ALIGN/mp_bits_per_limb)-1);
#endif
}

static void set_hbs(void)
{
	mp_limb_t t, *s;

	bw_scalar_alloc(zero,BW_SCALAR_SHORT);
	bw_scalar_set_zero(zero);

	bw_scalar_alloc(xxmodulus_hbs,BW_SCALAR_LONG);
	*xxmodulus_hbs=*xmodulus_hbs=0;
	bw_scalar_set_zero(modulus_hbs);	/* mostly useless... */
	memcpy(modulus_hbs,modulus->_mp_d,bw_allocsize*sizeof(mp_limb_t));
	for(shift_amount=0 ; ; shift_amount++) {
		t=modulus_hbs[bw_allocsize-1];
		if (mpn_lshift(modulus_hbs,modulus_hbs,bw_allocsize,1)!=0)
			break;
	}
	mpn_rshift(modulus_hbs,modulus_hbs,bw_allocsize,1); /* come back */
	modulus_hbs[bw_allocsize-1]=t;
	/* We must do so to recover the most significant bit (otherwise
	 * it is cleared. A nicer way would be to use mpzs... */

	two_n_plus_2_mod_p=malloc(bw_allocsize*sizeof(mp_limb_t));
	s=malloc((bw_allocsize+7)*sizeof(mp_limb_t));
	memset(s,0,(bw_allocsize+7)*sizeof(mp_limb_t));
	s[bw_allocsize+2]=1;
	mpn_tdiv_qr(s+bw_allocsize+3,two_n_plus_2_mod_p,0,
			s,bw_allocsize+3,
			modulus_plain,bw_allocsize);
	free(s);

}

void bw_shorten_long_scalar(bw_scalar s)
{
	mp_limb_t q[3];

	mpn_tdiv_qr(q,s, 0, s, bw_allocsize+2, modulus_hbs, bw_allocsize);
	s[bw_allocsize]=s[bw_allocsize+1]=0;
}

void bw_reduce_short_scalar(bw_scalar s)
{
	mp_limb_t q;

	mpn_tdiv_qr(&q,s, 0, s, bw_allocsize, modulus_plain, bw_allocsize);
}

#if 0
CRAP void bw_totally_reduce(bw_scalar s)
CRAP {
CRAP #define n bw_allocsize
CRAP 	mp_limb_t mask1, mask2;
CRAP 	int i,j;
CRAP #ifdef HARDCODE_PARAMS
CRAP 	mp_limb_t buf[n];
CRAP #else
CRAP 	mp_limb_t * buf;
CRAP 
CRAP 	buf=malloc(n*sizeof(mp_limb_t));
CRAP #endif
CRAP 	
CRAP 	mask1=(1UL<<(mp_bits_per_limb-1));
CRAP 	for(i=0;i<=shift_amount;i++,mask1>>=1) {
CRAP 		if (s[n-1] & mask1) {
CRAP 			/* There's one bit to kill */
CRAP 			mpn_lshift(buf,modulus_plain,n,shift_amount-i);
CRAP 			if (mpn_cmp(s,buf,n)>=0) {
CRAP 				mpn_sub_n(s,s,buf,n);
CRAP 				continue;
CRAP 			} else {
CRAP 				mask2=(1UL<<(mp_bits_per_limb-1-i-1));
CRAP 				for(j=i+1;j<shift_amount && (s[n-1] & mask1);j++) {
CRAP 					if (s[n-1] & mask2) {
CRAP 						mpn_lshift(buf,modulus_plain,n,shift_amount-j);
CRAP 						mpn_sub_n(s,s,buf,n);
CRAP 					}
CRAP 					mask2>>=1;
CRAP 				}
CRAP 			}
CRAP 		}
CRAP 	}
CRAP #undef n
CRAP }
#endif

static int bw_read_raw_modulus(const char * modulusname)
{
	if (!init_done) {
		mpz_init(modulus);
		init_done=1;
	}

	if (mpz_set_str(modulus, modulusname, 0)<0) {
		fprintf(stderr,"%s is not an acceptable modulus.\n",
				modulusname);
		return -1;
	}

	set_allocsize();
	register_factor(modulus);
	set_hbs();

	return 0;
}
	

int bw_read_modulus_info(const char * modulusname, int allow)
{
	FILE *f;
	mpz_t tempo;
	char idcheck[7];
	int i;

	f = fopen(modulusname,"r");
	if (f==NULL) {
		if (allow)
			return bw_read_raw_modulus(modulusname);
		else
			return -1;
	}

	if (!init_done) {
		mpz_init(modulus);
		init_done=1;
	}
	
	mpz_init(tempo);

	fgets(idcheck,7,f);
	if (strcmp(idcheck,"VALUE:")!=0) goto fail;
	if (mpz_inp_str(tempo,f,10)==0) goto fail;

	mpz_set(modulus,tempo);

	set_allocsize();

	fgets(idcheck,7,f);	/* finish the line */
	fgets(idcheck,7,f);
	if (strcmp(idcheck,"FACTS:")!=0) goto fail;

	for(;;)
	{
		if (mpz_inp_str(tempo,f,10)==0)
			break;
		register_factor(tempo);
	}

	mpz_set_ui(tempo,1);

	for(i=0;i<modulus_nfactors;i++)
		mpz_mul(tempo,tempo,modulus_fact[i]);

	if (mpz_cmp(tempo,modulus)!=0)
		goto fail;
	
	mpz_clear(tempo);

	set_hbs();
	
	return 0;
fail:
	mpz_clear(tempo);
	fclose(f);
	return -1;
}

mp_size_t _bw_scalar_alloc(bw_scalar * px, int mode)
{
	mp_size_t s;
	switch(mode) {
		case BW_SCALAR_SHORT: s=bw_allocsize; break;
		case BW_SCALAR_LONG:  s=bw_longsize; break;
		default: /*BW_SCALAR_MEDIUM*/ s=bw_allocsize+1; break;
	}

	*px=malloc(s*sizeof(mp_limb_t));

	return s;
}

void bw_scalar_set_one(bw_scalar x)
{
	bw_scalar_set_zero(x);
	(*x)++;
}

/* This function is superseded in spirit by the new set of functions
 * mpz_fits_<blah>_p in GMP 3 ; However, what we really wanna check here
 * is coherence with a type32 object, which is not exactly the same.
 *
 * XXX BTW, This breaks if unsigned long < 16bits.
 */
int bw_scalar_fits_word(bw_scalar a, type32 *p, int *s)
{
	mpz_t x,mmx,w;
	int res;
	int sgn;

	mpz_init(w);
	mpz_set_ui(w,(1<<16));
	mpz_mul(w,w,w);

	mpz_init(mmx);

	x->_mp_alloc=x->_mp_size=bw_allocsize;
	x->_mp_d=malloc(x->_mp_alloc*sizeof(mp_limb_t));
	memcpy(x->_mp_d,a,x->_mp_alloc*sizeof(mp_limb_t));
	for(;x->_mp_d[x->_mp_size-1]==0 && x->_mp_size>1;x->_mp_size--);

	mpz_sub(mmx,modulus,x);

	if (mpz_cmp(x,mmx)<0) {
		sgn=1;
		res=mpz_cmp(x,w)<0;
	} else {
		sgn=-1;
		res=mpz_cmp(mmx,w)<0;
	}

	if (s!=NULL && p!=NULL) switch (sgn) {
		case  1:
			*p=mpz_get_ui(x);
			*s= 1;
			break;
		default: /* case -1: */
			*p=mpz_get_ui(mmx);
			*s=-1;
			break;
	}
		
	mpz_clear(x);
	mpz_clear(mmx);
	mpz_clear(w);

	return res;
}

void bw_scalar_set_random(bw_scalar x)
{
	mpz_t y;
	mpz_init(y);
	mpz_random(y,bw_allocsize);
	memcpy(x,y->_mp_d,bw_allocsize*sizeof(mp_limb_t));
	mpz_clear(y);
}

/* vectors */
mp_size_t _bw_vector_alloc(bw_vector *pv, coord_t n, int mode)
{
	mp_size_t s;

	switch(mode) {
		case BW_SCALAR_SHORT: s=bw_allocsize; break;
		case BW_SCALAR_LONG:  s=bw_longsize; break;
		default: /*BW_SCALAR_MEDIUM*/ s=bw_allocsize+1; break;
	}

	*pv=malloc(s*n*sizeof(mp_limb_t));

	return s*n;
}

void bw_vector_set_random(bw_vector v, int n)
{
	bw_lemming l;
	int i;
	for(i=0,bw_lemming_reset(l,v);i<n;bw_lemming_step(l,1),i++)
		bw_scalar_set_random(l);
	bw_lemming_clear(l);
}

#ifndef NDEBUG
void pscalar(mp_limb_t * p)
{
	char * foo;
	mp_limb_t * bar;
	int i,l;

	bar=malloc((bw_allocsize+1)*sizeof(mp_limb_t));
	foo=malloc((bw_allocsize+1)*mp_bits_per_limb);
	/* overly exaggerated... */
	memcpy(bar,p,bw_allocsize*sizeof(mp_limb_t));
	l=mpn_get_str((unsigned char*)foo,10,bar,bw_allocsize);
	for(i=0;i<l;i++) {
		if (foo[i]<10)
			foo[i]+='0';
		else if (foo[i]<36)
			foo[i]+='A';
		else
			foo[i]='!';
		printf("%c",foo[i]);
	}
	printf("\n");
	free(foo);
	free(bar);
}
void fpscalar_n(FILE *f, mp_limb_t * p, mp_size_t n)
{
	char * foo;
	mp_limb_t * bar;
	int i,l;

	bar=malloc((n+1)*sizeof(mp_limb_t));
	foo=malloc((n+1)*mp_bits_per_limb);
	/* overly exaggerated... */
	memcpy(bar,p,n*sizeof(mp_limb_t));
	l=mpn_get_str((unsigned char*)foo,10,bar,n);
	for(i=0;i<l;i++) {
		if (foo[i]<10)
			foo[i]+='0';
		else if (foo[i]<36)
			foo[i]+='A';
		else
			foo[i]='!';
		fprintf(f,"%c",foo[i]);
	}
	free(foo);
	free(bar);
}
#endif
/* vim:set sw=8: */
