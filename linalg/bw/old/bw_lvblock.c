#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "slave_params.h"
#include "params.h"
#include <assert.h>
#include "macros.h"
#include "types.h"
#include "endian.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "bw_lvblock.h"
#include "slave.h"
#include "variables.h"
#include "bw_lvblock_steps.h"
#include "threaded.h"
#include "addmul.h"

#if 0
OLD 	/* obsolete (was intended for mksol) */
OLD void _bw_lvblock_alloc_n(bw_vector_block * px, coord_t cols)
OLD {
OLD 	assert(!is_multithreaded);
OLD 	*px=mymalloc(cols*bw_longsize*nbys*sizeof(mp_limb_t));
OLD }
#endif

void _bw_lvblock_alloc(bw_vector_block * px)
{
	assert(!is_multithreaded);
	*px=mymalloc(ncols*bw_longsize*nbys*sizeof(mp_limb_t));
}

void _bw_lvblock_free(bw_vector_block * px)
{
	assert(!is_multithreaded);
	free(*px);
	*px=NULL;
}

void bw_lvblock_set_zero_separated(bw_vector_block x)
{
	bw_lvblock_step_n00(x,tsv()->i0);
	memset(x,0,(tsv()->i1-tsv()->i0)*bw_longsize*nbys*sizeof(mp_limb_t));
}

void bw_lvblock_set_zero_full(bw_vector_block x)
{
	bw_lvblock_step_n00(x,tsv()->i0);
	memset(x,0,ncols*bw_longsize*nbys*sizeof(mp_limb_t));
}

#if 0
OLD 	/* obsolete (was intended for mksol) */
OLD void bw_lvblock_set_zero_n(bw_vector_block x, coord_t cols)
OLD {
OLD 	assert(!is_multithreaded);
OLD 	memset(x,0,cols*bw_longsize*nbys*sizeof(mp_limb_t));
OLD }
#endif


#if (defined(HARDCODE_PARAMS)&&(HARD_nbys==1))||!defined(BLOCKS_TOGETHER)
size_t bw_lvblock_read_1(bw_vector_block x, FILE *f)
{
	assert(!is_multithreaded);
	return bw_scalar_read(x,f);
}
size_t bw_lvblock_read(bw_vector_block x, FILE *f)
{
	size_t res=0;
	coord_t i;
	assert(!is_multithreaded);
	for(i=0;i<ncols;i++) {
		res+=bw_scalar_read(x,f);
		bw_lvblock_step_p00(x);
	}
	return res;
}

size_t bw_lvblock_write(FILE *f, bw_vector_block x)
{
	size_t res=0;
	coord_t i;
	assert(!is_multithreaded || is_master_thread());
	for(i=0;i<ncols;i++) {
		res+=bw_scalar_write(f,x);
		bw_lvblock_step_p00(x);
	}
	return res;
}
#else
#error "This piece of code seems to mess up between be_allocsize and filesize"
size_t bw_lvblock_read_1(bw_vector_block x, FILE *f)
{
	size_t res=0;
	mp_size_t k;
	assert(!is_multithreaded);
	memset(x,0,bw_filesize*sizeof(mp_size_t));
	for(k=0;k<bw_filesize;k++) {
		res+=fread(x,sizeof(mp_limb_t),1,f);
		if (M_BIG_ENDIAN) mswap32(*x);
		bw_lvblock_step_0p0(x);
	}
	return res;
}
size_t bw_lvblock_read(bw_vector_block x, FILE *f)
{
	size_t res=0;
	coord_t i;
	assert(!is_multithreaded);
	for(i=0;i<ncols;i++) {
		for(k=0;k<bw_filesize;k++) {
			res+=fread(x,sizeof(mp_limb_t),1,f);
			if (M_BIG_ENDIAN) mswap32(*x);
			bw_lvblock_step_0p0(x);
		}
		memset(x+bw_allocsize,0,(bw_allocsize+4)*sizeof(mp_size_t));
		bw_lvblock_step_0n0(x,bw_allocsize+4);
		bw_lvblock_step_pz0(x);
	}
	return res;
}

size_t bw_lvblock_write(bw_vector_block x, FILE *f, coord_t cols)
{
	size_t res=0;
	coord_t i;
	mp_limb_t buf;
	assert(!is_multithreaded);
	for(i=0;i<ncols;i++) {
		for(k=0;k<bw_filesize;k++) {
			buf=*x;
			if (M_BIG_ENDIAN) mswap32(buf);
			res+=fwrite(&buf,sizeof(mp_limb_t),1,f);
			bw_lvblock_step_0p0(x);
		}
		bw_lvblock_step_0n0(x,bw_allocsize+4);
		bw_lvblock_step_pz0(x);
	}
	return res;
}
#endif

#if (defined(HARDCODE_PARAMS)&&(HARD_nbys==1))||!defined(BLOCKS_TOGETHER)
/* We can do it in place. */
void bw_lvblock_reduce_separated(bw_vector_block x)
{
	int j,l;
	bw_vector_block y;
	y=x+bw_allocsize+2;

	bw_lvblock_step_n00(x,tsv()->i0);
	bw_lvblock_step_n00(y,tsv()->i0);
	for(j=tsv()->i0;j<tsv()->i1;j++) {
		for(l=0;l<nbys;l++) {
			/* x[bw_allocsize+2]==y[0] is writable, OK for
			 * routine */
			reduce_p_minus_q(x,x,y);
			memset(x+bw_allocsize,0,
				(bw_allocsize+4)*sizeof(mp_size_t));
			bw_lvblock_step_00p(x);
			/* should be 0zp but gmp stepped */
			bw_lvblock_step_00p(y);	/* idem */
		}
		bw_lvblock_step_p0z(x);
		bw_lvblock_step_p0z(y);
	}
}
#else
/* must be reasonably fast. */
void bw_lvblock_reduce_separated(bw_vector_block x)
{
	int j,k,l;
#ifndef NDEBUG
	bw_vector_block y=x;
#endif
#ifdef HARDCODE_PARAMS
	mp_limb_t buf1[bw_longsize];
	mp_limb_t buf3[bw_allocsize+2];
#else
	mp_limb_t *buf1,*buf3;
	buf1=malloc(bw_longsize*sizeof(mp_limb_t));
	buf3=malloc(bw_allocsize*sizeof(mp_limb_t));
#endif
#define buf2 (buf1+bw_allocsize+2)
	bw_lvblock_step_n00(x,tsv()->i0);
	for(j=tsv()->i0;j<tsv()->i1;j++) {
		for(l=0;l<nbys;l++) {
			for(k=0;k<bw_longsize;k++) {
				buf1[k]=*x;
				bw_lvblock_step_0p0(x);
			}
			/* buf1[bw_allocsize+2] is writable, OK for
			 * routine */
			reduce_p_minus_q(buf1,buf1,buf1+bw_allocsize+2);
			bw_lvblock_step_0z0(x);
#if 0
			mpn_sub_n(	buf3,
					xxmodulus_hbs,
					buf2,bw_allocsize+2);
			/* no result since high bit is set. */
			*buf2=mpn_add_n(buf1,
					buf3,
					buf1,bw_allocsize+2);
			mpn_divrem(	buf1+bw_allocsize,0,
					buf1,bw_allocsize+3,
					modulus_hbs,bw_allocsize)
#endif
			for(k=0;k<bw_allocsize;k++) {
				*x=buf1[k];
				bw_lvblock_step_0p0(x);
			}
			for(;k<bw_longsize;k++) {
				*x=0;
				bw_lvblock_step_0p0(x);
			}
			bw_lvblock_step_0zp(x);
		}
		bw_lvblock_step_p0z(x);
	}
#ifndef NDEBUG
	bw_lvblock_step_n00(y,tsv()->i1-tsv()->i0);
	assert(y==x);
#endif
#undef buf2
#ifndef HARDCODE_PARAMS
	free(buf1);
	free(buf3);
#endif
}
#endif

/*
 * multiproduct has a hard job to do.
 *
 * Given vectors v[0]...v[k-1] and scalars s[0]..s[k-1], compute :
 *
 * w[i] += v[i] * s[i];
 *
 * The complicated part is that v[i]*s[i] is too large to fit into the
 * standard buffer, so that we have to use (almost) the full length which
 * is dedicated to p - q . The necessary reduction is performed
 * afterwards.
 *
 * w must be reduced on input
 * v must be reduced on input.
 *
 */

#if (defined(HARDCODE_PARAMS)&&(HARD_nbys==1))||!defined(BLOCKS_TOGETHER)

#ifndef NDEBUG
static int buffer_is_zero(mp_limb_t * buf, mp_size_t len)
{
	int i;
	for(i=0;i<len && buf[i]==0UL;i++);
	return i==len;
}
#endif

void
bw_multiproduct_separated(bw_vector_block w, bw_vector_block v, mp_limb_t * s)
{
	int j,l;
	mp_limb_t * sl;
	mp_limb_t * buf;
	mp_limb_t carry;

	buf=malloc((1+(bw_allocsize<<1))*sizeof(mp_limb_t));
	
	bw_lvblock_step_n00(w,tsv()->i0);
	bw_lvblock_step_n00(v,tsv()->i0);
	for(j=tsv()->i0;j<tsv()->i1;j++) {
		sl=s;
		for(l=0;l<nbys;l++) {
			mpn_mul_n(buf,v,sl,bw_allocsize);
			carry=mpn_add_n(buf,w,buf,bw_allocsize);
			carry=mpn_add_1(buf+bw_allocsize,buf+bw_allocsize,
					bw_allocsize,carry);
			assert(carry==0);	/* guaranteed, really */
			mpn_tdiv_qr(buf+bw_allocsize,w, 0,
					buf, bw_allocsize<<1,
					modulus_hbs, bw_allocsize);
			/*
			 * Does the same thing actually 
			mpn_divrem(w+bw_allocsize,0,
				w,bw_allocsize<<1,
				modulus_hbs,bw_allocsize);
				*/

			assert(buffer_is_zero(w+bw_allocsize,
						bw_longsize-bw_allocsize));
#if 0
			/* I guess that's plain useless */
			memset(w+bw_allocsize,0,
				(bw_longsize-bw_allocsize)*sizeof(mp_limb_t));
#endif

			bw_lvblock_step_00p(w);
			bw_lvblock_step_00p(v);
			sl+=bw_allocsize;
		}
		bw_lvblock_step_p0z(w);
		bw_lvblock_step_p0z(v);
	}

	free(buf);
}
#else	/* BLOCKS_TOGETHER */
void
bw_multiproduct(bw_vector_block w, bw_vector_block v, bw_vector_block s)
{
	int j,k,l;

	for(j=0;j<ncols;j++) {
#ifdef HARDCODE_PARAMS
	mp_limb_t buf1[bw_longsize];
	mp_limb_t buf3[bw_allocsize+2];
#else
	mp_limb_t *buf1,*buf3;
	buf1=malloc(bw_longsize*sizeof(mp_limb_t));
	buf3=malloc(bw_allocsize*sizeof(mp_limb_t));
#endif

	/* TODO */
	choke me of course...;

}
#endif	/* BLOCKS_TOGETHER */




#ifndef NDEBUG
void plvblock(bw_vector_block x)
{
	int j,k,l;
	bw_vector_block y;
	mp_limb_t *buf1;
	char * str;
	int cn, i;
	int maxcn;
	int z;
	
	buf1=malloc(bw_longsize*sizeof(mp_limb_t));
#define buf2 (buf1+bw_allocsize+2)
	maxcn = (bw_allocsize+2)*10+2;
	str = malloc(maxcn);
	y=x;
	printf("\n");
	for(j=0;j<ncols;j++) {
		if (j==0)
			printf("[");
		else
			printf(" ");
		for(l=0;l<nbys;l++) {
			for(k=0;k<bw_longsize;k++) {
				buf1[k]=*x;
				bw_lvblock_step_0p0(x);
			}
			bw_lvblock_step_0zp(x);
			cn=mpn_get_str((unsigned char *)str,10,buf1,bw_allocsize+2);
			if (cn > maxcn) {
				fprintf(stderr,"Overflow in debug !\n");
			}
			for(i=0;i<cn;i++) {
				str[i]+='0';
			}
			if (i==0) str[i++]='0';
			str[i]=0;
			z=strspn(str,"0");
			if (z==strlen(str)) z--; /* We know that length >0 */
			if (l>=1) printf(", "); else printf(" ");
			printf("(%s - ", str+z);
			cn=mpn_get_str((unsigned char *)str,10,buf2,bw_allocsize+2);
			if (cn > maxcn) {
				fprintf(stderr,"Overflow in debug !\n");
			}
			for(i=0;i<cn;i++) {
				str[i]+='0';
			}
			if (i==0) str[i++]='0';
			str[i]=0;
			z=strspn(str,"0");
			if (z==strlen(str)) z--; /* We know that length >0 */
			printf("%s)", str+z);
		
		}
		if (j!=ncols-1)
			printf(",\n");
		else
			printf(" ]\n");
		bw_lvblock_step_p0z(x);
	}
	bw_lvblock_step_z00(x);
	if (y!=x) {
		fprintf(stderr,"Bug in debug !\n");
	}
#undef buf2
	free(buf1);
	free(str);
}
#endif

/* XXX Watch out, this fails for p = 2, and also with BLOCKS_TOGETHER.
 * If p is not prime, I don't know. */
/*
int
bw_lvblock_is_zero(bw_vector_block v)
{
	int j,l;
	int res;
	int cmp;
	mp_limb_t * buf[3];
	mp_size_t bufs;

	bw_lvblock_reduce(v);

	buf[0] = malloc(bw_allocsize * sizeof(mp_limb_t));
	buf[1] = malloc(bw_allocsize * sizeof(mp_limb_t));
	buf[2] = malloc(bw_allocsize * sizeof(mp_limb_t));

	res = 1;

	for(j=0;res && j<ncols;j++) {
	    for(l=0;res && l<nbys;l++) {
		if (!bw_scalar_is_zero(v)) {
			cmp  = mpn_cmp(v,modulus_plain,bw_allocsize);
			if (cmp < 0) { res = 0 ; break; }
			if (cmp) {
				memcpy(buf[1],v,bw_allocsize);
				memcpy(buf[2],modulus_plain,bw_allocsize);
				bufs = mpn_gcd(buf[0],
						buf[1],bw_allocsize,
						buf[2],bw_allocsize);
				res &= (bufs != 1 || buf[0][0] == 1UL);
			}
		}
		bw_lvblock_step_00p(v);
	    }
	    bw_lvblock_step_p0z(v);
	}

	free(buf[0]);
	free(buf[1]);
	free(buf[2]);
	
	return res;
}
*/

/* Much cleaner anyway. No BLOCKS_TOGETHER though */
int
bw_lvblock_is_zero_separated(bw_vector_block v)
{
	int j,l;
	int res;
	mp_limb_t c;

	assert(!is_multithreaded);
	bw_lvblock_reduce_separated(v);
	res = 1;
	for(j = 0 ; j < ncols ; j++) {
		for(l = 0 ; l < nbys ; l++) {
			c=mpn_lshift(v,v,bw_allocsize+1,shift_amount);
			assert(c==0);
			mpn_divrem(v+bw_allocsize,0,
					v,bw_allocsize+1,
					modulus_hbs,bw_allocsize);
			c=mpn_rshift(v,v,bw_allocsize,shift_amount);
			assert(c==0);
			memset(v+bw_allocsize,0,bw_allocsize+4);
			res = res && bw_scalar_is_zero(v);
			bw_lvblock_step_00p(v);
		}
		bw_lvblock_step_p0z(v);
	}

	return res;
}

#if 0
OLD 	/* obsolete (was intended for mksol) */
OLD int bw_lvblock_add(bw_vector_block a, bw_vector_block b)
OLD {
OLD 	int j,l;
OLD 	int res;
OLD 	mp_limb_t c;
OLD #define n bw_allocsize
OLD 	assert(!is_multithreaded);
OLD 	res = 1;
OLD 	for(j = 0 ; j < ncols ; j++) {
OLD 		for(l = 0 ; l < nbys ; l++) {
OLD 			c=mpn_add_n(a,a,b,n+2);
OLD 			assert(c==0);
OLD 			c=mpn_add_n(a+n+2,a+n+2,b+n+2,n+2);
OLD 			assert(c==0);
OLD 			bw_lvblock_step_00p(a);
OLD 			bw_lvblock_step_00p(b);
OLD 		}
OLD 		bw_lvblock_step_p0z(a);
OLD 		bw_lvblock_step_p0z(b);
OLD 	}
OLD #undef n
OLD 	return res;
OLD }
#endif
