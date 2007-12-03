
#define	USE_GMP_LONGLONG

#ifdef	USE_GMP_LONGLONG
#include <gmp/gmp-impl.h>
#include <gmp/longlong.h>
#endif

#ifndef INCLUDING_ADDMUL_INLINES_H_
#error "This file is not intended to be included directly\n"
#endif

#if 0
#ifdef GROSBUG
#include <sys/time.h>
#include "threaded.h"
extern void fpscalar_n(FILE * f, mp_limb_t * p, mp_size_t n);
#endif
#endif

#ifdef	__cplusplus
extern "C" {
#endif

INLINING_PREFIX void reduce_p_minus_q(mp_limb_t * dst, mp_limb_t * p, mp_limb_t *q)
#define mn	bw_allocsize
{
	int res;
	mp_limb_t buf[3];

#if 0
#ifdef GROSBUG
	FILE * dbstr=tsv()->dbstr;

	fprintf(dbstr,"K!(");
	fpscalar_n(dbstr,p,mn+2);
	fprintf(dbstr,"-");
	fpscalar_n(dbstr,q,mn+2);
	fprintf(dbstr,") eq ");
#endif
#endif

	res=(mpn_cmp(p,q,mn+2));
	if (res<0) {
#ifndef NDEBUG
		mp_limb_t cy;
		/* Guess that q<xxmodulus_hbs */
		cy=mpn_sub_n(q,xxmodulus_hbs,q,mn+2);
		assert(cy==0);
		cy=mpn_add_n(dst, p, q, mn+2);
		assert(cy==0);
#else
		mpn_sub_n(q,xxmodulus_hbs,q,mn+2);
		mpn_add_n(dst, p, q, mn+2);
#endif
		/* We now do the full reduction, since it's better to do
		 * it all in one run rather than do with msb set first,
		 * and then with modulus_plain */
		mpn_tdiv_qr(buf,dst, 0, dst, mn+2, modulus_plain, mn);
	} else if (res) {
		mpn_sub_n(dst, p, q, mn+2);
		mpn_tdiv_qr(buf,dst, 0, dst, mn+2, modulus_plain, mn);
	} else /* if (res==0) */
		memset(dst, 0, mn*sizeof(mp_limb_t));
#if 0
#ifdef GROSBUG
	fpscalar_n(dbstr,dst,mn);
	fprintf(dbstr,";\n");
#endif
#endif
}
#undef mn

/* XXX Important: other ways to do this computation, somewhat less
 * XXX satisfactory, are in rpmq.c
 */

#if (defined(HARDCODE_PARAMS)&&(HARD_nbys==1))||!defined(BLOCKS_TOGETHER)
#define	STREAMLINE
#if defined(STREAMLINE) && defined(HARDCODE_PARAMS) && (HARD_nbys==1)
INLINING_PREFIX
void
addmultiply(bw_vector_block w, bw_vector_block v, stype32 x)
{
#if 1
	/* 3% faster because there's no jump ! */
	int t=(x<0);
	bw_lvblock_step_0n0(w,t*(bw_allocsize+2));
	x-=(t*x)<<1;
#else
	if (x<0) {
		bw_lvblock_step_0n0(w,(bw_allocsize+2));
		x=-x;
	}
#endif
#ifdef USE_GMP_LONGLONG
	/* This version is slightly faster, and I presume this can even
	 * be improved a tiny bit (some loads and stores are avoidable).
	 */
	{
#if 0
		mp_limb_t in, middle, out;
		int i;
		in=0;
		/* out is at most (2^32-1)^2 div 2^32=(2^32-2) */
		/* in is at most 2^32-1 */
		/* middle is at most 2^32-1 */

		/*
		(out,middle)+in+*w is at most:
		(2^64-2^33+1) + (2^32-1) + (2^32-1) = 2^64-1.
		therefore we have no carry popping out.
		*w+=(middle+in);
		in=out+ carry;
		*/
		for(i=0;i<bw_allocsize;i++) {
			umul_ppmm(out,middle,*v,x);
			add_ssaaaa(out,middle,out,middle,0,in);
			add_ssaaaa(out,*w,out,middle,0,*w);
			in=out;
			v++;
			w++;
		}
		add_ssaaaa(w[1],w[0],w[1],w[0],0,in);
#endif
		/* This assembly code is 10% better than the generic
		 * version above
		 */
		type32 ecx_out;
		__asm__ __volatile__(
			"\t pxor	%%mm0, %%mm0\n"
			"\t movd	%[val], %%mm7\n"
			"\t 0:\n"	/* remove to explicitly unroll */
			"\t movd	(%[src]), %%mm1\n"
			"\t leal	4(%[src]), %[src]\n"
			"\t movd	(%[dst]), %%mm2\n"
			"\t pmuludq	%%mm7, %%mm1\n"
			"\t paddq	%%mm1, %%mm2\n"
			"\t paddq	%%mm2, %%mm0\n"
			"\t subl	$1, %[cnt]\n"
			"\t movd	%%mm0, (%[dst])\n"
			"\t psrlq	$32, %%mm0\n"
			"\t leal	4(%[dst]), %[dst]\n"
			"\t jnz 0b\n"	/* remove to explicitly unroll */

			/* This is the part for explicit unrolling of the
			 * loop. For 3 limbs, it wins very marginally,
			 * therefore it's disabled.
			 
			"\t movd	(%[src]), %%mm1\n"
			"\t leal	4(%[src]), %[src]\n"
			"\t movd	(%[dst]), %%mm2\n"
			"\t pmuludq	%%mm7, %%mm1\n"
			"\t paddq	%%mm1, %%mm2\n"
			"\t paddq	%%mm2, %%mm0\n"
			"\t movd	%%mm0, (%[dst])\n"
			"\t psrlq	$32, %%mm0\n"
			"\t leal	4(%[dst]), %[dst]\n"

			"\t movd	(%[src]), %%mm1\n"
			"\t movd	(%[dst]), %%mm2\n"
			"\t pmuludq	%%mm7, %%mm1\n"
			"\t paddq	%%mm1, %%mm2\n"
			"\t paddq	%%mm2, %%mm0\n"
			"\t movd	%%mm0, (%[dst])\n"
			"\t psrlq	$32, %%mm0\n"
			"\t leal	4(%[dst]), %[dst]\n"
			*/

			"\t movd	(%[dst]), %%mm2\n"
			"\t paddq	%%mm2, %%mm0\n"
			"\t movd	%%mm0, (%[dst])\n"
			"\t psrlq	$32, %%mm0\n"
			"\t leal	4(%[dst]), %[dst]\n"

			"\t movd	(%[dst]), %%mm2\n"
			"\t paddq	%%mm2, %%mm0\n"
			"\t movd	%%mm0, (%[dst])\n"


			: "=a" (v), "=d" (w), "=c" (ecx_out)
			: [src] "a" (v), [dst] "d" (w),
			  [cnt] "c" (bw_allocsize), [val] "g" (x)
			: "%mm0", "%mm1", "%mm2", "%mm7", "memory");
	}
#else
	{
		mp_limb_t c;
		c=mpn_addmul_1(w,v,bw_allocsize,x);
		bw_lvblock_step_0n0(w,bw_allocsize);
		/* Since we have either HARDCODE_PARAMS && HARD_nbys==1
		 * || !BLOCKS_TOGETHER, we know that
		 * bw_lvblock_step_0p0(w) means w++
		 */
		mpn_add_1(w,w,2,c);
		/*
		c=mpn_add_1(w,w,1,c);
		bw_lvblock_step_0p0(w);
		c=mpn_add_1(w,w,1,c);
		assert(c==0);
		*/
	}
#endif
}
#else
INLINING_PREFIX
void
addmultiply(bw_vector_block w, bw_vector_block v, stype32 x)
{
#define mn	bw_allocsize
	mp_limb_t c;
#if (HARD_nbys!=1)
	int l;
#endif
#ifdef SWITCH
	mp_limb_t * buf UNUSED_VARIABLE;

	switch(x) {
#if 0
	    case  0:
		BUG(); break;
#endif
	    case -1:
		bw_lvblock_step_0n0(w,bw_allocsize+2);
	    case  1:
		c=mpn_add_n(w,w,v,bw_allocsize+2);
		assert(c==0);
		break;
#if 0
	    case -4:
		bw_lvblock_step_0n0(w,bw_allocsize+2);
	    case  4:
		buf=FAST_ALLOC((bw_allocsize+2)*sizeof(mp_limb_t));
		c=mpn_lshift(buf,v,bw_allocsize+2,2);
		assert(c==0);
		c=mpn_add_n(w,w,buf,bw_allocsize+2);
		assert(c==0);
		FAST_FREE(buf);
		break;
	    case -16:
		bw_lvblock_step_0n0(w,bw_allocsize+2);
	    case  16:
		buf=FAST_ALLOC((bw_allocsize+2)*sizeof(mp_limb_t));
		c=mpn_lshift(buf,v,bw_allocsize+2,4);
		assert(c==0);
		c=mpn_add_n(w,w,buf,bw_allocsize+2);
		assert(c==0);
		FAST_FREE(buf);
		break;
	    case -64:
		bw_lvblock_step_0n0(w,bw_allocsize+2);
	    case  64:
		buf=FAST_ALLOC((bw_allocsize+2)*sizeof(mp_limb_t));
		c=mpn_lshift(buf,v,bw_allocsize+2,6);
		assert(c==0);
		c=mpn_add_n(w,w,buf,bw_allocsize+2);
		assert(c==0);
		FAST_FREE(buf);
		break;
	    case -256:
		bw_lvblock_step_0n0(w,bw_allocsize+2);
	    case  256:
		buf=FAST_ALLOC((bw_allocsize+2)*sizeof(mp_limb_t));
		c=mpn_lshift(buf,v,bw_allocsize+2,8);
		assert(c==0);
		c=mpn_add_n(w,w,buf,bw_allocsize+2);
		assert(c==0);
		FAST_FREE(buf);
		break;
#endif
	    default:
#endif	/* SWITCH */
		if (x<0) {
			bw_lvblock_step_0n0(w,bw_allocsize+2);
#if (HARD_nbys!=1)
			for(l=0;l<nbys;l++) {
#endif
				c=mpn_addmul_1(w,v,bw_allocsize,-x);
				bw_lvblock_step_0n0(w,bw_allocsize);
				c=mpn_add_1(w,w,1,c);
				bw_lvblock_step_0p0(w);
				c=mpn_add_1(w,w,1,c);
#if defined(HARDCODE_PARAMS)&&(HARD_nbys==1)&&defined(NDEBUG)
				assert(c==0);
				bw_lvblock_step_0n0(w,bw_allocsize+3);
				bw_lvblock_step_0zp(w);
				bw_lvblock_step_00p(v);
#endif
#if (HARD_nbys!=1)
			}
#endif
		} else if (x) {
#if (HARD_nbys!=1)
			for(l=0;l<nbys;l++) {
#endif
				c=mpn_addmul_1(w,v,bw_allocsize,x);
				bw_lvblock_step_0n0(w,bw_allocsize);
				c=mpn_add_1(w,w,1,c);
				bw_lvblock_step_0p0(w);
				c=mpn_add_1(w,w,1,c);
#if defined(HARDCODE_PARAMS)&&(HARD_nbys==1)&&defined(NDEBUG)
				assert(c==0);
				bw_lvblock_step_0n0(w,bw_allocsize+3);
				bw_lvblock_step_0zp(w);
				bw_lvblock_step_00p(v);
#endif
#if (HARD_nbys!=1)
			}
#endif
		}
#ifdef SWITCH
	}
#endif
#undef mn
}
#endif
#else
INLINING_PREFIX
void
addmultiply(bw_vector_block w, bw_vector_block v, stype32 x)
{
	int k,l;
#ifndef HARDCODE_PARAMS
	mp_limb_t *lbuf;
	lbuf=my_malloc(nbys*sizeof(mp_limb_t));
#else
	mp_limb_t lbuf[nbys];
#endif

	if (x<0) {
		memset(lbuf,0,nbys*sizeof(mp_limb_t));
		bw_lvblock_step_0n0(w,bw_allocsize+2)
		for(k=0;k<bw_allocsize;k++) {
			/* optimize this loop with MMX */ 
			for(l=0;l<nbys;l++) {
				lbuf[l]=mpn_add_1(w,w,1,lbuf[l]);
				lbuf[l]+=mpn_addmul_1(w,v,1,-x);
				bw_lvblock_step_00p(w);
				bw_lvblock_step_00p(v);
			}
			bw_lvblock_step_0pz(w);
			bw_lvblock_step_0pz(v);
		}
		for(l=0;l<nbys;l++) {
			lbuf[l]=mpn_add_1(w,w,1,lbuf[l]);
			bw_lvblock_step_00p(v);
			bw_lvblock_step_00p(w);
		}
		bw_lvblock_step_0pz(w);
		bw_lvblock_step_0pz(v);
		for(l=0;l<nbys;l++) {
			lbuf[l]=mpn_add_1(w,w,1,lbuf[l]);
			assert(lbuf[l]==0);
			bw_lvblock_step_00p(v);
			bw_lvblock_step_00p(w);
		}
	} else if (x) {
		memset(lbuf,0,nbys*sizeof(mp_limb_t));
		bw_lvblock_step_0n0(w,bw_allocsize+2)
		for(k=0;k<bw_allocsize;k++) {
			/* XXX optimize this loop with MMX */ 
			for(l=0;l<nbys;l++) {
				lbuf[l]=mpn_add_1(w,w,1,lbuf[l]);
				lbuf[l]+=mpn_addmul_1(w,v,1,x);
				bw_lvblock_step_00p(w);
				bw_lvblock_step_00p(v);
			}
			bw_lvblock_step_0pz(w);
			bw_lvblock_step_0pz(v);
		}
		for(l=0;l<nbys;l++) {
			lbuf[l]=mpn_add_1(w,w,1,lbuf[l]);
			bw_lvblock_step_00p(v);
			bw_lvblock_step_00p(w);
		}
		bw_lvblock_step_0pz(w);
		bw_lvblock_step_0pz(v);
		for(l=0;l<nbys;l++) {
			lbuf[l]=mpn_add_1(w,w,1,lbuf[l]);
			assert(lbuf[l]==0);
			bw_lvblock_step_00p(v);
			bw_lvblock_step_00p(w);
		}
	}
#ifndef HARDCODE_PARMS
	free(lbuf);
#endif
}
#endif

#ifdef	__cplusplus
}
#endif
