#include <sys/time.h>
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
#include "structure.h"
#include "gmp-hacks.h"
#include "modulus_hacks.h"
#include "twisting_polynomials.h"
#include "e_polynomial.h"
#include "fft_on_matrices.h"
#include "field_def.h"
#include "field_prime.h"
#include "field_quad.h"
#include "field_usage.h"
#include "timer.h"

#ifndef NDEBUG
#include <stdio.h>
/* ... */
#endif

unsigned int ceil_log2(unsigned int d)
{
	unsigned int k;
	for(k=0;d;k++) d>>=1;
	return k;
}

/* fft engine. */

/* compute s1 <-- s1 + omega * s2
 *         s2 <-- s1 - omega * s2
 */
void butterfly(mp_limb_t * s1, mp_limb_t * s2, mp_limb_t * omega)
{
	mp_limb_t * tmp;

	tmp = FAST_ALLOC(l_size * sizeof(mp_limb_t));
	l_mul(tmp,s2,omega);
	l_sub(s2,s1,tmp);
	l_add(s1,s1,tmp);
	FAST_FREE(tmp);
}

/* compute s1 <--  s1 + s2
 *         s2 <-- (s1 - s2) * omega
 */
void butterfly_inverse(mp_limb_t * s1, mp_limb_t * s2, mp_limb_t * omega)
{
	mp_limb_t * tmp;

	tmp = FAST_ALLOC(l_size * sizeof(mp_limb_t));
	l_sub(tmp,s1,s2);
	l_add(s1,s1,s2);
	l_mul(s2,tmp,omega);
	FAST_FREE(tmp);
}

#if 0
/*
 * compute s1 <-- s1 + s2
 *         s2 <-- s1 - s2
 */
void basic_butterfly(mp_limb_t * s1, mp_limb_t * s2)
{
	mp_limb_t * tmp;

	tmp=FAST_ALLOC(l_size * sizeof(mp_limb_t));
	l_add(tmp,s1,s2);
	l_sub(s2,s1,s2);
	l_set(s1,tmp);
	FAST_FREE(tmp);
}
#endif

static int * phi_tab;
static int tabulated_order;

static void compute_phi(int d)
{
	int i,j,l;
	phi_tab[0]=0;
	for(j=1;j<(1<<d);phi_tab[j++]=123124);	/* on fout la merde... */
	for(l=1;l<=d;l++) {
		for(i=0;i<(1<<(l-1));i++) {
			phi_tab[(i<<(d-l+1))+(1<<(d-l))]=phi_tab[i<<(d-l+1)]+(1<<(l-1));
		}
	}
}

static unsigned int t_func(unsigned int n)
{
	unsigned int k;
	if (n==0) return 0;
	for(k=1;k<=n;k<<=1);
	return n^((k>>1)-1);
}

/* this is the core of the FFT. Is works for any l with the good
 * properties */

typedef void (*fft_capable_function)(mp_limb_t *,mp_size_t,mp_limb_t *,unsigned int);

void direct_fft_recursive(mp_limb_t * buf, mp_size_t gap,
			mp_limb_t * roots, unsigned int d)
{
	unsigned int i;
	mp_limb_t *p, *q, *r;
	if (d>1) {
		direct_fft_recursive(buf,gap<<1,roots,d-1);
		direct_fft_recursive(buf+gap,gap<<1,roots,d-1);
	}
	p=buf;
	q=buf+gap;
	r=roots;
	for(i=0;i<(1<<(d-1));i++) {
		butterfly(p,q,r);
		p+=gap<<1;
		q+=gap<<1;
		r+=l_size<<1;
	}
}

void inverse_fft_recursive(mp_limb_t * buf, mp_size_t gap,
			mp_limb_t * roots, unsigned int d)
{
	unsigned int i;
	mp_limb_t *p, *q;
	p=buf;
	q=buf+gap;
	for(i=0;i<(1<<(d-1));i++) {
		butterfly_inverse(p,q,roots+t_func(i<<1)*l_size);
		p+=gap<<1;
		q+=gap<<1;
	}
	if (d>1) {
		inverse_fft_recursive(buf,gap<<1,roots,d-1);
		inverse_fft_recursive(buf+gap,gap<<1,roots,d-1);
	}
}

void direct_fft_iterative(mp_limb_t * buf, mp_size_t gap,
			mp_limb_t * roots, unsigned int d)
{
	unsigned int l,offset,k;
	for(l=1;l<=d;l++)
	for(offset=0;offset<(1<<(d-l));offset++)
	for(k=0;k<(1<<(l-1));k++) {
		butterfly(buf+(offset+(k<<(d-l+1)))*gap,
			buf+(offset+(((k<<1)+1)<<(d-l))*gap),
			roots+(k<<1)*l_size);
	}
}

void inverse_fft_iterative(mp_limb_t * buf, mp_size_t gap,
			mp_limb_t * roots, unsigned int d)
{
	unsigned int l,offset,k;
	for(l=d;l>=1;l--)
	for(offset=0;offset<(1<<(d-l));offset++)
	for(k=0;k<(1<<(l-1));k++) {
		butterfly_inverse(buf+(offset+(k<<(d-l+1)))*gap,
			buf+(offset+(((k<<1)+1)<<(d-l))*gap),
			roots+t_func(k<<1)*l_size);
	}
}


fft_capable_function direct_fft,inverse_fft;

mp_limb_t * l_roots;

/* XXX : locality is certainly important, unless multiplications take a
 * huge time anyway. So it might be better to keep polynomials gathered
 * in a single area.
 */
struct dft_mb *	fft_ec_dft(struct e_coeff * ec, unsigned int order, double *tm)
{
	struct dft_mb * res;
	unsigned int i,j;
	int	k,l;
	struct timeval tv;

	timer_r(&tv,TIMER_SET);

	res=my_malloc(sizeof(struct dft_mb));
	if (res==NULL) return res;
	res->degree=ec->degree;
	res->order=order;
	if (!mbdft_alloc(order,res->p)) {
		free(res);
		return NULL;
	}

	mbdft_zero(order,res->p);

	/*
	 * We want to compute the order n DFT, n==res->order. ec->degree
	 * might actually be bigger than 1<<n, so we reduce e(X)
	 * transparently modulo X^{1<<n} in order to get the proper DFT
	 */
	for(i=0;	i<m_param;	i++)
	for(j=0;	j<bigdim;	j++) 
	for(k=0;	k<(1<<order);	k++) 
	for(l=k;	l<=ec->degree;	l+=(1<<order)) {
		/* The proper function to use is extend_scalars, but it
		 * doesn't do the addup thing... */
		k_addto(	mbdft_scal(order,res->p,i,j,k),
				mbmat_scal(mbpoly_coeff(ec->p,l),i,j));
	
	}
	for(i=0;i<m_param;i++) for(j=0;j<bigdim;j++) {
		(*direct_fft)(mbdft_poly(order,res->p,i,j),
			      l_size,
			      l_roots,
			      res->order);
		printf(".");
	}
	printf("\n");
	
	*tm = timer_r(&tv,TIMER_ASK);
	
	return res;
}


struct dft_bb *	fft_tp_dft(struct t_poly * tp, unsigned int order, double *tm)
{
	struct dft_bb * res;
	unsigned int i,j,jr;
	int k,l;
	struct timeval tv;

	timer_r(&tv,TIMER_SET);

	res=my_malloc(sizeof(struct dft_bb));
	if (res==NULL) return res;
	res->degree=tp->degree;
	res->order=order;
	if (!bbdft_alloc(order,res->p)) {
		free(res);
		return NULL;
	}

	bbdft_zero(order,res->p);

	/* Same thing as above */
	for(i=0;	i<bigdim;	i++)
	for(j=0;	j<bigdim;	j++) { jr=tp->clist[j];
	for(k=0;	k<(1<<order);	k++)
	for(l=k;	l<=tp->degnom[jr];l+=(1<<order)) {
			k_addto(bbdft_scal(order,res->p,i,j,k),
				bbmat_scal(bbpoly_coeff(tp->p,l),i,jr));
	}}

	for(i=0;i<bigdim;i++) for(j=0;j<bigdim;j++) {
		(*direct_fft)(bbdft_poly(order,res->p,i,j),
			      l_size,
			      l_roots,
			      res->order);
		printf(".");
	}
	printf("\n");

	*tm = timer_r(&tv,TIMER_ASK);

	return res;
}
	
struct dft_mb *	fft_mbb_conv_sp(struct dft_mb * p,
		struct dft_bb * q,
		unsigned int dg_kill,
		double *tm)
{
	struct dft_mb * res;
	unsigned int i,j,k,l;
	struct timeval tv;
	mp_limb_t * tmp;

	timer_r(&tv,TIMER_SET);

	res=my_malloc(sizeof(struct dft_mb));
	if (res==NULL) return res;
	res->degree=p->degree + q->degree - dg_kill;
	assert(p->order==q->order);
	res->order=p->order;
	if (!mbdft_alloc(res->order,res->p)) {
		free(res);
		return NULL;
	}

	mbdft_zero(res->order,res->p);
	tmp=my_malloc(l_size*sizeof(mp_limb_t));

	for(i=0;i<m_param;i++)
	for(j=0;j<bigdim;j++) {
		printf(".");
	for(l=0;l<bigdim;l++)
	for(k=0;k<(1<<(res->order));k++) {
		int kr;
		kr=phi_tab[(phi_tab[k]*(-dg_kill)) & ((1<<tabulated_order)-1)];
		l_mul(tmp,	mbdft_scal(res->order,p->p,i,l,k),
				bbdft_scal(res->order,q->p,l,j,k));

		l_addmul(mbdft_scal(res->order,res->p,i,j,k),
				tmp, l_roots+kr*l_size);
	}}
	printf("\n");

	free(tmp);

	*tm = timer_r(&tv,TIMER_ASK);
	
	return res;
}

struct dft_bb *	fft_bbb_conv(struct dft_bb * p, struct dft_bb * q, double *tm)
{
	struct dft_bb * res;
	unsigned int i,j,k,l;
	struct timeval tv;
	
	timer_r(&tv,TIMER_SET);

	res=my_malloc(sizeof(struct dft_bb));
	if (res==NULL) return res;
	res->degree=0;	/* we don't care, this is meaningless anyway */
	assert(p->order==q->order);
	res->order=p->order;
	if (!bbdft_alloc(res->order,res->p)) {
		free(res);
		return NULL;
	}

	bbdft_zero(res->order,res->p);

	for(i=0;i<bigdim;i++)
	for(j=0;j<bigdim;j++) {
		printf(".");
	for(l=0;l<bigdim;l++)
	for(k=0;k<(1<<(res->order));k++) {
		l_addmul(	bbdft_scal(res->order,res->p,i,j,k),
				bbdft_scal(res->order,p->p,i,l,k),
				bbdft_scal(res->order,q->p,l,j,k));
	}}
	printf("\n");
	
	*tm = timer_r(&tv,TIMER_ASK);
	
	return res;
}

void fft_mb_invdft(bw_mbpoly dest,
		struct dft_mb * p,
		unsigned int deg,
		double *tm)
{
	unsigned int i,j,k;
	mp_limb_t * one_over_n;
	struct timeval tv;

	timer_r(&tv,TIMER_SET);

	one_over_n=my_malloc(k_size*sizeof(mp_limb_t));
	k_set_int(one_over_n,1<<(p->order));
	k_inv(one_over_n,one_over_n);

	assert(deg<(1<<p->order));

	for(i=0;i<m_param;i++) for(j=0;j<bigdim;j++) {
		(*inverse_fft)(mbdft_poly(p->order,p->p,i,j),
			       l_size,
			       l_roots,
			       p->order);
		for(k=0;k<=deg;k++) {
			/*
			 * The proper function to use is
			 * restrict_scalars. But there's some more to do.
			 */
#ifndef HAS_NATIVE_FFT
			assert(k_is_zero(mbdft_scal(p->order,p->p,i,j,k)+k_size));
#endif
			k_mul(	mbmat_scal(mbpoly_coeff(dest,k),i,j),
				mbdft_scal(p->order,p->p,i,j,k),
				one_over_n);
		}
		/* Watch out: The other coefficients might very well be
		 * non-zero. That's not incorrect, but these are
		 * meaningless */
		printf(".");
	}
	printf("\n");
	free(one_over_n);

	*tm = timer_r(&tv,TIMER_ASK);
}

void fft_tp_invdft(struct t_poly * tp, struct dft_bb * p, double *tm)
{
	unsigned int i,j,k;
	mp_limb_t * one_over_n;
	struct timeval tv;

	timer_r(&tv,TIMER_SET);

	one_over_n=my_malloc(k_size*sizeof(mp_limb_t));
	k_set_int(one_over_n,1<<(p->order));
	k_inv(one_over_n,one_over_n);
	/*
	 * Les degrés de destinations sont déjà dans tp.
	 *
	 * Les colonnes ont déjà été réordonnées par tp_comp_alloc, donc
	 * la permutation associée à tp est triviale : clist[j]==j.
	 *
	 * Par ailleurs, tp_comp_alloc a déjà tout mis à zéro.
	 */

	for(i=0;i<bigdim;i++) for(j=0;j<bigdim;j++) {
		(*inverse_fft)(bbdft_poly(p->order,p->p,i,j),
			      l_size,
			      l_roots,
			      p->order);
		for(k=0;k<=tp->degnom[j];k++) {
			/*
			 * The proper function to use is
			 * restrict_scalars. But there's some more to do.
			 */
#ifndef HAS_NATIVE_FFT
			assert(k_is_zero(bbdft_scal(p->order,p->p,i,j,k)+k_size));
#endif
			k_mul(	bbmat_scal(bbpoly_coeff(tp->p,k),i,j),
				bbdft_scal(p->order,p->p,i,j,k),
				one_over_n);
		}
		printf(".");
	}
	printf("\n");
	free(one_over_n);

	*tm = timer_r(&tv,TIMER_ASK);
}

struct dft_bb * fft_bb_dft_init_one(unsigned int deg)
{
	struct dft_bb *res;
	unsigned int j,k;

	res=my_malloc(sizeof(struct dft_bb));
	if (res==NULL) return res;
	res->degree=0;
	res->order=ceil_log2(deg+1);
	if (!bbdft_alloc(res->order,res->p)) {
		free(res);
		return NULL;
	}

	bbdft_zero(res->order,res->p);

	for(k=0;k<(1<<(res->order));k++) for(j=0;j<bigdim;j++)  {
		l_set_one(bbdft_scal(res->order,res->p,j,j,k));
	}
	return res;
}

void dft_mb_free(struct dft_mb * p)
{
	mbdft_free(p->order,p->p);
	free(p);
}

void dft_bb_free(struct dft_bb * p)
{
	bbdft_free(p->order,p->p);
	free(p);
}


void prepare_fft_engine(unsigned int max_order, int enable_cplx)
{
	mp_limb_t * r;
	int i;
	mp_limb_t * ptr, *nptr;

	tabulated_order=max_order;
	phi_tab=my_malloc((1<<max_order)*sizeof(int));
	compute_phi(max_order);

	prepare_fields_for_fft(max_order,enable_cplx);
	r=fetch_primitive_root(max_order);

	l_roots=my_malloc((1<<max_order)*l_size*sizeof(mp_limb_t));
	ptr=l_roots;	/* phi_tab[0]=0. Always. */
	l_set_one(ptr);
	for(i=1;i<(1<<max_order);i++) {
		nptr=l_roots+phi_tab[i]*l_size;
		l_mul(nptr,ptr,r);
		ptr=nptr;
	}
	free(r);

	/* Just a matter of taste, since the difference in terms of
	 * efficiency is hardly noticeable */

	direct_fft	= &direct_fft_recursive;
	inverse_fft	= &inverse_fft_recursive;
}

void cleanup_fft_engine(void)
{
	free(phi_tab);
	free(l_roots);
}

