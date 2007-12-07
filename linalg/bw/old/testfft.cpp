#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "master_params.h"
#include "params.h"
#include <assert.h>
#include "types.h"
#include "old-endian.h"
#include "macros.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "filenames.h"
#include "options.h"
#include "tagfile.h"
#include "timer.h"
#include "variables.h"
#include "structure.h"
#include "gmp-hacks.h"
#include "modulus_hacks.h"
#include "right_action.h"
#include "e_polynomial.h"
#include "twisting_polynomials.h"
#include "fft_on_matrices.hpp"
#include "field_def.h"
#include "field_prime.h"
#include "field_quad.h"
#include "field_usage.h"
#include "fft_core.h"

#define USE_MACROS

#define	ELT	mp_limb_t *
extern "C" {
extern void complex_field_addmul(ELT, ELT, ELT, struct field *);
extern void quad_field_addmul(ELT, ELT, ELT, struct field *);
extern void prime_field_mul(ELT, ELT, ELT, struct field *);
extern void prime_field_addmul(ELT, ELT, ELT, struct field *);
extern void prime_field_set_mpn(ELT, mp_limb_t *, mp_size_t, struct field *);
}

extern mp_limb_t * l_roots;

#define fdeg	bw_allocsize

static int	enable_cplx;
static int	order;
mp_limb_t * one_over_n;

struct multi_array {
	mp_limb_t * raw;
	int	dim[3];
	int	per[3];	/* normal storage yields {0,1,2} */
	int	ext;
	size_t	nlimbs;
};

struct multi_array * declare_multi_array(const int dim[3], const int per[3], int ext)
{
	struct multi_array * res;
	res         = (multi_array *) malloc(sizeof(struct multi_array));
	res->nlimbs = dim[0]*dim[1]*dim[2]*ext;
	res->raw    = (mp_limb_t *) malloc(res->nlimbs*sizeof(mp_limb_t));
	memcpy(res->dim,dim,3*sizeof(int));
	memcpy(res->per,per,3*sizeof(int));
	res->ext    = ext;
	return res;
}

void delete_multi_array(struct multi_array * blah)
{
	free(blah->raw);
	free(blah);
}

int multi_array_offset_tab(struct multi_array * x, const int tab[3])
{
	/* La formule magique à la mors-moi-le-noeud... */

	return x->ext*(tab[x->per[2]]+x->dim[x->per[2]]*(
			tab[x->per[1]]+x->dim[x->per[1]]*(
				tab[x->per[0]])));

}

int multi_array_offset_ijk(struct multi_array * x, int i, int j, int k)
{
	union {
		int	t[3];
		struct {
			int i;
			int j;
			int k;
		} s;
	} idx;

	idx.s.i=i;
	idx.s.j=j;
	idx.s.k=k;
	
	return multi_array_offset_tab(x,idx.t);
}

void pr_ma(struct multi_array * x)
{
	int i,j,k;
	mp_limb_t * p;
	struct field * bf;
	char ind;
	if (x->ext==fdeg) {
		bf=field_k;
		ind='X';
	} else if (x->ext==(fdeg<<1)) {
		bf=field_l;
		ind='Y';
	} else {
		fprintf(stderr,"Bad base field\n");
		return;
	}
		
	for(i=0;i<x->dim[0];i++) {
		printf("%c",i?' ':'[');
		for(j=0;j<x->dim[1];j++) {
			int z=1;
			printf("%c ",j?',':'[');
			for(k=0;k<x->dim[2];k++) {
				p=x->raw+multi_array_offset_ijk(x,i,j,k);
				if (bf->is_int(p,0,bf))
					continue;
				if (!z) printf("+");
				bf->print(p,bf);
				switch (k) {
				    case 0: break;
				    case 1: printf("*%c",ind); break;
				    default: printf("*%c^%d",ind,k); break;
				}
				z=0;
			}
			if (z) printf("0");
		}
		printf(" ]%c\n",(i==(x->dim[0]-1))?']':',');
	}
}

#if 0
/* These are just the dimensions. One of the dimensions is allowed to be
 * -1, in which case n can be anything. Otherwise the allowed values for n
 * are 0...d-1, and -d, where d is the dimension
 */

struct virtual_ma {
	int depth;
	int *dim;
};

struct virtual_ma * new_vma(int n)
{
	struct virtual_ma * res;
	res=malloc(sizeof(struct virtual_ma));
	res->depth=n;
	res->dim=malloc(n*sizeof(int));
	return res;
}

void delete_vma(struct virtual_ma * vma)
{
	free(vma->dim);
	free(vma);
}

int vma_getoffset(struct virtual_ma * vma, int i, int n)
{
	int res;

	if (vma->dim[i]==-1)
		return 0;

	assert(n==1 || n==-vma->dim[i] || (n>=0 && n<vma->dim[i]));

	res=n;
	for(i++;i<vma->depth;i++) res*=abs(vma->dim[i]);
	return res;
}
#endif

int dim_p[3];	/* Always specified as lines, columns, coeffs */
int dim_q[3];
int dim_r[3];

int dim_dp[3];
int dim_dq[3];
int dim_dr[3];

#define	O_ROWS		0
#define	O_COLUMNS	1
#define O_POLYS		2
#define O_IDX_R		2
#define O_INNER		3
#define O_IDX_P		4
#define O_IDX_Q		5

/* {O_ROWS,O_COLUMNS,O_POLYS} means: coeffs in polys are gathered, rows
 * are cast away. resp: small/large distance between consecutive elts.
 *
 * IOW: This is the ordering, by increasing order of locality */

int perm_p[3]=	{O_ROWS, O_COLUMNS, O_POLYS};
int perm_q[3]=	{O_ROWS, O_COLUMNS, O_POLYS};
int perm_r[3]=	{O_ROWS, O_COLUMNS, O_POLYS};

int perm_dp[3]=	{O_ROWS, O_COLUMNS, O_POLYS};
int perm_dq[3]=	{O_ROWS, O_COLUMNS, O_POLYS};
int perm_dr[3]=	{O_ROWS, O_COLUMNS, O_POLYS};

/* These ones *must* be so. Bangs otherwise */
const int perm_hp[3]=	{O_ROWS, O_COLUMNS, O_POLYS};
const int perm_hq[3]=	{O_ROWS, O_COLUMNS, O_POLYS};
const int perm_hr[3]=	{O_ROWS, O_COLUMNS, O_POLYS};

int comp_dp[2]=	{O_ROWS, O_COLUMNS};
int comp_dq[2]=	{O_ROWS, O_COLUMNS};
int comp_dr[2]=	{O_ROWS, O_COLUMNS};

/* These ones not worth being changed */
int trans_p[3]=	{O_ROWS, O_COLUMNS, O_POLYS};
int trans_q[3]=	{O_ROWS, O_COLUMNS, O_POLYS};
int trans_r[3]=	{O_ROWS, O_COLUMNS, O_POLYS};

int perm_conv[4]= {O_ROWS, O_COLUMNS, O_INNER, O_POLYS};
int perm_naive[5]= {O_ROWS, O_COLUMNS, O_INNER, O_IDX_R, O_IDX_Q};

int perm_mix1[3]= {O_ROWS, O_COLUMNS, O_INNER};
int perm_mix2[3]= {O_ROWS, O_COLUMNS, O_INNER};

void checkzero(mp_limb_t * x, size_t n)
{
	size_t i;
	for(i=0;i<n;i++) {
		if (x[i]) {
			die("Warning, got unexpected non-zero\n",1);
		}
	}
}

#define TIME_STACK_DEPTH	32
static double	time_stack[TIME_STACK_DEPTH];
static int	time_stack_pos;

void tm_push(const char * fmt, ...)
{
	va_list ap;

	if (time_stack_pos==TIME_STACK_DEPTH)	BUG();

	time_stack[time_stack_pos++]=timer(TIMER_ASK);

	va_start(ap,fmt);

	vprintf(fmt, ap);

	va_end(ap);
}

void tm_pop_fmt(const char * fmt)
{
	double tt;

	if (time_stack_pos==0)	BUG();

	tt=timer(TIMER_ASK)-time_stack[--time_stack_pos];

	printf(fmt,tt);
	fflush(stdout);
}

#define tm_pop_colon()	tm_pop_fmt("	: %.2fs\n")
#define tm_pop_cnum(n)	tm_pop_fmt("[" #n "G: %.2fs\n")
#define tm_pop_done()	tm_pop_fmt("	done in %.2fs\n")
#define tm_pop_table(n)	tm_pop_fmt("[" #n "Gdone in %.2fs\n")
	
/* Watch out, these macros are rather unsafe... They're here only to
 * avoid a huge number of function calls (and even then, the interest is
 * not obvious)
 */

#define  generic_operation_1op_macro(da,tt,blabla)		\
do {								\
	int z[3];						\
	mp_limb_t * x;						\
	int off[2];						\
	int gap;						\
								\
	x=da->raw;						\
								\
	z[tt[0]]=0;z[tt[1]]=1;z[2]=0;				\
	off[1]=multi_array_offset_tab(da,z);			\
								\
	z[tt[0]]=1;z[tt[1]]=-da->dim[tt[1]];z[2]=0;		\
	off[0]=multi_array_offset_tab(da,z);			\
								\
	z[0]=z[1]=0;						\
	z[2]=1;							\
	gap = multi_array_offset_tab(da,z);			\
								\
	for(z[tt[0]]=0;z[tt[0]]<da->dim[tt[0]];z[tt[0]]++) {	\
	for(z[tt[1]]=0;z[tt[1]]<da->dim[tt[1]];z[tt[1]]++) {	\
		blabla;						\
		x+=off[1];					\
	}							\
	x+=off[0];						\
	}							\
} while (0);

void generic_operation_1op_function(
	struct multi_array * da,
	int tt[2],
	void (*op)(mp_limb_t *, int))
{
	int z[3];
	mp_limb_t * x;
	int off[2];
	int gap;

	x=da->raw;

	z[tt[0]]=0;z[tt[1]]=1;z[2]=0;
	off[1]=multi_array_offset_tab(da,z);

	z[tt[0]]=1;z[tt[1]]=-da->dim[tt[1]];z[2]=0;
	off[0]=multi_array_offset_tab(da,z);

	z[0]=z[1]=0;
	z[2]=1;
	gap = multi_array_offset_tab(da,z);

	for(z[tt[0]]=0;z[tt[0]]<da->dim[tt[0]];z[tt[0]]++) {
	for(z[tt[1]]=0;z[tt[1]]<da->dim[tt[1]];z[tt[1]]++) {
		(*op)(x,gap);
		x+=off[1];
	}
	x+=off[0];
	}
}

#define DIMTEST(v,w,i)	\
	assert(v->dim[i]==w->dim[i] || (i==2 && v->dim[i]>=w->dim[i]))

#define  generic_operation_2op_macro(da,sa,tt,blabla)		\
do {								\
	int z[3];						\
	mp_limb_t * dst, * src;					\
	int doff[3];						\
	int soff[3];						\
	DIMTEST(da,sa,tt[0]);					\
	DIMTEST(da,sa,tt[1]);					\
	DIMTEST(da,sa,tt[2]);					\
	src=sa->raw;						\
	dst=da->raw;						\
								\
	z[tt[0]]=0;z[tt[1]]=0;z[tt[2]]=1;			\
	soff[2]=multi_array_offset_tab(sa,z);			\
	doff[2]=multi_array_offset_tab(da,z);			\
								\
	z[tt[0]]=0;z[tt[1]]=1;z[tt[2]]=-sa->dim[tt[2]];		\
	soff[1]=multi_array_offset_tab(sa,z);			\
	doff[1]=multi_array_offset_tab(da,z);			\
								\
	z[tt[0]]=1;z[tt[1]]=-sa->dim[tt[1]];z[tt[2]]=0;		\
	soff[0]=multi_array_offset_tab(sa,z);			\
	doff[0]=multi_array_offset_tab(da,z);			\
								\
	for(z[tt[0]]=0;z[tt[0]]<sa->dim[tt[0]];z[tt[0]]++) {	\
	for(z[tt[1]]=0;z[tt[1]]<sa->dim[tt[1]];z[tt[1]]++) {	\
	for(z[tt[2]]=0;z[tt[2]]<sa->dim[tt[2]];z[tt[2]]++) {	\
		blabla;						\
		src+=soff[2];					\
		dst+=doff[2];					\
	}							\
	src+=soff[1];						\
	dst+=doff[1];						\
	}							\
	src+=soff[0];						\
	dst+=doff[0];						\
	}							\
} while (0);

void generic_operation_2op_function(
	struct multi_array *da,
	struct multi_array *sa,
	int tt[3],
	void (*op)(mp_limb_t *, mp_limb_t *))
{
	int z[3];
	mp_limb_t * dst, * src;
	int doff[3];
	int soff[3];
	DIMTEST(da,sa,tt[0]);
	DIMTEST(da,sa,tt[1]);
	DIMTEST(da,sa,tt[2]);
	src=sa->raw;
	dst=da->raw;

	z[tt[0]]=0;z[tt[1]]=0;z[tt[2]]=1;
	soff[2]=multi_array_offset_tab(sa,z);
	doff[2]=multi_array_offset_tab(da,z);

	z[tt[0]]=0;z[tt[1]]=1;z[tt[2]]=-sa->dim[tt[2]];
	soff[1]=multi_array_offset_tab(sa,z);
	doff[1]=multi_array_offset_tab(da,z);

	z[tt[0]]=1;z[tt[1]]=-sa->dim[tt[1]];z[tt[2]]=0;
	soff[0]=multi_array_offset_tab(sa,z);
	doff[0]=multi_array_offset_tab(da,z);

	for(z[tt[0]]=0;z[tt[0]]<sa->dim[tt[0]];z[tt[0]]++) {
	for(z[tt[1]]=0;z[tt[1]]<sa->dim[tt[1]];z[tt[1]]++) {
	for(z[tt[2]]=0;z[tt[2]]<sa->dim[tt[2]];z[tt[2]]++) {
		(*op)(dst,src);
		src+=soff[2];
		dst+=doff[2];
	}
	src+=soff[1];
	dst+=doff[1];
	}
	src+=soff[0];
	dst+=doff[0];
	}
} 

#define  generic_operation_3op_macro(zr,zp,zq,tt,blabla)		\
{									\
	int z[4];							\
	int p_pos[4];							\
	int q_pos[4];							\
	int r_pos[4];							\
	int d[4];							\
	int p_off[4];							\
	int q_off[4];							\
	int r_off[4];							\
	int k;								\
	mp_limb_t * yr, * yp, * yq;					\
	assert(zr->dim[O_ROWS]==zp->dim[O_ROWS]);			\
	assert(zr->dim[O_COLUMNS]==zq->dim[O_COLUMNS]);			\
	assert(zp->dim[O_COLUMNS]==zq->dim[O_ROWS]);			\
	assert(zr->dim[O_POLYS]==zp->dim[O_POLYS]);			\
	assert(zr->dim[O_POLYS]==zq->dim[O_POLYS]);			\
	for(k=0;k<4;k++) {						\
		switch(tt[k]) {						\
		    case O_ROWS:					\
			p_pos[k]=O_ROWS;				\
			q_pos[k]=-1;					\
			r_pos[k]=O_ROWS;				\
			d[k]=zp->dim[O_ROWS];				\
			break;						\
		    case O_COLUMNS:					\
			p_pos[k]=-1;					\
			q_pos[k]=O_COLUMNS;				\
			r_pos[k]=O_COLUMNS;				\
			d[k]=zq->dim[O_COLUMNS];			\
			break;						\
		    case O_INNER:					\
			p_pos[k]=O_COLUMNS;				\
			q_pos[k]=O_ROWS;				\
			r_pos[k]=-1;					\
			d[k]=zq->dim[O_ROWS];				\
			break;						\
		    case O_POLYS:					\
			p_pos[k]=O_POLYS;				\
			q_pos[k]=O_POLYS;				\
			r_pos[k]=O_POLYS;				\
			d[k]=zp->dim[O_POLYS];				\
			break;						\
		    default:						\
			BUG();						\
		}							\
	}								\
	for(k=0;k<4;k++) {						\
		int zz[3]={0,0,0};					\
		int  r;							\
		if ((r=p_pos[k])!=-1) {					\
			zz[r]=1;					\
		}							\
		p_off[k]=multi_array_offset_tab(zp,zz);			\
	}								\
	for(k=0;k<4;k++) {						\
		int zz[3]={0,0,0};					\
		int  r;							\
		if ((r=q_pos[k])!=-1) {					\
			zz[r]=1;					\
		}							\
		q_off[k]=multi_array_offset_tab(zq,zz);			\
	}								\
	for(k=0;k<4;k++) {						\
		int zz[3]={0,0,0};					\
		int  r;							\
		if ((r=r_pos[k])!=-1) {					\
			zz[r]=1;					\
		}							\
		r_off[k]=multi_array_offset_tab(zr,zz);			\
	}								\
	for(k=0;k<3;k++) {						\
		p_off[k]+=-d[k+1]*p_off[k+1];				\
		q_off[k]+=-d[k+1]*q_off[k+1];				\
		r_off[k]+=-d[k+1]*r_off[k+1];				\
	}								\
	yr = zr->raw;							\
	yp = zp->raw;							\
	yq = zq->raw;							\
									\
	for(z[0]=0;z[0]<d[0];z[0]++) {					\
	for(z[1]=0;z[1]<d[1];z[1]++) {					\
	for(z[2]=0;z[2]<d[2];z[2]++) {					\
	for(z[3]=0;z[3]<d[3];z[3]++) {					\
		blabla;							\
		yr+=r_off[3];						\
		yp+=p_off[3];						\
		yq+=q_off[3];						\
	}								\
	yr+=r_off[2];							\
	yp+=p_off[2];							\
	yq+=q_off[2];							\
	}								\
	yr+=r_off[1];							\
	yp+=p_off[1];							\
	yq+=q_off[1];							\
	}								\
	yr+=r_off[0];							\
	yp+=p_off[0];							\
	yq+=q_off[0];							\
	}								\
}

void generic_operation_3op_function(
	struct multi_array *zr,
	struct multi_array *zp,
	struct multi_array *zq,
	int tt[4],
	void (*op)(mp_limb_t *, mp_limb_t *, mp_limb_t *))
{
	int z[4];
	int p_pos[4];
	int q_pos[4];
	int r_pos[4];
	int d[4];
	int p_off[4];
	int q_off[4];
	int r_off[4];
	int k;
	mp_limb_t * yr, * yp, * yq;
	assert(zr->dim[O_ROWS]==zp->dim[O_ROWS]);
	assert(zr->dim[O_COLUMNS]==zq->dim[O_COLUMNS]);
	assert(zp->dim[O_COLUMNS]==zq->dim[O_ROWS]);
	assert(zr->dim[O_POLYS]==zp->dim[O_POLYS]);
	assert(zr->dim[O_POLYS]==zq->dim[O_POLYS]);
	for(k=0;k<4;k++) {
		switch(tt[k]) {
		    case O_ROWS:
			p_pos[k]=O_ROWS;
			q_pos[k]=-1;
			r_pos[k]=O_ROWS;
			d[k]=zp->dim[O_ROWS];
			break;
		    case O_COLUMNS:
			p_pos[k]=-1;
			q_pos[k]=O_COLUMNS;
			r_pos[k]=O_COLUMNS;
			d[k]=zq->dim[O_COLUMNS];
			break;
		    case O_INNER:
			p_pos[k]=O_COLUMNS;
			q_pos[k]=O_ROWS;
			r_pos[k]=-1;
			d[k]=zq->dim[O_ROWS];
			break;
		    case O_POLYS:
			p_pos[k]=O_POLYS;
			q_pos[k]=O_POLYS;
			r_pos[k]=O_POLYS;
			d[k]=zp->dim[O_POLYS];
			break;
		    default:
			BUG();
		}
	}
	for(k=0;k<4;k++) {
		int zz[3]={0,0,0};
		int  r;
		if ((r=p_pos[k])!=-1) {
			zz[r]=1;
		}
		p_off[k]=multi_array_offset_tab(zp,zz);
	}
	for(k=0;k<4;k++) {
		int zz[3]={0,0,0};
		int  r;
		if ((r=q_pos[k])!=-1) {
			zz[r]=1;
		}
		q_off[k]=multi_array_offset_tab(zq,zz);
	}
	for(k=0;k<4;k++) {
		int zz[3]={0,0,0};
		int  r;
		if ((r=r_pos[k])!=-1) {
			zz[r]=1;
		}
		r_off[k]=multi_array_offset_tab(zr,zz);
	}
	for(k=0;k<3;k++) {
		p_off[k]+=-d[k+1]*p_off[k+1];
		q_off[k]+=-d[k+1]*q_off[k+1];
		r_off[k]+=-d[k+1]*r_off[k+1];
	}
	yr = zr->raw;
	yp = zp->raw;
	yq = zq->raw;

	for(z[0]=0;z[0]<d[0];z[0]++) {
	for(z[1]=0;z[1]<d[1];z[1]++) {
	for(z[2]=0;z[2]<d[2];z[2]++) {
	for(z[3]=0;z[3]<d[3];z[3]++) {
			(*op)(yr,yp,yq);
		yr+=r_off[3];
		yp+=p_off[3];
		yq+=q_off[3];
	}
	yr+=r_off[2];
	yp+=p_off[2];
	yq+=q_off[2];
	}
	yr+=r_off[1];
	yp+=p_off[1];
	yq+=q_off[1];
	}
	yr+=r_off[0];
	yp+=p_off[0];
	yq+=q_off[0];
	}
}

void extend(mp_limb_t * dst, mp_limb_t * src)
{
	memcpy(dst,src,fdeg*sizeof(mp_limb_t));
	memset(dst+fdeg,0,fdeg*sizeof(mp_limb_t));
}

/* Watch out, the names are misleading here. It's here only because we
 * want the caller to be driven by the length of the table which yields
 * the src argument to this function. This argument is not necessarily
 * the source argument.
 */
void do_restrict(mp_limb_t * dst, mp_limb_t * src)
{
	/*	memcpy(src,dst,fdeg*sizeof(mp_limb_t)); */
	prime_field_mul(src,dst,one_over_n,field_k);
	checkzero(dst+fdeg,fdeg);
}

void do_dft(mp_limb_t * x, int gap)
{
	(*direct_fft)(x,gap,l_roots,order);
}

void do_ift(mp_limb_t * x, int gap)
{
	(*inverse_fft)(x,gap,l_roots,order);
}

void do_conv(mp_limb_t * r, mp_limb_t * p, mp_limb_t * q)
{
	if (enable_cplx) {
		complex_field_addmul(r,p,q,field_l);
	} else {
		quad_field_addmul(r,p,q,field_l);
	}
}

void do_naive_mult(
	struct multi_array *zr,
	struct multi_array *zp,
	struct multi_array *zq,
	int tt[5])
{
	int z[5];
	int p_pos[5];
	int q_pos[5];
	int r_pos[5];
	int d[5];
	int p_off[5];
	int q_off[5];
	int r_off[5];
	int emul=0;
	int emul_off[5];
	int d_emul=0;
	int k;
	mp_limb_t * yr, * yp, * yq;
	int *ps,*pt,*pu;	/* maintained such that s+t=u */
	memset(zr->raw,0,zr->nlimbs*sizeof(mp_limb_t));
	ps=pt=pu=NULL;
	assert(zr->dim[0]==zp->dim[0]);	/* rows in p & r */
	assert(zr->dim[1]==zq->dim[1]);	/* cols in q & r */
	assert(zp->dim[1]==zq->dim[0]);	/* cols in p, rows in r */
	for(k=0;k<5;k++) {
		switch(tt[k]) {
		    case O_ROWS:
			p_pos[k]=O_ROWS;
			q_pos[k]=-1;
			r_pos[k]=O_ROWS;
			d[k]=zp->dim[O_ROWS];
			break;
		    case O_COLUMNS:
			p_pos[k]=-1;
			q_pos[k]=O_COLUMNS;
			r_pos[k]=O_COLUMNS;
			d[k]=zq->dim[O_COLUMNS];
			break;
		    case O_INNER:
			p_pos[k]=O_COLUMNS;
			q_pos[k]=O_ROWS;
			r_pos[k]=-1;
			d[k]=zq->dim[O_ROWS];
			break;
		    case O_IDX_P:
			p_pos[k]=O_POLYS;
			q_pos[k]=-1;
			r_pos[k]=-1;
			d[k]=zp->dim[O_POLYS];
			ps=z+k;
			break;
		    case O_IDX_Q:
			p_pos[k]=-1;
			q_pos[k]=O_POLYS;
			r_pos[k]=-1;
			d[k]=zq->dim[O_POLYS];
			pt=z+k;
			break;
		    case O_IDX_R:
			p_pos[k]=-1;
			q_pos[k]=-1;
			r_pos[k]=O_POLYS;
			d[k]=zr->dim[O_POLYS];
			pu=z+k;
			break;
		    default:
			BUG();
		}
	}

	if (!ps) {emul=O_IDX_P; d_emul=zp->dim[O_POLYS];}
	if (!pt) {emul=O_IDX_Q; d_emul=zq->dim[O_POLYS];}
	if (!pu) {emul=O_IDX_R; d_emul=zr->dim[O_POLYS];}

	for(k=0;k<5;k++) emul_off[k]=0;

	for(k=0;k<5;k++) {
		int zz[3]={0,0,0};
		int  r;
		if (emul==O_IDX_P && tt[k]==O_IDX_Q) {
			emul_off[k]=zz[O_POLYS]=-1;
		} else if  (emul==O_IDX_P && tt[k]==O_IDX_R) {
			emul_off[k]=zz[O_POLYS]=1;
		} else if ((r=p_pos[k])!=-1) {
			zz[r]=1;
		}
		p_off[k]=multi_array_offset_tab(zp,zz);
	}
	for(k=0;k<5;k++) {
		int zz[3]={0,0,0};
		int  r;
		if (emul==O_IDX_Q && tt[k]==O_IDX_P) {
			emul_off[k]=zz[O_POLYS]=-1;
		} else if  (emul==O_IDX_Q && tt[k]==O_IDX_R) {
			emul_off[k]=zz[O_POLYS]=1;
		} else if ((r=q_pos[k])!=-1) {
			zz[r]=1;
		}
		q_off[k]=multi_array_offset_tab(zq,zz);
	}
	for(k=0;k<5;k++) {
		int zz[3]={0,0,0};
		int  r;
		if (emul==O_IDX_R && tt[k]==O_IDX_P) {
			emul_off[k]=zz[O_POLYS]=1;
		} else if  (emul==O_IDX_R && tt[k]==O_IDX_Q) {
			emul_off[k]=zz[O_POLYS]=1;
		} else if ((r=r_pos[k])!=-1) {
			zz[r]=1;
		}
		r_off[k]=multi_array_offset_tab(zr,zz);
	}
	for(k=0;k<4;k++) {
		p_off[k]+=-d[k+1]*p_off[k+1];
		q_off[k]+=-d[k+1]*q_off[k+1];
		r_off[k]+=-d[k+1]*r_off[k+1];
		emul_off[k]+=-d[k+1]*emul_off[k+1];
	}
	yr = zr->raw;
	yp = zp->raw;
	yq = zq->raw;
	emul=0;	/* We reuse the counter */

	for(z[0]=0;z[0]<d[0];z[0]++) {
	for(z[1]=0;z[1]<d[1];z[1]++) {
	for(z[2]=0;z[2]<d[2];z[2]++) {
	for(z[3]=0;z[3]<d[3];z[3]++) {
	for(z[4]=0;z[4]<d[4];z[4]++) {
		if (emul>=0 && emul<d_emul)
			prime_field_addmul(yr,yp,yq,field_k);
		yr+=r_off[4];
		yp+=p_off[4];
		yq+=q_off[4];
		emul+=emul_off[4];
	}
	yr+=r_off[3];
	yp+=p_off[3];
	yq+=q_off[3];
	emul+=emul_off[3];
	}
	yr+=r_off[2];
	yp+=p_off[2];
	yq+=q_off[2];
	emul+=emul_off[2];
	}
	yr+=r_off[1];
	yp+=p_off[1];
	yq+=q_off[1];
	emul+=emul_off[1];
	}
	yr+=r_off[0];
	yp+=p_off[0];
	yq+=q_off[0];
	emul+=emul_off[0];
	}
}

static void do_my_fft_mul(
	struct multi_array *pr,
	struct multi_array *pp,
	struct multi_array *pq)
{
	struct multi_array * dp;
	struct multi_array * dq;
	struct multi_array * dr;

	dp = declare_multi_array(dim_dp,perm_dp,fdeg<<1);
	dq = declare_multi_array(dim_dq,perm_dq,fdeg<<1);
	dr = declare_multi_array(dim_dr,perm_dr,fdeg<<1);

	memset(dp->raw,0,dp->nlimbs*sizeof(mp_limb_t));
	memset(dq->raw,0,dq->nlimbs*sizeof(mp_limb_t));
	memset(dr->raw,0,dr->nlimbs*sizeof(mp_limb_t));

	tm_push("Filling in the DFT areas.");
#ifdef USE_MACROS
	generic_operation_2op_macro(dp,pp,trans_p,
		memcpy(dst,src,fdeg*sizeof(mp_limb_t));
		memset(dst+fdeg,0,fdeg*sizeof(mp_limb_t)));

	generic_operation_2op_macro(dq,pq,trans_q,
		memcpy(dst,src,fdeg*sizeof(mp_limb_t));
		memset(dst+fdeg,0,fdeg*sizeof(mp_limb_t)));
#else	/* USE_MACROS */
	generic_operation_2op_function(dp,pp,trans_p,&extend);
	generic_operation_2op_function(dq,pq,trans_q,&extend);
#endif	/* USE_MACROS */
	tm_pop_table(40);

	tm_push("Computing the DFTs.");
#ifdef USE_MACROS
	generic_operation_1op_macro(dp,perm_dp,
			(*direct_fft)(x,gap,l_roots,order));

	generic_operation_1op_macro(dq,perm_dq,
			(*direct_fft)(x,gap,l_roots,order));
#else	/* USE_MACROS */
	generic_operation_1op_function(dp,perm_dp,&do_dft);
	generic_operation_1op_function(dq,perm_dq,&do_dft);
#endif	/* USE_MACROS */
	tm_pop_table(40);

	tm_push("Computing the convolution product.");
#ifdef USE_MACROS
	if (enable_cplx) {
		generic_operation_3op_macro(dr,dp,dq,perm_conv,
				complex_field_addmul(yr,yp,yq,field_l));
	} else {
		generic_operation_3op_macro(dr,dp,dq,perm_conv,
				quad_field_addmul(yr,yp,yq,field_l));
	}
#else	/* USE_MACROS */
	generic_operation_3op_function(dr,dp,dq,perm_conv,&do_conv);
#endif	/* USE_MACROS */
	tm_pop_table(40);
			
	tm_push("Computing the inverse DFT.");
#ifdef USE_MACROS
	generic_operation_1op_macro(dr,perm_dr,
			(*inverse_fft)(x,gap,l_roots,order));
#else	/* USE_MACROS */
	generic_operation_1op_function(dr,perm_dr,&do_ift);
#endif	/* USE_MACROS */
	tm_pop_table(40);

	tm_push("Storing the result.");
#ifdef USE_MACROS
	generic_operation_2op_macro(dr,pr,trans_r,
		/*memcpy(src,dst,fdeg*sizeof(mp_limb_t));*/
		prime_field_mul(src,dst,one_over_n,field_k);
		checkzero(dst+fdeg,fdeg));
#else	/* USE_MACROS */
	generic_operation_2op_function(dr,pr,trans_r,&do_restrict);
#endif	/* USE_MACROS */
	tm_pop_table(40);

	delete_multi_array(dp);
	delete_multi_array(dq);
	delete_multi_array(dr);
}

void do_mix_gmp_mult(
	struct multi_array *zr,
	struct multi_array *zp,
	struct multi_array *zq,
	int tt[3])
{
	int z[3];
	int p_pos[3];
	int q_pos[3];
	int r_pos[3];
	int d[3];
	int p_off[3];
	int q_off[3];
	int r_off[3];
	int k;
	mp_limb_t * yr, * yp, * yq;
	struct multi_array * hp;
	struct multi_array * hq;
	struct multi_array * hr;
	int nlimbs_p, nlimbs_q, nlimbs_r;
	mp_limb_t * buf;
	int padded;
	int mou_hr;

	padded=1+(fdeg<<1);	/* must be at least fdeg<<1 */

	/* The good size is fdeg*2+ceil(log2(min(deg p, deg q)+1)) */

	/* One should also add ceil(log2(pcols)) */

	hp = declare_multi_array(dim_p,perm_hp,padded);
	hq = declare_multi_array(dim_q,perm_hq,padded);
	hr = declare_multi_array(dim_r,perm_hr,padded);

	memset(hp->raw,0,hp->nlimbs*sizeof(mp_limb_t));
	memset(hq->raw,0,hq->nlimbs*sizeof(mp_limb_t));
	memset(hr->raw,0,hr->nlimbs*sizeof(mp_limb_t));

	generic_operation_2op_macro(hp,zp,trans_p,
		memcpy(dst,src,fdeg*sizeof(mp_limb_t));
		memset(dst+fdeg,0,(padded-fdeg)*sizeof(mp_limb_t)));
	generic_operation_2op_macro(hq,zq,trans_q,
		memcpy(dst,src,fdeg*sizeof(mp_limb_t));
		memset(dst+fdeg,0,(padded-fdeg)*sizeof(mp_limb_t)));
	
	nlimbs_p=padded*(dim_p[2]-1)+fdeg;
	nlimbs_q=padded*(dim_q[2]-1)+fdeg;
	nlimbs_r=padded*(dim_r[2]-1)+(fdeg<<1);

	assert(nlimbs_p+nlimbs_q==nlimbs_r);
	mou_hr=hr->dim[2]*hr->ext-nlimbs_r;
	assert(mou_hr>=0);

	buf=(mp_limb_t *) malloc(MAX(nlimbs_r,1+padded)*sizeof(mp_limb_t));
	
	memset(zr->raw,0,zr->nlimbs*sizeof(mp_limb_t));
	assert(zr->dim[0]==zp->dim[0]);	/* rows in p & r */
	assert(zr->dim[1]==zq->dim[1]);	/* cols in q & r */
	assert(zp->dim[1]==zq->dim[0]);	/* cols in p, rows in r */
	for(k=0;k<3;k++) {
		switch(tt[k]) {
		    case O_ROWS:
			p_pos[k]=O_ROWS;
			q_pos[k]=-1;
			r_pos[k]=O_ROWS;
			d[k]=zp->dim[O_ROWS];
			break;
		    case O_COLUMNS:
			p_pos[k]=-1;
			q_pos[k]=O_COLUMNS;
			r_pos[k]=O_COLUMNS;
			d[k]=zq->dim[O_COLUMNS];
			break;
		    case O_INNER:
			p_pos[k]=O_COLUMNS;
			q_pos[k]=O_ROWS;
			r_pos[k]=-1;
			d[k]=zq->dim[O_ROWS];
			break;
		    default:
			BUG();
		}
	}

	for(k=0;k<3;k++) {
		int zz[3]={0,0,0};
		int  r;
		if ((r=p_pos[k])!=-1) {
			zz[r]=1;
		}
		p_off[k]=multi_array_offset_tab(hp,zz);
	}
	for(k=0;k<3;k++) {
		int zz[3]={0,0,0};
		int  r;
		if ((r=q_pos[k])!=-1) {
			zz[r]=1;
		}
		q_off[k]=multi_array_offset_tab(hq,zz);
	}
	for(k=0;k<3;k++) {
		int zz[3]={0,0,0};
		int  r;
		if ((r=r_pos[k])!=-1) {
			zz[r]=1;
		}
		r_off[k]=multi_array_offset_tab(hr,zz);
	}
	for(k=0;k<2;k++) {
		p_off[k]+=-d[k+1]*p_off[k+1];
		q_off[k]+=-d[k+1]*q_off[k+1];
		r_off[k]+=-d[k+1]*r_off[k+1];
	}
	yr = hr->raw;
	yp = hp->raw;
	yq = hq->raw;
	
	for(z[0]=0;z[0]<d[0];z[0]++) {
	for(z[1]=0;z[1]<d[1];z[1]++) {
	for(z[2]=0;z[2]<d[2];z[2]++) {
		mp_limb_t c;
		mpn_mul(buf,yp,nlimbs_p,yq,nlimbs_q);
		if ((c=mpn_add_n(yr,yr,buf,nlimbs_r))) {
			assert(mou_hr);
			c=mpn_add_1(yr+nlimbs_r,yr+nlimbs_r,c,mou_hr);
			assert(!c);
		}
		yr+=r_off[2];
		yp+=p_off[2];
		yq+=q_off[2];
	}
	yr+=r_off[1];
	yp+=p_off[1];
	yq+=q_off[1];
	}
	yr+=r_off[0];
	yp+=p_off[0];
	yq+=q_off[0];
	}

	assert(zr->nlimbs%fdeg==0);
	for(k=0;k< (int) zr->nlimbs/fdeg;k++) {
		memcpy(buf,hr->raw+k*padded,padded*sizeof(mp_limb_t));
		prime_field_set_mpn(zr->raw+k*fdeg,buf,padded,field_k);
	}

	free(buf);
}

static void	do_tests(int maxdeg)
{
	struct multi_array * pp;
	struct multi_array * pq;
	struct multi_array * pr;
	
	struct multi_array * pr_saved[4];
	int i;

	mp_limb_t * buf;

	pp = declare_multi_array(dim_p,perm_p,fdeg);
	pq = declare_multi_array(dim_q,perm_q,fdeg);
	pr = declare_multi_array(dim_r,perm_r,fdeg);

	printf("\n");

	tm_push("Filling the input with random data.");
	buf=(mp_limb_t *) malloc((fdeg+1)*sizeof(mp_limb_t));
	for(i=0;i<(int) pp->nlimbs/fdeg;i++) {
		int j;
		for(j=0;j<fdeg;j++) {
			buf[j]=rand();
		}
		prime_field_set_mpn(pp->raw+i*fdeg,buf,fdeg,field_k);
	}
	for(i=0;i<(int) pq->nlimbs/fdeg;i++) {
		int j;
		for(j=0;j<fdeg;j++) {
			buf[j]=rand();
		}
		prime_field_set_mpn(pq->raw+i*fdeg,buf,fdeg,field_k);
	}
	free(buf);
	tm_pop_table(40);

	tm_push("\n---------------------------------------------------\n"
		"*** METHOD 1 : my FFT\n\n");


	do_my_fft_mul(pr,pp,pq);
	tm_pop_fmt("my FFT: done in %.2fs\n\n");

	pr_saved[0]=pr;
	pr = declare_multi_array(dim_r,perm_r,fdeg);

	/**************************************************************/
	tm_push("\n---------------------------------------------------\n"
		"*** METHOD 2 : naive, with Gmp's fft\n");

	do_mix_gmp_mult(pr,pp,pq,perm_mix1);
	tm_pop_fmt("naive/GMP 1: done in %.2fs\n\n");

	pr_saved[1]=pr;
	pr = declare_multi_array(dim_r,perm_r,fdeg);

#if 0
	/**************************************************************/
	tm_push("\n---------------------------------------------------\n"
		"*** METHOD 3 : completely naive\n");

	do_naive_mult(pr,pp,pq,perm_naive);
	tm_pop_fmt("completely naive: done in %.2fs\n\n");

	pr_saved[2]=pr;
	pr = declare_multi_array(dim_r,perm_r,fdeg);
#else
	pr_saved[2]=pr_saved[1];
#endif

	/**************************************************************/

	if (	memcmp(pr_saved[1]->raw,pr_saved[0]->raw,
				pr->nlimbs*sizeof(mp_limb_t))!=0 ||
		memcmp(pr_saved[2]->raw,pr_saved[0]->raw,
				pr->nlimbs*sizeof(mp_limb_t))!=0)
	{
		printf("Got different results !!!\n");
		pr_ma(pr_saved[0]);
		pr_ma(pr_saved[1]);
		pr_ma(pr_saved[2]);
	} else {
		printf("OK, results match\n");
	}
}

int main(int argc, char *argv[])
{
	struct opt_desc * opts = NULL;
	int n_opts=0;
	int dm,dn,db,degree,tmp;
	char * modul = NULL;
	
	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	check_endianness(stdout);

	new_option(&opts,&n_opts,
			OPT_FILENAME(NULL, &modul,OPT_REQUIRED|OPT_IGNFURTHER));
	new_option(&opts,&n_opts,
			OPT_INTPRM(NULL,&dm,OPT_INDEP|OPT_REQUIRED));
	new_option(&opts,&n_opts,
			OPT_INTPRM(NULL,&dn,OPT_INDEP|OPT_REQUIRED));
	new_option(&opts,&n_opts,
			OPT_INTPRM(NULL,&db,OPT_INDEP|OPT_REQUIRED));
	new_option(&opts,&n_opts,
			OPT_INTPRM(NULL,&degree,OPT_INDEP|OPT_REQUIRED));
	new_option(&opts,&n_opts,
			OPT_FLAG("--enable-complex-field",&enable_cplx));
	
	process_options(argc,argv,n_opts,opts);

	coredump_limit(0);

	// set_all_filenames(bank_num);

	/* This is necessary since we want the modulus info. */
	if (bw_read_modulus_info(modul, 1)<0)
		die("Unable to read info for modulus %s\n",1,modul);

	printf("Testing (%d,%d) x (%d,%d) -> (%d,%d)\n",dm,dn,dn,db,dm,db);

	tmp=2*degree;
	for(order=0;tmp>>order;order++);

	dim_p[0]=dm; dim_p[1]=dn; dim_p[2]=degree+1;
	dim_q[0]=dn; dim_q[1]=db; dim_q[2]=degree+1;
	dim_r[0]=dm; dim_r[1]=db; dim_r[2]=tmp+1;

	dim_dp[0]=dm; dim_dp[1]=dn; dim_dp[2]=(1<<order);
	dim_dq[0]=dn; dim_dq[1]=db; dim_dq[2]=(1<<order);
	dim_dr[0]=dm; dim_dr[1]=db; dim_dr[2]=(1<<order);

	printf("Testing at degree %d, order = %d\n",degree,order);

	prepare_fft_engine(order,enable_cplx);

	one_over_n=(mp_limb_t *) malloc(k_size*sizeof(mp_limb_t));
	k_set_int(one_over_n,1<<order);
	k_inv(one_over_n,one_over_n);

	do_tests(degree);

	free(one_over_n);

	return 0;
}
