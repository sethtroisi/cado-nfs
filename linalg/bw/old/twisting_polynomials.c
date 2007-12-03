#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include "field_def.h"
#include "field_usage.h"

static int * degtable_entry(struct t_poly * tp, int i, int j)
{
	return tp->degree_table+i*bigdim+j;
}

static struct t_poly * tp_pre_alloc(void)
{
	struct t_poly * res;

	if ((res = malloc(sizeof(struct t_poly)))==NULL)
		goto fail_0;

	if ((res->degree_table=malloc(bigdim*bigdim*sizeof(int)))==NULL)
		goto fail_1;

	if ((res->clist=malloc(bigdim*sizeof(unsigned int)))==NULL)
		goto fail_2;

	if ((res->degnom=malloc(bigdim*sizeof(int)))==NULL)
		goto fail_3;

	res->alloc=0;
	res->degree=-1;
	return res;

	/* free(res->degnom); -- never done... */
fail_3:	free(res->clist);
fail_2:	free(res->degree_table);
fail_1:	free(res);
fail_0:	return NULL;
}

static struct t_poly * tp_post_alloc(struct t_poly *res, unsigned int degree)
{
	bbpoly_alloc(res->p,degree);
	if (STRICTTYPE_VAL(res->p) == NULL) {
		free(res->degnom);
		free(res->clist);
		free(res->degree_table);
		free(res);
		return NULL;
	}
	res->alloc=degree+1;
	res->degree=degree;
	return res;
}

struct t_poly * tp_alloc(unsigned int degree)
{
	return tp_post_alloc(tp_pre_alloc(),degree);
}

void tp_free(struct t_poly * pi)
{
	bbpoly_free(pi->p);
	free(pi->degnom);
	free(pi->clist);
	free(pi->degree_table);
	free(pi);
}
	
static int nbar_add(int a, int b)
{
	if (a==-1 || b==-1)
		return -1;
	else
		return a+b;
}

struct t_poly * tp_comp_alloc(struct t_poly * left,struct t_poly * right)
{
	int maxdeg,deg;
	int i,j,k;
	struct t_poly *res;
	int *d;

	res=tp_pre_alloc();

	maxdeg=-1;

	for(j=0;j<bigdim;j++) {
		res->clist[j]=j;
		res->degnom[j]=-1;
		for(i=0;i<bigdim;i++) {
			d=degtable_entry(res,i,j);
			*d=-1;
			for(k=0;k<bigdim;k++) {
				deg=nbar_add(
				    *degtable_entry(left,i,left->clist[k]),
				    *degtable_entry(right,k,right->clist[j]));
				if (deg > *d)
					*d=deg;
			}
			if (*d > res->degnom[j]) {
				res->degnom[j]=*d;
				if (*d > maxdeg)
					maxdeg=*d;
			}
		}
	}
		
	res=tp_post_alloc(res,maxdeg);

	if (res==NULL)
		return NULL;

	bbpoly_zero(res->p,maxdeg);
	return res;
}

void tp_act_on_delta(struct t_poly *tp, int * delta)
{
	int * tmp_delta, d, i, j;
	tmp_delta=malloc(bigdim*sizeof(int));
	for(j=0;j<bigdim;j++) {
		tmp_delta[j]=-1;
		for(i=0;i<bigdim;i++) {
			d=nbar_add(delta[i],*degtable_entry(tp,i,tp->clist[j]));
			if (d>tmp_delta[j])
				tmp_delta[j]=d;
		}
	}
	memcpy(delta,tmp_delta,bigdim*sizeof(int));
	free(tmp_delta);
}

void tp_set_ident(struct t_poly * tp)
{
	int i,j;

	tp->degree=0;

	for(j=0;j<bigdim;j++) {
		tp->degnom[j]=0;
		tp->clist[j]=j;
		for(i=0;i<tp->alloc;i++) {
			bcol_zero(bbmat_col(bbpoly_coeff(tp->p,i),j));
		}
		bw_scalar_set_one(bbmat_scal(bbpoly_coeff(tp->p,0),j,j));
		for(i=0;i<bigdim;i++) {
			*degtable_entry(tp,i,j)=(i==j)?0:-1;
		}
	}
}

void
tp_exchange_cols(struct t_poly * tp, int j1, int j2)
{
	int j0;
	j0=tp->clist[j1];
	tp->clist[j1]=tp->clist[j2];
	tp->clist[j2]=j0;
}

/* Operation : C_j <- C_{\sigma(j)} */
void
tp_apply_perm(struct t_poly * tp, unsigned int * sigma)
{
	unsigned int * tmp;
	int i;

	tmp=malloc(bigdim * sizeof(unsigned int));
	for(i=0;i<bigdim;i++) {
		tmp[i]=tp->clist[sigma[i]];
	}
	memcpy(tp->clist,tmp,bigdim*sizeof(unsigned int));
	free(tmp);
}

/* Operation : C2<-C2+lambda*C1 */
void
tp_transvec(struct t_poly * tp,
	       	int o_j1,
	       	int o_j2,
		int ikill UNUSED_VARIABLE,
	       	bw_scalar lambda)
{
	int i,t;
	int j1,j2;
	int d0;
	int *nd;

	j1=tp->clist[o_j1];
	j2=tp->clist[o_j2];

#if 0
	assert(tp->degnom[j1] <= tp->degnom[j2]);
#endif

	/* Apparently this might happen, and rightfully. One thing is
	 * certain, it's that the _deltas_ corresponding to these columns
	 * are in order. As for the local nominal degrees, anything can
	 * happen (particularly for small moduli.
	 */
	if (tp->degnom[j1] > tp->degnom[j2]) {
		printf("wild degree increase in tp_transvec\n");
		tp->degnom[j2]=tp->degnom[j1];
	}

	if (bw_scalar_is_zero(lambda))
		return;

	for(i=0;i<bigdim;i++) {
		d0=*degtable_entry(tp,i,j1);
		for(t=0;t<=tp->degnom[j1];t++) {
			nd=degtable_entry(tp,i,j2);
			if (d0>*nd)
				*nd=d0;
			addmul( bbmat_scal(bbpoly_coeff(tp->p,t),i,j2),
				bbmat_scal(bbpoly_coeff(tp->p,t),i,j1),
				lambda);
			k_reduce(bbmat_scal(bbpoly_coeff(tp->p,t),i,j2));

		}
	}
}

void
tp_x_multiply(struct t_poly * tp, int o_j)
{
	int t;
	int j;
	int i, *d;

	j=tp->clist[o_j];

	assert(tp->degnom[j] + 1 < tp->alloc);

	for(i=0;i<bigdim;i++) {
		d=degtable_entry(tp,i,j);
		*d=nbar_add(*d,1);
	}

	for ( t = tp->degnom[j] ; t >= 0 ; t-- ) {
		bcol_copy(	bbmat_col(bbpoly_coeff(tp->p,t+1),j),
				bbmat_col(bbpoly_coeff(tp->p,t),j));
	}

	tp->degnom[j]++;
	if (tp->degnom[j] > tp->degree)
		tp->degree++;
	
	bcol_zero(bbmat_col(bbpoly_coeff(tp->p,0),j));
}

void tp_print(struct t_poly * tp)
{
	int i, j, t, jr;

	printf("[ ");
	for(i = 0 ; i < bigdim ; i++) {
		for(j = 0 ; j < bigdim ; j++) {
			jr=tp->clist[j];
			for(t=0;t<=*degtable_entry(tp,i,jr);t++) {
				if (t) printf(" + Y^%d * ",t);
				k_print(bbmat_scal(bbpoly_coeff(tp->p,t),i,jr));
			}
			if (t==0)
				printf("0");
			if (i+1<bigdim || j+1<bigdim)
				printf(",	");
		}
		if (i+1<bigdim)
			printf("\n  ");
	}
	printf("	]\n");
}

int tp_write(FILE *f, struct t_poly *tp)
{
	int i,j,k,jr;
	mpz_t blah;

	mpz_init(blah);
	
	for(k=0;k<=tp->degree;k++) {
		for(i=0;i<bigdim;i++) {
			for(j=0;j<bigdim;j++) {
				jr=tp->clist[j];
				if (k <= tp->degnom[jr]) {
					MPZ_SET_MPN(blah, bbmat_scal( bbpoly_coeff(tp->p,k),i,jr), bw_allocsize);
					gmp_fprintf(f, "%Zd%c", blah,
							(j==bigdim-1)?'\n':' ');
				} else {
					gmp_fprintf(f, "0%c",
							(j==bigdim-1)?'\n':' ');
				}
			}
		}
		fprintf(f, "\n");
	}

	mpz_clear(blah);
	return 0;
}

struct t_poly * tp_read(FILE * f, unsigned int n)
{
	int i,j,k,dmax;
	struct t_poly *res;
	mpz_t blah;

	dmax = -1;

	res=tp_alloc(n);
	mpz_init(blah);
	
	for(j=0;j<bigdim;j++) {
		res->clist[j]=j;
		for(i=0;i<bigdim;i++) {
			*degtable_entry(res, i, j) = -1;
		}
	}

	int rc = 0;
	for(k=0;k<=res->degree;k++) {
		rc = 0;
		for(i=0;i<bigdim;i++) {
			for(j=0;j<bigdim;j++) {
				if (gmp_fscanf(f, "%Zd", blah) != 1)
					continue;
				rc++;
				MPN_SET_MPZ(bbmat_scal(bbpoly_coeff(res->p,k),i,j), bw_allocsize, blah);
				if (mpz_cmp_ui(blah, 0) != 0) {
					*degtable_entry(res, i, j) = k;
				}
			}
		}
		if (rc == 0) {
			res->degree = k-1;
			break;
		} else if (rc != bigdim * bigdim) {
			die("Reading pi file : not in sync\n",1);
		}
	}
	if (rc != 0) {
		die("Estimation of pi degree was wrong (more than %d)\n", 1, n);
	}
	for(j=0;j<bigdim;j++) {
		int d = -1;
		for(i=0;i<bigdim;i++) {
			int d2 = *degtable_entry(res, i, j);
			if (d2 > d)
				d = d2;
		}
		res->degnom[j]=d;
	}
	mpz_clear(blah);
	return res;
}
