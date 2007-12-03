#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include "e_polynomial.h"
#include "field_def.h"
#include "field_usage.h"

struct e_coeff * ec_encapsulate(bw_mbpoly p,
	       unsigned int * clist,
	       int degree)
{
	struct e_coeff *res;
	res=malloc(sizeof(struct e_coeff));
	if (res==NULL) return NULL;
	res->p=res->original_p=p;
	res->clist=clist;
	res->degree=degree;
	return res;
}

void ec_unlink(struct e_coeff *ec)
{
	free(ec);
}
	
/* Operation : C_j <- C_{\sigma(j)} */
void
ec_apply_perm(struct e_coeff * ec, unsigned int * sigma)
{
	unsigned int * tmp;
	unsigned int i;

	tmp=malloc(bigdim * sizeof(unsigned int));
	for(i=0;i<bigdim;i++) {
		tmp[i]=ec->clist[sigma[i]];
	}
	memcpy(ec->clist,tmp,bigdim*sizeof(unsigned int));
	free(tmp);
}

void
ec_exchange_cols(struct e_coeff * ec, int j1, int j2)
{
	unsigned int foo;
	memcpy(&foo, &(ec->clist[j1]), sizeof(unsigned int));
	memcpy(&(ec->clist[j1]), &(ec->clist[j2]), sizeof(unsigned int));
	memcpy(&(ec->clist[j2]), &foo, sizeof(unsigned int));
}

/*
 * Cancel a column.
 * The reference column is j1. The column j2 is added lambda times the
 * column j1 so that the coefficient at (ikill,j2) cancels out.
 *
 * On input, (ikill,j1) must be the first nonzero coefficient on the j1
 * column. Thus, the coefficients for i<ikill in column j1 and j2 are left
 * untouched. Coefficient (ikill,j2) is zeroed.
 *
 */
void
ec_transvec(struct e_coeff * ec,
		int o_j1,
		int o_j2,
		int ikill,
		bw_scalar lambda)
{
	int i,k;
	int j1,j2;

	j1=ec->clist[o_j1];
	j2=ec->clist[o_j2];

	if (bw_scalar_is_zero(lambda))
		return;

	bw_scalar_set_zero(mbmat_scal(mbpoly_coeff(ec->p,0),ikill,j2));
	

	for(i=ikill+1;i<m_param;i++) {
		addmul( mbmat_scal(mbpoly_coeff(ec->p,0),i,j2),
			mbmat_scal(mbpoly_coeff(ec->p,0),i,j1),
			lambda);
	}
	for(k=1;k<=ec->degree;k++) for(i=0;i<m_param;i++) {
		addmul( mbmat_scal(mbpoly_coeff(ec->p,k),i,j2),
			mbmat_scal(mbpoly_coeff(ec->p,k),i,j1),
			lambda);
	}
}

void
ec_x_multiply(struct e_coeff * ec, int o_j)
{
	int j,k;

	j=ec->clist[o_j];
	for(k=ec->degree;k>0;k--) {
		mcol_copy(mbmat_col(mbpoly_coeff(ec->p,k),j),
				mbmat_col(mbpoly_coeff(ec->p,k-1),j));
	}
	mcol_zero(mbmat_col(mbpoly_coeff(ec->p,0),j));
}

/* XXX This is rather ill-named...
 */
void ec_advance(struct e_coeff *ec, int n)
{
	ec->p=mbpoly_subpoly(ec->p,n);
	ec->degree-=n;
}

void ec_park(struct e_coeff * ec)
{
	ec->p=ec->original_p;
}

int ec_is_twisted(struct e_coeff * ec)
{
	unsigned int j;
	for(j=0;j<bigdim;j++) if (ec->clist[j]!=j)
		return 1;
	return 0;
}

void ec_untwist(struct e_coeff * ec)
{
	unsigned int j;
	for(j=0;j<bigdim;j++) ec->clist[j]=j;
}

void ec_print(struct e_coeff * ec)
{
	int i, j, t, jr;

	printf("[ ");
	for(i = 0 ; i < m_param ; i++) {
		for(j = 0 ; j < bigdim ; j++) {
			jr=ec->clist[j];
			for(t=0;t<=ec->degree;t++) {
				if (t) printf(" + Y^%d * ",t);
				k_print(mbmat_scal(mbpoly_coeff(ec->p,t),i,jr));
			}
			if (t==0)
				printf("0");
			if (i+1<m_param || j+1<bigdim)
				printf(",	");
		}
		if (i+1<m_param)
			printf("\n  ");
	}
	printf("	]\n");
}
