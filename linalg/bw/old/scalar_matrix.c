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
#include "scalar_matrix.h"

/* Operation : C_j <- C_{\sigma(j)} */
void
ec_apply_perm(struct e_coeff * ec, unsigned int * sigma)
{
	unsigned int * tmp;
	unsigned int i;

	tmp=malloc(bigdim * sizeof(unsigned int));
	for(i=0;i<bigdim;i++) tmp[i]=ec->clist[sigma[i]];
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
	int i;
	int j1,j2;

	j1=ec->clist[o_j1];
	j2=ec->clist[o_j2];


	bw_scalar_set_zero(mbmat_scal(ec->mat,ikill,j2));
	

	for(i=ikill+1;i<m_param;i++) {
		addmul( mbmat_scal(ec->mat,i,j2),
			mbmat_scal(ec->mat,i,j1),
			lambda);
	}
}

void
ec_x_multiply(struct e_coeff * ec UNUSED_VARIABLE, int o_j UNUSED_VARIABLE)
{
#ifndef NDEBUG
	int j;

	j=ec->clist[o_j];
	mcol_zero(mbmat_col(ec->mat,j));
#endif	/* NDEBUG */
}

	/* We zero it out because our computations are done ``modulo X''
	 * at this point.
	 *
	 * This has no real importance however, so we might as well make
	 * this function a no-op to save (very) few cycles. Hence the
	 * #ifdef / #endif pair....
	 */
