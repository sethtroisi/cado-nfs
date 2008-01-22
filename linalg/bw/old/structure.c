#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "lingen_params.h"
#include "params.h"
#include <assert.h>
#include "types.h"
#include "macros.h"
#include "bw_scalar.h"
#include "variables.h"
#include "auxfuncs.h"
#include "structure.h"
#include "gmp-hacks.h"

#ifdef __GNUC__
#undef  FASTFUNC
#define FASTFUNC
#define IMPLEMENT_INLINES
#define IMPLEMENT_FUNCS
#else
#ifndef INCLUDING_STRUCTURE_INLINES_H_
#undef	FASTFUNC
#define FASTFUNC
#define IMPLEMENT_INLINES
#define IMPLEMENT_FUNCS
#else
#undef	IMPLEMENT_INLINES
#undef	IMPLEMENT_FUNCS
#endif
#endif

#undef STRUCTURE_INLINES_AUTOMATIC_H_
#define INCLUDING_STRUCTURE_INLINES_AUTOMATIC_H_

#include "structure_inlines_automatic.h"


int nbpoly_write(FILE *f, bw_nbpoly p, int deg)
{
	int i,j,k;
	mpz_t blah;
	mpz_init(blah);
	
	for(k=0;k<=deg;k++) {
		for(i=0;i<n_param;i++) for(j=0;j<bigdim;j++) {
			MPZ_SET_MPN(blah, nbmat_scal(nbpoly_coeff(p,k),i,j), bw_allocsize);
			gmp_fprintf(f, "%Zd%c", blah, (j==bigdim-1)?'\n':' ');
		}
	}

	mpz_clear(blah);
	return 0;
}

bw_nbpoly nbpoly_read(FILE *f, int deg)
{
	int i,j,k;
	bw_nbpoly res;
	mpz_t blah;
	mpz_init(blah);

	nbpoly_alloc(res,deg);
	for(k=0;k<=deg;k++) {
		int rc = 0;
		for(i=0;i<n_param;i++) for(j=0;j<bigdim;j++) {
			if (gmp_fscanf(f, "%Zd", blah) != 1) 
				continue;
			rc++;
			MPN_SET_MPZ(nbmat_scal(nbpoly_coeff(res,k),i,j), bw_allocsize, blah);
		}
		if (rc == 0) {
			die("degree is too small (%d < %d)\n", 1, k-1, deg);
		} else if (rc != n_param * bigdim) {
			die("file not in sync\n", 1);
			/* XXX speak something more telling... */
		}
	}

	return res;
}
#if 0
#include <sys/types.h>
#include <sys/stat.h>
int nbpoly_write(FILE *f, bw_nbpoly p, int deg)
{
	int i,j,k,res;
	
	for(j=0;j<bigdim;j++)
	for(i=0;i<n_param;i++)
	for(k=0;k<=deg;k++) {
		res=bw_scalar_write(f,nbmat_scal(nbpoly_coeff(p,k),i,j));
		if (res != bw_filesize)
			return -1;
	}

	return 0;
}

bw_nbpoly nbpoly_read(int *deg, FILE *f)
{
	int i,j,k,n;
	struct stat sbuf;
	bw_nbpoly res;

	if (fstat(fileno(f),&sbuf)<0)
		return STRICTTYPE_CAST(bw_nbpoly,NULL);
	*deg=((sbuf.st_size)/((nbmat_size/bw_allocsize)*bw_filesize*
				sizeof(mp_limb_t)))-1;
	if (*deg==-1) {
		return STRICTTYPE_CAST(bw_nbpoly,NULL);
	}
	nbpoly_alloc(res,*deg);
	
	for(j=0;j<bigdim;j++)
	for(i=0;i<n_param;i++)
	for(k=0;k<=*deg;k++) {
		n=bw_scalar_read(nbmat_scal(nbpoly_coeff(res,k),i,j),f);
		if (n != bw_filesize) {
			nbpoly_free(res);
			return STRICTTYPE_CAST(bw_nbpoly,NULL);
		}
	}

	return res;
}
#endif
