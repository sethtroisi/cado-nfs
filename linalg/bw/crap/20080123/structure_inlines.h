#ifndef INCLUDING_STRUCTURE_INLINES_H_
#error "This file is not intended to be included directly\n"
#endif
#define INCLUDED_STRUCTURE_INLINES_H_

#ifdef	__cplusplus
extern "C" {
#endif

FASTFUNC void* _mnmat_alloc(bw_mnmat * px)
{
	return mnmat_alloc_m(STRICTTYPE_VAL(*px));
}
FASTFUNC void  mnmat_free(bw_mnmat x)
{
	mnmat_free_m(STRICTTYPE_VAL(x));
}
FASTFUNC mp_limb_t mnmat_pos(bw_mnmat x, int i, int j, int k)
{
	return mnmat_pos_m(STRICTTYPE_VAL(x),i,j,k);
}
FASTFUNC bw_scalar mnmat_scal(bw_mnmat x, int i, int j)
{
	return mnmat_scal_m(STRICTTYPE_VAL(x),i,j);
}
FASTFUNC void * mnmat_head(bw_mnmat x)
{
	return mnmat_head_m(STRICTTYPE_VAL(x));
}
FASTFUNC mp_limb_t * mnmat_step(mp_limb_t * x, int i, int j, int k)
{
	return mnmat_step_m(x,i,j,k);
}
FASTFUNC void mnmat_zero(bw_mnmat x)
{
	mnmat_zero_m(STRICTTYPE_VAL(x));
}
FASTFUNC bw_mnmat mnpoly_coeff(bw_mnpoly x, int t)
{
	return STRICTTYPE_CAST(bw_mnmat,mnpoly_coeff_m(STRICTTYPE_VAL(x),t));
}
FASTFUNC void * _mnpoly_alloc(bw_mnpoly * px, int d)
{
	return mnpoly_alloc_m(STRICTTYPE_VAL(*px),d);
}
FASTFUNC void mnpoly_zero(bw_mnpoly x, int t)
{
	mnpoly_zero_m(STRICTTYPE_VAL(x),t);
}
FASTFUNC void mnpoly_free(bw_mnpoly x)
{
	mnpoly_free_m(STRICTTYPE_VAL(x));
}
FASTFUNC void * _nbmat_alloc(bw_nbmat * px)
{
	return nbmat_alloc_m(STRICTTYPE_VAL(*px));
}
FASTFUNC void nbmat_free(bw_nbmat x)
{
	nbmat_free_m(STRICTTYPE_VAL(x));
}
FASTFUNC mp_limb_t nbmat_pos(bw_nbmat x, int i, int j, int k)
{
	return nbmat_pos_m(STRICTTYPE_VAL(x),i,j,k);
}
FASTFUNC bw_scalar nbmat_scal(bw_nbmat x, int i, int j)
{
	return nbmat_scal_m(STRICTTYPE_VAL(x),i,j);
}
FASTFUNC bw_ncol nbmat_col(bw_nbmat x, int j)
{
	return STRICTTYPE_CAST(bw_ncol,nbmat_col_m(STRICTTYPE_VAL(x),j));
}
FASTFUNC void * nbmat_head(bw_nbmat x)
{
	return nbmat_head_m(STRICTTYPE_VAL(x));
}
FASTFUNC mp_limb_t * nbmat_step(mp_limb_t * x, int i, int j, int k)
{
	return nbmat_step_m(x,i,j,k);
}
FASTFUNC void nbmat_zero(bw_nbmat x)
{
	nbmat_zero_m(STRICTTYPE_VAL(x));
}
FASTFUNC void ncol_copy(bw_ncol y, bw_ncol x)
{
	ncol_copy_m(STRICTTYPE_VAL(y),STRICTTYPE_VAL(x));
}
FASTFUNC void ncol_zero(bw_ncol x)
{
	ncol_zero_m(STRICTTYPE_VAL(x));
}
FASTFUNC bw_scalar ncol_scal(bw_ncol c, int i)
{
	return ncol_scal_m(STRICTTYPE_VAL(c),i);
}
FASTFUNC bw_nbmat nbpoly_coeff(bw_nbpoly x, int t)
{
	return STRICTTYPE_CAST(bw_nbmat,nbpoly_coeff_m(STRICTTYPE_VAL(x),t));
}
FASTFUNC void * _nbpoly_alloc(bw_nbpoly * px, int d)
{
	return nbpoly_alloc_m(STRICTTYPE_VAL(*px),d);
}
FASTFUNC void nbpoly_zero(bw_nbpoly x, int t)
{
	nbpoly_zero_m(STRICTTYPE_VAL(x),t);
}
FASTFUNC void nbpoly_free(bw_nbpoly x)
{
	nbpoly_free_m(STRICTTYPE_VAL(x));
}
FASTFUNC int ncol_is_zero(bw_ncol c)
{
	int i;
	for(i=0;i<n_param;i++) {
		if (!bw_scalar_is_zero(ncol_scal(c,i))) {
			return 0;
		}
	}
	return 1;
}
FASTFUNC void * _mbmat_alloc(bw_mbmat * px)
{
	return mbmat_alloc_m(STRICTTYPE_VAL(*px));
}
FASTFUNC void mbmat_free(bw_mbmat x)
{
	mbmat_free_m(STRICTTYPE_VAL(x));
}
FASTFUNC mp_limb_t mbmat_pos(bw_mbmat x, int i, int j, int k)
{
	return mbmat_pos_m(STRICTTYPE_VAL(x),i,j,k);
}
FASTFUNC bw_scalar mbmat_scal(bw_mbmat x, int i, int j)
{
	return mbmat_scal_m(STRICTTYPE_VAL(x),i,j);
}
FASTFUNC bw_mcol mbmat_col(bw_mbmat x, int j)
{
	return STRICTTYPE_CAST(bw_mcol,mbmat_col_m(STRICTTYPE_VAL(x),j));
}
FASTFUNC void * mbmat_head(bw_mbmat x)
{
	return mbmat_head_m(STRICTTYPE_VAL(x));
}
FASTFUNC void mbmat_zero(bw_mbmat x)
{
	mbmat_zero_m(STRICTTYPE_VAL(x));
}
FASTFUNC void mbmat_copy(bw_mbmat y, bw_mbmat x)
{
	mbmat_copy_m(STRICTTYPE_VAL(y),STRICTTYPE_VAL(x));
}
FASTFUNC void mcol_copy(bw_mcol y, bw_mcol x)
{
	mcol_copy_m(STRICTTYPE_VAL(y),STRICTTYPE_VAL(x));
}
FASTFUNC void mcol_zero(bw_mcol x)
{
	mcol_zero_m(STRICTTYPE_VAL(x));
}
FASTFUNC bw_scalar mcol_scal(bw_mcol c, int i)
{
	return mcol_scal_m(STRICTTYPE_VAL(c),i);
}
FASTFUNC bw_mbmat mbpoly_coeff(bw_mbpoly x, int t)
{
	return STRICTTYPE_CAST(bw_mbmat,mbpoly_coeff_m(STRICTTYPE_VAL(x),t));
}
FASTFUNC bw_mbpoly mbpoly_subpoly(bw_mbpoly x, int t)
{
	return STRICTTYPE_CAST(bw_mbpoly,mbpoly_coeff_m(STRICTTYPE_VAL(x),t));
}
FASTFUNC void * _mbpoly_alloc(bw_mbpoly * px, int d)
{
	return mbpoly_alloc_m(STRICTTYPE_VAL(*px),d);
}
FASTFUNC void mbpoly_zero(bw_mbpoly x, int t)
{
	mbpoly_zero_m(STRICTTYPE_VAL(x),t);
}
FASTFUNC void mbpoly_free(bw_mbpoly x)
{
	mbpoly_free_m(STRICTTYPE_VAL(x));
}
FASTFUNC int mcol_is_zero(bw_mcol c)
{
	int i;
	for(i=0;i<m_param;i++) {
		if (!bw_scalar_is_zero(mcol_scal(c,i))) {
			return 0;
		}
	}
	return 1;
}
FASTFUNC void * _bbmat_alloc(bw_bbmat * px)
{
	return bbmat_alloc_m(STRICTTYPE_VAL(*px));
}
FASTFUNC void bbmat_free(bw_bbmat x)
{
	bbmat_free_m(STRICTTYPE_VAL(x));
}
FASTFUNC mp_limb_t bbmat_pos(bw_bbmat x, int i, int j, int k)
{
	return bbmat_pos_m(STRICTTYPE_VAL(x),i,j,k);
}
FASTFUNC bw_scalar bbmat_scal(bw_bbmat x, int i, int j)
{
	return bbmat_scal_m(STRICTTYPE_VAL(x),i,j);
}
FASTFUNC bw_bcol bbmat_col(bw_bbmat x, int j)
{
	return STRICTTYPE_CAST(bw_bcol,bbmat_col_m(STRICTTYPE_VAL(x),j));
}
FASTFUNC void * bbmat_head(bw_bbmat x)
{
	return bbmat_head_m(STRICTTYPE_VAL(x));
}
FASTFUNC mp_limb_t * bbmat_step(mp_limb_t * x, int i, int j, int k)
{
	return bbmat_step_m(x,i,j,k);
}
FASTFUNC void bbmat_zero(bw_bbmat x)
{
	bbmat_zero_m(STRICTTYPE_VAL(x));
}
FASTFUNC void bbmat_copy(bw_bbmat y, bw_bbmat x)
{
	bbmat_copy_m(STRICTTYPE_VAL(y),STRICTTYPE_VAL(x));
}
FASTFUNC void bcol_copy(bw_bcol y, bw_bcol x)
{
	bcol_copy_m(STRICTTYPE_VAL(y),STRICTTYPE_VAL(x));
}
FASTFUNC void bcol_zero(bw_bcol x)
{
	bcol_zero_m(STRICTTYPE_VAL(x));
}
FASTFUNC bw_scalar bcol_scal(bw_bcol c, int i)
{
	return bcol_scal_m(STRICTTYPE_VAL(c),i);
}
FASTFUNC bw_bbmat bbpoly_coeff(bw_bbpoly x, int t)
{
	return STRICTTYPE_CAST(bw_bbmat,bbpoly_coeff_m(STRICTTYPE_VAL(x),t));
}
FASTFUNC void * _bbpoly_alloc(bw_bbpoly * px, int d)
{
	return bbpoly_alloc_m(STRICTTYPE_VAL(*px),d);
}
FASTFUNC void bbpoly_zero(bw_bbpoly x, int t)
{
	bbpoly_zero_m(STRICTTYPE_VAL(x),t);
}
FASTFUNC void bbpoly_free(bw_bbpoly x)
{
	bbpoly_free_m(STRICTTYPE_VAL(x));
}
FASTFUNC int bcol_is_zero(bw_bcol c)
{
	int i;
	for(i=0;i<bigdim;i++) {
		if (!bw_scalar_is_zero(bcol_scal(c,i))) {
			return 0;
		}
	}
	return 1;
}

#ifdef	__cplusplus
}
#endif
