divert(-1)

This file defines the necessary stuff to obtain header files for the
complicated structure definitions like foo x bar (polynomials of)
matrices over field gee, with locality emphasis on lines/columns...

dnl	First, matrices.

define(`matrix_byline_specific_defs',``
#ifdef	HARDCODE_PARAMS
typedef mp_limb_t		bw_$1mat_r[$2][$3][$5];
#define $1mat_pos_m(x,i,j,k)	x[i][j][k]
#define $1mat_scal_m(x,i,j)	((mp_limb_t *) (x[i][j]))
#else	/* HARDCODE_PARAMS */
typedef mp_limb_t		* bw_$1mat_r;
#define $1mat_pos_m(x,i,j,k)    *((x)+(k)+$5*((j)+$3*(i)))
#define $1mat_scal_m(x,i,j)     ((mp_limb_t *) ((x)+$5*((j)+$3*(i))))
#endif	/* HARDCODE_PARAMS */'')

define(`matrix_bycol_specific_defs',``
#ifdef HARDCODE_PARAMS
typedef mp_limb_t		bw_$1mat_r[$3][$2][$5];
#define $1mat_pos_m(x,i,j,k)	x[j][i][k]
#define $1mat_scal_m(x,i,j)	((mp_limb_t *) (x[j][i]))
#define $1mat_col_m(x,j)	((mp_limb_t *) (x[j]))
#else	/* HARDCODE_PARAMS */
typedef mp_limb_t		* bw_$1mat_r;
#define $1mat_pos_m(x,i,j,k)    *((x)+(k)+$5*((i)+$2*(j)))
#define $1mat_scal_m(x,i,j)     ((mp_limb_t *) ((x)+$5*((i)+$2*(j))))
#define $1mat_col_m(x,j)	((mp_limb_t *) ((x)+$5*$2*(j)))
#endif	/* HARDCODE_PARAMS */'')


define(`matrix_common_defs',``
#define $1mat_size		($2*$3*$5)
#ifdef HARDCODE_PARAMS
#define $1mat_alloc_m(x)	(x)
#define $1mat_free_m(x)
#define $1mat_head_m(x)		(&($1mat_pos_m(x,0,0,0)))
#else	/* HARDCODE_PARAMS */
#define $1mat_alloc_m(x)	(x = my_malloc($1mat_size*sizeof(mp_limb_t)))
#define $1mat_free_m(x)		free(x)
#define $1mat_head_m(x)		(x)
#endif	/* HARDCODE_PARAMS */
#define $1mat_zero_m(x)		(void) memset($1mat_head_m(x),0,$1mat_size*sizeof(mp_limb_t))
#define $1mat_copy_m(y,x)	(void) memcpy(y,x,$1mat_size*sizeof(mp_limb_t))
STRICTTYPE_DECL(bw_$1mat_r,bw_$1mat);'')

define(`matrix_ops_macros',``
#define $1mat_alloc(x)		$1mat_alloc_m(STRICTTYPE_VAL(x))
#define $1mat_free(x)		$1mat_free_m(STRICTTYPE_VAL(x))
#define $1mat_pos(x,i,j,k)	$1mat_pos_m(STRICTTYPE_VAL(x),i,j,k)
#define $1mat_scal(x,i,j)	$1mat_scal_m(STRICTTYPE_VAL(x),i,j)
#define $1mat_head(x)		$1mat_head_m(STRICTTYPE_VAL(x))
#define $1mat_zero(x)		$1mat_zero_m(STRICTTYPE_VAL(x))
#define $1mat_copy(y,x)		$1mat_copy_m(STRICTTYPE_VAL(y),STRICTTYPE_VAL(x))'')

define(`matrix_ops_inlines_protos',``
#ifdef	COMPILER_KNOWS_PRAGMA
#pragma inline(_$1mat_alloc)
#pragma inline($1mat_free)
#pragma inline($1mat_pos)
#pragma inline($1mat_scal)
#pragma inline($1mat_head)
#pragma inline($1mat_zero)
#endif	/* COMPILER_KNOWS_PRAGMA */

#define $1mat_alloc(x)		_$1mat_alloc(&(x))
FASTFUNC void*			_$1mat_alloc(bw_$1mat *);
FASTFUNC void			$1mat_free(bw_$1mat);
FASTFUNC mp_limb_t		$1mat_pos(bw_$1mat,  int, int, int);
FASTFUNC mp_limb_t *		$1mat_scal(bw_$1mat, int, int);
FASTFUNC void *			$1mat_head(bw_$1mat);
FASTFUNC void			$1mat_zero(bw_$1mat);'')

define(`matrix_ops_inlines_implementations',``
FASTFUNC void* _$1mat_alloc(bw_$1mat * px)
{
	return $1mat_alloc_m(STRICTTYPE_VAL(*px));
}
FASTFUNC void	$1mat_free(bw_$1mat x)
{
	$1mat_free_m(STRICTTYPE_VAL(x));
}
FASTFUNC mp_limb_t $1mat_pos(bw_$1mat x, int i, int j, int k)
{
	return $1mat_pos_m(STRICTTYPE_VAL(x),i,j,k);
}
FASTFUNC mp_limb_t * $1mat_scal(bw_$1mat x, int i, int j)
{
	return $1mat_scal_m(STRICTTYPE_VAL(x),i,j);
}
FASTFUNC void * $1mat_head(bw_$1mat x)
{
	return $1mat_head_m(STRICTTYPE_VAL(x));
}
FASTFUNC void $1mat_zero(bw_$1mat x)
{
	$1mat_zero_m(STRICTTYPE_VAL(x));
}
FASTFUNC void $1mat_copy(bw_$1mat y, bw_$1mat x)
{
	$1mat_copy_m(STRICTTYPE_VAL(y),STRICTTYPE_VAL(x));
}'')

dnl For column-oriented matrices, define some column opoerations.

define(`col_defs',``
#define $4col_size			($2*$5)
#define	$4col_copy_m(y,x)	(void) memcpy(y,x,$4col_size*sizeof(mp_limb_t))
#define	$4col_zero_m(x)		(void) memset(x,0,$4col_size*sizeof(mp_limb_t))
#define	$4col_scal_m(c,i)	((mp_limb_t *) ((c)+i*$5))
#define $4col_alloc_m(x)	(x=(bw_$4col_r) my_malloc($4col_size*sizeof(mp_limb_t)))
#define $4col_free_m(x)		free(x)
typedef mp_limb_t * bw_$4col_r;
STRICTTYPE_DECL(bw_$4col_r,bw_$4col);'')

define(`col_ops_macros',``
#define $1mat_col(x,j)		STRICTTYPE_CAST(bw_$4col,$1mat_col_m(STRICTTYPE_VAL(x),j))
#define $4col_copy(y,x)		$4col_copy_m(STRICTTYPE_VAL(y),STRICTTYPE_VAL(x))
#define $4col_zero(x)		$4col_zero_m(STRICTTYPE_VAL(x))
#define $4col_scal(c,i)		$4col_scal_m(STRICTTYPE_VAL(c),i)
#define $4col_alloc(x)		$4col_alloc_m(STRICTTYPE_VAL(x))
#define $4col_free(x)		$4col_free_m(STRICTTYPE_VAL(x))'')

define(`col_ops_funcs_protos',``
FASTFUNC int $4col_is_zero(bw_$4col);'')

define(`col_ops_inlines_protos',``
#ifdef	COMPILER_KNOWS_PRAGMA
#pragma inline($1mat_col)
#pragma inline($4col_copy)
#pragma inline($4col_zero)
#pragma inline($4col_scal)
#pragma inline(_$4col_alloc)
#pragma inline($4col_free)
#endif	/* COMPILER_KNOWS_PRAGMA */

#define	$4col_alloc(x)		_$4col_alloc(&(x))
FASTFUNC bw_$4col $1mat_col(bw_$1mat, int);
FASTFUNC void $4col_copy(bw_$4col, bw_$4col);
FASTFUNC void $4col_zero(bw_$4col);
FASTFUNC mp_limb_t * $4col_scal(bw_$4col, int);
FASTFUNC void * _$4col_alloc(bw_$4col *);
FASTFUNC void $4col_free(bw_$4col);'')

define(`col_ops_inlines_implementations',``
FASTFUNC bw_$4col $1mat_col(bw_$1mat x, int j)
{
	return STRICTTYPE_CAST(bw_$4col,$1mat_col_m(STRICTTYPE_VAL(x),j));
}
FASTFUNC void $4col_copy(bw_$4col y, bw_$4col x)
{
	$4col_copy_m(STRICTTYPE_VAL(y),STRICTTYPE_VAL(x));
}
FASTFUNC void $4col_zero(bw_$4col x)
{
	$4col_zero_m(STRICTTYPE_VAL(x));
}
FASTFUNC bw_scalar $4col_scal(bw_$4col c, int i)
{
	return $4col_scal_m(STRICTTYPE_VAL(c),i);
}
FASTFUNC void * _$4col_alloc(bw_$4col *c)
{
	return $4col_alloc_m(STRICTTYPE_VAL(*c));
}
FASTFUNC void $4col_free(bw_$4col c)
{
	$4col_free_m(STRICTTYPE_VAL(c));
}'')

define(`col_ops_funcs_implementations',``
FASTFUNC int $4col_is_zero(bw_$4col c)
{
	int i;
	for(i=0;i<$2;i++) {
		if (!bw_scalar_is_zero($4col_scal(c,i))) {
			return 0;
		}
	}
	return 1;
}'')


dnl	Now, implement polynomials

poly	: ARG 1 : name of the structure (suffixed by ``poly'')
	  ARG 2	: full reference type
	  ARG 3	: reference type size, in limbs

The reference type must admit a raw, _r suffixed version. typedef it if
it doesn't exist. The reference type must also be a pointer type.
integral types won't work.

define(`poly_defs',``
#ifdef	HARDCODE_PARAMS
typedef	$2_r		      * bw_$1poly_r;
#define $1poly_coeff_m(x,t)	x[t]
#else	/* HARDCODE_PARAMS */
typedef $2_r			bw_$1poly_r;
#define $1poly_coeff_m(x,t)	((x)+(t)*($3))
#endif	/* HARDCODE_PARAMS */
#define $1poly_alloc_m(x,d)	(x=(bw_$1poly_r) my_malloc($3 * ((d)+1) * sizeof(mp_limb_t)))
#define $1poly_zero_m(x,d)	(void) memset(x,0,$3 * ((d)+1) * sizeof(mp_limb_t))
#define $1poly_free_m(x)	free(x)
STRICTTYPE_DECL(bw_$1poly_r, bw_$1poly);'')

define(`poly_ops_macros',``
#define $1poly_coeff(x,t)	STRICTTYPE_CAST($2,$1poly_coeff_m(STRICTTYPE_VAL(x),t))
#define $1poly_subpoly(x,t)	STRICTTYPE_CAST(bw_$1poly,$1poly_coeff_m(STRICTTYPE_VAL(x),t))
#define $1poly_alloc(x,d)	$1poly_alloc_m(STRICTTYPE_VAL(x),d)
#define $1poly_zero(x,t)	$1poly_zero_m(STRICTTYPE_VAL(x),t)
#define $1poly_free(x)		$1poly_free_m(STRICTTYPE_VAL(x))'')

define(`poly_ops_inlines_protos',``
#ifdef	COMPILER_KNOWS_PRAGMA
#pragma inline($1poly_coeff)
#pragma inline($1poly_subpoly)
#pragma inline(_$1poly_alloc)
#pragma inline($1poly_zero)
#pragma inline($1poly_free)
#endif	/* COMPILER_KNOWS_PRAGMA */

#define $1poly_alloc(x,d)	_$1poly_alloc(&(x),d)
FASTFUNC $2			$1poly_coeff(bw_$1poly, int);
FASTFUNC bw_$1poly		$1poly_subpoly(bw_$1poly x, int t);
FASTFUNC void *			_$1poly_alloc(bw_$1poly *, int);
FASTFUNC void			$1poly_zero(bw_$1poly, int);
FASTFUNC void			$1poly_free(bw_$1poly);'')

define(`poly_ops_inlines_implementations',``
FASTFUNC $2 $1poly_coeff(bw_$1poly x, int t)
{
	return STRICTTYPE_CAST($2,$1poly_coeff_m(STRICTTYPE_VAL(x),t));
}
FASTFUNC bw_$1poly $1poly_subpoly(bw_$1poly x, int t)
{
	return STRICTTYPE_CAST(bw_$1poly,$1poly_coeff_m(STRICTTYPE_VAL(x),t));
}
FASTFUNC void * _$1poly_alloc(bw_$1poly * px, int d)
{
	return $1poly_alloc_m(STRICTTYPE_VAL(*px),d);
}
FASTFUNC void $1poly_zero(bw_$1poly x, int t)
{
	$1poly_zero_m(STRICTTYPE_VAL(x),t);
}
FASTFUNC void $1poly_free(bw_$1poly x)
{
	$1poly_free_m(STRICTTYPE_VAL(x));
}'')


dnl
dnl	Now, DFTs
dnl

dnl $1 name
dnl $2 lines
dnl $3 columns
dnl $4 scalar size, in limbs

define(`dft_defs',``
typedef mp_limb_t		* 	bw_$1dft_r;
#define $1dft_poly_m(s,x,i,j)	((x)+(($4)<<(s))*((i)+$2*(j)))
#define $1dft_scal_m(s,x,i,j,k)	($1dft_poly_m((s),(x),(i),(j))+(($4)*k))
#define $1dft_size(s)		(($2*$3*$4)<<(s))
#define $1dft_alloc_m(s,x)	(x = my_malloc($1dft_size((s))*sizeof(mp_limb_t)))
#define $1dft_free_m(s,x)	free(x)
#define $1dft_head_m(s,x)	(x)
#define $1dft_zero_m(s,x)	(void) memset($1dft_head_m((s),(x)),0,$1dft_size((s))*sizeof(mp_limb_t))
#define $1dft_copy_m(s,y,x)	(void) memcpy((y),(x),$1dft_size((s))*sizeof(mp_limb_t))
STRICTTYPE_DECL(bw_$1dft_r,bw_$1dft);'')

define(`dft_ops_macros',``
#define $1dft_alloc(s,x)	$1dft_alloc_m((s),STRICTTYPE_VAL(x))
#define $1dft_free(s,x)		$1dft_free_m((s),STRICTTYPE_VAL(x))
#define $1dft_poly(s,x,i,j)	$1dft_poly_m((s),STRICTTYPE_VAL(x),i,j)
#define $1dft_scal(s,x,i,j,k)	$1dft_scal_m((s),STRICTTYPE_VAL(x),i,j,k)
#define $1dft_head(s,x)		$1dft_head_m((s),STRICTTYPE_VAL(x))
#define $1dft_zero(s,x)		$1dft_zero_m((s),STRICTTYPE_VAL(x))
#define $1dft_copy(s,y,x)	$1dft_copy_m((s),STRICTTYPE_VAL(y),STRICTTYPE_VAL(x))'')

define(`dft_ops_inlines_protos',``
#ifdef	COMPILER_KNOWS_PRAGMA
#pragma inline(_$1dft_alloc)
#pragma inline($1dft_free)
#pragma inline($1dft_poly)
#pragma inline($1dft_scal)
#pragma inline($1dft_head)
#pragma inline($1dft_zero)
#pragma inline($1dft_copy)
#endif	/* COMPILER_KNOWS_PRAGMA */

#define $1dft_alloc(s, x)	_$1dft_alloc((s), &(x))
FASTFUNC void*			_$1dft_alloc(int, bw_$1dft *);
FASTFUNC void			$1dft_free(int, bw_$1dft);
FASTFUNC mp_limb_t *		$1dft_poly(int, bw_$1dft, int, int);
FASTFUNC mp_limb_t *		$1dft_scal(int, bw_$1dft, int, int, int);
FASTFUNC void *			$1dft_head(int, bw_$1dft);
FASTFUNC void			$1dft_zero(int, bw_$1dft);'')

define(`dft_ops_inlines_implementations',``
FASTFUNC void* _$1dft_alloc(int s,bw_$1dft * px)
{
	return $1dft_alloc_m(s,STRICTTYPE_VAL(*px));
}
FASTFUNC void	$1dft_free(int s,bw_$1dft x)
{
	$1dft_free_m(s,STRICTTYPE_VAL(x));
}
FASTFUNC mp_limb_t * $1dft_poly(int s,bw_$1dft x, int i, int j)
{
	return $1dft_poly_m(s,STRICTTYPE_VAL(x),i,j);
}
FASTFUNC mp_limb_t * $1dft_scal(int s,bw_$1dft x, int i, int j, int k)
{
	return $1dft_scal_m(s,STRICTTYPE_VAL(x),i,j,k);
}
FASTFUNC void * $1dft_head(int s,bw_$1dft x)
{
	return $1dft_head_m(s,STRICTTYPE_VAL(x));
}
FASTFUNC void $1dft_zero(int s,bw_$1dft x)
{
	$1dft_zero_m(s,STRICTTYPE_VAL(x));
}
FASTFUNC void $1dft_copy(int s,bw_$1dft y, bw_$1dft x)
{
	$1dft_copy_m(s,STRICTTYPE_VAL(y),STRICTTYPE_VAL(x));
}'')

dnl
dnl *******************************************************************
dnl

matrix_byline :	ARG 1 : name of the structure (suffixed by ``mat'')
		ARG 2 : number of lines
		ARG 3 : number of columns
		ARG 4 : unused
		ARG 5 : field size, in limbs

define(`matrix_byline_header',`
/*****************************************************************************/
/* Here are the definitions for bw_$1mat : $2 x $3 matrices, line local */
/*****************************************************************************/
matrix_byline_specific_defs($@)
matrix_common_defs($@)

#ifdef	PREFER_INLINES
matrix_ops_inlines_protos($@)

#else	/* PREFER_INLINES */
matrix_ops_macros($@)

#endif	/* PREFER_INLINES */
')

matrix_bycol :	ARG 1 : name of the structure (suffixed by ``mat'')
		ARG 2 : number of lines
		ARG 3 : number of columns
		ARG 4 : name of the column type (suffixed by ``col'')
		ARG 5 : field size, in limbs

define(`matrix_bycol_header',`
/*****************************************************************************/
/* Here are the definitions for bw_$1mat : $2 x $3 matrices, column local */
/*****************************************************************************/
matrix_bycol_specific_defs($@)
matrix_common_defs($@)
col_defs($@)

#ifdef	PREFER_INLINES
matrix_ops_inlines_protos($@)
col_ops_inlines_protos($@)

#else	/* PREFER_INLINES */
matrix_ops_macros($@)
col_ops_macros($@)

#endif	/* PREFER_INLINES */
col_ops_funcs_protos($@)
')

define(`matrix_companion_header',`
/*****************************************************************************/
/* Here are the implementations for bw_$1mat : $2 x $3 matrices */
/*****************************************************************************/

#ifdef IMPLEMENT_INLINES
matrix_ops_inlines_implementations($@)

#endif	/* IMPLEMENT_INLINES */

ifelse(`$4',`',,`

#ifdef IMPLEMENT_INLINES
col_ops_inlines_implementations($@)

#endif	/* IMPLEMENT_INLINES */

#if defined(IMPLEMENT_FUNCS) || !defined(PREFER_INLINES)
col_ops_funcs_implementations($@)

#endif	/* IMPLEMENT_FUNCS || !PREFER_INLINES */
')')

define(`poly_header',`
/*****************************************************************************/
/* Here are the definitions for bw_$1poly : $2 - polynomials */
/*****************************************************************************/
poly_defs($@)

#ifdef PREFER_INLINES
poly_ops_inlines_protos($@)

#else	/* PREFER_INLINES */
poly_ops_macros($@)

#endif	/* PREFER_INLINES */
')

define(`poly_companion_header',`
/*****************************************************************************/
/* Here are the implementations for bw_$1poly : $2 - polynomials */
/*****************************************************************************/

#ifdef IMPLEMENT_INLINES
poly_ops_inlines_implementations($@)

#endif	/* IMPLEMENT_INLINES */
')

define(`dft_header',`
/*****************************************************************************/
/* Here are the definitions for bw_$1dft : $2 x $3 - DFTs */
/*****************************************************************************/
dft_defs($@)

#ifdef PREFER_INLINES
dft_ops_inlines_protos($@)

#else	/* PREFER_INLINES */
dft_ops_macros($@)

#endif	/* PREFER_INLINES */
')

define(`dft_companion_header',`
/*****************************************************************************/
/* Here are the implementations for bw_$1dft : $2 x $3 - DFTs */
/*****************************************************************************/

#ifdef IMPLEMENT_INLINES
dft_ops_inlines_implementations($@)

#endif	/* IMPLEMENT_INLINES */
')


divert(0)dnl
dnl
dnl
dnl structure_automatic.h
dnl
ifelse(outputfile,`MAIN',`
/* This file was generated by m4 on esyscmd(echo -n $(date)) */
/* Do not edit this file. Re-generate instead from the source (using m4) */

#ifndef INCLUDING_STRUCTURE_AUTOMATIC_H_
#error "Do not include this file directly\n"
#endif

#ifndef STRUCTURE_AUTOMATIC_H_
#define STRUCTURE_AUTOMATIC_H_

#ifdef __cplusplus
extern "C" {
#endif
matrix_byline_header(mn,	m_param,	n_param,,	bw_allocsize)
matrix_bycol_header(nb,		n_param,	bigdim,	n,	bw_allocsize)
matrix_bycol_header(mb,		m_param,	bigdim,	m,	bw_allocsize)
matrix_bycol_header(bb,		bigdim,		bigdim,	b,	bw_allocsize)
poly_header(mn,			bw_mnmat,	mnmat_size)
poly_header(nb,			bw_nbmat,	nbmat_size)
poly_header(mb,			bw_mbmat,	mbmat_size)
poly_header(bb,			bw_bbmat,	bbmat_size)
#ifndef EXTFIELD_SIZE
#ifdef HAS_NATIVE_FFT
#define EXTFIELD_SIZE	bw_allocsize
#else	/* HAS_NATIVE_FFT */
#define EXTFIELD_SIZE   (bw_allocsize<<1)
#endif	/* HAS_NATIVE_FFT */
#endif
matrix_byline_header(x_mn,	m_param,	n_param,,	EXTFIELD_SIZE)
matrix_bycol_header(x_nb,	n_param,	bigdim,	x_n,	EXTFIELD_SIZE)
matrix_bycol_header(x_mb,	m_param,	bigdim,	x_m,	EXTFIELD_SIZE)
matrix_bycol_header(x_bb,	bigdim,		bigdim,	x_b,	EXTFIELD_SIZE)
poly_header(x_mn,		bw_x_mnmat,	x_mnmat_size)
poly_header(x_nb,		bw_x_nbmat,	x_nbmat_size)
poly_header(x_mb,		bw_x_mbmat,	x_mbmat_size)
poly_header(x_bb,		bw_x_bbmat,	x_bbmat_size)
dft_header(mn,	m_param,	n_param, EXTFIELD_SIZE)
dft_header(nb,	n_param,	bigdim, EXTFIELD_SIZE)
dft_header(mb,	m_param,	bigdim, EXTFIELD_SIZE)
dft_header(bb,	bigdim,		bigdim, EXTFIELD_SIZE)
#ifdef __cplusplus
}
#endif

#endif	/* STRUCTURE_AUTOMATIC_H_ */')
dnl
dnl
dnl structure_inlines_automatic.h
dnl
dnl
ifelse(outputfile,`COMPANION',`dnl
/* This file was generated by m4 on esyscmd(echo -n $(date)) */
/* Do not edit this file. Re-generate instead from the source (using m4) */

#ifndef INCLUDING_STRUCTURE_INLINES_AUTOMATIC_H_
#error "Do not include this file directly\n"
#endif

#ifndef STRUCTURE_INLINES_AUTOMATIC_H_
#define STRUCTURE_INLINES_AUTOMATIC_H_

#ifdef __cplusplus
extern "C" {
#endif
matrix_companion_header(mn,	m_param,	n_param,,	bw_allocsize)
matrix_companion_header(nb,	n_param,	bigdim,	n,	bw_allocsize)
matrix_companion_header(mb,	m_param,	bigdim,	m,	bw_allocsize)
matrix_companion_header(bb,	bigdim,		bigdim,	b,	bw_allocsize)
poly_companion_header(mn,	bw_mnmat,	mnmat_size)
poly_companion_header(nb,	bw_nbmat,	nbmat_size)
poly_companion_header(mb,	bw_mbmat,	mbmat_size)
poly_companion_header(bb,	bw_bbmat,	bbmat_size)
#ifndef EXTFIELD_SIZE
#ifdef HAS_NATIVE_FFT
#define EXTFIELD_SIZE	bw_allocsize
#else	/* HAS_NATIVE_FFT */
#define EXTFIELD_SIZE   (bw_allocsize<<1)
#endif	/* HAS_NATIVE_FFT */
#endif
matrix_companion_header(x_mn,	m_param,	n_param,,	EXTFIELD_SIZE)
matrix_companion_header(x_nb,	n_param,	bigdim,	x_n,	EXTFIELD_SIZE)
matrix_companion_header(x_mb,	m_param,	bigdim,	x_m,	EXTFIELD_SIZE)
matrix_companion_header(x_bb,	bigdim,		bigdim,	x_b,	EXTFIELD_SIZE)
poly_companion_header(x_mn,	bw_x_mnmat,	x_mnmat_size)
poly_companion_header(x_nb,	bw_x_nbmat,	x_nbmat_size)
poly_companion_header(x_mb,	bw_x_mbmat,	x_mbmat_size)
poly_companion_header(x_bb,	bw_x_bbmat,	x_bbmat_size)
dft_companion_header(mn,	m_param,	n_param,	EXTFIELD_SIZE)
dft_companion_header(nb,	n_param,	bigdim,		EXTFIELD_SIZE)
dft_companion_header(mb,	m_param,	bigdim, 	EXTFIELD_SIZE)
dft_companion_header(bb,	bigdim,		bigdim, 	EXTFIELD_SIZE)
#ifdef __cplusplus
}
#endif

#endif	/* STRUCTURE_INLINES_AUTOMATIC_H_ */')
