#ifndef OPS_POLY_FFT_H_
#define OPS_POLY_FFT_H_

#include <gmp.h>
#include "ops_poly.h"
#include "structure.h"
#include "auxfuncs.h"

#define FTF(X)  ft_ ## X ## _fft

#ifdef __cplusplus
extern "C" {
#endif

#define ft_mat_bb FTF(mat_bb)
#define ft_mat_mb FTF(mat_mb)
#define ft_order_t FTF(order_t)
#define ft_t FTF(t)
#define ft_order_eq FTF(order_eq)
#define ft_order_fits FTF(order_fits)
#define ft_order_int FTF(order_int)
#define ft_order_copy FTF(order_copy)

typedef bw_bbdft ft_mat_bb;
typedef bw_mbdft ft_mat_mb;
typedef mp_limb_t * ft_t;
typedef int ft_order_t;

static inline int ft_order_eq(ft_order_t a, ft_order_t b) { return a == b; }
static inline int ft_order_fits(ft_order_t s, int n) { return n<=(1 << s); }
static inline int ft_order_int(ft_order_t a) { return a; }
#define	ft_order_copy_fft(dst, src)	(dst) = (src)

#define ft_mat_mb_xalloc FTF(mat_mb_xalloc)
static inline void * ft_mat_mb_xalloc(ft_order_t s, ft_mat_mb * x)
{
	return mbdft_alloc(s,(*x));
}

#define ft_mat_bb_xalloc FTF(mat_bb_xalloc)
static inline void * ft_mat_bb_xalloc(ft_order_t s, ft_mat_bb * x)
{
	return bbdft_alloc(s,(*x));
}

#define ft_mat_mb_free FTF(mat_mb_free)
static inline void ft_mat_mb_free(ft_order_t s, ft_mat_mb x)
{
	mbdft_free(s,x);
}

#define ft_mat_bb_free FTF(mat_bb_free)
static inline void ft_mat_bb_free(ft_order_t s, ft_mat_bb x)
{
	bbdft_free(s,x);
}

#define ft_mat_mb_get FTF(mat_mb_get)
static inline ft_t ft_mat_mb_get(ft_order_t s, ft_mat_mb x, int i, int j)
{
	return mbdft_poly(s,x,i,j);
}

#define ft_mat_bb_get FTF(mat_bb_get)
static inline ft_t ft_mat_bb_get(ft_order_t s, ft_mat_bb x, int i, int j)
{
	return bbdft_poly(s,x,i,j);
}


#define	HAS_CONVOLUTION_SPECIAL

#ifdef __cplusplus
}
#endif

#endif	/* OPS_POLY_FFT_H_ */
