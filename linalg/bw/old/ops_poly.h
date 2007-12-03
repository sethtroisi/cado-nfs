#ifndef OPS_POLY_H_
#define OPS_POLY_H_

#include <stddef.h>
#include <gmp.h>


#define	OPS_POLY_FFT	1
#define	OPS_POLY_SCBK	2
#define	OPS_POLY_IFFT	3
#define	OPS_POLY_WRAP	4

#ifndef	OPS_POLY
#define	OPS_POLY	OPS_POLY_FFT
#endif

/* The ops_poly_XXX.h file has to define the things which might vary with
 * the inner interface. These are the following types and functions
 * (which may be macros).

typedef SOMETHING	ft_t;
typedef SOMETHING	ft_mat_mb;
typedef SOMETHING	ft_mat_bb;

ft_mat_mb ft_mat_mb_alloc(int order);
void ft_mat_mb_free(int order, ft_mat_mb x);
ft_t ft_mat_mb_get(int order, ft_mat_mb x, int i, int j);

ft_mat_bb ft_mat_bb_alloc(int order);
void ft_mat_bb_free(int order, ft_mat_bb x);
ft_t ft_mat_bb_get(int order, ft_mat_bb x, int i, int j);

 */

#if (OPS_POLY == OPS_POLY_FFT)
#include "ops_poly_fft.h"
#elif (OPS_POLY == OPS_POLY_SCBK)
#include "ops_poly_schoolbook.h"
#elif (OPS_POLY == OPS_POLY_IFFT)
#include "ops_poly_integerfft.h"
#elif (OPS_POLY == OPS_POLY_WRAP)
#include "ops_poly_wrap.h"
#endif

#ifndef	ft_mat_mb_alloc
#define	ft_mat_mb_alloc(s,x) ft_mat_mb_xalloc  (s,&(x))
#define	ft_mat_bb_alloc(s,x) ft_mat_bb_xalloc  (s,&(x))
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define ft_zero FTF(zero)  
void ft_zero  (ft_order_t order, ft_t p);

#define ft_one FTF(one)  
void ft_one  (ft_order_t order, ft_t p);

#define ft_convolution FTF(convolution)  
void ft_convolution  (ft_order_t order, ft_t r, ft_t p, ft_t q);

#ifdef	HAS_CONVOLUTION_SPECIAL
#define ft_convolution_special FTF(convolution_special)  
void ft_convolution_special  (ft_order_t order,

		ft_t r, ft_t p, ft_t q,
		unsigned int dg_kill);
#endif
#define ft_itransform FTF(itransform)  
void ft_itransform  (ft_order_t order,

		mp_limb_t * dst, ptrdiff_t stride, int deg,
		ft_t p);
#define ft_transform FTF(transform)  
void ft_transform  (ft_order_t order,

		ft_t p,
		mp_limb_t * src, ptrdiff_t stride, int deg);
#define ft_order FTF(order)  
void ft_order  (ft_order_t *, int ncoeffs);


#ifndef	charptr_list_t_defined
#define	charptr_list_t_defined
struct charptr_list_t {
	const char * p;
	struct charptr_list_t * next;
};
#endif	/* charptr_list_t_defined */

/* nmax : max # of coefficients we have to handle.
 * args : all the command-line arguments that are relevant to the
 * underlying interface.
 */
#define ft_ops_init FTF(ops_init)  
void ft_ops_init  (unsigned int nmax, struct charptr_list_t * args);

#define ft_ops_cleanup FTF(ops_cleanup)  
void ft_ops_cleanup  ();


#ifdef __cplusplus
}
#endif

#endif	/* OPS_POLY_H_ */
