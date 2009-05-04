#ifndef FFT_ON_MATRICES_H_
#define FFT_ON_MATRICES_H_

#include "ops_poly.hpp"

#ifdef	__cplusplus
extern "C" {
#endif

struct dft_mb {
	ft_order_t::mat_mb p;
	unsigned int degree;
	ft_order_t order;
};

struct dft_bb {
	ft_order_t::mat_bb p;
	unsigned int degree;
	ft_order_t order;
};

extern struct dft_mb *
fft_ec_dft(     struct e_coeff *, ft_order_t const&,    double *);

extern struct dft_bb *
fft_tp_dft(     struct t_poly *,  ft_order_t const&,    double *);

#ifdef	HAS_CONVOLUTION_SPECIAL
extern struct dft_mb *
fft_mbb_conv_sp(struct dft_mb *,  struct dft_bb *, unsigned int, double *);
#endif

extern struct dft_mb *
fft_mbb_conv(struct dft_mb *,  struct dft_bb *, double *);

extern struct dft_bb *
fft_bbb_conv(   struct dft_bb *,  struct dft_bb *, double *);

extern void 
fft_mb_invdft(  bw_mbpoly,        struct dft_mb *, unsigned int, double *);

extern void 
fft_tp_invdft(  struct t_poly *,  struct dft_bb *, double *);

extern struct dft_bb *
fft_bb_dft_init_one(unsigned int);

extern void
dft_mb_free(struct dft_mb *);

extern void
dft_bb_free(struct dft_bb *);

#ifdef	__cplusplus
}
#endif

#endif	/* FFT_ON_MATRICES_H_ */
