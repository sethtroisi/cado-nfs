#ifndef FFT_ON_MATRICES_H_
#define FFT_ON_MATRICES_H_

#ifdef	__cplusplus
extern "C" {
#endif

struct dft_mb {
	bw_mbdft p;
	unsigned int degree;
	unsigned int order;
};

struct dft_bb {
	bw_bbdft p;
	unsigned int degree;
	unsigned int order;
};

extern unsigned int 
ceil_log2(unsigned int);

extern struct dft_mb *
fft_ec_dft(     struct e_coeff *, unsigned int,    double *);

extern struct dft_bb *
fft_tp_dft(     struct t_poly *,  unsigned int,    double *);

extern struct dft_mb *
fft_mbb_conv_sp(struct dft_mb *,  struct dft_bb *, unsigned int, double *);

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

extern void
prepare_fft_engine(unsigned int, int);

#ifdef	__cplusplus
}
#endif

#endif	/* FFT_ON_MATRICES_H_ */
