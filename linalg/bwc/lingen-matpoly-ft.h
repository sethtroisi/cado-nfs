#ifndef LINGEN_MATPOLY_FT_H_
#define LINGEN_MATPOLY_FT_H_

#ifndef HAVE_MPIR
#error "This interface is only available with MPIR"
#endif

#include "lingen-matpoly.h"
#include "flint-fft/fft.h"

struct matpoly_ft_s {
    unsigned int m;
    unsigned int n;
    void * data;
};

typedef struct matpoly_ft_s matpoly_ft[1];
typedef struct matpoly_ft_s * matpoly_ft_ptr;
typedef const struct matpoly_ft_s * matpoly_ft_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

void matpoly_ft_init(abdst_field ab, matpoly_ft_ptr t, unsigned int m, unsigned int n, struct fft_transform_info * fti);
void matpoly_ft_zero(abdst_field ab, matpoly_ft_ptr t, struct fft_transform_info * fti);
void matpoly_ft_export(abdst_field ab, matpoly_ft_ptr t, struct fft_transform_info * fti);
void matpoly_ft_import(abdst_field ab, matpoly_ft_ptr t, struct fft_transform_info * fti);
void matpoly_ft_clear(abdst_field ab, matpoly_ft_ptr t, struct fft_transform_info * fti);
void matpoly_ft_dft(abdst_field ab, matpoly_ft_ptr t, matpoly_ptr p, struct fft_transform_info * fti);
void matpoly_ft_add(abdst_field ab, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti);
/*
void matpoly_ft_sub(abdst_field ab, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti);
*/
void matpoly_ft_mul(abdst_field ab, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti);
void matpoly_ft_addmul(abdst_field ab, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti);
void matpoly_ft_ift(abdst_field ab, matpoly_ptr p, matpoly_ft_ptr t, struct fft_transform_info * fti);
void matpoly_ft_ift_mp(abdst_field ab, matpoly_ptr p, matpoly_ft_ptr t, unsigned int shift, struct fft_transform_info * fti);

#ifdef __cplusplus
}
#endif


#endif	/* LINGEN_MATPOLY_FT_H_ */
