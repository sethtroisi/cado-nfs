#ifndef FAKE_FFT_H_
#define FAKE_FFT_H_

#include <stdlib.h>
#include <string.h>
/* This file is a placeholder for the typical requirements of an FFT
 * interface. Of course, there is nothing interesting being done here.
 * It's just an E-X-A-M-P-L-E. See also the .c file.
 *
 *
 * A unified adapter for the CXX code is in fft-adapter.hpp
 */


#ifdef __cplusplus
extern "C" {
#endif

struct fake_info_struct {
    int d1,d2,d3;
    unsigned int acc;
    unsigned int size;
};
typedef struct fake_info_struct fake_info_t[1];

typedef ulong * fake_t;
typedef const ulong * fake_src_t;

extern void fake_setup(fake_info_t p, int dF, int dG); 
extern inline fake_t fake_alloc(const fake_info_t p, int n) {
    return (fake_t) malloc(n * p->size * sizeof(ulong));
}
extern inline void fake_free(const fake_info_t p, fake_t x, int n) {
    free(x);
}
extern inline fake_t fake_get(const fake_info_t p, fake_t x, int k) {
        return x + (k * p->size);
}
extern inline void fake_zero(const fake_info_t p, fake_t x, int n)
{
        memset(x, 0, n * p->size * sizeof(unsigned long));
}
extern void fake_dft(const fake_info_t p, fake_t x, unsigned long * F, int dF);
extern void fake_compose(const fake_info_t p, fake_t y, fake_src_t x1, fake_src_t x2);
extern void fake_add(const fake_info_t p, fake_t y, fake_src_t x1, fake_src_t x2);
extern void fake_ift(const fake_info_t p, unsigned long * H, int Hl, fake_src_t h);


#ifdef __cplusplus
}
#endif

/* vim: set sw=4 sta et: */
#endif	/* FAKE_FFT_H_ */
