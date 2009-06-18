#ifndef FAKE_FFT_H_
#define FAKE_FFT_H_

#include <stdlib.h>
#include <string.h>

#include "macros.h"
#include "alloc_proxy.h"

/* This file is a placeholder for the typical requirements of an FFT
 * interface. Of course, there is nothing interesting being done here.
 * It's just an E-X-A-M-P-L-E. See also the .c file.
 *
 *
 * A unified adapter for the CXX code is in fft-adapter.hpp
 */

/* TODO
 *
 * There is some space wasted here. We're only copying the data when
 * doing the transforms, so we might as well keep the data as is. This
 * implies that it would perhaps make more sense not to have a separate
 * structure for this case -- although of course this could make the
 * calling interface more complicated.
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

struct fake_info_struct {
    size_t n1,n2,n3;
    size_t acc;
    size_t size;
};
typedef struct fake_info_struct fake_info_t[1];
typedef const struct fake_info_struct * fake_info_ptr;
typedef struct fake_info_struct * fake_info_srcptr;

typedef unsigned long fake_t;
typedef fake_t * fake_ptr;
typedef const fake_t * fake_srcptr;

extern void fake_init(fake_info_t p, size_t nF, size_t nG, ...); 
static inline void fake_clear(fake_info_t p MAYBE_UNUSED) {}
static inline fake_ptr fake_alloc(const fake_info_t p, size_t n) {
    return (fake_ptr) mymalloc(n * p->size * sizeof(unsigned long));
}
static inline void fake_free(const fake_info_t p MAYBE_UNUSED,
        fake_ptr x,
        size_t n MAYBE_UNUSED)
{
    myfree(x, n * p->size * sizeof(unsigned long));
}
static inline fake_ptr fake_get(const fake_info_t p, fake_ptr x, size_t k) {
    return x + (k * p->size);
}
static inline void fake_zero(const fake_info_t p, fake_ptr x, size_t n)
{
    memset(x, 0, n * p->size * sizeof(unsigned long));
}
extern void fake_dft(const fake_info_t p, fake_ptr x, unsigned long * F, size_t nF);
extern void fake_compose(const fake_info_t p, fake_ptr y, fake_srcptr x1, fake_srcptr x2);
extern void fake_addcompose(const fake_info_t p, fake_ptr y, fake_srcptr x1, fake_srcptr x2);
extern void fake_add(const fake_info_t p, fake_ptr y, fake_srcptr x1, fake_srcptr x2);
extern void fake_cpy(const fake_info_t p, fake_ptr y, fake_srcptr x);
extern void fake_ift(const fake_info_t p, unsigned long * H, size_t Hl, fake_srcptr h);

extern size_t fake_size(const fake_info_t p);

#ifdef __cplusplus
}
#endif

/* vim: set sw=4 sta et: */
#endif	/* FAKE_FFT_H_ */
