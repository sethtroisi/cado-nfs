#ifndef CANTOR128_H_
#define CANTOR128_H_

#include <stdint.h>
#include <stdlib.h>
#include "mpfq_2_128.h"

#ifdef __cplusplus
extern "C" {
#endif

struct c128_info_struct {
    int k;
    int n;
};
typedef struct c128_info_struct c128_info_t[1];

typedef mpfq_2_128_elt * c128_t;
/* aliasing is not our friend here ;-( */
typedef mpfq_2_128_elt * c128_src_t;

extern void c128_setup(c128_info_t p, int dF, int dG);
extern c128_t c128_alloc(const c128_info_t p, int n);
extern void c128_free(
        const c128_info_t p MAYBE_UNUSED,
        c128_t x,
        int n MAYBE_UNUSED);
extern c128_t c128_get(const c128_info_t p, c128_t x, int k);
extern void c128_zero(const c128_info_t p, c128_t x, int n);

extern void c128_dft(const c128_info_t p, c128_t x, unsigned long * F, int dF);
extern void c128_compose(const c128_info_t p,
		c128_t y, c128_src_t x1, c128_src_t x2);
extern void c128_add(const c128_info_t p,
		c128_t y, c128_src_t x1, c128_src_t x2);
extern void c128_ift(const c128_info_t p,
		unsigned long * H, int Hl, c128_src_t h);

extern void mulCantor128(unsigned long *H, unsigned long *F, int Fl, unsigned long *G, int Gl);
#ifdef __cplusplus
}
#endif

#endif	/* CANTOR128_H_ */
