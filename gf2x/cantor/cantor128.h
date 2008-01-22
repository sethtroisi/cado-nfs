#ifndef CANTOR128_H_
#define CANTOR128_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include "mpfq_2_128.h"

struct cantor_info_struct {
    int k;
    int n;
};
typedef struct cantor_info_struct cantor_info_t[1];

typedef mpfq_2_128_elt * cantor_transform_t;

extern void cantor_setup_info(cantor_info_t p, int Fl, int Gl);
extern cantor_transform_t cantor_transform_alloc(const cantor_info_t p, int n);
extern void cantor_transform_free(const cantor_info_t p, cantor_transform_t x, int n);
extern inline cantor_transform_t cantor_transform_get(const cantor_info_t p, cantor_transform_t x, int k) {
	return x + k << p->k;
}
extern void cantor_transform(const cantor_info_t p, cantor_transform_t x, uint64_t * F, int Fl);
extern void cantor_compose(const cantor_info_t p, cantor_transform_t y, cantor_transform_t x1, cantor_transform_t x2);
extern void cantor_add(const cantor_info_t p, cantor_transform_t y, cantor_transform_t x1, cantor_transform_t x2);
extern void cantor_itransform(const cantor_info_t p, uint64_t * H, int Hl, cantor_transform_t h);


extern void mulCantor128(uint64_t *H, uint64_t *F, int Fl, uint64_t *G, int Gl);
#ifdef __cplusplus
}
#endif

#endif	/* CANTOR128_H_ */
