#ifndef gf2x_tfft_H_
#define gf2x_tfft_H_

#ifdef __cplusplus
extern "C" {
#endif

struct gf2x_tfft_info_s {
    size_t bits_a;  // number of bits of operand1
    size_t bits_b;  // number of bits of operand2
    size_t K;       // 0 indicates fallback. negative means split.
    size_t M;
    unsigned long * tmp;
    size_t * perm;
    int split;  // boolean
};

typedef struct gf2x_tfft_info_s gf2x_tfft_info_t[1];
typedef struct gf2x_tfft_info_s * gf2x_tfft_info_ptr;
typedef const struct gf2x_tfft_info_s * gf2x_tfft_info_srcptr;

typedef unsigned long gf2x_tfft_t;
typedef gf2x_tfft_t * gf2x_tfft_ptr;
typedef const gf2x_tfft_t * gf2x_tfft_srcptr;

extern size_t gf2x_tfft_size(gf2x_tfft_info_srcptr o);
extern gf2x_tfft_ptr gf2x_tfft_get(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr ptr, size_t k);
extern void gf2x_tfft_zero(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr ptr, size_t n);
extern void gf2x_tfft_cpy(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr y, gf2x_tfft_srcptr x);
extern gf2x_tfft_ptr gf2x_tfft_alloc(gf2x_tfft_info_srcptr o, size_t n);
extern void gf2x_tfft_free(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr ptr, size_t n);
extern void gf2x_tfft_dft(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tr, const unsigned long * a, size_t bits_a);
extern void gf2x_tfft_compose(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tc, gf2x_tfft_srcptr ta, gf2x_tfft_srcptr tb);
extern void gf2x_tfft_addcompose(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tc, gf2x_tfft_srcptr ta, gf2x_tfft_srcptr tb);
extern void gf2x_tfft_add(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tc, gf2x_tfft_srcptr ta, gf2x_tfft_srcptr tb);
extern void gf2x_tfft_ift(gf2x_tfft_info_srcptr o, unsigned long * c, size_t bits_c, gf2x_tfft_ptr tr);
extern void gf2x_tfft_init(gf2x_tfft_info_ptr o, size_t bits_a, size_t bits_b, ...);
extern void gf2x_tfft_clear(gf2x_tfft_info_ptr o);
extern void gf2x_tfft_init_similar(gf2x_tfft_info_ptr o, size_t bits_a, size_t bits_b, gf2x_tfft_info_srcptr other);
extern int gf2x_tfft_compatible(gf2x_tfft_info_srcptr o1, gf2x_tfft_info_srcptr o2);

#ifdef __cplusplus
}
#endif

#endif	/* gf2x_tfft_H_ */
