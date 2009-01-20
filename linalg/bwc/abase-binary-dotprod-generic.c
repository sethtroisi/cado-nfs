#include <emmintrin.h>

/* This must be included with the proper macros set */

/* Given the number K of entries represented by the abase type (which is
 * P(nbits)(x), here 1 * 64 = 64), write a K times K matrix as K elements
 * in the area pointed to by w.  */
ABASE_F(void,dotprod,(P(obj_srcptr) x MAYBE_UNUSED,
        P(base_type) * w,
        const P(base_type) * u,
        const P(base_type) * v,
        unsigned int n))
{
    memset(w, 0, P(bytes)(x,P(nbits)(x)));
    for(unsigned int i = 0 ; i < n ; i++) {
        __v2di * w0 = (__v2di*) w;
        for(unsigned int l = 0 ; l < P(repeat)(x) ; l++) {
            // TODO: It's possible to expand more, and use a __v2di
            // mb[4][2], or even [4]. This wouldn't change the code much
            // (see the u128 version), and is likely to speed things up a
            // wee bit maybe.
            __v2di mb[4] = {
                (__v2di) {0, 0},
                (__v2di) {*v, 0},
                (__v2di) {0, *v},
                (__v2di) {*v, *v},
            };
            v++;
            __v2di *sw = w0;
            const P(base_type) * ut = u;
            for(unsigned int k = 0 ; k < P(repeat)(x) ; k++) {
                uint64_t a = *ut++;
                for (unsigned int j = 0; j < 64; j += 2) {
                    *sw ^= mb[a & 3];
                    a >>= 2;
                    sw += P(repeat)(x);
                }
            }
            w0++;
        }
        u += P(repeat)(x);
    }
}
