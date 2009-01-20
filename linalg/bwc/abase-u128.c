#include "abase-u128.h"

#include "pad.h"
#define P(X)    PAD(abase_u128,X)
#define ABASE_F(t,n,a) t P(n) a

ABASE_F(void,dotprod,(P(obj_srcptr) x,
        P(base_type) * w,
        const P(base_type) * u,
        const P(base_type) * v,
        unsigned int n))
{
    memset(w, 0, P(bytes)(x,P(nbits)(x)));
    for(unsigned int i = 0 ; i < n ; i++) {
        __v2di * w0 = (__v2di*) w;
        for(unsigned int l = 0 ; l < P(repeat)(x) ; l++) {
#define SSE_0   ((__v2di) {0,0} )
            __v2di mb[4][2] = { 
                {SSE_0, SSE_0},
                {*v, SSE_0},
                {SSE_0, *v},
                {*v, *v},
            };
#undef SSE_0
            v++;
            __v2di *sw = w0;
            const uint64_t * ut = (const uint64_t *) u;
            for(unsigned int k = 0 ; k < 2 * P(repeat)(x) ; k++) {
                uint64_t a = *ut++;
                for (unsigned int j = 0; j < 64; j += 2) {
                    *sw ^= mb[a & 3][0];
                    sw += P(repeat)(x);
                    *sw ^= mb[a & 3][1];
                    sw += P(repeat)(x);
                    a >>= 2;
                }
            }
            w0++;
        }
        u += P(repeat)(x);
    }
}
