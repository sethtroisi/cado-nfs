#include "abase-common.h"
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
ABASE_F(void,vdotprod,(
            P(obj_srcptr) x,
            PV(obj_srcptr) y,
        P(base_type) * w,
        const PV(base_type) * u,
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
            const PV(base_type) * ut = u;
            for(unsigned int k = 0 ; k < PV(repeat)(y) ; k++) {
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
        u += PV(repeat)(y);
    }
}

/* XXX This one is copied from the generic version. I'm not completely
 * sure it does the right thing, in particular w.r.t endianness and so
 * on.
 */
/* Transpose a matrix of aboj_t's -- the number of rows of the matrix to
 * be transposed is a runtime variable. Therefore the output consists of
 * abvt's
 */
ABASE_F(void,vtranspose,(PV(obj_srcptr) y MAYBE_UNUSED,
        P(obj_srcptr) x MAYBE_UNUSED,
        PV(base_type) * w,
        const P(base_type) * u))
{
    ASSERT(P(nbits)(x) % 64 == 0);
    ASSERT(PV(nbits)(y) % 64 == 0);
    PV(zero)(y, w, P(nbits)(u));
    unsigned int i,j;
    uint64_t md = 1UL;
    unsigned int od = 0;
    uint64_t ms = 1UL;
    unsigned int os = 0;
    uint64_t * dst = (uint64_t *) w;
    const uint64_t * src = (const uint64_t *) u;
    for(i = 0 ; i < P(nbits)(x) ; i++) {
        for(j = 0 ; j < PV(nbits)(y) ; j++) {
            dst[od] |= md & -((src[os+j*P(nbits)(x)/64] & ms) != 0);
            md <<= 1;
            od += md == 0;
            md += md == 0;
        }
        ms <<= 1;
        os += ms == 0;
        ms += ms == 0;
    }
}

/* XXX This one is copied from the generic version. I'm not completely
 * sure it does the right thing, in particular w.r.t endianness and so
 * on.
 */
ABASE_F(void,vaddmul_tiny,(
            P(obj_srcptr) x MAYBE_UNUSED,
            PV(obj_srcptr) y MAYBE_UNUSED,
            PV(base_type) * w,
            const PV(base_type) * u,
            const P(base_type) * v,
            unsigned int n))
{
    ASSERT(P(nbits)(x) % 64 == 0);
    ASSERT(PV(nbits)(y) % 64 == 0);
    // PV(zero)(y,w,n);
    uint64_t * u0 = (uint64_t *) u;
    uint64_t * u1 = (uint64_t *) (u + P(offset)(x,1));
    uint64_t * w0 = (uint64_t *) w;
    uint64_t * w1 = (uint64_t *) (w + PV(offset)(y,1));
    for (unsigned int j = 0; j < n; j += 2 ) {
        uint64_t * v0 = (uint64_t *) v;
        for(unsigned int l = 0 ; l < PV(repeat)(y) ; l++) {
            union { __v2di s; uint64_t x[2]; } r;
            r.s = (__v2di) {0,0};
            uint64_t * vv = v0;
            for(unsigned int b = 0 ; b < 2 * P(repeat)(x) ; b++) {
                __v2di a = (__v2di) { u0[b], u1[b] };
                __v2di one = (__v2di) { 1, 1, };
                for (unsigned int i = 0; i < 64; i++) {
                    __v2di zw = { *vv, *vv, };
                    r.s ^= (zw & -(a & one));
                    a = _mm_srli_epi64(a, 1);
                    vv += PV(repeat)(y);
                }
            }
            w0[l] ^= r.x[0];
            w1[l] ^= r.x[1];
            v0++;
        }
        w0 += 2 * PV(repeat)(y);
        w1 += 2 * PV(repeat)(y);
        u0 += 4 * P(repeat)(x);
        u1 += 4 * P(repeat)(x);
    }
}

