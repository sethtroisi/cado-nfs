#include "abase-common.h"
#include "abase-u64.h"

#include "pad.h"
#define P(X)    PAD(abase_u64,X)
#define ABASE_F(t,n,a) t P(n) a

#define need_dotprod_64K_64
#define need_vaddmul_tiny_64K_64L
#define need_vtranspose_64K_64L

#include "abase-binary-dotprod-backends.h"

ABASE_F(void,dotprod,(P(obj_srcptr) x MAYBE_UNUSED,
        P(base_type) * w,
        const P(base_type) * u,
        const P(base_type) * v,
        unsigned int n))
{
    dotprod_64K_64(w,u,v,n,1);
}
ABASE_F(void,vdotprod,(
            P(obj_srcptr) x MAYBE_UNUSED,
            PV(obj_srcptr) y MAYBE_UNUSED,
        P(base_type) * w,
        const PV(base_type) * u,
        const P(base_type) * v,
        unsigned int n))
{
    dotprod_64K_64(w,u,v,n,PV(repeat)(y));
}
ABASE_F(void,vaddmul_tiny,(
            P(obj_srcptr) x MAYBE_UNUSED,
            PV(obj_srcptr) y MAYBE_UNUSED,
            PV(base_type) * w,
            const P(base_type) * u,
            const PV(base_type) * v,
            unsigned int n))
{
    vaddmul_tiny_64K_64L(w,u,v,n,P(repeat)(x),PV(repeat)(y));
}
ABASE_F(void,vtranspose,(PV(obj_srcptr) y MAYBE_UNUSED,
            P(obj_srcptr) x MAYBE_UNUSED,
            PV(base_type) * w,
            const P(base_type) * u))
{
    vtranspose_64K_64L(w,u,PV(repeat)(y),P(repeat)(x));
}
