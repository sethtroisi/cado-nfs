/* This header is intentionally NOT #ifdef-protected, because it may be
 * included several times.
 */

/* This is most of the meat for binary abase implementations.
 * Must be included with ABASE_F and P properly defined.
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef AVOID_GENERIC_obj_set_nbys
/* This generic version assumes compile-time-constant stride */
ABASE_F(void,obj_set_nbys,(P(obj_ptr) x MAYBE_UNUSED, unsigned int nbys))
{
    ASSERT_ALWAYS(nbys == P(nbits)(x));
}
#endif

#ifndef AVOID_GENERIC_offset
ABASE_F(ptrdiff_t, offset, (P(obj_srcptr) x MAYBE_UNUSED, unsigned int k))
{
    return k * P(repeat)(x);
}
#endif

#ifndef AVOID_GENERIC_bytes
ABASE_F(size_t, bytes, (P(obj_srcptr) x MAYBE_UNUSED, unsigned int k))
{
    return P(offset)(x, k) * sizeof(P(base_type));
}
#endif

#ifndef AVOID_GENERIC_init
ABASE_F(P(base_type) *,init,(P(obj_srcptr) x MAYBE_UNUSED, unsigned int n))
{
    return (P(base_type) *) malloc(P(bytes)(x,n));
    // return (P(base_type) *) electric_alloc(n * P(repeat)(x) * sizeof(P(base_type)));
}
#endif

#ifndef AVOID_GENERIC_clear
ABASE_F(P(base_type) *, clear, (P(obj_srcptr) x MAYBE_UNUSED,
        P(base_type) * p,
        unsigned int n MAYBE_UNUSED))
{
    free(p);
    // electric_free(p, P(bytes)(x,n));
    return NULL;
}
#endif

#ifndef AVOID_GENERIC_zero
ABASE_F(void, zero, (P(obj_srcptr) x MAYBE_UNUSED,
        P(base_type) * p,
        unsigned int n MAYBE_UNUSED))
{
    memset(p, 0, P(bytes)(x,n));
}
#endif

#ifndef AVOID_GENERIC_random
ABASE_F(void,random, (P(obj_srcptr) x MAYBE_UNUSED,
        P(base_type) * p,
        unsigned int n MAYBE_UNUSED))
{
    myrand_area(p, P(bytes)(x,n));
}
#endif

#ifndef AVOID_GENERIC_set_ui
ABASE_F(void,set_ui, (P(obj_srcptr) x MAYBE_UNUSED,
        P(base_type) * p,
        unsigned int k, unsigned int v))
{
    ASSERT(k < P(nbits)(x));
    uint64_t * xp = (uint64_t *) p;
    uint64_t mask = ((uint64_t)1) << (k%64);
    xp[k/64] = (xp[k/64] & ~mask) | ((v << (k%64))&mask);
}
#endif

#ifndef AVOID_GENERIC_copy
ABASE_F(void, copy, (P(obj_srcptr) x MAYBE_UNUSED,
        P(base_type) * q,
        const P(base_type) * p,
        unsigned int n MAYBE_UNUSED))
{
    memcpy(q, p, P(bytes)(x,n));
}
#endif

#ifndef AVOID_GENERIC_is_zero
ABASE_F(int, is_zero, (P(obj_srcptr) x MAYBE_UNUSED,
        const P(base_type) * p,
        unsigned int n MAYBE_UNUSED))
{
    for(unsigned int i = 0 ; i < n * P(repeat)(x) ; i++) {
        if (p[i])
            return 0;
    }
    return 1;
}
#endif

#ifndef AVOID_GENERIC_add
ABASE_F(void,add,
        (P(obj_srcptr) x MAYBE_UNUSED,
        P(base_type) * dst, const P(base_type) * src))
{
    for(unsigned int i = 0 ; i < P(repeat)(x) ; i++) {
        dst[i] ^= src[i];
    }
}
#endif

/* The initf versions do alloca -- it's meant for temporaries
 */
#ifndef AVOID_GENERIC_initf
ABASE_F(P(base_type) *, initf, (P(obj_srcptr) x MAYBE_UNUSED))
{
    return (P(base_type) *) alloca(P(bytes)(x,1));
}
#endif

#ifndef AVOID_GENERIC_clearf
ABASE_F(P(base_type) *, clearf, (P(obj_srcptr) x MAYBE_UNUSED,
        P(base_type) * p MAYBE_UNUSED))
{
    return NULL;
}
#endif

#ifndef AVOID_GENERIC_read
ABASE_F(unsigned int, read, (P(obj_srcptr) x MAYBE_UNUSED,
        FILE * f, P(base_type) * p, unsigned int n))
{
    return fread(p, P(bytes)(x,1), n, f);
}
#endif

#ifndef AVOID_GENERIC_write
ABASE_F(unsigned int, write, (P(obj_srcptr) x MAYBE_UNUSED,
        FILE * f, P(base_type) * p, unsigned int n))
{
    return fread(p, P(bytes)(x,1), n, f);
}
#endif

#ifdef __cplusplus
}
#endif
