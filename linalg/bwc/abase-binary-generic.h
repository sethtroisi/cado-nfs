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

/* The is_zero thing is sort of lacking from the libc interface. It makes
 * sense using only the striding value, but may be slower than this loop
 * in that byte-oriented scanning could be slow. */
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

#ifdef __cplusplus
}
#endif
