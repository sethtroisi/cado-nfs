#include <stdint.h>

/*
 * This file contains classical typedefs for SSE, and some inline macros
 * for instructions that are not given by gcc builtins.
 *
 * We put only what is needed by Cado. 
 */

typedef double v2df __attribute__ ((vector_size (16)));
typedef uint64_t v2di __attribute__ ((vector_size (16)));


/*
 * In gcc-4.3.0, the psrlq128 has been (vaguely) disabled. We have to 
 * use our own macro instead of __builtin_ia32_psrlq128().
 * Hopefully this gives the same efficiency.
 */

static inline v2di cado_psrlq128(v2di x, v2di sh)
{
    // In theory, the "x,m" constraint is appropriate, but for some
    // reason, gcc-4.3 does not like it (internal compiler error!).
    // Let's stick to "x,x" constraint, which might be slightly worse
    // in cases where the shift value changes frequently (this is not the
    // case in our application).
    // __asm__("psrlq %2,%0" : "+x,x"(x) : "0,0"(x), "x,m"(sh));

    __asm__("psrlq %1,%0" : "+x,x"(x) : "x,m"(sh));
    
    //__asm__("psrlq %2,%0" : "=x"(x) : "0"(x), "x"(sh));
    return x;
}
