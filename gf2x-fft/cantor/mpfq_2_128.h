#include <gmp.h>

#if GMP_LIMB_BITS == 32
#include "mpfq_2_128.h.32"
#elif GMP_LIMB_BITS == 64
#include "mpfq_2_128.h.64"
#else
#error "GMP_LIMB_BITS is neither 32 nor 64 ?"
#endif
