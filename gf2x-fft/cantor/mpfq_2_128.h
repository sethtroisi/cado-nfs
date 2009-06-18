#include <gmp.h>

#if GMP_NUMB_BITS == 32
#include "mpfq_2_128.h.32"
#else
#include "mpfq_2_128.h.64"
#endif
