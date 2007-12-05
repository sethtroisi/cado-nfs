#include <gmp.h>

#if GMP_NUMB_BITS == 32
#include "mpfq_2_128.c.32"
#else
#include "mpfq_2_128.c.64"
#endif
