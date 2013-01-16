#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

#include "macros.h"
#include "random_generation.h"

mp_limb_t cado_random(cado_random_state_ptr r MAYBE_UNUSED)
{
#ifdef  USE_GMP_RANDOM
    return gmp_urandomb_ui(r->g, GMP_LIMB_BITS);
#else
    return rand ();
#endif
}

void cado_random_area(cado_random_state_ptr r, void * p, size_t s)
{
    for(size_t i = 0 ; i < s ; i++) {
        ((char*)p)[i] = cado_random(r);
    }
#if 0
    /* FIXME This does not work with seeding ! */
    FILE * f = fopen("/dev/urandom", "r");
    BUG_ON(f == NULL);
    size_t r = fread(p, 1, s, f);
    BUG_ON(r != s);
    fclose(f);
#endif
}

void cado_random_init(cado_random_state_ptr r MAYBE_UNUSED, unsigned int s)
{
    if (s == 0)
        s = time(NULL);
#ifdef  USE_GMP_RANDOM
    gmp_randinit_mt(r->g);
    gmp_randseed_ui(r->g, s);
#endif
    srand(s);
}

void cado_random_clear(cado_random_state_ptr r MAYBE_UNUSED)
{
#ifdef  USE_GMP_RANDOM
    gmp_randclear_mt(r->g);
#endif
}
