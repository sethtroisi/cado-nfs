#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <gmp.h>

#include "macros.h"
#include "rootfinder.h"
#include "mpz_poly.h"
#include "modul_poly.h"
#include "gmp_aux.h"

/* put in r[0], ..., r[n-1] the roots of f (of degree d) modulo p,
   and the return value n is the number of roots (without multiplicities) */
int
poly_roots_ulong (unsigned long *r, mpz_t *f, int d, unsigned long p)
{
    int n;

    residueul_t * rr;
    modulusul_t pp;
    modul_initmod_ul(pp, p);
    int i;
    
    if (r == NULL)
      return modul_poly_roots(NULL, f, d, pp);
    
    rr = (residueul_t *) malloc(d * sizeof(residueul_t));
    for(i = 0 ; i < d ; i++) {
      modul_init_noset0(rr[i], pp);
    }
    n = modul_poly_roots(rr, f, d, pp);
    for(int i = 0 ; i < n ; i++) {
      /* The assumption is that p fits within an unsigned long
       * anyway. So the roots do as well.
       */
      r[i] = modul_get_ul(rr[i], pp);
    }
    for(i = 0 ; i < d ; i++) {
      modul_clear(rr[i], pp);
    }
    free(rr);
    modul_clearmod(pp);
    
    return n;
}

int
poly_roots_uint64 (uint64_t * r, mpz_t * f, int d, uint64_t p)
{
    /* This is glue around poly_roots_ulong, nothing more. When uint64
       is larger than ulong, we call mpz_poly_roots_mpz as a fallback */
    unsigned long *rr;
    int i, n;

    if (p > (uint64_t) ULONG_MAX)
      {
        mpz_t pp;

        if (sizeof (long) == 4) printf ("poly_roots_uint64\n");
        mpz_init (pp);
        mpz_set_uint64 (pp, p);
        if (sizeof (long) == 4) gmp_printf ("pp=%Zd\n", pp);
        if (r == NULL)
          {
            if (sizeof (long) == 4) printf ("enter mpz_poly_roots_mpz\n");
            n = mpz_poly_roots_mpz (NULL, f, d, pp);
            if (sizeof (long) == 4) printf ("exit, n=%d\n", n);
          }
        else
          {
            mpz_t *rr;
            rr = malloc ((d + 1) * sizeof (mpz_t));
            for (i = 0; i <= d; i++)
              mpz_init (rr[i]);
            n = mpz_poly_roots_mpz (NULL, f, d, pp);
            for (i = 0; i <= d; i++)
              {
                if (i < n)
                  r[i] = mpz_get_uint64 (rr[i]);
                mpz_clear (rr[i]);
              }
            free (rr);
          }
        mpz_clear (pp);
        return n;
      }

    if (r == NULL)
      return poly_roots_ulong (NULL, f, d, p);

    if (sizeof (unsigned long) != sizeof (uint64_t))
      rr = (unsigned long *) malloc(d * sizeof(unsigned long));
    else
      rr = (unsigned long *) r;
    n = poly_roots_ulong (rr, f, d, p);
    if (sizeof (unsigned long) != sizeof (uint64_t))
      {
        for(i = 0 ; i < n ; i++)
          r[i] = rr[i];
        free (rr);
      }
    return n;
}

int
poly_roots (mpz_t * r, mpz_t * f, int d, mpz_t p)
{
    if (mpz_cmp_ui(p, ULONG_MAX) <= 0) {
        /* There's a chance of using one of our layers. */
        unsigned long pp = mpz_get_ui(p);
        unsigned long * rr;
        int i;
        int n;

        if (r == NULL)
            return poly_roots_ulong(NULL, f, d, pp);

        rr = (unsigned long *) malloc(d * sizeof(unsigned long));
        n = poly_roots_ulong(rr,  f, d, pp);

        for(i = 0 ; i < n ; i++) {
            /* The assumption is that p fits within an unsigned long
             * anyway. So the roots do as well.
             */
            mpz_set_ui(r[i], rr[i]);
        }
        free(rr);
        return n;
    } else {
      int n;
      n = mpz_poly_roots_mpz (r, f, d, p);
      return n;
    }
}
