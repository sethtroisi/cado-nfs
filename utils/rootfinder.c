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


/* Entry point for rootfind routines */


int
mpz_poly_roots (mpz_t * r, mpz_poly_t F, mpz_t p)
{
    int d = F->deg;

    if (mpz_cmp_ui(p, ULONG_MAX) <= 0) {
        /* There's a chance of using one of our layers. */
        unsigned long pp = mpz_get_ui(p);
        unsigned long * rr;
        int i;
        int n;
	
        if (r == NULL)
            return mpz_poly_roots_ulong (NULL, F, pp);

        rr = (unsigned long *) malloc(d * sizeof(unsigned long));
        n = mpz_poly_roots_ulong (rr, F, pp);

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
      n = mpz_poly_roots_mpz (r, F, p);
      return n;
    }
}


/* put in r[0], ..., r[n-1] the roots of F modulo p,
   and the return value n is the number of roots (without multiplicities)     */
int
mpz_poly_roots_ulong (unsigned long *r, mpz_poly_t F, unsigned long p)
{
    int n;

    residueul_t * rr;
    modulusul_t pp;
    modul_initmod_ul(pp, p);
    int i;
    int d = F->deg;
        
    if (r == NULL)
      return modul_poly_roots(NULL, F, pp);
    
    rr = (residueul_t *) malloc(d * sizeof(residueul_t));
    for(i = 0 ; i < d ; i++) {
      modul_init_noset0(rr[i], pp);
    }
    n = modul_poly_roots(rr, F, pp);
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
mpz_poly_roots_uint64 (uint64_t * r, mpz_poly_t F, uint64_t p)
{
    /* This is glue around poly_roots_ulong, nothing more. When uint64
       is larger than ulong, we call mpz_poly_roots_mpz as a fallback */
    unsigned long *rr;
    int i, n;
    int d = F->deg;

    if (p > (uint64_t) ULONG_MAX)
      {
        mpz_t pp;
        mpz_init (pp);
        mpz_set_uint64 (pp, p);
        if (r == NULL)
          n = mpz_poly_roots_mpz (NULL, F, pp);
        else
          {
            mpz_t *rr;
            rr = malloc ((d + 1) * sizeof (mpz_t));
            for (i = 0; i <= d; i++)
              mpz_init (rr[i]);
            n = mpz_poly_roots_mpz (NULL, F, pp);
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
      return mpz_poly_roots_ulong (NULL, F, p);

    if (sizeof (unsigned long) != sizeof (uint64_t))
      rr = (unsigned long *) malloc(d * sizeof(unsigned long));
    else
      rr = (unsigned long *) r;
    n = mpz_poly_roots_ulong (rr, F, p);
    if (sizeof (unsigned long) != sizeof (uint64_t))
      {
        for(i = 0 ; i < n ; i++)
          r[i] = rr[i];
        free (rr);
      }
    return n;
}




typedef int (*sortfunc_t) (const void *, const void *);

static int mpz_poly_coeff_cmp(const mpz_t *a, const mpz_t *b) {
  return mpz_cmp(*a, *b) < 0 ? -1 : 1;
}

/* Assuming f is a (squarefree) product of linear factors mod p, splits it
   and put the corresponding roots mod p in r[]. Return number of roots
   which should be degree of f. Assumes p is odd, and deg(f) >= 1. */
static int
mpz_poly_cantor_zassenhaus (mpz_t *r, mpz_poly_t f, const mpz_t p, int depth)
{
  mpz_t a, aux;
  mpz_poly_t q, h, ff;
  int d = f->deg, dq, n, m;

  mpz_init (a);
  mpz_init (aux);

  /* linear polynomial */
  if (d == 1) {
    mpz_neg (aux, f->coeff[1]);
    mpz_invert (a, aux, p);
    mpz_mul (r[0], a, f->coeff[0]);
    mpz_fdiv_r (r[0], r[0], p);
    n = 1;
    goto clear_a;
  }

  /* if f has degree d, then q,h may have up to degree 2d-1 in the
     powering algorithm */
  mpz_poly_init (q, 2 * d - 1);
  mpz_poly_init (h, 2 * d - 1);
  mpz_poly_init (ff, d);

  /* random polynomial by a */
  mpz_set_ui (a, lrand48());
  mpz_fdiv_r (a, a, p);
  for (;;)
  {
    /* q=x+a */
    mpz_set_ui (aux, 1);
    mpz_poly_setcoeff (q, 1, aux);
    mpz_poly_setcoeff (q, 0, a);

    /* h=(x+a)^((p-1)/2) mod (f, p) */
    mpz_sub_ui (aux, p, 1);
    mpz_divexact_ui (aux, aux, 2);
    mpz_poly_power_mod_f_mod_mpz (h, q, f, aux, p);
    mpz_poly_sub_ui (h, 1);

    /* q = gcd(f,h) */
    mpz_poly_copy (q, f);
    mpz_poly_gcd_mpz (q, h, p);
    dq = q->deg;
    ASSERT (dq >= 0);

    /* recursion-split */
    if (0 < dq && dq < d)
    {
      n = mpz_poly_cantor_zassenhaus (r, q, p, depth+1);
      ASSERT (n == dq);

      mpz_poly_copy (ff, f);
      /* mpz_poly_divexact clobbers its 2nd arg */
      mpz_poly_divexact (h, ff, q, p);
      m = mpz_poly_cantor_zassenhaus (r + n, h, p, depth + 1);
      ASSERT (m == h->deg);
      n += m;
      break;
    }

    mpz_add_ui (a, a, 1); /* no need to reduce mod p, since it will be done
                             in mpz_poly_power_mod_f_mod_mpz */
  }

  mpz_poly_clear (q);
  mpz_poly_clear (h);
  mpz_poly_clear (ff);

clear_a:
  mpz_clear (a);
  mpz_clear (aux);
  return n;
}


/* Solve f(x)=0 (mod p), where p is a prime. Return the number of roots.
   Assume d (the degree of f) is at least 1.
 */
int
mpz_poly_roots_mpz (mpz_t *r, mpz_poly_t f, const mpz_t p)
{
  int nr = 0;
  mpz_t tmp;
  mpz_poly_t fp, g, h;
  int d = f->deg;

  ASSERT(d >= 1);

  mpz_init (tmp);

  mpz_poly_init (fp, d);
  mpz_poly_init (g, 2*d-1);
  mpz_poly_init (h, 2*d-1);

  /* reduce f to monic and modulo p */
  /* mpz_poly_set (mpz_poly_f, f, d); */
  mpz_poly_reduce_makemonic_mod_mpz (fp, f, p);
  if (fp->deg <= 0)
    goto clear_and_exit;
  /* h=x^p-x (mod mpz_poly_fp) */
  mpz_set_ui (tmp, 1UL);
  mpz_poly_setcoeff (g, 1, tmp);
  mpz_poly_power_mod_f_mod_mpz (h, g, fp, p, p);
  mpz_poly_sub (h, h, g);
  /* g = gcd (mpz_poly_fp, h) */
  mpz_poly_gcd_mpz (fp, h, p);
  /* mpz_poly_fp contains gcd(x^p-x, f) */
  nr = fp->deg;
  ASSERT (nr >= 0);

  /* If r is NULL, we only return the number of roots. */
  if (r != NULL && nr > 0)
  {
    int n MAYBE_UNUSED = mpz_poly_cantor_zassenhaus (r, fp, p, 0);
    ASSERT (n == nr);
  }

 clear_and_exit:
  mpz_poly_clear(fp);
  mpz_poly_clear(g);
  mpz_poly_clear(h);
  mpz_clear (tmp);

  /* Sort the roots */
  if (r && nr)
    qsort(r, nr, sizeof(mpz_t), (sortfunc_t) &mpz_poly_coeff_cmp);

  return nr;
}

