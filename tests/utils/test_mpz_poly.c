#include "cado.h"
#include "macros.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpz_poly.h"
#include "tests_common.h"

/* put random coefficients of k bits in a polynomial (already initialized) */
static void
mpz_poly_random (mpz_poly_t f, int d, int k)
{
  int i;
  mpz_t u;

  ASSERT_ALWAYS (k > 0);
  mpz_poly_realloc (f, d + 1);
  mpz_init_set_ui (u, 1);
  mpz_mul_2exp (u, u, k - 1); /* u = 2^(k-1) */
  for (i = 0; i <= d; i++)
    {
      mpz_rrandomb (f->coeff[i], state, k); /* 0 to 2^k-1 */
      mpz_sub (f->coeff[i], f->coeff[i], u); /* -2^(k-1) to 2^(k-1)-1 */
    }
  mpz_clear (u);
  while (d >= 0 && mpz_cmp_ui (f->coeff[d], 0) == 0)
    d--;
  f->deg = d;
}

static void
mpz_poly_mul_xk (mpz_poly_t f, int k)
{
  int i;

  ASSERT_ALWAYS (k >= 0);

  if (k == 0)
    return;

  mpz_poly_realloc (f, f->deg + k + 1);
  for (i = f->deg; i >= 0; i--)
    mpz_set (f->coeff[i + k], f->coeff[i]);
  for (i = 0; i < k; i++)
    mpz_set_ui (f->coeff[i], 0);
  f->deg += k;
}

static void mpz_poly_setcoeffs_ui_var(mpz_poly_t f, int d, ...)
{
    va_list ap;
    va_start(ap, d);
    mpz_poly_realloc(f, d + 1);
    for(int i = 0 ; i <= d ; i++) {
        mpz_set_ui(f->coeff[i], va_arg(ap, int));
    }
    mpz_poly_cleandeg(f, d);
    va_end(ap);
}

void
test_polymodF_mul ()
{
  int d1, d2, d;
  mpz_poly_t F, T, U;
  polymodF_t P1, P2, Q, P1_saved;
  int k = 2 + (lrand48 () % 128), count = 0;
  mpz_t c;

  mpz_poly_init (T, -1);
  mpz_poly_init (U, -1);
  mpz_init (c);
  for (d = 1; d <= 10; d++)
    {
      mpz_poly_init (F, d);
      mpz_poly_init (Q->p, d-1);
      do mpz_poly_random (F, d, k); while (F->deg == -1);
      for (d1 = 1; d1 <= 10; d1++)
        {
          mpz_poly_init (P1->p, d1);
          mpz_poly_init (P1_saved->p, d1);
          mpz_poly_random (P1->p, d1, k);
          mpz_poly_set (P1_saved->p, P1->p);
          P1->v = 0;
          for (d2 = 1; d2 <= 10; d2++)
            {
              mpz_poly_init (P2->p, d2);
              mpz_poly_random (P2->p, d2, k);
              P2->v = 0;
              if ((++count % 3) == 0)
                polymodF_mul (Q, P1, P2, F);
              else if ((count % 3) == 1)
                {
                  polymodF_mul (P1, P1, P2, F);
                  mpz_poly_set (Q->p, P1->p);
                  Q->v = P1->v;
                  mpz_poly_set (P1->p, P1_saved->p);
                  P1->v = 0;
                }
              else
                {
                  polymodF_mul (P1, P2, P1, F);
                  mpz_poly_set (Q->p, P1->p);
                  Q->v = P1->v;
                  mpz_poly_set (P1->p, P1_saved->p);
                  P1->v = 0;
                }
              /* check that Q->p = lc(F)^Q->v * P1 * P1 mod F */
              ASSERT_ALWAYS (Q->p->deg < F->deg);
              mpz_poly_mul (T, P1->p, P2->p);
              mpz_pow_ui (c, F->coeff[F->deg], Q->v);
              mpz_poly_mul_mpz (T, T, c);
              mpz_poly_sub (T, T, Q->p);
              /* T should be a multiple of F */
              while (T->deg >= F->deg)
                {
                  int oldd = T->deg;
                  if (!mpz_divisible_p (T->coeff[T->deg], F->coeff[F->deg]))
                    {
                      printf ("Error in test_polymodF_mul\n");
                      printf ("F="); mpz_poly_fprintf (stdout, F);
                      printf ("P1="); mpz_poly_fprintf (stdout, P1->p);
                      printf ("P2="); mpz_poly_fprintf (stdout, P2->p);
                      printf ("Q="); mpz_poly_fprintf (stdout, Q->p);
                      exit (1);
                    }
                  mpz_divexact (c, T->coeff[T->deg], F->coeff[F->deg]);
                  mpz_poly_mul_mpz (U, F, c);
                  /* multiply U by x^(T->deg - F->deg) */
                  mpz_poly_mul_xk (U, T->deg - F->deg);
                  mpz_poly_sub (T, T, U);
                  ASSERT_ALWAYS (T->deg < oldd);
                }
              if (T->deg != -1)
                {
                  printf ("count=%d\n", count);
                }
              ASSERT_ALWAYS (T->deg == -1);
              mpz_poly_clear (P2->p);
            }
          mpz_poly_clear (P1->p);
          mpz_poly_clear (P1_saved->p);
        }
      mpz_poly_clear (F);
      mpz_poly_clear (Q->p);
    }
  mpz_poly_clear (T);
  mpz_poly_clear (U);
  mpz_clear (c);
}

/* void */
/* test_mpz_poly_roots_mpz (unsigned long iter) */
/* { */
/*   mpz_t r[10], f[10], p, res; */
/*   unsigned long i, n, d; */
/*   mpz_poly_t F; */

/*   for (i = 0; i < 10; i++) */
/*     { */
/*       mpz_init (r[i]); */
/*       mpz_init (f[i]); */
/*     } */
/*   mpz_init (p); */
/*   mpz_init (res); */

/*   /\* -16*x^2 - x - 2 mod 17 *\/ */
/*   mpz_set_si (f[2], -16); */
/*   mpz_set_si (f[1], -1); */
/*   mpz_set_si (f[0], -2); */
/*   F->coeff = f; */
/*   F->deg = 2; */
/*   mpz_set_ui (p, 17); */
/*   n = mpz_poly_roots_mpz (r, F, p); */
/*   ASSERT_ALWAYS(n == 2); */
/*   ASSERT_ALWAYS(mpz_cmp_ui (r[0], 2) == 0); */
/*   ASSERT_ALWAYS(mpz_cmp_ui (r[1], 16) == 0); */

/*   /\* 9*x^2 + 6*x + 3 mod 3 *\/ */
/*   mpz_set_si (f[2], 9); */
/*   mpz_set_si (f[1], 6); */
/*   mpz_set_si (f[0], 3); */
/*   mpz_set_ui (p, 3); */
/*   n = mpz_poly_roots_mpz (r, F, p); */
/*   ASSERT_ALWAYS(n == 0); */

/*   /\* try random polynomials *\/ */
/*   for (i = 0; i < iter; i++) */
/*     { */
/*       d = 1 + (lrand48 () % 7); */
/*       for (n = 0; n <= d; n++) */
/*         mpz_set_si (f[n], mrand48 ()); */
/*       mpz_urandomb (p, state, 128); */
/*       mpz_nextprime (p, p); */
/*       while (mpz_divisible_p (f[d], p)) */
/*         mpz_set_si (f[d], mrand48 ()); */
/*       F->coeff = f; */
/*       F->deg = d; */
/*       n = mpz_poly_roots_mpz (r, F, p); */
/*       ASSERT_ALWAYS (n <= d); */
/*       while (n-- > 0) */
/*         { */
/*           mpz_poly_eval (res, F, r[n]); */
/*           ASSERT_ALWAYS (mpz_divisible_p (res, p)); */
/*         } */
/*     } */

/*   for (i = 0; i < 10; i++) */
/*     { */
/*       mpz_clear (r[i]); */
/*       mpz_clear (f[i]); */
/*     } */
/*   mpz_clear (p); */
/*   mpz_clear (res); */
/* } */

/* also exercises mpz_poly_mul */
void
test_mpz_poly_sqr_mod_f_mod_mpz (unsigned long iter)
{
  while (iter--)
    {
      mpz_poly_t Q, P, f;
      mpz_t m, invm;
      int k = 2 + (lrand48 () % 127);
      int d = 1 + (lrand48 () % 7);

      mpz_init (m);
      do mpz_urandomb (m, state, k); while (mpz_tstbit (m, 0) == 0);
      mpz_poly_init (f, d);
      mpz_init (invm);
      while (1)
        {
          mpz_poly_random (f, d, k);
          if (f->deg < d)
            continue;
          mpz_gcd (invm, m, f->coeff[d]);
          if (mpz_cmp_ui (invm, 1) == 0)
            break;
        }
      barrett_init (invm, m);
      mpz_poly_init (P, d - 1);
      if (iter)
        mpz_poly_random (P, d - 1, k);
      else
        P->deg = -1; /* P=0 */
      mpz_poly_init (Q, d - 1);
      mpz_poly_sqr_mod_f_mod_mpz (Q, P, f, m, invm);
      if (iter == 0)
        ASSERT_ALWAYS(Q->deg == -1);
      mpz_poly_mul_mod_f_mod_mpz (Q, P, P, f, m, invm);
      if (iter == 0)
        ASSERT_ALWAYS(Q->deg == -1);
      mpz_poly_clear (f);
      mpz_poly_clear (P);
      mpz_poly_clear (Q);
      mpz_clear (m);
      mpz_clear (invm);
    }
}

/* Also exercises mpz_poly_getcoeff, mpz_poly_setcoeff_int64,
   mpz_poly_setcoeff_si, mpz_poly_cmp, mpz_poly_eval,
   mpz_poly_eval_mod_mpz_barrett and mpz_poly_eval_several_mod_mpz_barrett */
void
test_mpz_poly_fprintf (void)
{
  mpz_poly_t f, g;
  mpz_t c, v[2], m, invm;
  int res;
  mpz_poly_srcptr F[2];
  mpz_ptr V[2];

  mpz_poly_init (f, 1);
  mpz_poly_init (g, 1);
  F[0] = f;
  F[1] = g;
  V[0] = (mpz_ptr) v[0];
  V[1] = (mpz_ptr) v[1];
  mpz_init (c);
  mpz_init (v[0]);
  mpz_init (v[1]);
  mpz_init_set_ui (m, 11);
  mpz_init (invm);

  barrett_init (invm, m);

  f->deg = -1;
  mpz_poly_fprintf (stdout, f);
  mpz_poly_getcoeff (c, 0, f);
  ASSERT_ALWAYS (mpz_cmp_ui (c, 0) == 0);
  mpz_set_ui (c, 17);
  mpz_poly_eval (v[0], f, c);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 0) == 0);
  mpz_poly_eval_mod_mpz_barrett (v[0], f, c, m, invm);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 0) == 0);

  f->deg = 0;
  mpz_set_ui (f->coeff[0], 17); /* f = 17 */
  mpz_poly_fprintf (stdout, f);
  mpz_set_ui (c, 42);
  mpz_poly_eval (v[0], f, c);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 17) == 0);
  mpz_poly_eval_mod_mpz_barrett (v[0], f, c, m, invm);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 6) == 0);

  mpz_poly_setcoeff_int64 (f, 1, 42); /* f = 42*x+17 */
  mpz_poly_fprintf (stdout, f);
  mpz_set_ui (c, 1);
  mpz_poly_eval (v[0], f, c);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 59) == 0);
  mpz_set_si (c, -1);
  mpz_poly_eval (v[0], f, c);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], -25) == 0);
  mpz_poly_eval_mod_mpz_barrett (v[0], f, c, m, invm);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 8) == 0);

  mpz_poly_setcoeff_si (f, 2, -3); /* f = -3*x^2+42*x+17 */
  mpz_poly_fprintf (stdout, f);

  mpz_poly_set (g, f);
  res = mpz_poly_cmp (f, g);
  ASSERT_ALWAYS (res == 0);
  mpz_add_ui (g->coeff[g->deg], g->coeff[g->deg], 1); /* g = -2*x^2+42*x+17 */
  res = mpz_poly_cmp (f, g);
  ASSERT_ALWAYS (res != 0);
  mpz_set_si (c, 3);
  mpz_poly_eval_several_mod_mpz_barrett (V, F, 1, c, m, invm);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], 6) == 0);
  mpz_poly_eval_several_mod_mpz_barrett (V, F, 2, c, m, invm);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], 6) == 0);
  ASSERT_ALWAYS (mpz_cmp_si (v[1], 4) == 0);
  mpz_poly_setcoeff_si (g, g->deg + 1, 1); /* g = x^3-2*x^2+42*x+17 */
  res = mpz_poly_cmp (f, g);
  ASSERT_ALWAYS (res != 0);
  mpz_set_si (c, -3);
  mpz_poly_eval_several_mod_mpz_barrett (V, F, 1, c, m, invm);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], 7) == 0);
  mpz_poly_eval_several_mod_mpz_barrett (V, F, 2, c, m, invm);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], 7) == 0);
  ASSERT_ALWAYS (mpz_cmp_si (v[1], 0) == 0);
  /* test with one zero polynomial */
  g->deg = -1;
  mpz_set_si (c, 3);
  mpz_poly_eval_several_mod_mpz_barrett (V, F, 2, c, m, invm);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], 6) == 0);
  ASSERT_ALWAYS (mpz_cmp_si (v[1], 0) == 0);

  mpz_poly_clear (f);
  mpz_poly_clear (g);
  mpz_clear (c);
  mpz_clear (v[0]);
  mpz_clear (v[1]);
  mpz_clear (m);
  mpz_clear (invm);
}

void
test_mpz_poly_div_2_mod_mpz (void)
{
  mpz_poly_t f;
  mpz_t m;

  mpz_init_set_ui (m, 17);
  mpz_poly_init (f, -1);
  mpz_poly_setcoeff_si (f, 0, 1);
  mpz_poly_setcoeff_si (f, 1, -2);
  mpz_poly_setcoeff_si (f, 2, -3);
  mpz_poly_setcoeff_si (f, 3, 4);
  mpz_poly_div_2_mod_mpz (f, f, m);
  ASSERT_ALWAYS(mpz_cmp_si (f->coeff[0], 9) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (f->coeff[1], -1) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (f->coeff[2], 7) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (f->coeff[3], 2) == 0);
  mpz_poly_clear (f);
  mpz_clear (m);
}

void
test_mpz_poly_derivative (void)
{
  mpz_poly_t f, df;

  mpz_poly_init (f, -1);
  mpz_poly_init (df, 1);

  mpz_poly_derivative (df, f);
  ASSERT_ALWAYS(df->deg == -1);

  mpz_poly_setcoeff_si (f, 0, 17); /* f = 17 */
  mpz_poly_derivative (df, f);
  ASSERT_ALWAYS(df->deg == -1);

  mpz_poly_setcoeff_si (f, 1, 42); /* f = 42*x + 17 */
  mpz_poly_derivative (df, f);
  ASSERT_ALWAYS(df->deg == 0);
  ASSERT_ALWAYS(mpz_cmp_si (df->coeff[0], 42) == 0);

  mpz_poly_setcoeff_si (f, 2, -3); /* f = -3*x^2 + 42*x + 17 */
  mpz_poly_derivative (df, f);
  ASSERT_ALWAYS(df->deg == 1);
  ASSERT_ALWAYS(mpz_cmp_si (df->coeff[0], 42) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (df->coeff[1], -6) == 0);

  mpz_poly_clear (f);
  mpz_poly_clear (df);
}

/* also exercises mpz_poly_power_mod_f_mod_mpz and
   mpz_poly_power_mod_f_mod_mpz_Barrett */
void
test_mpz_poly_power_mod_f_mod_ui (void)
{
  mpz_poly_t Q, P, f;
  mpz_t a, pp, invp;
  unsigned long p = 4294967291UL;

  mpz_poly_init (Q, -1);
  mpz_poly_init (P, -1);
  mpz_poly_init (f, -1);
  mpz_init (a);
  mpz_init_set_ui (pp, p);
  mpz_init (invp);
  barrett_init (invp, pp);
  mpz_poly_setcoeff_si (f, 4, 60);
  mpz_poly_setcoeff_si (f, 3, 165063);
  mpz_poly_setcoeff_int64 (f, 2, (int64_t) 2561596016);
  mpz_poly_setcoeff_int64 (f, 1, (int64_t) -4867193837504);
  mpz_poly_setcoeff_int64 (f, 0, (int64_t) -9292909378109715);
  mpz_poly_setcoeff_si (P, 0, 0);
  mpz_poly_setcoeff_si (P, 1, 1); /* P = x */

  mpz_set_ui (a, 0);
  mpz_poly_power_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(Q->deg == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 1) == 0);
  mpz_poly_power_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(Q->deg == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 1) == 0);
  mpz_poly_power_mod_f_mod_mpz_Barrett (Q, P, f, a, pp, invp);
  ASSERT_ALWAYS(Q->deg == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 1) == 0);

  mpz_set_ui (a, 1);
  mpz_poly_power_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(mpz_poly_cmp (Q, P) == 0);
  mpz_poly_power_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(mpz_poly_cmp (Q, P) == 0);
  mpz_poly_power_mod_f_mod_mpz_Barrett (Q, P, f, a, pp, invp);
  ASSERT_ALWAYS(mpz_poly_cmp (Q, P) == 0);

  mpz_set_ui (a, 2);
  mpz_poly_power_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(Q->deg == 2);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 1) == 0);
  mpz_poly_power_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(Q->deg == 2);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 1) == 0);
  mpz_poly_power_mod_f_mod_mpz_Barrett (Q, P, f, a, pp, invp);
  ASSERT_ALWAYS(Q->deg == 2);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 1) == 0);

  mpz_set_ui (a, 3);
  mpz_poly_power_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[3], 1) == 0);
  mpz_poly_power_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[3], 1) == 0);
  mpz_poly_power_mod_f_mod_mpz_Barrett (Q, P, f, a, pp, invp);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[3], 1) == 0);

  mpz_set_ui (a, 4);
  mpz_poly_power_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 2081229567) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 3524154901) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 1102631344) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[3], 2362229259) == 0);
  mpz_poly_power_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 2081229567) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 3524154901) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 1102631344) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[3], 2362229259) == 0);
  mpz_poly_power_mod_f_mod_mpz_Barrett (Q, P, f, a, pp, invp);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 2081229567) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 3524154901) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 1102631344) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[3], 2362229259) == 0);

  mpz_set_ui (a, 999999);
  mpz_poly_power_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 4223801964) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 502704799) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 3358125388) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[3], 1722383279) == 0);
  mpz_poly_power_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 4223801964) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 502704799) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 3358125388) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[3], 1722383279) == 0);
  mpz_poly_power_mod_f_mod_mpz_Barrett (Q, P, f, a, pp, invp);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[0], 4223801964) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[1], 502704799) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[2], 3358125388) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (Q->coeff[3], 1722383279) == 0);

  mpz_clear (a);
  mpz_clear (pp);
  mpz_clear (invp);
  mpz_poly_clear (Q);
  mpz_poly_clear (P);
  mpz_poly_clear (f);
}

void
test_barrett_mod (unsigned long iter)
{
  mpz_t m, invm, a, b, c;
  unsigned long k;

  mpz_init (m);
  mpz_init (invm);
  mpz_init (a);
  mpz_init (b);
  mpz_init (c);
  while (iter--)
    {
      k = 2 + (lrand48 () % 128);
      do mpz_urandomb (m, state, k); while (mpz_sizeinbase (m, 2) < k);
      mpz_setbit (m, 0); /* make m odd */
      barrett_init (invm, m);

      mpz_urandomb (a, state, 2 * k);
      barrett_mod (b, a, m, invm);
      mpz_mod (c, a, m);
      ASSERT_ALWAYS (mpz_cmp (b, c) == 0);

      mpz_urandomb (a, state, 3 * k);
      barrett_mod (b, a, m, invm);
      mpz_mod (c, a, m);
      ASSERT_ALWAYS (mpz_cmp (b, c) == 0);

      mpz_urandomb (a, state, 4 * k);
      barrett_mod (b, a, m, invm);
      mpz_mod (c, a, m);
      ASSERT_ALWAYS (mpz_cmp (b, c) == 0);
    }
  mpz_clear (m);
  mpz_clear (invm);
  mpz_clear (a);
  mpz_clear (b);
  mpz_clear (c);
}

/* also test mpz_poly_base_modp_lift and mpz_poly_sizeinbase */
void
test_mpz_poly_base_modp_init (unsigned long iter)
{
  int p, l, k, K[32], d, i;
  mpz_poly_t f, *P, g;
  mpz_t pk;
  size_t s;

  mpz_poly_init (f, -1);
  mpz_poly_init (g, -1);
  mpz_init (pk);
  while (iter--)
    {
      p = lrand48 ();
      if (p < 2)
        p = 2;
      k = lrand48 () % 1024;
      if (k < 2)
        k = 2; /* ensures l > 0 */
      for (K[0] = k, l = 0; K[l] > 1; K[l+1] = (K[l] + 1) >> 1, l++);
      d = lrand48 () % 10;
      int m = (1 + (lrand48 () % 9)) * k;
      if (iter == 0) /* exercise bug found on 32-bit MinGW */
        {
          p = 2048;
          k = 119;
          d = 1;
          m = 833;
        }
      mpz_poly_random (f, d, m);
      s = mpz_poly_sizeinbase (f, d, 2);
      for (i = 0; i <= d; i++)
        ASSERT_ALWAYS(mpz_sizeinbase (f->coeff[i], 2) <= s);
      P = mpz_poly_base_modp_init (f, p, K, l);
      /* check f = P[0] + p^K[l]*P[1] + p^K[l-1]*P[2] + ... + p^K[1]*P[l] */
      for (i = 1; i <= l; i++)
        {
          mpz_ui_pow_ui (pk, p, K[l + 1 - i]);
          mpz_poly_base_modp_lift (P[0], P, i, pk);
        }
      ASSERT_ALWAYS(mpz_poly_cmp (f, P[0]) == 0);
      mpz_poly_base_modp_clear (P);
    }
  mpz_poly_clear (f);
  mpz_poly_clear (g);
  mpz_clear (pk);
}

void test_mpz_poly_is_root(unsigned long iter)
{
    mpz_t p, r;
    mpz_poly_t f, ell;
    mpz_init(p);
    mpz_init(r);
    mpz_poly_init(f, 10);
    mpz_poly_init(ell, 1);

    for( ; iter--; ) {
        mpz_poly_random(f, 10, 100);
        mpz_urandomb(p, state, 100);
        mpz_rrandomb(r, state, 100);
        mpz_poly_setcoeff_si(ell, 1, 1);
        mpz_neg(ell->coeff[0], r);
        mpz_poly_mul(f, f, ell);
        mpz_poly_mod_mpz(f, f, p, NULL);
        mpz_mod(r, r, p);
    }

    ASSERT_ALWAYS(mpz_poly_is_root(f, r, p));
    mpz_poly_clear(f);
    mpz_poly_clear(ell);
    mpz_clear(r);
    mpz_clear(p);
}

void test_mpz_poly_factor(unsigned long iter)
{
    mpz_t p;
    mpz_poly_factor_list lf;
    mpz_poly_t f;

    mpz_init(p);
    mpz_poly_init(f, -1);
    mpz_poly_factor_list_init(lf);

    mpz_set_ui(p, 13);
    mpz_poly_setcoeffs_ui_var(f, 21, 3, 8, 5, 5, 6, 1, 9, 4, 3, 3, 8, 7, 7, 7, 0, 12, 5, 11, 11, 1, 7, 10);
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 2);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 20);

    ASSERT_ALWAYS(mpz_poly_is_irreducible(lf->factors[1]->f, p));

    mpz_set_ui(p, 13);
    mpz_poly_setcoeffs_ui_var(f, 25, 0, 0, 0, 0, 5, 1, 1, 3, 7, 9, 7, 4, 11, 8, 2, 10, 8, 11, 0, 5, 12, 12, 1, 11, 2, 3);
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 6);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[0]->m == 4);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[2]->f->deg == 2);
    ASSERT_ALWAYS(lf->factors[3]->f->deg == 3);
    ASSERT_ALWAYS(lf->factors[4]->f->deg == 6);
    ASSERT_ALWAYS(lf->factors[5]->f->deg == 9);

    mpz_set_ui(p, 13);
    mpz_poly_setcoeffs_ui_var(f, 21, 0, 7, 4, 3, 6, 11, 3, 6, 12, 2, 5, 4, 7, 6, 4, 6, 12, 8, 0, 3, 3, 6);
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 5);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[2]->f->deg == 3);
    ASSERT_ALWAYS(lf->factors[2]->m == 2);
    ASSERT_ALWAYS(lf->factors[3]->f->deg == 5);
    ASSERT_ALWAYS(lf->factors[4]->f->deg == 8);

    mpz_set_ui(p, 3);
    mpz_poly_setcoeffs_ui_var(f, 20, 1, 1, 1, 2, 1, 0, 0, 1, 0, 0, 0, 1, 2, 0, 2, 0, 2, 0, 1, 1, 2);
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 4);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[1]->m == 3);
    ASSERT_ALWAYS(lf->factors[2]->f->deg == 6);
    ASSERT_ALWAYS(lf->factors[3]->f->deg == 10);

    mpz_set_ui(p, 3);
    mpz_poly_setcoeffs_ui_var(f, 20, 2, 0, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 0, 2, 2, 1, 1, 2, 0, 2, 2);
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 5);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[0]->m == 6);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 3);
    ASSERT_ALWAYS(lf->factors[2]->f->deg == 3);
    ASSERT_ALWAYS(lf->factors[3]->f->deg == 4);
    ASSERT_ALWAYS(lf->factors[4]->f->deg == 4);

    for( ; iter-- ; ) {
        // fprintf(stderr, "%lu ", iter);
        mpz_rrandomb(p, state, 20);
        mpz_nextprime(p, p);
        mpz_poly_random(f, 10, 10);
        mpz_poly_mod_mpz(f, f, p, NULL);
        // mpz_poly_fprintf(stderr, f);
        mpz_poly_factor(lf, f, p, state);
        mpz_poly_t g;
        mpz_poly_init(g, f->deg);
        mpz_poly_set_xi(g, 0);
        for(int i = 0 ; i < lf->size ; i++) {
            mpz_poly_with_m_ptr fx = lf->factors[i];
            mpz_poly_power_ui_mod_f_mod_mpz(fx->f, fx->f, NULL, fx->m, p); 
            mpz_poly_mul(g, g, fx->f);
            mpz_poly_mod_mpz(g, g, p, NULL);
        }
        mpz_poly_reduce_makemonic_mod_mpz(f, f, p);
        ASSERT_ALWAYS(mpz_poly_cmp(f, g) == 0);
        mpz_poly_clear(g);
    }

    mpz_poly_factor_list_clear(lf);
    mpz_clear(p);
    mpz_poly_clear(f);
}

void test_mpz_poly_trivialities()
{
    mpz_t a[5], p;
    mpz_poly_t f, g, q, r;
    mpz_poly_factor_list lf;
    int rc;

    mpz_init_set_ui(a[0], 1);
    mpz_init_set_ui(a[1], 4);
    mpz_init_set_ui(a[2], 6);
    mpz_init_set_ui(a[3], 4);
    mpz_init_set_ui(a[4], 1);
    mpz_init_set_ui(p, 13);
    mpz_poly_init(f, -1);
    mpz_poly_init(g, -1);
    mpz_poly_init(q, -1);
    mpz_poly_init(r, -1);
    mpz_poly_setcoeffs(f, a, 4);
    {
        mpz_poly_factor_list_init(lf);
        mpz_poly_factor(lf, f, p, state);
        ASSERT_ALWAYS(lf->factors[0]->m == 4);
        mpz_poly_swap(f, lf->factors[0]->f);
        mpz_poly_factor_list_clear(lf);
    }
    /* we expect to have f == x + 1 */
    mpz_poly_setcoeffs_ui_var(g, 1, 1, 1);
    mpz_poly_getcoeff(a[0], 0, f);
    ASSERT_ALWAYS(mpz_cmp_ui(a[0], 1) == 0);
    mpz_poly_getcoeff(a[1], 1, f);
    ASSERT_ALWAYS(mpz_cmp_ui(a[1], 1) == 0);
    mpz_poly_sub_mod_mpz(f, f, g, p);
    mpz_poly_set_zero(g);
    ASSERT_ALWAYS(mpz_poly_cmp(f, g) == 0);

    /* check that (x+1)^5 div x is irreducible mod 13 */
    mpz_poly_setcoeffs_ui_var(g, 1, 1, 1);
    mpz_poly_power_ui_mod_f_mod_mpz(g, g, NULL, 5, p);
    mpz_poly_div_xi(f, g, 6);
    ASSERT_ALWAYS(f->deg == -1);
    mpz_poly_div_xi(f, g, 0);
    ASSERT_ALWAYS(mpz_poly_cmp(f, g) == 0);
    mpz_poly_div_xi(g, g, 1);
    ASSERT_ALWAYS(mpz_poly_is_irreducible(g, p));

    /* check that x+1 - (x+1)^2 = -x^2-x */
    mpz_poly_set_zero(f);
    mpz_poly_add_ui(f, f, 1);
    mpz_poly_set_xi(g, 1);
    mpz_poly_add(f, f, g);
    mpz_poly_mul(g, f, f);
    mpz_poly_sub(f, f, g);
    mpz_poly_set_zero(g);
    mpz_poly_sub_ui(g, g, 1);
    mpz_poly_mul(f, f, g);
    mpz_poly_set_zero(g);
    mpz_set_ui(a[0], 1);
    mpz_poly_set_xi(g, 2);
    mpz_poly_setcoeff(g, 1, a[0]);
    ASSERT_ALWAYS(mpz_poly_cmp(f, g) == 0);

    /* multiply by zero */
    mpz_poly_random(f, 10, 10);
    mpz_poly_set_zero(g);
    mpz_poly_mul(f, f, g);
    ASSERT_ALWAYS(mpz_poly_cmp(f, g) == 0);

    /* div modulo non-prime N should get factor */
    mpz_set_ui(p, 143);
    mpz_poly_setcoeffs_ui_var(f, 2, 1, 1, 3);
    mpz_poly_setcoeffs_ui_var(g, 1, 1, 11);
    rc = mpz_poly_div_r(f, g, p);
    ASSERT_ALWAYS(rc == 0);
    mpz_poly_setcoeffs_ui_var(f, 2, 1, 1, 3);
    mpz_poly_setcoeffs_ui_var(g, 1, 1, 11);
    /* same test. Of course it doesn't exactly divide, but we do expect
     * failure nonetheless */
    rc = mpz_poly_divexact(f, f, g, p);
    ASSERT_ALWAYS(rc == 0);

    /* test div_qr */
    mpz_poly_setcoeffs_ui_var(f, 4, 1, 1, 1, 1, 3);
    mpz_poly_setcoeffs_ui_var(g, 2, 1, 2, 11);
    rc = mpz_poly_div_qr(q, r, f, g, p);
    ASSERT_ALWAYS(rc == 0);
    mpz_set_ui(p, 13);
    rc = mpz_poly_div_qr(q, r, f, g, p);
    mpz_poly_setcoeffs_ui_var(f, 2, 0, 11, 5);
    mpz_poly_setcoeffs_ui_var(g, 1, 1, 3);
    ASSERT_ALWAYS(rc);
    ASSERT_ALWAYS(mpz_poly_cmp(f, q) == 0);
    ASSERT_ALWAYS(mpz_poly_cmp(g, r) == 0);

    /* multiply by p, then reduce mod p */
    mpz_poly_random(f, 10, 10);
    mpz_poly_mul_mpz (f, f, p);
    mpz_poly_reduce_makemonic_mod_mpz(f, f, p);
    ASSERT_ALWAYS(f->deg < 0);

    /* a non-irreducible polynomial */
    mpz_poly_setcoeffs_ui_var(g, 2, 1, 2, 1);
    rc = mpz_poly_is_irreducible(g, p);
    ASSERT_ALWAYS(rc == 0);


    mpz_poly_clear(f);
    mpz_poly_clear(g);
    mpz_poly_clear(q);
    mpz_poly_clear(r);
    mpz_clear(a[4]);
    mpz_clear(a[3]);
    mpz_clear(a[2]);
    mpz_clear(a[1]);
    mpz_clear(a[0]);
    mpz_clear(p);
}

int
main (int argc, const char *argv[])
{
  unsigned long iter = 500;
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);

  test_polymodF_mul ();
  /* test_mpz_poly_roots_mpz (iter); */
  test_mpz_poly_sqr_mod_f_mod_mpz (iter);
  test_mpz_poly_fprintf ();
  test_mpz_poly_div_2_mod_mpz ();
  test_mpz_poly_derivative ();
  test_mpz_poly_power_mod_f_mod_ui ();
  test_barrett_mod (iter);
  test_mpz_poly_base_modp_init (iter);
  test_mpz_poly_is_root(iter);
  test_mpz_poly_factor(iter / 5);
  test_mpz_poly_trivialities ();
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}
