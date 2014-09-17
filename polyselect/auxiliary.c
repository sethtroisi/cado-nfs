/* Auxiliary routines for polynomial selection

Copyright 2008, 2009, 2010, 2013 Emmanuel Thome, Paul Zimmermann

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h> /* for DBL_MAX */
#include <math.h>
#include "gmp.h"
#include "macros.h" /* for ASSERT_ALWAYS */
#include "portability.h"
#include "utils.h"
#include "auxiliary.h"
#include "rootsieve.h"
#include "murphyE.h"

/* define OPTIMIZE_MP to perform computations in multiple-precision */
//#define OPTIMIZE_MP

//#define OPTIMIZE_LLL_LIST
#ifdef OPTIMIZE_LLL_LIST
typedef struct polylll_pq_t {
  mpz_t **f;
  mpz_t **g;
  mpz_t *fd;
  double *lognorm;
  double *alpha;
  double *skew;
  int used;
  int len;
  int degree;
} polylll_pq;
void new_polylll_pq (polylll_pq **queue, int len, int d);
void insert_polylll_pq (polylll_pq *queue, mpz_t *f, mpz_t *g,
                        double lognorm, double alpha, double skew);
void free_polylll_pq (polylll_pq **queue);
#endif

//#define DEBUG_OPTIMIZE_AUX

/* for the rotation, we try (j*x+k) for |k| <= 2^MAX_k */
int MAX_k = 16;

/************************* norm and skewness *********************************/

/* Same as L2_lognorm, but takes 'double' instead of 'mpz_t' as coefficients.
   Returns 1/2*log(int(int(F(r*cos(t)*s,r*sin(t))^2*r/s^d, r=0..1), t=0..2*Pi))
   (circular method). Cf Remark 3.2 in Kleinjung paper, Math. of Comp., 2006.

   Maple code for degree 2:
   f:=x->a2*x^2+a1*x+a0: d:=degree(f(x),x):
   int(int((f(s*cos(t)/sin(t))*(r*sin(t))^d)^2*r/s^d, r=0..1), t=0..2*Pi);
 */
static double
L2_lognorm_d (double_poly_srcptr p, double s)
{
  double n;
  double *a = p->coeff;
  unsigned int d = p->deg;

  ASSERT_ALWAYS(1 <= d && d <= 7);

  if (d == 1)
  {
    double a1, a0;
    a1 = a[1] * s;
    a0 = a[0];
    /* use circular integral (Sage code):
       var('a1,a0,x,y,r,s,t')
       f = a1*x+a0
       F = expand(f(x=x/y)*y)
       F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
       v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
       (s*v).expand().collect(pi)
    */
    n = a0 * a0 + a1 * a1;
    n = n * 0.785398163397448310; /* Pi/4 */
    return 0.5 * log (n / s);
  }
  else if (d == 2)
  {
    double a2, a1, a0;
    a2 = a[2] * s * s;
    a1 = a[1] * s;
    a0 = a[0];
    /* use circular integral (Sage code):
       var('a2,a1,a0,x,y,r,s,t')
       f = a2*x^2+a1*x+a0
       F = expand(f(x=x/y)*y^2)
       F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
       v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
       (s^2*v).expand().collect(pi)
    */
    n = 3.0 * (a2 * a2 + a0 * a0) + 2.0 * a0 * a2 + a1 * a1;
    n = n * 0.130899693899574704; /* Pi/24 */
    return 0.5 * log(n / (s * s));
  }
  else if (d == 3)
    {
      double a3, a2, a1, a0;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      /* use circular integral (Sage code):
         var('a3,a2,a1,a0,x,y,r,s,t')
         f = a3*x^3+a2*x^2+a1*x+a0
         F = expand(f(x=x/y)*y^3)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         (s^3*v).expand().collect(pi)
      */
      n = 5.0 * (a3 * a3 + a0 * a0) + 2.0 * (a3 * a1 + a0 * a2)
        + a1 * a1 + a2 * a2;
      n = n * 0.049087385212340519352; /* Pi/64 */
      return 0.5 * log(n / (s * s * s));
    }
  else if (d == 4)
    {
      double a4, a3, a2, a1, a0;

      a4 = a[4] * s * s * s * s;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      /* use circular integral (Sage code):
         var('a4,a3,a2,a1,a0,x,r,s,t')
         f = a4*x^4+a3*x^3+a2*x^2+a1*x+a0
         F = expand(f(x=x/y)*y^4)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         (s^4*v).expand().collect(pi)
      */
      n = 35.0 * (a4 * a4 + a0 * a0) + 10.0 * (a4 * a2 + a2 * a0)
        + 5.0 * (a3 * a3 + a1 * a1) + 6.0 * (a4 * a0 + a3 * a1)
        + 3.0 * a2 * a2;
      n = n * 0.0049087385212340519352; /* Pi/640 */
      return 0.5 * log(n / (s * s * s * s));
    }
  else if (d == 5)
    {
      double a5, a4, a3, a2, a1, a0, s2 = s * s, s3 = s2 * s, s4 = s3 * s,
        s5 = s4 * s;

      /*
        f := a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0:
        F := expand(y^5*subs(x=x/y,f));
        int(int(subs(x=x*s,F)^2/s^5, x=-1..1), y=-1..1);
       */
      a0 = a[0];
      a1 = a[1] * s;
      a2 = a[2] * s2;
      a3 = a[3] * s3;
      a4 = a[4] * s4;
      a5 = a[5] * s5;
      /* use circular integral (Sage code):
         var('a5,a4,a3,a2,a1,a0,x,r,s,t')
         f = a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0
         F = expand(f(x=x/y)*y^5)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         (s^5*v).expand().collect(pi)
      */
      n = 6.0 * (a3 * a1 + a1 * a5 + a4 * a2 + a0 * a4)
        + 14.0 * (a0 * a2 + a3 * a5) + 63.0 * (a0 * a0 + a5 * a5)
        + 7.0 * (a4 * a4 + a1 * a1) + 3.0 * (a3 * a3 + a2 * a2);
      n = n * 0.0020453077171808549730; /* Pi/1536 */
      return 0.5 * log(n / s5);
    }
  else if (d == 6)
    {
      double a6, a5, a4, a3, a2, a1, a0;

      a6 = a[6] * s * s * s * s * s * s;
      a5 = a[5] * s * s * s * s * s;
      a4 = a[4] * s * s * s * s;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      /* use circular integral (Sage code):
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=6; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         var('r,s,t')
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         (s^d*v).expand().collect(pi)
      */
      n = 231.0 * (a6 * a6 + a0 * a0) + 42.0 * (a6 * a4 + a2 * a0)
        + 21.0 * (a5 * a5 + a1 * a1) + 7.0 * (a4 * a4 + a2 * a2)
        + 14.0 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1)
        + 10.0 * (a6 * a0 + a5 * a1 + a4 * a2) + 5.0 * a3 * a3;
      n = n * 0.00043828022511018320850; /* Pi/7168 */
      return 0.5 * log(n / (s * s * s * s * s * s));
    }
  else /* d == 7 */
    {
      double a7, a6, a5, a4, a3, a2, a1, a0;

      a7 = a[7] * s * s * s * s * s * s * s;
      a6 = a[6] * s * s * s * s * s * s;
      a5 = a[5] * s * s * s * s * s;
      a4 = a[4] * s * s * s * s;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      /* use circular integral (Sage code):
         var('r,s,t,y')
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=7; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         (s^d*v).expand().collect(pi)
      */
      n = 429.0*(a0*a0+a7*a7) + 33.0*(a1*a1+a6*a6) + 66.0*(a0*a2+a5*a7)
        + 9*(a2*a2+a5*a5) + 18*(a1*a3+a0*a4+a4*a6+a3*a7) + 5*(a3*a3+a4*a4)
        + 10*(a2*a4+a1*a5+a3*a5+a0*a6+a2*a6+a1*a7);
      n = n * 0.000191747598485705154; /* Pi/16384 */
      return 0.5 * log(n / (s * s * s * s * s * s * s));
    }
}

/* Returns the logarithm of the L2-norm as defined by Kleinjung, i.e.,
   log(1/2 sqrt(int(int((F(sx,y)/s^(d/2))^2, x=-1..1), y=-1..1))).
   Since we only want to compare norms, we don't consider the log(1/2) term,
   and compute only 1/2 log(int(int(...))) [here the 1/2 factor is important,
   since it is added to the alpha root property term].

   Circular method: integrate over the unit circle.
*/
double
L2_lognorm (mpz_poly_ptr f, double s)
{
  double res;
  double_poly_t a;
  double a_coeffs[MAXDEGREE];

  a->coeff = a_coeffs;
  double_poly_set_mpz_poly (a, f);
  res = L2_lognorm_d (a, s);
  return res;
}

#ifdef OPTIMIZE_MP
/* The name is misleading. It returns the before-log part in
   1/2 log(int(int(...)))
*/
void
L2_lognorm_mp (mpz_poly_ptr f, mpz_t s, mpz_t norm)
{
  mpz_t n, tmp, tmp1, tmpsum;
  unsigned long i;
  unsigned int d = f->deg;

  mpz_init_set_ui (n, 1);
  mpz_init_set_ui (tmp, 1);
  mpz_init_set_ui (tmp1, 1);
  mpz_init_set_ui (tmpsum, 1);

  if (d != 6)
    {
      fprintf (stderr, "not yet implemented for degree %u\n", d);
      exit (1);
    }
  else {
    mpz_t a[d+1];
    mpz_init_set_ui (a[0], 1);
    for (i=1; i<=d; i++) {
      mpz_init (a[i]);
      mpz_mul (a[i], a[i-1], s);
    }
    for (i=0; i<=d; i++)
      mpz_mul (a[i], a[i], f->coeff[i]);

    // n = 231.0 * (a6 * a6 + a0 * a0)
    mpz_mul (tmp, a[6], a[6]);
    mpz_mul (tmp1, a[0], a[0]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_mul_ui(n, tmpsum, 231);
    //   + 42.0 * (a6 * a4 + a2 * a0)
    mpz_mul (tmp, a[6], a[4]);
    mpz_mul (tmp1, a[2], a[0]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_addmul_ui (n, tmpsum, 42);
    //   + 21.0 * (a5 * a5 + a1 * a1)
    mpz_mul (tmp, a[5], a[5]);
    mpz_mul (tmp1, a[1], a[1]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_addmul_ui (n, tmpsum, 21);
    // + 7.0 * (a4 * a4 + a2 * a2)
    mpz_mul (tmp, a[4], a[4]);
    mpz_mul (tmp1, a[2], a[2]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_addmul_ui (n, tmpsum, 7);
    // +  14.0 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1)
    mpz_mul (tmp, a[6], a[2]);
    mpz_mul (tmp1, a[5], a[3]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_mul (tmp, a[4], a[0]);
    mpz_mul (tmp1, a[3], a[1]);
    mpz_add (tmp, tmp, tmp1);
    mpz_add (tmpsum, tmp, tmpsum);
    mpz_addmul_ui (n, tmpsum, 14);
    // + 10.0 * (a6 * a0 + a5 * a1 + a4 * a2)
    mpz_mul (tmp, a[6], a[0]);
    mpz_mul (tmp1, a[5], a[1]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_mul (tmp, a[4], a[2]);
    mpz_add (tmpsum, tmpsum, tmp);
    mpz_addmul_ui (n, tmpsum, 10);
    //  + 5.0 * a3 * a3;
    mpz_mul (tmp, a[3], a[3]);
    mpz_addmul_ui (n, tmp, 5);

    for (i=0; i<=d; i++)
      mpz_clear (a[i]);
  }

  mpz_set (norm, n);
  mpz_clear (tmp);
  mpz_clear (tmp1);
  mpz_clear (tmpsum);
  mpz_clear (n);
}
#endif

static double
L2_skewness_deg6_approx (mpz_poly_ptr f MAYBE_UNUSED, double_poly_ptr ff,
                         double_poly_ptr dff, int prec)
{
  double *dfd = dff->coeff;
  double *fd = ff->coeff;
  double s, nc, s2, s4, s5, s6, a, b, c, smin, smax;
  int sign_changes = 0;
  double q[7], logmu, best_logmu = DBL_MAX, best_s = DBL_MAX;

  dfd[6] = 99.0 * fd[6] * fd[6];
  dfd[5] = 6.0 * (2.0 * fd[4] * fd[6] + fd[5] * fd[5]);
  dfd[4] = 2.0 * (fd[2] * fd[6] + fd[3] * fd[5]) + fd[4] * fd[4];
  dfd[2] = -2.0 * (fd[0] * fd[4] + fd[1] * fd[3]) - fd[2] * fd[2];
  dfd[1] = -6.0 * (2.0 * fd[0] * fd[2] + fd[1] * fd[1]);
  dfd[0] = -99.0 * fd[0] * fd[0];
  if (dfd[1] > 0)
    sign_changes ++; /* since dfd[0] < 0 */
  if (dfd[1] * dfd[2] < 0)
    sign_changes ++;
  if (dfd[2] * dfd[4] < 0)
    sign_changes ++; /* since dfd[3] = 0 */
  if (dfd[4] * dfd[5] < 0)
    sign_changes ++;
  if (dfd[5] < 0)
    sign_changes ++; /* since dfd[6] > 0 */
  /* since dfd[6] and dfd[0] have opposite signs, we have an odd number of
     roots on [0,+inf[. Moreover since dfd[3]=0, we can't have 5 positive
     roots, thus we have either 1 or 3. */

  q[6] = dfd[6];
  q[5] = (dfd[5] < 0) ? dfd[5] : 0.0;
  q[4] = (dfd[4] < 0) ? dfd[4] : 0.0;
  q[2] = (dfd[2] < 0) ? dfd[2] : 0.0;
  q[1] = (dfd[1] < 0) ? dfd[1] : 0.0;
  q[0] = dfd[0]; /* always negative */
  s = 1.0;
  while ((((((q[6]*s)+q[5])*s+q[4])*s*s+q[2])*s+q[1])*s+q[0] < 0)
    s = s + s;
  if (s == 1.0)
    {
      while ((((((q[6]*s)+q[5])*s+q[4])*s*s+q[2])*s+q[1])*s+q[0] > 0)
        s = s * 0.5;
      s = s + s;
    }
  smax = s;

  q[6] = dfd[6]; /* always positive */
  q[5] = (dfd[5] > 0) ? dfd[5] : 0.0;
  q[4] = (dfd[4] > 0) ? dfd[4] : 0.0;
  q[2] = (dfd[2] > 0) ? dfd[2] : 0.0;
  q[1] = (dfd[1] > 0) ? dfd[1] : 0.0;
  q[0] = dfd[0];
  s = smax;
  while ((((((q[6]*s)+q[5])*s+q[4])*s*s+q[2])*s+q[1])*s+q[0] > 0)
    s = s * 0.5;
  smin = s;

  /* positive roots are in [smin, smax] */

#define MULTS 2.0
  double v = -1.0;
  for (double t = MULTS * smin; t <= smax; t = MULTS * t)
    {
      /* invariant: q(smin) < 0 */
      v = (((((dfd[6]*t)+dfd[5])*t+dfd[4])*t*t+dfd[2])*t+dfd[1])*t+dfd[0];
      if (v <= 0)
        smin = t;
      else
        {
          /* the derivative has a root in [smin,b] and is increasing, thus
             the norm has a minimum in [smin,b], we refine it by dichotomy */
          a = smin;
          b = t;
          for (int i = 0; i < prec; i++)
            {
              c = (a + b) * 0.5;
              s2 = c * c;
              s4 = s2 * s2;
              s5 = s4 * c;
              s6 = s5 * c;

              nc = dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
                + dfd[2] * s2 + dfd[1] * c + dfd[0];
              if (nc > 0)
                b = c;
              else
                a = c;
            }
          s = sqrt ((a + b) * 0.5);
          if (sign_changes == 1)
            {
              best_s = s;
              break;
            }
          logmu = L2_lognorm_d (ff, s);
          if (logmu < best_logmu)
            {
              best_logmu = logmu;
              best_s = s;
            }
        }
    }

  ASSERT_ALWAYS (best_s < DBL_MAX);

  return best_s;
}

double
L2_skewness_deg6 (mpz_poly_ptr f MAYBE_UNUSED, double_poly_srcptr ff,
                  double_poly_srcptr dff MAYBE_UNUSED, int prec MAYBE_UNUSED)
{
  double s, logmu, logmu_min = DBL_MAX, s_min = 1.0;
  mpz_poly_t df;
  root_struct Roots[6];
  int i, k;

  mpz_poly_init (df, 6);
  for (i = 0; i < 6; i++)
    root_struct_init (Roots + i);

  /* Sage code:
     var('r,s,t,y')
     R.<x> = PolynomialRing(ZZ)
     S.<a> = InfinitePolynomialRing(R)
     d=6; f = SR(sum(a[i]*x^i for i in range(d+1)))
     F = expand(f(x=x/y)*y^d)
     F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
     v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
     v = (7168*v/pi).expand()
     dv = v.diff(s)
     dv = (dv*s^7/14).expand().collect(s)
     99*a6^2 * s^12 +
     6*(2*a4*a6 + a5^2) * s^10 +
     (2*a2*a6 + 2*a3*a5 + a4^2) * s^8 -
     (2*a0*a4 + 2*a1*a3 + a2^2) * s^4 -
     6*(2*a0*a2 + a1^2) * s^2 -
     99*a0^2
  */
#if 1 /* using numberOfRealRoots */
  mpz_mul (df->coeff[6], f->coeff[6], f->coeff[6]);
  mpz_mul_ui (df->coeff[6], df->coeff[6], 99); /* 99*a6^2 */
  mpz_mul (df->coeff[5], f->coeff[4], f->coeff[6]);
  mpz_mul_2exp (df->coeff[5], df->coeff[5], 1);
  mpz_addmul (df->coeff[5], f->coeff[5], f->coeff[5]);
  mpz_mul_ui (df->coeff[5], df->coeff[5], 6); /* 6*(2*a4*a6 + a5^2) */
  mpz_mul (df->coeff[4], f->coeff[2], f->coeff[6]);
  mpz_addmul (df->coeff[4], f->coeff[3], f->coeff[5]);
  mpz_mul_2exp (df->coeff[4], df->coeff[4], 1);
  mpz_addmul (df->coeff[4], f->coeff[4], f->coeff[4]); /*2*a2*a6+2*a3*a5+a4^2*/
  mpz_set_ui (df->coeff[3], 0);
  mpz_mul (df->coeff[2], f->coeff[0], f->coeff[4]);
  mpz_addmul (df->coeff[2], f->coeff[1], f->coeff[3]);
  mpz_mul_2exp (df->coeff[2], df->coeff[2], 1);
  mpz_addmul (df->coeff[2], f->coeff[2], f->coeff[2]);
  mpz_neg (df->coeff[2], df->coeff[2]); /* -(2*a0*a4+2*a1*a3+a2^2) */
  mpz_mul (df->coeff[1], f->coeff[0], f->coeff[2]);
  mpz_mul_2exp (df->coeff[1], df->coeff[1], 1);
  mpz_addmul (df->coeff[1], f->coeff[1], f->coeff[1]);
  mpz_mul_si (df->coeff[1], df->coeff[1], -6); /* -6*(2*a0*a2 + a1^2) */
  mpz_mul (df->coeff[0], f->coeff[0], f->coeff[0]);
  mpz_mul_si (df->coeff[0], df->coeff[0], -99); /* -99*a0^2 */

  df->deg = 6;
  k = numberOfRealRoots (df->coeff, 6, 0.0, 0, Roots);
  int kpos = 0;
  for (i = 0; i < k; i++)
    if (mpz_sgn (Roots[i].b) > 0)
      {
        kpos ++;
        s = rootRefine (Roots + i, df->coeff, 6, ldexp (1.0, -prec));
        s = sqrt (s);
        logmu = L2_lognorm_d (ff, s);
        if (logmu < logmu_min)
          {
            logmu_min = logmu;
            s_min = s;
          }
      }
#else /* using double_poly_compute_roots */
  double *dfd = dff->coeff;
  double *fd = ff->coeff;
  double roots[6];

  dfd[6] = 99.0 * fd[6] * fd[6];
  dfd[5] = 6.0 * ( 2.0 * fd[4] * fd[6] + fd[5] * fd[5] );
  dfd[4] = 2.0 * ( fd[2] * fd[6] + fd[3] * fd[5] ) + fd[4] * fd[4];
  dfd[3] = 0.0;
  dfd[2] = -2.0 * ( fd[0] * fd[4] + fd[1] * fd[3] ) - fd[2] * fd[2];
  dfd[1] = -6.0 * ( 2.0 * fd[0] * fd[2] + fd[1] * fd[1] );
  dfd[0] = -99.0 * fd[0] * fd[0];
  double B = double_poly_bound_roots (dff);
  k = double_poly_compute_roots (roots, dff, B);
  ASSERT_ALWAYS(k > 0);
  for (i = 0; i < k; i++)
    {
      s = sqrt (roots[i]);
      logmu = L2_lognorm_d (ff, s);
      if (logmu < logmu_min)
        {
          logmu_min = logmu;
          s_min = s;
        }
    }
#endif

  mpz_poly_clear (df);
  for (i = 0; i < 6; i++)
    root_struct_clear (Roots + i);

  return s_min;
}

/* Use derivative test, with ellipse regions */
double
L2_skewness (mpz_poly_ptr f, int prec)
{
  double_poly_t ff, df;
  double s = 0.0, a = 0.0, b = 0.0, c, nc, *fd, *dfd,
    s1, s2, s3, s4, s5, s6, s7;
  unsigned int d = f->deg;

  ASSERT_ALWAYS(1 <= d && d <= 7);

  double_poly_init (ff, d);
  double_poly_init (df, d);

  /* convert once for all to double's to avoid expensive mpz_get_d() */
  double_poly_set_mpz_poly (ff, f);
  fd = ff->coeff;
  dfd = df->coeff;
  if (d == 7)
    {
      /* Sage code:
         var('r,s,t,y')
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=7; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         v = (16384*v/pi).expand()
         dv = v.diff(s)
         dv = (dv*s^8).expand().collect(s)
         3003*a7^2*s^14 + 165*a6^2*s^12 + 330*a5*a7*s^12 + 27*a5^2*s^10
         + 54*a4*a6*s^10 + 54*a3*a7*s^10 + 5*a4^2*s^8 + 10*a3*a5*s^8
         + 10*a2*a6*s^8 + 10*a1*a7*s^8 - 5*a3^2*s^6 - 10*a2*a4*s^6
         - 10*a1*a5*s^6 - 10*a0*a6*s^6 - 27*a2^2*s^4 - 54*a1*a3*s^4
         - 54*a0*a4*s^4 - 165*a1^2*s^2 - 330*a0*a2*s^2 - 3003*a0^2
      */
      dfd[7] = 3003.0 * fd[7] * fd[7];
      dfd[6] = 165.0 * (fd[6] * fd[6] + 2.0 * fd[5] * fd[7]);
      dfd[5] = 27.0 * (fd[5]*fd[5] + 2.0*fd[4]*fd[6] + 2.0*fd[3]*fd[7]);
      dfd[4] = 5.0*(fd[4]*fd[4]+2.0*fd[3]*fd[5]+2.0*fd[2]*fd[6]+2.0*fd[1]*fd[7]);
      dfd[3] = 5.0*(fd[3]*fd[3]+2.0*fd[2]*fd[4]+2.0*fd[1]*fd[5]+2.0*fd[0]*fd[6]);
      dfd[2] = 27.0 * (fd[2]*fd[2] + 2.0*fd[1]*fd[3] + 2.0*fd[0]*fd[4]);
      dfd[1] = 165.0 * (fd[1]*fd[1] + 2.0*fd[0]*fd[2]);
      dfd[0] = 3003 * fd[0] * fd[0];
      s = 1.0;
      nc = dfd[7] + dfd[6] + dfd[5] + dfd[4] - dfd[3] - dfd[2] - dfd[1]
        - dfd[0];
      /* first isolate the minimum in an interval [s, 2s] by dichotomy */
      while (nc > 0)
        {
          s = 0.5 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          s6 = s3 * s3; /* s^12 */
          s7 = s6 * s1; /* s^14 */
          nc = dfd[7] * s7 + dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
            - dfd[3] * s3 - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      do
        {
          s = 2.0 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          s6 = s3 * s3; /* s^12 */
          s7 = s6 * s1; /* s^14 */
          nc = dfd[7] * s7 + dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
            - dfd[3] * s3 - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      while (nc < 0);

      /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
      a = (s == 2.0) ? 1.0 : 0.5 * s;
      b = s;
      /* use dichotomy to refine the root */
      while (prec--)
        {
          c = (a + b) * 0.5;
          s1 = c * c;
          s2 = s1 * s1;
          s3 = s2 * s1;
          s4 = s2 * s2;
          s5 = s4 * s1;
          s6 = s3 * s3;
          s7 = s6 * s1;

          nc = dfd[7] * s7 + dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
            - dfd[3] * s3 - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
          if (nc > 0)
            b = c;
          else
            a = c;
        }
    }
  else if (d == 6)
    {
      s = L2_skewness_deg6_approx (f, ff, df, prec);
      goto end;
    }
  else if (d == 5)
    {
      /* Sage code:
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=5; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         var('r,s,t')
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         v = (1536*v/pi).expand()
         dv = v.diff(s)
         dv = (dv*s^6/3).expand().collect(s)
      */
      dfd[5] = 105.0 * fd[5] * fd[5];
      dfd[4] = 7.0 * (2.0 * fd[3] * fd[5] + fd[4] * fd[4]);
      dfd[3] = 2.0 * (fd[1] * fd[5] + fd[2] * fd[4]) + fd[3] * fd[3];
      dfd[2] = 2.0 * (fd[0] * fd[4] + fd[1] * fd[3]) + fd[2] * fd[2];
      dfd[1] = 7.0 * (2.0 * fd[0] * fd[2] + fd[1] * fd[1]);
      dfd[0] = 105.0 * fd[0] * fd[0];
      s = 1.0;
      nc = dfd[5] + dfd[4] + dfd[3] - dfd[2] - dfd[1] - dfd[0];
      /* first isolate the minimum in an interval [s, 2s] by dichotomy */
      while (nc > 0)
        {
          s = 0.5 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          nc = dfd[5] * s5 + dfd[4] * s4 + dfd[3] * s3
            - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      do
        {
          s = 2.0 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          nc = dfd[5] * s5 + dfd[4] * s4 + dfd[3] * s3
            - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      while (nc < 0);

      /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
      a = (s == 2.0) ? 1.0 : 0.5 * s;
      b = s;
      /* use dichotomy to refine the root */
      while (prec--)
        {
          c = (a + b) * 0.5;
          s1 = c * c;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          nc = dfd[5] * s5 + dfd[4] * s4 + dfd[3] * s3
            - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
          if (nc > 0)
            b = c;
          else
            a = c;
        }
    }
  else if (d == 4)
    {
      /* Sage code:
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=4; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         var('r,s,t')
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         v = (640*v/pi).expand()
         dv = v.diff(s)
         dv = (dv*s^5/10).expand().collect(s)
      */
      dfd[4] = 14.0 * fd[4] * fd[4];
      dfd[3] = 2.0 * fd[2] * fd[4] + fd[3] * fd[3];
      dfd[1] = 2.0 * fd[0] * fd[2] + fd[1] * fd[1];
      dfd[0] = 14.0 * fd[0] * fd[0];
      s = 1.0;
      nc = dfd[4] + dfd[3] - dfd[1] - dfd[0];
      /* first isolate the minimum in an interval [s, 2s] by dichotomy */
      while (nc > 0)
        {
          s = 0.5 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          nc = dfd[4] * s4 + dfd[3] * s3 - dfd[1] * s1 - dfd[0];
        }
      do
        {
          s = 2.0 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          nc = dfd[4] * s4 + dfd[3] * s3 - dfd[1] * s1 - dfd[0];
        }
      while (nc < 0);

      /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
      a = (s == 2.0) ? 1.0 : 0.5 * s;
      b = s;
      /* use dichotomy to refine the root */
      while (prec--)
        {
          c = (a + b) * 0.5;
          s1 = c * c;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          nc = dfd[4] * s4 + dfd[3] * s3 - dfd[1] * s1 - dfd[0];
          if (nc > 0)
            b = c;
          else
            a = c;
        }
    }
  else if (d == 3)
    {
      /* Sage code:
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=3; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         var('r,s,t')
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         v = (64*v/pi).expand()
         dv = v.diff(s)
         dv = (dv*s^4).expand().collect(s)
      */
      dfd[3] = 15.0 * fd[3] * fd[3];
      dfd[2] = 2.0 * fd[1] * fd[3] + fd[2] * fd[2];
      dfd[1] = 2.0 * fd[0] * fd[2] + fd[1] * fd[1];
      dfd[0] = 15.0 * fd[0] * fd[0];
      s = 1.0;
      nc = dfd[3] + dfd[2] - dfd[1] - dfd[0];
      /* first isolate the minimum in an interval [s, 2s] by dichotomy */
      while (nc > 0)
        {
          s = 0.5 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          nc = dfd[3] * s3 + dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      do
        {
          s = 2.0 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          nc = dfd[3] * s3 + dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      while (nc < 0);

      /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
      a = (s == 2.0) ? 1.0 : 0.5 * s;
      b = s;
      /* use dichotomy to refine the root */
      while (prec--)
        {
          c = (a + b) * 0.5;
          s1 = c * c;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          nc = dfd[3] * s3 + dfd[2] * s2 - dfd[1] * s1 - dfd[0];
          if (nc > 0)
            b = c;
          else
            a = c;
        }
    }
  else if (d == 2)
    {
      /* Sage code:
         var('r,s,t,y')
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=2; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         v = (24*v/pi).expand()
         dv = v.diff(s)
         dv = (dv*s^3).expand().collect(s)
      */
      dfd[2] = 6.0 * fd[2] * fd[2];
      dfd[1] = 0.0;
      dfd[0] = 6.0 * fd[0] * fd[0];
      s = 1.0;
      nc = dfd[2] + dfd[1] - dfd[0];
      /* first isolate the minimum in an interval [s, 2s] by dichotomy */
      while (nc > 0)
        {
          s = 0.5 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          nc = dfd[2] * s2 + dfd[1] - dfd[0];
        }
      do
        {
          s = 2.0 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          nc = dfd[2] * s2 + dfd[1] - dfd[0];
        }
      while (nc < 0);
      /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
      a = (s == 2.0) ? 1.0 : 0.5 * s;
      b = s;
      /* use dichotomy to refine the root */
      while (prec--)
        {
          c = (a + b) * 0.5;
          s1 = c * c;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          nc = dfd[2] * s2 + dfd[1] - dfd[0];
          if (nc > 0)
            b = c;
          else
            a = c;
        }
    }
  else /* d == 1 */
    a = b = (fd[0] / fd[1] >= 0) ? fd[0] / fd[1] : -fd[0] / fd[1];

  s = (a + b) * 0.5;

 end:
  double_poly_clear (ff);
  double_poly_clear (df);

  return s;
}

double L2_skew_lognorm (mpz_poly_ptr f, int prec)
{
  return L2_lognorm (f, L2_skewness (f, prec));
}

#ifdef OPTIMIZE_MP

/* Use derivative test, with ellipse regions */
void
L2_skewness_derivative_mp (mpz_poly_ptr F, int prec, mpz_t skewness)
{
  mpz_t s, s1, s2, s3, s4, s5, s6, a, b, c, nc;
  mpz_init (s);
  mpz_init (s1);
  mpz_init (s2);
  mpz_init (s3);
  mpz_init (s4);
  mpz_init (s5);
  mpz_init (s6);
  mpz_init (a);
  mpz_init (b);
  mpz_init (c);
  mpz_init (nc);

  int i, d = F->deg;
  mpz_t *f = F->coeff;

  if (d == 6) {
    mpz_t df[d+1];
    mpz_t tmp;
    mpz_init_set_ui (tmp, 1);
    for (i=0; i<=d; i++)
      mpz_init (df[i]);

    // dfd[6] = 99.0 * fd[6] * fd[6];
    mpz_mul (df[6], f[6], f[6]);
    mpz_mul_ui(df[6], df[6], 99);
    // dfd[5] = 6.0 * ( 2.0 * fd[4] * fd[6] + fd[5] * fd[5] );
    mpz_mul (df[5], f[5], f[5]);
    mpz_mul (tmp, f[4], f[6]);
    mpz_add (df[5], df[5], tmp);
    mpz_add (df[5], df[5], tmp);
    mpz_mul_ui (df[5], df[5], 6);
    // dfd[4] = 2.0 * ( fd[2] * fd[6] + fd[3] * fd[5] ) + fd[4] * fd[4];
    mpz_mul (df[4], f[4], f[4]);
    mpz_mul (tmp, f[3], f[5]);
    mpz_add (df[4], df[4], tmp);
    mpz_add (df[4], df[4], tmp);
    mpz_mul (tmp, f[2], f[6]);
    mpz_add (df[4], df[4], tmp);
    mpz_add (df[4], df[4], tmp);
    //  dfd[2] = 2.0 * ( fd[0] * fd[4] + fd[1] * fd[3] ) + fd[2] * fd[2];
    mpz_mul (df[2], f[2], f[2]);
    mpz_mul (tmp, f[1], f[3]);
    mpz_add (df[2], df[2], tmp);
    mpz_add (df[2], df[2], tmp);
    mpz_mul (tmp, f[0], f[4]);
    mpz_add (df[2], df[2], tmp);
    mpz_add (df[2], df[2], tmp);
    // dfd[1] = 6.0 * ( 2.0 * fd[0] * fd[2] + fd[1] * fd[1] );
    mpz_mul (df[1], f[1], f[1]);
    mpz_mul (tmp, f[0], f[2]);
    mpz_add (df[1], df[1], tmp);
    mpz_add (df[1], df[1], tmp);
    mpz_mul_ui (df[1], df[1], 6);
    // dfd[0] = 99.0 * fd[0] * fd[0] ;
    mpz_mul (df[0], f[0], f[0]);
    mpz_mul_ui(df[0], df[0], 99);
    /*
    gmp_fprintf (stderr, "df[6]: %Zd\n", df[6]);
    gmp_fprintf (stderr, "df[5]: %Zd\n", df[5]);
    gmp_fprintf (stderr, "df[4]: %Zd\n", df[4]);
    gmp_fprintf (stderr, "df[2]: %Zd\n", df[2]);
    gmp_fprintf (stderr, "df[1]: %Zd\n", df[1]);
    gmp_fprintf (stderr, "df[0]: %Zd\n", df[0]);
    */

    mpz_set_si (nc, -1);
    mpz_set_ui (s, 1);

    /* first isolate the minimum in an interval [s, 2s] by dichotomy */
    while ( mpz_cmp_ui(nc, 0) < 0 ) {

      mpz_add (s, s, s); /* s = 2.0 * s */
      mpz_mul (s1, s, s); /* s^2 */
      mpz_mul (s2, s1, s1); /* s^4 */
      mpz_mul (s4, s2, s2); /* s^8 */
      mpz_mul (s5, s4, s1); /* s^10 */
      mpz_mul (s6, s5, s1); /* s^12 */

      /* nc = dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
         - dfd[2] * s2 - dfd[1] * s1 - dfd[0];         */
      mpz_mul (s6, s6, df[6]);
      mpz_mul (s5, s5, df[5]);
      mpz_mul (s4, s4, df[4]);
      mpz_mul (s2, s2, df[2]);
      mpz_mul (s1, s1, df[1]);
      mpz_add (nc, s6, s5);
      mpz_add (nc, nc, s4);
      mpz_sub (nc, nc, s2);
      mpz_sub (nc, nc, s1);
      mpz_sub (nc, nc, df[0]);

    }

    /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
    mpz_cdiv_q_2exp (a, s, 1);
    mpz_set (b, s);
    /* use dichotomy to refine the root */
    while (prec--)
    {
      mpz_add (tmp, a, b);
      mpz_cdiv_q_2exp (c, tmp, 1);
      mpz_mul (s1, c, c); //s1 = c * c;
      mpz_mul (s2, s1, s1); //s2 = s1 * s1;
      mpz_mul (s4, s2, s2); //s4 = s2 * s2;
      mpz_mul (s5, s4, s1); //s5 = s4 * s1;
      mpz_mul (s6, s5, s1); //s6 = s5 * s1;

      /* nc = dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
         - dfd[2] * s2 - dfd[1] * s1 - dfd[0]; */
      mpz_mul (s6, s6, df[6]);
      mpz_mul (s5, s5, df[5]);
      mpz_mul (s4, s4, df[4]);
      mpz_mul (s2, s2, df[2]);
      mpz_mul (s1, s1, df[1]);
      mpz_add (nc, s6, s5);
      mpz_add (nc, nc, s4);
      mpz_sub (nc, nc, s2);
      mpz_sub (nc, nc, s1);
      mpz_sub (nc, nc, df[0]);

      if (mpz_cmp_ui (nc, 0) > 0)
        mpz_set (b, c);
      else
        mpz_set (a, c);
    }

    mpz_clear (tmp);
    for (i=0; i<=d; i++)
      mpz_clear (df[i]);

  } // end
  else  {
    fprintf (stderr, "L2_skewness_derivative_mp not yet implemented for degree %d\n", d);
    exit (1);
  }

  mpz_add (s, a, b);
  mpz_cdiv_q_2exp (skewness, s, 1);

  mpz_clear (s);
  mpz_clear (s1);
  mpz_clear (s2);
  mpz_clear (s3);
  mpz_clear (s4);
  mpz_clear (s5);
  mpz_clear (s6);
  mpz_clear (a);
  mpz_clear (b);
  mpz_clear (c);
  mpz_clear (nc);
}
#endif

/************************** polynomial arithmetic ****************************/

/* g <- content(f) where deg(f)=d */
void
content_poly (mpz_t g, mpz_poly_ptr f)
{
  unsigned int i;
  unsigned int d = f->deg;

  ASSERT(d >= 1);

  mpz_gcd (g, f->coeff[0], f->coeff[1]);
  for (i = 2; i <= d; i++)
    mpz_gcd (g, g, f->coeff[i]);
}

/* v <- f(r), where f is of degree d */
void
eval_poly_ui (mpz_t v, mpz_t *f, unsigned int d, unsigned long r)
{
  int i;

  mpz_set (v, f[d]);
  for (i = d - 1; i >= 0; i--)
    {
      mpz_mul_ui (v, v, r);
      mpz_add (v, v, f[i]);
    }
}

/* v <- f'(r), where f is of degree d */
void
eval_poly_diff_ui (mpz_t v, mpz_t *f, unsigned int d, unsigned long r)
{
  int i;

  mpz_mul_ui (v, f[d], d);
  for (i = d - 1; i >= 1; i--)
    {
      mpz_mul_ui (v, v, r);
      mpz_addmul_ui (v, f[i], i); /* v <- v + i*f[i] */
    }
}

/* h(x) <- h(x + r/p), where the coefficients of h(x + r/p) are known to
   be integers */
static void
poly_shift_divp (mpz_t *h, unsigned int d, unsigned long r, unsigned long p)
{
  unsigned int i, k;
  mpz_t t;

  mpz_init (t);
  for (i = 1; i <= d; i++)
    for (k = d - i; k < d; k++)
      { /* h[k] <- h[k] + r/p h[k+1] */
        ASSERT (mpz_divisible_ui_p (h[k+1], p) != 0);
        mpz_divexact_ui (t, h[k+1], p);
        mpz_addmul_ui (h[k], t, r);
      }
  mpz_clear (t);
}

/********************* computation of alpha **********************************/

/* Auxiliary routine for special_valuation(), see below. It returns the
   average p-valuation of the polynomial f. Works recursively. */
double
special_val0 (mpz_poly_ptr f, unsigned long p)
{
  double v;
  mpz_t c,  *h;
  int v0;
  unsigned long *roots, r, r0;
  int i, d = f->deg, nroots;
  mpz_poly_t g, H;

  mpz_init (c);
  content_poly (c, f);
  for (v = 0.0; mpz_divisible_ui_p (c, p); v++, mpz_divexact_ui (c, c, p));
  v0 = (int) v;

  mpz_poly_init(g, d);
  g->deg = d;

  /* g <- f/p^v */
  if (v != 0.0)
    {
      mpz_ui_pow_ui (c, p, (unsigned long) v); /* p^v */
      for (i = 0; i <= d; i++)
        mpz_divexact (g->coeff[i], f->coeff[i], c);
    }
  else
    mpz_poly_set (g, f);

  mpz_poly_init (H, d);
  H->deg = d;
  h = H->coeff;
  /* first compute h(x) = g(px) */
  mpz_set_ui (c, 1);
  for (i = 0; i <= d; i++)
    {
      mpz_mul (h[i], g->coeff[i], c);
      mpz_mul_ui (c, c, p);
    }
  /* Search for roots of g mod p */
  ASSERT (d > 0);
  roots = (unsigned long*) malloc (d * sizeof (unsigned long));
  FATAL_ERROR_CHECK(roots == NULL, "not enough memory");

  nroots = mpz_poly_roots_ulong (roots, g, p);
  ASSERT (nroots <= d);
  for (r0 = 0, i = 0; i < nroots; i++)
    {
      r = roots[i];
      eval_poly_diff_ui (c, g->coeff, d, r);
      if (mpz_divisible_ui_p (c, p) == 0) /* g'(r) <> 0 mod p */
  v += 1.0 / (double) (p - 1);
      else /* hard case */
  {
    /* g(px+r) = h(x + r/p), thus we can go from h0(x)=g(px+r0)
       to h1(x)=g(px+r1) by computing h0(x + (r1-r0)/p).
       Warning: we can have h = f, and thus an infinite loop, when
       the p-valuation of f is d, and f has a single root r/(1-p) of
       multiplicity d.
       Moreover if f(x) = c*p^d*(x-r+b*p)^d, where c is coprime to p,
       then h(x) = f(p*x+r)/p^d = c*p^d*(x+b)^d, and most likely after
       at most p iterations we'll go back to f(x), thus we should avoid
       all cases where f(x) has a root of multiplicity d, but how to
       check that efficiently? And which value to return in such a case?
    */
    ASSERT_ALWAYS (r >= r0); /* the roots are sorted */
    poly_shift_divp (h, d, r - r0, p);
    r0 = r;
    if (v0 != d) /* this should catch all the cases where f(x) has a
        root of multiplicity d, but also more cases.
        In those cases we avoid an infinite loop, but the
        result is probably wrong. */
      v += special_val0 (H, p) / (double) p;
  }
    }
  free (roots);
  mpz_poly_clear (H);
  mpz_poly_clear (g);
  mpz_clear (c);

  return v;
}

/* Compute the average valuation of F(a,b) for gcd(a,b)=1, for a prime p
   dividing the discriminant of f, using the following algorithm from
   Guillaume Hanrot (which is some kind of p-adic variant of Uspensky's
   algorithm):

   val(f, p)
     return val0(f, p) * p / (p+1) + val0(f(1/(p*x))*(p*x)^d, p) * 1/(p+1)

   val0(f, p).
     v <- valuation (content(f), p);
     f <- f/p^v

     r <- roots mod p(f, p)

     for r_i in r do
         if f'(r_i) <> 0 mod p then v +=  1/(p-1).
         else
              f2 <- f(p*x + r_i)
              v += val0(f2, p) / p.
         endif
     endfor
     Return v.

A special case when:
(a) p^2 does not divide disc(f),
(b) p does not divide lc(f),
then the average valuation is (p q_p - 1)/(p^2 - 1), where q_p is the number
of roots of f mod p. When q_p=1, we get 1/(p+1).

Note: when p does not divide lc(f), the val0(f(1/(p*x))*(p*x)^d, p) call
always returns 0 in val(f,p).

Assumes p divides disc = disc(f), d is the degree of f.
*/
double
special_valuation (mpz_poly_ptr f, unsigned long p, mpz_t disc)
{
    double v;
    int p_divides_lc;
    int pvaluation_disc = 0;
    double pd = (double) p;
    int d = f->deg;

    if (mpz_divisible_ui_p(disc, p)) {
  mpz_t t;
  pvaluation_disc++;
  mpz_init(t);
  mpz_divexact_ui(t, disc, p);
  if (mpz_divisible_ui_p(t, p))
      pvaluation_disc++;
  mpz_clear(t);
    }

    p_divides_lc = mpz_divisible_ui_p(f->coeff[d], p);

    if (pvaluation_disc == 0) {
  /* easy ! */
  int e;
  e = mpz_poly_roots_ulong (NULL, f, p);
  if (p_divides_lc) {
      /* Or the discriminant would have valuation 1 at least */
      ASSERT(mpz_divisible_ui_p(f->coeff[d - 1], p) == 0);
      e++;
  }
  return (pd * e) / (pd * pd - 1);
    } else if (pvaluation_disc == 1) {
      /* special case where p^2 does not divide disc */
  int e;
  e = mpz_poly_roots_ulong (NULL, f, p);
        if (p_divides_lc)
          e ++;
  /* something special here. */
  return (pd * e - 1) / (pd * pd - 1);
    } else {
  v = special_val0(f, p) * pd;
  if (p_divides_lc) {
      /* compute g(x) = f(1/(px))*(px)^d, i.e., g[i] = f[d-i]*p^i */
      /* IOW, the reciprocal polynomial evaluated at px */
      mpz_poly_t G;
      mpz_t *g;
      mpz_t t;
      int i;

      mpz_poly_init (G, d);
      G->deg = d;
      g = G->coeff;
      mpz_init_set_ui(t, 1);  /* will contains p^i */
      for (i = 0; i <= d; i++) {
        mpz_mul(g[i], f->coeff[d - i], t);
        mpz_mul_ui(t, t, p);
      }
      v += special_val0(G, p);
      mpz_poly_clear (G);
      mpz_clear(t);
  }
  v /= pd + 1.0;
  return v;
    }
}

/* Compute the value alpha(F) from Murphy's thesis, page 49:
   alpha(F) = sum(prime p <= B, (1 - q_p*p/(p+1)) log(p)/(p-1))
   where q_p is the number of roots of F mod p, including the number of
   projective roots (i.e., the zeros of the reciprocal polynomial mod p).

   alpha(F) is an estimate of the average logarithm of the part removed
   from sieving, compared to a random integer.

   We want alpha as small as possible, i.e., alpha negative with a large
   absolute value. Typical good values are alpha=-4, -5, ...
*/
double
get_alpha (mpz_poly_ptr f, unsigned long B)
{
  double alpha, e;
  unsigned long p;
  mpz_t disc;
  unsigned int d = f->deg;

  /* for F linear, we have q_p = 1 for all p, thus
     alpha(F) = sum(prime p <= B, log(p)/(p^2-1)) ~ 0.569959993064325 */
  if (d == 1)
    return 0.569959993064325;

  mpz_init (disc);
  discriminant (disc, f->coeff, d);

  /* special_valuation returns the expected average exponent of p in F(a,b)
     for coprime a, b, i.e., e = q_p*p/(p^2-1), thus the contribution for p
     is (1/(p-1) - e) * log(p) */

  /* prime p=2 */
  e = special_valuation (f, 2, disc);
  alpha = (1.0 - e) * log (2.0);

  /* FIXME: generate all primes up to B and pass them to get_alpha */
  for (p = 3; p <= B; p += 2)
    if (ulong_isprime (p))
      {
        e = special_valuation (f, p, disc);
        alpha += (1.0 / (double) (p - 1) - e) * log ((double) p);
      }
  mpz_clear (disc);
  return alpha;
}

/* affine part of the special valution for polynomial f over p. */
double
special_valuation_affine (mpz_poly_ptr f, unsigned long p, mpz_t disc )
{
   double v;
   int pvaluation_disc = 0;
   double pd = (double) p;

   if (mpz_divisible_ui_p(disc, p)) {
      mpz_t t;
      pvaluation_disc++;
      mpz_init(t);
      mpz_divexact_ui(t, disc, p);
      if (mpz_divisible_ui_p(t, p))
         pvaluation_disc++;
      mpz_clear(t);
   }

   if (pvaluation_disc == 0) {
      /* case 1: root must be simple*/
      int e = 0;
      e = mpz_poly_roots_ulong (NULL, f, p);

      return (pd * e) / (pd * pd - 1);
   }
   /* else if (pvaluation_disc == 1) { */
   /*     /\* case 2: special case where p^2 does not divide disc *\/ */
   /*     int e = 0; */
   /*     e = mpz_poly_roots_ulong (NULL, f, p); */

   /*     /\* something special here. *\/ */
   /*     return (pd * e - 1) / (pd * pd - 1); */

   /* } */
   else {
      v = special_val0(f, p) * pd;
      v /= pd + 1.0;
      return v;
   }
}


/*
  Find biased alpha_projective for a poly f. It uses
  some hacks here which need to be changed in future.
  Until now, since this will only be done several
  times, hence the speed is not critical.

  Note that, the returned alpha is the  -val * log(p)
  biased part in the alpha. Hence, we can just add
  this to our affine part.
*/
double
get_biased_alpha_projective (mpz_poly_ptr f, unsigned long B)
{
   double alpha, e;
   unsigned long p;
   mpz_t disc;
   int d = f->deg;

   mpz_init (disc);
   discriminant (disc, f->coeff, d);

   /* prime p=2 */
   e = special_valuation (f, 2, disc) - special_valuation_affine (f, 2, disc);

   /* 1/(p-1) is counted in the affine part */
   alpha =  (- e) * log (2.0);

   /* FIXME: generate all primes up to B and pass them to get_alpha */
   for (p = 3; p <= B; p += 2)
      if (ulong_isprime (p)) {
         e = special_valuation(f, p, disc) - special_valuation_affine (f, p, disc);
         alpha += (- e) * log ((double) p);
      }

   mpz_clear (disc);

   return alpha;
}

#if 0
/*
  Similar to above, but for affine part.
*/
double
get_biased_alpha_affine (mpz_poly_ptr f, unsigned long B)
{
   double alpha, e;
   unsigned long p;
   mpz_t disc;

   mpz_init (disc);
   discriminant (disc, f, d);

   /* prime p=2 */
   e = special_valuation_affine (f, 2, disc);
   alpha =  (1.0 - e) * log (2.0);

   //printf ("\np: %u, val: %f, alpha: %f\n", 2, e, alpha);

   /* FIXME: generate all primes up to B and pass them to get_alpha */
   for (p = 3; p <= B; p += 2)
      if (ulong_isprime (p)) {
         e = special_valuation_affine (f, p, disc);
         alpha += (1.0 / (double) (p - 1) - e) * log ((double) p);
         //printf ("\np: %u, val: %f, alpha: %f\n", p, e, alpha);

      }
   mpz_clear (disc);
   return alpha;
}


/*
  Contribution from a particular multiple root r of the polynomial f
  over p. Note, r must also be a double root of f mod p.
*/
static double
average_valuation_affine_root (mpz_poly_ptr f, unsigned long p, unsigned long r )
{
   unsigned long v = 0UL;
   int i, j;
   mpz_t c, *fv;
   double val;

   mpz_init (c);

   /* init fv */
   fv = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
   if (fv == NULL) {
      fprintf (stderr, "Error, cannot allocate memory in average_valuation_affine_root.\n");
      exit (1);
   }

   for (i = 0; i <= d; i++)
      mpz_init_set (fv[i], f[i]);

   /* remove the p-valuations from fv */
   content_poly (c, f);
   while (mpz_divisible_ui_p(c, p)) {
      v += 1;
      for (i = 0; i <= d; i ++) {
         mpz_fdiv_q_ui (fv[i], f[i], p);
      }
   }

   /* first translate, then scale */
   for (i = d - 1; i >= 0; i--)
      for (j = i; j < d; j++)
         mpz_addmul_ui (fv[j], fv[j+1], r);
   /* t is p^i */
   mpz_set_ui(c, 1);
   for (i = 0; i <= d; i++) {
      mpz_mul(fv[i], fv[i], c);
      mpz_mul_ui(c, c, p);
   }

   /* now c is disc. */
   discriminant (c, fv, d);
   val = special_valuation_affine (fv, d, p, c);
   val = val / (double) p;

   /* clear */
   for (i = 0; i <= d; i++) {
      mpz_clear (fv[i]);
   }

   /* !!! REMEMBER THIS !!! */
   free (fv);
   mpz_clear(c);
   return val;
}
#endif


/**************************** rotation ***************************************/

/* D <- discriminant (f+k*g), which has degree d */
static void
discriminant_k (mpz_t *D, mpz_poly_ptr f, mpz_t m, mpz_t b)
{
  int i, j, k;
  uint32_t **M, pivot;
  int d = f->deg;

  ASSERT_ALWAYS(d <= 9);

  /* we first put in D[i] the value of disc(f + i*g) for 0 <= i <= d,
     thus if disc(f + k*g) = a[d]*k^d + ... + a[0], then
          D[0] = a[0]
          D[1] = a[0] + a[1] + ... + a[d]
          ...
          D[d] = a[0] + a[1]*d + ... + a[d]*d^d */

  discriminant (D[0], f->coeff, d); /* a[0] */
  for (i = 1; i <= d; i++)
    {
      /* add b*x - m */
      mpz_add (f->coeff[1], f->coeff[1], b);
      mpz_sub (f->coeff[0], f->coeff[0], m);
      discriminant (D[i], f->coeff, d);
    }

  /* initialize matrix coefficients */
  M = (uint32_t**) malloc ((d + 1) * sizeof(uint32_t*));
  for (i = 0; i <= d; i++)
    {
      M[i] = (uint32_t*) malloc ((d + 1) * sizeof(uint32_t));
      M[i][0] = 1;
      for (j = 1; j <= d; j++)
        M[i][j] = i * M[i][j-1];
    }

  for (j = 0; j < d; j++)
    {
      /* invariant: D[i] = M[i][0] * a[0] + ... + M[i][d] * a[d]
         with M[i][k] = 0 for k < j and k < i */
      for (i = j + 1; i <= d; i++)
        {
          /* eliminate M[i][j] */
          pivot = M[i][j] / M[j][j];
          mpz_submul_ui (D[i], D[j], pivot);
          for (k = j; k <= d; k++)
            M[i][k] -= pivot * M[j][k];
        }
    }

  /* now we have an upper triangular matrix */
  for (j = d; j > 0; j--)
    {
      for (k = j + 1; k <= d; k++)
        mpz_submul_ui (D[j], D[k], M[j][k]);
      ASSERT_ALWAYS(mpz_divisible_ui_p (D[j], M[j][j]));
      mpz_divexact_ui (D[j], D[j], M[j][j]);
    }

  /* restore the original f[] */
  mpz_submul_ui (f->coeff[1], b, d);
  mpz_addmul_ui (f->coeff[0], m, d);

  for (i = 0; i <= d; i++)
    free (M[i]);
  free (M);
}

/* replace f + k0 * x^t * (b*x - m) by f + k * x^t * (b*x - m), and return k */
long
rotate_aux (mpz_t *f, mpz_t b, mpz_t m, long k0, long k, unsigned int t)
{
  mpz_addmul_si (f[t + 1], b, k - k0);
  mpz_submul_si (f[t], m, k - k0);
  return k;
}

/* add k*x^t*g to f */
static void
rotate_auxg_si (mpz_t *f, mpz_t *g, long k, unsigned int t)
{
  mpz_addmul_si (f[t + 1], g[1], k);
  mpz_addmul_si (f[t], g[0], k);
}

/* replace f by f + k * x^t * (b*x + g0) */
void
rotate_auxg_z (mpz_t *f, const mpz_t b, const mpz_t g0, const mpz_t k, unsigned int t)
{
  mpz_addmul (f[t + 1], b, k);
  mpz_addmul (f[t], g0, k);
}

/* replace f by f - k * x^t * (b*x + g0) */
void
derotate_auxg_z (mpz_t *f, const mpz_t b, const mpz_t g0, const mpz_t k, unsigned int t)
{
  mpz_submul (f[t + 1], b, k);
  mpz_submul (f[t], g0, k);
}

unsigned long
rotate_area (long K0, long K1, long J0, long J1)
{
  return (unsigned long) (J1 - J0 + 1) * (unsigned long) (K1 - K0 + 1);
}

/* Use Emmanuel Thome's idea: assuming alpha(f + k*g) admits a Gaussian
   distribution with mean 'm' and standard deviation 's', then the probability
   that one polynomial has a value of alpha >= A is
   1/2*(1 - erf((A-m)/s/sqrt(2))), thus the probability that K polynomials
   have all their alpha values >= A is [1/2*(1 - erf((A-m)/s/sqrt(2)))]^K.
   For 'm' and 's' fixed, and a given K, we define the value of A for which
   this probability is 1/2 to be the expected value of A for K polynomials.

   We assume lognorm(f + k*g) + E(alpha(f + k*g)) is first decreasing, then
   increasing, then the optimal K corresponds to the minimum of that function.
*/
void
rotate_bounds (mpz_poly_ptr f, mpz_t b, mpz_t m, long *K0, long *K1,
               long *J0, long *J1, int verbose)
{
  /* exp_alpha[i] corresponds to K=2^i polynomials for m=0, s=1
     f := x -> 1/2*(1 - erf(x/sqrt(2)));
     seq(fsolve(f(x)^(2^k) = 1/2, x=-log[2](1.0+k)), k=0..30);
   */
  double exp_alpha[] = {0.0 /* K=1 */, -0.5449521356 /* 2 */,
                        -0.9981488825 /* 4 */, -1.385198061 /* 8 */,
                        -1.723526050 /* 16 */, -2.025111894 /* 32 */,
                        -2.298313131 /* 64 */, -2.549067054 /* 128 */,
                        -2.781676726 /* 256 */, -2.999326227 /* 512 */,
                        -3.204420841 /* 1024 */, -3.398814100 /* 2048 */,
                        -3.583961388 /* 4096 */, -3.761025425 /* 8192 */,
                        -3.930949902 /* 16384 */, -4.094511733 /* 32768 */,
                        -4.252358774 /* 65536 */, -4.405037486 /* 131072 */,
                        -4.553013560 /* 262144 */, -4.696687518 /* 524288 */,
                        -4.836406692 /* 2^20 */, -4.972474538 /* 2^21 */,
                        -5.105157963 /* 2^22 */, -5.234693169 /* 2^23 */,
                        -5.361290351 /* 2^24 */, -5.485137511 /* 2^25 */,
                        -5.606403590 /* 2^26 */, -5.725241052 /* 2^27 */,
                        -5.841788041 /* 2^27 */, -5.956170181 /* 2^29 */,
                        DBL_MAX};
  int i;
  long k0 = 0, j0 = 0;
  double lognorm, alpha, E0, E, best_E;
  double skewness = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
  long jmax = (long) ((double) (1L << MAX_k) / skewness);
  unsigned long max_area = 1UL << MAX_k;

#define MARGIN 0.12 /* we allow a small error margin in the expected lognorm
                       + alpha values, to get a larger search range */

  E0 = L2_lognorm (f, skewness);

  *K0 = -2;
  *J0 = -16;
  *J1 = 16;

  /* look for negative k: -2, -4, -8, ... */
  best_E = E0;
  for (i = 1; rotate_area (*K0, -*K0, *J0, *J1) < max_area; i++, *K0 *= 2)
    {
      k0 = rotate_aux (f->coeff, b, m, k0, *K0, 0);
      lognorm = L2_skew_lognorm (f, SKEWNESS_DEFAULT_PREC);
      alpha = exp_alpha[i];
      E = lognorm + alpha;
      if (E < best_E + MARGIN)
        {
          if (E < best_E)
            best_E = E;
        }
      else
        break;
    }
  /* go back to k=0 */
  k0 = rotate_aux (f->coeff, b, m, k0, 0, 0);
  *K1 = -*K0;

  /* now try negative j: -1, -3, -7, ... */
  for (i++; exp_alpha[i] != DBL_MAX && rotate_area (*K0, *K1, *J0, -*J0) < max_area; i++, *J0 = 2 * *J0 - 1)
    {
      j0 = rotate_aux (f->coeff, b, m, j0, *J0, 1);
      lognorm = L2_skew_lognorm (f, SKEWNESS_DEFAULT_PREC);
      alpha = exp_alpha[i];
      E = lognorm + alpha;
      if (E < best_E + MARGIN)
        {
          if (E < best_E)
            best_E = E;
        }
      else
        break;
      if (1 - 2 * *J0 > jmax)
        break;
    }
  *J1 = -*J0;

  if (verbose > 0)
    fprintf (stderr, "# Rotate bounds: %ld <= j <= %ld, %ld <= k <= %ld\n",
             *J0, *J1, *K0, *K1);

  /* rotate back to j=k=0 */
  rotate_aux (f->coeff, b, m, k0, 0, 0);
  rotate_aux (f->coeff, b, m, j0, 0, 1);
}

/* Return the smallest value of lognorm + alpha(f + (j*x+k)*(b*x-m)) for
   j and k small enough such that the norm does not increase too much, and
   modify f[] accordingly.
   The parameter "multi" means that several polynomials are wanted. If
   multi=0 or 1, then only 1 polynomial is returned (classical behavior).
   Otherwise, multi polynomials are stored in jmin and kmin (that
   must be initialized arrays with at least multi elements). This option
   might be useful for Coppersmith variant (a.k.a. MNFS).
   In the multi case, the smallest of the returned values of lognorm + alpha
   is returned (and f[] accordingly).
   Warning: the caller is responsible to update the skewness if needed.

   FIXME: it would be better to pass the polynomial g to rotate, instead of
   b and m.
   */
double
rotate (mpz_poly_ptr f, unsigned long alim, mpz_t m, mpz_t b,
        long *jmin, long *kmin, int multi, int verbose)
{
  mpz_array_t *D;
  long K0, K1, J0, J1, k0, k, i, j, j0, bestk;
  double *A, alpha, lognorm, best_alpha = DBL_MAX, best_lognorm = DBL_MAX;
  double corr = 0.0;
  double alpha0;
  unsigned long p;
  double *best_E = NULL; /* set to NULL to avoid warning... */
  double time_alpha = 0.0, time_norm = 0.0;
  int d = f->deg;
  unsigned long pp;
  double one_over_pm1, logp, average_alpha = 0.0;

  /* allocate best_E, to store the best (lognorm+alpha) in multi mode */
  if (multi > 1)
    {
      best_E = (double *) malloc (multi * sizeof (double));
      for (i = 0; i < multi; ++i)
        best_E[i] = DBL_MAX;
    }

  /* allocate D(k) = disc(f + (j*x+k)*g, x) */
  D = alloc_mpz_array (d + 1);

  /* compute range for k */
  rotate_bounds (f, b, m, &K0, &K1, &J0, &J1, verbose);
  ASSERT_ALWAYS(K0 <= 0 && 0 <= K1);

  /* allocate sieving zone for computing alpha */
  A = (double*) malloc ((K1 + 1 - K0) * sizeof (double));
  j0 = k0 = 0; /* the current coefficients f[] correspond to f+(j*x+k)*g */

  *jmin = *kmin = 0;

  alpha0 = get_alpha (f, alim); /* value of alpha without rotation */

  ASSERT_ALWAYS(J0 < 0 && 0 < J1);
  for (j = 0;;)
    {
      /* we consider j=0, 1, ..., J1, then J0, J0+1, ..., -1 */
      j0 = rotate_aux (f->coeff, b, m, j0, j, 1);
      /* go back to k=0 for the discriminant */
      k0 = rotate_aux (f->coeff, b, m, k0, 0, 0);
      /* D(k) = disc(f + (j*x+k)*g, x) (j is now fixed) */
      discriminant_k (D->data, f, m, b);

      for (k = K0; k <= K1; k++)
  A[k - K0] = 0.0; /* A[k - K0] will store the value alpha(f + k*g) */

  for (p = 2; p <= alim; p += 1 + (p & 1))
    if (ulong_isprime (p))
      {
        int i;
        /* We skip primes which divide all coefficients of f, since then
           f mod p is zero. This can only happen when p divides N, which is
           silly in GNFS, but maybe the user is stupid. */
        for (i = d; i >= 0 && mpz_divisible_ui_p (f->coeff[i], p); i--);
        if (i < 0)
          continue;

  if (k0 != 0)
    k0 = rotate_aux (f->coeff, b, m, k0, 0, 0);

  time_alpha -= seconds ();

  one_over_pm1 = 1.0 / (double) (p - 1);
  logp = log ((double) p);
  for (pp = p; pp <= alim; pp *= p)
    {
      /* Murphy (page 48) defines cont_p(F) = q_p*p/(p^2-1)
         = q_p*p/(p+1)*(1/p+1/p^2+...)
         The contribution for p^k is thus q_p*p/(p+1)/p^k. */
      alpha = logp / (double) (p + 1) * (double) p / (double) pp;
      /* the following is the average contribution for a prime not
         dividing the discriminant, cf alpha.sage, function alpha_p.
         We take it into account only for p, not for p^2, p^3, ... */
      if (p == pp)
        average_alpha += logp * one_over_pm1;
      /* we do not consider the projective roots here, since the
         corresponding correction will be considered separately for each
         row below */
      /* + alpha_p_projective (f, d, (D->data)[0], p); */
      update_table (f->coeff, d, m, b, A, K0, K1, pp, alpha);
    }

  time_alpha += seconds ();
      } /* end of loop on primes p */

  /* determine the best alpha in each row */
  bestk = K0;
  for (k = K0 + 1; k <= K1; k++)
    if (A[k - K0] < A[bestk - K0])
      bestk = k;

  /* Correction to apply to the current row (takes into account the projective
     roots). FIXME: we are lazy here, we should only consider the contribution
     from the projective roots. */
  k0 = rotate_aux (f->coeff, b, m, k0, bestk, 0);
  corr = get_alpha (f, alim) - A[bestk - K0];

  if (verbose > 1)
    fprintf (stderr, "# best alpha for j=%ld: k=%ld with %f\n",
             j, bestk, A[bestk - K0] + corr);

  /* now finds the best lognorm+alpha */
  time_norm -= seconds ();
  for (k = K0; k <= K1; k++)
    {
      alpha = A[k - K0] + corr;
      if (alpha < best_alpha + 2.0)
        {
          /* check lognorm + alpha < best_lognorm + best_alpha */

          /* translate from k0 to k */
          k0 = rotate_aux (f->coeff, b, m, k0, k, 0);
          lognorm = L2_skew_lognorm (f, SKEWNESS_DEFAULT_PREC);
          if (multi <= 1) {
              if (lognorm + alpha < best_lognorm + best_alpha) {
                  best_lognorm = lognorm;
                  best_alpha = alpha;
                  *kmin = k;
                  *jmin = j;
              }
          } else { /* multi mode */
              /* Rem: best_lognorm + best_alpha is the worse of the
                 preselected */
              double newE = lognorm + alpha;
              if (newE < best_E[multi-1]) {
                  int ii;
                  /* find position; assume list of preselected is sorted */
                  for (ii = 0; ii < multi; ++ii) {
                      if (best_E[ii] > newE)
                          break;
                  }
                  ASSERT_ALWAYS(ii < multi);
                  /* insert */
                  for (i = multi - 1; i > ii; --i) {
                      kmin[i] = kmin[i-1];
                      jmin[i] = jmin[i-1];
                      best_E[i] = best_E[i-1];
                  }
                  kmin[ii] = k;
                  jmin[ii] = j;
                  best_E[ii] = newE;
              }
          }
        }
    }
  time_norm += seconds ();

  j++;
  if (j > J1)
    j = J0;
  else if (j == 0)
    break;

    } /* end of loop on j */

  /* we now have f + (j0*x+k0)*(bx-m) and we want f + (jmin*x+kmin)*(bx-m),
     thus we have to add ((jmin-j0)*x+(kmin-k0)*(bx-m) */
  /* if you are in multi-mode, we use the best polynomial */
  rotate_aux (f->coeff, b, m, k0, *kmin, 0);
  rotate_aux (f->coeff, b, m, j0, *jmin, 1);

  if ((verbose > 0) && (multi <= 1))
    {
      fprintf (stderr, "# Rotate by ");
      if (*jmin != 0)
        {
          if (*jmin == -1)
            fprintf (stderr, "-");
          else if (*jmin != 1)
            fprintf (stderr, "%ld*", *jmin);
          fprintf (stderr, "x");
          if (*kmin >= 0)
            fprintf (stderr, "+");
        }
      fprintf (stderr, "%ld: alpha improved from %1.2f to %1.2f (alpha %1.2fs, norm %1.2fs)\n", *kmin, alpha0, best_alpha, time_alpha, time_norm);
    }

  if (verbose && (multi > 1)) {
      fprintf(stderr, "Found the following polynomials  (j, k, E):\n");
      for (i = 0; i < multi; ++i) {
          fprintf(stderr, "  %ld\t%ld\t%1.2f\n", jmin[i], kmin[i], best_E[i]);
      }
  }

  free (A);

  clear_mpz_array (D);

  {
      double ret_val = best_lognorm + best_alpha;
      if (multi>1) {
          ret_val = best_E[0]; /* we return the smallest */
          free(best_E);
      }
      return ret_val;
  }
}


/*
   Print f, g only.
   Note: it's a backend for print_cadopoly().
*/
void
print_cadopoly_fg (FILE *fp, mpz_t *f, int df, mpz_t *g, int dg, mpz_t n )
{
   int i;

   /* n */
   gmp_fprintf (fp, "\nn: %Zd\n", n);

   /* Y[i] */
   for (i = dg; i >= 0; i--)
     gmp_fprintf (fp, "Y%d: %Zd\n", i, g[i]);

   /* c[i] */
   for (i = df; i >= 0; i--)
     gmp_fprintf (fp, "c%d: %Zd\n", i, f[i]);
}


/*
   Print f, g only, lognorm, skew, alpha, MurphyE.
   Note:  it's a backend for print_cadopoly_extra().
*/
double
print_cadopoly (FILE *fp, cado_poly p)
{
   unsigned int nroots = 0;
   double alpha, alpha_proj, logmu, e;
   mpz_poly_t F, G;

   F->coeff = p->alg->coeff;
   F->deg = p->alg->deg;
   G->coeff = p->rat->coeff;
   G->deg = p->rat->deg;

   /* print f, g only*/
   print_cadopoly_fg (fp, F->coeff, F->deg, G->coeff, G->deg, p->n);

#ifdef DEBUG
   fprintf (fp, "# ");
   fprint_polynomial (fp, F->coeff, F->deg);
   fprintf (fp, "# ");
   fprint_polynomial (fp, G->coeff, G->deg);
#endif

   fprintf (fp, "skew: %1.3f\n", p->skew);

   if (G->deg > 1)
   {
    logmu = L2_lognorm (G, p->skew);
    alpha = get_alpha (G, ALPHA_BOUND);
    alpha_proj = get_biased_alpha_projective (G, ALPHA_BOUND);
    nroots = numberOfRealRoots (G->coeff, G->deg, 0, 0, NULL);
    fprintf (fp, "# lognorm: %1.2f, alpha: %1.2f (proj: %1.2f), E: %1.2f, "
                 "nr: %u\n", logmu, alpha, alpha_proj, logmu + alpha, nroots);
   }

   logmu = L2_lognorm (F, p->skew);
   alpha = get_alpha (F, ALPHA_BOUND);
   alpha_proj = get_biased_alpha_projective (F, ALPHA_BOUND);
   nroots = numberOfRealRoots (F->coeff, F->deg, 0, 0, NULL);
   fprintf (fp, "# lognorm: %1.2f, alpha: %1.2f (proj: %1.2f), E: %1.2f, "
                "nr: %u\n", logmu, alpha, alpha_proj, logmu + alpha, nroots);

   e = MurphyE (p, bound_f, bound_g, area, MURPHY_K);
   fprintf (fp, "# MurphyE(Bf=%.2e,Bg=%.2e,area=%.2e)=%.2e\n",
        bound_f, bound_g, area, e);

   return e;
}


/*
   Print f, g, lognorm, skew, alpha, MurphyE, REV, time duration.
*/
void
print_cadopoly_extra (FILE *fp, cado_poly cpoly, int argc, char *argv[], double st)
{
   int i;

   print_cadopoly (fp, cpoly);
   /* extra info */
   fprintf (fp, "# generated by %s: %s", CADO_REV, argv[0]);
   for (i = 1; i < argc; i++)
      fprintf (fp, " %s", argv[i]);
   fprintf (fp, " in %.2fs\n", (seconds () - st));
}


/*
  Call print_cadopoly, given f, g and return MurphyE.
*/
double
print_poly_fg (mpz_poly_ptr f, mpz_t *g, mpz_t N, int mode)
{
   double e;
   int i;
   int d = f->deg;

   cado_poly cpoly;
   cado_poly_init(cpoly);
   for (i = 0; i < (d + 1); i++)
      mpz_set(cpoly->alg->coeff[i], f->coeff[i]);
   for (i = 0; i < 2; i++)
      mpz_set(cpoly->rat->coeff[i], g[i]);
   mpz_set(cpoly->n, N);
   cpoly->skew = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
   cpoly->alg->deg = d;
   cpoly->rat->deg = 1;

   if (mode == 1)
     {
       e = print_cadopoly (stdout, cpoly);
       fflush(stdout);
     }
   else
     e = MurphyE (cpoly, bound_f, bound_g, area, MURPHY_K);

   cado_poly_clear (cpoly);
   return e;
}

/* f <- f(x+k), g <- g(x+k) */
void
do_translate_z (mpz_poly_ptr f, mpz_t *g, const mpz_t k)
{
  int i, j;
  int d = f->deg;

  for (i = d - 1; i >= 0; i--)
    for (j = i; j < d; j++)
      mpz_addmul (f->coeff[j], f->coeff[j+1], k);
  mpz_addmul (g[0], g[1], k);
}

/* f <- f(x-k), g <- g(x-k) */
void
do_detranslate_z (mpz_poly_ptr f, mpz_t *g, const mpz_t k)
{
  int i, j;
  int d = f->deg;

  for (i = d - 1; i >= 0; i--)
    for (j = i; j < d; j++)
      mpz_submul (f->coeff[j], f->coeff[j+1], k);
  mpz_submul (g[0], g[1], k);
}

#define MAX_ROT 4 /* we consider rotation by x^k for 0 <= k < MAX_ROT */

/* Use rotation and translation to find a polynomial with smaller norm
   (local minimum). Modify f and g accordingly.
   If flag = 0, use only translation.
   If flag = 1, use rotation and translation.
   If flag = 2, use only rotation.
*/
void
optimize_aux (mpz_poly_ptr f, mpz_t *g, int verbose, int flag, int max_iter)
{
  mpz_t kt, k[MAX_ROT], l; /* current offset */
  /* compute total translation and rotation, such that current polynomials are
     [f+(khitot*x^2+lamtot*x+mutot)*g](x+k) and g(x+k) */
  mpz_t tmp;
  int changed[MAX_ROT], changedt;
  double logmu00, logmu0, logmu, skew;
  int prec = SKEWNESS_DEFAULT_PREC;
  int count = 0;
  int d = f->deg, i, maxi;
  mpz_poly_t G;

  maxi = 1 + (d >= 5);
  ASSERT_ALWAYS(maxi < MAX_ROT);

  G->coeff = g;
  G->deg = 1;
  skew = L2_skewness (f, prec);
  logmu00 = logmu0 = L2_lognorm (f, skew);
  for (i = 0; i <= maxi; i++)
    mpz_init_set_ui (k[i], 1);
  mpz_init_set_ui (kt, 1);
  mpz_init (l);
  mpz_init (tmp);

#define GUARD 0.001
  /* initialize k[i] to the smallest power of two that increases the lognorm
     by GUARD with respect to the initial polynomial, to avoid being stuck
     in a very flat region */
  for (i = 0; i <= maxi; i++)
    {
      while (1)
        {
          rotate_auxg_z (f->coeff, g[1], g[0], k[i], i);
          logmu = L2_lognorm (f, skew);
          derotate_auxg_z (f->coeff, g[1], g[0], k[i], i); /* restore f */
          if (logmu > logmu0 + GUARD)
	    {
	      mpz_neg (l, k[i]);
	      rotate_auxg_z (f->coeff, g[1], g[0], l, i);
	      logmu = L2_lognorm (f, skew);
	      derotate_auxg_z (f->coeff, g[1], g[0], l, i);
	      if (logmu > logmu0 + GUARD)
		break;
	    }
          mpz_mul_2exp (k[i], k[i], 1);
        }
    }

  /* initialize kt likewise */
  while (1)
    {
      do_translate_z (f, g, kt);
      logmu = L2_lognorm (f, skew);
      do_detranslate_z (f, g, kt);
      if (logmu > logmu0 + GUARD)
	{
	  mpz_neg (l, kt);
	  do_translate_z (f, g, l);
	  logmu = L2_lognorm (f, skew);
	  do_detranslate_z (f, g, l);
	  if (logmu > logmu0 + GUARD)
	    break;
	}
      mpz_mul_2exp (kt, kt, 1);
    }

  while (1)
    {
      changed[0] = changedt = changed[2] = changed[1] = 0;

      if (flag != 2)
        {
      /* first try translation by kt */
      do_translate_z (f, g, kt); /* f(x+kt) */
      logmu = L2_skew_lognorm (f, prec);
      if (logmu < logmu0)
        {
          changedt = 1;
          logmu0 = logmu;
        }
      else
        {
          mpz_mul_si (l, kt, -2); /* l = -2*kt */
          do_translate_z (f, g, l); /* f(x-kt) */
          logmu = L2_skew_lognorm (f, prec);
          if (logmu < logmu0)
          {
              changedt = 1;
              logmu0 = logmu;
            }
          else
            do_translate_z (f, g, kt); /* original f */
        }
        }

      for (i = 0; (i <= maxi) && flag; i++)
        {
          /* do rotation by k[i]*x^i*g */
          rotate_auxg_z (f->coeff, g[1], g[0], k[i], i);
          logmu = L2_skew_lognorm (f, prec);
          if (logmu < logmu0)
            {
              changed[i] = 1;
              logmu0 = logmu;
            }
          else
            {
              mpz_mul_si (l, k[i], -2); /* l = -2*k2 */
              rotate_auxg_z (f->coeff, g[1], g[0], l, i); /* f - k[i]*x^i*g */
              logmu = L2_skew_lognorm (f, prec);
              if (logmu < logmu0)
                {
                  changed[i] = 1;
                  logmu0 = logmu;
                }
              else
                rotate_auxg_z (f->coeff, g[1], g[0], k[i], i); /*previous f*/
            }
        }

      if (changedt == 0 &&
          changed[0] == 0 && changed[2] == 0 && changed[1] == 0 &&
          mpz_cmp_ui (kt, 1) == 0 && mpz_cmp_ui (k[0], 1) == 0 &&
          (maxi < 2 || mpz_cmp_ui (k[2], 1) == 0) && mpz_cmp_ui (k[1], 1) == 0)
        break;

      if (changedt == 1)
        mpz_mul_2exp (kt, kt, 1);       /* kt <- 2*kt */
      else if (mpz_cmp_ui (kt, 1) > 0)
        mpz_div_2exp (kt, kt, 1);       /* kt <- kt/2 */

      for (i = 0; i <= maxi; i++)
        {
          if (changed[i] == 1)
            mpz_mul_2exp (k[i], k[i], 1);       /* k[i] <- 2*k[i] */
          else if (mpz_cmp_ui (k[i], 1) > 0)
            mpz_div_2exp (k[i], k[i], 1);       /* k[i] <- k[i]/2 */
        }

      if (count++ > max_iter) /* abort in cases where the descent procedure
                                 takes too long to converge; 200 iterations
                                 seems largely enough in most cases */
        break;
    }

  if (verbose > 1)
    {
      gmp_fprintf (stderr, "# ad=%Zd: optimized lognorm from %.2f to %.2f in %d iterations\n",
                   f->coeff[d], logmu00, logmu0, count);
      if (verbose > 2)
        {
          fprintf (stderr, "# "); mpz_poly_fprintf (stderr, f);
          fprintf (stderr, "# "); mpz_poly_fprintf (stderr, G);
        }
    }

  for (i = 0; i <= maxi; i++)
    mpz_clear (k[i]);
  mpz_clear (kt);
  mpz_clear (l);
  mpz_clear (tmp);
}

#ifdef OPTIMIZE_MP
/*
   use mpz
*/
void
optimize_aux_mp (mpz_poly_ptr f, mpz_t *g, int verbose, int use_rotation)
{
  mpz_t kt, k2, k1, k, l; /* current offset */
  mpz_t ktot, khitot, lamtot, mutot; /* compute total translation and rotation,
                                        such that current polynomial are
                                        [f+(khitot*x^2+lamtot*x+mutot)*g](x+k)
                                        and g(x+k) */
  mpz_t tmp, tmp0, logmu00, logmu0, logmu, skew, skew0;
  int changed, changedt, changed2, changed1;
  int prec = SKEWNESS_DEFAULT_PREC;
  int count = 0;
  int d = f->deg;
  mpz_poly_t G;

  G->coeff = g;
  G->deg = 1;
  mpz_init (tmp);
  mpz_init (tmp0);
  mpz_init (skew);
  mpz_init (skew0);
  mpz_init (logmu);
  mpz_init (logmu0);
  mpz_init (logmu00);
  mpz_init (l);
  mpz_init (ktot);
  mpz_init (lamtot);
  mpz_init (khitot);
  mpz_init (mutot);
  mpz_init_set_ui (k, 1);
  mpz_init_set_ui (k2, 1);
  mpz_init_set_ui (k1, 1);
  mpz_init_set_ui (kt, 1);

  // init norm
  L2_skewness_derivative_mp (f, prec, skew0);
  L2_lognorm_mp (f, skew0, logmu0);
  mpz_set (logmu00, logmu0);

  while (1)
  {
    changed = changedt = changed2 = changed1 = 0;

    /* first try translation by kt */
    do_translate_z (f, g, kt); /* f(x+kt) */
    mpz_add (ktot, ktot, kt);
    L2_skewness_derivative_mp (f, prec, skew);
    L2_lognorm_mp (f, skew, logmu);

    mpz_pow_ui (tmp, skew0, d); // s_old^6
    mpz_mul (tmp, tmp, logmu);
    mpz_pow_ui (tmp0, skew, d); // s^6
    mpz_mul (tmp0, tmp0, logmu0);
    if (mpz_cmp (tmp, tmp0) < 0)
    {
      changedt = 1;
      mpz_set (logmu0, logmu);
      mpz_set (skew0, skew);
    }
    else
    {
      mpz_mul_si (l, kt, -2); /* l = -2*kt */
      do_translate_z (f, g, l); /* f(x-kt) */
      mpz_add (ktot, ktot, l);
      L2_skewness_derivative_mp (f, prec, skew);
      L2_lognorm_mp (f, skew, logmu);

      mpz_pow_ui (tmp, skew0, d); // s_old^6
      mpz_mul (tmp, tmp, logmu);
      mpz_pow_ui (tmp0, skew, d); // s^6
      mpz_mul (tmp0, tmp0, logmu0);
      if (mpz_cmp (tmp, tmp0) < 0)
      {
        changedt = 1;
        mpz_set (logmu0, logmu);
        mpz_set (skew0, skew);
      }
      else
      {
        do_translate_z (f, g, kt); /* original f */
        mpz_add (ktot, ktot, kt);
      }
    }

    /* then do rotation by k2*x^2*g if d >= 6 */
    if (d >= 6 && use_rotation)
    {
      rotate_auxg_z (f->coeff, g[1], g[0], k2, 2);
      mpz_mul (tmp, ktot, k2);
      mpz_submul_ui (lamtot, tmp, 2);
      mpz_addmul (mutot, tmp, ktot);
      mpz_add (khitot, khitot, k2);
      L2_skewness_derivative_mp (f, prec, skew);
      L2_lognorm_mp (f, skew, logmu);

      mpz_pow_ui (tmp, skew0, d); // s_old^6
      mpz_mul (tmp, tmp, logmu);
      mpz_pow_ui (tmp0, skew, d); // s^6
      mpz_mul (tmp0, tmp0, logmu0);
      if (mpz_cmp (tmp, tmp0) < 0)
      {
        changed2 = 1;
        mpz_set (logmu0, logmu);
        mpz_set (skew0, skew);
      }
      else
      {
        mpz_mul_si (l, k2, -2); /* l = -2*k2 */
        rotate_auxg_z (f->coeff, g[1], g[0], l, 2); /* f - k2*x^2*g */
        mpz_mul (tmp, ktot, l);
        mpz_submul_ui (lamtot, tmp, 2);
        mpz_addmul (mutot, tmp, ktot);
        mpz_add (khitot, khitot, l);
        L2_skewness_derivative_mp (f, prec, skew);
        L2_lognorm_mp (f, skew, logmu);

        mpz_pow_ui (tmp, skew0, d); // s_old^6
        mpz_mul (tmp, tmp, logmu);
        mpz_pow_ui (tmp0, skew, d); // s^6
        mpz_mul (tmp0, tmp0, logmu0);
        if (mpz_cmp (tmp, tmp0) < 0)
        {
          changed2 = 1;
          mpz_set (logmu0, logmu);
          mpz_set (skew0, skew);
        }
        else
        {
          rotate_auxg_z (f->coeff, g[1], g[0], k2, 2); /* previous f */
          mpz_mul (tmp, ktot, k2);
          mpz_submul_ui (lamtot, tmp, 2);
          mpz_addmul (mutot, tmp, ktot);
          mpz_add (khitot, khitot, k2);
        }
      }
    }

    if (use_rotation)
    {
      /* then do rotation by k1*x*g */
      rotate_auxg_z (f->coeff, g[1], g[0], k1, 1);
      mpz_submul (mutot, ktot, k1);
      mpz_add (lamtot, lamtot, k1);
      L2_skewness_derivative_mp (f, prec, skew);
      L2_lognorm_mp (f, skew, logmu);

      mpz_pow_ui (tmp, skew0, d); // s_old^6
      mpz_mul (tmp, tmp, logmu);
      mpz_pow_ui (tmp0, skew, d); // s^6
      mpz_mul (tmp0, tmp0, logmu0);
      if (mpz_cmp (tmp, tmp0) < 0)
      {
        changed1 = 1;
        mpz_set (logmu0, logmu);
        mpz_set (skew0, skew);
      }
      else
      {
        mpz_mul_si (l, k1, -2); /* l = -2*k1 */
        rotate_auxg_z (f->coeff, g[1], g[0], l, 1); /* f - k1*x*g */
        mpz_submul (mutot, ktot, l);
        mpz_add (lamtot, lamtot, l);
        L2_skewness_derivative_mp (f, prec, skew);
        L2_lognorm_mp (f, skew, logmu);

        mpz_pow_ui (tmp, skew0, d); // s_old^6
        mpz_mul (tmp, tmp, logmu);
        mpz_pow_ui (tmp0, skew, d); // s^6
        mpz_mul (tmp0, tmp0, logmu0);
        if (mpz_cmp (tmp, tmp0) < 0)
        {
          changed1 = 1;
          mpz_set (logmu0, logmu);
          mpz_set (skew0, skew);
        }
        else
        {
          rotate_auxg_z (f->coeff, g[1], g[0], k1, 1); /* previous f */
          mpz_submul (mutot, ktot, k1);
          mpz_add (lamtot, lamtot, k1);
        }
      }

      /* then do rotation by k*g */
      rotate_auxg_z (f->coeff, g[1], g[0], k, 0);
      mpz_add (mutot, mutot, k);
      L2_skewness_derivative_mp (f, prec, skew);
      L2_lognorm_mp (f, skew, logmu);

      mpz_pow_ui (tmp, skew0, d); // s_old^6
      mpz_mul (tmp, tmp, logmu);
      mpz_pow_ui (tmp0, skew, d); // s^6
      mpz_mul (tmp0, tmp0, logmu0);
      if (mpz_cmp (tmp, tmp0) < 0)
      {
        changed = 1;
        mpz_set (logmu0, logmu);
        mpz_set (skew0, skew);
      }
      else
      {
        mpz_mul_si (l, k, -2); /* l = -2*k */
        rotate_auxg_z (f->coeff, g[1], g[0], l, 0); /* f - k*g */
        mpz_add (mutot, mutot, l);
        L2_skewness_derivative_mp (f, prec, skew);
        L2_lognorm_mp (f, skew, logmu);

        mpz_pow_ui (tmp, skew0, d); // s_old^6
        mpz_mul (tmp, tmp, logmu);
        mpz_pow_ui (tmp0, skew, d); // s^6
        mpz_mul (tmp0, tmp0, logmu0);
        if (mpz_cmp (tmp, tmp0) < 0)
        {
          changed = 1;
          mpz_set (logmu0, logmu);
          mpz_set (skew0, skew);
        }
        else
        {
          rotate_auxg_z (f->coeff, g[1], g[0], k, 0); /* previous f */
          mpz_add (mutot, mutot, k);
        }
      }
    } /* use_rotation */

    if (changedt == 1)
      mpz_mul_2exp (kt, kt, 1);       /* kt <- 2*kt */
    else if (mpz_cmp_ui (kt, 1) > 0)
      mpz_div_2exp (kt, kt, 1);       /* kt <- kt/2 */
    if (changed == 1)
      mpz_mul_2exp (k, k, 1);       /* k <- 2*k */
    else if (mpz_cmp_ui (k, 1) > 0)
      mpz_div_2exp (k, k, 1);       /* k <- k/2 */
    if (changed2 == 1)
      mpz_mul_2exp (k2, k2, 1);       /* k2 <- 2*k2 */
    else if (mpz_cmp_ui (k2, 1) > 0)
      mpz_div_2exp (k2, k2, 1);       /* k2 <- k2/2 */
    if (changed1 == 1)
      mpz_mul_2exp (k1, k1, 1);       /* k1 <- 2*k1 */
    else if (mpz_cmp_ui (k1, 1) > 0)
      mpz_div_2exp (k1, k1, 1);       /* k1 <- k1/2 */
    if (changedt == 0 && changed == 0 && changed2 == 0 && changed1 == 0 &&
        mpz_cmp_ui (kt, 1) == 0 && mpz_cmp_ui (k, 1) == 0 &&
        mpz_cmp_ui (k2, 1) == 0 && mpz_cmp_ui (k1, 1) == 0)
      break;
    if (count++ > 10000) /* avoid an infinite loop due to the random
                            choices when logmu=logmu0 */
      break;
  }

  if (verbose > 0)
  {
    gmp_fprintf (stderr, "# ad=%Zd: optimized lognorm from %Zd to %Zd\n",
                 f[d], logmu00, logmu0);
    gmp_fprintf (stderr, "# (rotation %Zd*x^2+%Zd*x+%Zd, translation %Zd)\n",
                 khitot, lamtot, mutot, ktot);
    if (verbose > 1)
    {
      fprintf (stderr, "# ");
      mpz_poly_fprintf (stderr, f);
      fprintf (stderr, "# ");
      mpz_poly_fprintf (stderr, G);
    }
  }
  mpz_clear (tmp);
  mpz_clear (tmp0);
  mpz_clear (logmu);
  mpz_clear (logmu0);
  mpz_clear (logmu00);
  mpz_clear (skew);
  mpz_clear (skew0);
  mpz_clear (k2);
  mpz_clear (k1);
  mpz_clear (k);
  mpz_clear (kt);
  mpz_clear (l);
  mpz_clear (ktot);
  mpz_clear (khitot);
  mpz_clear (lamtot);
  mpz_clear (mutot);
}
#endif

/* puts in h[2], ..., h[0] the coefficients (in k) of x^(d-2) in f(x+k) */
static void
fdminus2_translated (mpz_poly_ptr h, mpz_poly_ptr f)
{
  unsigned int d = f->deg;

  mpz_mul_ui (h->coeff[2], f->coeff[d], (d * (d-1)) / 2);
  mpz_mul_ui (h->coeff[1], f->coeff[d-1], d-1);
  mpz_set (h->coeff[0], f->coeff[d-2]);
  h->deg = 2;
}

/* given f(x) of degree d, the coefficient of x^(d-3) in f(x+k) is a
   polynomial of degree 3 in k, say h[3]*k^3 + h[2]*k^2 + h[1]*k + h[0].
   Put in h[3], ..., h[0] those coefficients. */
static void
fdminus3_translated (mpz_poly_ptr h, mpz_poly_ptr f)
{
  unsigned int d = f->deg;

  mpz_mul_ui (h->coeff[3], f->coeff[d], (d * (d-1) * (d-2)) / 6);
  mpz_mul_ui (h->coeff[2], f->coeff[d-1], ((d-1) * (d-2)) / 2);
  mpz_mul_ui (h->coeff[1], f->coeff[d-2], d-2);
  mpz_set (h->coeff[0], f->coeff[d-3]);
  h->deg = 3;
}

#if 0
/* puts in h[4], ..., h[0] the coefficients (in k) of x^(d-4) in f(x+k) */
static void
fdminus4_translated (mpz_poly_ptr h, mpz_poly_ptr f)
{
  unsigned int d = f->deg;

  mpz_mul_ui (h->coeff[4], f->coeff[d], (d * (d-1) * (d-2) * (d-3)) / 24);
  mpz_mul_ui (h->coeff[3], f->coeff[d-1], ((d-1) * (d-2) * (d-3)) / 6);
  mpz_mul_ui (h->coeff[2], f->coeff[d-2], ((d-2) * (d-3)) / 2);
  mpz_mul_ui (h->coeff[1], f->coeff[d-3], d-3);
  mpz_set (h->coeff[0], f->coeff[d-4]);
  h->deg = 4;
}
#endif

#define MAX_EFFORT 1000

int optimize_effort = 1; /* should be <= MAX_EFFORT */

/* LMAX is the number of l-values kept for which the coefficient
   of degree d-2 in [f(x)+l*x^(d-3)*g(x)](x+k) is small for some k.
   The optimization time is roughly linear in LMAX. */
int LMAX;

/* MAXL is such that we consider only values of l in [-MAXL,MAXL].
   The optimization time does not depend much on MAXL, but taking a
   too large value will find values of k which are large, and might give
   worse results. Taking MAXL = 2*LMAX seems to work well. */
int MAXL;

void
set_optimize_effort (int e)
{
  ASSERT_ALWAYS(e >= 1);
  optimize_effort = e;
}

#define LMAX_MAX MAX_EFFORT
#define MAXY 0
/* we need at least LMAX+1 entries in each array (for a sentinel and to
   store l=0 at the end) */
long best_l[LMAX_MAX+1], best_l2[LMAX_MAX+1];
double best_v[LMAX_MAX+1], best_v2[LMAX_MAX+1];

void
insert_l (long *L, double *V, int *n, long l, double v)
{
  int i = *n;

  while (0 < i && v < V[i-1])
    {
      V[i] = V[i-1];
      L[i] = L[i-1];
      i--;
    }
  if (i == LMAX)
    return;
  V[i] = v;
  L[i] = l;
  if (*n < LMAX)
    (*n)++;
}

/* Put in best_l[0..n-1] the best values of l such that the coefficient
   of degree d-2 of (f+l*x^(d-3)*g)(x=x+k) is small, i.e.,
   d*(d-1)/2*c[d]*k^2 + (d-1)*c[d-1]*k + c[d-2] + l*g1.
   To ensure we have at two real roots in k, we should have:
   (d-1)^2*c[d-1]^2 - 2*d*(d-1)*c[d]*(c[d-2] + l*g1) >= 0
   l*g1 <= (d-1)*c[d-1]^2/(2*d*c[d]) - c[d-2].
   Return the number n of values stored, with n <= LMAX.
*/
int
find_good_l (mpz_poly_ptr f, mpz_t *g)
{
  int d = f->deg, n = 0;
  mpz_t lmax;
  long l, l0, l1;
  double v0, v1;

  mpz_init (lmax);
  mpz_mul (lmax, f->coeff[d-1], f->coeff[d-1]);
  mpz_mul_ui (lmax, lmax, d-1);
  mpz_fdiv_q (lmax, lmax, f->coeff[d]);
  mpz_fdiv_q_ui (lmax, lmax, 2*d);
  mpz_sub (lmax, lmax, f->coeff[d-2]);
  mpz_fdiv_q (lmax, lmax, g[1]);
  if (mpz_fits_slong_p (lmax))
    {
      mpz_t az, bz, cz, delta, t, sdelta, k, v;
      mpz_init (az);
      mpz_init (bz);
      mpz_init (cz);
      mpz_init (delta);
      mpz_init (sdelta);
      mpz_init (t);
      mpz_init (k);
      mpz_init (v);
      l0 = mpz_get_si (lmax);
      if (l0 > MAXL)
        {
          l0 = MAXL;
          l1 = -MAXL;
        }
      else if (-MAXL <= l0)
        l1 = -MAXL;
      else
        l1 = l0 - MAXL / 2;
      mpz_mul_ui (az, f->coeff[d], (d * (d-1)) / 2);
      mpz_mul_ui (bz, f->coeff[d-1], d-1);
      mpz_mul_si (cz, g[1], l0);
      mpz_add (cz, cz, f->coeff[d-2]);
      mpz_mul (delta, az, cz);
      mpz_mul_si (delta, delta, -4);
      mpz_addmul (delta, bz, bz);
      ASSERT_ALWAYS(mpz_cmp_ui (delta, 0) >= 0);
      mpz_mul (t, az, g[1]);
      mpz_mul_ui (t, t, 4);
      for (l = l0; l >= l1; l--)
        {
          /* invariant: cz = c[d-2] + l * g1
                        delta = bz^2 - 4*az*cz */
          mpz_sqrt (sdelta, delta);

          /* first root, rounded downwards */
          mpz_sub (k, sdelta, bz);
          mpz_fdiv_q (k, k, az);
          mpz_fdiv_q_ui (k, k, 2);
          mpz_mul (v, az, k);
          mpz_add (v, v, bz);
          mpz_mul (v, v, k);
          mpz_add (v, v, cz);
          v0 = fabs (mpz_get_d (v));

          /* second root, rounded downwards */
          mpz_add (k, bz, sdelta);
          mpz_neg (k, k);
          mpz_fdiv_q (k, k, az);
          mpz_fdiv_q_ui (k, k, 2);
          mpz_mul (v, az, k);
          mpz_add (v, v, bz);
          mpz_mul (v, v, k);
          mpz_add (v, v, cz);
          v1 = fabs (mpz_get_d (v));

          /* insert only the smallest absolute value */
          insert_l (best_l, best_v, &n, l, (v0 < v1) ? v0 : v1);

          /* update invariants */
          mpz_sub (cz, cz, g[1]);
          mpz_add (delta, delta, t);
        }
      mpz_clear (az);
      mpz_clear (bz);
      mpz_clear (cz);
      mpz_clear (delta);
      mpz_clear (sdelta);
      mpz_clear (t);
      mpz_clear (k);
      mpz_clear (v);
    }
  mpz_clear (lmax);
  return n;
}

/* Put in best_l2[0..n-1] the best values of l such that the coefficient
   of degree d-3 of (f+l*x^(d-3)*g)(x+k) is small, i.e.,
   d*(d-1)*(d-2)/6*c[d]*k^3 + (d-1)*(d-2)/2*c[d-1]*k^2 +
   (d-2)*(c[d-2]+l*g1)*k + c[d-3] + l*g0.
   Return the number n of values stored, with n <= LMAX.
*/
int
find_good_l2 (mpz_poly_ptr f, mpz_t *g)
{
  int d = f->deg, n = 0, nr, j;
  long l;
  mpz_poly_t c3k;
  mpz_t r[3], v;

  mpz_poly_init (c3k, 3);
  c3k->deg = 3;
  l = -MAXL;
  mpz_mul_ui (c3k->coeff[3], f->coeff[d], (d * (d-1) * (d-2)) / 6);
  mpz_mul_ui (c3k->coeff[2], f->coeff[d-1], ((d-1) * (d-2)) / 2);
  mpz_mul_ui (c3k->coeff[1], f->coeff[d-2], d-2);
  mpz_addmul_si (c3k->coeff[1], g[1], l * (d-2));
  mpz_set (c3k->coeff[0], f->coeff[d-3]);
  mpz_addmul_si (c3k->coeff[0], g[0], l);
  mpz_init (r[0]);
  mpz_init (r[1]);
  mpz_init (r[2]);
  mpz_init (v);
  for (; l <= MAXL; l++)
    {
      nr = mpz_poly_mpz_roots (r, c3k);
      for (j = 0; j < nr; j++)
        {
          mpz_poly_eval (v, c3k, r[j]);
          insert_l (best_l2, best_v2, &n, l, fabs (mpz_get_d (v)));
        }
      mpz_add (c3k->coeff[0], c3k->coeff[0], g[0]);
      mpz_addmul_ui (c3k->coeff[1], g[1], d-2);
    }
  mpz_poly_clear (c3k);
  mpz_clear (r[0]);
  mpz_clear (r[1]);
  mpz_clear (r[2]);
  mpz_clear (v);
  return n;
}

/* assuming h(a)*h(b) < 0 where h(k) = h[n]*k^n+...+h[0] and a < b,
   refines the root such that a + 1 = b.
   The refined root is in a.
*/
static void
root_refine (mpz_t a, mpz_t b, mpz_poly_ptr h)
{
  mpz_t c, v;
  int sa, i;
  unsigned int n = h->deg;

  mpz_init (c);
  mpz_init (v);
  mpz_poly_eval (v, h, a);
  sa = mpz_sgn (v);
  while (1)
    {
      mpz_add (c, a, b);
      mpz_fdiv_q_2exp (c, c, 1);
      if (mpz_cmp (c, a) == 0)
        break;
      mpz_poly_eval (v, h, c);
      if (mpz_sgn (v) == sa)
        mpz_swap (a, c);
      else
        mpz_swap (b, c);
    }

  /* now b = a+1, we do one more iteration to round to nearest, i.e.,
     we check the sign of 2^n*h(a+1/2) = h[n]*(2a+1)^n + ... + 2^n*h[0] */

  mpz_mul_2exp (c, a, 1);
  mpz_add_ui (c, c, 1); /* c = 2a+1 */
  mpz_mul (v, h->coeff[n], c);
  mpz_addmul_ui (v, h->coeff[n-1], 2);
  for (i = n-2; i >= 0; i--)
    {
      mpz_mul (v, v, c);
      mpz_addmul_ui (v, h->coeff[i], 1UL << (n - i));
    }
  if (mpz_sgn (v) == sa) /* root is in [a+1/2,b], round to b */
    mpz_add_ui (a, a, 1);

  mpz_clear (c);
  mpz_clear (v);
}

/* put in r[0], r[1], r[2] integer approximations of the real roots of
   h[3]*x^3+...+h[0], and return the number of real roots */
static int
roots3 (mpz_t *r, mpz_poly_ptr H)
{
  int n = 0;
  mpz_t v, k;
  int s1, s2;
  mpz_t *h = H->coeff;

  ASSERT_ALWAYS(H->deg == 3);
  /* if h2^2-3*h1*h3 >= 0, the derivative of h has two real roots */
  mpz_mul (r[2], h[2], h[2]);
  mpz_mul (r[1], h[1], h[3]);
  mpz_mul_ui (r[1], r[1], 3);
  mpz_sub (r[2], r[2], r[1]);
  if (mpz_sgn (r[2]) >= 0) /* 1 or 3 real roots */
    {
      /* the roots of h' are (-h2 +/- sqrt(h2^2-3*h1*h3))/(3*h3) */
      mpz_sqrt (r[2], r[2]);
      mpz_neg (r[1], h[2]);
      mpz_sub (r[1], r[1], r[2]);
      mpz_add (r[2], r[2], h[2]);
      mpz_tdiv_q_ui (r[1], r[1], 3);
      mpz_tdiv_q (r[1], r[1], h[3]);
      mpz_tdiv_q_ui (r[2], r[2], 3);
      mpz_tdiv_q (r[2], r[2], h[3]);
      if (mpz_cmp (r[1], r[2]) > 0)
        mpz_swap (r[1], r[2]);
      /* now r[1] < r[2] are approximations of the two real roots of h' */
    }
  else /* only 1 real root */
    {
      mpz_set_ui (r[1], 0);
      mpz_set_ui (r[2], 0);
    }
  /* now we have four control points -Inf < r[1] < r[2] < +Inf */
  mpz_init (v);
  mpz_init (k);
  mpz_poly_eval (v, H, r[1]);
  s1 = mpz_sgn (v);
  if (-mpz_sgn (h[3]) * s1 < 0) /* one root in -Inf..r[1] */
    {
      mpz_set_si (k, -1);
      while (mpz_cmp (r[1], k) <= 0)
        mpz_mul_2exp (k, k, 1);
      while (1)
        {
          mpz_poly_eval (v, H, k);
          if (mpz_sgn (v) * s1 < 0)
            break;
          mpz_mul_2exp (k, k, 1);
        }
      root_refine (k, r[1], H);
      mpz_swap (r[n++], k);
    }
  mpz_poly_eval (v, H, r[2]);

  s2 = mpz_sgn (v);
  if (s1 * s2 < 0) /* one root in r[1]..r[2] */
    {
      root_refine (r[1], r[2], H);
      mpz_swap (r[n++], r[1]);
    }
  if (s2 * mpz_sgn (h[3]) < 0) /* one root in r[2]..+Inf */
    {
      mpz_set_ui (k, 1);
      while (mpz_cmp (k, r[2]) <= 0)
        mpz_mul_2exp (k, k, 1);
      while (1)
        {
          mpz_poly_eval (v, H, k);
          if (mpz_sgn (v) * s2 < 0)
            break;
          mpz_mul_2exp (k, k, 1);
        }
      root_refine (r[2], k, H);
      mpz_swap (r[n++], r[2]);
    }
  mpz_clear (v);
  mpz_clear (k);

  return n;
}

#if 0
/* put in r[0..n-1] the real roots of x^3 + a*x^2 + b*x + c and return n */
static int
cubic_roots (double *r, double a, double b, double c)
{
  double Q = (a * a - 3.0 * b) / 9.0;
  double R = ((2.0 * a * a - 9.0 * b) * a + 27.0 * c) / 54.0;
  double R2 = R * R;
  double Q3 = Q * Q * Q;
  double theta;

  a = a / 3.0;
  if (R2 < Q3)
    {
      double twopioverthree = 2.0943951023931954923;

      theta = acos (R / sqrt (Q3)) / 3.0;
      Q = -2.0 * sqrt (Q);
      r[0] = Q * cos (theta) - a;
      r[1] = Q * cos (theta + twopioverthree) - a;
      r[2] = Q * cos (theta - twopioverthree) - a;
      return 3;
    }
  else
    {
      double A = cbrt (fabs (R) + sqrt (R2 - Q3)), B;
      A = (R >= 0) ? -A : A;
      B = Q / A;
      r[0] = (A + B) - a;
      return 1;
    }
}

/* for degree 6, and a given l, estimate the smallest degree-3 coefficient
   of f + l*x^3*g */
static double
optimize_c3 (mpz_t *f, mpz_t *g, long l)
{
  double c6 = mpz_get_d (f[6]);
  double c5 = mpz_get_d (f[5]);
  double c4 = mpz_get_d (f[4]);
  double c3 = mpz_get_d (f[3]);
  double g1 = mpz_get_d (g[1]);
  double g0 = mpz_get_d (g[0]);
  double roots[3], k, v, vmin = DBL_MAX;
  int n;

  c6 = 20.0 * c6;
  c5 = 10.0 * c5;
  c4 = 4.0 * (c4 + (double) l * g1);
  c3 = c3 + (double) l * g0;
  n = cubic_roots (roots, c5 / c6, c4 / c6, c3 / c6);
  while (n)
    {
      k = round (roots[--n]);
      v = ((c6 * k + c5) * k + c4) * k + c3;
      v = fabs (v);
      if (v < vmin)
        vmin = v;
    }
  return vmin;
}
#endif

/* return 1 if new, 0 if already present in table */
int
insert_one (unsigned long *H, unsigned long alloc, unsigned long h)
{
  unsigned long i = h % alloc;

  while (H[i] != 0 && H[i] != h)
    if (++i == alloc)
      i = 0;
  if (H[i] == h) /* already present */
    return 0;
  else
  {
    ASSERT_ALWAYS(H[i] == 0);
    H[i] = h;
    return 1;
  }
}

int
insert_hash (unsigned long h)
{
  static unsigned long *H = NULL, size = 0, alloc = 0, inserted = 0;

  inserted++;

  if (h == 0) /* free table */
  {
    //printf ("inserted %lu elements, among which %lu new\n", inserted, size);
    free (H);
    H = NULL;
    size = alloc = inserted = 0;
    return 0;
  }

  if (2 * size >= alloc) /* realloc */
  {
    unsigned long *newH, newalloc, i, newsize = 0;

    newalloc = ulong_nextprime (2 * alloc + 1);
    newH = malloc (newalloc * sizeof (unsigned long));
    memset (newH, 0, newalloc * sizeof (unsigned long));
    for (i = 0; i < alloc; i++)
      if (H[i])
        newsize += insert_one (newH, newalloc, H[i]);
    free (H);
    H = newH;
    alloc = newalloc;
    ASSERT_ALWAYS(newsize == size);
  }

  if (insert_one (H, alloc, h))
  {
    size ++;
    return 1;
  }

  return 0;
}


static void
pre_reduce(mpz_poly_ptr f, mpz_t *g, mpz_t temp)
{
  for (int i = 0; i < f->deg - 3; i++) {
    mpz_ndiv_q (temp, f->coeff[i], g[0]);
    mpz_neg (temp, temp);
    rotate_auxg_z (f->coeff, g[1], g[0], temp, i);
  }
}

#define MAXQ 8

/* store in Q[0..MAXQ-1] the MAXQ best rational approximations of q
   with denominator <= bound. Q should have at least MAXQ+1 entries.
   Return the number of found approximations. */
int
make_rational (double *Q, double q, double bound)
{
  double i, p, e;
  int n = 0, j;
  double E[MAXQ+1];

  for (i = 2.0; i <= bound; i += 1.0)
    {
      p = floor (i * q + 0.5);
      e = fabs (q - p / i);
      /* search for duplicate before inserting */
      for (j = 0; j < n && e > E[j]; j++);
      if (j < n && e == E[j])
        continue;
      for (j = n; j > 0 && e < E[j-1]; j--)
	{
	  Q[j] = Q[j-1];
	  E[j] = E[j-1];
	}
      /* now j=0 or e > E[j-1] */
      Q[j] = p / i;
      E[j] = e;
      n += (n < MAXQ - 1); /* saturate at MAXQ */
    }
  /* force approximation with denominator 1 at the end */
  Q[n++] = floor (q + 0.5);
  ASSERT_ALWAYS(n <= MAXQ);
  return n;
}

/* Keep only those extrema which are closer to 0 than the points in their
   neighbourhood, i.e., where f(x) * f''(x) > 0. E.g., if f'(x) = 0 and
   f''(x) < 0, it's a local maximum, and if f(x) < 0, then this maximum is
   closer to 0 than its neighbourhood.
   The list of extrma is updated in-place, and the number of surviving
   extrema is returned. */
static size_t
keep_extrema_close_to_0 (double *extrema, const size_t n, double_poly_srcptr f)
{
  const int verbose = 0;
  if (f->deg < 2) {
    /* f''(x) = 0, so no point satisfies f(x)*f''(x) > 0 */
    return 0;
  }

  /* Compute f''(x) */
  double_poly_t f_deriv2;
  double_poly_init (f_deriv2, f->deg - 1);
  double_poly_derivative (f_deriv2, f);
  double_poly_derivative (f_deriv2, f_deriv2);

  size_t kept = 0;

  for (size_t i = 0; i < n; i++) {
    const double x = extrema[i];
    const double fx = double_poly_eval(f, x);
    const double fddx = double_poly_eval(f_deriv2, x);
    if (fx * fddx > 0.) {
      extrema[kept++] = x;
      if (verbose)
        printf ("Keeping x = %f, f(x) = %f, f''(x) = %f\n", x, fx, fddx);
    } else {
      if (verbose)
        printf ("Not keeping x = %f, f(x) = %f, f''(x) = %f\n", x, fx, fddx);
    }
  }
  double_poly_clear (f_deriv2);
  return kept;
}

/* Computes the roots of R, and the extrema of R which are closer to 0 than
   their neighbourhood. Stores those roots and extrema in roots_and_extrema,
   and returns their number. Note: roots_and_extrema must have enough storage
   for all roots and ALL extrema, whether close to 0 or not; i.e., 2*deg - 1
   entries is enough. */
static int
compute_roots_and_extrema_close_to_0(double *roots_and_extrema, double_poly_srcptr R)
{
  int n;
  n = double_poly_compute_all_roots (roots_and_extrema, R);

  /* add roots of the derivative */
  double_poly_t R_deriv;
  double_poly_init (R_deriv, R->deg - 1);
  double_poly_derivative (R_deriv, R);
  const int n_deriv = double_poly_compute_all_roots (roots_and_extrema + n, R_deriv);
  double_poly_clear (R_deriv);

  /* Keep only those extrema which are closer to 0 than points in their
     neighbourhood */
  n += keep_extrema_close_to_0 (roots_and_extrema + n, n_deriv, R);
  return n;
}

/* values of q greater than 1e10 in absolute value do not help */
#define MAX_Q 1e10

static int
find_best_k_deg5 (mpz_t *K, mpz_poly_srcptr f, mpz_t *g)
{
  int n, ret = 0;
  /* Let f(x) = a5 x^5 + ... + a0 and g(x) = g1 x + g0.
     Let f~(x,k,q) = f(x+k)+q*x^2*g(x+k), then the coefficient of x^2 in f~ is
     c_2(k,q) = 10*k^3*a5 + k*q*g1 + 6*k^2*a4 + q*g0 + 3*k*a3 + a2
     and the coefficient of x^3 in f~ is
     c_3(k,q) = 10*k^2*a5 + q*g1 + 4*k*a4 + a3
     We want to choose q such that these two polynomials have a common root,
     such that the resultant is 0.

     Sage code:
     R.<x,k,q,g0,g1,a0,a1,a2,a3,a4,a5> = ZZ[]
     f = a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0
     g = g1*x+g0
     ff = f(x=x+k)+q*x^2*g(x=x+k)
     c2 = ff.coefficient({x:2})
     c3 = ff.coefficient({x:3})
     r = c2.resultant(c3,k)//(40*a5)
     SR(r.coefficient({q:2})).factor()
  */
  double_poly_t R;
  double *r, roots_q[5], roots_k[5];
  double_poly_init (R, 3);
  r = R->coeff;

  const double a5 = mpz_get_d (f->coeff[5]),
               a4 = mpz_get_d (f->coeff[4]),
               a3 = mpz_get_d (f->coeff[3]),
               a2 = mpz_get_d (f->coeff[2]),
               g1 = mpz_get_d (g[1]),
               g0 = mpz_get_d (g[0]);

  R->deg = 2;
  /* r.coefficient({q:2}) */
  r[2] = g1*g1*a4*a4 - 10.*g0*g1*a4*a5 + 25.*g0*g0*a5*a5;

  /* r.coefficient({q:1}) */
  r[1] = -2.*g1*a3*a4*a4 + 8.*g0*a4*a4*a4 + 10.*g1*a3*a3*a5
         - 10.*g1*a2*a4*a5 - 30.*g0*a3*a4*a5 + 50.*g0*a2*a5*a5;

  /* r.coefficient({q:0}) */
  r[0] = -3.*a3*a3*a4*a4 + 8.*a2*a4*a4*a4 + 10.*a3*a3*a3*a5 - 30.*a2*a3*a4*a5
         + 25.*a2*a2*a5*a5;

  n = compute_roots_and_extrema_close_to_0(roots_q, R);

  for (int i = 0; i < n; i++)
    {
      double Q[MAXQ+1];
      int m, t;
      /* it looks like values of q larger than 1 in absolute value do
         not give any good translation */
      if (fabs (roots_q[i]) > MAX_Q)
        continue;
      t = make_rational (Q, roots_q[i], 100.0);
      while (t) {
        const double q = roots_q[i] = Q[--t];
        /* find roots k of c3(k): 10*k^2*a5 + q*g1 + 4*k*a4 + a3 */
        R->deg = 2;
        R->coeff[0] = q*g1 + a3;
        R->coeff[1] = 4*a4;
        R->coeff[2] = 10*a5;
        m = double_poly_compute_all_roots (roots_k, R);
        ASSERT_ALWAYS(m <= 2);

        /* find roots k of c2(k): 10*k^3*a5 + k*q*g1 + 6*k^2*a4 + q*g0 + 3*k*a3 + a2 */

        R->deg = 3;
        R->coeff[0] = q*g0 + a2;
        R->coeff[1] = q*g1 + 3*a3;
        R->coeff[2] = 6*a4;
        R->coeff[3] = 10*a5;
        m += double_poly_compute_all_roots (roots_k + m, R);
        ASSERT_ALWAYS(m <= 5);

        for (int j = 0; j < m; j++, ret++)
          mpz_set_d (K[ret], roots_k[j] > 0 ? roots_k[j] + 0.5
                     : roots_k[j] - 0.5);
      }
    }
  double_poly_clear (R);

  return ret;
}

static int
find_best_k_deg6 (mpz_t *K, mpz_poly_ptr f, mpz_t *g)
{
  int n, i, j, t, m, ret = 0;
  /* Let f(x) = a6 x^6 + ... + a0 and g(x) = g1 x + g0.
     Let f~(x,k,q) = f(x+k)+q*x^3*g(x+k), then the coefficient of x^3 in f~ is
     c_3(k,q) = 20 a6 k^3 + 10 a5 k^2 + 4 a4 k + q g1 k + a3 - q g0
     We want to choose q such that this polynomial in k has a root
     of multiplicity 2, i.e., such that the discriminant is 0.
     The discriminant of c_3(k,q) w.r.t. k is a polynomial of degree 3 in q,
     we call it r(q) below. We need the roots of this polynomial.

     Sage code:
     R.<x,a6,k,a5,g1,q,g0,a3,a4> = ZZ[]
     f = a6*x^6+a5*x^5+a4*x^4+a3*x^3
     g = g1*x+g0
     ff=f(x=x+k)+q*x^3*g(x=x+k)
     c3=ff.coefficient({x:3})
     c4=ff.coefficient({x:4})
     r = c3.resultant(c4,k)//(25*a6)
     SR(r.coefficient({q:3})).factor() */
  double_poly_t R, C;
#define MAX_ROOTSK 40
  double roots_k[MAX_ROOTSK];
  double *r, a6, a5, a4, a3, a2, g1, g0, roots_q[5];
  double_poly_init (R, 3);
  R->deg = 3;
  r = R->coeff;

  a6 = mpz_get_d (f->coeff[6]);
  a5 = mpz_get_d (f->coeff[5]);
  a4 = mpz_get_d (f->coeff[4]);
  a3 = mpz_get_d (f->coeff[3]);
  a2 = mpz_get_d (f->coeff[2]);
  g1 = mpz_get_d (g[1]);
  g0 = mpz_get_d (g[0]);

  /* degree 3: a6*g1^3 */
  r[3] = a6 * g1 * g1 * g1;

 /* 2: 135*a6^2*g0^2 - 45*a5*a6*g0*g1 + 10*a5^2*g1^2 - 15*a4*a6*g1^2 */
  r[2] = 135.0 * a6 * a6 * g0 * g0 - 45.0 * a5 * a6 * g0 * g1
    + 10.0 * a5 * a5 * g1 * g1 - 15.0 * a4 * a6 * g1 * g1;

  /* 1: 50*a5^3*g0 - 180*a4*a5*a6*g0 + 270*a3*a6^2*g0 - 10*a4*a5^2*g1
     + 48*a4^2*a6*g1 - 45*a3*a5*a6*g1 */
  r[1] = 50.0 * a5 * a5 * a5 * g0 - 180.0 * a4 * a5 * a6 * g0
    + 270.0 * a3 * a6 * a6 * g0 - 10.0 * a4 * a5 * a5 * g1
    + 48.0 * a4 * a4 * a6 * g1 - 45.0 * a3 * a5 * a6 * g1;

  /* 0: -20*a4^2*a5^2 + 50*a3*a5^3 + 64*a4^3*a6 - 180*a3*a4*a5*a6
     + 135*a3^2*a6^2 */
  r[0] = -20.0 * a4 * a4 * a5 * a5 + 50.0 * a3 * a5 * a5 * a5
    + 64.0 * a4 * a4 * a4 * a6 - 180.0 * a3 * a4 * a5 * a6
    + 135.0 * a3 * a3 * a6 * a6;

  n = double_poly_compute_all_roots_with_bound (roots_q, R, MAX_Q);

  /* add roots of the derivative */
  r[0] = r[1];
  r[1] = r[2] * 2.0;
  r[2] = r[3] * 3.0;
  R->deg = 2;

  n += double_poly_compute_all_roots_with_bound (roots_q + n, R, MAX_Q);

  double_poly_init (C, 4);
  for (i = 0; i < n; i++)
    {
      double Q[MAXQ+1];
      t = make_rational (Q, roots_q[i], 100.0);
      while (t) {
      roots_q[i] = Q[--t];
      /* find roots k of c4(k): 15*a6*k^2 + 5*k*a5 + g1*q + a4 = 0 */
      C->deg = 2;
      C->coeff[0] = g1 * roots_q[i] + a4;
      C->coeff[1] = 5.0 * a5;
      C->coeff[2] = 15.0 * a6;

      m = double_poly_compute_all_roots (roots_k, C);
      ASSERT_ALWAYS(m <= 2);

      /* find roots k of c3(k):
         20*a6*k^3 + 10*k^2*a5 + k*g1*q + q*g0 + 4*k*a4 + a3 */
      C->deg = 3;
      C->coeff[3] = 20.0 * a6;
      C->coeff[2] = 10.0 * a5;
      C->coeff[1] = g1 * roots_q[i] + 4.0 * a4;
      C->coeff[0] = roots_q[i] * g0 + a3;

      m += double_poly_compute_all_roots (roots_k + m, C);
      ASSERT_ALWAYS(m <= 5);

      /* add roots of Res(c3,c2) with 2*x^2*g(x+k):
         R.<a6,k,a5,g1,q,q2,g0,a3,a4,a2,x> = ZZ[]
         f = a6*x^6+a5*x^5+a4*x^4+a3*x^3+a2*x^2
         g = g1*x+g0
         ff = f(x=x+k)+q*x^3*g(x=x+k)+q2*x^2*g(x=x+k)
         c3 = ff.coefficient({x:3})
         c2 = ff.coefficient({x:2})
         c3.resultant(c2,q2) */
      C->deg = 4;
      C->coeff[4] = -5.0 * a6 * g1;
      C->coeff[3] =  -20.0 * a6 * g0;
      C->coeff[2] =  -g1 * g1 * roots_q[i] - 10.0 * a5 * g0 + 2.0 * g1 * a4;
      C->coeff[1] =  2.0 * g1 * a3 - 2.0 * g1 * roots_q[i] * g0 - 4.0 * g0 * a4;
      C->coeff[0] = -roots_q[i] * g0 * g0 - g0 * a3 + g1 * a2;

      m += double_poly_compute_all_roots (roots_k + m, C);
      ASSERT_ALWAYS(m <= 9);

      /* roots of the derivative of Res(c3,c2) */
      C->deg = 3;
      C->coeff[0] = C->coeff[1];
      C->coeff[1] = C->coeff[2] * 2.0;
      C->coeff[2] = C->coeff[3] * 3.0;
      C->coeff[3] = C->coeff[4] * 4.0;

      m += double_poly_compute_all_roots (roots_k + m, C);
      ASSERT_ALWAYS(m <= 12);


      /* add roots of Res(c3(q3=roots_q[i](k),c2(k)) (eliminate k) */
      double m0 = -g0;
      double roots_q2[7];
      int ii, nb_roots_q2 = 0;

      C->deg = 4;
      C->coeff[4] = -125*a6*a6*a6*g1*g1*g1*g1;
      C->coeff[3] = (3000*a6*a6*a6*m0*g1*g1*g1 + 500*a5*a6*a6*g1*g1*g1*g1)*roots_q[i] - 160000*a6*a6*a6*a6*m0*m0*m0 - 80000*a5*a6*a6*a6*m0*m0*g1 - 5000*a5*a5*a6*a6*m0*g1*g1 - 20000*a4*a6*a6*a6*m0*g1*g1 - 1000*a4*a5*a6*a6*g1*g1*g1 - 3500*a3*a6*a6*a6*g1*g1*g1;
      C->coeff[2] = -225*a6*a6*g1*g1*g1*g1*g1*roots_q[i]*roots_q[i]*roots_q[i] + (15750*a6*a6*a6*m0*m0*g1*g1 + 5250*a5*a6*a6*m0*g1*g1*g1 - 500*a5*a5*a6*g1*g1*g1*g1 + 2250*a4*a6*a6*g1*g1*g1*g1)*roots_q[i]*roots_q[i] + (20000*a5*a5*a6*a6*m0*m0*g1 - 48000*a4*a6*a6*a6*m0*m0*g1 + 10000*a5*a5*a5*a6*m0*g1*g1 - 28000*a4*a5*a6*a6*m0*g1*g1 + 18000*a3*a6*a6*a6*m0*g1*g1 + 2000*a4*a5*a5*a6*g1*g1*g1 - 7200*a4*a4*a6*a6*g1*g1*g1 + 5000*a3*a5*a6*a6*g1*g1*g1 - 4000*a2*a6*a6*a6*g1*g1*g1)*roots_q[i] - 50000*a5*a5*a5*a5*a6*m0*m0 + 240000*a4*a5*a5*a6*a6*m0*m0 - 192000*a4*a4*a6*a6*a6*m0*m0 - 240000*a3*a5*a6*a6*a6*m0*m0 + 480000*a2*a6*a6*a6*a6*m0*m0 - 20000*a4*a5*a5*a5*a6*m0*g1 + 80000*a4*a4*a5*a6*a6*m0*g1 + 10000*a3*a5*a5*a6*a6*m0*g1 - 216000*a3*a4*a6*a6*a6*m0*g1 + 160000*a2*a5*a6*a6*a6*m0*g1 - 2000*a4*a4*a5*a5*a6*g1*g1 + 7200*a4*a4*a4*a6*a6*g1*g1 - 1000*a3*a4*a5*a6*a6*g1*g1 + 5000*a2*a5*a5*a6*a6*g1*g1 - 33750*a3*a3*a6*a6*a6*g1*g1 + 20000*a2*a4*a6*a6*a6*g1*g1;
      C->coeff[1] = (-18000*a6*a6*a6*m0*m0*m0*g1 - 9000*a5*a6*a6*m0*m0*g1*g1 - 3600*a4*a6*a6*m0*g1*g1*g1 - 900*a3*a6*a6*g1*g1*g1*g1)*roots_q[i]*roots_q[i]*roots_q[i] + (-15000*a5*a5*a6*a6*m0*m0*m0 + 36000*a4*a6*a6*a6*m0*m0*m0 - 20000*a5*a5*a5*a6*m0*m0*g1 + 63000*a4*a5*a6*a6*m0*m0*g1 - 67500*a3*a6*a6*a6*m0*m0*g1 - 8000*a4*a5*a5*a6*m0*g1*g1 + 25200*a4*a4*a6*a6*m0*g1*g1 - 7500*a3*a5*a6*a6*m0*g1*g1 - 30000*a2*a6*a6*a6*m0*g1*g1 - 2000*a3*a5*a5*a6*g1*g1*g1 + 6300*a3*a4*a6*a6*g1*g1*g1 - 5000*a2*a5*a6*a6*g1*g1*g1)*roots_q[i]*roots_q[i] + (8000*a4*a4*a5*a5*a6*m0*g1 - 20000*a3*a5*a5*a5*a6*m0*g1 - 28800*a4*a4*a4*a6*a6*m0*g1 + 84000*a3*a4*a5*a6*a6*m0*g1 - 20000*a2*a5*a5*a6*a6*m0*g1 - 81000*a3*a3*a6*a6*a6*m0*g1 + 48000*a2*a4*a6*a6*a6*m0*g1 + 2000*a3*a4*a5*a5*a6*g1*g1 - 10000*a2*a5*a5*a5*a6*g1*g1 - 7200*a3*a4*a4*a6*a6*g1*g1 + 4500*a3*a3*a5*a6*a6*g1*g1 + 32000*a2*a4*a5*a6*a6*g1*g1 - 36000*a2*a3*a6*a6*a6*g1*g1)*roots_q[i] + 16000*a4*a4*a4*a5*a5*a6*m0- 60000*a3*a4*a5*a5*a5*a6*m0+ 100000*a2*a5*a5*a5*a5*a6*m0- 57600*a4*a4*a4*a4*a6*a6*m0+ 240000*a3*a4*a4*a5*a6*a6*m0+ 15000*a3*a3*a5*a5*a6*a6*m0- 480000*a2*a4*a5*a5*a6*a6*m0- 324000*a3*a3*a4*a6*a6*a6*m0+ 384000*a2*a4*a4*a6*a6*a6*m0+ 480000*a2*a3*a5*a6*a6*a6*m0- 480000*a2*a2*a6*a6*a6*a6*m0+ 4000*a3*a4*a4*a5*a5*a6*g1 - 20000*a3*a3*a5*a5*a5*a6*g1 + 20000*a2*a4*a5*a5*a5*a6*g1 - 14400*a3*a4*a4*a4*a6*a6*g1 + 81000*a3*a3*a4*a5*a6*a6*g1 - 80000*a2*a4*a4*a5*a6*a6*g1 - 10000*a2*a3*a5*a5*a6*a6*g1 - 121500*a3*a3*a3*a6*a6*a6*g1 + 216000*a2*a3*a4*a6*a6*a6*g1 - 80000*a2*a2*a5*a6*a6*a6*g1;
      C->coeff[0] = (3375*a6*a6*a6*m0*m0*m0*m0 + 2250*a5*a6*a6*m0*m0*m0*g1 + 1350*a4*a6*a6*m0*m0*g1*g1 + 675*a3*a6*a6*m0*g1*g1*g1 + 225*a2*a6*a6*g1*g1*g1*g1)*roots_q[i]*roots_q[i]*roots_q[i]*roots_q[i] + (5000*a5*a5*a5*a6*m0*m0*m0 - 18000*a4*a5*a6*a6*m0*m0*m0 + 27000*a3*a6*a6*a6*m0*m0*m0 + 3000*a4*a5*a5*a6*m0*m0*g1 - 10800*a4*a4*a6*a6*m0*m0*g1 + 4500*a3*a5*a6*a6*m0*m0*g1 + 18000*a2*a6*a6*a6*m0*m0*g1 + 1500*a3*a5*a5*a6*m0*g1*g1 - 5400*a3*a4*a6*a6*m0*g1*g1 + 6000*a2*a5*a6*a6*m0*g1*g1 + 500*a2*a5*a5*a6*g1*g1*g1 - 675*a3*a3*a6*a6*g1*g1*g1)*roots_q[i]*roots_q[i]*roots_q[i] + (-6000*a4*a4*a5*a5*a6*m0*m0 + 15000*a3*a5*a5*a5*a6*m0*m0 + 21600*a4*a4*a4*a6*a6*m0*m0 - 63000*a3*a4*a5*a6*a6*m0*m0 + 15000*a2*a5*a5*a6*a6*m0*m0 + 60750*a3*a3*a6*a6*a6*m0*m0 - 36000*a2*a4*a6*a6*a6*m0*m0 - 3000*a3*a4*a5*a5*a6*m0*g1 + 15000*a2*a5*a5*a5*a6*m0*g1 + 10800*a3*a4*a4*a6*a6*m0*g1 - 6750*a3*a3*a5*a6*a6*m0*g1 - 48000*a2*a4*a5*a6*a6*m0*g1 + 54000*a2*a3*a6*a6*a6*m0*g1 - 1500*a3*a3*a5*a5*a6*g1*g1 + 3000*a2*a4*a5*a5*a6*g1*g1 + 4050*a3*a3*a4*a6*a6*g1*g1 - 7200*a2*a4*a4*a6*a6*g1*g1 - 3000*a2*a3*a5*a6*a6*g1*g1 + 12000*a2*a2*a6*a6*a6*g1*g1)*roots_q[i]*roots_q[i] + 6000*a3*a3*a4*a4*a5*a5*a6 - 16000*a2*a4*a4*a4*a5*a5*a6 - 20000*a3*a3*a3*a5*a5*a5*a6 + 60000*a2*a3*a4*a5*a5*a5*a6 - 50000*a2*a2*a5*a5*a5*a5*a6 - 21600*a3*a3*a4*a4*a4*a6*a6 + 57600*a2*a4*a4*a4*a4*a6*a6 + 81000*a3*a3*a3*a4*a5*a6*a6 - 240000*a2*a3*a4*a4*a5*a6*a6 - 15000*a2*a3*a3*a5*a5*a6*a6 + 240000*a2*a2*a4*a5*a5*a6*a6 - 91125*a3*a3*a3*a3*a6*a6*a6 + 324000*a2*a3*a3*a4*a6*a6*a6 - 192000*a2*a2*a4*a4*a6*a6*a6 - 240000*a2*a2*a3*a5*a6*a6*a6 + 160000*a2*a2*a2*a6*a6*a6*a6;


      nb_roots_q2 += double_poly_compute_all_roots_with_bound (roots_q2, C,
                                                               MAX_Q);
      ASSERT_ALWAYS(nb_roots_q2 <= 4);

      /* roots of the derivative of Res(c3,c2) */
      C->deg = 3;
      C->coeff[0] = C->coeff[1];
      C->coeff[1] = C->coeff[2] * 2.0;
      C->coeff[2] = C->coeff[3] * 3.0;
      C->coeff[3] = C->coeff[4] * 4.0;
      nb_roots_q2 += double_poly_compute_all_roots_with_bound (roots_q2
                                                  + nb_roots_q2, C, MAX_Q);
      ASSERT_ALWAYS(nb_roots_q2 <= 7);

      for (ii = 0; ii < nb_roots_q2; ii++)
      {
        /* Compute roots of c2(q2=roots_q2[ii]) */
        C->deg = 4;
        C->coeff[4] = 15.0 * a6;
        C->coeff[3] = 10.0 * a5;
        C->coeff[2] = 6.0 * a4;
        C->coeff[1] = g1 * roots_q2[ii] + 3.0 * a3;
        C->coeff[0] = a2 + g0 * roots_q2[ii];

        m += double_poly_compute_all_roots (roots_k + m, C);
        ASSERT_ALWAYS(nb_roots_q2 <= 12+4*(ii+1));
      }

#define MAX_ROOTS (MAXQ*5*MAX_ROOTSK+1) /* q-roots: 3 for Delta and 2 for its derivative
                       k-roots: 2 for c4, 3 for c3, 4 for Res(c3,c2),
                       3 for its derivative, +1 for k=0 */
      for (j = 0; j < m; j++, ret++)
        mpz_set_d (K[ret], roots_k[j] > 0 ? roots_k[j] + 0.5
                   : roots_k[j] - 0.5);
      }
    }
  double_poly_clear (C);

  ASSERT_ALWAYS(ret < MAX_ROOTS);

  double_poly_clear (R);
  return ret;
}

/* put in K[0], K[1], ... good values of k such that f(x+k) has small
   coefficients of degree d-2 and d-3.
   For degree 6, we should have space for at least 6 coefficients in K. */
int
find_best_k (mpz_t *K, mpz_poly_ptr f, mpz_t *g)
{
  int d = f->deg, i, j, ret = 0;

  if (d == 6)
    {
      ret = find_best_k_deg6 (K, f, g);
    }
  else if (d == 5)
    {
      ret = find_best_k_deg5 (K, f, g);
    }
  else
    ASSERT_ALWAYS(0);

  mpz_set_ui (K[ret++], 0); /* add k=0 */

  /* remove duplicates */
  int uniq = 0;
  for (i = 0; i < ret; i++)
    {
      for (j = 0; j < i; j++)
	if (mpz_cmp (K[j], K[i]) == 0)
	  break;
      if (j == i)
	mpz_set (K[uniq++], K[i]);
    }

  return uniq;
}

#ifdef OPTIMIZE_LLL_LIST

void new_polylll_pq (polylll_pq **pqueue, int len, int d) {

  int i, j;

  if (len < 2) {
    fprintf(stderr,"Error: len too short.\n");
    exit(1);
  }

  (*pqueue) = (polylll_pq *) malloc (sizeof (polylll_pq));
  if ((*pqueue) == NULL) {
    fprintf(stderr,"Error: malloc failed.\n");
    exit(1);
  }

  /* len */
  (*pqueue)->len = len;
  (*pqueue)->degree = d;

  /* f */
  (*pqueue)->f = (mpz_t **) malloc (len*sizeof(mpz_t*));
  if ((*pqueue)->f != NULL) {
    for (i = 0; i < (*pqueue)->len; i++) {
      (*pqueue)->f[i] = (mpz_t *) malloc ((d+1)*sizeof(mpz_t));
      for (j = 0; j <= d; j++) {
        mpz_init ((*pqueue)->f[i][j]);
      }
    }
  }

  /* g */
  (*pqueue)->g = (mpz_t **) malloc (len*sizeof(mpz_t*));
  if ((*pqueue)->g != NULL) {
    for (i = 0; i < (*pqueue)->len; i++) {
      (*pqueue)->g[i] = (mpz_t *) malloc (2*sizeof(mpz_t));
      mpz_init ((*pqueue)->g[i][0]);
      mpz_init ((*pqueue)->g[i][1]);
    }
  }

  /* others */
  (*pqueue)->fd = (mpz_t *) malloc (len*sizeof(mpz_t));
  (*pqueue)->lognorm = (double *) malloc (len*sizeof(double));
  (*pqueue)->alpha = (double *) malloc (len*sizeof(double));
  (*pqueue)->skew = (double *) malloc (len*sizeof(double));
  if ((*pqueue)->lognorm != NULL && (*pqueue)->alpha != NULL &&
      (*pqueue)->skew != NULL && (*pqueue)->fd != NULL) {
    for (i = 0; i < (*pqueue)->len; i++) {
      mpz_init_set_ui ((*pqueue)->fd[i], 0);
      (*pqueue)->lognorm[i] = 0;
      (*pqueue)->alpha[i] = 0;
      (*pqueue)->skew[i] = 0;
    }
  }

  (*pqueue)->lognorm[0] = DBL_MAX;
  (*pqueue)->used = 1;
}


static inline int pq_parent ( int i ) {return (i>>1);}


static inline void
insert_polylll_pq_up (polylll_pq *pqueue, mpz_t *f, mpz_t *g,
                      double lognorm, double alpha, double skew)
{
  int i, j, d;
  d = pqueue->degree;
  for (i=pqueue->used; lognorm>pqueue->lognorm[pq_parent(i)]; i/=2) {
    /* only record polynomials with different ad */
    if (mpz_cmp(f[d], pqueue->f[pq_parent(i)][d])==0)
      return;
    mpz_set (pqueue->g[i][0], pqueue->g[pq_parent(i)][0]);
    mpz_set (pqueue->g[i][1], pqueue->g[pq_parent(i)][1]);
    for (j=0;j<=pqueue->degree;j++)
      mpz_set (pqueue->f[i][j], pqueue->f[pq_parent(i)][j]);
    pqueue->lognorm[i] = pqueue->lognorm[pq_parent(i)];
    pqueue->alpha[i] = pqueue->alpha[pq_parent(i)];
    pqueue->skew[i] = pqueue->skew[pq_parent(i)];
  }
  mpz_set (pqueue->g[i][0], g[0]);
  mpz_set (pqueue->g[i][1], g[1]);
  for (j=0;j<=pqueue->degree;j++)
    mpz_set (pqueue->f[i][j], f[j]);
  pqueue->lognorm[i] = lognorm;
  pqueue->alpha[i] = alpha;
  pqueue->skew[i] = skew;
  pqueue->used ++;
}


static inline void
insert_polylll_pq_down (polylll_pq *pqueue, mpz_t *f, mpz_t *g,
                        double lognorm, double alpha, double skew)
{
  int i, j, l, d;
  d = pqueue->degree;

  for (i = 1; i*2 < pqueue->used; i = l) {

    l = (i << 1);

    if (mpz_cmp(f[d], pqueue->f[i][d])==0 ||
        mpz_cmp(f[d], pqueue->f[l][d])==0)
      return;

    if ((l+1) < pqueue->used && pqueue->lognorm[l+1]
        > pqueue->lognorm[l])
      l ++;

    /* only record polynomials with different ad */
    if (mpz_cmp(f[d], pqueue->f[l][d])==0)
      return;

    if(pqueue->lognorm[l] > lognorm) {
      mpz_set (pqueue->g[i][1], pqueue->g[l][1]);
      mpz_set (pqueue->g[i][0], pqueue->g[l][0]);
      for (j=0;j<=pqueue->degree;j++)
        mpz_set (pqueue->f[i][j], pqueue->f[l][j]);
      pqueue->lognorm[i] = pqueue->lognorm[l];
      pqueue->alpha[i] = pqueue->alpha[l];
      pqueue->skew[i] = pqueue->skew[l];
    }
    else
      break;
  }
  mpz_set (pqueue->g[i][0], g[0]);
  mpz_set (pqueue->g[i][1], g[1]);
  for (j=0;j<=pqueue->degree;j++)
    mpz_set (pqueue->f[i][j], f[j]);
  pqueue->lognorm[i] = lognorm;
  pqueue->alpha[i] = alpha;
  pqueue->skew[i] = skew;
}


void insert_polylll_pq (polylll_pq *pqueue, mpz_t *f, mpz_t *g,
                        double lognorm, double alpha, double skew)
{

  /* only allow non-repeated cd */
  int i, d;
  d = pqueue->degree;
  for (i=0; i<pqueue->used; i++)
    if (mpz_cmp (pqueue->f[i][d], f[d]) == 0)
      return;
  /**/

  if (pqueue->len == pqueue->used) {
    if (lognorm < pqueue->lognorm[1]) {
      insert_polylll_pq_down (pqueue, f, g, lognorm, alpha, skew);
    }
  }
  else if (pqueue->len > pqueue->used) {
    insert_polylll_pq_up (pqueue, f, g, lognorm, alpha, skew);
  }
  else {
    fprintf(stderr,"Error: error (pqueue->len < pqueue->used) "
            "in insert_sublattice_pq()\n");
    exit(1);
  }
}


void extract_polylll_pq (polylll_pq *pqueue, mpz_t *f, mpz_t *g,
                         double *lognorm)
{
  /* don't extract u[0] since it is just a placeholder. */
  pqueue->used --;
  mpz_set (g[0], pqueue->g[1][0]);
  mpz_set (g[1], pqueue->g[1][1]);
  *lognorm = pqueue->lognorm[1];
  for (int j=0;j<=pqueue->degree;j++)
    mpz_set (f[j], pqueue->f[1][j]);

  insert_polylll_pq_down (pqueue,
                          pqueue->f[pqueue->used],
                          pqueue->g[pqueue->used],
                          pqueue->lognorm[pqueue->used],
                          pqueue->alpha[pqueue->used],
                          pqueue->skew[pqueue->used]);
}


void free_polylll_pq (polylll_pq **pqueue)
{
  int i, j;
  for (i = 0; i < (*pqueue)->len; i++) {
    mpz_clear ((*pqueue)->g[i][0]);
    mpz_clear ((*pqueue)->g[i][1]);
    mpz_clear ((*pqueue)->fd[i]);
    for (j = 0; j <= (*pqueue)->degree; j++) {
      mpz_clear ((*pqueue)->f[i][j]);
    }
    free ((*pqueue)->g[i]);
    free ((*pqueue)->f[i]);
  }
  free ((*pqueue)->f);
  free ((*pqueue)->g);
  free ((*pqueue)->lognorm);
  free ((*pqueue)->fd);
  free ((*pqueue)->alpha);
  free ((*pqueue)->skew);
  free (*pqueue);
}

void
optimize_lll_list (mpz_poly_ptr f, mpz_t *g, polylll_pq *pqueue, double lognorm0, int verbose)
{
  double skew, skew2, best_lognorm, lognorm, smax;
  mpz_t s, det, a, b;
  int d = f->deg, i, j, ss;
  mat_Z m;
  mpz_poly_t best_f, copy_f;
  mpz_t best_g0, copy_g0, k, copy_g0k;

  mpz_init (s);
  mpz_init (det);

  /* 1/4 < a/b < 1: the closer a/b is from 1, the better the reduction is.
     Some experiments suggest that increasing a/b does not yield better
     polynomials, thus we take the smallest possible value. */
  mpz_init_set_ui (a, 1);
  mpz_init_set_ui (b, 4);
  mpz_poly_init (best_f, d);
  mpz_poly_init (copy_f, d);
  mpz_poly_set (copy_f, f);
  mpz_init_set (copy_g0, g[0]);
  mpz_poly_set (best_f, f);
  mpz_init_set (best_g0, g[0]);
  best_lognorm = L2_skew_lognorm (f, SKEWNESS_DEFAULT_PREC);
  mpz_init (k);
  mpz_init (copy_g0k);

  /* m has d-1 (row) vectors of d+1 coefficients each */
  m.NumRows = d - 1;
  m.NumCols = d + 1;
  m.coeff = (mpz_t **) malloc (d * sizeof(mpz_t*));
  for (j = 1; j <= d-1; j++) {
    m.coeff[j] = (mpz_t *) malloc ((d + 2) * sizeof(mpz_t));
    for (i = 0; i <= d + 1; i++)
      mpz_init (m.coeff[j][i]);
  }
  smax = pow (fabs (mpz_get_d (g[0]) / mpz_get_d (f->coeff[d])),
              1.0 / (double) d);

#define K ((16*optimize_effort)/10)
#define MULT (1000*optimize_effort)
//#define NSKEW (1+(5*optimize_effort)/10)
#define NSKEW (1+(10*optimize_effort))

  /* for each translation k */
  for (mpz_set_si (k, (-K) * MULT); mpz_cmp_ui (k, K * MULT) <= 0; mpz_add_ui (k, k, MULT)) {

        /* for each skewness */
    for (ss = NSKEW + 1; ss <= 2 * NSKEW ; ss++) {

      skew = pow (smax, (double) ss / (double) (2*NSKEW));

      mpz_poly_set (f, copy_f);
      mpz_set (g[0], copy_g0);
      do_translate_z (f, g, k);
      mpz_set (copy_g0k, g[0]);
      mpz_set_d (s, skew + 0.5); /* round to nearest */

      //gmp_fprintf (stdout, "# ==> skew is %Zd\n", s);
      //print_poly_fg (f, g, g[0], 1);

      for (j = 1; j <= d-1; j++)
        for (i = 0; i <= d; i++)
          mpz_set_ui (m.coeff[j][i+1], 0);
      mpz_set_ui (det, 1);

      for (i = 0; i <= d; i++) {
        if (i > 0)
          mpz_mul (det, det, s); /* s^i */
        if (i <= d - 3)
          mpz_set (m.coeff[i+2][i+1], det);
        mpz_mul (m.coeff[1][i+1], det, f->coeff[i]);
      }

      for (i = 0; i <= d - 3; i++) {
        /* m.coeff[i+2][i+1] is already s^i */
        mpz_mul (m.coeff[i+2][i+2], m.coeff[i+2][i+1], s); /* s^(i+1) */
        mpz_mul (m.coeff[i+2][i+1], m.coeff[i+2][i+1], g[0]);
        mpz_mul (m.coeff[i+2][i+2], m.coeff[i+2][i+2], g[1]);
      }

      /*
      for (j = 1; j <= d-1; j++) {
        for (i = 1; i <= d+1; i++) {
          gmp_fprintf (stdout, " %Zd ", m.coeff[j][i]);
        }
      fprintf (stdout, "\n");
      }
      */
      LLL (det, m, NULL, a, b);

      /*
      fprintf (stdout, "--- AFTER --\n");
      for (j = 1; j <= d-1; j++) {
        for (i = 1; i <= d+1; i++) {
          gmp_fprintf (stdout, " %Zd ", m.coeff[j][i]);
        }
      fprintf (stdout, "\n");
      }
      */

      /* for each row */
      for (j = 1; j <= d-1; j++) {

        if (mpz_sgn (m.coeff[j][d]) != 0) {
          mpz_set_ui (det, 1);

          for (i = 0; i <= d; i++) {
            if (i > 0)
              mpz_mul (det, det, s); /* det = s^i */
            ASSERT_ALWAYS(mpz_divisible_p (m.coeff[j][i+1], det));
            mpz_divexact (f->coeff[i], m.coeff[j][i+1], det);
          }

          if (insert_hash (mpz_getlimbn (f->coeff[0], 0))) {
            mpz_set (g[0], copy_g0k);
            optimize_aux (f, g, verbose, 1, OPT_STEPS);

            skew2 = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
            lognorm = L2_lognorm (f, skew2);

            if (lognorm <= lognorm0 + 1) {
              insert_polylll_pq (pqueue, f->coeff, g, lognorm, 0.0, skew2);
              //gmp_fprintf (stdout, "c6: %Zd\n", f->coeff[6]);
              //gmp_fprintf (stdout, " lognorm: %f\n", lognorm);
            }

            if (lognorm < best_lognorm) {
              best_lognorm = lognorm;
              mpz_poly_set (best_f, f);
              mpz_set (best_g0, g[0]);
            }
          }
        }
      }
    }
  }
  for (j = 1; j <= d-1; j++)
  {
    for (i = 0; i <= d + 1; i++)
      mpz_clear (m.coeff[j][i]);
    free (m.coeff[j]);
  }
  free (m.coeff);
  mpz_poly_set (f, best_f);
  mpz_set (g[0], best_g0);
  mpz_clear (s);
  mpz_clear (det);
  mpz_clear (a);
  mpz_clear (b);
  mpz_poly_clear (best_f);
  mpz_poly_clear (copy_f);
  mpz_clear (best_g0);
  mpz_clear (copy_g0);
  mpz_clear (k);
  mpz_clear (copy_g0k);
}


#define POLYLLL_NUM 1024

/* Pre-optimize routine: we assume that f[d], f[d-1] and f[d-2] are small,
   and we want to force f[d-3] to be small. */
static void
optimize_deg6_list ( mpz_poly_ptr f, mpz_t *g, const int verbose,
                     const int use_rotation )
{
  int d = f->deg;
#define MAXR 5 /* max. 3 roots of c3 and 2 roots of c4 */
  mpz_t k, r[MAXR], g0_copy, best_g0;
  long l;
  mpz_poly_t h, best_f, f_copy;
  LMAX = optimize_effort;
  ASSERT_ALWAYS(LMAX <= LMAX_MAX);
  MAXL = 2 * optimize_effort;
  mpz_init (k);
  mpz_init (best_g0);
  double skew, logmu, best_logmu = DBL_MAX, lognorm0;

  /* original lognorm */
  skew = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
  lognorm0 = L2_lognorm (f, skew);

  /* record top polynomials */
  polylll_pq *pqueue;
  new_polylll_pq (&pqueue, POLYLLL_NUM, d);

  /* holders, g[1] is not changed below, thus we
     only save g[0] */
  mpz_poly_init (f_copy, d);
  mpz_poly_init (best_f, d);
  mpz_poly_init (h, 3);
  for (int i = 0; i < MAXR; i++)
    mpz_init (r[i]);
  mpz_poly_set (f_copy, f);
  mpz_init_set (g0_copy, g[0]);

  /* add y*x^(d-2)*g and reduce coeff. of degree d-1 by
     translation. This is 0 for now. */
  for (int y = 0; y <= MAXY; y++) {

    /* recover original polynomial */
    mpz_poly_set (f, f_copy);
    mpz_set (g[0], g0_copy);

    /* add y*x^(d-2)*g and reduce */
    rotate_auxg_si (f->coeff, g, y, f->deg - 2);
    mpz_mul_si (r[0], f->coeff[d], - (long) d);
    mpz_ndiv_q (r[1], f->coeff[d-1], r[0]);
    do_translate_z (f, g, r[1]);

    /* idea by Thorsten Kleinjung: rotating by l*x^(d-3)*g
       for several values of l, keep best l */
    int nl = find_good_l (f, g), il;
    // if (nl > 0) printf ("nl=%d last v=%e\n", nl, best_v[nl-1]);
    best_l[nl++] = 0; /* add 0 as last element */
    int nl2 = find_good_l2 (f, g);

    for (il = 0; il < nl + nl2; il++) {

      l = (il < nl) ? best_l[il] : best_l2[il-nl];

      /* we consider f + l*x^(d-3)*g */
      mpz_poly_set (f, f_copy);
      mpz_set (g[0], g0_copy);

      /* f(x) := f(x) + l*x^3*g(x) */
      rotate_auxg_si (f->coeff, g, l, 3);

#if 1 /* consider both roots of c_{d-2}(k) and c_{d-3}(k) */
      /* h(k) = coefficient of x^3 in  f(x+k), which is a degree-3
         polynomial in k */
      fdminus3_translated (h, f);
      /* Store in r the integer approximation of the real roots of this
         cubic h(k), and the number of roots in nr_roots */
      int nr_roots = roots3 (r, h);
      fdminus2_translated (h, f);
      nr_roots += mpz_poly_mpz_roots (r + nr_roots, h);
#else /* only consider roots of c_{d-2}(k) */
      fdminus2_translated (h, f);
      int nr_roots = mpz_poly_mpz_roots (r, h);
      if (nr_roots == 0) {
        fdminus3_translated (h, f);
        nr_roots = roots3 (r, h);
      }
#endif

      /* For each root R, rounded to an integer, optimize f(x + R) */
      for (int j = 0; j < nr_roots; j++) {
        mpz_poly_set (f, f_copy);
        mpz_set (g[0], g0_copy);
        rotate_auxg_si (f->coeff, g, l, 3);

        do_translate_z (f, g, r[j]);

        /* after translation, c4/c3 are reduced */
        // print_poly_fg (f, g, g[0], 1);
        //fprintf (stdout, "# l: %ld, index: %d, nl1: %d, nl2: %d, j: %d\n", l, il, nl, nl2, j);

        /* pre-reduce smallest coefficients using rotation: this does not
           change the results of LLL reduction, but makes it faster */
        pre_reduce (f, g, k);

        optimize_lll_list (f, g, pqueue, lognorm0, verbose & use_rotation);
        skew = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
        logmu = L2_lognorm (f, skew);

        if (logmu < best_logmu) {
          //gmp_printf ("l=%ld k=%Zd lognorm=%.2f\n", l, r[j], logmu);
          best_logmu = logmu;
          mpz_poly_set (best_f, f);
          mpz_set (best_g0, g[0]);
        }

      } // for each root r of a translation k
    }  // for each translations k
  } // currently null

  /* output */
  mpz_set_ui (k, 0);
  for (int i=0;i<pqueue->used;i++) {
    extract_polylll_pq (pqueue, f->coeff, g, &logmu);
    fprintf (stdout, "\n# Size-optimized polynomial.\n");
    print_poly_fg (f, g, k, 1);
    //fprintf (stdout, "i: %d, lognorm: %f, used: %d\n", i, logmu, pqueue->used);
  }

  fprintf (stdout, "\n# used: %d\n", pqueue->used);

  /* set f and g to the best polynomials found */
  mpz_set (g[0], best_g0);
  mpz_poly_set (f, best_f);

  mpz_clear (k);
  mpz_clear (best_g0);
  mpz_poly_clear (h);
  mpz_poly_clear (best_f);
  mpz_poly_clear (f_copy);
  for (int i = 0; i < MAXR; i++)
    mpz_clear (r[i]);
  mpz_clear (g0_copy);
  free_polylll_pq (&pqueue);
  insert_hash (0);
}

#else

void
optimize_lll (mpz_poly_ptr f, mpz_t *g, int verbose)
{
  double skew, best_lognorm, lognorm, smax;
  mpz_t s, det, a, b;
  int d = f->deg, i, j, ss;
  mat_Z m;
  mpz_poly_t best_f, copy_f;
  mpz_t best_g0, copy_g0, k, copy_g0k;

  mpz_init (s);
  mpz_init (det);
  /* 1/4 < a/b < 1: the closer a/b is from 1, the better the reduction is.
     Some experiments suggest that increasing a/b does not yield better
     polynomials, thus we take the smallest possible value. */
  mpz_init_set_ui (a, 1);
  mpz_init_set_ui (b, 4);
  mpz_poly_init (best_f, d);
  mpz_poly_init (copy_f, d);
  mpz_poly_set (copy_f, f);
  mpz_init_set (copy_g0, g[0]);
  mpz_poly_set (best_f, f);
  mpz_init_set (best_g0, g[0]);
  best_lognorm = L2_skew_lognorm (f, SKEWNESS_DEFAULT_PREC);
  mpz_init (k);
  mpz_init (copy_g0k);

  /* m has d-1 (row) vectors of d+1 coefficients each */
  m.NumRows = d - 1;
  m.NumCols = d + 1;
  m.coeff = (mpz_t **) malloc (d * sizeof(mpz_t*));
  for (j = 1; j <= d-1; j++)
  {
    m.coeff[j] = (mpz_t *) malloc ((d + 2) * sizeof(mpz_t));
    for (i = 0; i <= d + 1; i++)
      mpz_init (m.coeff[j][i]);
  }
  smax = pow (fabs (mpz_get_d (g[0]) / mpz_get_d (f->coeff[d])),
              1.0 / (double) d);

  mpz_t best_k[MAX_ROOTS];
  int nk;
  for (i = 0; i < MAX_ROOTS; i++)
    mpz_init (best_k[i]);
  nk = find_best_k (best_k, f, g);
  ASSERT_ALWAYS(nk <= MAX_ROOTS);

  double best_lognorm_k = 0.0, old_best_lognorm_k = DBL_MAX;
  long kmin, kinc = 1;
  long kmax = kinc * ((optimize_effort + 1) / 2);
#define NSKEW 2 /* seems to be quite good */
  int optimize_effort2 = (int) sqrt ((double) optimize_effort);
  for (int ik = 0;; ik++) {
    if (ik < nk * optimize_effort2)
      {
	mpz_set (k, best_k[ik / optimize_effort2]);
	mpz_add_si (k, k, (ik % optimize_effort2) - (optimize_effort2 / 2));
      }
    else
      {
	/* We try values around 1, 10, 100, ... until we get a too large
	   average lognorm. The number of tried values depends on the
	   optimize effort. */
	if (ik == nk * optimize_effort2) /* initialize */
	  {
	    kmin = -kinc * (optimize_effort / 2);
	    mpz_set_si (k, kmin);
	    best_lognorm_k = DBL_MAX;
	  }
	else
	  {
	    mpz_add_ui (k, k, kinc);
	    if (mpz_cmp_si (k, kmax) > 0)
	      {
		if (best_lognorm_k > old_best_lognorm_k + 1.0)
		  break;
		old_best_lognorm_k = best_lognorm_k;
                if (kinc > LONG_MAX / 10)
                  break;
		kinc *= 10;
		kmin = -kinc * (optimize_effort / 2);
		kmax = kinc * ((optimize_effort + 1) / 2);
		mpz_set_si (k, kmin);
		best_lognorm_k = DBL_MAX;
	      }
	  }
	if (mpz_cmp_ui (k, 0) == 0)
	  continue; /* skip 0, already tried in best_k[] */
      }
    for (ss = NSKEW; ss <= 3 * NSKEW; ss++) {
      skew = pow (smax, (double) ss / (double) (2*NSKEW));
      mpz_poly_set (f, copy_f);
      mpz_set (g[0], copy_g0);
      do_translate_z (f, g, k);
      mpz_set (copy_g0k, g[0]);
      mpz_set_d (s, skew + 0.5); /* round to nearest */
      for (j = 1; j <= d-1; j++)
        for (i = 0; i <= d; i++)
          mpz_set_ui (m.coeff[j][i+1], 0);
      mpz_set_ui (det, 1);
      for (i = 0; i <= d; i++)
      {
        if (i > 0)
          mpz_mul (det, det, s); /* s^i */
        if (i <= d - 3)
          mpz_set (m.coeff[i+2][i+1], det);
        mpz_mul (m.coeff[1][i+1], det, f->coeff[i]);
      }
      for (i = 0; i <= d - 3; i++)
      {
        /* m.coeff[i+2][i+1] is already s^i */
        mpz_mul (m.coeff[i+2][i+2], m.coeff[i+2][i+1], s); /* s^(i+1) */
        mpz_mul (m.coeff[i+2][i+1], m.coeff[i+2][i+1], g[0]);
        mpz_mul (m.coeff[i+2][i+2], m.coeff[i+2][i+2], g[1]);
      }
      LLL (det, m, NULL, a, b);
      for (j = 1; j <= d-1; j++)
      {
        if (mpz_sgn (m.coeff[j][d]) != 0)
        {
          mpz_set_ui (det, 1);
          for (i = 0; i <= d; i++)
          {
            if (i > 0)
              mpz_mul (det, det, s); /* det = s^i */
            ASSERT_ALWAYS(mpz_divisible_p (m.coeff[j][i+1], det));
            mpz_divexact (f->coeff[i], m.coeff[j][i+1], det);
          }

	  /* force leading coefficient of f to be positive */
          if (mpz_sgn (f->coeff[d]) < 0)
            for (i = 0; i <= d; i++)
              mpz_neg (f->coeff[i], f->coeff[i]);

          if (insert_hash (mpz_getlimbn (f->coeff[0], 0)))
          {
            mpz_set (g[0], copy_g0k);
            optimize_aux (f, g, verbose, 1, OPT_STEPS);
            lognorm = L2_skew_lognorm (f, SKEWNESS_DEFAULT_PREC);
	    if (lognorm < best_lognorm_k)
	      best_lognorm_k = lognorm;
            if (lognorm < best_lognorm)
            {
              best_lognorm = lognorm;
              mpz_poly_set (best_f, f);
              mpz_set (best_g0, g[0]);
            }
          }
        }
      }
    }
  }

  for (i = 0; i < MAX_ROOTS; i++)
    mpz_clear (best_k[i]);

  for (j = 1; j <= d-1; j++)
  {
    for (i = 0; i <= d + 1; i++)
      mpz_clear (m.coeff[j][i]);
    free (m.coeff[j]);
  }
  free (m.coeff);
  mpz_poly_set (f, best_f);
  mpz_set (g[0], best_g0);
  mpz_clear (s);
  mpz_clear (det);
  mpz_clear (a);
  mpz_clear (b);
  mpz_poly_clear (best_f);
  mpz_poly_clear (copy_f);
  mpz_clear (best_g0);
  mpz_clear (copy_g0);
  mpz_clear (k);
  mpz_clear (copy_g0k);
}

/* Pre-optimize routine: we assume that f[d], f[d-1] and f[d-2] are small,
   and we want to force f[d-3] to be small. */
static void
optimize_deg6 (mpz_poly_ptr f, mpz_t *g, const int verbose,
               const int use_rotation)
{
  int d = f->deg;
#define MAXR 5 /* max. 3 roots of c3 and 2 roots of c4 */
  mpz_t k, r[MAXR], g0_copy, best_g0;
  long l;
  mpz_poly_t h, best_f, f_copy;

  LMAX = optimize_effort;
  //  ASSERT_ALWAYS(LMAX <= LMAX_MAX);

  MAXL = 2 * optimize_effort;

  mpz_init (k);
  mpz_init (best_g0);

#ifdef OPTIMIZE_MP
  mpz_t skew, logmu, best_logmu;
  mpz_init (skew);
  mpz_init_set_ui (logmu, 1);
  mpz_init_set_ui (best_logmu, 1);
  mpz_ui_pow_ui (best_logmu, 10, 100);
  double best_logmud = DBL_MAX;
#else
  double skew, logmu, best_logmu = DBL_MAX;
#endif

  mpz_poly_init (f_copy, d);
  mpz_poly_init (best_f, d);
  mpz_poly_init (h, 3);
  for (int i = 0; i < MAXR; i++)
    mpz_init (r[i]);

  /* g[1] is not changed below, thus we only save g[0] */
  mpz_poly_set (f_copy, f);
  mpz_init_set (g0_copy, g[0]);

  for (int y = 0; y <= MAXY; y++) {
    /* add y*x^(d-2)*g and reduce coeff. of degree d-1 by translation */
    mpz_poly_set (f, f_copy);
    mpz_set (g[0], g0_copy);
    rotate_auxg_si (f->coeff, g, y, f->deg - 2);
    mpz_mul_si (r[0], f->coeff[d], - (long) d);
    mpz_ndiv_q (r[1], f->coeff[d-1], r[0]);
    do_translate_z (f, g, r[1]);

    {
      l = 0;

      /* we consider f + l*x^(d-3)*g */
      mpz_poly_set (f, f_copy);
      mpz_set (g[0], g0_copy);

      /* f(x) := f(x) + l*x^3*g(x) */
      rotate_auxg_si (f->coeff, g, l, 3);

#if 0 /* consider both roots of c_{d-2}(k) and c_{d-3}(k) */
      /* h(k) = coefficient of x^3 in  f(x+k), which is a degree-3
         polynomial in k */
      fdminus3_translated (h, f);
      /* Store in r the integer approximation of the real roots of this
         cubic h(k), and the number of roots in nr_roots */
      int nr_roots = roots3 (r, h);
      fdminus2_translated (h, f);
      nr_roots += mpz_poly_mpz_roots (r + nr_roots, h);
#else /* only consider roots of c_{d-2}(k) */
      fdminus2_translated (h, f);
      int nr_roots = mpz_poly_mpz_roots (r, h);
      if (nr_roots == 0)
      {
        fdminus3_translated (h, f);
        nr_roots = roots3 (r, h);
      }
#endif

      /* For each root R, rounded to an integer, optimize f(x + R) */
      nr_roots = 1;
      mpz_set_ui (r[0], 0);
      for (int j = 0; j < nr_roots; j++)
      {
        mpz_poly_set (f, f_copy);
        mpz_set (g[0], g0_copy);
        rotate_auxg_si (f->coeff, g, l, 3);

        do_translate_z (f, g, r[j]);

        /* pre-reduce smallest coefficients using rotation: this does not
           change the results of LLL reduction, but makes it faster */
	pre_reduce (f, g, k);

#ifdef OPTIMIZE_MP
        optimize_aux_mp (f, g, verbose, use_rotation);
        L2_skewness_derivative_mp (f, SKEWNESS_DEFAULT_PREC, skew);
        L2_lognorm_mp (f, skew, logmu);

        double logmud = 0.0, skewd = 0.0;
        logmud = mpz_get_d (logmu);
        mpz_pow_ui (skew, skew, d); // s^6
        skewd = mpz_get_d (skew);
        logmud = logmud / skewd * 0.00043828022511018320850; /* Pi/7168 */
        logmud = 0.5 * log(logmud);

        if (logmud < best_logmud) {
          best_logmud = logmud;
          mpz_poly_set (best_f, f);
          mpz_set (best_g0, g[0]);
        }
#else
        // optimize_aux (f, g, verbose, use_rotation);
        optimize_lll (f, g, verbose & use_rotation);
        skew = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
        logmu = L2_lognorm (f, skew);

        if (logmu < best_logmu)
        {
          // gmp_printf ("l=%ld k=%Zd lognorm=%.2f\n", l, r[j], logmu);
          best_logmu = logmu;
          mpz_poly_set (best_f, f);
          mpz_set (best_g0, g[0]);
        }
#endif

      } // #roots
    }  // #translations
  }

  /* set f and g to the best polynomials found */
  mpz_set (g[0], best_g0);
  mpz_poly_set (f, best_f);

  mpz_clear (k);
  mpz_clear (best_g0);
  mpz_poly_clear (h);
  mpz_poly_clear (best_f);
  mpz_poly_clear (f_copy);
  for (int i = 0; i < MAXR; i++)
    mpz_clear (r[i]);
  mpz_clear (g0_copy);

#ifdef OPTIMIZE_MP
  mpz_clear (logmu);
  mpz_clear (skew);
  mpz_clear (best_logmu);
#endif

  insert_hash (0);
}
#endif


#if 0
/* try all translation -- this is old -- only for test */
static void
optimize_deg6_translate ( mpz_poly_ptr f, mpz_t *g, const int verbose,
                          const int use_rotation )
{
  mpz_t k0, k, g0_copy, best_g0, best_k, startk, endk, step;
  mpz_poly_t best_f, f_copy;
  double skew, logmu, best_logmu = DBL_MAX;
  int d = f->deg;
  mpz_init (k0);
  mpz_init (step);
  mpz_init (best_g0);
  mpz_init (best_k);
  mpz_init (startk);
  mpz_init (endk);
  mpz_poly_init (f_copy, d);
  mpz_poly_init (best_f, d);
  mpz_poly_set (f_copy, f);
  mpz_init_set (g0_copy, g[0]);

  // assuming monotonic, range 2^k
  mpz_tdiv_q (k0, g[0], f->coeff[d]);
  mpz_abs (k0, k0);
  mpz_root (k0, k0, d);
  mpz_mul (k0, k0, k0);

  mpz_init_set_si (k, 1);
  for (int i = 0; i < 30; i++, mpz_mul_si (k, k, 2)) {
    mpz_poly_set (f, f_copy);
    mpz_set (g[0], g0_copy);
    do_translate_z (f, g, k);
    optimize_aux (f, g, verbose, use_rotation, OPT_STEPS_FINAL);
    skew = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
    logmu = L2_lognorm (f, skew);
    gmp_fprintf (stdout, "# k: %Zd, logmu %f\n", k, logmu);
    if (logmu < best_logmu) {
      best_logmu = logmu;
      mpz_poly_set (best_f, f);
      mpz_set (best_g0, g[0]);
      mpz_set (best_k, k);
    }
  }

  // range -2^k
  mpz_init_set_si (k, -1);
  for (int i = 0; i < 30; i++, mpz_mul_si (k, k, 2)) {
    mpz_poly_set (f, f_copy);
    mpz_set (g[0], g0_copy);
    do_translate_z (f, g, k);
    optimize_aux (f, g, verbose, use_rotation, OPT_STEPS_FINAL);
    skew = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
    logmu = L2_lognorm (f, skew);
    gmp_fprintf (stdout, "# k: %Zd, logmu %f\n", k, logmu);
    if (logmu < best_logmu) {
      best_logmu = logmu;
      mpz_poly_set (best_f, f);
      mpz_set (best_g0, g[0]);
      mpz_set (best_k, k);
    }
  }

  /* more */
  mpz_init_set_si (k, 1);
  if (mpz_cmp_ui(best_k, 0) >= 0) {
    mpz_tdiv_q_ui (startk, best_k, 2);
    mpz_mul_ui (endk, best_k, 2);
  }
  else {
    mpz_tdiv_q_ui (endk, best_k, 2);
    mpz_mul_ui (startk, best_k, 2);
  }
  mpz_sub (step, startk, endk);
  mpz_abs (step, step);
  mpz_tdiv_q_ui (step, step, 1000);
  if (mpz_cmp_ui(step, 1) < 0)
    mpz_set_ui (step, 2);
  mpz_set (k, startk);
  while ( mpz_cmp(k, endk) < 0 ) {
    mpz_poly_set (f, f_copy);
    mpz_set (g[0], g0_copy);
    do_translate_z (f, g, k);
    optimize_aux (f, g, verbose, use_rotation, OPT_STEPS_FINAL);
    skew = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
    logmu = L2_lognorm (f, skew);
    //gmp_fprintf (stdout, "# k: %Zd, logmu %f\n", k, logmu);
    if (logmu < best_logmu) {
      gmp_fprintf (stdout, "# more k: %Zd, logmu %f\n", k, logmu);
      best_logmu = logmu;
      mpz_poly_set (best_f, f);
      mpz_set (best_g0, g[0]);
      mpz_set (best_k, k);
    }
    mpz_add (k, k, step);
  }

  // end
  mpz_poly_set (f, f_copy);
  mpz_set (g[0], g0_copy);
  do_translate_z (f, g, best_k);

  mpz_clear (k);
  mpz_clear (startk);
  mpz_clear (endk);
  mpz_clear (step);
  mpz_clear (best_k);
  mpz_clear (k0);
  mpz_clear (g0_copy);
  mpz_clear (best_g0);
  mpz_poly_clear (best_f);
  mpz_poly_clear (f_copy);
}
#endif

/* if use_rotation is non-zero, also use rotation */
void
optimize (mpz_poly_ptr f, mpz_t *g, const int verbose, const int use_rotation)
{
  const int d = f->deg;

  if (d == 6 || d == 5) {
#ifdef OPTIMIZE_LLL_LIST
    optimize_deg6_list (f, g, verbose, use_rotation);
#else
    optimize_deg6 (f, g, verbose, use_rotation);
#endif
    return; /* optimize_deg6 already performs optimize_aux */
  }

#ifdef OPTIMIZE_MP
  optimize_aux_mp (f, g, verbose, use_rotation);
#else
  optimize_aux (f, g, verbose, use_rotation, OPT_STEPS_FINAL);
#endif
}
