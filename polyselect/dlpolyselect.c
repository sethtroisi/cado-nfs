/* 
   Test data
    P  = 31081938120519680804196101011964261019661412191103091971180537759
    (P - 1)/2 = 2 * Q where Q is a prime
    Q = 15540969060259840402098050505982130509830706095551545985590268879

    --run--
    ./dlpolyselect -df 3 -N 31081938120519680804196101011964261019661412191103091971180537759
    
    ./dlpolyselect -df 4 -dg 3 -N 191147927718986609689229466631454649812986246276667354864188503638807260703436799058776201365135161278134258296128109200046702912984568752800330221777752773957404540495707851421851

    ./dlpolyselect -df 3 -dg 2 -N 191147927718986609689229466631454649812986246276667354864188503638807260703436799058776201365135161278134258296128109200046702912984568752800330221777752773957404540495707851421851

*/

#include "cado.h"
#include "auxiliary.h"
#include "area.h"
#include "utils.h"
#include "portability.h"
#include "murphyE.h"
#include <ctype.h>
#include <stdlib.h>
#include <time.h> 

const double exp_rot[] = {0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0};
#define SIZE(p) (((long *) (p))[1])
#define DATA(p) ((mp_limb_t *) (((long *) (p)) + 2))
#define swap(x, y) { long _tmp = (x); (x) = (y); (y) = _tmp; }
#define mpz_swap_n(x, y) { mpz_t *_tmp = (x); (x) = (y); (y) = _tmp; }

double area=AREA;
double bound_f=BOUND_F;
double bound_g=BOUND_G;

/* This uses Paul Zimmermann implementaion. I kept the header here:


   LLL using exact multiprecision arithmetic.
   Translated from NTL 4.1a <http://www.shoup.net/>
   into GMP <http://www.swox.se/gmp/> 
   by Paul Zimmermann, July 2000.

   Revised April 4, 2002 (bug found by Jens Franke <franke (at) math
   (dot) uni-bonn (dot) de>).

   This program is open-source software distributed under the terms 
   of the GNU General Public License <http://www.fsf.org/copyleft/gpl.html>.

   Usage: lll <mat_size> [a] [b] < file

   mat_size - size of the input matrix
   a, b - positive integer coefficients used for swapping vectors. 
   We should have 1/4 < delta=a/b <= 1. The closest delta is from 1,
   the shortest are the output vectors, but the computation takes longer.
 */

#define ROW
#define LEN_QQ 11

typedef struct {
  mpz_t **coeff;
  int NumRows, NumCols;
} mat_Z;

const unsigned int QQ[LEN_QQ] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};

static void
ident ( mat_Z X, long n )
{
   long i, j;
   X.NumRows = X.NumCols = n;
  
   for (i = 1; i <= n; i++)  
      for (j = 1; j <= n; j++)  
         if (i == j)
	    mpz_set_ui(X.coeff[i][j], 1);
         else  
            mpz_set_ui(X.coeff[i][j], 0);
} 


static void
InnerProduct ( mpz_t x, mpz_t *a, mpz_t *b, long n, mpz_t t1 )
{
   long i;

   mpz_set_ui(x, 0);
   for (i = 1; i <= n; i++) {
      mpz_mul(t1, a[i], b[i]);
      mpz_add(x, x, t1);
   }
}


static void
IncrementalGS ( mat_Z B, long *P, mpz_t *D, mpz_t **lam, long *s, long k )
{
   long n = B.NumCols;
   mpz_t u, t1, t2;
   long i, j, posj;

   mpz_init(u);
   mpz_init(t1);
   mpz_init(t2);

   for (j = 1; j <= k-1; j++) {
      posj = P[j];
      if (posj == 0) continue;

      InnerProduct(u, B.coeff[k], B.coeff[j], n, t1);
      for (i = 1; i <= posj-1; i++) {
         mpz_mul(t1, D[i], u);
         mpz_mul(t2, lam[k][i], lam[j][i]);
         mpz_sub(t1, t1, t2);
         mpz_div(t1, t1, D[i-1]);
         mpz_set(u, t1);
      }

      mpz_set(lam[k][posj], u);
   }

   InnerProduct(u, B.coeff[k], B.coeff[k], n, t1);

   for (i = 1; i <= *s; i++) {
      mpz_mul(t1, D[i], u);
      mpz_mul(t2, lam[k][i], lam[k][i]);
      mpz_sub(t1, t1, t2);
      mpz_div(t1, t1, D[i-1]);
      mpz_set(u, t1);
   }

   if (mpz_cmp_ui(u, 0) == 0)
     {
       P[k] = 0;
     }
   else
     {
       (*s)++;
       P[k] = *s;
       mpz_set (D[*s], u);
     }

   mpz_clear(u);
   mpz_clear(t1);
   mpz_clear(t2);
}


static void
BalDiv ( mpz_t q, mpz_t a, mpz_t d, mpz_t r )
/*  rounds a/d to nearest integer, breaking ties
    by rounding towards zero.  Assumes d > 0. */
{
   long cmp;

   mpz_fdiv_qr(q, r, a, d);

   mpz_mul_2exp(r, r, 1);

   cmp = mpz_cmp(r, d);
   if (cmp > 0 || (cmp == 0 && mpz_cmp_ui(q, 0) < 0))
      mpz_add_ui(q, q, 1);
}


static void
MulSub( mpz_t c, mpz_t c1, mpz_t c2, mpz_t x, mpz_t tmp )
/* c = c1 - x*c2 */
{
   mpz_mul(tmp, x, c2);
   mpz_sub(c, c1, tmp);
}


static void
MulSubN ( mpz_t *c, mpz_t *c2, mpz_t x, long n, mpz_t tmp )
/* c = c - x*c2 */
{
   long i;
   signed long int x0;

   x0 = mpz_get_si (x);
   if (mpz_cmp_si (x, x0) == 0 && 
       x0 != ((signed long int) 1 << (mp_bits_per_limb - 1))) {
     if (x0 > 0)
       for (i = 1; i <= n; i++) {
	 mpz_mul_ui (tmp, c2[i], x0);
	 mpz_sub (c[i], c[i], tmp);
       }
     else if (x0 < 0) {
       x0 = -x0;
       for (i = 1; i <= n; i++)
	 mpz_addmul_ui(c[i], c2[i], x0);
     }
   }
   else
     {
       for (i = 1; i <= n; i++)
         {
           mpz_mul (tmp, c2[i], x);
           mpz_sub (c[i], c[i], tmp);
         }
     }
}


static void
reduce ( long k, long l, mat_Z B, long *P, mpz_t *D, 
         mpz_t **lam, mat_Z* U, mpz_t t1, mpz_t r )
/* t1 and r are temporary variables */
{
    long j;

   if (P[l] == 0) return;

   mpz_mul_2exp (t1, lam[k][P[l]], 1);
   mpz_abs (t1, t1);
   if (mpz_cmp(t1, D[P[l]]) <= 0)
     return;

   BalDiv (r, lam[k][P[l]], D[P[l]], t1);
   MulSubN (B.coeff[k], B.coeff[l], r, B.NumRows, t1);

   if (U)
     MulSubN (U->coeff[k], U->coeff[l], r, B.NumRows, t1);

   for (j = 1; j <= l-1; j++)
     if (P[j] != 0)
       MulSub(lam[k][P[j]], lam[k][P[j]], lam[l][P[j]], r, t1);

   MulSub(lam[k][P[l]], lam[k][P[l]], D[P[l]], r, t1);
}


static long
SwapTest ( mpz_t d0, mpz_t d1, mpz_t d2, mpz_t lam,
           mpz_t a, mpz_t b, mpz_t t1, mpz_t t2)
/* test if a*d1^2 > b*(d0*d2 + lam^2)
   t1 and t2 are temporary variables */
{
    mpz_mul(t1, d0, d2);
    mpz_mul(t2, lam, lam);
    mpz_add(t1, t1, t2);
    mpz_mul(t1, t1, b);

    mpz_mul(t2, d1, d1);
    mpz_mul(t2, t2, a);

    return (mpz_cmp(t2, t1) > 0);
}


static void
MulAddDiv( mpz_t c, mpz_t c1, mpz_t c2, 
           mpz_t x, mpz_t y, mpz_t z, mpz_t t1, mpz_t t2 )

/* c = (x*c1 + y*c2)/z
   warning: c and z can be the same variable
   t1 and t2 are temporary variables */
{
    mpz_mul(t1, x, c1);
    mpz_mul(t2, y, c2);
    mpz_add(t1, t1, t2);
    mpz_divexact(c, t1, z);
}


static void
MulSubDiv ( mpz_t c, mpz_t c1, mpz_t c2, 
            mpz_t x, mpz_t y, mpz_t z, mpz_t t1 )
/* c = (x*c1 - y*c2)/z
   t1 is a temporary variable */
{
    mpz_mul(t1, x, c1);
    mpz_mul(c, y, c2);
    mpz_sub(t1, t1, c);
    mpz_divexact(c, t1, z);
}


static void
RowTransform ( mpz_t c1,
               mpz_t c2,
               mpz_t x,
               mpz_t y,
               mpz_t u,
               mpz_t v )
/* (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2) */
{
    mpz_t t1, t2;

    mpz_init(t1);
    mpz_init(t2);

    mpz_mul(t1, x, c1);
    mpz_mul(t2, y, c2);
    mpz_add(t1, t1, t2);

    mpz_mul(t2, u, c1);
    mpz_set(c1, t1);
    mpz_mul(t1, v, c2);
    mpz_add(c2, t1, t2);

    mpz_clear(t1);
    mpz_clear(t2);
}


static void
RowTransformN ( mpz_t *c1,
                mpz_t *c2,
                mpz_t x,
                mpz_t y,
                mpz_t u,
                mpz_t v,
                long n )
/* (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2) */
{
   mpz_t t1, t2;
   long i;

   mpz_init(t1);
   mpz_init(t2);

   for (i = 1; i <= n; i++)
     {
       mpz_mul(t1, x, c1[i]);
       mpz_mul(t2, y, c2[i]);
       mpz_add(t1, t1, t2);

       mpz_mul(t2, u, c1[i]);
       mpz_set(c1[i], t1);
       mpz_mul(t1, v, c2[i]);
       mpz_add(c2[i], t1, t2);
     }

   mpz_clear(t1);
   mpz_clear(t2);
}


static void
swapLLL ( long k, mat_Z B, long *P, mpz_t *D, 
          mpz_t **lam, mat_Z* U, long m, long verbose )
/* swaps vectors k-1 and k;  assumes P(k-1) != 0 */
{
   long i, j;
   mpz_t t1, t2, t3, e, x, y;

   mpz_init(t1);
   mpz_init(t2);
   mpz_init(t3);
   mpz_init(e);
   mpz_init(x);
   mpz_init(y);

   if (P[k] != 0) {
      if (verbose) fprintf(stderr, "swap case 1: %ld\n", k);

      mpz_swap_n(B.coeff[k-1], B.coeff[k]);
      if (U)
        mpz_swap_n(U->coeff[k-1], U->coeff[k]);
      
      for (j = 1; j <= k-2; j++)
         if (P[j] != 0)
            mpz_swap(lam[k-1][P[j]], lam[k][P[j]]);

      for (i = k+1; i <= m; i++) {
         MulAddDiv(t1, lam[i][P[k]-1], lam[i][P[k]],
                   lam[k][P[k]-1], D[P[k]-2], D[P[k]-1], t2, t3);
         MulSubDiv(lam[i][P[k]], lam[i][P[k]-1], lam[i][P[k]], 
                   D[P[k]], lam[k][P[k]-1], D[P[k]-1], t2);
         mpz_set(lam[i][P[k]-1], t1);
      }

      MulAddDiv(D[P[k]-1], D[P[k]], lam[k][P[k]-1],
                D[P[k]-2], lam[k][P[k]-1], D[P[k]-1], t2, t3);
   }
   else if (mpz_cmp_ui(lam[k][P[k-1]], 0) != 0) {
      if (verbose) fprintf(stderr, "swap case 2: %ld\n", k);
      mpz_gcdext(e, x, y, lam[k][P[k-1]], D[P[k-1]]);

      mpz_divexact(t1, lam[k][P[k-1]], e);
      mpz_divexact(t2, D[P[k-1]], e);

      mpz_set(t3, t2);
      mpz_neg(t2, t2);
      RowTransformN(B.coeff[k-1], B.coeff[k], t1, t2, y, x, B.NumRows);
      if (U)
        RowTransformN(U->coeff[k-1], U->coeff[k], t1, t2, y, x, B.NumRows);
      for (j = 1; j <= k-2; j++)
         if (P[j] != 0)
            RowTransform(lam[k-1][P[j]], lam[k][P[j]], t1, t2, y, x);

      mpz_mul(t2, t2, t2);
      mpz_divexact(D[P[k-1]], D[P[k-1]], t2);

      for (i = k+1; i <= m; i++)
         if (P[i] != 0) {
            mpz_divexact(D[P[i]], D[P[i]], t2);
            for (j = i+1; j <= m; j++) {
               mpz_divexact(lam[j][P[i]], lam[j][P[i]], t2);
            }
         }

      for (i = k+1; i <= m; i++) {
         mpz_divexact(lam[i][P[k-1]], lam[i][P[k-1]], t3);
      }

      swap(P[k-1], P[k]);
   }
   else {
      if (verbose) fprintf(stderr, "swap case 3: %ld\n", k);

      mpz_swap_n(B.coeff[k-1], B.coeff[k]);
      if (U)
        mpz_swap_n(U->coeff[k-1], U->coeff[k]);
   
      for (j = 1; j <= k-2; j++)
         if (P[j] != 0)
            mpz_swap(lam[k-1][P[j]], lam[k][P[j]]);

      swap(P[k-1], P[k]);
   }

   mpz_clear(t1);
   mpz_clear(t2);
   mpz_clear(t3);
   mpz_clear(e);
   mpz_clear(x);
   mpz_clear(y);
}


static long LLL ( mpz_t det, mat_Z B, mat_Z* U, mpz_t a,
                  mpz_t b, long verbose )
{
   long m, *P, j, s, k, max_k;
   mpz_t *D, **lam, tmp1, tmp2;

   mpz_init (tmp1);
   mpz_init (tmp2);

   m = B.NumRows;

   P = (long*) malloc((m+1) * sizeof(long));

   D = (mpz_t*) malloc((m+1) * sizeof(mpz_t));
   for (j=0; j<=m; j++)
     mpz_init_set_ui(D[j], j==0);

   lam = (mpz_t**) malloc((m+1) * sizeof(mpz_t*));
   for (j = 0; j <= m; j++)
     {
       lam[j] = (mpz_t*) malloc((m + 1) * sizeof(mpz_t));
       for (k = 0; k <= m; k++) mpz_init_set_ui(lam[j][k], 0);
     }

   if (U) ident(*U, m);

   s = 0;

   k = 1;
   max_k = 0;

   while (k <= m) {
      if (k > max_k)
        {
          IncrementalGS (B, P, D, lam, &s, k);
          max_k = k;
        }

      if (k == 1) {
         k++;
         continue;
      }

      reduce (k, k-1, B, P, D, lam, U, tmp1, tmp2);

      if (P[k-1] != 0 && 
          (P[k] == 0 || 
           SwapTest(D[P[k]], D[P[k]-1], D[P[k]-2], lam[k][P[k]-1], a, b, tmp1, tmp2))) {
         swapLLL (k, B, P, D, lam, U, max_k, verbose);
         k--;
      }
      else {	
         for (j = k-2; j >= 1; j--) 
            reduce(k, j, B, P, D, lam, U, tmp1, tmp2);
         k++;
      }
   }

   mpz_set(det, D[s]);
   for (j=0; j<=m; j++) mpz_clear(D[j]); free(D);
   for (j = 0; j <= m; j++) {
      for (k = 0; k <= m; k++) mpz_clear(lam[j][k]);
      free (lam[j]);
   }
   free (lam);

   mpz_clear(tmp1);
   mpz_clear(tmp2);

   free(P);
   return s;
}


/*
  Print two nonlinear poly info
*/
void
print_nonlinear_poly_info ( mpz_t *f,
                            mpz_t *g,
                            unsigned int df,
                            unsigned int dg,
                            int format )
{
    unsigned int i;
    double skew[2], logmu[2], alpha[2];
    mpz_poly_t ff;
    ff->deg = df;
    ff->coeff = f;
    mpz_poly_t gg;
    gg->deg = dg;
    gg->coeff = g;

    skew[0] = L2_skewness (ff, SKEWNESS_DEFAULT_PREC);
    logmu[0] = L2_lognorm (ff, skew[0]);
    alpha[0] = get_alpha (ff, ALPHA_BOUND);
    skew[1] = L2_skewness (gg, SKEWNESS_DEFAULT_PREC);
    logmu[1] = L2_lognorm (gg, skew[0]);
    alpha[1] = get_alpha (gg, ALPHA_BOUND);

    if (format == 1) {
        for (i = df + 1; i -- != 0; )
            gmp_printf ("c%u: %Zd\n", i, f[i]);
    }
    else {
        for (i = df + 1; i -- != 0; )
            gmp_printf ("X%u %Zd\n", i, f[i]);
    }
    if (format == 1) {
        for (i = dg + 1; i -- != 0; )
            gmp_printf ("Y%u: %Zd\n", i, g[i]);
    }
    else {
        for (i = dg + 1; i -- != 0; )
            gmp_printf ("Y%u %Zd\n", i, g[i]);
    }
    printf ("# f lognorm %1.2f, skew %1.2f, alpha %1.2f, E %1.2f, " \
            "exp_E %1.2f\n",
            logmu[0], skew[0], alpha[0], logmu[0] + alpha[0],
            logmu[0] - 0.824 * sqrt (2.0 * exp_rot[df] * log (skew[0])));
    printf ("# g lognorm %1.2f, skew %1.2f, alpha %1.2f, E %1.2f, " \
            "exp_E %1.2f\n",
            logmu[1], skew[1], alpha[1], logmu[1] + alpha[1],
            logmu[1] - 0.824 * sqrt (2.0 * exp_rot[dg] * log (skew[1])));
    printf ("# f+g score %1.2f\n", logmu[1] + alpha[1] + logmu[0] + alpha[0]);

}


/*
  Generate polynomial f(x) of degree d.
*/
static int
polygen_JL_f ( mpz_t n,
               unsigned int d,
               unsigned int bound,
               mpz_t *f,
               mpz_t *rf,
               int ad )
{
    unsigned int i;
    unsigned long *rq;
    int nr, *fint;
    mpz_t t;
    mpz_init (t);
    fint = (int *) malloc ((d + 1)*sizeof(int));
    rq = (unsigned long *) malloc ((d + 1)*sizeof(unsigned long));

    /* find irreducible polynomial f */
    while (1)
    {
        fint[d] = ad;
        mpz_set_ui (f[d], ad);
        for (i = 0; i < d; i ++) {
            fint[i] = (2*rand()- RAND_MAX) % bound;
            mpz_set_si (f[i], fint[i]);
        }

        /* content test (not necessary for now) */
        mpz_poly_t ff;
        ff->deg = d;
        ff->coeff = f;
        content_poly (t, ff);
        if (mpz_cmp_ui(t, 1) != 0 ) {
            for (i = 0; i < d+1; i ++) {
                mpz_divexact (f[i], f[i], t);
            }
        }

        /* irreducibility test */
        /*
        f_ZZ.SetLength(d+1);
        for (i = 0; i <= d; i ++)
            SetCoeff(f_ZZ, i, fint[i]);
        if (is_irreducible(f_ZZ) == 0)
            continue;
        */
        int test = 0;
        for (i = 0; i < LEN_QQ; i ++) {
            test = mpz_poly_roots_ulong (rq, ff, QQ[i]);
            if (test == 0)
                break;
        }
        if (test != 0)
            continue;

        /* find roots */
        nr = mpz_poly_roots_mpz (rf, f, d, n);
        if (nr > 0)
            break;
    }

    mpz_clear (t);
    free(fint);
    free(rq);
    return nr;
}


/*
  Generate polynomial g(x) of degree dg, given root
*/
static void
polygen_JL_g ( mpz_t N,
               int dg,
               mat_Z g,
               mpz_t root )
{
    int i, j;
    mpz_t a, b, det, r;
    mpz_init (det);
    mpz_init_set_ui(a, 3);
    mpz_init_set_ui(b, 4);
    mpz_init_set(r, root);
    for (i = 0; i <= dg + 1; i ++) {
        for (j = 0; j <= dg + 1; j ++) {
            mpz_set_ui (g.coeff[i][j], 0);
        }
    }

    for (i = 1;  i <= dg + 1; i++) {
        for (j = 1; j <= dg + 1; j++) {
            if (i == 1) {
                if (j == 1) {
                    mpz_set (g.coeff[j][i], N);
                }
                else {
                    mpz_set (g.coeff[j][i], r);
                    mpz_neg (g.coeff[j][i], g.coeff[j][i]);
                    mpz_mul (r, r, root);
                }
            }
            else {
                mpz_set_ui (g.coeff[j][i], i==j);
            }
        }
    }

    LLL(det, g, NULL, a, b, 0);

    mpz_clear (det);
    mpz_clear (a);
    mpz_clear (b);
    mpz_clear (r);
}


/* JL method to generate d and d-1 polynomial */
static void
polygen_JL ( mpz_t n,
             unsigned int df,
             unsigned int dg,
             unsigned int bound,
             int ad )
{
    ASSERT_ALWAYS (df >= 3);
    unsigned int i, j, nr, format=0;
    mpz_t *f, *rf;
    mat_Z g;
    f = (mpz_t *) malloc ((df + 1)*sizeof(mpz_t));
    rf = (mpz_t *) malloc ((df + 1)*sizeof(mpz_t));
    for (i = 0; i <= df; i ++) {
        mpz_init (f[i]);
        mpz_init (rf[i]);
    }
    g.coeff = (mpz_t **) malloc ((dg + 2)*sizeof(mpz_t*));
    g.NumRows = g.NumCols = dg + 1;
    for (i = 0; i <= dg + 1; i ++) {
        g.coeff[i] = (mpz_t *) malloc ((dg + 2)*sizeof(mpz_t));
        for (j = 0; j <= dg + 1; j ++) {
            mpz_init (g.coeff[i][j]);
        }
    }
    
    /* generate f of degree d of small coefficients */
    nr = polygen_JL_f (n, df, bound, f, rf, ad);

    for (i = 0; i < nr; i ++) {
        /* generate g of degree dg */
        polygen_JL_g (n, dg, g, rf[i]);

        for (j = 1; j <= dg + 1; j ++) {
            if (format == 1)
                gmp_printf ("n: %Zd\n", n);
            else
                gmp_printf ("N %Zd\n", n);

            print_nonlinear_poly_info (f, &((g.coeff[j])[1]), df, dg, format);
            if (format == 1)
                gmp_printf ("m: %Zd\n", n);
            else
                gmp_printf ("M %Zd\n\n", rf[i]);
        }
    }
    
    /* clear */
    for (i = 0; i <= df; i ++) {
        mpz_clear (f[i]);
        mpz_clear (rf[i]);
    }
    for (i = 0; i <= dg + 1; i ++) {
        for (j = 0; j <= dg + 1; j ++)
            mpz_clear (g.coeff[i][j]);
        free(g.coeff[i]);
    }
    free (f);
    free (g.coeff);
    free (rf);
}


static void
usage ()
{
    fprintf (stderr, "./dlpolyselect -N xxx -df xxx -dg xxx -ad xxx\n");
    exit (1);
}


int
main (int argc, char *argv[])
{
    int i, ad = 1;
    mpz_t N;
    unsigned int df = 0, dg = 0;
    mpz_init (N);

    /* printf command-line */
    printf ("#");
    for (i = 0; i < argc; i++)
        printf (" %s", argv[i]);
    printf ("\n");
    fflush (stdout);

    /* parsing */
    while (argc >= 3 && argv[1][0] == '-')
    {
        if (argc >= 3 && strcmp (argv[1], "-N") == 0) {
            mpz_set_str (N, argv[2], 10);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-df") == 0) {
            df = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-dg") == 0) {
            dg = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-ad") == 0) {
            ad = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else {
            fprintf (stderr, "Invalid option: %s\n", argv[1]);
            usage();
            exit (1);
        }
    }

    if (mpz_cmp_ui (N, 0) <= 0) {
        fprintf (stderr, "Error, missing input number (-N option)\n");
        usage ();
    }

    if (df == 0) {
        fprintf (stderr, "Error, error degree (-df option)\n");
        usage ();
    }

    if (dg == 0 || dg >= df) {
        fprintf (stderr, "Error, missing or error degree (-dg option)\n");
        fprintf (stderr, "       only support dg < df.\n");
        usage ();
    }

    if (ad <= 0) {
        fprintf (stderr, "Error, need ad > 0 (-ad option)\n");
        usage (); 
    }

    srand(time(NULL));
    unsigned int bound = 1000;
    unsigned int c = 0;
    while (c < bound*bound) {
        polygen_JL(N, df, dg, bound, ad);
        c++;
    }
    mpz_clear (N);
}
