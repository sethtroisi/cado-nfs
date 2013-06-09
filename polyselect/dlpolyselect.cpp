/* 
   Test data
    P  = 31081938120519680804196101011964261019661412191103091971180537759
    (P - 1)/2 = 2 * Q where Q is a prime
    Q = 15540969060259840402098050505982130509830706095551545985590268879

    --run--
    ./dlpolyselect -n 31081938120519680804196101011964261019661412191103091971180537759 -d 4
    
    Large ones with deg 4 and 3,

    ./dlpolyselect -d 4 -N 191147927718986609689229466631454649812986246276667354864188503638807260703436799058776201365135161278134258296128109200046702912984568752800330221777752773957404540495707851421851

    ./dlpolyselect -d 3 -N 191147927718986609689229466631454649812986246276667354864188503638807260703436799058776201365135161278134258296128109200046702912984568752800330221777752773957404540495707851421851

*/

#include "cado.h"
#include "auxiliary.h"
#include "utils.h"
#include "portability.h"
#include "murphyE.h"
#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZVec.h>
#include <NTL/LLL.h>
#include <NTL/mat_ZZ.h>
#include <NTL/pair_ZZ_pX_long.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>

using namespace std;
using namespace NTL;

const double exp_rot[] = {0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0};

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

    skew[0] = L2_skewness (f, df, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    logmu[0] = L2_lognorm (f, df, skew[0], DEFAULT_L2_METHOD);
    alpha[0] = get_alpha (f, df, ALPHA_BOUND);
    skew[1] = L2_skewness (g, dg, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    logmu[1] = L2_lognorm (g, dg, skew[0], DEFAULT_L2_METHOD);
    alpha[1] = get_alpha (g, dg, ALPHA_BOUND);

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
  Test irreducibility over ZZx using NTL.
*/
static int
is_irreducible (ZZX f_ZZ) {
    Vec<Pair<ZZX, long> > factors;
    ZZ c = to_ZZ(1);
    factor (c, factors, f_ZZ, 0, 0);
    if (factors.length() != 1 )
        return 0;
    else
        return 1;
}


/*
  Find roots: square free and cantorzass.
  static int
  findRoot (ZZ_pX f_Fp, vec_pair_ZZ_pX_long& factors) {

  SquareFreeDecomp (factors, f_Fp);

  if (factors.length() != 1) {
  return 0;
  }
  else {
  CanZass (factors, f_Fp, 0);
  ZZ_pX ff;
  mul(ff, factors);
  if (f_Fp != ff)
  return 0;
  if (deg(factors[0].a) != 1)
  {
  return 0;
  }
  else
  return 1;
  }
  return 1;
  }
*/


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
    int nr, *fint;
    mpz_t t;
    mpz_init (t);
    fint = (int *) malloc ((d + 1)*sizeof(int));
    /* ntl */
    ZZX f_ZZ;

    /* find irreducible polynomial f */
    srand(time(NULL));
    while (1)
    {
        mpz_set_ui (f[d], ad);
        fint[0] = ad;
        for (i = 0; i < d; i ++) {
            fint[i] = (2*rand()- RAND_MAX) % bound;
            mpz_set_si (f[i], fint[i]);
        }

        /* content test (not necessary for now) */
        content_poly (t, f, d);
        if (mpz_cmp_ui(t, 1) != 0 ) {
            for (i = 0; i < d+1; i ++) {
                mpz_divexact (f[i], f[i], t);
            }
        }

        /* irreducibility test */
        f_ZZ.SetLength(d+1);
        for (i = 0; i <= d; i ++)
            SetCoeff(f_ZZ, i, fint[i]);
        if (is_irreducible(f_ZZ) == 0)
            continue;

        /* find roots */
        nr = poly_roots_mpz (rf, f, d, n);
        if (nr > 0)
            break;
    }

    mpz_clear (t);
    free(fint);
    return nr;
}


/*
  Generate polynomial g(x) of degree dg, given root
*/
static void
polygen_JL_g (mpz_t N, int dg, mpz_t **g, mpz_t root)
{
    char *stmp1;
    char *stmp2;
    ZZ zX, zY, zP;
    stmp1 = (char *) malloc (mpz_sizeinbase(N, 10) + 2);
    mpz_get_str (stmp1, 10, N);
    zP = to_ZZ (stmp1);
    stmp2 = (char *) malloc (mpz_sizeinbase(root, 10) + 2);
    mpz_get_str (stmp2, 10, root);
    zX = to_ZZ (stmp2);

    /* coeff matrix of dimention (dg+1)x(dg+1) */
    mat_ZZ Matrix;
    Matrix.SetDims (dg+1, dg+1);
    zY = zX;
    for (int i = 0;  i <= dg; i++) {
        for (int j = 0; j <= dg; j++) {
            if (i == 0) {
                if (j == 0)
                    Matrix [j][i] = zP;
                else {
                    Matrix [j][i] = zY;
                    Matrix [j][i] = zY * (-1);
                    zY = zY*zX;
                }
            }
            else
                Matrix [j][i] = to_ZZ( (i == j) );
        }
    }
		
    LLL (zY, Matrix, 0);
	
    /* only use first shortest coefficient */
    for (int k = 0; k <= dg; k ++ ) {
        for (int i = 0; i <= dg; i ++) {
            for (int j = 0; j < NumBits (Matrix [k] [i]); j ++)
                if (bit(Matrix[k][i], j) == 1)
                    mpz_setbit (g[k][i], j);
            if (Matrix[k][i] < 0)
                mpz_mul_si (g[k][i], g[k][i], -1);
        }
    }
    
    free(stmp1);
    free(stmp2);
}


/* JL method to generate d and d-1 polynomial */
static void
polygen_JL ( mpz_t n,
             unsigned int d,
             unsigned int bound,
             int ad )
{
    ASSERT_ALWAYS (d >= 3);
    unsigned int i, j, nr, dg, format=0;
    dg = d - 1;
    mpz_t *f, *rf, **g;
    f = (mpz_t *) malloc ((d + 1)*sizeof(mpz_t));
    rf = (mpz_t *) malloc ((d + 1)*sizeof(mpz_t));
    g = (mpz_t **) malloc ((dg + 1)*sizeof(mpz_t*));
    for (i = 0; i <= d; i ++) {
        mpz_init (f[i]);
        mpz_init (rf[i]);
    }
    for (i = 0; i <= dg; i ++) {
        g[i] = (mpz_t *) malloc ((dg+1)*sizeof(mpz_t));
        for (j = 0; j <= dg; j ++)
            mpz_init_set_ui (g[i][j], 0);
    }
    
    /* generate f of degree d of small coefficients */
    nr = polygen_JL_f (n, d, bound, f, rf, ad);

    for (i = 0; i < nr; i ++) {
        
        /* generate g of degree dg */
        polygen_JL_g (n, dg, g, rf[i]);

        for (j = 0; j < dg; j ++) {
            if (format == 1)
                gmp_printf ("n: %Zd\n", n);
            else
                gmp_printf ("N %Zd\n", n);
            print_nonlinear_poly_info (f, g[j], d, dg, format);
            if (format == 1)
                gmp_printf ("m: %Zd\n", n);
            else
                gmp_printf ("M %Zd\n\n", n);
        }
    }
    
    /* clear */
    for (i = 0; i <= d; i ++) {
        mpz_clear (f[i]);
        mpz_clear (rf[i]);
    }
    for (i = 0; i <= dg; i ++) {
        for (j = 0; j <= dg; j ++)
            mpz_clear (g[i][j]);
        free(g[i]);
    }
    free (f);
    free (g);
    free (rf);
}


static void
usage ()
{
    fprintf (stderr, "./dlpolyselect -n xxx -d xxx -ad xxx\n");
    exit (1);
}


int
main (int argc, char *argv[])
{
    int i, ad = 1;
    mpz_t N;
    unsigned int d = 0;
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
        else if (argc >= 3 && strcmp (argv[1], "-d") == 0) {
            d = atoi (argv[2]);
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
        exit (1);
    }

    if (d == 0) {
        fprintf (stderr, "Error, missing degree (-d option)\n");
        exit (1);
    }

    if (ad <= 0) {
        fprintf (stderr, "Error, need ad > 0 (-ad option)\n");
        exit (1);
    }

    unsigned int bound = 1000;
    unsigned int c = 0;
    while (c < bound*bound) {
        polygen_JL(N, d, bound, ad);
        c++;
    }
    mpz_clear (N);
}
