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
#define LEN_QQ 11

const unsigned int QQ[LEN_QQ] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};

double area=AREA;
double bound_f=BOUND_F;
double bound_g=BOUND_G;

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
        mpz_poly_content (t, ff);
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
        nr = mpz_poly_roots_mpz (rf, ff, n);
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

    LLL (det, g, NULL, a, b);

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
