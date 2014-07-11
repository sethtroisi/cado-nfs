#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include "bwc_config.h"
#include "portability.h"
#include "macros.h"
#include "utils.h"

int verbose = 0;

int nrhs = 0;
mpz_t prime;
char * rhsname;
FILE * rhsfile;
gmp_randstate_t rstate;

int * colweights;

double random_uniform(gmp_randstate_t rstate)
{
    mpf_t x;
    mpf_init2(x, 53);
    mpf_urandomb(x, rstate, 53);
    double y = mpf_get_d(x);
    mpf_clear(x);
    return y;
}

double random_normal_standard(gmp_randstate_t rstate)
{
    static int last = 0;
    static double vlast = 0;
    if (last) { --last; return vlast; }
    double u = random_uniform(rstate);
    double v = random_uniform(rstate);
    double rho = sqrt(-2*log(u));
    double theta = 2*M_PI*v;
    double x = rho*cos(theta);
    double y = rho*sin(theta);
    vlast=y;
    last++;
    return x;
}

double random_normal(gmp_randstate_t rstate, double mean, double sdev)
{
    return mean + sdev * random_normal_standard(rstate);
}

/* return the expected maximum among n picks for a normal low with given
 * mean and sdev */
double extreme_normal(double n, double mean, double sdev)
{
    /* See Cramer, Mathematical methods of statistics, p. 575 (7th
     * printing).
     */
    double x = sqrt(2*log(n))+(log(log(n))+log(4*M_PI)-2*0.5772)/(2*sqrt(2*log(n)));
    return mean + x * sdev;
}

double random_normal_constrained(gmp_randstate_t rstate, double mean, double sdev, double a, double b)
{
    for(;;) {
        double x = floor(random_normal(rstate, mean, sdev));
        if (x >= a && x < b) return x;
    }
}

/* This is the random variable associated to the *size* of the sample */
double random_binomial(gmp_randstate_t rstate, unsigned long n, double p)
{
    double mean = n * p;
    double sdev = sqrt(n * p * (1-p));
    return random_normal_constrained(rstate, mean, sdev, 0, n);
}

struct distribution_data_s {
    double alpha;
    int offset;         /* this controls the peakedness for the leftmost
                           columns. It is difficult to make this much
                           smaller than 32 presently. Quite unsafe to
                           change. */
    double scale;       /* event function is scale/(x+offset)^alpha */
    double mean;        /* computed */
    double sdev;        /* computed */
    unsigned long ncols;        /* only for constraint correction */
    unsigned long nrows;        /* informational */
};
typedef struct distribution_data_s distribution_data[1];
typedef struct distribution_data_s * distribution_data_ptr;

distribution_data F = {{.alpha=0.94, .offset=32, .scale=1}};


/* the probability function for value i is scale/(i+offset)^alpha.
 * the cumulative function (sum on [0,j[ ) is thus:
 *
 * scale * ((j+offset)^beta - offset^beta)/beta
 *
 * with beta = 1-alpha
 *
 * the reverse cumulative function is:
 *
 * (x*beta/scale+offset^beta)^(1/beta)-offset
 *
 * the variance is \sum p(i) - \sum p(i)^2  [proof left to reader]
 *
 * and based on this, the last term is:
 *
 * scale^2*((j+offset)^gamma- offset^gamma)/gamma
 *
 * with gamma = 1-2*alpha
 */

double dist_p(distribution_data_ptr f, double x)
{
    return f->scale*pow(x+f->offset,-f->alpha);
}

double dist_q(distribution_data_ptr f, double x)
{
    double beta = 1 - f->alpha;
    double u = f->scale / beta;
    return u * (pow(x + f->offset, beta) - pow(f->offset, beta));
}

double dist_qrev(distribution_data_ptr f, double y)
{
    double beta = 1 - f->alpha;
    double u = f->scale / beta;
    double r = pow(y / u + pow(f->offset, beta), 1 / beta) - f->offset;
    return r;
}

double dist_qq(distribution_data_ptr f, double x)
{
    double gamma = 1 - 2 * f->alpha;
    double v = f->scale * f->scale / gamma;
    return v * (pow(x + f->offset, gamma) - pow(f->offset, gamma));
}

void adjust_distribution_parameters(distribution_data_ptr f, unsigned long nrows, unsigned long ncols, double density)
{
    /* sets the scale parameter so that the expected row weight matches
     * our desired density target */
    f->scale = density / dist_q(f, nrows);
    f->mean = dist_q(f, nrows);
    f->sdev = sqrt(f->mean * f->mean - dist_qq(f, nrows));
    f->nrows = nrows;
    f->ncols = ncols;

    /* some checking and info */
    double p0 = dist_p(f, 0);
    if (p0 >= 1.0) {
        fprintf(stderr, "Error: this density is not acceptable for the current distribution equation. Please adjust the internal offset parameter to something larger.\nrows");
        exit(1);
    }
    double mean0 = nrows * p0;
    double sdev0 = sqrt(nrows * p0 * (1-p0));
    double pn = dist_p(f, ncols-1);
    double mean_n = nrows * pn;
    double sdev_n = sqrt(nrows * pn * (1-pn));
    printf("Expected row weight: %.3f, sdev %.3f\nrows", f->mean, f->sdev);
    printf("Expected weight for first column is %.3f (sdev %.3f, m/sdev=%.1f)\n",
            mean0, sdev0, mean0 / sdev0);
    printf("Expected weight for last column is %.3f (sdev %.3f, m/sdev=%.1f)\n",
            mean_n, sdev_n, mean_n / sdev_n);
    printf("Tail size for rejected normal picks is 2^%.2f\n",
            log2(1-erf(mean_n / sdev_n)));
    printf("Worst-case expectation for last column weight by normal approximation: %.3f\n",
            extreme_normal(nrows, mean_n, -sdev_n));
}

struct punched_interval_s;
struct punched_interval_s {
    double b0, b1;
    double holes;
    struct punched_interval_s * left;
    struct punched_interval_s * right;
};
typedef struct punched_interval_s * punched_interval_ptr;

void punched_interval_free(punched_interval_ptr c)
{
    if (!c) return;
    punched_interval_free(c->left);
    punched_interval_free(c->right);
    free(c);
}

punched_interval_ptr punched_interval_alloc(double b0, double b1)
{
    punched_interval_ptr x = malloc(sizeof(struct punched_interval_s));
    memset(x, 0, sizeof(struct punched_interval_s));
    x->b0 = b0;
    x->b1 = b1;
    return x;
}

unsigned long pick_and_punch(distribution_data_ptr f, punched_interval_ptr c, double x)
{
    /* x should be within [c->b0, c->b1 - c->holes] */
    if (c->left == NULL) {
        /* no holes ! */
        double r = dist_qrev(f, x);
        unsigned long i;
        if (r < 0) {
            i = 0;
        } else if (r >= f->ncols) {
            i = f->ncols - 1;
        } else {
            i = floor(r);
        }
        double x0 = dist_q(f, i);
        double x1 = dist_q(f, i + 1);
        c->holes += x1 - x0;
        c->left = punched_interval_alloc(c->b0, x0);
        c->right = punched_interval_alloc(x1, c->b1);
        return i;
    }
    /* try to correct x with all left holes */
    double xc = x + c->left->holes;
    if (xc < c->left->b1) {
        double h = c->left->holes;
        unsigned long i = pick_and_punch(f, c->left, x);
        c->holes += c->left->holes - h;
        return i;
    } else {
        /* modify x. It's more than just xc ! */
        x = xc + c->right->b0 - c->left->b1;
        double h = c->right->holes;
        unsigned long i = pick_and_punch(f, c->right, x);
        c->holes += c->left->holes - h;
        return i;
    }
}

typedef int (*sortfunc_t)(const void *, const void *);

int cmp_ulong(unsigned long * a, unsigned long * b)
{
    return (*a > *b) - (*b > *a);
}

unsigned long generate_row(gmp_randstate_t rstate, distribution_data_ptr f, unsigned long * ptr)
{
    /* pick a row weight */
    unsigned long weight = random_normal_constrained(rstate, f->mean, f->sdev, 0, f->ncols);
    punched_interval_ptr range = punched_interval_alloc(0, f->mean);
    for(unsigned long i = 0 ; i < weight ; i++) {
        double x = random_uniform(rstate) * (range->b1 - range->holes);
        unsigned long k = pick_and_punch(f, range, x);
        ptr[i] = k;
    }
    qsort(ptr, weight, sizeof(unsigned long), (sortfunc_t) &cmp_ulong);
    punched_interval_free(range);
    return weight;
}


void printrows(gmp_randstate_t rstate, distribution_data_ptr f, unsigned long nrows, unsigned long ncols, int maxc)
{
    unsigned long * ptr = malloc(ncols * sizeof(unsigned long));
    int total_coeffs = 0;
    double tot_sq = 0;
    for(unsigned long i = 0 ; i < nrows ; i++) {
        long v = 0;
        unsigned long c = generate_row(rstate, f, ptr);
        printf("%lu", c);
        for(unsigned long j = 0 ; j < c ; j++) {
            printf(" %lu", ptr[j]);
            if (maxc) {
                int co = gmp_urandomm_ui(rstate, 2 * maxc + 1) - maxc;
                printf(":%d", co);
                if (nrhs) v += co * (1+ptr[j]);
            }
        }
        printf("\n");
        if (nrhs) {
            mpz_t x,s;
            mpz_init(x);
            mpz_init(s);
            mpz_set_si(s, v);

            for(int j = 0 ; j < nrhs - 1 ; j++) {
                mpz_urandomm(x, rstate, prime);
                mpz_addmul_ui(s, x, ncols + j + 1);
                gmp_fprintf(rhsfile, "%Zd ", x);
            }
            mpz_set_si(x, -(ncols + nrhs));
            mpz_invert(x, x, prime);
            mpz_mul(s, s, x);
            mpz_mod(s, s, prime);
            gmp_fprintf(rhsfile, "%Zd\n", s);
            mpz_clear(x);
            mpz_clear(s);
        }

        total_coeffs += c;
        tot_sq += (double) c * (double) c;
    }
    free(ptr);
    if (verbose) {
        double e = (double) total_coeffs / nrows;
        double s = (double) tot_sq / nrows;
        double sdev = sqrt(s - e*e);
        fprintf(stderr, "Actual density per row avg %.2f sdev %.2f\n",
                e, sdev);
    }
}


#if 0
int gen_row(double lambda, int n, int * ptr)
{
    int c = 0;
    for(int i = 0 ; ; ) {
        int e = 1 + (int) (lambda * dist_func() * (i+1));
        i += e;
        c += 1;
        if (i < 0 || i > n)
            break;
        if (ptr) *ptr++ = i - 1;
    }
    return c-1;
}

void printrows(int nrows, int ncols, double density, int maxc)
{
    int * ptr = malloc((ncols + 1) * sizeof(int));
    int total_coeffs = 0;
    double tot_sq = 0;
    for(int i = 0 ; i < nrows ; i++) {
        long v = 0;
        int c = gen_row(lambda, ncols, ptr);
        printf("%d", c);
        for(int j = 0 ; j < c ; j++) {
            printf(" %d", ptr[j]);
            if (maxc) {
                int co = gmp_urandomm_ui(rstate, 2 * maxc + 1) - maxc;
                printf(":%d", co);
                if (nrhs) v += co * (1+ptr[j]);
            }
        }
        printf("\n");
        if (nrhs) {
            mpz_t x,s;
            mpz_init(x);
            mpz_init(s);
            mpz_set_si(s, v);

            for(int j = 0 ; j < nrhs - 1 ; j++) {
                mpz_urandomm(x, rstate, prime);
                mpz_addmul_ui(s, x, ncols + j + 1);
                gmp_fprintf(rhsfile, "%Zd ", x);
            }
            mpz_set_si(x, -(ncols + nrhs));
            mpz_invert(x, x, prime);
            mpz_mul(s, s, x);
            mpz_mod(s, s, prime);
            gmp_fprintf(rhsfile, "%Zd\n", s);
            mpz_clear(x);
            mpz_clear(s);
        }

        total_coeffs += c;
        tot_sq += (double) c * (double) c;
    }
    free(ptr);
    if (verbose) {
        double e = (double) total_coeffs / nrows;
        double s = (double) tot_sq / nrows;
        double sdev = sqrt(s - e*e);
        fprintf(stderr, "Actual density per row avg %.2f sdev %.2f\n",
                e, sdev);
    }
}
#endif

void usage()
{
    fprintf(stderr, "Usage: ./random <nrows> [<ncols>] [<density>] [options]\n"
            "Options:\n"
            "\t-d <density> : desired density per row\n"
            "\t-s <seed> : seed\n"
            "\t-c <maxc> : add coefficients\n"
            "\t-v : turn verbosity on\n"
            "\t--kleft <d>: ensure at least a left kernel of dimension d\n"
            "\t--kright <d>: ditto for right kernel\n"
            "\t--rhs <nrhs>,<prime>,<filename>: rhs output\n"
            );
    exit(1);
}


int main(int argc, char * argv[])
{
    param_list pl;

    param_list_init(pl);

    argv++, argc--;
    int wild = 0;

    int maxcoeff = 0;

    int wild_args[3] = { -1, -1, -1 }; // nrows ncols coeffs_per_row

    param_list_configure_alias(pl, "density", "-d");
    param_list_configure_alias(pl, "seed", "-s");
    param_list_configure_switch(pl, "-v", &verbose);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

        if (argv[0][0] != '-' && wild < 3) {
            char * tmp;
            wild_args[wild] = strtoul(argv[0], &tmp, 0);
            if (*tmp != '\0') {
                fprintf(stderr, "Parse error for parameter %s\n", argv[0]);
                exit(1);
            }
            argv++, argc--;
            wild++;
            continue;
        }

        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage();
    }

    int seed = 0;
    param_list_parse_int(pl, "seed", &seed);

    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, seed ? seed : time(NULL));


    param_list_parse_int(pl, "c", &maxcoeff);

    int nrows = wild_args[0];
    if (nrows < 0) {
        fprintf(stderr, "Please specify nrows\n");
        exit(1);
    }
    int ncols = wild_args[1];
    if (ncols < 0) {
        ncols = nrows;
    }
    int density = wild_args[2];      // average coeffs per row
    if (density < 0) {
        param_list_parse_int(pl, "density", &density);
    }

    ASSERT_ALWAYS(ncols > 10 && nrows > 10);

    /* default is to make the matrix invertible */
    int kernel_left = 0;
    param_list_parse_int(pl, "kleft", &kernel_left);

    int kernel_right = 0;
    param_list_parse_int(pl, "kright", &kernel_right);
    
    /* we've been given dimensions for the target matrix, together with a
     * constraint on the kernel size (to be understood as "at least that
     * large").
     *
     * given the row/col unbalance, we may or may not have to generate
     * fewer rows/cols.
     */
    if (ncols > nrows) kernel_right -= ncols - nrows;
    if (kernel_right < 0) kernel_right = 0;

    if (nrows > ncols) kernel_left -= nrows - ncols;
    if (kernel_left < 0) kernel_left = 0;

    if (kernel_right > nrows / 4) {
        fprintf(stderr, "Warning, right kernel is large."
                " Could trigger misbehaviours\n");
    }

    if (kernel_left > ncols / 4) {
        fprintf(stderr, "Warning, left kernel is large."
                " Could trigger misbehaviours\n");
    }

    const char * tmp = param_list_lookup_string(pl, "rhs");

    if (param_list_warn_unused(pl)) usage();

    if (tmp) {
        /* try to parse the rhs info */
        ASSERT_ALWAYS(maxcoeff > 0);
        ASSERT_ALWAYS(kernel_left == 0);
        kernel_right = 0;
        char * rhsname = malloc(1 + strlen(tmp));
        mpz_init(prime);
        int rc = gmp_sscanf(tmp, "%d,%Zd,%s", &nrhs, prime, rhsname);
        ASSERT_ALWAYS(rc == 3);
        if (nrhs == 0) {
            fprintf(stderr, "--rhs argument requires setting more than 0 vectors !\n");
            exit(1);
        }
        rhsfile = fopen(rhsname, "w");
        fprintf(rhsfile, "%d %d\n", nrows, nrhs);
    }


    colweights = malloc(ncols * sizeof(int));
    memset(colweights, 0, ncols * sizeof(int));

    adjust_distribution_parameters(F, nrows - kernel_right, ncols - kernel_left, density);

    printf("%d %d\n", nrows, ncols);
    printrows(rstate, F, nrows - kernel_right, ncols - kernel_left, maxcoeff);
    for(int i = 0 ; i < kernel_right ; i++) {
        printf("0\n");
    }

    free(colweights);

    gmp_randclear(rstate);

    return 0;
}
