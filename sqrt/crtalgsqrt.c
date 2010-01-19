#define _BSD_SOURCE     /* M_LN2 */
/*
 * Program: crtalgsqrt
 * Authors: E. Thomé.
 * Purpose: computing the squareroots and finishing the factorization
 *
 */

/* TODO list.
 *
 * This program seems to work, but is not complete.
 *
 * For correctness:
 *
 * - finish binding with the rest: check also the rational square root,
 *   and produce the factorization. Not necessarily independent from the
 *   above, since the elementary check is whether we get x^2=a mod N.
 * - update this TODO list.
 *
 * For speed / memory:
 *
 * - The program is probably overzealous with mpz_mod's sometimes. Some
 *   can be s(h)aved.
 * - I havent' check the peak memory usage. Notwithstanding the extra
 *   memory amount used by full-length multiplications, it should be 3
 *   times the ram_gb parameter.
 *
 * Important in terms of functionality, but not critical:
 *
 * - Use fpLLL for solving the knapsack. Should make it possible to go to
 *   12 primes or so in degree 6 (at least).
 * - MPI-ify the relevant bits. The two main loops on pgnum and apnum can
 *   be split in an s times t process grid.
 *
 * Not important:
 *
 * - understand and, if relevant, fix the memory leak diagnosed by
 *   valgrind. I don't understand.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <complex.h>
#include <float.h>
#include <time.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include "cado.h"
#include "utils.h"
#include "modul_poly.h"
#include "powers_of_p.h"
#include "polyroots.h"

static int verbose = 0;
double print_delay = 1;
double ram_gb = 3.0;    // Number of gigabytes. Note that this is the
                        // maximum size of product that are obtained.
                        // The actual memory footprint can be quite
                        // considerably larger due to FFT allocation (by
                        // a constant factor, though).

static void usage()
{
    fprintf(stderr, "usage: crtalgsqrt algdepfile ratdepfile polyfile\n");
    exit(1);
}

// {{{ interface for reading the list of (a,b)'s, with sort of a random
// access (for the starting point only).

#define ABFILE_MAX_LINE_LENGTH  256

struct ab_source_s {
    const char * fname0;
    char * sname;
    size_t sname_len;
    size_t nab;
    char prefix[100];
    int depnum;
    // int cado; // use !numfiles instead.

    /* this relates to rough estimations based on the size(s) of the
     * file(s). Note however that apparently the guess logic isn't too
     * good at taking the leading coefficient into account, so it's not
     * meant to be used before it gets fixed. Doing the accurate
     * evaluation via complex embeddings is ultra-fast anyway */
    size_t nab_estim;
    size_t digitbytes_estim;

    /* This relates to the different files (if several), their sizes, and
     * the current position. */
    size_t * file_bases;
    int nfiles; // 0 for cado format.
    size_t totalsize;

    FILE * f;
    int c;
    size_t cpos;        // position within current file.
    size_t tpos;        // position within totality.
};

typedef struct ab_source_s ab_source[1];
typedef struct ab_source_s * ab_source_ptr;

void ab_source_init(ab_source_ptr ab, const char * fname)
{
    memset(ab, 0, sizeof(ab_source));
    ab->fname0 = fname;
    char * magic;
    if ((magic = strstr(fname, ".prep.")) != NULL) {
        // then assume kleinjung format.
        strncpy(ab->prefix, fname, magic-fname);
        ab->prefix[magic-fname]='\0';
        magic++;
        int fnum;
        if (sscanf(magic, "prep.%d.rel.%d", &ab->depnum, &fnum) == 2) {
            ab->nfiles = -1;  // to be determined later on.
        } else {
            FATAL_ERROR_CHECK(1, "error in parsing filename");
        }
    } else if ((magic = strstr(fname, ".dep.alg.")) != NULL) {
        // assume cado format (means only one file, so we don't need to
        // parse, really.
        strncpy(ab->prefix, fname, magic-fname);
        ab->prefix[magic-fname]='\0';
        magic++;
        if (sscanf(magic, "dep.alg.%d", &ab->depnum) == 1) {
            ab->nfiles = 0;
        } else {
            FATAL_ERROR_CHECK(1, "error in parsing filename");
        }
    } else {
        FATAL_ERROR_CHECK(1, "error in parsing filename");
    }

    // do some size estimations;
    size_t tsize = 0;
    struct stat sbuf[1];
    int rc;
    if (ab->nfiles == 0) {
        rc = stat(fname, sbuf);
        ASSERT_ALWAYS(rc == 0);
        tsize=sbuf->st_size;
        // we have 2.5 non-digit bytes per file line. However we
        // don't know the line count, so we can't subtract. As a
        // guess, we read the first 16kb, and count the number of
        // lines in there.
        char buf[16384];
        FILE * f = fopen(fname, "r");
        fread(buf, 1, sizeof(buf), f);
        fclose(f);
        int nrows_16k = 0;
        for(unsigned int i = 0 ; i < sizeof(buf) ; i++) {
            nrows_16k += buf[i] == '\n';
        }
        ab->nab_estim = (double) tsize * nrows_16k / sizeof(buf);
        ab->digitbytes_estim = tsize - 2.5 * ab->nab_estim;
        ab->file_bases = malloc(2 * sizeof(size_t));
        ab->file_bases[0] = 0;
        ab->file_bases[1] = tsize;
        ab->totalsize = tsize;
    } else {
        FILE * f = fopen(fname, "r");
        size_t dummy;
        size_t hdrbytes;
        char line[ABFILE_MAX_LINE_LENGTH];
        char * xx = fgets(line, sizeof(line), ab->f);
        DIE_ERRNO_DIAG(xx == NULL, "fgets", fname);
        rc = sscanf(line, "AB %zu %zu", &ab->nab, &dummy);
        DIE_ERRNO_DIAG(rc != 2, "parse", fname);
        hdrbytes = ftell(f);
        fclose(f);
        ab->sname_len = strlen(fname) + 10;
        ab->sname = malloc(ab->sname_len);
        for(ab->nfiles = 0 ; ; ab->nfiles++) {
            snprintf(ab->sname, ab->sname_len, "%s.prep.%d.rel.%d",
                    ab->prefix, ab->depnum, ab->nfiles);
            rc = stat(ab->sname, sbuf);
            ASSERT_ALWAYS(rc == 0 || errno == ENOENT);
            if (rc < 0) break;
            tsize += sbuf->st_size - hdrbytes;
        }
        ASSERT_ALWAYS(ab->nfiles > 0);
        ab->nab_estim = ab->nab;
        ab->digitbytes_estim = tsize - 5 * ab->nab_estim;
        ab->file_bases = malloc((ab->nfiles+1) * sizeof(size_t));
        ab->file_bases[0] = 0;
        for(int i = 0 ; i < ab->nfiles ; i++) {
            snprintf(ab->sname, ab->sname_len, "%s.prep.%d.rel.%d",
                    ab->prefix, ab->depnum, ab->nfiles);
            rc = stat(ab->sname, sbuf);
            ASSERT_ALWAYS(rc == 0);
            ab->file_bases[i+1]=ab->file_bases[i] + sbuf->st_size;
        }
        ab->totalsize = ab->file_bases[ab->nfiles];
    }
    fprintf(stderr, "# [%2.2lf] %s: roughly %zu rows,"
            " %zu digits (%zu bits/c)\n",
            seconds(),
            ab->nfiles ? "kleinjung" : "cado",
            ab->nab_estim, ab->digitbytes_estim,
            (size_t) (ab->digitbytes_estim * M_LN10 / M_LN2));
}

void ab_source_rewind(ab_source_ptr ab)
{
    if (ab->f) fclose(ab->f);
    ab->f = NULL;
    ab->c = 0;
    ab->nab = 0;
    ab->cpos = 0;
    ab->tpos = 0;
}

void ab_source_clear(ab_source_ptr ab)
{
    ab_source_rewind(ab);
    free(ab->sname);
    free(ab->file_bases);
    memset(ab, 0, sizeof(ab_source));
}

int ab_openfile_internal(ab_source_ptr ab)
{
    const char * s;
    if (ab->nfiles == 0) {
        ab->f = fopen(s=ab->fname0, "r");
        ab->tpos = 0;
    } else {
        snprintf(ab->sname, ab->sname_len, "%s.prep.%d.rel.%d",
                ab->prefix, ab->depnum, ab->c);
        ab->f = fopen(s=ab->sname, "r");
        if (ab->f == NULL && errno == ENOENT)
            return 0;
        ab->tpos = ab->file_bases[ab->c];
        ab->cpos = 0;
        char header[80];
        fgets(header, sizeof(header), ab->f);
    }
    ab->cpos = ftell(ab->f);
    ab->tpos += ab->cpos;
    DIE_ERRNO_DIAG(ab->f == NULL, "fopen", s);
    return 1;
}

int ab_source_next(ab_source_ptr ab, int64_t * a, uint64_t * b)
{
    if (ab->f) {
        int rc;
        char line[ABFILE_MAX_LINE_LENGTH];
        char * xx = fgets(line, sizeof(line), ab->f);
        size_t cpos = ftell(ab->f);
        if (xx) {
            if (ab->nfiles == 0) {
                rc = sscanf(line, "%" SCNd64 " %" SCNu64, a, b);
                ASSERT(rc == 2);
            } else {
                int dummy;
                rc = sscanf(line, "%d %" SCNd64 " %" SCNu64, &dummy, a, b);
                ASSERT(rc == 3);
            }
            ab->tpos += cpos - ab->cpos;
            ab->cpos = cpos;
            ab->nab++;
            return 1;
        }
        fclose(ab->f); ab->f = NULL;
        ab->tpos += cpos - ab->tpos;
        // don't update cpos, it is not defined in this situation.
        if (ab->nfiles == 0)
            return 0;
        ab->c++;
    }
    if (ab_openfile_internal(ab) == 0)
        return 0;
    return ab_source_next(ab, a, b);
}

void ab_source_move_afterpos(ab_source_ptr ab, size_t offset)
{
    /* move the file pointer to the earliest non-header line starting at
     * an offset which is greater than or equal to offset.
     *
     * IF the current file position is already >= offset, do nothing.
     *
     * IF the current file position is < offset, then the returned offset
     * is >= offset. This is achieved by seeking to some pre_offset <
     * offset, and
     * advancing to the next line starting at >= offset.
     *
     * This is used to read a collection of files in chunks whose size is
     * governed by the file size.
     */

    /* note that the test below always succeeds for offset == 0 */
    if (ab->tpos >= offset)
        return;
    // otherwise it's really, really a can of worms.
    ASSERT_ALWAYS(ab->f == NULL);

    // which file ?
    FATAL_ERROR_CHECK(offset >= ab->totalsize,
            "attempt to seek beyond end of files");

    size_t pre_offset = MAX(offset, 10) - 10;
    ASSERT_ALWAYS(pre_offset < offset);

    if (ab->nfiles == 0) {
        // well, does not really make a lot of sense here, but anyway.
        ab_openfile_internal(ab);
        fseek(ab->f, pre_offset, SEEK_SET);
    } else {
        for( ; ab->file_bases[ab->c+1] <= pre_offset ; ab->c++) ;
        ab_openfile_internal(ab);
        // so we know that
        // ab->file_bases[ab->c] <= pre_offset < ab->file_bases[ab->c+1]
        fseek(ab->f, pre_offset - ab->file_bases[ab->c], SEEK_SET);
        // note that it is possible that tpos, as obtained after
        // openfile, is >= offset. but we know that we'll be able to seek
        // to pre_offset, and that one is appropriately < offset, so no
        // special case.
    }
    char line[ABFILE_MAX_LINE_LENGTH];
    char * xx = fgets(line, sizeof(line), ab->f);
    DIE_ERRNO_DIAG(xx == NULL, "fgets", ab->nfiles ? ab->sname : ab->fname0);
    size_t cpos = ftell(ab->f);
    ab->tpos += cpos - ab->cpos;
    ab->cpos = cpos;
    for(int n_adjust = 0 ; ab->tpos < offset ; n_adjust++) {
        FATAL_ERROR_CHECK(n_adjust > 10, "adjustment on the runaway");
        int64_t a;
        uint64_t b;
        int r = ab_source_next(ab, &a, &b);
        FATAL_ERROR_CHECK(r == 0, "adjustment failed");
    }
    fprintf(stderr, "# [%2.2lf] (a,b) rewind to %s, pos %zu\n", seconds(),
            ab->nfiles ? ab->sname : ab->fname0, ab->cpos);
}
// }}}

// some global variables (sheeh !)

struct sqrt_globals {
    int m;      // nprimes
    int n;      // degree
    int prec;
    mpz_t P;    // prime product (not to the power prec)
    size_t nbits_sqrt;
    size_t nbits_a;
    mpz_t * f_hat;
    mpz_t * f_hat_diff;
    double f_hat_coeffs;
    cado_poly pol;
    poly_t t_abpoly;
    poly_t F;
    ab_source ab;
    int lll_maxdim;
};

struct sqrt_globals glob = { .m=4, .lll_maxdim=50 };

// {{{ TODO: Now that the v field is gone, replace the polymodF layer.
// Here's the only fragments which need to remain.
    static int
poly_normalized_p (const poly_t f)
{
    return (f->deg == -1) || mpz_cmp_ui (f->coeff[f->deg], 0) != 0;
}

static void
poly_from_ab_monic(poly_t tmp, long a, unsigned long b) {
    tmp->deg = b != 0;
    mpz_set_ui (tmp->coeff[1], b);
    mpz_neg (tmp->coeff[1], tmp->coeff[1]);
    mpz_set_si (tmp->coeff[0], a);
    mpz_mul(tmp->coeff[0], tmp->coeff[0], glob.pol->f[glob.n]);
}

static void
poly_reducemodF_monic(poly_t P, poly_t p, const poly_t F)
{
    if (p->deg < F->deg) {
        poly_copy(P, p);
        return;
    }
    const int d = F->deg;
    while (p->deg >= d) {
        const int k = p->deg;
        for (int i = 0; i < d; ++i) 
            mpz_submul (p->coeff[k-d+i], p->coeff[k], F->coeff[i]);

        cleandeg (p, k-1);
    }

    poly_copy(P, p);
}

void
polymodF_mul_monic (poly_t Q, const poly_t P1, const poly_t P2,
        const poly_t F)
{
    poly_t prd;
    poly_alloc(prd, P1->deg+P2->deg);
    ASSERT_ALWAYS(poly_normalized_p (P1));
    ASSERT_ALWAYS(poly_normalized_p (P2));
    poly_mul(prd, P1, P2);
    poly_reducemodF_monic(Q, prd, F);
    poly_free(prd);
}

void poly_swap(poly_t a, poly_t b)
{
    ASSERT_ALWAYS(a->deg + 1 <= b->alloc);
    ASSERT_ALWAYS(b->deg + 1 <= a->alloc);
    for(int i = 0 ; i <= a->deg || i<=b->deg ; i++) {
        mpz_swap(a->coeff[i], b->coeff[i]);
    }
    int d = a->deg;
    a->deg = b->deg;
    b->deg = d;
}

// }}}

// {{{ floating point stuff
// {{{ getting the coefficients of the lagrange interpolation matrix.
long double complex lagrange_polynomial(long double complex * res, double * f, int deg, long double complex r)
{
    long double complex z = f[deg];
    long double complex y = deg * f[deg];
    res[deg-1] = z;
    for(int i = deg-1 ; i > 0 ; i--) {
        z *= r; z += f[i];
        y *= r; y += i * f[i];
        res[i-1] = z;
    }
    for(int i = deg-1 ; i >= 0 ; i--) {
        res[i] /= y;
    }
    return z;
}

double lagrange_polynomial_abs(double * res, double * f, int deg, long double complex r)
{
    long double complex * cres = malloc(deg * sizeof(long double complex));
    long double complex z = lagrange_polynomial(cres, f, deg, r);
    for(int i = deg-1 ; i >= 0 ; i--) {
        res[i] = cabsl(cres[i]);
    }
    free(cres);
    return cabsl(z);
}
// }}}

void estimate_nbits_sqrt(size_t * sbits, size_t * abits, ab_source_ptr ab) // , int guess)
{
    /*
    if (guess) {
        *abits = ab->digitbytes_estim * M_LN10 / M_LN2;
        *sbits = *abits / 2;
        // when doing gross estimates like this, we can hardly avoid
        // taking a safety margin.
        *abits += *abits / 10;
        *sbits += *sbits / 10;
        fprintf(stderr, "# [%2.2lf] coefficients of A"
                " have at most %zu bits (ESTIMATED)\n", seconds(), *abits);
        fprintf(stderr, "# [%2.2lf] square root coefficients"
                " have at most %zu bits (ESTIMATED)\n", seconds(), *sbits);
        return;
    }
    */

    double t1,tt;

    long double complex * eval_points = malloc(glob.n * sizeof(long double complex));
    double * double_coeffs = malloc((glob.n+1) * sizeof(double));
    double * evaluations = malloc(glob.n * sizeof(double));

    // take the roots of f, and multiply later on to obtain the roots of
    // f_hat. Otherwise we encounter precision issues.
    for(int i = 0 ; i <= glob.n ; i++) {
        double_coeffs[i] = mpz_get_d(glob.pol->f[i]);
    }

    int rc = poly_roots_longdouble(double_coeffs, glob.n, eval_points);
    if (rc) {
        fprintf(stderr, "# [%2.2lf] Warning: rootfinder had accuracy problem with %d roots\n", seconds(), rc);
    }


    // {{{ compress the list of roots.
    int nreal = 0, ncomplex = 0, rs = 0;
    for(int i = 0 ; i < glob.n ; i++) {
        if (cimagl(eval_points[i]) > 0) {
            eval_points[rs] = eval_points[i];
            rs++;
            ncomplex++;
        } else if (cimagl(eval_points[i]) < 0) {
            continue;
        } else {
            eval_points[rs] = creall(eval_points[i]);
            // eval_points[rs] = eval_points[i];
            rs++;
            nreal++;
        }
    }
    // }}}

    // {{{ post-scale to roots of f_hat, and store f_hat instead of f
    for(int i = 0 ; i < rs ; i++) {
        eval_points[i] *= mpz_get_d(glob.pol->f[glob.n]);
    }
    for(int i = 0 ; i <= glob.n ; i++) {
        double_coeffs[i] = mpz_get_d(glob.f_hat[i]);
    }
    // }}}

    // {{{ print the roots.
    if (nreal) {
        fprintf(stderr, "# [%2.2lf]", seconds());
        fprintf(stderr, " real");
        for(int i = 0 ; i < rs ; i++) {
            double r = cimagl(eval_points[i]);
            if (r == 0) {
                fprintf(stderr, " %.4Lg", creall(eval_points[i]));
            }
        }
    }
    if (ncomplex) {
        fprintf(stderr, " complex");
        for(int i = 0 ; i < rs ; i++) {
            if (cimagl(eval_points[i]) > 0) {
                fprintf(stderr, " %.4Lg+i*%.4Lg", creall(eval_points[i]), cimagl(eval_points[i]));
            }
        }
    }
    fprintf(stderr, "\n");
    // }}}

    // {{{ now evaluate the product.
    for(int i = 0 ; i < glob.n ; i++) {
        evaluations[i] = 0;
    }

    // Consider the product A of all a-b\alpha's, which is a
    // polynomial in \alpha. This polynomial has _rational_ coefficients,
    // since each reduction involves dividing out by f_d.
    //
    // Easier to handle is f_d^nab*A (product of all f_d*a-b*f_d*\alpha),
    // which is an element of the order Z[f_d alpha] (f_d
    // alpha is an algebraic integer). It can be expressed with integer
    // coefficients in the powers of f_d\alpha. Its square root, however,
    // is not necessarily in this sub-order of the ring of integers.
    // Therefore we multiply by f_hat'(f_d\alpha), where f_hat is the
    // minimal polynomial of f_d\alpha.
    //
    // We ensure that we've rounded up nab to the next even multiple.

    int64_t a;
    uint64_t b;
    t1 = seconds();
    ab_source_rewind(ab);
    for( ; ab_source_next(ab, &a, &b) ; ) {
        for(int i = 0 ; i < rs ; i++) {
            long double complex y = a * mpz_get_d(glob.pol->f[glob.n]);   
            long double complex w = eval_points[i] * b;
            y = y - w;
            evaluations[i] += log(cabsl(y));
        }
        tt = seconds();
        if (tt > t1 + print_delay || !(ab->nab % 10000000)) {
            t1 = tt;
            fprintf(stderr,
                    "# [%2.2lf] floating point evaluation: %zu (%.1f%%)\n",
                    t1, ab->nab, 100.0*(double)ab->nab/ab->nab_estim);
        }
    }
    // }}}
    // note that now that we've read everything, we know the precise
    // number of (a,b)'s. Thus we can replace the estimation.
    ab->nab_estim = ab->nab;

    // {{{ post-process evaluation: f'(alpha), and even nab. print.

    if (ab->nab & 1) {
        fprintf(stderr, "# [%2.2lf] odd number of pairs !\n", seconds());
        for(int i = 0 ; i < rs ; i++) {
            evaluations[i] += log(fabs(mpz_get_d(glob.pol->f[glob.n])));
        }
    }

    // multiply by the square of f_hat'(f_d\alpha).
    for(int i = 0 ; i < rs ; i++) {
        long complex double s = glob.n;
        for(int j = glob.n - 1 ; j >= 0 ; j--) {
            s *= eval_points[i];
            s += double_coeffs[j] * j;
        }
        evaluations[i] += 2 * clog(s);
    }
    fprintf(stderr, "# [%2.2lf] Log_2(A)", seconds());
    for(int i = 0 ; i < rs ; i++) {
        fprintf(stderr, " %.4Lg", creall(evaluations[i]) / M_LN2);
        if (cimagl(eval_points[i]) > 0) {
            fprintf(stderr, "*2");
        }
    }
    fprintf(stderr, "\n");
    // }}}

    // {{{ deduce the lognorm. print.
    double lognorm = 0;
    for(int i = 0 ; i < rs ; i++) {
        if (cimagl(eval_points[i]) > 0) { 
            lognorm += 2*creall(evaluations[i]);
        } else {
            lognorm += creall(evaluations[i]);
        }
    }
    fprintf(stderr, "# [%2.2lf] log_2(norm(A)) %.4g\n",
            seconds(), lognorm / M_LN2);
    // }}}

    // {{{ now multiply this by the Lagrange matrix.
    double * a_bounds = malloc(glob.n * sizeof(double));
    double * sqrt_bounds = malloc(glob.n * sizeof(double));

    for(int j = 0 ; j < glob.n ; j++) {
        a_bounds[j] = 0;
        sqrt_bounds[j] = 0;
    }
    for(int i = 0 ; i < rs ; i++) {
        double * lmat = malloc(glob.n * sizeof(double ));
        lagrange_polynomial_abs(lmat, double_coeffs, glob.n, eval_points[i]);
        for(int j = 0 ; j < glob.n ; j++) {
            double za, zs;
            za = zs = log(lmat[j]);
            za += evaluations[i];
            zs += evaluations[i] / 2;
            if (cimagl(eval_points[i]) > 0) {
                za += log(2);
                zs += log(2);
            }
            if (za > a_bounds[j]) a_bounds[j] = za;
            if (zs > sqrt_bounds[j]) sqrt_bounds[j] = zs;
        }
        free(lmat);
    }
    // }}}

    // {{{ get global logbounds
    double logbound_sqrt = 0;
    double logbound_a = 0;
    for(int j = 0 ; j < glob.n ; j++) {
        // note that we might have added up to n times the same thing.
        // (inequality a+b < 2max(a,b) )
        sqrt_bounds[j] += log(glob.n);
        a_bounds[j] += log(glob.n);

        // safety margin for inaccuracies ?
        sqrt_bounds[j] += 100 * M_LN2;
        a_bounds[j] += 100 * M_LN2;
        if (sqrt_bounds[j] > logbound_sqrt) logbound_sqrt = sqrt_bounds[j];
        if (a_bounds[j] > logbound_a) logbound_a = a_bounds[j];
    }
    // }}}

#if 0
    // {{{ print logbounds
    fprintf(stderr, "# [%2.2lf] logbounds per coeff of A", seconds());
    for(int j = 0 ; j < glob.n ; j++) {
        fprintf(stderr, " %.4Lg", a_bounds[j] / M_LN2);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "# [%2.2lf] logbounds per coeff of sqrt(A)", seconds());
    for(int j = 0 ; j < glob.n ; j++) {
        fprintf(stderr, " %.4Lg", sqrt_bounds[j] / M_LN2);
    }
    fprintf(stderr, "\n");
    // }}}
#endif

    // Since we're being very gross, we use a stupid bound on the last
    // coefficient.

    *sbits = ceil(logbound_sqrt / M_LN2);
    *abits = ceil(logbound_a / M_LN2);
    fprintf(stderr, "# [%2.2lf] coefficients of A"
            " have at most %zu bits\n", seconds(), *abits);
    fprintf(stderr, "# [%2.2lf] square root coefficients"
            " have at most %zu bits\n", seconds(), *sbits);

    free(a_bounds);
    free(sqrt_bounds);
    free(eval_points);
    free(double_coeffs);
    free(evaluations);
}
/* }}} */

#define ABPOLY_OFFSET_THRESHOLD        16
// NOTE: This does not depend on p (nor r of course).
// the accumulation is done for all data between:
// the first data line starting at offset >= off0 (inclusive)
// the first data line starting at offset >= off1 (exclusive)
void accumulate_ab_poly(poly_t P, ab_source_ptr ab, size_t off0, size_t off1)
{
    mpz_set_ui(P->coeff[0], 1);
    P->deg = 0;
    if (off1 - off0 < ABPOLY_OFFSET_THRESHOLD) {
        ab_source_move_afterpos(ab, off0);
        for( ; ab->tpos < off1 ; ) {
            int64_t a;
            uint64_t b;
            int r = ab_source_next(ab, &a, &b);
            FATAL_ERROR_CHECK(!r, "dep file ended prematurely\n");
            poly_from_ab_monic(glob.t_abpoly, a, b);
            polymodF_mul_monic(P, P, glob.t_abpoly, glob.F);
        }
        return;
    }
    size_t d = (off1 - off0) / 2;
    poly_t Pl, Pr;
    poly_alloc(Pl, glob.n);
    poly_alloc(Pr, glob.n);
    accumulate_ab_poly(Pl, ab, off0, off0 + d);
    accumulate_ab_poly(Pr, ab, off0 + d, off1);
    polymodF_mul_monic(P, Pl, Pr, glob.F);
    poly_free(Pl);
    poly_free(Pr);
}

/*{{{ have to pre-declare prime_data */
struct alg_ptree_s;

struct individual_contribution {
    uint64_t ratio;
    mpz_t modN;
};

struct prime_data {
    unsigned long p;
    unsigned long * r;
    // unsigned long rj;
    void * powers;      // see .cpp file.

    // computed somewhat late.
    mpz_t iHx;

    // these are temporary.
    mpz_t invdev_rx;
    mpz_t ta, tb;

    struct alg_ptree_s * T;

    // after the rational ptree reduction, this contains the share of A
    // with coefficients reduced. Short-lived.
    poly_t A;

    // this is the set of evaluations.
    poly_t evals;

    // roots.
    poly_t lroots;
    mpz_ptr rx;

    // square roots of A(x)     (only in the end !)
    poly_t sqrts;
    mpz_ptr sx;

};/* }}} */

/* {{{ product tree (rational) */
struct rat_ptree_s {
    struct rat_ptree_s * t0;
    struct rat_ptree_s * t1;
    mpz_t z;
    mpz_srcptr zx;
    struct prime_data * p;
};
typedef struct rat_ptree_s rat_ptree_t;

rat_ptree_t * rat_ptree_build(struct prime_data * p, int i0, int i1)
{
    ASSERT_ALWAYS(i0 < i1);
    rat_ptree_t * res = malloc(sizeof(rat_ptree_t));
    memset(res, 0, sizeof(rat_ptree_t));
    if (i1-i0 == 1) {
        res->p = p + i0;
        res->zx = power_lookup(res->p->powers, glob.prec);
        return res;
    }
    int d = (i1-i0)/2;
    res->t0 = rat_ptree_build(p, i0, i0+d);
    res->t1 = rat_ptree_build(p, i0+d, i1);
    mpz_init(res->z);
    res->zx = res->z;
    mpz_mul(res->z, res->t0->zx, res->t1->zx);
    return res;
}

void rat_ptree_clear(rat_ptree_t * t)
{
    if (t == NULL) return;
    rat_ptree_clear(t->t0);
    rat_ptree_clear(t->t1);
    if (t->p == NULL) mpz_clear(t->z);
    free(t);
}

#if 0
unsigned int rat_ptree_nleaves(rat_ptree_t * t)
{
    if (t->p)
        return 1;
    unsigned int n0 = rat_ptree_nleaves(t->t0);
    unsigned int n1 = rat_ptree_nleaves(t->t1);
    return n0 + n1;
}
#endif

/* }}} */

/* {{{ product tree (algebraic) */
struct alg_ptree_s {
    struct alg_ptree_s * t0;
    struct alg_ptree_s * t1;
    poly_t s;
};
typedef struct alg_ptree_s alg_ptree_t;

alg_ptree_t * alg_ptree_build(struct prime_data * p, int i0, int i1)
{
    /* Everything being done the naive way, the algebraic ptree wins
     * nothing: the count of multiplications is exactly the same (n^2-n),
     * and the algebraic ptree incurs an extra storage for n
     * coefficients. It's thus disabled.
     */
    return NULL;

    mpz_srcptr px = power_lookup(p->powers, glob.prec);
    ASSERT_ALWAYS(i0 < i1);
    alg_ptree_t * res = malloc(sizeof(alg_ptree_t));
    memset(res, 0, sizeof(alg_ptree_t));
    poly_alloc(res->s, i1-i0);
    if (i1-i0 == 1) {
        res->s->deg = 1;
        mpz_set_ui(res->s->coeff[1], 1);
        mpz_neg(res->s->coeff[0], p->rx);
        return res;
    }
    int d = (i1-i0)/2;
    res->t0 = alg_ptree_build(p, i0, i0+d);
    res->t1 = alg_ptree_build(p, i0+d, i1);
    poly_mul(res->s, res->t0->s, res->t1->s);
    poly_reducemodF_monic(res->s, res->s, glob.F);
    poly_reduce_mod_mpz (res->s, res->s, px);
    return res;
}

void alg_ptree_clear(alg_ptree_t * T)
{
    if (T == NULL) return;
    alg_ptree_clear(T->t0);
    alg_ptree_clear(T->t1);
    free(T);
}

/* }}} */

// of course this destroys P.
// we assume that P is reduced mod the top-level T.
void reduce_poly_mod_rat_ptree(poly_t P, rat_ptree_t * T)
{
    if (T->p) {
        poly_swap(T->p->A, P);
        return;
    }
    poly_t temp;
    poly_alloc(temp, glob.n);
    temp->deg = glob.n - 1;
    for(int i = 0 ; i < glob.n ; i++) {
        mpz_mod(temp->coeff[i], P->coeff[i], T->t0->zx);
        mpz_mod(P->coeff[i], P->coeff[i], T->t1->zx);
    }
    cleandeg(temp, glob.n - 1);
    cleandeg(P, glob.n - 1);
    reduce_poly_mod_rat_ptree(temp, T->t0);
    reduce_poly_mod_rat_ptree(P, T->t1);
    poly_free(temp);
}

void reduce_poly_mod_alg_ptree(alg_ptree_t * T, struct prime_data * p)
{
    FATAL_ERROR_CHECK(T != NULL, "not implemented");
    mpz_srcptr px = power_lookup(p->powers, glob.prec);
    for(int j = 0 ; j < glob.n ; j++) {
        mpz_set(p->ta, p->A->coeff[glob.n - 1]);
        for(int k = glob.n - 2 ; k >= 0 ; k--) {
            mpz_mul(p->ta, p->ta, p->lroots->coeff[j]);
            mpz_add(p->ta, p->ta, p->A->coeff[k]);
            mpz_mod(p->ta, p->ta, px);
        }
        mpz_mul(p->evals->coeff[j], p->evals->coeff[j], p->ta);
        if (mpz_size(p->evals->coeff[j]) >= mpz_size(px) * 3/2)
            mpz_mod(p->evals->coeff[j], p->evals->coeff[j], px);
    }
    for(int k = glob.n - 1 ; k >= 0 ; k--) {
        mpz_realloc(p->A->coeff[k], 0);
    }
    mpz_realloc(p->ta, 0);
}

/* {{{ finding the CRT primes */

typedef int (*sortfunc_t) (const void *, const void *);
int modul_cmp_sortfunc(const unsigned long * a, const unsigned long * b)
{
    return (*a>*b) - (*b>*a);
}

struct prime_data * suitable_crt_primes()
{
    unsigned int m = glob.m;
    struct prime_data * res = malloc(m * sizeof(struct prime_data));
    unsigned int i = 0;

    // Note that p0 has to be well above the factor base bound. Thus on
    // 32-bit machines, it's possibly problematic.
    unsigned long p0 = ((~0UL)>>1)+1;    // 2^31 or 2^63
    unsigned long p = p0;
    unsigned long * roots = malloc(glob.n * sizeof(unsigned long));

    fprintf(stderr, "# [%2.2lf] Searching for CRT primes\n", seconds());
    // fprintf(stderr, "# [%2.2lf] p0=%lu\n", seconds(), p0);

    for( ; i < m ; ) {
        p = ulong_nextprime(p);
        modulusul_t q;
        modul_initmod_ul(q, p);
        memset(roots, 0, glob.n * sizeof(unsigned long));
        int nr = modul_poly_roots_ulong(roots, glob.pol->f, glob.n, q);
        if (nr != glob.n) continue;
        memset(&(res[i]), 0, sizeof(struct prime_data));
        res[i].r = malloc(glob.n * sizeof(unsigned long));
        memcpy(res[i].r, roots, glob.n * sizeof(unsigned long));
        res[i].p = p;
        // res[i].log2_p = log(p)/M_LN2;
        // fprintf(stderr, "# [%2.2lf] p0+%lu", seconds(), p-p0);
        residueul_t fd;
        modul_init(fd, q);
        modul_set_ul_reduced(fd, mpz_fdiv_ui(glob.pol->f[glob.n], modul_getmod_ul (q)), q);
        for(int j = 0 ; j < glob.n ; j++) {
            residueul_t r;
            modul_init(r, q);
            modul_set_ul_reduced(r, res[i].r[j], q);
            modul_mul(r,r,fd,q);
            res[i].r[j] = modul_get_ul(r, q);
            // fprintf(stderr, " %lu", res[i].r[j]);
            modul_clear(r, q);
        }
        // sort with respect to the values of the scaled roots.
        qsort(res[i].r, glob.n, sizeof(unsigned long), (sortfunc_t) &modul_cmp_sortfunc);
        modul_clear(fd, q);
        modul_clearmod(q);
        // fprintf(stderr, "\n");
        i++;
    }
    fprintf(stderr, "# [%2.2lf] Found all CRT primes\n", seconds());
    free(roots);
    getprime(0);

    return res;
}
/* }}} */

/* {{{ everything that happens only modulo one prime */

void inversion_lift(struct prime_data * p, mpz_t Hx, int precision)/* {{{ */
{
    double t0 = seconds();

    mpz_srcptr pk = power_lookup(p->powers, precision);
    assert(precision > 0);

    if (precision == 1) {
        mpz_invert(p->iHx, Hx, pk);
        return;
    }
    int lower = precision - precision / 2;

    mpz_srcptr pl = power_lookup(p->powers, lower);

    // we're going to recurse. change the long Hx_mod value stored by a
    // temporary small one. The problem with this approach is that we
    // store many reductions in memory.
    mpz_t Hx_save;
    mpz_init(Hx_save);
    mpz_mod(Hx_save, Hx, pl);
    // recurse.
    inversion_lift(p, Hx_save, lower);
    mpz_clear(Hx_save);

    if (seconds() > t0 + print_delay)
        fprintf(stderr, "# [%2.2lf] precision %d\n", seconds(), precision);

    mpz_mul(p->ta, p->iHx, Hx);
    mpz_mod(p->ta, p->ta, pk);

    mpz_sub_ui(p->ta, p->ta, 1);
    mpz_mul(p->ta, p->ta, p->iHx);
    mpz_sub(p->iHx, p->iHx, p->ta);
    mpz_mod(p->iHx, p->iHx, pk);
    // gmp_fprintf(stderr, "# [%2.2lf] %Zd\n", seconds(), p->iHx_mod);
}/* }}} */

/* {{{ tonelli-shanks */
void modul_find_ts_gen(residueul_t z, modulusul_t p)
{
    unsigned long pp = modul_getmod_ul(p)-1;
    int e = ctzl(pp);
    pp >>= e;
    unsigned long s = 1UL << (e-1);
    residueul_t r;
    modul_init(r, p);
    do {
        modul_set_ul(z, random(), p);
        modul_pow_ul(z, z, pp, p);
        modul_pow_ul(r, z, s, p);
        modul_add_ul(r, r, 1, p);
    } while (!modul_is0(r, p));
    modul_clear(r, p);
}

int modul_field_sqrt(residueul_t z, residueul_t a, residueul_t g, modulusul_t p)
{
    unsigned long pp = modul_getmod_ul(p)-1;
    int e = ctzl(pp);
    pp >>= e+1;
    if (modul_is0(a, p)) {
        modul_set0(z, p);
        return 1;
    }
    if (modul_is0(g, p)) {
        modul_find_ts_gen(g, p);
    }
    residueul_t b, x, y, t;
    modul_init(b, p);
    modul_init(x, p);
    modul_init(y, p);
    modul_init(t, p);
    int r = e;
    unsigned long s = 1UL << (e-1);

    // modul_set(x, a, p);
    modul_set(y, g, p);

    modul_pow_ul(x, a, pp, p);
    modul_sqr(b, x, p);
    modul_mul(x, x, a, p);
    modul_mul(b, b, a, p);

    int m;
    for(;;) {
        modul_set(t, b, p);
        for(m=0; !modul_is1(t, p); m++)
            modul_sqr(t, t, p);
        assert(m<=r);

        if (m==0 || m==r)
            break;

        s = 1UL << (r-m-1);
        r = m;

        modul_pow_ul(t, y, s, p);
        modul_sqr(y, t, p);
        modul_mul(x, x, t, p);
        modul_mul(b, b, y, p);
    }

    modul_set(z, x, p);
    modul_clear(t, p);
    modul_clear(x, p);
    modul_clear(y, p);
    modul_clear(b, p);
    return (m==0);
}
/* }}} */
/* {{{ sqrt / invsqrt */
void invsqrt_lift(struct prime_data * p, mpz_t A, int precision)
{
    double t0 = seconds();
    assert(precision > 0);

    if (precision == 1) {
        residueul_t z, a;
        modulusul_t q;
        modul_initmod_ul(q, p->p);
        modul_init(z, q);
        modul_init(a, q);
        modul_set_ul_reduced(a, mpz_get_ui(A), q);
        // fprintf(stderr, "A mod (%lu, alpha-r)=%lu\n", p->p, modul_get_ul(a, q));
        int issquare = modul_field_sqrt(a, a, z, q);
        ASSERT_ALWAYS(issquare);
        modul_inv(a, a, q);
        if (modul_get_ul(a, q) & 1) modul_neg(a, a, q);
        mpz_set_ui(p->sx, modul_get_ul(a, q));
        modul_clear(a, q);
        modul_clear(z, q);
        modul_clearmod(q);
        return;
    }
    int lower = precision - precision / 2;

    mpz_srcptr pl = power_lookup(p->powers, lower);

    // we're going to recurse ; arrange so that we recurse on something
    // acceptably small.  The problem with this approach is that we store
    // many reductions in memory.
    mpz_t A_save;
    mpz_init(A_save);
    mpz_mod(A_save, A, pl);
    // recurse.
    invsqrt_lift(p, A_save, lower);
    mpz_clear(A_save);

    if (seconds() > t0 + print_delay)
        fprintf(stderr, "# [%2.2lf] precision %d\n", seconds(), precision);

    mpz_srcptr pk = power_lookup(p->powers, precision);

    mpz_mul(p->ta, p->sx, p->sx);
    mpz_mod(p->ta, p->ta, pk);
    mpz_mul(p->ta, p->ta, A);
    mpz_mod(p->ta, p->ta, pk);
    mpz_sub_ui(p->ta, p->ta, 1);
    if (mpz_odd_p(p->ta)) mpz_add(p->ta, p->ta, pk);
    mpz_div_2exp(p->ta, p->ta, 1);
    mpz_submul(p->sx, p->ta, p->sx);
    mpz_mod(p->sx, p->sx, pk);
    // gmp_fprintf(stderr, "# [%2.2lf] %Zd\n", seconds(), p->iH);
}

void sqrt_lift(struct prime_data * p, mpz_t A, int precision)
{
    double t0 = seconds();
    int lower = precision - precision / 2;
    mpz_srcptr pk = power_lookup(p->powers, precision);
    mpz_srcptr pl = power_lookup(p->powers, lower);
    mpz_t A_save;
    mpz_init(A_save);
    mpz_mod(A_save, A, pl);
    invsqrt_lift(p, A_save, lower);
    mpz_clear(A_save);

    if (seconds() > t0 + print_delay)
        fprintf(stderr, "# [%2.2lf] precision %d\n", seconds(), precision);
    // inverse square root now in sx.
    mpz_mul(p->ta, A, p->sx);
    mpz_mod(p->ta, p->ta, pl);
    // XXX This destroys A !!!
    mpz_submul(A, p->ta, p->ta);
    if (mpz_odd_p(p->sx)) mpz_add(p->sx, p->sx, pk);
    mpz_div_2exp(p->sx, p->sx, 1);
    mpz_addmul(p->ta, p->sx, A);
    mpz_mod(p->sx, p->ta, pk);
}
/* }}} */

/* {{{ lifting the roots */
// assume that q is very considerably larger than all coefficients of f
// -- this will be the case for the iterations that matter.
static void mp_poly_eval_mod(mpz_t r, mpz_t * poly, int deg, mpz_srcptr a, mpz_srcptr q)
{
    int i;

    mpz_set(r, poly[deg]);
    for (i = deg - 1; i >= 0; i--) {
        mpz_mul(r, r, a);
        mpz_mod(r, r, q);
        mpz_add(r, r, poly[i]);
    }
}
/* This implements the following iteration */
/*
p:=goodprimes[1];
r1:=GF(p)!goodprimes_lroots[1][1];
z1:=(1/Evaluate(Derivative(PolynomialRing(GF(p))!f),r1));
r:=GF(p)!r1;
z:=GF(p)!z1;
k:=1;
while k lt lift_prec do
k*:=2;
R:=Integers(p^k);
r:=R!Z!r;
z:=R!Z!z;
fr:=Evaluate(PolynomialRing(R)!f, r);
r:=r-fr*z;
fdr:=Evaluate(Derivative(PolynomialRing(R)!f), r);
z:=z-z*(fdr*z-1);
end while;
*/
void root_lift_innerlevels(struct prime_data * p, int precision)/* {{{ */
{
    double t0 = seconds();
    assert(precision > 0);

    if (precision == 1)
        return;
    int lower = precision - precision / 2;

    // recurse.
    root_lift_innerlevels(p, lower);

    if (seconds() > t0 + print_delay)
        fprintf(stderr, "# [%2.2lf] precision %d\n", seconds(), precision);

    mpz_srcptr pk = power_lookup(p->powers, precision);

    // we know r to half-precision, 1/f'(r) to half-precision.
    mp_poly_eval_mod(p->tb, glob.f_hat, glob.n, p->rx, pk);
    mpz_mul(p->ta, p->tb, p->invdev_rx);
    // printf("%zu %zu %zu %zu\n", mpz_size(pk), mpz_size(p->ta), mpz_size(p->tb), mpz_size(p->invdev_rx));
    mpz_sub(p->rx, p->rx, p->ta);
    mpz_mod(p->rx, p->rx, pk);

    mp_poly_eval_mod(p->tb, glob.f_hat_diff, glob.n-1, p->rx, pk);
    mpz_mul(p->ta, p->invdev_rx, p->tb);
    mpz_sub_ui(p->ta, p->ta, 1);
    mpz_mod(p->ta, p->ta, pk);
    mpz_mul(p->tb, p->ta, p->invdev_rx);
    mpz_sub(p->ta, p->invdev_rx, p->tb);
    mpz_mod(p->invdev_rx, p->ta, pk);
    /*
       if (precision < 40) {
       gmp_fprintf(stderr, "# [%2.2lf] %Zd\n", seconds(), p->rx);
       gmp_fprintf(stderr, "# [%2.2lf] %Zd\n", seconds(), p->invdev_rx);
       }
       */
}
void root_lift(struct prime_data * p, int precision)
{
    double t0 = seconds();
    assert(precision > 0);

    if (precision == 1)
        return;
    int lower = precision - precision / 2;

    // recurse.
    root_lift_innerlevels(p, lower);

    if (seconds() > t0 + print_delay)
        fprintf(stderr, "# [%2.2lf] precision %d\n", seconds(), precision);

    mpz_srcptr pk = power_lookup(p->powers, precision);

    // we know r to half-precision, 1/f'(r) to half-precision.
    mp_poly_eval_mod(p->tb, glob.f_hat, glob.n, p->rx, pk);
    mpz_mul(p->ta, p->tb, p->invdev_rx);
    // printf("%zu %zu %zu %zu\n", mpz_size(pk), mpz_size(p->ta), mpz_size(p->tb), mpz_size(p->invdev_rx));
    mpz_sub(p->rx, p->rx, p->ta);
    mpz_mod(p->rx, p->rx, pk);
    // don't do the lower part.
}/* }}} */
/* }}} */

void prime_initialization(struct prime_data * p)/* {{{ */
{
    // trigger computation of p^glob.prec
    p->powers = power_lookup_table_init(p->p);
    mpz_init(p->iHx);
    p->rx = NULL;
    p->sx = NULL;

    mpz_init(p->invdev_rx);
    mpz_init(p->ta);
    mpz_init(p->tb);

    poly_alloc(p->A, glob.n-1);
    poly_alloc(p->evals, glob.n-1);
    poly_alloc(p->lroots, glob.n-1);
    poly_alloc(p->sqrts, glob.n-1);
}
/* }}} */

void prime_precomputations(struct prime_data * p)/* {{{ */
{
    /*
    mpz_srcptr px = power_lookup(p->powers, glob.prec);

    fprintf(stderr, "# [%2.2lf] considering p=%lu\n",
            seconds(), p->p);
    fprintf(stderr, 
            "# [%2.2lf] size(p^l) %.1f MB\n",
            seconds(), mpz_sizeinbase(px, 256) * 1.0e-6);
            */

    for(int j = 0 ; j < glob.n ; j++) {
        mpz_srcptr p1 = power_lookup(p->powers, 1);
        // mpz_srcptr px = power_lookup(p->powers, glob.prec);
        // fprintf(stderr, "# [%2.2lf] considering r=%lu\n", seconds(), p->r[j]);

        // lift the root.
        // p->rj = p->r[j];
        mpz_set_ui(p->lroots->coeff[j], p->r[j]);
        p->rx = p->lroots->coeff[j];
        mp_poly_eval_mod(p->invdev_rx, glob.f_hat_diff, glob.n-1, p->rx, p1);
        mpz_invert(p->invdev_rx, p->invdev_rx, p1);
        fprintf(stderr, "# [%2.2lf] lifting p=%lu, r=%lu\n",
                seconds(), p->p, p->r[j]);
        root_lift(p, glob.prec);
        // fprintf(stderr, "# [%2.2lf] done\n", seconds());
        mpz_realloc(p->ta, 0);
        mpz_realloc(p->tb, 0);
        mpz_realloc(p->invdev_rx, 0);
    }
    p->T = alg_ptree_build(p, 0, glob.n);

}
/* }}} */

void prime_postcomputations(struct individual_contribution * cont, struct prime_data * p)/* {{{ */
{
    mpz_srcptr px = power_lookup(p->powers, glob.prec);

    mpz_t Hx;
    mpz_init(Hx);
    mpz_divexact_ui(Hx, glob.P, p->p);
    mpz_powm_ui(Hx, Hx, glob.prec, px);
    /*
    fprintf(stderr, "# [%2.2lf] size(H^l) %.1f MB\n",
            seconds(), mpz_sizeinbase(p->Hx, 256) * 1.0e-6);
            */
    // compute the inverse of H^prec modulo p^prec.
    // need a recursive function for computing the inverse.
    fprintf(stderr, "# [%2.2lf] lifting H^-l\n", seconds());
    inversion_lift(p, Hx, glob.prec);
    mpz_clear(Hx);
    // fprintf(stderr, "# [%2.2lf] done\n", seconds());
    mpz_realloc(p->ta, 0);
    mpz_realloc(p->tb, 0);
    
    // Lagrange reconstruction.
    //
    // Normally each coefficient has to be divided by the evaluation of
    // the derivative. However we skip this division, effectively
    // reconstructing the polynomial multiplied by the square of the
    // derivative -- which is exactly what we're looking for, in fact.

    // recall that we're working with the number field sieve in mind. So
    // we don't really care about the whole reconstruction in the number
    // field, and from here on we are going to take wild shortcuts.
    // Indeed, even though the ``magical sign combination'' is not known
    // at this point, we do know that the eventual reconstruction will be
    // linear. Thus instead of storing n^2 full length modular integers
    // (n for each root),  and do this for each prime, we store only the
    // pair (quotient mod p^x, residue mod N).

    mpf_t pxf, ratio;
    mpz_t z;

    mpf_init2(pxf, 256);
    mpf_init2(ratio, 256);
    mpz_init(z);

    mpf_set_z(pxf, px);

    mpz_t Hxm;
    mpz_init(Hxm);
    mpz_set_ui(Hxm, p->p);
    mpz_invert(Hxm, Hxm, glob.pol->n);
    mpz_mul(Hxm, Hxm, glob.P);
    mpz_mod(Hxm, Hxm, glob.pol->n);
    mpz_powm_ui(Hxm, Hxm, glob.prec, glob.pol->n);

    for(int j = 0 ; j < glob.n ; j++) {
        p->rx = p->lroots->coeff[j];
        p->sx = p->sqrts->coeff[j];

        // so we have this nice square root. The first thing we do on our
        // list is to scramble it by multiplying it with the inverse of
        // H^x...
        mpz_mul(p->sx, p->sx, p->iHx);
        mpz_mod(p->sx, p->sx, px);

        // Now use the evaluation of f_hat mod rx to obtain the lagrange
        // coefficients.
        mpz_set_ui(p->ta, 1);
        for(int k = glob.n - 1 ; k >= 0 ; k--) {
            if (k < glob.n - 1) {
                mpz_mul(p->ta, p->ta, p->rx);
                mpz_add(p->ta, p->ta, glob.f_hat[k+1]);
                mpz_mod(p->ta, p->ta, px);
            }
            // multiply directly with H^-x * sqrt
            mpz_mul(p->tb, p->ta, p->sx);
            if (k < glob.n - 1) {
                mpz_mod(p->tb, p->tb, px);
            }
            ASSERT_ALWAYS(mpz_cmp_ui(p->tb, 0) >= 0);
            ASSERT_ALWAYS(mpz_cmp(p->tb, px) < 0);

            // now the shortcuts.
            mpf_set_z(ratio, p->tb);
            mpf_div(ratio, ratio, pxf);
            mpf_mul_2exp(ratio, ratio, 64);
            mpz_set_f(z, ratio);

#if GMP_LIMB_BITS == 64
            cont[j*glob.n+k].ratio = mpz_get_ui(z);
#else
            cont[j*glob.n+k].ratio   = (uint64_t) mpz_getlimb(z,1);
            cont[j*glob.n+k].ratio <<= 32;
            cont[j*glob.n+k].ratio  |= (uint64_t) mpz_getlimb(z,0);
#endif
            mpz_mul(p->tb, p->tb, Hxm);
            mpz_mod(cont[j*glob.n+k].modN, p->tb, glob.pol->n);
        }
    }
    mpz_clear(Hxm);
    mpz_clear(z);
    mpf_clear(pxf);
    mpf_clear(ratio);
}
/* }}} */

void prime_sqrt(struct prime_data * p)/* {{{ */
{
    for(int j = 0 ; j < glob.n ; j++) {
        p->sx = p->sqrts->coeff[j];
        fprintf(stderr, "# [%2.2lf] lifting sqrt(A) mod (%lu, alpha-%lu)\n", seconds(), p->p, p->r[j]);
        sqrt_lift(p, p->evals->coeff[j], glob.prec);
        // fprintf(stderr, "# [%2.2lf] done\n", seconds());
        mpz_realloc(p->evals->coeff[j], 0);
    }
    mpz_realloc(p->ta, 0);
    mpz_realloc(p->tb, 0);
}/* }}} */

void prime_cleanup(struct prime_data * p)/* {{{ */
{
    power_lookup_table_clear(p->powers);
    mpz_clear(p->iHx);
    mpz_clear(p->ta);
    mpz_clear(p->tb);

    mpz_clear(p->invdev_rx);

    free(p->r);
    alg_ptree_clear(p->T);

    poly_free(p->A);
    poly_free(p->evals);
    poly_free(p->lroots);
    poly_free(p->sqrts);
    // pointers p->rx and p->sx untouched of course.

    memset(p, 0, sizeof(struct prime_data));
}
/* }}} */

/* }}} */

struct kns_sum {
    uint64_t x;
    uint64_t v;
};

int kns_sum_cmp(const struct kns_sum * s, const struct kns_sum * t)
{
    int d = (s->x > t->x) - (s->x < t->x);
    if (d) return d;
    else return (s->v > t->v) - (s->v < t->v);
}

struct kns_sum * all_sums(const struct individual_contribution * t, unsigned int t_stride, unsigned long k, size_t extra_alloc, uint64_t offset)
{
    struct kns_sum * r = malloc((extra_alloc + (1UL << k)) * sizeof(struct kns_sum));
    uint64_t * t2 = malloc(k * sizeof(uint64_t));
    uint64_t x = offset;
    for(uint64_t i = 0 ; i < k ; i++, t += t_stride) {
        x -= t->ratio;
        t2[i] = t->ratio * 2;
    }
    r[0].x=x;
    uint64_t v = 0;
    r[0].v=0;
    for(uint64_t i = 1 ; i < (1UL << k) ; i++) {
        unsigned int s = __builtin_ffs(i)-1;
        uint64_t sh = 1UL << s;
        unsigned int down = v & sh;
        v ^= sh;
        x += down ? -t2[s] : t2[s];
        r[i].x = x; r[i].v = v;
    }
    free(t2);

    return r;
}

#define KNAPSACK_ALLOC_EXTRA     1024

void solve_knapsack(mpz_t sqrt_modN, struct individual_contribution * contribs, size_t lc_exp)
{
    fprintf(stderr, "# [%2.2lf] beginning recombination\n", seconds());

    /* our array has length glob.m * glob.n * glob.n */
    uint64_t bound = 1+ceil(log(glob.m*glob.n)/log(2));
    /* purpose: find a combination of the 64-bit integers in tab[], with
     * coefficients -1 or +1, which is within the interval
     * [-bound..bound].
     */

    unsigned int nelems = glob.m * glob.n;
    unsigned int k1 = nelems / 2;
    unsigned int k2 = nelems - k1;
    struct kns_sum * s1, * s2;
    /* set a stride value. */
    s1 = all_sums(contribs, glob.n, k1, KNAPSACK_ALLOC_EXTRA, bound);
    s2 = all_sums(contribs + k1 * glob.n, glob.n, k2, 0, 0);

    uint64_t n1 = (1UL << k1);
    uint64_t n2 = (1UL << k2);

    qsort(s1, n1, sizeof(struct kns_sum), (sortfunc_t) & kns_sum_cmp);
    qsort(s2, n2, sizeof(struct kns_sum), (sortfunc_t) & kns_sum_cmp);

    uint64_t ebound = 2 * bound;
    // how many elements in the head of s1 are below ebound-s2[0] ?
    unsigned int wrap = 0;
    for( ; s1[wrap].x + s2[0].x < ebound ; wrap++);
    wrap++;
    assert(wrap < KNAPSACK_ALLOC_EXTRA);
    memcpy(s1 + n1, s1, wrap * sizeof(struct kns_sum));

    unsigned int u1 = n1 + wrap;
    unsigned int pd = 0;
    assert(s2[0].x + s1[u1-1].x >= ebound);
    char * signs = malloc(nelems+1);
    memset(signs, 0, nelems+1);

    mpz_t Px;
    mpz_init(Px);
    mpz_powm_ui(Px, glob.P, glob.prec, glob.pol->n);

    mpz_t lcx;
    mpz_init(lcx);
    mpz_set(lcx, glob.pol->f[glob.n]);
    mpz_powm_ui(lcx, glob.pol->f[glob.n], lc_exp, glob.pol->n);
    mpz_invert(lcx, lcx, glob.pol->n);

    mpz_t fhdiff_modN;
    mpz_init(fhdiff_modN);
    // evaluate the derivative of f_hat in alpha_hat mod N, that is lc*m.
    mpz_set_ui(fhdiff_modN, 0);
    for(int k = glob.n ; k > 0 ; k--) {
        mpz_mul(fhdiff_modN, fhdiff_modN, glob.pol->m);
        mpz_mul(fhdiff_modN, fhdiff_modN, glob.pol->f[glob.n]);
        mpz_addmul_ui(fhdiff_modN, glob.f_hat[k], k);
        mpz_mod(fhdiff_modN, fhdiff_modN, glob.pol->n);
    }
    mpz_invert(fhdiff_modN, fhdiff_modN, glob.pol->n);

    for(unsigned int u2 = 0 ; u2 < n2 ; u2++) {
        unsigned int k;
        for( ; u1 ; u1--, pd++) {
            if ((int64_t) (s2[u2].x + s1[u1-1].x) < 0)
                break;
        }
        for( ; pd && s2[u2].x + s1[u1+pd-1].x >= ebound ; pd--);
        for(k = 0; k < pd ; k++) {
            unsigned int s;
            uint64_t v = s2[u2].v << k1 | s1[u1+k].v;
            uint64_t x = s2[u2].x + s1[u1+k].x - bound;
            for(s = 0 ; s < nelems ; s++) {
                signs[s]=(v & (1UL << s)) ? '+' : '-';
                /*
                fprintf(stderr, "%c", signs[s]);
                for(int k = 0 ; k < glob.n ; k++) {
                    fprintf(stderr, " %"PRIu64, contribs[s*glob.n+k].ratio);
                    // gmp_fprintf(stderr, "%Zd %c\n", contribs[s*glob.n+k].modN, signs[s]);
                }
                fprintf(stderr, "\n");
                */
            }
            // now try to look at this solution more closely */
            mpz_t e;
            mpz_init_set_ui(e, 0);
            int spurious = 0;
            for(int k = glob.n - 1 ; k >= 0 ; k--) {
                mpz_t z, t, zN, qz, rz;
                mpz_init_set_ui(z, bound);
                mpz_init_set_ui(zN, 0);
                mpz_init_set_ui(t, 0);
                mpz_init(qz);
                mpz_init(rz);
                uint64_t sk = bound;
                for(int j = 0 ; j < glob.m * glob.n ; j++) {
                    mpz_set_uint64(t, contribs[j * glob.n + k].ratio);
                    if (v & (1UL << j)) {
                        sk+=contribs[j * glob.n + k].ratio;
                        mpz_add(z, z, t);
                        mpz_add(zN, zN, contribs[j * glob.n + k].modN);
                    } else {
                        sk-=contribs[j * glob.n + k].ratio;
                        mpz_sub(z, z, t);
                        mpz_sub(zN, zN, contribs[j * glob.n + k].modN);
                    }
                }
                if (sk >= 2 * bound) {
                    fprintf(stderr, "recombination of coeff in X^%d yields noise (%"PRIu64")\n",k, sk);
                    spurious++;
                    // break;
                }
                // so we have (sum of r's) mod p^k = something * p^k + small
                // the ``something'' is in the quotient of the division.
                // FIXME
                // The problem is that the ``small'' thing might be small
                // enough that we won't be able to tell apart s and p-s.
                // This implies that we will have to try 2^degree
                // combinations: those with the quotients as given, and
                // the other combinations with one or several quotients
                // lowered by one unit.
                mpz_set(qz,z);
                if (mpz_cmp_ui(qz,0) >= 0) {
                    mpz_fdiv_q_2exp(qz, qz, 63);
                    mpz_add_ui(qz,qz,mpz_odd_p(qz) != 0);
                    mpz_fdiv_q_2exp(qz, qz, 1);
                } else {
                    mpz_neg(qz,qz);
                    mpz_fdiv_q_2exp(qz, qz, 63);
                    mpz_add_ui(qz,qz,mpz_odd_p(qz) != 0);
                    mpz_fdiv_q_2exp(qz, qz, 1);
                    mpz_neg(qz,qz);
                }
                mpz_mul_2exp(rz, qz, 64);
                mpz_sub(rz, z, rz);
                // gmp_printf("[%d] %Zd = %Zd * 2^64 + %Zd\n", k, z, qz, rz);
                mpz_submul(zN, qz, Px);
                mpz_mod(zN, zN, glob.pol->n);
                // gmp_printf("[X^%d] %Zd\n", k, zN);
                // good. we have the coefficient !
                mpz_mul(e, e, glob.pol->m);
                mpz_mul(e, e, glob.pol->f[glob.n]);
                mpz_add(e, e, zN);
                mpz_mod(e, e, glob.pol->n);
                mpz_clear(z);
                mpz_clear(t);
                mpz_clear(zN);
                mpz_clear(qz);
                mpz_clear(rz);

            }

            mpz_mul(e,e,lcx);
            mpz_mod(e,e,glob.pol->n);

            mpz_mul(e,e,fhdiff_modN);
            mpz_mod(e,e,glob.pol->n);

            mpz_set(sqrt_modN, e);

            mpz_clear(e);
            if (spurious) {
                fprintf(stderr, "# [%2.2lf] %"PRIx64" (%s) %"PRId64" [SPURIOUS]\n",
                        seconds(), v, signs, (int64_t) x);
            } else {
                gmp_fprintf(stderr, "# [%2.2lf] %"PRIx64" (%s) %"PRId64" [%Zd]\n",
                        seconds(), v, signs, (int64_t) x, sqrt_modN);
            }
        }
    }
    mpz_clear(Px);
    mpz_clear(lcx);
    mpz_clear(fhdiff_modN);
    free(signs);
    free(s1);
    free(s2);
    fprintf(stderr, "# [%2.2lf] recombination finished \n", seconds());

}

int main(int argc, char **argv)
{
    int ret, i;
    // int size_guess = 0;

    /* {{{ parameter parsing */
    /* print the command line */
    fprintf(stderr, "%s.r%s", argv[0], CADO_REV);
    for (i = 1; i < argc; i++)
        fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n");

    param_list pl;
    param_list_init(pl);
    param_list_configure_knob(pl, "-v", &verbose);
    // param_list_configure_knob(pl, "--size-guess", &size_guess);
    int wild = 0;
    argv++, argc--;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        if (wild == 0) {
            param_list_add_key(pl, "depfile", argv[0],
                    PARAMETER_FROM_CMDLINE);
            wild++;
            argv++, argc--;
            continue;
        }
        if (wild == 1) {
            param_list_add_key(pl, "ratdepfile", argv[0],
                    PARAMETER_FROM_CMDLINE);
            wild++;
            argv++, argc--;
            continue;
        }
        if (wild == 2) {
            param_list_add_key(pl, "polyfile", argv[0],
                    PARAMETER_FROM_CMDLINE);
            wild++;
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage();
    }

    param_list_parse_double(pl, "print_delay", &print_delay);
    param_list_parse_double(pl, "ram", &ram_gb);

    if (param_list_lookup_string(pl, "depfile") == NULL)
        usage();
    if (param_list_lookup_string(pl, "ratdepfile") == NULL)
        usage();
    if (param_list_lookup_string(pl, "polyfile") == NULL)
        usage();
    /* }}} */

    cado_poly_init(glob.pol);
    ret = cado_poly_read(glob.pol, param_list_lookup_string(pl, "polyfile"));
    glob.n = glob.pol->degree;
    ASSERT_ALWAYS(ret == 1);

    /* {{{ create f_hat, the minimal polynomial of alpha_hat = lc(f) *
     * alpha */
    {
        mpz_t tmp;
        mpz_init_set_ui(tmp, 1);

        glob.f_hat = malloc((glob.n + 1) * sizeof(mpz_t));
        glob.f_hat_diff = malloc(glob.n * sizeof(mpz_t));
        mpz_init_set_ui(glob.f_hat[glob.n], 1);
        for(int i = glob.n - 1 ; i >= 0 ; i--) {
            mpz_init(glob.f_hat[i]);
            mpz_init(glob.f_hat_diff[i]);
            mpz_mul(glob.f_hat[i], tmp, glob.pol->f[i]);
            mpz_mul(tmp, tmp, glob.pol->f[glob.n]);
            mpz_mul_ui(glob.f_hat_diff[i], glob.f_hat[i+1], i+1);
        }
        mpz_clear(tmp);
        fprintf(stderr, "# [%2.2lf] Note: all computations done with polynomial f_hat\n", seconds());
    }
    /* }}} */

    // FIXME: merge these two instances of the algebraic polynomial in
    // one type only.
    poly_alloc(glob.F, glob.n);
    for (int i = glob.pol->degree; i >= 0; --i)
        poly_setcoeff(glob.F, i, glob.f_hat[i]);


    // fprintf(stderr, "# [%2.2lf] A is f_d^%zu*f_hat'(alpha_hat)*prod(f_d a - b alpha_hat)\n", seconds(), nab + (nab &1));

    ab_source_init(glob.ab, param_list_lookup_string(pl, "depfile"));

    // note that for rsa768, this estimation takes only 10 minutes, so
    // it's not a big trouble.
    estimate_nbits_sqrt(&glob.nbits_sqrt, &glob.nbits_a, glob.ab); //, size_guess);

    // how many parts of A shall we consider ?
    // now how do we choose the appropriate number of primes ? Let:
    //
    // R: the number of bits of the ram guide amount given.
    // n: the degree
    // T: the number of bits of one coefficient of the square root.
    //
    // Note thus that A occupies a memory space of 2Tn bits.
    //
    // and
    //
    // s: the number of parts of A considered
    // B: the number of bits of lifted primes p^lambda
    // r: the number of primes stored together in a product tree.
    // t: the number of product trees considered.
    //
    // We have:
    //
    // [1] sR >= 2Tn (parts of A must form the totality of A)
    // [2] nrB <= R  (sets of ``reduced'' A's must fit in RAM)
    // [3] rtB >= T  (lifted primes must be large enough for the reconstruction)
    //
    // A further condition is given by the maximum knapsack lattice
    // dimension which can be efficiently handled. This is somewhat
    // arbitrary, let's say:
    //
    // [4] rtn <= 80
    //
    // Note that [1] readily gives s, and that [2+3] give t. Then per
    // [2,3], the product rB is placed within an interval [T/t, R/n].
    // Obeying [3], either by choosing the largest admissible power of two
    // for r, or simply the largest admissible value, we deduce r and B.

    int n,s,t,r;
    {
        size_t R = ram_gb * (1UL << 30) * 8UL;
        n = glob.n;
        size_t T = glob.nbits_sqrt + 80;    // some margin for knapsack.
        double ramfraction = 2*T*n / (double) R;
        s = ceil(ramfraction);
        t = ceil((double)n*T / R);
        int dmax = glob.lll_maxdim - n;
        double rmax = (double) dmax / n / t;
        r = rmax;
        for(int mask = ~0 ; r & mask ; mask <<= 1) r &= mask;
        size_t B = ceil((double)T/t/r);

        fprintf(stderr, "# [%2.2lf] A is %.1f%% of %.1fGB."
                " Splitting in %d parts\n",
                seconds(), 100.0 * ramfraction, ram_gb, s);

        fprintf(stderr, "# [%2.2lf] %d groups of %d prime powers,"
                " each of %zu bits\n",
                seconds(), t, r, B);

        ASSERT_ALWAYS((double)T/t < (double) R/n);
        ASSERT_ALWAYS(B*n <= R);

        glob.m = r * t;

        fprintf(stderr, "# [%2.2lf] number of pairs is %zu\n", seconds(),
                glob.ab->nab);
    }
    int odd_ab = glob.ab->nab & 1;
    size_t global_nab = glob.ab->nab;

    // TODO MPI --> agree on all of the above (perhaps do it only on the
    // master node at first).
    
    struct prime_data * primes = suitable_crt_primes();

    // TODO MPI --> collect complete list.
    
    mpz_init_set_ui(glob.P, 1);
    double log2_P = 0;
    for(int i = 0 ; i < glob.m ; i++) {
        log2_P += log(primes[i].p)/M_LN2;
        mpz_mul_ui(glob.P, glob.P, primes[i].p);
    }
    glob.prec = ceil((glob.nbits_sqrt + 128)/ log2_P);
    
    // TODO MPI --> agree on prec
    
    fprintf(stderr, "# [%2.2lf] Lifting to precision l=%d\n", seconds(), glob.prec);

    struct individual_contribution * contribs;

    contribs = malloc(glob.m * glob.n * glob.n * sizeof(struct individual_contribution));
    for(int i = 0 ; i < glob.m * glob.n * glob.n ; i++) {
        mpz_init(contribs[i].modN);
    }

    // XXX THE TWO LOOPS ON pgnum AND apnum CAN BE SPLIT WITH MPI in a grid-like
    // manner.
    for(int pgnum = 0 ; pgnum < t ; pgnum++) {
        fprintf(stderr, "# [%2.2lf] prime group %d\n", seconds(), pgnum);

        int i0 = pgnum*r;
        int i1 = i0 + r;

        for(int i = i0 ; i < i1 ; i++) {
            // note that the data computed locally is not read by other
            // nodes, so it needs not be duplicated everywhere if the
            // calculation is processed in parallel over several nodes.
            prime_initialization(&(primes[i]));
            prime_precomputations(&(primes[i]));
            // XXX This is a hack. we are evaluating the products
            // f_d^\epsilon*A*f_hat'(\hat\alpha)
            //
            // where A is the value denoted by the variable named A. It
            // is defined as:
            //
            // A(\hat\alpha) = (\prod_{(a,b)}(f_da-b\hat\alpha)
            //
            // Instead of fixing what's missing in A, we use shortcuts.
            // Extra f_d coefficients are added to the evaluation array,
            // and the derivative is eliminated later on by avoiding the
            // normalization in the Lagrange step.
            for(int j =  0 ; j < glob.n ; j++) {
                if (odd_ab) {
                    mpz_set(primes[i].evals->coeff[j], glob.pol->f[glob.n]);
                } else {
                    mpz_set_ui(primes[i].evals->coeff[j], 1);
                }
            }
        }

        // rational product tree: all prime powers.
        rat_ptree_t * rat_ptree = rat_ptree_build(primes, i0, i1);
        fprintf(stderr, "# [%2.2lf] product tree %d ready\n", seconds(), pgnum);

        for(int apnum = 0 ; apnum < s ; apnum++) {
            fprintf(stderr, "# [%2.2lf] A part %d\n", seconds(), apnum);

            size_t off0 = apnum * glob.ab->totalsize / s;
            size_t off1 = (apnum+1) * glob.ab->totalsize / s;

            poly_t P;
            poly_alloc(P, glob.n);
            poly_alloc(glob.t_abpoly, 1);
            ab_source_rewind(glob.ab);
            accumulate_ab_poly(P, glob.ab, off0, off1);
            fprintf(stderr, "# [%2.2lf] A%d product ready\n", seconds(), apnum);
            poly_free(glob.t_abpoly);

            // reduce modulo the top level (in place).
            for(int i = 0 ; i < glob.n ; i++) {
                mpz_mod(P->coeff[i], P->coeff[i], rat_ptree->zx);
                mpz_realloc(P->coeff[i], mpz_size(P->coeff[i]));
            }
            cleandeg(P, glob.n);
            reduce_poly_mod_rat_ptree(P, rat_ptree);
            poly_free(P);

            for(int i = i0 ; i < i1 ; i++) {
                reduce_poly_mod_alg_ptree(primes[i].T, &(primes[i]));
            }
        }

        rat_ptree_clear(rat_ptree);
        for(int i = i0 ; i < i1 ; i++) {
            alg_ptree_clear(primes[i].T);
            primes[i].T = NULL;
        }

        for(int i = i0 ; i < i1 ; i++) {
            prime_sqrt(&primes[i]);
            prime_postcomputations(contribs + i * glob.n * glob.n, &primes[i]);
        }

        for(int i = i0 ; i < i1 ; i++) {
            prime_cleanup(&(primes[i]));
        }
    }
    free(primes);

#if 0
    {
        for(int i = 0 ; i < glob.m ; i++) {
            for(int j =  0 ; j < glob.n ; j++) {
                for(int k = 0 ; k < glob.n ; k++) {
                    struct individual_contribution * c;
                    c = contribs + ((i * glob.n) + j) * glob.n + k;
                    printf(" %"PRIu64, c->ratio);
                }
                for(int s = 0 ; s < glob.m * glob.n ; s++) {
                    printf(" %d", s == (i * glob.n + j));
                }
                printf("\n");
            }
        }
        mpz_t z;
        mpz_init(z);
        mpz_ui_pow_ui(z, 2, 64);
        for(int k = 0 ; k < glob.n ; k++) {
            for(int s = 0 ; s < glob.n + glob.m * glob.n ; s++) {
                if (s == k) {
                    gmp_printf(" %Zd", z);
                } else {
                    printf(" 0");
                }
            }
            printf("\n");
        }
        mpz_clear(z);
    }
#endif

    {
        mpz_t sqrt_modN;
        mpz_init(sqrt_modN);
        solve_knapsack(sqrt_modN, contribs, (global_nab+odd_ab) / 2);
        mpz_clear(sqrt_modN);
    }

    /****************************************************************/

    for(int i = 0 ; i < glob.m * glob.n * glob.n ; i++) {
        mpz_clear(contribs[i].modN);
    }
    free(contribs);

    {
        for(int i = glob.n ; i >= 0 ; i--) {
            mpz_clear(glob.f_hat[i]);
        }
        for(int i = glob.n - 1 ; i >= 0 ; i--) {
            mpz_clear(glob.f_hat_diff[i]);
        }
        free(glob.f_hat);
        free(glob.f_hat_diff);
    }

    mpz_clear(glob.P);
    cado_poly_clear(glob.pol);
    poly_free(glob.F);

    ab_source_clear(glob.ab);
    param_list_clear(pl);
    return 0;
}

#if 0

    // Init F to be the algebraic polynomial
    poly_alloc(F, degree);
    for (i = pol->degree; i >= 0; --i)
        poly_setcoeff(F, i, pol->f[i]);

    // Init prd to 1.
    poly_alloc(prd->p, pol->degree);
    mpz_set_ui(prd->p->coeff[0], 1);
    prd->p->deg = 0;
    prd->v = 0;

    // Allocate tmp
    poly_alloc(tmp->p, 1);

    // Accumulate product
    int nab = 0, nfree = 0;
#if 0
    // Naive version, without subproduct tree
    while (fscanf(depfile, "%ld %lu", &a, &b) != EOF) {
        if (!(nab % 100000))
            fprintf(stderr, "# Reading ab pair #%d at %2.2lf\n", nab,
                    seconds());
        if ((a == 0) && (b == 0))
            break;
        polymodF_from_ab(tmp, a, b);
        polymodF_mul(prd, prd, tmp, F);
        nab++;
    }
#else
    // With a subproduct tree
    {
        polymodF_t *prd_tab;
        unsigned long lprd = 1;	/* number of elements in prd_tab[] */
        unsigned long nprd = 0;	/* number of accumulated products in prd_tab[] */
        prd_tab = (polymodF_t *) malloc(lprd * sizeof(polymodF_t));
        poly_alloc(prd_tab[0]->p, F->deg);
        mpz_set_ui(prd_tab[0]->p->coeff[0], 1);
        prd_tab[0]->p->deg = 0;
        prd_tab[0]->v = 0;
        while (fscanf(depfile, "%ld %lu", &a, &b) != EOF) {
            if (!(nab % 100000))
                fprintf(stderr, "# Reading ab pair #%d at %2.2lf\n", nab,
                        seconds());
            if ((a == 0) && (b == 0))
                break;
            polymodF_from_ab(tmp, a, b);
            prd_tab = accumulate_fast(prd_tab, tmp, F, &lprd, nprd++);
            nab++;
            if (b == 0)
                nfree++;
        }
        fprintf(stderr, "# Read %d including %d free relations\n", nab,
                nfree);
        ASSERT_ALWAYS((nab & 1) == 0);
        ASSERT_ALWAYS((nfree & 1) == 0);
        accumulate_fast_end(prd_tab, F, lprd);
        fclose(depfile);

        poly_copy(prd->p, prd_tab[0]->p);
        prd->v = prd_tab[0]->v;
        for (i = 0; i < (long) lprd; ++i)
            poly_free(prd_tab[i]->p);
        free(prd_tab);
    }
#endif

    fprintf(stderr, "Finished accumulating the product at %2.2lf\n",
            seconds());
    fprintf(stderr, "nab = %d, nfree = %d, v = %d\n", nab, nfree, prd->v);
    fprintf(stderr, "maximal polynomial bit-size = %lu\n",
            (unsigned long) poly_sizeinbase(prd->p, pol->degree - 1, 2));

    p = FindSuitableModP(F);
    fprintf(stderr, "Using p=%lu for lifting\n", p);

    double tm = seconds();
    polymodF_sqrt(prd, prd, F, p);


    fprintf(stderr, "Square root lifted in %2.2lf\n", seconds() - tm);
    mpz_t algsqrt, aux;
    mpz_init(algsqrt);
    mpz_init(aux);
    poly_eval_mod_mpz(algsqrt, prd->p, pol->m, pol->n);
    mpz_invert(aux, F->coeff[F->deg], pol->n);	// 1/fd mod n
    mpz_powm_ui(aux, aux, prd->v, pol->n);	// 1/fd^v mod n
    mpz_mul(algsqrt, algsqrt, aux);
    mpz_mod(algsqrt, algsqrt, pol->n);

    gmp_fprintf(stderr, "Algebraic square root is: %Zd\n", algsqrt);

    depfile = fopen(param_list_lookup_string(pl, "ratdepfile"), "r");
    ASSERT_ALWAYS(depfile != NULL);
    gmp_fscanf(depfile, "%Zd", aux);
    fclose(depfile);

    gmp_fprintf(stderr, "Rational square root is: %Zd\n", aux);
    fprintf(stderr, "Total square root time is %2.2lf\n", seconds());

    int found = 0;
    {
        mpz_t g1, g2;
        mpz_init(g1);
        mpz_init(g2);

        // First check that the squares agree
        mpz_mul(g1, aux, aux);
        mpz_mod(g1, g1, pol->n);
        if (mpz_cmp_ui(pol->g[1], 1) != 0) {
            // case g(X)=m1*X+m2 with m1 != 1
            // we should have prod (a+b*m2/m1) = A^2 = R^2/m1^(nab-nfree)
            // and therefore nab should be even
            if (nab & 1) {
                fprintf(stderr, "Sorry, but #(a, b) is odd\n");
                fprintf(stderr,
                        "Bug: this should be patched! Please report your buggy input\n");
                printf("Failed\n");
                return 0;
            }
            mpz_powm_ui(g2, pol->g[1], (nab - nfree) >> 1, pol->n);
            mpz_mul(algsqrt, algsqrt, g2);
            mpz_mod(algsqrt, algsqrt, pol->n);
        }
        mpz_mul(g2, algsqrt, algsqrt);
        mpz_mod(g2, g2, pol->n);
        if (mpz_cmp(g1, g2) != 0) {
            fprintf(stderr, "Bug: the squares do not agree modulo n!\n");
            //      gmp_printf("g1:=%Zd;\ng2:=%Zd;\n", g1, g2);
        }
        mpz_sub(g1, aux, algsqrt);
        mpz_gcd(g1, g1, pol->n);

        if (mpz_cmp(g1, pol->n)) {
            if (mpz_cmp_ui(g1, 1)) {
                found = 1;
                gmp_printf("%Zd\n", g1);
            }
        }
        mpz_add(g2, aux, algsqrt);
        mpz_gcd(g2, g2, pol->n);
        if (mpz_cmp(g2, pol->n)) {
            if (mpz_cmp_ui(g2, 1)) {
                found = 1;
                gmp_printf("%Zd\n", g2);
            }
        }
        mpz_clear(g1);
        mpz_clear(g2);
    }
    if (!found)
        printf("Failed\n");

    cado_poly_clear(pol);
    mpz_clear(aux);
    mpz_clear(algsqrt);
    poly_free(F);
    poly_free(prd->p);
    poly_free(tmp->p);

    return 0;
#endif
