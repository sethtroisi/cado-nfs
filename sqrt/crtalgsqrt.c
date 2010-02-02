#define _BSD_SOURCE     /* M_LN2 */
#define _POSIX_C_SOURCE 200112L
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
#include "knapsack.h"

#include "select_mpi.h"
#include "gmp-hacks.h"          // TODO: REMOVE !

static int verbose = 0;
int ncores = 2;
double print_delay = 1;
double ram_gb = 3.0;    // Number of gigabytes. Note that this is the
// maximum size of product that are obtained.
// The actual memory footprint can be quite
// considerably larger due to FFT allocation (by
// a constant factor, though).

#define MPI_MY_SIZE_T   MPI_UNSIGNED_LONG
#define MPI_MY_UINT64_T   MPI_UNSIGNED_LONG
#define MPI_MY_INT64_T   MPI_LONG
#define MPI_MY_MP_SIZE_T        MPI_LONG
#define MPI_MY_MP_LIMB_T        MPI_UNSIGNED_LONG
#define MPI_MY_GMP_INTERNAL_SIZE_FIELD_T MPI_INT

static void usage()
{
    fprintf(stderr, "usage: crtalgsqrt algdepfile ratdepfile polyfile\n");
    exit(1);
}

/* {{{ mpi-gmp helpers */
static int all_agree(void *buffer, int count, MPI_Datatype datatype,/*{{{*/
        MPI_Comm comm)
{
    int s;
    MPI_Type_size(datatype, &s);
    void * b_and = malloc(count * s);
    void * b_or = malloc(count * s);
    MPI_Allreduce(buffer, b_and, count, datatype, MPI_BAND, comm);
    MPI_Allreduce(buffer, b_or, count, datatype, MPI_BOR, comm);
    int ok = memcmp(b_and, b_or, count * s) == 0;
    free(b_and);
    free(b_or);
    return ok;
}/*}}}*/
static void broadcast_mpz(mpz_ptr z, int root, MPI_Comm comm) /*{{{*/
{
    // that's much, much uglier than it should.
    int k;
    MPI_Comm_rank(comm, &k);
    mp_size_t nlimbs = mpz_size(z);
    MPI_Bcast(&nlimbs, 1, MPI_MY_MP_SIZE_T, root, comm);
    if (k != root) _mpz_realloc(z, nlimbs);
    MPI_Bcast(z->_mp_d, nlimbs, MPI_MY_MP_LIMB_T, root, comm);
    MPI_Bcast(&z->_mp_size, 1, MPI_MY_GMP_INTERNAL_SIZE_FIELD_T, root, comm);
}/*}}}*/
static void send_reduce_mul_mpz(mpz_ptr z, int r, int s, MPI_Comm comm)/*{{{*/
{
    // r <-- s
    int k;
    MPI_Comm_rank(comm, &k);
    mp_size_t nlimbs = mpz_size(z);
    if (k == r) {
        MPI_Recv(&nlimbs, 1, MPI_MY_MP_SIZE_T, s, (r<<4), comm, MPI_STATUS_IGNORE);
        if (nlimbs == 0) {
            /* It's not fundamentally wrong, but in our case, it
             * indicates a bug, for sure. */
            fprintf(stderr, "Warning: received 0 from job %d\n", s);
            mpz_set_ui(z,0);
            return;
        }

        mpz_t t;
        mpz_init(t);
        _mpz_realloc(t, nlimbs);
        MPI_Recv(t->_mp_d, nlimbs, MPI_MY_MP_LIMB_T, s, 1+(r<<4), comm, MPI_STATUS_IGNORE);
        MPI_Recv(&(t->_mp_size), 1, MPI_MY_GMP_INTERNAL_SIZE_FIELD_T, s, 2+(r<<4), comm, MPI_STATUS_IGNORE);
        mpz_mul(z,z,t);
        mpz_clear(t);
    } else {
        MPI_Send(&nlimbs, 1, MPI_MY_MP_SIZE_T, r, (r<<4), comm);
        if (nlimbs == 0)
            return;

        MPI_Send(z->_mp_d, nlimbs, MPI_MY_MP_LIMB_T, r, 1+(r<<4), comm);
        MPI_Send(&(z->_mp_size), 1, MPI_MY_GMP_INTERNAL_SIZE_FIELD_T, r, 2+(r<<4), comm);
    }
}/*}}}*/
/* }}} */

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
    // size_t digitbytes_estim;

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

void ab_source_init(ab_source_ptr ab, const char * fname, int rank, int root, MPI_Comm comm)
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
        ASSERT_ALWAYS(all_agree(&rc, 1, MPI_INT, comm));
        ASSERT_ALWAYS(rc == 0);
        tsize=sbuf->st_size;
        // we have 2.5 non-digit bytes per file line. However we
        // don't know the line count, so we can't subtract. As a
        // guess, we read the first 16kb, and count the number of
        // lines in there.
        if (rank == root) {
            char buf[16384];
            FILE * f = fopen(fname, "r");
            fread(buf, 1, sizeof(buf), f);
            fclose(f);
            int nrows_16k = 0;
            for(unsigned int i = 0 ; i < sizeof(buf) ; i++) {
                nrows_16k += buf[i] == '\n';
            }
            ab->nab_estim = (double) tsize * nrows_16k / sizeof(buf);
            // ab->digitbytes_estim = tsize - 2.5 * ab->nab_estim;
        }
        // XXX OK, this requires endianness consistency.
        MPI_Bcast(&ab->nab_estim, 1, MPI_MY_SIZE_T, root, comm);
        ab->file_bases = malloc(2 * sizeof(size_t));
        ab->file_bases[0] = 0;
        ab->file_bases[1] = tsize;
        ab->totalsize = tsize;
    } else {
        size_t dummy;
        size_t hdrbytes;
        char line[ABFILE_MAX_LINE_LENGTH];
        if (rank == root) {
            FILE * f = fopen(fname, "r");
            char * xx = fgets(line, sizeof(line), f);
            DIE_ERRNO_DIAG(xx == NULL, "fgets", fname);
            rc = sscanf(line, "AB %zu %zu", &ab->nab, &dummy);
            DIE_ERRNO_DIAG(rc != 2, "parse", fname);
            hdrbytes = ftell(f);
            fclose(f);
        }
        MPI_Bcast(&hdrbytes, 1, MPI_MY_SIZE_T, root, comm);

        ab->sname_len = strlen(fname) + 10;
        ab->sname = malloc(ab->sname_len);
        for(ab->nfiles = 0 ; ; ab->nfiles++) {
            snprintf(ab->sname, ab->sname_len, "%s.prep.%d.rel.%d",
                    ab->prefix, ab->depnum, ab->nfiles);
            rc = stat(ab->sname, sbuf);
            ASSERT_ALWAYS(rc == 0 || errno == ENOENT);
            ASSERT_ALWAYS(all_agree(&rc, 1, MPI_INT, comm));
            if (rc < 0) break;
            tsize += sbuf->st_size - hdrbytes;
        }
        ASSERT_ALWAYS(ab->nfiles > 0);
        ab->nab_estim = ab->nab;
        // ab->digitbytes_estim = tsize - 5 * ab->nab_estim;
        ab->file_bases = malloc((ab->nfiles+1) * sizeof(size_t));
        ab->file_bases[0] = 0;
        for(int i = 0 ; i < ab->nfiles ; i++) {
            snprintf(ab->sname, ab->sname_len, "%s.prep.%d.rel.%d",
                    ab->prefix, ab->depnum, i);
            rc = stat(ab->sname, sbuf);
            ASSERT_ALWAYS(rc == 0);
            ASSERT_ALWAYS(all_agree(&rc, 1, MPI_INT, comm));
            ab->file_bases[i+1]=ab->file_bases[i] + sbuf->st_size;
        }
        ab->totalsize = ab->file_bases[ab->nfiles];
    }
    /*
       fprintf(stderr, "# [%2.2lf] %s: roughly %zu rows,"
       " %zu digits (%zu bits/c)\n",
       seconds(),
       ab->nfiles ? "kleinjung" : "cado",
       ab->nab_estim, ab->digitbytes_estim,
       (size_t) (ab->digitbytes_estim * M_LN10 / M_LN2));
       */
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
    fprintf(stderr, "# [%2.2lf] [----] (a,b) rewind to %s, pos %zu\n",
            seconds(),
            ab->nfiles ? ab->sname : ab->fname0, ab->cpos);
}
// }}}

/* {{{ POSIX threads stuff : work queues */

enum wq_task_code { LIFT_ROOT, LIFT_SQRT, };
struct wq_task {
    // TODO: it's somewhat ugly to resort to this, but anyway.
    enum wq_task_code what;
    int i, j;
    void * ptr;

    // these three are reserved fields
    int done;
    void * res; /* not necessarily meaningful, depends on the task */
    pthread_mutex_t m_[1];
    pthread_cond_t c_[1];
};

struct wq_tasklist_item {
    struct wq_task * t;
    struct wq_tasklist_item * next;
    struct wq_tasklist_item * prev;
};

struct work_queue {
    pthread_mutex_t m[1];
    pthread_cond_t c[1];   // used by waiters.
    pthread_t * clients;
    // TODO: Use a single tasklist item for a doubly linked list.
    struct wq_tasklist_item * head;
};

void wq_push(struct work_queue * wq, struct wq_task * t)
{
    if (t) {
        /* initialize the reserved fields */
        pthread_mutex_init(t->m_, NULL);
        pthread_cond_init(t->c_, NULL);
        t->done = 0;
    }

    struct wq_tasklist_item * t1 = malloc(sizeof(struct wq_tasklist_item));
    t1->t = t;
    pthread_mutex_lock(wq->m);
    t1->next = wq->head;
    wq->head = t1;
    pthread_cond_signal(wq->c);
    pthread_mutex_unlock(wq->m);
}

struct wq_task * wq_pop_wait(struct work_queue * wq)
{
    struct wq_tasklist_item * t1;
    pthread_mutex_lock(wq->m);
    for( ; wq->head == NULL ; ) {
        pthread_cond_wait(wq->c, wq->m);
    }
    t1 = wq->head;
    wq->head = t1->next;
    pthread_mutex_unlock(wq->m);
    struct wq_task * t = t1->t;
    free(t1);
    return t;
}

void * wq_subtask(struct work_queue * wq, struct wq_task * t);

void * wq_waiter(void * wq)
{
    for( ; ; ) {
        struct wq_task * t = wq_pop_wait((struct work_queue *) wq);
        /* t == NULL is an indication that we're reaching the end-of-work
         * signal */
        if (t == NULL)
            break;
        void * res = wq_subtask((struct work_queue *) wq, t);
        pthread_mutex_lock(t->m_);
        t->done = 1;
        t->res = res;
        pthread_cond_signal(t->c_);
        pthread_mutex_unlock(t->m_);
    }
    return NULL;
}

void wq_join(struct work_queue * wq MAYBE_UNUSED, struct wq_task * t)
{
    pthread_mutex_lock(t->m_);
    for( ; t->done == 0 ; ) {
        pthread_cond_wait(t->c_, t->m_);
    }
    pthread_mutex_unlock(t->m_);
    pthread_cond_destroy(t->c_);
    pthread_mutex_destroy(t->m_);
}

void wq_init(struct work_queue * wq, unsigned int n)
{
    pthread_mutex_init(wq->m, NULL);
    pthread_cond_init(wq->c, NULL);
    wq->head = NULL;
    wq->clients = malloc(n * sizeof(pthread_t));
    for(unsigned int i = 0 ; i < n ; i++) {
        pthread_create(wq->clients + i, NULL, &wq_waiter, wq);
    }
}

void wq_clear(struct work_queue * wq, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        // schedule end-of-work for everybody.
        wq_push(wq, NULL);
    }
    for(unsigned int i = 0 ; i < n ; i++) {
        pthread_join(wq->clients[i], NULL);
    }

    ASSERT_ALWAYS(wq->head == NULL);
    pthread_mutex_destroy(wq->m);
    pthread_cond_destroy(wq->c);
}

/* }}} */

// some global variables (sheeh !)

struct sqrt_globals {
    int m;      // nprimes
    int n;      // degree
    int s;      // number of shares of A
    int t;      // number of prime groups
    int r;      // number of primes in each prime group
    int prec;
    mpz_t P;    // prime product (not to the power prec)
    size_t nbits_sqrt;
    // size_t nbits_a;
    mpz_t * f_hat;
    mpz_t * f_hat_diff;
    double f_hat_coeffs;
    cado_poly pol;
    poly_t t_abpoly;
    poly_t F;
    ab_source ab;
    int lll_maxdim;
    int rank;
    int nprocs;
    struct work_queue wq[1];
    MPI_Comm acomm;     // same share of A
    MPI_Comm pcomm;     // same sub-product tree
    int arank, asize;
    int prank, psize;
};

struct sqrt_globals glob = { .lll_maxdim=50 };

// {{{ trivial utility
static const char * size_disp(size_t s, char buf[16])
{
    char * prefixes = "bkMGT";
    double ds = s;
    const char * px = prefixes;
    for( ; px[1] && ds > 500.0 ; ) {
        ds /= 1024.0;
        px++;
    }
    snprintf(buf, 10, "%.1f%c", ds, *px);
    return buf;
}
// }}}

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

void estimate_nbits_sqrt(size_t * sbits, ab_source_ptr ab) // , int guess)
{
    size_t abits[1];
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
    int n = glob.pol->degree;

    long double complex * eval_points = malloc(n * sizeof(long double complex));
    double * double_coeffs = malloc((n+1) * sizeof(double));
    double * evaluations = malloc(n * sizeof(double));

    // take the roots of f, and multiply later on to obtain the roots of
    // f_hat. Otherwise we encounter precision issues.
    for(int i = 0 ; i <= n ; i++) {
        double_coeffs[i] = mpz_get_d(glob.pol->f[i]);
    }

    int rc = poly_roots_longdouble(double_coeffs, n, eval_points);
    if (rc) {
        fprintf(stderr, "# [%2.2lf] Warning: rootfinder had accuracy problem with %d roots\n", seconds(), rc);
    }


    // {{{ compress the list of roots.
    int nreal = 0, ncomplex = 0, rs = 0;
    for(int i = 0 ; i < n ; i++) {
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
        eval_points[i] *= mpz_get_d(glob.pol->f[n]);
    }
    for(int i = 0 ; i <= n ; i++) {
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
    for(int i = 0 ; i < n ; i++) {
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
            long double complex y = a * mpz_get_d(glob.pol->f[n]);   
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
            evaluations[i] += log(fabs(mpz_get_d(glob.pol->f[n])));
        }
    }

    // multiply by the square of f_hat'(f_d\alpha).
    for(int i = 0 ; i < rs ; i++) {
        long complex double s = n;
        for(int j = n - 1 ; j >= 0 ; j--) {
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
    double * a_bounds = malloc(n * sizeof(double));
    double * sqrt_bounds = malloc(n * sizeof(double));

    for(int j = 0 ; j < n ; j++) {
        a_bounds[j] = 0;
        sqrt_bounds[j] = 0;
    }
    for(int i = 0 ; i < rs ; i++) {
        double * lmat = malloc(n * sizeof(double ));
        lagrange_polynomial_abs(lmat, double_coeffs, n, eval_points[i]);
        for(int j = 0 ; j < n ; j++) {
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
    for(int j = 0 ; j < n ; j++) {
        // note that we might have added up to n times the same thing.
        // (inequality a+b < 2max(a,b) )
        sqrt_bounds[j] += log(n);
        a_bounds[j] += log(n);

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
    for(int j = 0 ; j < n ; j++) {
        fprintf(stderr, " %.4Lg", a_bounds[j] / M_LN2);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "# [%2.2lf] logbounds per coeff of sqrt(A)", seconds());
    for(int j = 0 ; j < n ; j++) {
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

/*{{{ have to pre-declare prime_data */
struct alg_ptree_s;

#if 0
struct individual_contribution {
    uint64_t ratio;
    mpz_t modN;
};
#endif

struct prime_data {
    unsigned long p;
    unsigned long * r;
    // unsigned long rj;
    void * powers;      // see .cpp file.

    // computed somewhat late.
    mpz_t iHx;

    struct alg_ptree_s * T;     // always NULL.

    // after the rational ptree reduction, this contains the share of A
    // with coefficients reduced. Short-lived.
    poly_t A;

    poly_t evals;       // set of evaluations.
    poly_t lroots;      // (lifted) roots.
    poly_t sqrts;       // square roots of A(x)     (only in the end !)

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

//  rational product tree: all prime powers.
// it's quite easy to set up.
//
// As such the code does not use the worker threads at all.
// It is possible to do it in a distributed manner, however it is not
// clear that a lot is to be gained.

rat_ptree_t * rat_ptree_build_inner(struct prime_data * p, int i0, int i1)
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
    res->t0 = rat_ptree_build_inner(p, i0, i0+d);
    res->t1 = rat_ptree_build_inner(p, i0+d, i1);
    mpz_init(res->z);
    res->zx = res->z;
    mpz_mul(res->z, res->t0->zx, res->t1->zx);
    return res;
}

rat_ptree_t * rat_ptree_build(struct prime_data * p, int i0, int i1)
{
    if (glob.prank == 0) {
        fprintf(stderr, "# [%2.2lf] [P%d  ] rational ptree starts\n",
                seconds(), glob.arank);
    }
    rat_ptree_t * res = rat_ptree_build_inner(p, i0, i1);
    if (glob.prank == 0) {
        fprintf(stderr, "# [%2.2lf] [P%d  ] rational ptree done\n",
                seconds(), glob.arank);
    }
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

/* }}} */

/* {{{ product tree (algebraic) */
struct alg_ptree_s {
    struct alg_ptree_s * t0;
    struct alg_ptree_s * t1;
    poly_t s;
};
typedef struct alg_ptree_s alg_ptree_t;

#if 0
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
#endif

void alg_ptree_clear(alg_ptree_t * T)
{
    if (T == NULL) return;
    alg_ptree_clear(T->t0);
    alg_ptree_clear(T->t1);
    free(T);
}

void reduce_poly_mod_alg_ptree(alg_ptree_t * T, struct prime_data * p)
{
    FATAL_ERROR_CHECK(T != NULL, "not implemented");
    mpz_t ta;
    mpz_init(ta);

    mpz_srcptr px = power_lookup_const(p->powers, glob.prec);
    for(int j = 0 ; j < glob.n ; j++) {
        mpz_set(ta, p->A->coeff[glob.n - 1]);
        for(int k = glob.n - 2 ; k >= 0 ; k--) {
            mpz_mul(ta, ta, p->lroots->coeff[j]);
            mpz_add(ta, ta, p->A->coeff[k]);
            mpz_mod(ta, ta, px);
        }
        mpz_mul(p->evals->coeff[j], p->evals->coeff[j], ta);
        if (mpz_size(p->evals->coeff[j]) >= mpz_size(px) * 3/2)
            mpz_mod(p->evals->coeff[j], p->evals->coeff[j], px);
    }
    for(int k = glob.n - 1 ; k >= 0 ; k--) {
        mpz_realloc(p->A->coeff[k], 0);
    }
    mpz_clear(ta);
}

/* }}} */

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

    if (glob.rank == 0)
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
        // if (glob.rank == 0) fprintf(stderr, "#[%2.2lf] %lu\n", seconds(), p);
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
    if (glob.rank == 0)
        fprintf(stderr, "# [%2.2lf] Found all CRT primes\n", seconds());
    free(roots);

    return res;
}
/* }}} */

/* {{{ everything that happens only modulo one prime */

void inversion_lift(struct prime_data * p, mpz_ptr iHx, mpz_srcptr Hx, int precision)/* {{{ */
{
    double t0 = seconds();

    mpz_srcptr pk = power_lookup_const(p->powers, precision);
    assert(precision > 0);

    if (precision == 1) {
        mpz_invert(iHx, Hx, pk);
        return;
    }
    int lower = precision - precision / 2;

    mpz_srcptr pl = power_lookup_const(p->powers, lower);

    // we're going to recurse. change the long Hx_mod value stored by a
    // temporary small one. The problem with this approach is that we
    // store many reductions in memory.
    mpz_t Hx_save;
    mpz_init(Hx_save);
    mpz_mod(Hx_save, Hx, pl);
    // recurse.
    inversion_lift(p, iHx, Hx_save, lower);
    mpz_clear(Hx_save);

    if (seconds() > t0 + print_delay)
        fprintf(stderr, "# [%2.2lf] [P%dA%d] precision %d\n",
                seconds(),
                glob.arank, glob.prank,
                precision);

    mpz_t ta;
    mpz_init(ta);

    mpz_mul(ta, iHx, Hx);
    mpz_mod(ta, ta, pk);

    mpz_sub_ui(ta, ta, 1);
    mpz_mul(ta, ta, iHx);
    mpz_sub(iHx, iHx, ta);
    mpz_mod(iHx, iHx, pk);

    mpz_clear(ta);
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
void invsqrt_lift(struct prime_data * p, mpz_ptr A, mpz_ptr sx, int precision)
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
        int issquare = modul_field_sqrt(a, a, z, q);
        ASSERT_ALWAYS(issquare);
        modul_inv(a, a, q);
        if (modul_get_ul(a, q) & 1) modul_neg(a, a, q);
        mpz_set_ui(sx, modul_get_ul(a, q));
        modul_clear(a, q);
        modul_clear(z, q);
        modul_clearmod(q);
        return;
    }
    int lower = precision - precision / 2;

    mpz_srcptr pl = power_lookup_const(p->powers, lower);

    // we're going to recurse ; arrange so that we recurse on something
    // acceptably small.  The problem with this approach is that we store
    // many reductions in memory.
    mpz_t A_save;
    mpz_init(A_save);
    mpz_mod(A_save, A, pl);
    // recurse.
    invsqrt_lift(p, A_save, sx, lower);
    mpz_clear(A_save);

    if (seconds() > t0 + print_delay)
        fprintf(stderr, "# [%2.2lf] [P%dA%d] precision %d\n",
                seconds(),
                glob.arank, glob.prank,
                precision);

    mpz_srcptr pk = power_lookup_const(p->powers, precision);

    mpz_t ta;
    mpz_init(ta);

    mpz_mul(ta, sx, sx);
    mpz_mod(ta, ta, pk);
    mpz_mul(ta, ta, A);
    mpz_mod(ta, ta, pk);
    mpz_sub_ui(ta, ta, 1);
    if (mpz_odd_p(ta)) mpz_add(ta, ta, pk);
    mpz_div_2exp(ta, ta, 1);
    mpz_submul(sx, ta, sx);
    mpz_mod(sx, sx, pk);

    mpz_clear(ta);
}

void sqrt_lift(struct prime_data * p, mpz_ptr A, mpz_ptr sx, int precision)
{
    double t0 = seconds();
    int lower = precision - precision / 2;
    mpz_srcptr pk = power_lookup_const(p->powers, precision);
    mpz_srcptr pl = power_lookup_const(p->powers, lower);
    mpz_t A_save;
    mpz_init(A_save);
    mpz_mod(A_save, A, pl);
    invsqrt_lift(p, A_save, sx, lower);
    mpz_clear(A_save);

    if (seconds() > t0 + print_delay)
        fprintf(stderr, "# [%2.2lf] [P%dA%d] precision %d\n",
                seconds(),
                glob.arank, glob.prank,
                precision);
    // inverse square root now in sx.

    mpz_t tmp;
    mpz_init(tmp);

    mpz_mul(tmp, A, sx);
    mpz_mod(tmp, tmp, pl);
    // XXX This destroys A !!!
    mpz_submul(A, tmp, tmp);
    if (mpz_odd_p(sx)) mpz_add(sx, sx, pk);
    mpz_div_2exp(sx, sx, 1);
    mpz_addmul(tmp, sx, A);
    mpz_mod(sx, tmp, pk);

    mpz_clear(tmp);
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
// improves the approximation of rx and 1/f'(rx)
void root_lift_innerlevels(struct prime_data * p, mpz_ptr rx, mpz_ptr irx, int precision)/* {{{ */
{
    double t0 = seconds();
    assert(precision > 0);

    if (precision == 1)
        return;
    int lower = precision - precision / 2;

    // recurse.
    root_lift_innerlevels(p, rx, irx, lower);

    if (seconds() > t0 + print_delay)
        fprintf(stderr, "# [%2.2lf] [P%dA%d] precision %d\n",
                seconds(),
                glob.arank, glob.prank,
                precision);

    mpz_srcptr pk = power_lookup_const(p->powers, precision);

    mpz_t ta, tb;
    mpz_init(ta);
    mpz_init(tb);

    // we know r to half-precision, 1/f'(r) to half-precision.
    mp_poly_eval_mod(tb, glob.f_hat, glob.n, rx, pk);
    mpz_mul(ta, tb, irx);
    // printf("%zu %zu %zu %zu\n", mpz_size(pk), mpz_size(ta), mpz_size(tb), mpz_size(irx));
    mpz_sub(rx, rx, ta);
    mpz_mod(rx, rx, pk);

    mp_poly_eval_mod(tb, glob.f_hat_diff, glob.n-1, rx, pk);
    mpz_mul(ta, irx, tb);
    mpz_sub_ui(ta, ta, 1);
    mpz_mod(ta, ta, pk);
    mpz_mul(tb, ta, irx);
    mpz_sub(ta, irx, tb);
    mpz_mod(irx, ta, pk);

    mpz_clear(ta);
    mpz_clear(tb);
    /*
       if (precision < 40) {
       gmp_fprintf(stderr, "# [%2.2lf] %Zd\n", seconds(), rx);
       gmp_fprintf(stderr, "# [%2.2lf] %Zd\n", seconds(), p->invdev_rx);
       }
       */
}
void root_lift(struct prime_data * p, mpz_ptr rx, mpz_ptr irx, int precision)
{
    double t0 = seconds();
    assert(precision > 0);

    if (precision == 1)
        return;
    int lower = precision - precision / 2;

    // recurse.
    root_lift_innerlevels(p, rx, irx, lower);

    if (seconds() > t0 + print_delay)
        fprintf(stderr, "# [%2.2lf] [P%dA%d] precision %d\n",
                seconds(),
                glob.arank, glob.prank,
                precision);

    mpz_srcptr pk = power_lookup_const(p->powers, precision);

    mpz_t ta;
    mpz_init(ta);

    // we know r to half-precision, 1/f'(r) to half-precision.
    mp_poly_eval_mod(ta, glob.f_hat, glob.n, rx, pk);
    mpz_mul(ta, ta, irx);
    // printf("%zu %zu %zu %zu\n", mpz_size(pk), mpz_size(ta), mpz_size(tb), mpz_size(p->invdev_rx));
    mpz_sub(rx, rx, ta);
    mpz_mod(rx, rx, pk);
    // don't do the lower part (thus irx is not relevant).

    mpz_clear(ta);
}/* }}} */
/* }}} */

void prime_initialization(struct prime_data * p)/* {{{ */
{
    // trigger computation of p^glob.prec
    p->powers = power_lookup_table_init(p->p);
    mpz_init(p->iHx);

    poly_alloc(p->A, glob.n-1);
    poly_alloc(p->evals, glob.n-1);
    poly_alloc(p->lroots, glob.n-1);
    poly_alloc(p->sqrts, glob.n-1);
}
/* }}} */

void * prime_precomputations_doit(struct prime_data * p, int j)/* {{{ */
{
    /*
       mpz_srcptr px = power_lookup(p->powers, glob.prec);

       fprintf(stderr, "# [%2.2lf] considering p=%lu\n",
       seconds(), p->p);
       fprintf(stderr, 
       "# [%2.2lf] size(p^l) %.1f MB\n",
       seconds(), mpz_sizeinbase(px, 256) * 1.0e-6);
       */

    mpz_srcptr p1 = power_lookup_const(p->powers, 1);
    // mpz_srcptr px = power_lookup(p->powers, glob.prec);
    // fprintf(stderr, "# [%2.2lf] considering r=%lu\n", seconds(), p->r[j]);

    // lift the root.
    // p->rj = p->r[j];
    mpz_ptr rx = p->lroots->coeff[j];
    mpz_set_ui(rx, p->r[j]);

    mpz_t irx;
    mpz_init(irx);

    mp_poly_eval_mod(irx, glob.f_hat_diff, glob.n-1, rx, p1);
    mpz_invert(irx, irx, p1);

    fprintf(stderr, "# [%2.2lf] [P%dA%d] lifting p=%lu, r=%lu\n",
            seconds(), glob.arank, glob.prank, p->p, p->r[j]);
    root_lift(p, rx, irx, glob.prec);
    // fprintf(stderr, "# [%2.2lf] done\n", seconds());

    mpz_clear(irx);

    return NULL;
}

/* }}} */

void prime_postcomputations(int64_t * c64, mp_limb_t * cN, struct prime_data * p)/* {{{ */
{
    mpz_srcptr px = power_lookup_const(p->powers, glob.prec);

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
    fprintf(stderr, "# [%2.2lf] [P%dA%d] lifting H^-l\n", seconds(), glob.arank, glob.prank);
    inversion_lift(p, p->iHx, Hx, glob.prec);
    mpz_clear(Hx);
    // fprintf(stderr, "# [%2.2lf] done\n", seconds());

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

    mpz_t ta, tb;
    mpz_init(ta);
    mpz_init(tb);

    // XXX Eh ! mpi-me !
    for(int j = 0 ; j < glob.n ; j++) {
        mpz_srcptr rx = p->lroots->coeff[j];
        mpz_ptr sx = p->sqrts->coeff[j];

        // so we have this nice square root. The first thing we do on our
        // list is to scramble it by multiplying it with the inverse of
        // H^x...
        mpz_mul(sx, sx, p->iHx);
        mpz_mod(sx, sx, px);

        // Now use the evaluation of f_hat mod rx to obtain the lagrange
        // coefficients.
        mpz_set_ui(ta, 1);
        for(int k = glob.n - 1 ; k >= 0 ; k--) {
            if (k < glob.n - 1) {
                mpz_mul(ta, ta, rx);
                mpz_add(ta, ta, glob.f_hat[k+1]);
                mpz_mod(ta, ta, px);
            }
            // multiply directly with H^-x * sqrt
            mpz_mul(tb, ta, sx);
            if (k < glob.n - 1) {
                mpz_mod(tb, tb, px);
            }
            ASSERT_ALWAYS(mpz_cmp_ui(tb, 0) >= 0);
            ASSERT_ALWAYS(mpz_cmp(tb, px) < 0);

            // now the shortcuts.
            mpf_set_z(ratio, tb);
            mpf_div(ratio, ratio, pxf);
            mpf_mul_2exp(ratio, ratio, 64);
            mpz_set_f(z, ratio);

            uint64_t u;
#if GMP_LIMB_BITS == 64
            u = mpz_get_ui(z);
#else
            u = (uint64_t) mpz_getlimb(z,1);
            u <<= 32;
            u |= (uint64_t) mpz_getlimb(z,0);
#endif
            c64[j*glob.n+k] = (int64_t) u;

            mpz_mul(tb, tb, Hxm);
            mpz_mod(tb, tb, glob.pol->n);
            mp_size_t sN = mpz_size(glob.pol->n);
            ASSERT_ALWAYS(SIZ(tb) > 0);
            MPN_SET_MPZ(cN + (j*glob.n+k) * sN, sN, tb);
        }
    }
    mpz_clear(ta);
    mpz_clear(tb);
    mpz_clear(Hxm);
    mpz_clear(z);
    mpf_clear(pxf);
    mpf_clear(ratio);
}
/* }}} */

void * prime_sqrt_doit(struct prime_data * p, int j)/* {{{ */
{
    fprintf(stderr, "# [%2.2lf] [P%dA%d] lifting sqrt mod (%lu, x-%lu)\n", seconds(), glob.arank, glob.prank, p->p, p->r[j]);
    sqrt_lift(p, p->evals->coeff[j], p->sqrts->coeff[j], glob.prec);
    // fprintf(stderr, "# [%2.2lf] done\n", seconds());
    mpz_realloc(p->evals->coeff[j], 0);
    return NULL;
}/* }}} */

void prime_cleanup(struct prime_data * p)/* {{{ */
{
    power_lookup_table_clear(p->powers);
    mpz_clear(p->iHx);
    free(p->r);
    alg_ptree_clear(p->T);

    poly_free(p->A);
    poly_free(p->evals);
    poly_free(p->lroots);
    poly_free(p->sqrts);

    memset(p, 0, sizeof(struct prime_data));
}
/* }}} */

/* }}} */

void * wq_subtask(struct work_queue * wq MAYBE_UNUSED, struct wq_task * t)
{
    switch(t->what) {
        case LIFT_ROOT:
            return prime_precomputations_doit(t->ptr, t->j);
            break;
        case LIFT_SQRT:
            return prime_sqrt_doit(t->ptr, t->j);
            break;
        default: ASSERT_ALWAYS(0); break;
    }
}

/* {{{ knapsack stuff */
struct crtalgsqrt_knapsack {
    knapsack_object ks;
    mpz_t Px;
    mpz_t lcx;
    mpz_t fhdiff_modN;
    mpz_t sqrt_modN;
    const mp_limb_t * tabN;
};

void crtalgsqrt_knapsack_init(struct crtalgsqrt_knapsack * cks)
{
    knapsack_object_init(cks->ks);
    mpz_init(cks->Px);
    mpz_init(cks->lcx);
    mpz_ptr z = cks->fhdiff_modN;
    mpz_init(z);
    mpz_init(cks->sqrt_modN);
}

int crtalgsqrt_knapsack_callback(struct crtalgsqrt_knapsack * cks,
        unsigned long v, int64_t x)
{
    knapsack_object_srcptr ks = cks->ks;
    const int64_t * c64 = cks->ks->tab;
    const mp_limb_t * cN = cks->tabN;
    unsigned int s;
    char * signs = malloc(ks->nelems+1);

    memset(signs, 0, ks->nelems+1);
    for(s = 0 ; s < ks->nelems ; s++) {
        signs[s]=(v & (1UL << s)) ? '+' : '-';
    }

    // now try to look at the solution more closely */
    mpz_t e;
    mpz_init_set_ui(e, 0);
    int spurious = 0;
    for(int k = glob.n - 1 ; k >= 0 ; k--) {
        mpz_t z, t, zN, qz, rz;
        mpz_init_set_ui(z, cks->ks->bound);
        mpz_init_set_ui(zN, 0);
        mpz_init_set_ui(t, 0);
        mpz_init(qz);
        mpz_init(rz);
        int64_t sk = ks->bound;
        for(int j = 0 ; j < glob.m * glob.n ; j++) {
            mpz_set_uint64(t, c64[j * glob.n + k]);
            mpz_t w;
            mp_size_t sN = mpz_size(glob.pol->n);
            MPZ_INIT_SET_MPN(w, cN + sN * (j * glob.n + k), sN);
            if (v & (1UL << j)) {
                sk+=c64[j * glob.n + k];
                mpz_add(z, z, t);
                mpz_add(zN, zN, w);
            } else {
                sk-=c64[j * glob.n + k];
                mpz_sub(z, z, t);
                mpz_sub(zN, zN, w);
            }
            mpz_clear(w);
        }
        if (sk >= 2 * ks->bound || sk < 0) {
            fprintf(stderr, "[%s] recombination of coeff in X^%d yields noise (%"PRIu64")\n",signs, k, sk);
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
        mpz_submul(zN, qz, cks->Px);
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

    mpz_mul(e,e,cks->lcx);
    mpz_mod(e,e,glob.pol->n);

    mpz_mul(e,e,cks->fhdiff_modN);
    mpz_mod(e,e,glob.pol->n);

    mpz_set(cks->sqrt_modN, e);

    mpz_clear(e);

    if (spurious) {
        fprintf(stderr, "# [%2.2lf] %"PRIx64" (%s) %"PRId64" [SPURIOUS]\n",
                seconds(), v, signs, (int64_t) x);
    } else {
        gmp_fprintf(stderr, "# [%2.2lf] %"PRIx64" (%s) %"PRId64" [%Zd]\n",
                seconds(), v, signs, (int64_t) x, cks->sqrt_modN);
    }
    free(signs);

    return spurious == 0;
}

/* Compute all the stuff which will be needed for checking the solutions
*/
void crtalgsqrt_knapsack_prepare(struct crtalgsqrt_knapsack * cks, size_t lc_exp)
{
#if 0 /* {{{ display the knapsack contents */
    {
        for(int i = 0 ; i < glob.m ; i++) {
            for(int j =  0 ; j < glob.n ; j++) {
                for(int k = 0 ; k < glob.n ; k++) {
                    int64_t c = contribs64[ ((i * glob.n) + j) * glob.n + k];
                    printf(" %"PRId64, c);
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
#endif/*}}}*/

    mpz_powm_ui(cks->Px, glob.P, glob.prec, glob.pol->n);

    mpz_set(cks->lcx, glob.pol->f[glob.n]);
    mpz_powm_ui(cks->lcx, glob.pol->f[glob.n], lc_exp, glob.pol->n);
    mpz_invert(cks->lcx, cks->lcx, glob.pol->n);

    mpz_ptr z = cks->fhdiff_modN;
    // evaluate the derivative of f_hat in alpha_hat mod N, that is lc*m.
    mpz_set_ui(z, 0);
    for(int k = glob.n ; k > 0 ; k--) {
        mpz_mul(z, z, glob.pol->m);
        mpz_mul(z, z, glob.pol->f[glob.n]);
        mpz_addmul_ui(z, glob.f_hat[k], k);
        mpz_mod(z, z, glob.pol->n);
    }
    mpz_invert(z, z, glob.pol->n);

    cks->ks->cb_arg = cks;
    cks->ks->cb = (knapsack_object_callback_t) crtalgsqrt_knapsack_callback;

    cks->ks->bound = 1+ceil(log(glob.m*glob.n)/log(2));

    cks->ks->stride = glob.n;

    unsigned int nelems = cks->ks->nelems = glob.m * glob.n;
    unsigned int k1 = nelems / 2;
    unsigned int k2 = nelems - k1;
    uint64_t n1 = (1UL << k1);
    uint64_t n2 = (1UL << k2);
    char buf[16];

    fprintf(stderr,
            "# [%2.2lf] Recombination: dimension %u = %u + %u, %s needed\n",
            seconds(),
            nelems, k1, k2,
            size_disp((n1+n2) * sizeof(uint64_t), buf)
           );

}

void crtalgsqrt_knapsack_clear(struct crtalgsqrt_knapsack * cks)
{
    mpz_clear(cks->sqrt_modN);
    mpz_clear(cks->Px);
    mpz_clear(cks->lcx);
    mpz_clear(cks->fhdiff_modN);
    knapsack_object_clear(cks->ks);
}
/* }}} */

void mpi_custom_types_check()/*{{{*/
{
    int ret;
    /* If this ever fails, we have to do something to the defines above
     * (and maybe touch the cmake files) */
    MPI_Type_size(MPI_MY_SIZE_T,&ret);   ASSERT_ALWAYS(ret==sizeof(size_t));
    MPI_Type_size(MPI_MY_MP_SIZE_T,&ret);ASSERT_ALWAYS(ret==sizeof(mp_size_t));
    MPI_Type_size(MPI_MY_MP_LIMB_T,&ret);ASSERT_ALWAYS(ret==sizeof(mp_limb_t));
    MPI_Type_size(MPI_MY_UINT64_T,&ret);ASSERT_ALWAYS(ret==sizeof(uint64_t));
    MPI_Type_size(MPI_MY_INT64_T,&ret);ASSERT_ALWAYS(ret==sizeof(int64_t));
    {
        mpz_t z;
        MPI_Type_size(MPI_MY_GMP_INTERNAL_SIZE_FIELD_T,&ret);
        ASSERT_ALWAYS(ret==sizeof(z->_mp_size));
    }
}/*}}}*/

// {{{ parallelism parameters selection
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
// Obeying [3], either by choosing the smallest/largest admissible
// power of two for r, or simply the smallest/largest admissible
// value, we deduce r and B. (large r will yield more
// parallelism-multithreading, but OTOH small r will yield less
// time).

void get_parameters(int * pr, int * ps, int * pt)
{
    int n=0,r=0,s=0,t=0;
    size_t R = ram_gb * (1UL << 30) * 8UL;
    n = glob.n;
    size_t T = glob.nbits_sqrt + 80;    // some margin for knapsack.
    double ramfraction = 2*T*n / (double) R;
    s = ceil(ramfraction);
    t = ceil((double)n*T / R);
    int dmax = glob.lll_maxdim - n;
    double rmax = (double) dmax / n / t;
    r = rmax;
    if (r == 0) {
        r = 1;
        fprintf(stderr, "Warning: resorting to choosing r=1, which will give a reconstruction effort of dimension %d\n",
                n*t);
    }
    // for(int mask = ~0 ; r & mask ; mask <<= 1) r &= mask;
    // try r=1 always.
    r = 1;
    size_t B = ceil((double)T/t/r);

    char buf[16];
    char buf2[16];
    fprintf(stderr, "# [%2.2lf] A is %s, %.1f%% of %s."
            " Splitting in %d parts\n",
            seconds(),
            size_disp(2*T*n/8.0, buf),
            100.0 * ramfraction, size_disp(ram_gb*1024.0*1048576.0, buf2), s);

    fprintf(stderr, "# [%2.2lf] %d groups of %d prime powers,"
            " each of %zu bits\n",
            seconds(), t, r, B);

    /* If these asserts fail, it means that the constraints are
     * impossible to satisfy */
    ASSERT_ALWAYS((double)T/t < (double) R/n);
    ASSERT_ALWAYS(B*n <= R);

    fprintf(stderr, "# [%2.2lf] number of pairs is %zu\n", seconds(),
            glob.ab->nab);
    *pr = r;
    *ps = s;
    *pt = t;
}
/* }}} */

void lifting_roots(struct prime_data * primes, int i0, int i1)/*{{{*/
{
    int pgnum = glob.arank;
    int n = glob.n;

    //  lifting all roots mod all primes (per pgroup)
    /* Again, this one can be distributed t-fold. In this case, it's
     * quite likely that there actually _is_ something to gain, since the
     * tasks are so clearly independent from one another. The code below
     * tries to give it a simple go.
     */
    if (glob.prank == 0) {
        fprintf(stderr, "# [%2.2lf] [P%d  ] lifting roots starts\n",
                seconds(), glob.arank);

    }
    {
        struct wq_task * lift_tasks = malloc((i1-i0)*n*sizeof(struct wq_task));
        int pushed = 0;
        for(int j = 0 ; j < n ; j++) {
            for(int i = i0 ; i < i1 ; i++) {
                int k = (i-i0) * n + j;
                if (k % glob.psize != glob.prank)
                    continue;
                struct wq_task * task = lift_tasks + k;
                task->what = LIFT_ROOT;
                task->ptr = primes + i;
                task->j = j;
                wq_push(glob.wq, task);
                pushed++;
            }
        }
        /* we're doing nothing */
        for(int j = 0 ; j < n ; j++) {
            for(int i = i0 ; i < i1 ; i++) {
                int k = (i-i0) * n + j;
                if (k % glob.psize != glob.prank)
                    continue;
                struct wq_task * task = lift_tasks + (i-i0) * n + j;
                wq_join(glob.wq, task);
            }
        }
        free(lift_tasks);
    }

    if (glob.prank == 0) {
        fprintf(stderr, "# [%2.2lf] [P%d  ] lifting roots partially done, now sharing\n",
                seconds(), pgnum);
    }

    {
        for(int j = 0 ; j < n ; j++) {
            for(int i = i0 ; i < i1 ; i++) {
                int k = (i-i0) * n + j;
                struct prime_data * p = primes + i;
                mpz_ptr z = p->lroots->coeff[j];
                broadcast_mpz(z, k % glob.psize, glob.pcomm);
            }
        }
    }
    if (glob.prank == 0) {
        fprintf(stderr, "# [%2.2lf] [P%d  ] lifting roots ends\n",
                seconds(), pgnum);
    }

    for(int i = i0 ; i < i1 ; i++) {
        for(int j = 0 ; j < n ; j++) {
            ASSERT_ALWAYS(mpz_size(primes[i].lroots->coeff[j]));
        }
    }
}/*}}}*/

/*{{{ a_poly_read_share and companion */ 
#define ABPOLY_OFFSET_THRESHOLD        16
// NOTE: This does not depend on p (nor r of course).
// the accumulation is done for all data between:
// the first data line starting at offset >= off0 (inclusive)
// the first data line starting at offset >= off1 (exclusive)
size_t accumulate_ab_poly(poly_t P, ab_source_ptr ab, size_t off0, size_t off1)
{
    size_t res = 0;
    mpz_set_ui(P->coeff[0], 1);
    P->deg = 0;
    if (off1 - off0 < ABPOLY_OFFSET_THRESHOLD) {
        ab_source_move_afterpos(ab, off0);
        for( ; ab->tpos < off1 ; res++) {
            int64_t a;
            uint64_t b;
            int r = ab_source_next(ab, &a, &b);
            FATAL_ERROR_CHECK(!r, "dep file ended prematurely\n");
            poly_from_ab_monic(glob.t_abpoly, a, b);
            polymodF_mul_monic(P, P, glob.t_abpoly, glob.F);
        }
        return res;
    }
    size_t d = (off1 - off0) / 2;
    poly_t Pl, Pr;
    poly_alloc(Pl, glob.n);
    poly_alloc(Pr, glob.n);
    res += accumulate_ab_poly(Pl, ab, off0, off0 + d);
    res += accumulate_ab_poly(Pr, ab, off0 + d, off1);
    polymodF_mul_monic(P, Pl, Pr, glob.F);
    poly_free(Pl);
    poly_free(Pr);
    return res;
}

size_t a_poly_read_share(poly_t P, size_t off0, size_t off1)
{
    size_t nab_loc = 0;
    if (glob.arank == 0) {
        poly_alloc(glob.t_abpoly, 1);
        ab_source_rewind(glob.ab);
        nab_loc = accumulate_ab_poly(P, glob.ab, off0, off1);
        fprintf(stderr, "# [%2.2lf] [P%dA%d] product ready (%zu pairs)\n", seconds(), glob.arank, glob.prank, nab_loc);
        poly_free(glob.t_abpoly);
    }
    MPI_Allreduce(MPI_IN_PLACE, &nab_loc, 1, MPI_MY_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    for(int i = 0 ; i < glob.n ; i++) { broadcast_mpz(P->coeff[i], 0, glob.acomm); }
    if (glob.arank == 0) {
        fprintf(stderr, "# [%2.2lf] [  A%d] product broadcasted\n", seconds(), glob.prank);
    }
    return nab_loc;
}/*}}}*/

void complete_reduction(poly_t P, struct prime_data * primes, int i0, int i1, rat_ptree_t * rat_ptree)/*{{{*/
{
    // reduce modulo the top level (in place).
    for(int i = 0 ; i < glob.n ; i++) {
        mpz_mod(P->coeff[i], P->coeff[i], rat_ptree->zx);
        mpz_realloc(P->coeff[i], mpz_size(P->coeff[i]));
    }
    cleandeg(P, glob.n);
    reduce_poly_mod_rat_ptree(P, rat_ptree);

    fprintf(stderr, "# [%2.2lf] [P%dA%d] rational reduction done\n", seconds(), glob.arank, glob.prank);

    for(int i = i0 ; i < i1 ; i++) {
        reduce_poly_mod_alg_ptree(primes[i].T, &(primes[i]));
        fprintf(stderr, "# [%2.2lf] [P%dA%d] algebraic reduction done (prime %d)\n", seconds(), glob.arank, glob.prank, i);
    }
}/*}}}*/

// {{{ starting values of the evaluation coefficients.
void set_starting_values(struct prime_data * primes, int i0, int i1, int odd_ab)
{
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
    for(int i = i0 ; i < i1 ; i++) {
        for(int j =  0 ; j < glob.n ; j++) {
            if (odd_ab) {
                mpz_set(primes[i].evals->coeff[j], glob.pol->f[glob.n]);
            } else {
                mpz_set_ui(primes[i].evals->coeff[j], 1);
            }
        }
    }
}
// }}}

void multiply_all_shares(struct prime_data * primes, int i0, int i1)
{
    // ok. for all eval points, collect the reduced stuff. we'll do the
    // collection in a tree-like manner.
    // this all happens within pcomm, which has size s.
    int me = glob.prank;
    int s = glob.s;

    for(int done_size = 1 ; done_size < s ; done_size <<= 1) {
        MPI_Barrier(glob.pcomm);
        if (me % done_size) continue;
        int receiver = me & ~done_size;
        int sender = me | done_size;
        if (sender >= s)
            continue;
        for(int i = i0 ; i < i1 ; i++) {
            mpz_srcptr px = power_lookup_const(primes[i].powers, glob.prec);
            for(int j = 0 ; j < glob.n ; j++) {
                mpz_ptr z = primes[i].evals->coeff[j];
                send_reduce_mul_mpz(z, receiver, sender, glob.pcomm);
                // could be off-loaded to a working queue.
                if (me == receiver) {
                    mpz_mod(z,z, px);
                }
            }
        }
    }
}

void local_square_roots(struct prime_data * primes, int i0, int i1)
{
    /* Only one job per pgroup has collected the complete data. Shall we
     * dispatch it again ? We might, if a need is felt, but that's
     * unclear */
    if (glob.prank != 0)
        return;

    struct wq_task * sqrt_tasks = malloc((i1-i0)*glob.n*sizeof(struct wq_task));
    int pushed = 0;
    for(int j = 0 ; j < glob.n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * glob.n + j;
            // if (k % glob.psize != glob.prank) continue;
            struct wq_task * task = sqrt_tasks + k;
            task->what = LIFT_SQRT;
            task->ptr = primes + i;
            task->j = j;
            wq_push(glob.wq, task);
            pushed++;
        }
    }
    /* we're doing nothing */
    for(int j = 0 ; j < glob.n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * glob.n + j;
            // if (k % glob.psize != glob.prank) continue;
            struct wq_task * task = sqrt_tasks + (i-i0) * glob.n + j;
            wq_join(glob.wq, task);
        }
    }
    free(sqrt_tasks);
}

void mpi_set_communicators()/*{{{*/
{
    int s = glob.s;
    int t = glob.t;
    if (s * t > glob.nprocs) {
        if (glob.rank == 0) {
            fprintf(stderr, "Error: need at least %d*%d=%d jobs (or more RAM)\n",
                    s,t,s*t);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    } else if (s * t < glob.nprocs) {
        if (glob.rank == 0) {
            fprintf(stderr, "Warning: only %d*%d=%d jobs needed, have %d\n",
                    s,t,s*t,glob.nprocs);
        }
    }

    MPI_Comm_split(MPI_COMM_WORLD, glob.rank / t, glob.rank, &glob.acomm);
    MPI_Comm_rank(glob.acomm, &glob.arank);
    MPI_Comm_size(glob.acomm, &glob.asize);
    ASSERT_ALWAYS(glob.asize == t);

    MPI_Comm_split(MPI_COMM_WORLD, glob.rank % t, glob.rank, &glob.pcomm);
    MPI_Comm_rank(glob.pcomm, &glob.prank);
    MPI_Comm_size(glob.pcomm, &glob.psize);
    ASSERT_ALWAYS(glob.psize == s);
}
/*}}}*/

int main(int argc, char **argv)
{
    int ret, i;
    // int size_guess = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &glob.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &glob.nprocs);

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

    param_list_parse_int(pl, "ncores", &ncores);
    param_list_parse_int(pl, "lll_maxdim", &glob.lll_maxdim);

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

    ab_source_init(glob.ab, param_list_lookup_string(pl, "depfile"),
            glob.rank, 0, MPI_COMM_WORLD);

    // note that for rsa768, this estimation takes only 10 minutes, so
    // it's not a big trouble.
    if (!param_list_parse_ulong(pl, "sqrt_coeffs_bits", &glob.nbits_sqrt)) {
        if (glob.rank == 0) {
            estimate_nbits_sqrt(&glob.nbits_sqrt, glob.ab); //, size_guess);
        }
        MPI_Bcast(&glob.nbits_sqrt, 1, MPI_MY_SIZE_T, 0, MPI_COMM_WORLD);
    }
    // MPI_Bcast(&glob.nbits_a, 1, MPI_MY_SIZE_T, 0, MPI_COMM_WORLD);
    // we no longer need to know nab, so let's drop it as a proof !
    glob.ab->nab = 0;

    glob.r=0;
    glob.s=0;
    glob.t=0;
    if (glob.rank == 0) {
        get_parameters(&glob.r,&glob.s,&glob.t);
    }
    MPI_Bcast(&glob.s, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&glob.t, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&glob.r, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int r = glob.r;
    int s = glob.s;
    int n = glob.n;
    int t = glob.t;
    glob.m = r * t;


    // this is sufficiently trivial to do.
    struct prime_data * primes = suitable_crt_primes();

    mpz_init_set_ui(glob.P, 1);
    double log2_P = 0;
    for(int i = 0 ; i < glob.m ; i++) {
        log2_P += log(primes[i].p)/M_LN2;
        mpz_mul_ui(glob.P, glob.P, primes[i].p);
    }
    glob.prec = ceil((glob.nbits_sqrt + 128)/ log2_P);

    ASSERT_ALWAYS(all_agree(&glob.prec, 1, MPI_INT, MPI_COMM_WORLD));

    if (glob.rank == 0) {
        char sbuf[32];
        fprintf(stderr, "# [%2.2lf] Lifting to precision l=%d (p^l is approx %s)\n", seconds(), glob.prec, size_disp(glob.prec * log(primes[0].p)/M_LN2 / 8, sbuf));
    }

    mpi_set_communicators();

    if (glob.rank == 0) {
        fprintf(stderr, "# [%2.2lf] starting %d worker threads on each node\n", seconds(), ncores);
    }
    wq_init(glob.wq, ncores);

    int pgnum = glob.rank % t;
    ASSERT_ALWAYS(glob.arank == pgnum);

    int i0 = pgnum*r;
    int i1 = i0 + r;

    for(int i = i0 ; i < i1 ; i++) prime_initialization(&(primes[i]));

    rat_ptree_t * rat_ptree = rat_ptree_build(primes, i0, i1);

    lifting_roots(primes, i0, i1);

    int apnum = glob.rank / t;
    ASSERT_ALWAYS(glob.prank == apnum);

    size_t off0 = apnum * glob.ab->totalsize / s;
    size_t off1 = (apnum+1) * glob.ab->totalsize / s;

    poly_t P;
    poly_alloc(P, glob.n);

    size_t nab = a_poly_read_share(P, off0, off1);
    set_starting_values(primes, i0, i1, nab & 1);

    complete_reduction(P, primes, i0, i1, rat_ptree);

    poly_free(P);

    MPI_Barrier(MPI_COMM_WORLD);

    multiply_all_shares(primes, i0, i1);

    /* some housekeeping */
    rat_ptree_clear(rat_ptree);
    for(int i = i0 ; i < i1 ; i++) {
        alg_ptree_clear(primes[i].T);
        primes[i].T = NULL;
    }

    local_square_roots(primes, i0, i1);

    size_t nc = glob.m * glob.n * glob.n;
    int64_t * contribs64 = malloc(nc * sizeof(int64_t));
    mp_size_t sN = mpz_size(glob.pol->n);
    memset(contribs64,0,nc * sizeof(int64_t));
    mp_limb_t * contribsN = malloc(nc * sN * sizeof(mp_limb_t));
    memset(contribsN, 0,  nc * sN * sizeof(mp_limb_t)); 

    if (glob.prank == 0) {

        for(int i = i0 ; i < i1 ; i++) {
            int disp = i * glob.n * glob.n;
            prime_postcomputations(contribs64 + disp, contribsN + disp * sN, &primes[i]);
        }
    }

    if (glob.rank == 0) {
        fprintf(stderr, "# [%2.2lf] clearing work queues\n", seconds());
    }
    wq_clear(glob.wq, ncores);

    for(int i = i0 ; i < i1 ; i++) prime_cleanup(&(primes[i]));
    free(primes);

    // now share the contribs !
    if (glob.prank == 0) {
        for(int i = 0 ; i < glob.m ; i++) {
            int root = i / r;
            int d64 = glob.n * glob.n;
            int dN = d64 * mpz_size(glob.pol->n);
            MPI_Bcast(contribs64 + i * d64, d64, MPI_MY_INT64_T, root, glob.acomm);
            MPI_Bcast(contribsN + i * dN, dN, MPI_MY_MP_LIMB_T, root, glob.acomm);
        }
    }

    // and only the leader does the knapsack.
    if (glob.rank == 0) {
        struct crtalgsqrt_knapsack cks[1];
        crtalgsqrt_knapsack_init(cks);
        crtalgsqrt_knapsack_prepare(cks, (nab+(nab&1)) / 2);
        cks->ks->tab = contribs64;
        cks->tabN = contribsN;
        knapsack_solve(cks->ks);
        crtalgsqrt_knapsack_clear(cks);
    }

    /****************************************************************/

    free(contribs64);
    free(contribsN);

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

    MPI_Finalize();

    return 0;
}

#if 0/*{{{*/

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
#endif/*}}}*/
