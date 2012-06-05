#include "cado.h"

/* This program is the simplest interface to the bare matrix
 * multiplication routine. It's meant to provide an easy way of benching,
 * and comparing, different matrix product implementations.
 *
 * It must be used with the INNER matrix, not the .info one. So either use
 * the vanilla cado format (before balance), or one of the .h<i>.v<j>
 * matrices as output by balance.
 */

#include <stdio.h>
#include <time.h>
#include <errno.h>
#include <limits.h>
#include <inttypes.h>
#include <unistd.h>
#include <string.h>
#include "bwc_config.h"
#include "matmul.h"
#include "macros.h"
#include "params.h"
#include "worker-threads.h"
#include "utils.h"
#include "mpfq/abase_vbase.h"
#include "matmul-mf.h"
// #include "debug.h"

void usage()
{
    fprintf(stderr,
            "Usage: ./bench [--impl <implementation>] [--tmax <time>] [--nmax <n_iter>] [--nchecks <number> | --nocheck] [-r|--rebuild] [-t|--transpose] [--nthreads <number>] [--cycles <frequency>] -- <file0> [<file1> ... ]\n");
    exit(1);
}


struct private_args {
    matmul_t mm;
    void * src;
    void * dst;
};

struct bench_args {
    abase_vbase xx;
    param_list pl;
    const char * impl;
    int nthreads;
    int transpose;
    int nchecks;
    int rebuild;
    double freq;
    char ** mfiles;
    const char * source_vec;
    struct private_args * p;
    struct worker_threads_group * tg;
};

void init_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    int t0 = time(NULL);
    struct private_args * p = ba->p + tnum;

    p->mm = matmul_init(ba->xx, 0, 0, ba->mfiles[tnum], ba->impl, ba->pl, !ba->transpose);

    fprintf(stderr, "Expect to find or create cache file %s\n", p->mm->cachefile_name);

    if (!ba->rebuild && matmul_reload_cache(p->mm)) {
        pthread_mutex_lock(&tg->mu);
        fprintf(stderr, "T%d Reusing cache file for %s\n",
                tnum, ba->mfiles[tnum]);
        fprintf(stderr, "T%d Cache load time %ds wct\n",
                tnum, (int) time(NULL) - t0);
        pthread_mutex_unlock(&tg->mu);
    } else {
        clock_t ct0 = clock();
        pthread_mutex_lock(&tg->mu);
        fprintf(stderr, "T%d Building cache file for %s\n",
                tnum, ba->mfiles[tnum]);
        pthread_mutex_unlock(&tg->mu);
        ASSERT_ALWAYS(p->mm->store_transposed == ba->transpose);
        matrix_u32 m;
        mf_prepare_matrix_u32(p->mm, m, ba->mfiles[tnum]);
        matmul_build_cache(p->mm, m);
        pthread_mutex_lock(&tg->mu);
        fprintf(stderr, "T%d Cache build time %.2fs cpu\n",
                tnum, (double) (clock()-ct0) / CLOCKS_PER_SEC);
        fprintf(stderr, "T%d Saving cache file for %s\n",
                tnum, ba->mfiles[tnum]);
        pthread_mutex_unlock(&tg->mu);
        t0 = time(NULL);
        matmul_save_cache(p->mm);
        pthread_mutex_lock(&tg->mu);
        fprintf(stderr, "T%d Cache save time %ds wct\n",
                tnum, (int) time(NULL) - t0);
        pthread_mutex_unlock(&tg->mu);
    }

    pthread_mutex_lock(&tg->mu);
    fprintf(stderr, "T%u uses cache file %s\n",
            tnum,
            /* cache for mmt->locfile, */
            p->mm->cachefile_name);
    pthread_mutex_unlock(&tg->mu);

    unsigned int nr = p->mm->dim[0];
    unsigned int nc = p->mm->dim[1];
    matmul_aux(p->mm, MATMUL_AUX_GET_READAHEAD, &nr);
    matmul_aux(p->mm, MATMUL_AUX_GET_READAHEAD, &nc);

    ba->xx->vec_init(ba->xx, &p->src, nc);
    ba->xx->vec_init(ba->xx, &p->dst, nr);
    ba->xx->vec_set_zero(ba->xx, p->src, nc);
    ba->xx->vec_set_zero(ba->xx, p->dst, nr);
}/*}}}*/

void check_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    struct private_args * p = ba->p + tnum;
    abase_vbase_ptr A = ba->xx;

    unsigned int nr = p->mm->dim[0];
    unsigned int nc = p->mm->dim[1];
    unsigned int nr0 = nr;
    unsigned int nc0 = nc;
    matmul_aux(p->mm, MATMUL_AUX_GET_READAHEAD, &nr);
    matmul_aux(p->mm, MATMUL_AUX_GET_READAHEAD, &nc);

    void * dstT;
    void * srcT;
    A->vec_init(A, &dstT, nc);
    A->vec_init(A, &srcT, nr);
    A->vec_set_zero(A, dstT, nc);
    A->vec_set_zero(A, srcT, nr);

    void * check0;
    void * check1;
    A->vec_init(A, &check0, A->groupsize(A));
    A->vec_init(A, &check1, A->groupsize(A));

    A->vec_set(A, dstT, p->src, nc);
    A->vec_set(A, srcT, p->dst, nr);

    printf("T%d src(%u): %08" PRIx32 "\n", tnum,
            nc, crc32((unsigned long*) p->src, A->vec_elt_stride(A, nc0) / sizeof(unsigned long)));
    matmul_mul(p->mm, p->dst, p->src, 1);
    printf("T%d dst(%u): %08" PRIx32 "\n", tnum,
            nr, crc32((unsigned long*) p->dst, A->vec_elt_stride(A, nr0) / sizeof(unsigned long)));
    // debug_write(p->dst, abbytes(A, nr), "/tmp/Lmul");

    A->dotprod(A, check0, p->dst, srcT, nr0);

    printf("T%d srcT(%u): %08" PRIx32 "\n", tnum,
            nr, crc32((unsigned long*) srcT, A->vec_elt_stride(A, nr0) / sizeof(unsigned long)));
    matmul_mul(p->mm, dstT, srcT, 0);
    printf("T%d dstT(%u): %08" PRIx32 "\n", tnum,
            nc, crc32((unsigned long*) dstT, A->vec_elt_stride(A, nc0) / sizeof(unsigned long)));

    // debug_write(dstT, abbytes(A, nc), "/tmp/Rmul");

    A->dotprod(A, check1, p->src, dstT, nc0);

    if (A->vec_cmp(A, check0, check1, A->groupsize(A)) != 0) {
        pthread_mutex_lock(&tg->mu);
        fprintf(stderr, "T%d : Check failed\n", tnum);
        pthread_mutex_unlock(&tg->mu);
        abort();
    }

    A->vec_clear(A, &dstT, nc);
    A->vec_clear(A, &srcT, nr);
    A->vec_clear(A, &check0, A->groupsize(A));
    A->vec_clear(A, &check1, A->groupsize(A));


    matmul_aux(p->mm, MATMUL_AUX_ZERO_STATS);
}/*}}}*/

void mul_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    struct private_args * p = ba->p + tnum;
    matmul_mul(p->mm, ba->transpose ? p->src : p->dst,
            ba->transpose ? p->dst : p->src,
            !ba->transpose);
}/*}}}*/

void clear_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    struct private_args * p = ba->p + tnum;
    abase_vbase_ptr A = ba->xx;
    unsigned int nr = p->mm->dim[0];
    unsigned int nc = p->mm->dim[1];
    pthread_mutex_lock(&tg->mu);
    matmul_report(p->mm, ba->freq);
    printf("\n");
    pthread_mutex_unlock(&tg->mu);
    matmul_clear(p->mm);
    A->vec_clear(A, &p->src, nc);
    A->vec_clear(A, &p->dst, nr);
}/*}}}*/

void banner(int argc, char * argv[])
{
    /* print command line */
    fprintf (stderr, "# (%s) %s", CADO_REV, (argv)[0]);
    for (int i = 1; i < (argc); i++)
        fprintf (stderr, " %s", (argv)[i]);
    fprintf (stderr, "\n");

#ifdef  __GNUC__
    fprintf(stderr, "# Compiled with gcc " __VERSION__ "\n");
#endif
    fprintf(stderr, "# Compilation flags " CFLAGS "\n");
}

int main(int argc, char * argv[])
{
    struct bench_args ba[1];

    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    banner(argc, argv);

    ba->impl = "bucket";
    char * file = NULL;
    int nocheck = 0;
    ba->nchecks = 4;
    ba->nthreads = 1;
    ba->freq = 1;

    /* {{{ */
    param_list_init(ba->pl);
    argv++,argc--;
    param_list_configure_knob(ba->pl, "--transpose", &ba->transpose);
    param_list_configure_knob(ba->pl, "--rebuild", &ba->rebuild);
    param_list_configure_knob(ba->pl, "--nocheck", &nocheck);
    param_list_configure_alias(ba->pl, "--transpose", "-t");
    param_list_configure_alias(ba->pl, "--rebuild", "-r");

    int wild = 0;
    for( ; argc ; ) {
        if (param_list_update_cmdline(ba->pl, &argc, &argv)) { continue; }
        if (strcmp(argv[0],"--") == 0) {
            argv++, argc--;
            break;
        }
        if (argv[0][0] != '-' && wild == 0) {
            file = argv[0];
            argv++,argc--;
            wild++;
            continue;
        }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage();
    }

    unsigned int nmax = UINT_MAX;
    param_list_parse_uint(ba->pl, "nmax", &nmax);
    param_list_parse_int(ba->pl, "nthreads", &ba->nthreads);
    const char * unit;
    if (param_list_parse_double(ba->pl, "cycles", &ba->freq)) {
        unit = "cy/c";
    } else {
        unit = "ns/c";
    }

    if (ba->nthreads == 1) {
        if (file) {
            ba->mfiles = &file;
        } else {
            if (argc) {
                ba->mfiles = argv;
            } else {
                fprintf(stderr, "No file !\n");
                exit(1);
            }
        }
    } else {
        ba->mfiles = argv;
        if (file) {
            fprintf(stderr, "When using --nthreads, please specify files after a terminating --\n");
            exit(1);
        } else if (argc != ba->nthreads) {
            fprintf(stderr, "%u threads requested, but %u files given on the command line.\n", ba->nthreads, argc);
            exit(1);
        }
    }

    param_list_lookup_string(ba->pl, "matmul_bucket_methods");

    double tmax = 100.0;
    param_list_parse_double(ba->pl, "tmax", &tmax);
    param_list_parse_int(ba->pl, "nchecks", &ba->nchecks);
    if (nocheck) ba->nchecks = 0;

    const char * tmp;
    if ((tmp = param_list_lookup_string(ba->pl, "impl")) != NULL) {
        ba->impl = tmp;
    }

    unsigned int nbys = 64;

    abase_vbase_ptr A = ba->xx;

    /* may leave nbys unchanged ! */
    param_list_parse_uint(ba->pl, "nbys", &nbys);
    /* The api mandates that we set the desired value for nbys. Here,
     * that's not really our intent, since we really want to bench
     * the layer in its favorite working context. Most of the time,
     * setting nbys is pointless.
     */
    abase_vbase_oo_field_init_bygroupsize(A, nbys);
    A->set_groupsize(A, nbys);

    param_list_lookup_string(ba->pl, "srcvec");
    if (param_list_warn_unused(ba->pl)) {
        usage();
    }/*}}}*/

    ba->p = malloc(ba->nthreads * sizeof(struct private_args));
    ba->tg = worker_threads_init(ba->nthreads);
    fprintf(stderr, "Using implementation \"%s\"\n", ba->impl);
    worker_threads_do(ba->tg, (worker_func_t) &init_func, ba);

    /* {{{ display some info */
    uint64_t ncoeffs_total = 0;
    for(int tnum = 0 ; tnum < ba->nthreads ; tnum++) {
        struct private_args * p = ba->p + tnum;
        fprintf (stderr, "T%d %s: %u rows %u cols %" PRIu64 " coeffs\n",
                tnum, ba->mfiles[tnum], p->mm->dim[0], p->mm->dim[1], p->mm->ncoeffs);
        ncoeffs_total += p->mm->ncoeffs;
    }
    fprintf (stderr, "total %" PRIu64 " coeffs\n", ncoeffs_total);
    /* }}} */

    setup_seeding(1);
    /* {{{ Do checks if requested */
    for(int t = 0; t < ba->nchecks ; t++) {
        /* create deterministic test values */
        for(int tnum = 0 ; tnum < ba->nthreads ; tnum++) {
            struct private_args * p = ba->p + tnum;
            A->vec_random(A, p->src, p->mm->dim[1]);
            A->vec_random(A, p->dst, p->mm->dim[0]);
            /* If we want shared vectors, this is the way to go. */
            /* Note that for such a test, the clear_func must be skipped
             * or improved, since we don't really want to free() the same
             * pointer twice */
            /*
            if (ba->transpose)
                p->dst = ba->p[0].dst;
            else
                p->src = ba->p[0].src;
                */
        }
        worker_threads_do(ba->tg, (worker_func_t) &check_func, ba);
        fprintf(stderr, "Check %d ok\n", t);
    }
    if (ba->nchecks)
        printf("All %d checks passed\n", ba->nchecks);
    /* }}} */

    for(int tnum = 0 ; tnum < ba->nthreads ; tnum++) {
        struct private_args * p = ba->p + tnum;
        A->vec_random(A, p->src, p->mm->dim[1]);
        A->vec_random(A, p->dst, p->mm->dim[0]);
    }

    if ((tmp = param_list_lookup_string(ba->pl, "srcvec")) != NULL) {
        struct private_args * p = ba->p;
        if (ba->nthreads > 1) {
            fprintf(stderr, "srcvec incompatible with multithread\n");
            exit(1);
        }
        A->vec_set_zero(A, p->src, p->mm->dim[1]);
        A->vec_set_zero(A, p->dst, p->mm->dim[0]);
        FILE * f = fopen(tmp, "r");
        if (f == NULL) {
            fprintf(stderr, "fopen(%s): %s\n", tmp, strerror(errno));
            exit(1);
        }
        void * src = ba->transpose ? p->dst : p->src;
        void * dst = ba->transpose ? p->src : p->dst;
        int n = ba->transpose ? p->mm->dim[0] : p->mm->dim[1];
        fprintf(stderr, "reading %zu bytes from %s\n",
                n * sizeof(uint64_t), tmp);
        int nread = fread(src, sizeof(uint64_t), n, f);
        if (nread != n) {
            fprintf(stderr, "short read (%d < %d)\n", nread, n);
            exit(1);
        }
        fclose(f);
        worker_threads_do(ba->tg, (worker_func_t) &mul_func, ba);
        char * dstvec;
        int rc = asprintf(&dstvec, "%s.dst", tmp);
        ASSERT_ALWAYS(rc >= 0);
        f = fopen(dstvec, "w");
        int nw = ba->transpose ? p->mm->dim[1] : p->mm->dim[0];
        fprintf(stderr, "writing %zu bytes to %s\n",
                nw * sizeof(uint64_t), dstvec);
        int nwritten = fwrite(dst, sizeof(uint64_t), nw, f);
        if (nwritten != nw) {
            fprintf(stderr, "short write (%d < %d)\n", nwritten, nw);
            exit(1);
        }
        fclose(f);
        fprintf(stderr, "Saved [%s]%s * [%s] to [%s]\n",
                p->mm->cachefile_name,
                ba->transpose ? "^T" : "",
                tmp,
                dstvec);
        free(dstvec);
    }



#define NLAST   10
    double last[10]={0,};
    double sum_last = 0;

    clock_t t0 = clock();
    double next = 0.25 * CLOCKS_PER_SEC;
    double t1 = 0;
    double t, dt;
    if (next > tmax * CLOCKS_PER_SEC) { next = tmax * CLOCKS_PER_SEC; }
    unsigned int n;
    printf("Note: timings in seconds per iterations are given in cpu seconds ;\n");
    printf("Note: with %d threads, this means %.2f walltime seconds.\n", ba->nthreads, 1.0 / ba->nthreads);
    for(n = 0 ; n < nmax ; ) {
        worker_threads_do(ba->tg, (worker_func_t) &mul_func, ba);
        t = clock();
        dt = t - t0;
        sum_last -= last[n % NLAST];
        last[n % NLAST] = (t - t1) / CLOCKS_PER_SEC;
        sum_last += (t - t1) / CLOCKS_PER_SEC;
        t1 = t;
        int nlast = NLAST;
        n++;
        if (n < nlast) nlast = n;
        if (dt > next || n == nmax - 1) {
            do { next += 0.25 * CLOCKS_PER_SEC; } while (dt > next);
            if (next > tmax * CLOCKS_PER_SEC) { next = tmax * CLOCKS_PER_SEC; }
            dt /= CLOCKS_PER_SEC;
            printf("%d iters in %2.fs, %.3f/1, %.2f %s (last %u : %.3f/1, %.2f %s)        \r",
                    n, dt, dt/n, ba->freq * 1.0e9 * dt/n/ncoeffs_total, unit,
                    nlast, sum_last/nlast, ba->freq * 1.0e9 *sum_last/nlast/ncoeffs_total, unit);
            fflush(stdout);
            if (dt > tmax)
                break;
        }
    }
    printf("\n");
    // printf("Scanned %lu coeffs in total\n", n * ncoeffs_total);

    worker_threads_do(ba->tg, (worker_func_t) &clear_func, ba);
    worker_threads_clear(ba->tg);
    param_list_clear(ba->pl);
    free(ba->p);

    return 0;
}
