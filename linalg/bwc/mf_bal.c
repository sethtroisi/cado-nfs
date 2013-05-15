#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include "portability.h"
#include "utils.h"
#include "mf.h"
#include "balancing.h"
#include "rowset_heap.h"

/* This program computes how a matrix would have to be balanced for
 * fitting a job grid of given size. This does _not_ read the matrix,
 * only the row-weight and col-weight files are read.
 */

/* TODO: Now that all files input and output by bwc are permutation
 * independent, there is room for computing this permutation on the fly
 * if the balancing argument is not provided (and only in this case).
 * Essentially we would hook onto the mmt_fill_fields_from_balancing()
 * function called from mmt_finish_init().
 */
void usage(int rc) {
    fprintf(stderr,
            "Usage: ./mf_bal [options,flags] <nh> <nv> <matrix file>\n");
    fprintf(stderr,
            "Recognized options"
                " (<option_name>=<value>, or --<option_name> <value>:\n"
            " mfile       matrix file (can also be given freeform)\n"
            " rwfile      row weight file (defaults to <mfile>.rw)\n"
            " cwfile      col weight file (defaults to <mfile>.cw)\n"
            " like        balance in a way compatible with this file\n"
            " out         output file name (defaults to stdout)\n"
            "Recognized flags:"
            " --balance2d balance both dimensions\n"
            " --ascii     output in ascii\n"
            " --quiet     be quiet\n"
            " --rowperm   permute rows in priority (defaults to auto)\n"
            " --colperm   permute rows in priority (defaults to auto)\n"
            " --shuffled-product   prepare permutation for shuffled product\n"
           );
    exit(rc);
}

typedef int (*sortfunc_t) (const void *, const void *);

int revcmp_2u32(const uint32_t * a, const uint32_t * b)
{
    if (a[0] > b[0]) return -1;
    if (b[0] > a[0]) return 1;
    if (a[1] > b[1]) return -1;
    if (b[1] > a[1]) return 1;
    return 0;
}

/* {{{ Datatype containing description for slices of the matrix */
struct slice {
    uint32_t * r;
    uint32_t nrows;
    uint32_t coeffs;
    uint32_t i0;
};

int slice_finding_helper(const uint32_t * key, const struct slice * elem)/*{{{*/
{
    if (*key < elem->i0) {
        return -1;
    } else if (*key >= elem->i0 + elem->nrows) {
        return 1;
    }
    return 0;
}

unsigned int which_slice(const struct slice * slices, unsigned int ns, uint32_t k)
{
    const struct slice * found = (const struct slice *) bsearch(
            (const void *) &k, (const void *) slices,
            ns, sizeof(struct slice),
            (sortfunc_t) &slice_finding_helper);
    if (found == NULL) {
        return UINT_MAX;
    } else {
        return found - slices;
    }
}/*}}}*/

void free_slices(struct slice * slices, unsigned int n)
{
    unsigned int i;
    for(i = 0 ; i < n ; i++) {
        free(slices[i].r);
    }
    free(slices);
}

struct slice * alloc_slices(unsigned int water, unsigned int n)
{
    struct slice * res;
    unsigned int i;

    res = (struct slice*) malloc(n * sizeof(struct slice));

    uint32_t common_size = water / n;

    // we now require equal-sized blocks.
    ASSERT_ALWAYS(water % n == 0);

    for(i = 0 ; i < n ; i++) {
        res[i].nrows = common_size; // min_size + (i < (water % n));
        res[i].r = (uint32_t *) malloc(res[i].nrows * sizeof(uint32_t));
        res[i].coeffs = 0;
    }

    return res;
}
/* }}} */
/* {{{ shuffle_rtable: This is the basic procedure which dispatches rows
 * (or columns) in buckets according to our preferred strategy
 */
struct slice * shuffle_rtable(
        const char * text,
        uint32_t * rt,
        uint32_t n,
        unsigned int ns)
{
    uint32_t i;
    struct bucket * heap;
    struct slice * slices;
    clock_t t;

    t = -clock();
    qsort(rt, n, 2 * sizeof(uint32_t), (sortfunc_t) &revcmp_2u32);
    t += clock();
    fprintf(stderr, "sort time %.1f s\n", (double) t / CLOCKS_PER_SEC);

    
    t = -clock();

    slices = alloc_slices(n, ns);

    heap = (struct bucket *) malloc(ns * sizeof(struct bucket));
    for(i = 0 ; i < ns ; i++) {
        heap[i].s = 0;
        heap[i].i = i;
        heap[i].room = slices[i].nrows;
    }
    make_heap(heap, heap + ns);

    /* Then we're putting each row through the appropriate bucket, using
     * the heap to constantly update.
     * Eventually, we have constructed the set of rows going in the
     * bucket into slices[0]...slices[nslices-1]
     */
    for(i = 0 ; i < n ; i++) {
        int j = heap[0].i;
        int pos = slices[j].nrows-heap[0].room;
        assert(heap[0].room);
        slices[j].r[pos] = rt[2*i+1];
        heap[0].s += rt[2*i];
        heap[0].room--;
        pop_heap(heap, heap + ns);
        push_heap(heap, heap + ns);
    }

    t += clock();

    fprintf(stderr, "heap fill time %.1f\n", (double) t / CLOCKS_PER_SEC);

    qsort(heap, ns, sizeof(struct bucket), (sortfunc_t) &heap_index_compare);

    for(i = 0 ; i < ns ; i++) {
        int j = heap[i].i;
        assert(heap[i].i == (int) i);
        printf("%s slice %d, span=%ld, weight=%ld\n",
                text,
                i, slices[j].nrows - heap[i].room,
                heap[i].s);
        slices[j].coeffs = heap[i].s;
        if (heap[i].room != 0) {
            abort();
        }

#if 0
        // This is for giving the possibility of reading the matrix
        // linearly.
        if (sort_positionally) {
            qsort(slices[j].data, slices[j].nrows, sizeof(struct row),
                    (sortfunc_t) row_compare_index);
        }
#endif
    }
    free(heap);

    uint32_t i0 = 0;
    for(i = 0 ; i < ns ; i++) {
        slices[i].i0=i0;
        i0 += slices[i].nrows;
    }

    return slices;
}

/* }}} */

int main(int argc, char * argv[])
{
    param_list pl;
    const char * rwfile = NULL;
    const char * cwfile = NULL;
    const char * mfile = NULL;
    unsigned int wild =  0;
    int rc;
    int quiet =  0;
    int nh=0;
    int nv=0;
    int twodim=0;
    int withcoeffs=0;
    // int ascii_in = 0;
    // int ascii_out = 0;
    // int binary_in = 0;
    // int binary_out = 0;
    // int ascii_freq = 0;
    // int binary_freq = 0;
    /* Whether our priority is the balancing of rows or or columns. It's
     * a misleading name that should be changed. Anyway the default
     * should be right (max standard deviation to be balanced) */
    /* Note that we don't even have to read the .cw file if we're being
     * told that only rows matter (and conversely, of course) */
    // to be supported one day maybe. for now, it's not useful.
    int rowperm = 0;
    int colperm = 0;
    int shuffled_product = 0;
    // const char * ref_balance = NULL;
    int display_correlation = 0;

    param_list_init(pl);
    argv++,argc--;
    param_list_configure_switch(pl, "--quiet", &quiet);
    param_list_configure_switch(pl, "--display-correlation", &display_correlation);
    // param_list_configure_switch(pl, "--ascii-in", &ascii_in);
    // param_list_configure_switch(pl, "--binary-in", &binary_in);
    // param_list_configure_switch(pl, "--ascii-out", &ascii_out);
    // param_list_configure_switch(pl, "--binary-out", &binary_out);
    // param_list_configure_switch(pl, "--ascii-freq", &ascii_freq);
    // param_list_configure_switch(pl, "--binary-freq", &binary_freq);
    param_list_configure_switch(pl, "--balance2d", &twodim);
    param_list_configure_switch(pl, "--shuffled-product", &shuffled_product);
    param_list_configure_switch(pl, "--withcoeffs", &withcoeffs);

    for(;argc;) {
        char * q;
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;

        if (argv[0][0] != '-' && wild == 0 && (q = strchr(argv[0],'x')) != NULL) {
            nh = atoi(argv[0]);
            nv = atoi(q+1);
            wild+=2;
            argv++,argc--;
            continue;
        }

        if (argv[0][0] != '-' && wild == 0) { nh = atoi(argv[0]); wild++,argv++,argc--; continue; }
        if (argv[0][0] != '-' && wild == 1) { nv = atoi(argv[0]); wild++,argv++,argc--; continue; }
        if (argv[0][0] != '-' && wild == 2) {
            mfile = argv[0];
            wild++;
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "unknown option %s\n", argv[0]);
        exit(1);
    }

    if (nh == 0 || nv == 0)
        usage(1);

    if (rowperm && colperm) {
        fprintf(stderr, "--rowperm and --colperm are incompatible\n");
        usage(1);
    }

    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "mfile")) != NULL) {
        mfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "rwfile")) != NULL) {
        rwfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "cwfile")) != NULL) {
        cwfile = tmp;
    }

    if (mfile && !rwfile) {
        char * leakme;
        rwfile = leakme = build_mat_auxfile(mfile, "rw", ".bin");
    }

    if (mfile && !cwfile) {
        char * leakme;
        cwfile = leakme = build_mat_auxfile(mfile, "cw", ".bin");
    }

    if (rowperm) {
        if (!rwfile) { fprintf(stderr, "No rwfile given\n"); exit(1); }
    } else if (colperm) {
        if (!cwfile) { fprintf(stderr, "No cwfile given\n"); exit(1); }
    } else {
        /* If we have to decide automatically, we need everything ! */
        if (!rwfile) { fprintf(stderr, "No rwfile given\n"); exit(1); }
        if (!cwfile) { fprintf(stderr, "No cwfile given\n"); exit(1); }
    }

    if (!mfile) {
        fprintf(stderr, "Matrix file name (mfile) must be given, even though the file itself does not have to be present\n");
        exit(1);
    }

    balancing bal;
    balancing_init(bal);

    bal->h->nh = nh;
    bal->h->nv = nv;

    struct stat sbuf_rw[1];
    rc = stat(rwfile, sbuf_rw);
    if (rc < 0) { perror(rwfile); exit(1); }
    bal->h->nrows = sbuf_rw->st_size / sizeof(uint32_t);

    struct stat sbuf_cw[1];
    rc = stat(cwfile, sbuf_cw);
    if (rc < 0) { perror(cwfile); exit(1); }
    bal->h->ncols = sbuf_cw->st_size / sizeof(uint32_t);

    size_t maxdim = MAX(bal->h->nrows, bal->h->ncols);
    if (maxdim != bal->h->nrows) {
        fprintf(stderr, "Warning. More columns than rows. There could be bugs.\n");
    }

    // we are forcibly computing a symmetric permutation, according to a
    // situation where the _columns_ of the input matrix are the most
    // unbalanced.
    
    uint32_t elem_block = iceildiv(maxdim, nh * nv);
    // uint32_t rslice_size = nv * elem_block;
    uint32_t cslice_size = nh * elem_block;
    uint32_t padding = nv * nh * elem_block - maxdim;

    fprintf(stderr,
            "Using %" PRIu32" padding %ss to obtain %u blocks of %u*%"PRIu32"=%"PRIu32" %ss\n",
                padding, "col", nv, nh, elem_block, cslice_size, "col");


    struct stat sbuf_mat[1];

    rc = stat(mfile, sbuf_mat);
    if (rc < 0) {
        fprintf(stderr, "Reading %s: %s (not fatal)\n", mfile, strerror(errno));
        fprintf(stderr, "%s: %" PRIu32 " rows %" PRIu32 " cols\n",
                mfile, bal->h->nrows, bal->h->ncols);
        fprintf(stderr,
                "%s: main input file not present locally, total weight unknown\n", mfile);
        bal->h->ncoeffs = 0;
    } else {
        bal->h->ncoeffs = sbuf_mat->st_size / sizeof(uint32_t) - bal->h->nrows;
        if (withcoeffs) {
            if (bal->h->ncoeffs & 1) {
                fprintf(stderr, "Matrix with coefficient must have an even number of 32-bit entries for all (col index, coeff). Here, %"PRIu64" is odd.\n", bal->h->ncoeffs);
                abort();
            }
            bal->h->ncoeffs /= 2;
        }

        int extra = bal->h->ncols - bal->h->nrows;
        if (extra > 0) {
            fprintf(stderr,
                    "%s: %" PRIu32 " rows %" PRIu32 " cols"
                    " (%d extra cols)"
                    " weight %" PRIu64 "\n",
                    mfile, bal->h->nrows, bal->h->ncols,
                    extra,
                    bal->h->ncoeffs);
        } else if (extra < 0) {
            fprintf(stderr,
                    "%s: %" PRIu32 " rows %" PRIu32 " cols"
                    " (%d extra rows)"
                    " weight %" PRIu64 "\n",
                    mfile, bal->h->nrows, bal->h->ncols,
                    -extra,
                    bal->h->ncoeffs);
        } else {
            fprintf(stderr,
                    "%s: %" PRIu32 " rows %" PRIu32 " cols"
                    " weight %" PRIu64 "\n",
                    mfile, bal->h->nrows, bal->h->ncols, bal->h->ncoeffs);
        }
    }
    /* TODO: In case we rely on the rows, not columns for doing the
     * balancing, there is of course some code which must be modified
     * along the lines of the following.
     */

    /* read column data */
    uint32_t * colweights;
    {
        size_t n = maxdim + padding;
        bal->colperm = malloc(2 * n * sizeof(uint32_t));
        FILE * fcw = fopen(cwfile, "rb");
        if (fcw == NULL) { perror(cwfile); exit(1); }
        colweights = malloc(n * sizeof(uint32_t));
        memset(colweights, 0, n * sizeof(uint32_t));
        /* Padding cols have zero weight of course */
        double t_cw;
        t_cw = -wct_seconds();
        size_t nr = fread32_little(colweights, bal->h->ncols, fcw);
        t_cw += wct_seconds();
        fclose(fcw);
        if (nr < bal->h->ncols) {
            fprintf(stderr, "%s: short col count\n", cwfile);
            exit(1);
        }
        fprintf(stderr, "read %s in %.1f s (%.1f MB / s)\n", cwfile, t_cw,
                1.0e-6 * sbuf_cw->st_size / t_cw);
    }

    /* Compute the de-correlating permutation */
    {
        modulusul_t M;
        modul_initmod_ul(M, bal->h->ncols);
        residueul_t a,b;
        residueul_t ai,bi;
        modul_init(a, M);
        modul_init(b, M);
        modul_init(ai, M);
        modul_init(bi, M);
        modul_set_ul(a, (unsigned long) sqrt(bal->h->ncols), M);
        modul_set_ul(b, 42, M);

        for( ; modul_inv(ai, a, M) == 0 ; modul_add_ul(a,a,1,M)) ;
        modul_mul(bi, ai, b, M);
        modul_neg(bi, bi, M);

        bal->h->pshuf[0] = modul_get_ul(a, M);
        bal->h->pshuf[1] = modul_get_ul(b, M);
        bal->h->pshuf_inv[0] = modul_get_ul(ai, M);
        bal->h->pshuf_inv[1] = modul_get_ul(bi, M);

        modul_clear(a, M);
        modul_clear(b, M);
        modul_clear(ai, M);
        modul_clear(bi, M);
        modul_clearmod(M);
    }
    if (param_list_lookup_string(pl, "noshuffle")) { /* internal, for debugging. don't use */
        bal->h->pshuf[0] = 1;
        bal->h->pshuf[1] = 0;
        bal->h->pshuf_inv[0] = 1;
        bal->h->pshuf_inv[1] = 0;
    }

    /* prepare for qsort */
    double t_cw = -wct_seconds();
    double s1 = 0, s2 = 0;
    uint64_t tw = 0;
    for(size_t r = 0 ; r < maxdim + padding ; r++) {
        /* Column r in the matrix we work with is actually column rx in
         * the original matrix ! */
        size_t rx = r;
        if (r < bal->h->ncols) {
            rx = balancing_pre_unshuffle(bal, r);
            ASSERT(balancing_pre_shuffle(bal, rx) == r);
            ASSERT(rx < bal->h->ncols);
        }
        uint32_t w = colweights[rx];
#ifdef HAVE_MINGW
        fprintf (stderr, "read w=%u from dcwfile\n", w);
#endif
        tw += w;
        double x = w;
        bal->colperm[2*r]=w;
        bal->colperm[2*r+1]=r;
        s1 += x;
        s2 += x * x;
    }
    double avg = s1 / bal->h->ncols;
    double sdev = sqrt(s2 / bal->h->ncols - avg*avg);
    t_cw += wct_seconds();

    if (bal->h->ncoeffs) {
        if (tw != bal->h->ncoeffs) {
            fprintf(stderr, "Inconsistency in number of coefficients\n"
                    "From %s: %"PRIu64", from file sizes; %"PRIu64"\n",
                    cwfile, tw, bal->h->ncoeffs);
            fprintf(stderr, "Maybe use the --withcoeffs option for DL matrices ?\n");
            exit(1);
        }
    } else {
        bal->h->ncoeffs = tw;
        fprintf(stderr, "%"PRIu64" coefficients counted\n", tw);
    }

    fprintf(stderr, "%"PRIu32" cols ; avg %.1f sdev %.1f [scan time %.1f s]\n",
            bal->h->ncols, avg, sdev, t_cw);

    if (display_correlation) {
        size_t n = maxdim + padding;
        FILE * frw = fopen(rwfile, "rb");
        if (frw == NULL) { perror(rwfile); exit(1); }
        uint32_t * rowweights = malloc(n * sizeof(uint32_t));
        memset(rowweights, 0, n * sizeof(uint32_t));
        double t_rw;
        t_rw = -wct_seconds();
        size_t nr = fread32_little(rowweights, bal->h->nrows, frw);
        t_rw += wct_seconds();
        fclose(frw);
        if (nr < bal->h->nrows) {
            fprintf(stderr, "%s: short row count\n", rwfile);
            exit(1);
        }
        double rs1 = 0;
        double rs2 = 0;
        double rc_plain = 0;
        double rc_decorr = 0;
        for(size_t r = 0 ; r < n ; r++) {
            double x = rowweights[r];
            rs1 += x;
            rs2 += x * x;
            rc_plain += x * colweights[r];
            rc_decorr += x * bal->colperm[2*r];
        }
        double ravg = rs1 / n;
        double rsdev = sqrt(rs2 / n - ravg*ravg);
        double cavg = s1 / n;
        double csdev = sqrt(s2 / n - cavg*cavg);
        double pcov = rc_plain / n - ravg * cavg;
        double pcorr = pcov / csdev / rsdev;
        double dcov = rc_decorr / n - ravg * cavg;
        double dcorr = dcov / csdev / rsdev;

        fprintf(stderr, "%"PRIu32" rows ; avg %.1f sdev %.1f [scan time %.1f s]\n",
                bal->h->nrows, ravg, rsdev, t_rw);
        fprintf(stderr, "row-column correlation coefficient is %.4f\n",
                pcorr);
        fprintf(stderr, "row-column correlation coefficient after decorrelation is %.4f\n",
                dcorr);
    }
    free(colweights);

    struct slice * h = shuffle_rtable("col", bal->colperm, maxdim + padding, bal->h->nv);

    /* we can now make the column permutation tidier */
    bal->colperm = realloc(bal->colperm, (maxdim + padding) * sizeof(uint32_t));
    for(int ii = 0 ; ii < nv ; ii++) {
        const struct slice * r = &(h[ii]);
        memcpy(bal->colperm + r->i0, r->r, r->nrows * sizeof(uint32_t));
    }
    free_slices(h, bal->h->nv);

    bal->h->flags = FLAG_PADDING|FLAG_COLPERM;
    if (shuffled_product) bal->h->flags |= FLAG_SHUFFLED_MUL;

    bal->h->flags |= FLAG_REPLICATE;

    balancing_finalize(bal);

    balancing_write(bal, mfile, param_list_lookup_string(pl, "out"));
    balancing_clear(bal);

    param_list_clear(pl);

    return 0;
}
