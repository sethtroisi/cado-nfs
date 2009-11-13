#define _GNU_SOURCE     /* asprintf */
#define _DARWIN_C_SOURCE     /* asprintf */
#define _POSIX_C_SOURCE 200112L /* fileno */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "select_mpi.h"

#include "utils.h"
#include "mf.h"

int quiet;
int master_rank;
int extra_master;
int rank, size;
int gsize;
const char * mfile;
const char * bfile;

/* The number of send queues allocated on the master is 2 * njobs */
size_t queue_bytes = 1 << 20;
#define queue_words (queue_bytes / sizeof(uint32_t))

FILE * matrix;

/* filled completely only for the master job. Other jobs get only the
 * header. */
balancing bal;
uint32_t row_block0;
uint32_t col_block0;
uint32_t row_cellbase;
uint32_t col_cellbase;

int mcheck;     // for firing off only the master job.

struct send_queue_s {
    uint32_t * buf[2];
    int cur;
    time_t pending;
    time_t waited;
    time_t total_io_wct;        // including offloaded time.
    int tag;
    size_t pos;
    size_t total;
    size_t sent_rows;
    size_t exp_rows;
    MPI_Request req;
};
typedef struct send_queue_s send_queue[1];
typedef struct send_queue_s * send_queue_ptr;

/* slaves only */
int my_i, my_j;
uint32_t my_nrows;
uint32_t my_ncols;
size_t expected_localsize;

char * tmatfile = NULL;     /* This stores the transposed matrix file */
uint32_t * rq;
uint32_t * rq_ptr;

/* master only */
send_queue * sq;
uint32_t * sq_pool;     // pointer to be freed.
uint32_t * fw_rowperm;
uint32_t * fw_colperm;
uint32_t * ioq;

int which_node_has_row(uint32_t rnum)
{
    int b;
    if (rnum < row_block0) {
        b = rnum / (row_cellbase + 1);
    } else {
        int q = bal->h->nrows % bal->h->nh;
        b = (rnum - q) / row_cellbase;
    }
    b *= bal->h->nv;
    ASSERT(b < gsize);
    return b;
}
int which_node_has_col(uint32_t cnum)
{
    int b;
    if (cnum < col_block0) {
        b = cnum / (col_cellbase + 1);
    } else {
        int q = bal->h->ncols % bal->h->nv;
        b = (cnum - q) / col_cellbase;
    }
    ASSERT(b < gsize);
    return b;
}



void sq_send_mandatory(int i)
{
    send_queue_ptr s = sq[i];
    time_t t0 = s->pending;
    time_t t = time(NULL);
    if (t0 != 0) {
        time_t t1 = t;
        MPI_Wait(&s->req, MPI_STATUS_IGNORE);
        t = time(NULL);
        s->waited += t - t1;
        s->total_io_wct += t - t0;
    }
    s->pending = t;
    // fprintf(stderr, "send #%d -> %d\n", s->tag, i);
    MPI_Isend(s->buf[s->cur], queue_bytes, MPI_BYTE, i, s->tag, MPI_COMM_WORLD, &s->req);
    s->tag++;
    s->cur ^= 1;
    s->pos = 0;
}

int sq_send_if_possible(int i)
{
    send_queue_ptr s = sq[i];
    time_t t0 = s->pending;
    time_t t = time(NULL);
    if (t0 != 0) {
        /* If there's a request pending, then we'll do the send later. */
        return 0;
    }
    s->pending = t;
    // fprintf(stderr, "send #%d -> %d\n", s->tag, i);
    MPI_Isend(s->buf[s->cur], queue_bytes, MPI_BYTE, i, s->tag, MPI_COMM_WORLD, &s->req);
    s->tag++;
    s->cur ^= 1;
    s->pos = 0;
    return 1;
}


void read_bfile()
{
    balancing_init(bal);
    if (rank == master_rank || mcheck) {
        balancing_read(bal, bfile);
    }
    if (mcheck) {
        gsize = bal->h->nh * bal->h->nv;
    } else {
        if (rank == master_rank && !mcheck) {
            if (bal->h->nh * bal->h->nv != (uint32_t) gsize) {
                fprintf(stderr, "Error: to use %s, one needs %d processes, or %d with --extra-master (got %d here)\n", bfile, bal->h->nh * bal->h->nv, bal->h->nh * bal->h->nv+1, size);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }
    MPI_Bcast(bal->h, sizeof(balancing_header), MPI_BYTE, master_rank, MPI_COMM_WORLD);
    row_cellbase = bal->h->nrows / bal->h->nh;
    col_cellbase = bal->h->ncols / bal->h->nv;
    row_block0 = (row_cellbase + 1) * (bal->h->nrows % bal->h->nh);
    col_block0 = (col_cellbase + 1) * (bal->h->ncols % bal->h->nv);
}

void set_slave_variables()
{
    if (mcheck) return;
    if (rank >= gsize)
        return;

    /* Of course this does not make sense to come here if we have a
     * separate master job.  */
    my_j = rank % bal->h->nv;
    my_i = rank / bal->h->nv;
    my_ncols = bal->h->ncols / bal->h->nv + ((uint32_t) my_j < (bal->h->ncols % bal->h->nv));
    my_nrows = bal->h->nrows / bal->h->nh + ((uint32_t) my_i < (bal->h->nrows % bal->h->nh));
    char * suffix;
    int rc = asprintf(&suffix, "tr.h%dv%d", my_i, my_j);
    ASSERT_ALWAYS(rc>=0);
    tmatfile = build_mat_auxfile(mfile, suffix, ".bin");
    free(suffix);
    expected_localsize = my_nrows + my_ncols + bal->h->ncoeffs / gsize;
    expected_localsize += expected_localsize / 10;
    expected_localsize += queue_words;
    printf("Job %3d (%2d,%2d) expects"
            " %"PRIu32" rows %"PRIu32" cols."
            // " Datafile %s."
            " Allocating %zu MB\n",
            rank, my_i, my_j, my_nrows, my_ncols,
            // tmatfile,
            expected_localsize * sizeof(uint32_t) >> 20);
    rq = malloc(expected_localsize * sizeof(uint32_t));
    rq_ptr = rq;
}

void set_master_variables()
{
    if (rank != master_rank && !mcheck)
        return;

    sq = malloc(gsize * sizeof(send_queue));
    sq_pool = malloc(gsize * 2 * queue_bytes);
    for(int i = 0 ; i < gsize ; i++) {
        memset(sq[i], 0, sizeof(send_queue));
        sq[i]->buf[0] = sq_pool + queue_words * (2 * i);
        sq[i]->buf[1] = sq_pool + queue_words * (2 * i + 1);
        // sq[i]->pending = 0;
        // sq[i]->cur = 0;
        int h = bal->h->nh;
        int v = bal->h->nv;
        sq[i]->exp_rows = bal->h->nrows / h + ((uint32_t) (i / v) < (bal->h->nrows % h));
    }

    ioq = malloc(queue_bytes);

    uint32_t maxdim = MAX(bal->h->nrows, bal->h->ncols);

    if (bal->h->flags && FLAG_REPLICATE) {
        fw_rowperm = fw_colperm = malloc(maxdim * sizeof(uint32_t));
        memset(fw_rowperm, -1, maxdim * sizeof(uint32_t));
        for(uint32_t i = 0 ; i < maxdim ; i++) {
            uint32_t j = (bal->colperm ? bal->colperm : bal->rowperm)[i];
            ASSERT(fw_colperm[j] == UINT32_MAX);
            fw_colperm[j]=i;
        }
    } else if (bal->h->flags & FLAG_PADDING) {
        ASSERT_ALWAYS(bal->h->flags & FLAG_COLPERM);
        ASSERT_ALWAYS(bal->h->flags & FLAG_ROWPERM);
        fw_rowperm = malloc(maxdim * sizeof(uint32_t));
        fw_colperm = malloc(maxdim * sizeof(uint32_t));
        memset(fw_rowperm, -1, maxdim * sizeof(uint32_t));
        memset(fw_colperm, -1, maxdim * sizeof(uint32_t));
        for(uint32_t i = 0 ; i < maxdim ; i++) {
            ASSERT(fw_colperm[bal->colperm[i]] == UINT32_MAX);
            ASSERT(fw_rowperm[bal->rowperm[i]] == UINT32_MAX);
            fw_colperm[bal->colperm[i]]=i;
            fw_rowperm[bal->rowperm[i]]=i;
        }
    } else {
        ASSERT_ALWAYS(bal->h->flags & FLAG_COLPERM);
        fw_colperm = malloc(bal->h->ncols * sizeof(uint32_t));
        memset(fw_colperm, -1, bal->h->ncols * sizeof(uint32_t));
        for(uint32_t i = 0 ; i < bal->h->nrows ; i++) {
            ASSERT(fw_colperm[bal->colperm[i]] == UINT32_MAX);
            fw_colperm[bal->colperm[i]]=i;
        }

        ASSERT_ALWAYS(bal->h->flags & FLAG_ROWPERM);
        fw_rowperm = malloc(bal->h->nrows * sizeof(uint32_t));
        memset(fw_rowperm, -1, bal->h->nrows * sizeof(uint32_t));
        for(uint32_t i = 0 ; i < bal->h->ncols ; i++) {
            ASSERT(fw_rowperm[bal->colperm[i]] == UINT32_MAX);
            fw_rowperm[bal->rowperm[i]]=i;
        }
    }
    
    /* one more check. The cost is tiny compared to what we do in other
     * parts of the code. */
    printf("Consistency check..."); fflush(stdout);
    uint32_t * ttab = malloc(bal->h->nh * sizeof(uint32_t));
    memset(ttab, 0, bal->h->nh * sizeof(uint32_t));
    for(uint32_t j = 0 ; j < bal->h->nrows ; j++) {
        ttab[which_node_has_row(fw_rowperm[j]) / bal->h->nv]++;
    }
    for(uint32_t k = 0 ; k < bal->h->nh ; k++) {
        ASSERT_ALWAYS(ttab[k] == sq[k * bal->h->nv]->exp_rows);
    }
    printf("ok\n");
}

void clear_master_variables()
{
    if (rank != master_rank)
        return;
    if (bal->h->flags && FLAG_REPLICATE) {
        free(fw_rowperm);
    } else {
        free(fw_rowperm);
        free(fw_colperm);
    }
    fw_rowperm = fw_colperm = NULL;
    free(sq_pool);
    free(sq);
    free(ioq);
}

void clear_slave_variables()
{
    free(rq); rq = NULL;
    free(tmatfile); tmatfile = NULL;
}

void sq_add_value(int i, uint32_t v)
{
    if (i == master_rank) {
        ASSERT_ALWAYS((size_t) (rq_ptr - rq) < expected_localsize);
        *rq_ptr++ = v;
        return;
    }

    send_queue_ptr s = sq[i];
    if (s->pos == queue_words) {
        /* At this point, if a request was pending (thus on the other
         * buffer), it is time to wait for its completion, and once it
         * has completed, initiate a second one.
         */
        sq_send_mandatory(i);
    }
    ASSERT(s->pos < queue_words);
    s->buf[s->cur][s->pos] = v;
    ++s->total;
    if (++s->pos == queue_words) {
        /* if there's a request pending, we don't need to block yet
         * before being able to send our payload */
        sq_send_if_possible(i);
    }
}

void master_loop()
{
    FILE * fmat;
    fmat = fopen(mfile, "r");
    if (fmat == NULL) {
        perror(mfile);
        exit(1);
    }
    struct stat sbuf[1];
    fstat(fileno(fmat), sbuf);
    size_t esz = (bal->h->ncoeffs + bal->h->nrows) * sizeof(uint32_t);
    if ((size_t) sbuf->st_size != esz) {
        fprintf(stderr, "%s: expected size %zu, not %zu\n",
                mfile, esz, sbuf->st_size);
        exit(1);
    }

    size_t pos = 0;
    uint32_t crow_togo = 0;
    uint32_t r = 0;
    int noderow = 0;
    time_t t0 = time(NULL);
    time_t t1 = t0 + 1;
    uint32_t maxdim = MAX(bal->h->nrows, bal->h->ncols);

    for( ; ; ) {
        size_t n = fread(ioq, sizeof(uint32_t), queue_words, fmat);
        if (n == 0) {
            if (ferror(fmat)) {
                perror(mfile);
                exit(1);
            } else {
                if (pos != (size_t) sbuf->st_size) {
                    fprintf(stderr, "Ended at wrong position in %s\n", mfile);
                    exit(1);
                }
            }
            break;
        }
        pos += n * sizeof(uint32_t);
        for(size_t s = 0 ; s < n ; s++) {
            if (crow_togo == 0) {
                crow_togo = ioq[s];
                noderow = which_node_has_row(fw_rowperm[r]);
                /* Send a new row info _only_ to the nodes on the
                 * corresponding row in the process grid ! */
                for(int i = 0 ; i < (int) bal->h->nv ; i++) {
                    sq_add_value(noderow + i, UINT32_MAX);
                    sq[noderow + i]->sent_rows++;
                    if (sq[noderow + i]->sent_rows > sq[noderow + i]->exp_rows) {
                        fprintf(stderr, "Number of rows sent to node %d (%zu) exceeds expected value %zu\n", noderow+i, sq[noderow + i]->sent_rows, sq[noderow + i]->exp_rows);
                        // abort();
                    }
                }
                /* This advances, even if the row is not complete. */
                r++;
            } else {
                uint32_t c = fw_colperm[ioq[s]];
                int nodecol = which_node_has_col(c);
                sq_add_value(noderow + nodecol, c);
                --crow_togo;
            }
        }
        time_t t = time(NULL);
        if (t >= t1) {
            t1 = t + 1;
            fprintf(stderr,
                    "Read %zu MB, %"PRIu32" rows in %d s ; %.0f MB/s  \n",
                    pos>>20, r, (int) (t-t0), ((double)pos/(t1-t0))*1.0e-6);
        }
    }
    ASSERT(r == bal->h->nrows);
    fprintf(stderr, "Master loop finished\n");
    for(int i = 0 ; i < (int) gsize ; i++) {
        /* complete from nrows to maxdim if necessary ! */
        if (bal->h->flags && FLAG_PADDING) {
            for( ; r < maxdim ; r++) {
                sq_add_value(i, UINT32_MAX);
            }
        }
        sq_add_value(i, UINT32_MAX);
        sq_send_mandatory(i);
    }
    fprintf(stderr, "Final sends completed, everybody should finish.\n");
    for(int i = 0 ; i < (int) gsize ; i++) {
        printf("[%d] sent %zu rows, %zu coeffs\n",
                i, sq[i]->sent_rows, sq[i]->total - sq[i]->sent_rows);
    }
}

void slave_loop()
{
    /* If the master also plays the role of a slave, then everything
     * is handled as well from the master loop, not from here.
     */
    ASSERT(rank != master_rank);

    int tag = 0;
    uint32_t r = 0;
    uint32_t * rq_marker = NULL;
    uint64_t tw=0;
    double s1 = 0;
    double s2 = 0;
    for( ; r < my_nrows ; tag++ ) {
        // fprintf(stderr, "[%d] expects send #%d at offset %td\n",
                // rank, tag, rq_ptr - rq);
#if 0
        fprintf(stderr, "[%d] %"PRIu32" rows so far (%.1f%%),"
                " %"PRIu64" coeffs (%.1f%%)\n",
                rank,
                r,  100.0 * (double) r / my_nrows,
                tw, 100.0 * (double) (tw+r) / expected_localsize);
#endif
#if 0
        if (rq_marker) {
            ASSERT_ALWAYS((size_t) (rq_marker - rq) == tw + r);
        } else {
            ASSERT_ALWAYS(rq_ptr == rq);
        }
        if (!((size_t)(rq_ptr + queue_words - rq) < expected_localsize)) {
            fprintf(stderr, "argh on rank %d\n", rank);
            fprintf(stderr, "tw %"PRIu64" r %"PRIu32" e %zu qw %zu delta %td\n",
                    tw, r, expected_localsize, queue_words, rq_ptr - rq);
            ASSERT_ALWAYS(0);
        }
#endif
        ASSERT_ALWAYS((size_t)(rq_ptr + queue_words - rq) < expected_localsize);
        MPI_Recv(rq_ptr, queue_bytes, MPI_BYTE,
                master_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // fprintf(stderr, "Recv #%d -> %d\n", tag, rank);
        // fflush(stderr);
        // fflush(stdout);
        uint32_t * rq_tail = rq_ptr + queue_words;
        for( ; r < my_nrows && rq_ptr != rq_tail ; rq_ptr++) {
            uint32_t x = *rq_ptr;
            if (x != UINT32_MAX)
                continue;
            if (rq_marker) {
                uint32_t w = rq_ptr - rq_marker - 1;
                *rq_marker = w;
                s1 += w;
                s2 += ((double)w)*((double)w);
                tw += w;
                r++;
            }
            rq_marker = rq_ptr;
        }
    }
    printf("[%d] %"PRIu32" rows, %"PRIu64" coeffs.\n", rank, r, tw);
}

void all()
{
    read_bfile();

    if (rank == master_rank) {
        printf("Total %"PRIu32" rows %"PRIu32" cols %"PRIu64" coeffs\n",
                bal->h->nrows,
                bal->h->ncols,
                bal->h->ncoeffs);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    set_slave_variables();
    set_master_variables();

    //////////////////////////////////////////////////////////////////
    
    if (rank == master_rank) {
        master_loop();
    } else {
        slave_loop();
    }

    /* master must flush all pending send queues */

    MPI_Barrier(MPI_COMM_WORLD);

    /* slaves must now transpose their data and save it. */

    clear_slave_variables();
    clear_master_variables();
}

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    param_list pl;
    param_list_init(pl);
    int wild = 0;
    argv++,argc--;
    param_list_configure_knob(pl, "--quiet", &quiet);
    param_list_configure_knob(pl, "--mcheck", &mcheck);
    param_list_configure_knob(pl, "--extra-master", &extra_master);
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        if (wild == 0) { mfile = argv[0]; wild++,argv++,argc--; continue; }
        if (wild == 1) { bfile = argv[0]; wild++,argv++,argc--; continue; }
        fprintf(stderr, "Unknown option %s\n", argv[0]);
        exit(1);
    }

    param_list_parse_ulong(pl, "qsize", &queue_bytes);
    master_rank = extra_master ? size-1 : 0;
    gsize = size - extra_master;
    ASSERT_ALWAYS(queue_bytes % sizeof(uint32_t) == 0);

    if (mcheck) {
        ASSERT_ALWAYS(size == 1);
        ASSERT_ALWAYS(rank == 0);
    }

    all();

    param_list_clear(pl);
    MPI_Finalize();
}
