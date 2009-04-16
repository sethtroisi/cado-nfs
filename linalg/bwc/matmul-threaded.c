#define _GNU_SOURCE     /* asprintf */
#define _DARWIN_C_SOURCE    /* for asprintf. _ANSI_SOURCE must be undefined */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <pthread.h>

#include "matmul.h"
#include "matmul-threaded.h"
#include "abase.h"
#include "manu.h"

#include "readmat-easy.h"

/* This extension is used to distinguish between several possible
 * implementations of the product */
#define MM_EXTENSION   "-threaded"
#define MM_MAGIC_FAMILY        0xa002UL
#define MM_MAGIC_VERSION       0x1004UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)

#ifdef  __GNUC__
#define ASM_COMMENT(x)  __asm__("#\t" x "\n")
#else
#define ASM_COMMENT(x)  /**/
#endif

/* A fair part of these should become tunable */
#define MM_DENSIFY_TOLERANCE    0.5    /* Set to zero to skip dense blocks */
#define MM_DGRP_SIZE            32
#define MM_NTHREADS             4
#define MM_SGRP_SIZE            64
#define MM_DO_TRANSPOSE         1
#define MM_PREFETCH_OFFSET_1    32
#define MM_PREFETCH_OFFSET_2    16
#define MM_PREFETCH_OFFSET_3    16
#define MM_BLOCKSIZE            8


struct thread_info {
    int i;
    matmul_ptr mm;
};

struct matmul_threaded_data_s {
    /* repeat the fields from the public interface */
    // unsigned int nrows;
    // unsigned int ncols;
    unsigned int dim[2];        /* dim[0] is nrows, dim[1] is ncols. The
                                   organization of the matrix in memory
                                   also obeys the ``transposed'' flag */
    unsigned long ncoeffs;

    /* Now our private fields */
    abobj_t xab;

    struct {
        uint32_t n;
        uint32_t * q;
        uint32_t weight;
        /* This is protected by the mutex, and used for serializing the
         * writing of the data area */
        uint32_t processed;
        pthread_cond_t recheck_please;
    } dense[1];
    struct {
        uint32_t n;
        uint32_t * len;
        uint32_t tlen;
        uint32_t * q;
    } sparse[1];

    uint32_t sgrp_size;
    uint32_t nthreads;
    uint32_t blocksize;
    uint32_t transposed;
    uint32_t off1;
    uint32_t off2;
    uint32_t off3;

    pthread_mutex_t mu;
    pthread_cond_t cond_in;
    pthread_cond_t cond_out;

    /* Number of worker threads at work (not idle) */
    unsigned int working;
    unsigned int done;

    struct thread_info * thread_table;

    pthread_t * threads;
    
    void * dst;
    const void * src;
    int d;
    int iteration;
};

#define MM      ((struct matmul_threaded_data_s *) mm)

extern void matmul_threaded_mul_sub_dense(matmul_ptr mm, abt * dst, abt const * src, int d, int i);
extern void matmul_threaded_mul_sub_sparse(matmul_ptr mm, abt * dst, abt const * src, int d, int i);

void * thread_job(void * ti)
{
    matmul_ptr mm = ((struct thread_info *) ti) -> mm;
    int i = ((struct thread_info *) ti) -> i;

    pthread_mutex_lock(&MM->mu);

    MM->working++;
    pthread_cond_signal(&MM->cond_out);
    /* The mutex is released later by pthread_cond_wait */

    unsigned int res = 0;

    for( ; ; ) {
        pthread_cond_wait(&MM->cond_in,&MM->mu);
        MM->working++;
        pthread_mutex_unlock(&MM->mu);
        if (MM->iteration < 0)
            break;

        /* Do our job with mutexes released */
        matmul_threaded_mul_sub_dense(mm, MM->dst, MM->src, MM->d, i);
        matmul_threaded_mul_sub_sparse(mm, MM->dst, MM->src, MM->d, i);

        /* resume serialized behaviour */
        pthread_mutex_lock(&MM->mu);
        res = ++(MM->done);

        if (res == MM->nthreads) {
            pthread_cond_signal(&MM->cond_out);
        }
    }

    pthread_mutex_unlock(&MM->mu);

    return NULL;
}

void worker_threads_init(matmul_ptr mm)
{
    MM->threads = malloc(MM->nthreads * sizeof(pthread_t));
    MM->thread_table = malloc(MM->nthreads * sizeof(struct thread_info));
    pthread_cond_init(&MM->cond_in, NULL);
    pthread_cond_init(&MM->cond_out, NULL);
    pthread_cond_init(&MM->dense->recheck_please, NULL);
    pthread_mutex_init(&MM->mu, NULL);
    MM->working = 0;
    MM->done = 0;
    MM->iteration = 0;
    for(unsigned int i = 0 ; i < MM->nthreads ; i++) {
        int rc;
        MM->thread_table[i].i = i;
        MM->thread_table[i].mm = mm;
        void * ti = &(MM->thread_table[i]);
        rc = pthread_create(&(MM->threads[i]), NULL, &thread_job, ti);
        ASSERT_ALWAYS(rc == 0);
    }

    /* We must make sure that all threads have started properly. Since
     * we've got everything available for writing a barrier, we spare the
     * hassle of defining an extra data member for it. */
    pthread_mutex_lock(&MM->mu);
    for( ; MM->working != MM->nthreads ; ) {
        pthread_cond_wait(&MM->cond_out, &MM->mu);
    }
    pthread_mutex_unlock(&MM->mu);
}

void worker_threads_clear(matmul_ptr mm)
{
    MM->iteration = -1;
    pthread_mutex_lock(&MM->mu);
    pthread_cond_broadcast(&MM->cond_in);
    pthread_mutex_unlock(&MM->mu);
    for(unsigned int i = 0 ; i < MM->nthreads ; i++) {
        void * res;
        int rc = pthread_join(MM->threads[i], &res);
        ASSERT_ALWAYS(rc == 0);
    }
    free(MM->threads);
    free(MM->thread_table);
    pthread_mutex_destroy(&MM->mu);
    pthread_cond_destroy(&MM->cond_in);
    pthread_cond_destroy(&MM->cond_out);
}

void matmul_threaded_clear(matmul_ptr mm)
{
    worker_threads_clear(mm);

    free(MM->dense->q);
    free(MM->sparse->len);
    free(MM->sparse->q);
    free(MM);
}

static matmul_ptr matmul_threaded_init()
{
    matmul_ptr mm = malloc(sizeof(struct matmul_threaded_data_s));
    memset(mm, 0, sizeof(struct matmul_threaded_data_s));

    /* TODO: Turn these into parameters */
    /* TODO: remove it from the on-disk file, since it does not impact
     * the data itself. */
    MM->nthreads = MM_NTHREADS;


    /* TODO: Turn these into parameters */
    /* Transposing here is just an optimisation concern. The code does
     * not compute anything different depending on this value. But it
     * performs differently */
    MM->transposed = MM_DO_TRANSPOSE;
    MM->sgrp_size = MM_SGRP_SIZE;

    /* Must correspond to a L1 cache line. */
    /* XXX Must also correspond to the selected abase ! */
    MM->blocksize = MM_BLOCKSIZE;

    MM->off1 = MM_PREFETCH_OFFSET_1;
    MM->off2 = MM_PREFETCH_OFFSET_2;
    MM->off3 = MM_PREFETCH_OFFSET_3;

    /* TODO of course, once settled, generate assembly code */
    
    return mm;
}

typedef int (*sortfunc_t) (const void *, const void *);

static int uint_cmp(unsigned int * a, unsigned int * b)
{
    if (*a < *b) return -1;
    else if (*b < *a) return 1;
    return 0;
}

void matmul_threaded_blocks_info(matmul_ptr mm)
{
    printf("// %u dense blocks of %d rows\n",
            MM->dense->n, MM_DGRP_SIZE);
    unsigned long sparse_coeffs = MM->ncoeffs - MM->dense->weight;
    unsigned long padding = MM->sparse->tlen - sparse_coeffs;
    printf("// %u sparse blocks of %d rows ; %lu coeffs + %lu padding (+%.2f%%)\n",
            MM->sparse->n, MM->sgrp_size,
            sparse_coeffs, padding,
            (100.0 * padding / sparse_coeffs));
}

matmul_ptr matmul_threaded_build(abobj_ptr xx MAYBE_UNUSED, const char * filename)
{
    matmul_ptr mm = matmul_threaded_init();

    abobj_init_set(MM->xab, xx);

    /* All fields are initially set to zero */

    /* We give some notations and example values.  */

    // uint32_t n = MM->nthreads;      /* e.g. 4 */
    uint32_t g = MM->sgrp_size;    /* e.g. 32 */
    // uint32_t b = MM->blocksize;     /* e.g. 16 */
    ASSERT_ALWAYS(g % (MM->nthreads*MM->blocksize) == 0);
    size_t sparse_readahead = MM->nthreads * MM->off2 * MM->blocksize;

    /* Okay, this matrix reading stage is really a memory hog. But it's
     * done only once, and we can even consider doing it only once, so we
     * practically don't care. */
    /* dims of the transposed matrix */
    unsigned int nrows_t;
    unsigned int ncols_t;
    uint32_t * data;
    if (MM->transposed) {
        data = read_easy_transposed(filename, &nrows_t, &ncols_t);
    } else {
        data = read_easy(filename, &nrows_t, &ncols_t);
    }
    MM->dim[ MM->transposed] = nrows_t;
    MM->dim[!MM->transposed] = ncols_t;

    /* count dense and sparse blocks. Since one entry in the sparse part
     * occupies a 32-bit location, we consider that it's preferrable to
     * use a dense mode when it's no heavier. For 32-bit consecutive
     * rows, we then store a sequence of N 32-bit integers, whose bit
     * indicates entries in the N columns of the block.
     *
     * 32 is somehow a constant here. But we may allow some waste
     * storage, say up to a factor of 2. That's called MM_DENSIFY_TOLERANCE
     */

    uint32_t * ptr = data;


    unsigned int i = 0;
    /* Count dense blocks */
    for( ; i < nrows_t ; MM->dense->n++) {
        unsigned int di;
        unsigned int weight = 0;
        uint32_t * q = ptr;
        for(di = 0 ; di < MM_DGRP_SIZE && i+di < nrows_t ; di++) {
            weight += *q;
            q += 1 + *q;
        }
        int dense = weight * MM_DENSIFY_TOLERANCE >= ncols_t;
        const char * rc[2] = { "row","col", };
        printf("%ss %u..%u ; total %u coeffs = n%ss(%u) * %.2f ;%s dense\n",
                rc[MM->transposed],
                i, i+di-1, weight,
                rc[!MM->transposed],
                MM->dim[!MM->transposed],
                (double) weight / ncols_t,
                dense ? "" : " not");
        if (!dense)
            break;
        MM->ncoeffs += weight;
        MM->dense->weight += weight;
        ptr = q;
        i += di;
    }
    /* Count sparse blocks */
    for( ; i < nrows_t ; MM->sparse->n++) {
        unsigned long slice_coeffs = 0;
        uint32_t mlen = 0;
        unsigned int di;
        for(di = 0 ; di < g && i + di < nrows_t ; di++) {
            if (*ptr >= mlen)
                mlen = *ptr;
            slice_coeffs += *ptr;
            ptr += 1 + *ptr;
        }
        MM->ncoeffs += slice_coeffs;
        MM->sparse->tlen += g * mlen;
        i += di;
    }

    matmul_threaded_blocks_info(mm);

    size_t dblock_stride = next_multiple_of_powerof2(ncols_t, MM->blocksize);

    MM->dense->q = malloc(MM->dense->n * dblock_stride * sizeof(uint32_t));
    MM->sparse->q = malloc((MM->sparse->tlen + sparse_readahead) * sizeof(uint32_t));
    MM->sparse->len = malloc(MM->sparse->n * sizeof(uint32_t));

    /* reset pointer to beginning, start to fill dense blocks */
    ptr = data;
    i = 0;
    memset(MM->dense->q, 0, MM->dense->n * dblock_stride * sizeof(uint32_t));
    for(unsigned int blocknum = 0 ; blocknum < MM->dense->n ; blocknum++) {
        uint32_t * q = MM->dense->q + blocknum * dblock_stride;
        unsigned int di;
        for(di = 0 ; di < MM_DGRP_SIZE && i + di < nrows_t ; di++) {
            uint32_t rlen = *ptr++;
            for(unsigned int j = 0 ; j < rlen ; j++) {
                q[*ptr++] ^= ((uint32_t) 1) << di;
            }
        }
        i += di;
    }

    /* Now proceed with sparse blocks */
    uint32_t * sq = MM->sparse->q;
    /* First fill padding data */
    for(unsigned int j = 0 ; j < MM->sparse->tlen ; j++)
        sq[j] = ncols_t;
    for(unsigned int blocknum = 0 ; blocknum < MM->sparse->n ; blocknum++) {
        unsigned int di;
        /* We'll discover the max size in the strip once again */
        uint32_t mlen = 0;

        /* striding for the entries in one row is g==sgrp_size.
         * Entry g*l+k corresponds to index l in row k.  A
         * complete block of g indices is covered by the n threads in
         * g/(nb) blocks of b rows (g/nb blocks from each thread).
         */

        for(di = 0 ; di < g ; di++) {
            uint32_t rlen;
            
            if (i + di < nrows_t) {
                rlen = *ptr++;
                qsort(ptr, rlen, sizeof(uint32_t), (sortfunc_t) uint_cmp);
            } else {
                rlen = 0;
            }

            ASSERT(sq-MM->sparse->q+rlen*g <= MM->sparse->tlen);

            for(unsigned int j = 0 ; j < rlen ; j++) sq[di + j * g] = *ptr++;

            if (rlen >= mlen) mlen = rlen;
        }

        MM->sparse->len[blocknum] = mlen * g;
        i += di;
        sq += mlen * g;
    }

    memset(sq, 0, sparse_readahead * sizeof(uint32_t));

    free(data);

    worker_threads_init(mm);

    return mm;
}

matmul_ptr matmul_threaded_reload_cache(abobj_ptr xx MAYBE_UNUSED, const char * filename)
{
    char * base;
    FILE * f;

    {
        int rc = asprintf(&base, "%s" MM_EXTENSION ".bin", filename);
        FATAL_ERROR_CHECK(rc < 0, "out of memory");
    }

    f = fopen(base, "r");

    if (f == NULL) {
        // fprintf(stderr, "fopen(%s): %s\n", base, strerror(errno));
        fprintf(stderr, "no cache file %s\n", base);
        free(base);
        return NULL;
    }
    free(base);

    size_t rc;
    matmul_ptr mm = matmul_threaded_init();

    abobj_init_set(MM->xab, xx);

    unsigned long magic_check;
    rc = fread(&magic_check, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    FATAL_ERROR_CHECK(magic_check != MM_MAGIC,
            "Wrong magic in cached matrix file");
    rc = fread(MM->dim, sizeof(unsigned int), 2, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->ncoeffs, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->dense->weight, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->sgrp_size, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->transposed, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->dense->n, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->sparse->n, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->sparse->tlen, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");

    // unsigned int nrows_t = MM->dim[ MM->transposed];
    unsigned int ncols_t = MM->dim[!MM->transposed];
    size_t dblock_stride = next_multiple_of_powerof2(ncols_t, MM->blocksize);

    /* sparse_readahead depends on data which is initialized by us, not
     * recovered from the data file */
    size_t sparse_readahead = MM->nthreads * MM->off2 * MM->blocksize;

    /* Now read the interesting blocks */
    MM->dense->q = malloc(MM->dense->n * dblock_stride * sizeof(uint32_t));
    memset(MM->dense->q, 0, MM->dense->n * dblock_stride * sizeof(uint32_t));
    MM->sparse->q = malloc((MM->sparse->tlen + sparse_readahead) * sizeof(uint32_t));
    MM->sparse->len = malloc(MM->sparse->n * sizeof(uint32_t));

    size_t len;
    for(unsigned int i = 0 ; i < MM->dense->n ; i++) {
        len = ncols_t;
        rc = fread(MM->dense->q + i * dblock_stride, sizeof(uint32_t), len, f);
        FATAL_ERROR_CHECK(rc < len, "Short read from cached matrix file");
    }
    len = MM->sparse->n;
    rc = fread(MM->sparse->len, sizeof(uint32_t), len, f);
    FATAL_ERROR_CHECK(rc < len, "Short read from cached matrix file");
    len = MM->sparse->tlen;
    rc = fread(MM->sparse->q, sizeof(uint32_t), len, f);
    FATAL_ERROR_CHECK(rc < len, "Short read from cached matrix file");

    /* fill up the readahead zone properly */
    memset(MM->sparse->q + len, 0, sparse_readahead * sizeof(uint32_t));

    fclose(f);

    matmul_threaded_blocks_info(mm);

    worker_threads_init(mm);

    return mm;
}

void matmul_threaded_save_cache(matmul_ptr mm, const char * filename)
{
    char * base;
    FILE * f;

    {
        int rc = asprintf(&base, "%s" MM_EXTENSION ".bin", filename);
        FATAL_ERROR_CHECK(rc < 0, "out of memory");
    }

    f = fopen(base, "w");
    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", base, strerror(errno));
        free(base);
        exit(1);
    }
    free(base);

    unsigned long magic = MM_MAGIC;
    size_t rc;

    rc = fwrite(&magic, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(MM->dim, sizeof(unsigned int), 2, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->ncoeffs, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->dense->weight, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->sgrp_size, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->transposed, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->dense->n, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->sparse->n, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->sparse->tlen, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");

    // unsigned int nrows_t = MM->dim[ MM->transposed];
    unsigned int ncols_t = MM->dim[!MM->transposed];
    size_t dblock_stride = next_multiple_of_powerof2(ncols_t, MM->blocksize);

    size_t len;
    for(unsigned int i = 0 ; i < MM->dense->n ; i++) {
        len = ncols_t;
        rc = fwrite(MM->dense->q + i * dblock_stride, sizeof(uint32_t), len, f);
        FATAL_ERROR_CHECK(rc < len, "Short write to cached matrix file");
    }
    len = MM->sparse->n;
    rc = fwrite(MM->sparse->len, sizeof(uint32_t), len, f);
    FATAL_ERROR_CHECK(rc < len, "Short write to cached matrix file");
    len = MM->sparse->tlen;
    rc = fwrite(MM->sparse->q, sizeof(uint32_t), len, f);
    FATAL_ERROR_CHECK(rc < len, "Short write to cached matrix file");

    fclose(f);
}

void matmul_threaded_mul_sub_dense(matmul_ptr mm, abt * dst, abt const * src, int d, int throff)
{
    if (MM->dense->n == 0) {
        /* We have nothing to do. In this special case, the task of
         * clearing the destination area is handed over to the sparse
         * block case */
        return;
    }

    ASM_COMMENT("dense multiplication code");
    abobj_ptr x = MM->xab;
    uint32_t nthr = MM->nthreads;      /* e.g. 4 */
    uint32_t b = MM->blocksize;     /* e.g. 16 */
    uint32_t off3 = MM->off3;

    uint32_t * A = MM->dense->q + throff * b;
    uint32_t i0 = throff * b;

    abt * rx = abinitf(x, 2);
    abzero(x, rx, 1);

    abt * dst0 = dst;
    abt const * src0 = src;

    if (d == !MM->transposed) {
        /* t == 1 is well suited for d == 0. In this case, dst is
         * progressively filled. The same holds for t == 0 and d == 1 */

        /* We expect to have very few dense blocks, so it's acceptable to
         * have temporary storage for all of them on the stack */
        abt * tdst = abinitf(x, MM_DGRP_SIZE * MM->dense->n);
        abt * rb = tdst;
        for(unsigned int blocknum = 0 ; blocknum < MM->dense->n ; blocknum++) {
            src = src0 + aboffset(x, i0);
            abzero(x, rb, MM_DGRP_SIZE);
            ASM_COMMENT("critical loop");
            for(uint32_t i = i0 ; i < MM->dim[d] ; i += nthr * b) {
                for(uint32_t j = 0 ; j < b ; j++) {
                    __builtin_prefetch(A + nthr * b * off3, 0);
                    __builtin_prefetch(src + aboffset(x, nthr * b * off3), 0);
                    uint32_t rA = *A++;
                    abcopy(x, rx + aboffset(x, 1), src, 1);
                    src += aboffset(x, 1);
                    for(unsigned int k = 0 ; k < MM_DGRP_SIZE ; k++) {
                        abadd(x, rb + aboffset(x, k), rx + aboffset(x, rA&1));
                        rA >>= 1;
                    }
                }
                A += (nthr-1) * b;
                src += aboffset(x, (nthr-1) * b);
            }
            rb += aboffset(x, MM_DGRP_SIZE);
            ASM_COMMENT("end of critical loop");
        }

        /* Now the painful section. Experimentally, it does not seem
         * expensive to adopt this crude approach. */
        pthread_mutex_lock(&MM->mu); 
        if (MM->dense->processed++ == 0) {
            abcopy(x, dst, tdst, MM_DGRP_SIZE * MM->dense->n);
        } else {
            for(unsigned int k = 0 ; k < MM_DGRP_SIZE * MM->dense->n ; k++) {
                abadd(x, dst + k, tdst + k);
            }
        }
        pthread_mutex_unlock(&MM->mu);

#if 0
        /* There could be a way to seemingly avoid requiring a mutex
         * here, although in the end it probably wouldn't change much.
         * However, in any case we have to zero out the dst area in a
         * locked manner beforehand, and this implies serialization as
         * well.
         *
         * For reference, plan B uses the atomic xor operation
         * __sync_xor_and_fetch, which is a gcc builtin (corresponds to
         * using the lock prefix on x86*).  We're cheating with the abase
         * interface here. A lot.  Because depending on the active abase,
         * we may not be allowed to use __sync_xor_and_fetch. Thus we
         * have to cast everything to unsigned longs */

        unsigned int dsize_ulongs;
        dsize_ulongs = aboffset(x, MM_DGRP_SIZE);
        dsize_ulongs *= MM->dense->n;
        dsize_ulongs *= sizeof(abt);
        ASSERT(dsize_ulongs % sizeof(unsigned long) == 0);
        dsize_ulongs /= sizeof(unsigned long);

        for(unsigned int k = 0 ; k < dsize_ulongs ; k++) {
            __sync_xor_and_fetch(
                    ((unsigned long*)dst) + k,
                    ((unsigned long*)tdst) [k]);
        }
#endif
        abclearf(x, tdst, MM_DGRP_SIZE * MM->dense->n);
    } else {
        for(uint32_t i = i0 ; i < MM->dim[!d] ; i += nthr * b) {
            abzero(x, dst0 + aboffset(x, i), b);
        }
        for(unsigned int blocknum = 0 ; blocknum < MM->dense->n ; blocknum++) {
            dst = dst0 + aboffset(x, i0);
            ASM_COMMENT("critical loop");
            for(uint32_t i = i0 ; i < MM->dim[!d] ; i += nthr * b) {
                for(uint32_t j = 0 ; j < b ; j++) {
                    __builtin_prefetch(A + nthr * b * off3, 0);
                    __builtin_prefetch(dst + aboffset(x, nthr * b * off3), 1);
                    uint32_t rA = *A++;
                    abzero(x, rx + aboffset(x, 1), 1);
                    for(unsigned int k = 0 ; k < MM_DGRP_SIZE ; k++) {
                        abadd(x, rx + aboffset(x, rA&1), src + aboffset(x, k));
                        rA >>= 1;
                    }
                    abadd(x, dst, rx + aboffset(x, 1));
                    dst += aboffset(x, 1);
                }
                A += (nthr-1) * b;
                dst += aboffset(x, (nthr-1) * b);
            }
            src += aboffset(x, MM_DGRP_SIZE);
            ASM_COMMENT("end of critical loop");
        }
        pthread_mutex_lock(&MM->mu); 
        /* It's a pity, but if some other thread is already doing the
         * sparse blocks at this point (it can't be anything else than
         * thread 0), we must be sure that everybody here is done with
         * the dense blocks.
         */
        if (++MM->dense->processed == MM->nthreads) {
            pthread_cond_signal(&MM->dense->recheck_please);
        }
        pthread_mutex_unlock(&MM->mu);
    }
    abclearf(x, rx, 2);
    ASM_COMMENT("end of dense multiplication code");
}

void matmul_threaded_mul_sub_sparse(matmul_ptr mm, abt * dst, abt const * src, int d, int throff)
{
    ASM_COMMENT("sparse multiplication code");
    abobj_ptr x = MM->xab;
    uint32_t nthr MAYBE_UNUSED = MM->nthreads;
    uint32_t g = MM->sgrp_size;
    uint32_t b = MM->blocksize;
    uint32_t off1 = MM->off1;
    uint32_t off2 = MM->off2;

    uint32_t * A = MM->sparse->q + throff * b;
    uint32_t i0 = throff * b;


    /* If d == 0, output has ncols == dim[1] entries
     *             input has nrows == dim[0] entries
     */

    /*
    abt * dst0 = dst;
    abt const * src0 = src;
    unsigned int ra[2];
    ra[0] = MM->dim[0];
    matmul_threaded_aux(mm, MATMUL_AUX_GET_READAHEAD, &(ra[0]));
    ra[1] = MM->dim[1];
    matmul_threaded_aux(mm, MATMUL_AUX_GET_READAHEAD, &(ra[1]));
    */

    if (d == !MM->transposed) {
        dst += aboffset(x, MM->dense->n * MM_DGRP_SIZE);
#if 1
#define b       MM_BLOCKSIZE
#define g       MM_SGRP_SIZE
#define nthr    MM_NTHREADS
#define off1    MM_PREFETCH_OFFSET_1
#define off2    MM_PREFETCH_OFFSET_2
#endif
        for(unsigned int blocknum = 0 ; blocknum < MM->sparse->n ; blocknum++) {

            // number of blocks
            uint32_t rlen= MM->sparse->len[blocknum];
            ASM_COMMENT("critical loop");
            for(uint32_t i = i0 ; i < g ; i += nthr * b) {
                abzero(x, dst + aboffset(x, i), b);
            }
            for(uint32_t i = i0 ; i < rlen ; i+= nthr * b) {
                for(unsigned int j = 0 ; j < b ; j++) {
                    __builtin_prefetch(A + nthr * b * off1, 0);
                    __builtin_prefetch(src + aboffset(x, A[nthr * b * off2]), 0);
                    uint32_t z = (i + j) & (g-1);
                    uint32_t k = *A++;
                    /*
                    ASSERT((dst - dst0) < aboffset(x, ra[!d]-z));
                    ASSERT((src - src0) < aboffset(x, 1 + MM->dim[d]-k));
                    */
                    abadd(x, dst + aboffset(x, z), src + aboffset(x, k));
                }
                A += (nthr-1) * b;
            }
            ASM_COMMENT("end of critical loop");

            dst += aboffset(x, g);
        }
#if 1
#undef  b
#undef  g
#undef  nthr
#undef  off1
#undef  off2
#endif
    } else {
        /* Then it's not so fun. In fact it's even terribly slow, and
         * bug-prone (concurrent writes all over the place). Therefore we
         * prefer to stay single-threaded here.
         *
         * Note that there _could_ be a way around the problem, using
         * atomic xors for abadd. But it isn't worth the trouble.
         */
        if (throff != 0)
            return;

        if (MM->iteration == 10) {
            fprintf(stderr, "Warning: Doing many iterations with bad code\n");
        }

        src += aboffset(x, MM->dense->n * MM_DGRP_SIZE);

        if (MM->dense->n == 0) {
            abzero(x, dst, MM->dim[!d]);
        } else {
            /* We must test whether all threads are done with the dense
             * blocks */
            if (MM->dense->processed != MM->nthreads) {
                pthread_mutex_lock(&MM->mu);
                for( ; MM->dense->processed != MM->nthreads ; ) {
                    pthread_cond_wait(&MM->dense->recheck_please, &MM->mu);
                }
                pthread_mutex_unlock(&MM->mu);
            }
        }
        for(unsigned int blocknum = 0 ; blocknum < MM->sparse->n ; blocknum++) {
            ASM_COMMENT("critical loop (transposed mult, SLOW)");
            uint32_t rlen= MM->sparse->len[blocknum];
            for(uint32_t i = i0 ; i < rlen ; i+= b) {
                for(unsigned int j = 0 ; j < b ; j++) {
                    __builtin_prefetch(A + b * off1, 0);
                    __builtin_prefetch(dst + aboffset(x, A[b * off2]), 1);
                    uint32_t z = (i + j) % g;
                    uint32_t k = *A++;
                    /*
                    ASSERT((dst - dst0) < aboffset(x, 1 + MM->dim[!d]-k));
                    ASSERT((src - src0) < aboffset(x, ra[d]-z));
                    */
                    abadd(x, dst + aboffset(x, k), src + aboffset(x, z));
                }
            }
            ASM_COMMENT("end of critical loop (transposed mult, SLOW)");
            src += aboffset(x, g);
        }
    }
    ASM_COMMENT("end of sparse multiplication code");
}

void matmul_threaded_mul(matmul_ptr mm, abt * dst, abt const * src, int d)
{
    MM->working = 0;
    MM->done = 0;
    MM->dense->processed = 0;
    MM->dst = dst;
    MM->src = src;
    MM->d = d;

    /* This mutex is associated to the cond_in and cond_out conditions. */
    pthread_mutex_lock(&MM->mu);

    pthread_cond_broadcast(&MM->cond_in);
    /* All threads are now working */

    /* Now switch to the outgoing condition, wait for the last thread to
     * release us */
    pthread_cond_wait(&MM->cond_out,&MM->mu);
    pthread_mutex_unlock(&MM->mu);

    MM->iteration++;
}

void matmul_threaded_report(matmul_ptr mm MAYBE_UNUSED) {
}

void matmul_threaded_auxv(matmul_ptr mm, int op, va_list ap)
{
    if (op == MATMUL_AUX_GET_READAHEAD) {
        unsigned int * res = va_arg(ap, unsigned int *);
        // First, we allow at least one extra.
        // But more subtly, we need to make sure that irrespective of the
        // number of dense and sparse blocks, the last sparse block will
        // have plenty of space to end. For example if nrows=40 and there
        // is one dense block of 32 rows, with sgrp_size=64, we need to
        // reach 96 entries, not merely 64. So if s and d are the sparse
        // and dense group sizes, we need to reach at least (x&~(d-1)+s)
        unsigned int start = *res;
        *res = (start & ~(MM_DGRP_SIZE-1)) + MM->sgrp_size;
        if (*res > start) {
            return;
        } else {
            *res = next_multiple_of_powerof2(1+*res, MM->sgrp_size);
        }
    }
}

void matmul_threaded_aux(matmul_ptr mm, int op, ...)
{
    va_list ap;
    va_start(ap, op);
    matmul_threaded_auxv (mm, op, ap);
    va_end(ap);
}
/* vim: set sw=4: */
