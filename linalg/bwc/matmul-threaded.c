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
#define MM_MAGIC_VERSION       0x1000UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)


#define RLEN_ALIGN      1024

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
    size_t datasize;

    /* Now our private fields */
    abobj_t xab;
    uint32_t * q;

    uint32_t numrows_group;
    uint32_t nthreads;
    uint32_t blocksize;
    uint32_t transposed;
    uint32_t off1;
    uint32_t off2;

    pthread_mutex_t mu;
    pthread_cond_t cond_in;
    pthread_cond_t cond_out;

    /* Number of worker threads at work (not idle) */
    unsigned int working;

    struct thread_info * thread_table;

    pthread_t * threads;
    
    void * dst;
    const void * src;
    int d;
    int iteration;
};

#define MM      ((struct matmul_threaded_data_s *) mm)

extern void matmul_threaded_mul_sub(matmul_ptr mm, abt * dst, abt const * src, int d, int i);

void * thread_job(void * ti)
{
    matmul_ptr mm = ((struct thread_info *) ti) -> mm;
    int i = ((struct thread_info *) ti) -> i;

    pthread_mutex_lock(&MM->mu);

    MM->working++;
    pthread_cond_signal(&MM->cond_out);
    /* The mutex will soon be released */

    unsigned int res = 0;

    for( ; ; ) {
        pthread_cond_wait(&MM->cond_in,&MM->mu);
        MM->working++;
        pthread_mutex_unlock(&MM->mu);
        if (MM->iteration < 0)
            break;

        /* Do our job with mutexes released */
        matmul_threaded_mul_sub(mm, MM->dst, MM->src, MM->d, i);

        pthread_mutex_lock(&MM->mu);
        res = ++(MM->working);

        if (res == 2 * MM->nthreads) {
            pthread_cond_signal(&MM->cond_out);
        }
    }

    pthread_mutex_unlock(&MM->mu);

    printf("Thread %d exits\n", i);

    return NULL;
}


void worker_threads_init(matmul_ptr mm)
{
    MM->threads = malloc(MM->nthreads * sizeof(pthread_t));
    MM->thread_table = malloc(MM->nthreads * sizeof(struct thread_info));
    pthread_cond_init(&MM->cond_in, NULL);
    pthread_cond_init(&MM->cond_out, NULL);
    pthread_mutex_init(&MM->mu, NULL);
    MM->working = 0;
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

    free(MM->q);
    free(MM);
}

static matmul_ptr matmul_threaded_init()
{
    matmul_ptr mm = malloc(sizeof(struct matmul_threaded_data_s));

    /* TODO: Turn these into parameters */
    /* TODO: remove it from the on-disk file, since it does not impact
     * the data itself. */
    MM->nthreads = 4;

    /* TODO: Turn these into parameters */
    /* Transposing here is just an optimisation concern. The code does
     * not compute anything different depending on this value. But it
     * performs differently */
    MM->transposed = 1;
    MM->numrows_group = 32;

    /* Must correspond to a L1 cache line. */
    /* XXX Must also correspond to the selected abase ! */
    MM->blocksize = 8;

    MM->off1 = 32;
    MM->off2 = 16;

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

#define CHECK_REALLOC(ptr__,alloc__,size__,n__) do {    		\
        if (size__ + n__ >= alloc__) {					\
            alloc__ = size__ + n__;                                     \
            alloc__ += alloc__ / 2;                                     \
            ptr__ = realloc(ptr__, alloc__ * sizeof(uint32_t));		\
        }								\
    } while (0)


matmul_ptr matmul_threaded_build(abobj_ptr xx MAYBE_UNUSED, const char * filename)
{
    matmul_ptr mm = matmul_threaded_init();

    ASSERT_ALWAYS(MM->numrows_group % MM->blocksize == 0);

    abobj_init_set(MM->xab, xx);

    MM->ncoeffs = 0;
    MM->q = NULL;
    size_t alloc = 0;
    size_t size = 0;

    /* Okay, this matrix reading stage is really a memory hog. But it's
     * done only once, and we can even consider doing it only once, so we
     * practically don't care. */

    uint32_t * data;

    /* dimensions of the transposed matrix */
    unsigned int nrows_t;
    unsigned int ncols_t;

    if (MM->transposed) {
        data = read_easy_transposed(filename, &nrows_t, &ncols_t);
    } else {
        data = read_easy(filename, &nrows_t, &ncols_t);
    }

    MM->dim[ MM->transposed] = nrows_t;
    MM->dim[!MM->transposed] = ncols_t;

    /* We give some notations and example values.  */
    uint32_t n = MM->nthreads;      /* e.g. 4 */
    uint32_t g = MM->numrows_group; /* e.g. 32 */
    uint32_t b = MM->blocksize;     /* e.g. 16 */

    uint32_t number_of_big_blocks;
    number_of_big_blocks = iceildiv(nrows_t, g);
    printf("// %u blocks of %d rows\n", number_of_big_blocks, g);

    /* Make sure the area we allocate for lengths is properly aligned */
    unsigned int rlen_area = ((number_of_big_blocks-1) | (RLEN_ALIGN-1))+1;

    CHECK_REALLOC(MM->q,alloc,size,rlen_area);

    memset(MM->q, 0, rlen_area * sizeof(uint32_t));

    size= rlen_area;

    unsigned int bigblock_number = 0;

    unsigned int padding_coeffs = 0;

    uint32_t * ptr = data;

    for(unsigned int i = 0 ; i < nrows_t ; ) {
        unsigned int di;

        /* common length of the rows */
        uint32_t l = *ptr;

        /* striding for the entries in one row is g==numrows_group */
         
        CHECK_REALLOC(MM->q,alloc,size,l*g);
        
        /*
         * We have (ff == and following, for b indices).
         *
         * blk 0x00: idx 0 for row   0      ff --> thr 0
         * blk 0x01: idx 0 for row   b      ff --> thr 1
         * blk 0x02: idx 0 for row  2b      ff --> thr 2
         * blk 0x03: idx 0 for row  3b      ff --> thr 3
         * blk 0x04: idx 0 for row  4b      ff --> thr 0
         * 
         * IOW, it's easy. entry g*l+k : corresponds to index l in
         * row k, and we recover the striding g for row entries. A
         * complete block of g indices is covered by the n threads in
         * g/(nb) blocks of b rows.
         */
        ASSERT_ALWAYS(g % (n*b) == 0);

        /* We store the total number of entries in
         * MM->q[bigblock_number]. A given thread increments its counters
         * by steps of n*b, because it effectively uses only b coeffs
         * every n*b.
         */
        ASSERT_ALWAYS(bigblock_number < number_of_big_blocks);
        MM->q[bigblock_number++] = l * g;

        for(di = 0 ; di < g ; di++) {
            uint32_t rlen;
            
            if (i + di < nrows_t) {
                rlen = *ptr++;
                qsort(ptr, rlen, sizeof(uint32_t), (sortfunc_t) uint_cmp);
            } else {
                rlen = 0;
            }

            uint32_t * myrow = MM->q + size + di;

            unsigned int j;
            for(j = 0 ; j < rlen ; j++) {
                myrow[j * g] = *ptr++;
            }
            MM->ncoeffs += rlen;

            /* padding, trash data */
            for( ; j < l ; j++) {
                myrow[j * g] = ncols_t;
                padding_coeffs++;
            }
        }

        i += g;
        size += l * g;
    }

    /* fill up a small readahead zone */
    CHECK_REALLOC(MM->q,alloc,size,g+MM->off2*b);
    for(unsigned int i = g ; i ; i--) {
        MM->q[size++]=0;
    }

    free(data);

    printf("// %u padding coeffs, %lu real coeffs\n",
            padding_coeffs, MM->ncoeffs);

    MM->datasize = size * sizeof(uint32_t);

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
    rc = fread(&MM->datasize, sizeof(size_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");

    rc = fread(&MM->numrows_group, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->nthreads, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->transposed, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");

    MM->q = malloc(MM->datasize);
    rc = fread(MM->q, 1, MM->datasize, f);
    FATAL_ERROR_CHECK(rc < MM->datasize, "Short read from cached matrix file");
    fclose(f);

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
    rc = fwrite(&MM->datasize, sizeof(size_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");

    rc = fwrite(&MM->numrows_group, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->nthreads, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->transposed, sizeof(uint32_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");

    rc = fwrite(MM->q, 1, MM->datasize, f);
    FATAL_ERROR_CHECK(rc < MM->datasize, "Short write from cached matrix file");
    fclose(f);
}

void matmul_threaded_mul_sub(matmul_ptr mm, abt * dst, abt const * src, int d, int throff)
{
    __asm__("# multiplication code\n");
    uint32_t * A = MM->q;
    abobj_ptr x = MM->xab;

    uint32_t nthr = MM->nthreads;      /* e.g. 4 */
    uint32_t g = MM->numrows_group; /* e.g. 32 */
    uint32_t b = MM->blocksize;     /* e.g. 16 */

    uint32_t off1 = MM->off1;
    uint32_t off2 = MM->off2;

    uint32_t number_of_big_blocks;
    number_of_big_blocks = iceildiv(MM->dim[MM->transposed], g);

    unsigned int rlen_area = ((number_of_big_blocks-1) | (RLEN_ALIGN-1))+1;

    A += rlen_area;
    A += throff * b;

    /* If d == 0, output has ncols == dim[1] entries
     *             input has nrows == dim[0] entries
     */

    abt * dst0 = dst;
    abt const * src0 = src;

    unsigned int ra[2];
    ra[0] = MM->dim[0];
    matmul_threaded_aux(mm, MATMUL_AUX_GET_READAHEAD, &(ra[0]));
    ra[1] = MM->dim[1];
    matmul_threaded_aux(mm, MATMUL_AUX_GET_READAHEAD, &(ra[1]));

    if (d == !MM->transposed) {
        /* t == 1 is well suited for d == 0. In this case, dst is
         * progressively filled. The same holds for t == 0 and d == 1 */
#define b       8
#define g       32
#define nthr    4
#define off1    32
#define off2    16
        for(unsigned int bnum = 0 ; bnum < number_of_big_blocks ; bnum++) {

            // number of blocks
            uint32_t rlen= MM->q[bnum];

            uint32_t i0 = throff * b;

            __asm__("# critical loop\n");
            for(uint32_t i = i0 ; i < g ; i += nthr * b) {
                abzero(x, dst + aboffset(x, i), b);
            }
            for(uint32_t i = i0 ; i < rlen ; i+= nthr * b) {
                for(unsigned int j = 0 ; j < b ; j++) {
                    __builtin_prefetch(A + nthr * b * off1, 0);
                    __builtin_prefetch(src + aboffset(x, A[nthr * b * off2]), 0);
                    uint32_t z = (i + j) & (g-1);
                    uint32_t k = *A++;
                    ASSERT((dst - dst0) < aboffset(x, ra[!d]-z));
                    ASSERT((src - src0) < aboffset(x, 1 + MM->dim[d]-k));
                    abadd(x, dst + aboffset(x, z), src + aboffset(x, k));
                }
                A += (nthr-1) * b;
            }
            __asm__("# end of critical loop\n");

            dst += aboffset(x, g);
        }
#undef  b
#undef  g
#undef  nthr
#undef  off1
#undef  off2
    } else {
        /* Then it's not so fun. In fact it's even terribly slow, and
         * bug-prone (concurrent writes all over the place). Therefore we
         * prefer to stay single-threaded here.
         */
        if (throff != 0)
            return;
        if (MM->iteration == 10) {
            fprintf(stderr, "Warning: Doing many iterations with bad code\n");
        }

        abzero(x, dst, MM->dim[!d]);
        for(unsigned int bnum = 0 ; bnum < number_of_big_blocks ; bnum++) {
            __asm__("# critical loop (transposed mult, SLOW)\n");
            uint32_t rlen= MM->q[bnum];
            uint32_t i = throff * b;
            for(i = 0 ; i < rlen ; i+= b) {
                for(unsigned int j = 0 ; j < b ; j++) {
                    __builtin_prefetch(A + b * off1, 0);
                    __builtin_prefetch(dst + aboffset(x, A[b * off2]), 1);
                    uint32_t z = (i + j) % g;
                    uint32_t k = *A++;
                    ASSERT((dst - dst0) < aboffset(x, 1 + MM->dim[!d]-k));
                    ASSERT((src - src0) < aboffset(x, ra[d]-z));
                    abadd(x, dst + aboffset(x, k), src + aboffset(x, z));
                }
            }
            __asm__("# end of critical loop (transposed mult, SLOW)\n");
            src += aboffset(x, g);
        }
    }
    __asm__("# end of multiplication code\n");
}

void matmul_threaded_mul(matmul_ptr mm, abt * dst, abt const * src, int d)
{
    MM->working = 0;
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
        // allow at least one extra.
        *res += 1;
        // and pad to a multiple of numrows_group
        *res = ((*res-1) | (MM->numrows_group - 1)) + 1;
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
