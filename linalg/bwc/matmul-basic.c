#define _GNU_SOURCE     /* asprintf */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>

#include "matmul.h"
#include "readmat.h"
#include "abase.h"
#include "manu.h"

/* This extension is used to distinguish between several possible
 * implementations of the product */
#define MM_EXTENSION   "-basic"
#define MM_MAGIC_FAMILY        0xa001UL
#define MM_MAGIC_VERSION       0x1000UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)


struct matmul_data_s {
    /* repeat the fields from the public interface */
    unsigned int nrows;
    unsigned int ncols;
    unsigned long ncoeffs;
    size_t datasize;
    abobj_t xab;
    uint32_t * q;
};


#define MM      ((struct matmul_data_s *) mm)

void matmul_basic_clear(matmul_ptr mm)
{
    free(MM->q);
    free(MM);
}

matmul_ptr matmul_basic_init()
{
    return malloc(sizeof(struct matmul_data_s));
}

typedef int (*sortfunc_t) (const void *, const void *);

static int uint_cmp(unsigned int * a, unsigned int * b)
{
    if (*a < *b) return -1;
    else if (*b < *a) return 1;
    return 0;
}

matmul_ptr matmul_basic_build(abobj_ptr xx MAYBE_UNUSED, const char * filename)
{
    matmul_ptr mm = matmul_basic_init();

    abobj_init_set(MM->xab, xx);

    sparse_mat_t smat;
    FILE * f;

    sparse_mat_init(smat);

    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", filename, strerror(errno));
        exit(1);
    }
    read_matrix_header(f, smat);

    size_t alloc = 0;
    size_t size = 0;

#define CHECK_REALLOC(n)      do {					\
        if (size + n >= alloc) {					\
            alloc = size + n;                                           \
            alloc += alloc / 2;                                         \
            MM->q = realloc(MM->q, alloc * sizeof(uint32_t));		\
        }								\
    } while (0)

    MM->nrows = smat->nrows;
    MM->ncols = smat->ncols;
    MM->ncoeffs = 0;
    MM->q = NULL;

    for(unsigned int i = 0 ; i < MM->nrows ; i++) {
        /* Our dst pointer for read_matrix_row never changes, since
         * we write matrix rows over and over again in memory. */
        read_matrix_row(f,smat,smat->data,1);
        CHECK_REALLOC(1 + smat->data[0]);
        uint32_t last = 0;
        MM->q[size++] = smat->data[0];

        qsort(smat->data + 1, smat->data[0],
                sizeof(smat->data[0]), (sortfunc_t) uint_cmp);

        for(unsigned int j = 0 ; j < smat->data[0] ; j++) {
            MM->q[size++] = smat->data[1+j] - last;
            last = smat->data[1+j];
        }
        MM->ncoeffs += smat->data[0];
    }

    MM->datasize = size * sizeof(uint32_t);

    return mm;
}

matmul_ptr matmul_basic_reload_cache(abobj_ptr xx MAYBE_UNUSED, const char * filename)
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
    matmul_ptr mm = matmul_basic_init();

    abobj_init_set(MM->xab, xx);

    unsigned long magic_check;
    rc = fread(&magic_check, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    FATAL_ERROR_CHECK(magic_check != MM_MAGIC,
            "Wrong magic in cached matrix file");
    rc = fread(&MM->nrows, sizeof(unsigned int), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->ncols, sizeof(unsigned int), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->ncoeffs, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&MM->datasize, sizeof(size_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    MM->q = malloc(MM->datasize);
    rc = fread(MM->q, 1, MM->datasize, f);
    FATAL_ERROR_CHECK(rc < MM->datasize, "Short read from cached matrix file");
    fclose(f);

    return mm;
}

void matmul_basic_save_cache(matmul_ptr mm, const char * filename)
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
    rc = fwrite(&MM->nrows, sizeof(unsigned int), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->ncols, sizeof(unsigned int), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->ncoeffs, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&MM->datasize, sizeof(size_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(MM->q, 1, MM->datasize, f);
    FATAL_ERROR_CHECK(rc < MM->datasize, "Short write from cached matrix file");
    fclose(f);
}

void matmul_basic_mul(matmul_ptr mm, abt * dst, abt const * src, int d)
{
    __asm__("# multiplication code\n");
    uint32_t * q = MM->q;
    abobj_ptr x = MM->xab;

    if (d == 1) {
        abzero(x, dst, MM->nrows);
        __asm__("# critical loop\n");
        for(unsigned int i = 0 ; i < MM->nrows ; i++) {
            uint32_t len = *q++;
            unsigned int j = 0;
            for( ; len-- ; ) {
                j += *q++;
                ASSERT(j < mm->ncols);
                abadd(x, dst + aboffset(x, i), src + aboffset(x, j));
            }
        }
        __asm__("# end of critical loop\n");
    } else {
        abzero(x, dst, MM->ncols);
        __asm__("# critical loop (transposed mult)\n");
        for(unsigned int i = 0 ; i < MM->nrows ; i++) {
            uint32_t len = *q++;
            unsigned int j = 0;
            for( ; len-- ; ) {
                j += *q++;
                ASSERT(j < mm->ncols);
                abadd(x, dst + aboffset(x, j), src + aboffset(x, i));
            }
        }
        __asm__("# end of critical loop (transposed mult)\n");
    }
    __asm__("# end of multiplication code\n");
}

void matmul_basic_report(matmul_ptr mm MAYBE_UNUSED) {
}

/* vim: set sw=4: */
