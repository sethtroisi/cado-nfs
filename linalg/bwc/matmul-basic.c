#define _GNU_SOURCE     /* asprintf */
#define _DARWIN_C_SOURCE    /* for asprintf. _ANSI_SOURCE must be undefined */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>

#include "matmul.h"
#include "matmul-basic.h"
#include "readmat.h"
#include "abase.h"

/* This extension is used to distinguish between several possible
 * implementations of the product */
#define MM_EXTENSION   "-basic"
#define MM_MAGIC_FAMILY        0xa001UL
#define MM_MAGIC_VERSION       0x1000UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)


struct matmul_basic_data_s {
    /* repeat the fields from the public interface */
    struct matmul_public_s public_[1];
    /* now our private fields */
    size_t datasize;
    abobj_t xab;
    uint32_t * q;
};

void matmul_basic_clear(struct matmul_basic_data_s * mm)
{
    free(mm->q);
    free(mm);
}

static struct matmul_basic_data_s * matmul_basic_init()
{
    struct matmul_basic_data_s * mm;
    mm = malloc(sizeof(struct matmul_basic_data_s));
    return mm;
}

typedef int (*sortfunc_t) (const void *, const void *);

static int uint_cmp(unsigned int * a, unsigned int * b)
{
    if (*a < *b) return -1;
    else if (*b < *a) return 1;
    return 0;
}

struct matmul_basic_data_s * matmul_basic_build(abobj_ptr xx MAYBE_UNUSED, const char * filename, param_list pl MAYBE_UNUSED)
{
    struct matmul_basic_data_s * mm = matmul_basic_init();

    abobj_init_set(mm->xab, xx);

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
            mm->q = realloc(mm->q, alloc * sizeof(uint32_t));		\
        }								\
    } while (0)

    mm->public_->dim[0] = smat->nrows;
    mm->public_->dim[1] = smat->ncols;
    mm->public_->ncoeffs = 0;
    mm->q = NULL;

    for(unsigned int i = 0 ; i < mm->public_->dim[0] ; i++) {
        /* Our dst pointer for read_matrix_row never changes, since
         * we write matrix rows over and over again in memory. */
        read_matrix_row(f,smat,smat->data,1);
        CHECK_REALLOC(1 + smat->data[0]);
        uint32_t last = 0;
        mm->q[size++] = smat->data[0];

        qsort(smat->data + 1, smat->data[0],
                sizeof(smat->data[0]), (sortfunc_t) uint_cmp);

        for(unsigned int j = 0 ; j < smat->data[0] ; j++) {
            mm->q[size++] = smat->data[1+j] - last;
            last = smat->data[1+j];
        }
        mm->public_->ncoeffs += smat->data[0];
    }

    mm->datasize = size * sizeof(uint32_t);

    fclose(f);

    return mm;
}

struct matmul_basic_data_s * matmul_basic_reload_cache(abobj_ptr xx MAYBE_UNUSED, const char * filename, param_list pl MAYBE_UNUSED)
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
    struct matmul_basic_data_s * mm = matmul_basic_init();

    abobj_init_set(mm->xab, xx);

    unsigned long magic_check;
    rc = fread(&magic_check, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    FATAL_ERROR_CHECK(magic_check != MM_MAGIC,
            "Wrong magic in cached matrix file");
    rc = fread(&mm->public_->dim, sizeof(unsigned int), 2, f);
    FATAL_ERROR_CHECK(rc < 2, "No valid data in cached matrix file");
    rc = fread(&mm->public_->ncoeffs, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    rc = fread(&mm->datasize, sizeof(size_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "No valid data in cached matrix file");
    mm->q = malloc(mm->datasize);
    rc = fread(mm->q, 1, mm->datasize, f);
    FATAL_ERROR_CHECK(rc < mm->datasize, "Short read from cached matrix file");
    fclose(f);

    return mm;
}

void matmul_basic_save_cache(struct matmul_basic_data_s * mm, const char * filename)
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
    rc = fwrite(&mm->public_->dim, sizeof(unsigned int), 2, f);
    FATAL_ERROR_CHECK(rc < 2, "Cannot write to cached matrix file");
    rc = fwrite(&mm->public_->ncoeffs, sizeof(unsigned long), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(&mm->datasize, sizeof(size_t), 1, f);
    FATAL_ERROR_CHECK(rc < 1, "Cannot write to cached matrix file");
    rc = fwrite(mm->q, 1, mm->datasize, f);
    FATAL_ERROR_CHECK(rc < mm->datasize, "Short write from cached matrix file");
    fclose(f);
}

void matmul_basic_mul(struct matmul_basic_data_s * mm, abt * dst, abt const * src, int d)
{
    __asm__("# multiplication code\n");
    uint32_t * q = mm->q;
    abobj_ptr x = mm->xab;

    if (d == 1) {
        abzero(x, dst, mm->public_->dim[0]);
        __asm__("# critical loop\n");
        for(unsigned int i = 0 ; i < mm->public_->dim[0] ; i++) {
            uint32_t len = *q++;
            unsigned int j = 0;
            for( ; len-- ; ) {
                j += *q++;
                ASSERT(j < mm->public_->dim[1]);
                abadd(x, dst + aboffset(x, i), src + aboffset(x, j));
            }
        }
        __asm__("# end of critical loop\n");
    } else {
        abzero(x, dst, mm->public_->dim[1]);
        __asm__("# critical loop (transposed mult)\n");
        for(unsigned int i = 0 ; i < mm->public_->dim[0] ; i++) {
            uint32_t len = *q++;
            unsigned int j = 0;
            for( ; len-- ; ) {
                j += *q++;
                ASSERT(j < mm->public_->dim[1]);
                abadd(x, dst + aboffset(x, j), src + aboffset(x, i));
            }
        }
        __asm__("# end of critical loop (transposed mult)\n");
    }
    __asm__("# end of multiplication code\n");
}

void matmul_basic_report(struct matmul_basic_data_s * mm MAYBE_UNUSED) {
}

void matmul_basic_auxv(struct matmul_basic_data_s * mm MAYBE_UNUSED, int op MAYBE_UNUSED, va_list ap MAYBE_UNUSED)
{
}

void matmul_basic_aux(struct matmul_basic_data_s * mm, int op, ...)
{
    va_list ap;
    va_start(ap, op);
    matmul_basic_auxv (mm, op, ap);
    va_end(ap);
}

/* vim: set sw=4: */
