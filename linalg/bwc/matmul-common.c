#define _GNU_SOURCE         /* asprintf */
#define _DARWIN_C_SOURCE    /* for asprintf. _ANSI_SOURCE must be undefined */

#include <errno.h>
#include <string.h>
#include "matmul-common.h"
#include "params.h"
#include "readmat-easy.h"

#define MM_COMMON_MAGIC 0xb0010000UL

const char * rowcol[2] = { "row", "col", };

/* Factor out some stuff which turns out to appear fairly often */

static FILE * fopen_cache(struct matmul_public_s * mm, const char * filename, const char * ext, const char * mode)
{
    char * base;
    FILE * f;
    int rc;
   
    rc = asprintf(&base, "%s%s%s.bin", filename, ext, mm->store_transposed ? "T" : "");
    FATAL_ERROR_CHECK(rc < 0, "out of memory");

    f = fopen(base, mode);

    if (f == NULL) {
        // fprintf(stderr, "fopen(%s): %s\n", base, strerror(errno));
        fprintf(stderr, "no cache file %s\n", base);
        free(base);
        return NULL;
    }
    free(base);
    return f;
}

FILE * matmul_common_reload_cache_fopen(size_t stride, struct matmul_public_s * mm, const char * filename, const char * ext, uint32_t magic)
{
    FILE * f = fopen_cache(mm, filename, ext, "r");
    if (f == NULL) return NULL;

    uint32_t magic_check;
    MATMUL_COMMON_READ_ONE32(magic_check, f);
    
    if (magic_check != magic) {
        fprintf(stderr, "Wrong magic in cached matrix file\n");
        return NULL;
    }   
    
    MATMUL_COMMON_READ_ONE32(magic_check, f);
    if (magic_check != MM_COMMON_MAGIC) {
        fprintf(stderr, "Wrong magic in cached matrix file\n");
        return NULL;
    }   

    uint32_t nbytes_check;
    MATMUL_COMMON_READ_ONE32(nbytes_check, f);
    /* It's not fatal. It only deserves a warning */
    if (nbytes_check != stride) {
        fprintf(stderr, "Warning: cached matrix file fits data with different striding\n");
    }

    MATMUL_COMMON_READ_ONE32(mm->dim[0], f);
    MATMUL_COMMON_READ_ONE32(mm->dim[1], f);

    MATMUL_COMMON_READ_ONE64(mm->ncoeffs, f);

    return f;
}

FILE * matmul_common_save_cache_fopen(size_t stride, struct matmul_public_s * mm, const char * filename, const char * ext, uint32_t magic)
{
    FILE * f = fopen_cache(mm, filename, ext, "w");
    if (f == NULL) return NULL;

    MATMUL_COMMON_WRITE_ONE32(magic,f);
    MATMUL_COMMON_WRITE_ONE32(MM_COMMON_MAGIC,f);
    MATMUL_COMMON_WRITE_ONE32(stride,f);
    MATMUL_COMMON_WRITE_ONE32(mm->dim[0],f);
    MATMUL_COMMON_WRITE_ONE32(mm->dim[1],f);
    MATMUL_COMMON_WRITE_ONE64(mm->ncoeffs,f);

    return f;
}

void matmul_common_init_post(struct matmul_public_s * mm, param_list pl, int suggest)
{
    mm->store_transposed = suggest;
    
    if (pl) {
        param_list_parse_uint(pl, "bwc_mm_store_transposed", &mm->store_transposed);
        if (mm->store_transposed != (unsigned int) suggest) {
            fprintf(stderr, "Warning, bwc_mm_store_transposed"
                    " overrides suggested matrix storage ordering\n");
        }   
    }
}

/* Okay, this matrix reading stage is really a memory hog. But it's
 * done only once, so we practically don't care. */
uint32_t * matmul_common_read_stupid_data(struct matmul_public_s * mm, const char * filename)
{
    /* dims of the possibly transposed matrix */
    unsigned int nrows_t;
    unsigned int ncols_t;
    uint32_t * data;
    if (mm->store_transposed) {
        data = read_easy_transposed(filename, &nrows_t, &ncols_t);
    } else {
        data = read_easy(filename, &nrows_t, &ncols_t);
    }
    mm->dim[ mm->store_transposed] = nrows_t;
    mm->dim[!mm->store_transposed] = ncols_t;
    return data;
}
