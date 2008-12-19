/* Cutting a matrix file into small pieces for parallel linear algebra

Copyright 2008 Emmanuel Thom\'e

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.


The purpose of this program is to cut a big matrix file into several small
files, corresponding to the subsequent distribution of a linear algebra step
on several computers.

Features:
  - The output splitting is a grid n x m
  - The rows (resp. cols) are weight sorted during the process
  - If the sorting is done according to the weight of rows, then
    the same permutation is used for columns
  - The temporary files, the input files and so on, can be optionnally
    kept or removed asap
  - Memory usage is under control (tunable)
  - Can pad the matrix with zeros to get a square matrix

Disk usage:
  If the input matrix has size 100 (whatever the unit is), then the total
  size of the output will be 100. The top disk usage (if the input is kept)
  will be 100 + 100 + 100/n, where
    - the first 100 are for the input
    - the second 100 are for the ouput
    - the 100/n are for temporaries
    - n is the number of horizontal slices
  If the ram-limit is high enough to contain a whole horizontal slice, then
  the 100/n term is not present.

  In the case where the input matrix can be erased asap, then the top disk
  usage is 100 + 100.
  TODO: Would it be possible to have it in 100 + 100/n ???
        This is not at all clear that POSIX allows to do that with files.

List of options: (-in and -out are mandatory)
  -in <path>    input matrix filename
  -out <path>   output matrix filename 
                If there is a directory path prepending the name, then
                the same subdirectory is used for temp files.
  -nslices n    balance for n horizontal slices and 1 vertical slice 
           nxm  balance for n horizontal slices and m vertical slices
  -square       Pad the matrix with zeroes to get a square matrix
  -permute-rows-like-columns
                The weight-sorting is done on columns and the permutation
                is applied also to rows
  -permute-columns-like-rows
                The weight-sorting is done on rows and the permutation
                is applied also to columns

Other options that do not change the ouput:
  -remove-input Erase the input matrix file as soon as possible
  -ram-limit <x>
                Do not use more than x memory. x is of the form 14M, or 2G
                or 8000k ...
                The default is 256 MB.
  -keep-temps   Do not erase temp files (for debugging purpose)
  -legacy       Obsolete output format
  
*/

#define __STDC_FORMAT_MACROS    /* For C++ to have PRId64 and friends */
#define _GNU_SOURCE    /* for strndup and asprintf */

#include <stddef.h>     /* ptrdiff_t */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <inttypes.h>   /* SCNu32 etc */
#include <assert.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>     /* for sqrt -- we're computing the sdev of the weights */

/* platform headers */
#include <sys/stat.h>   /* mkdir */
#include <sys/types.h>  /* mkdir */
#include <unistd.h>     /* chdir */
#include <time.h>       /* ctime */

#include "readbuffer.h"
#include "rowset_heap.h"
#include "utils.h"

/*{{{ macros */
#ifdef  __cplusplus
using namespace std;
#endif

#define PARSE_CHECK(tst, kind, text) do {				\
    if ((tst)) {							\
        fprintf(stderr, "parse error while reading %s: %s",		\
                kind, text);						\
        exit(1);							\
    }									\
} while (0)
/*}}}*/

typedef int (*sortfunc_t)(const void *, const void *);

/*{{{ globals ; flags, mostly */

#define BIGCHUNK        16384

unsigned int nhslices = 1;
unsigned int nvslices = 1;

int keep_temps = 0;
int remove_input = 0;
int pad_to_square = 0;
int legacy = 0;

int permute_rows = 0;
int permute_cols = 0;

size_t ram_limit = 1 << 28;     /* 256 MB */

/* If we do this, then necessarily we'll have to copy the matrix,
 * unless... */
int weight_sort_in_cells = 0;

/* Upon reading the matrix, we'll most probably clear these flags */
/* If the flag stays up, then we won't have to copy the matrix */
/* Note also that if we slice the matrix, the way we populate slices
 * implies that sorting will be stay (we'll still have to do the copy though).
 */
int rows_are_weight_sorted = 1;
int cols_are_weight_sorted = 1;

int replicate_rows_perm_for_columns = 0;
int replicate_columns_perm_for_rows = 0;

/* full path to the input matrix */
char * pristine_filename;

/* basename of the resulting matrix filename */
char * working_filename;
/*}}}*/

void usage()
{
    fprintf(stderr, "Usage: ./bw-balance <options>\n");
    fprintf(stderr,
            "Typical options:\n"
            "--in <path>\tinput matrix filename\n"
            "--out <path>\toutput matrix filename\n"
            "--nslices <n1>[x<n2>]\toptimize for <n1>x<n2> strips\n"
            "--square\tpad matrix with zeroes to obtain square size\n"
        );
    fprintf(stderr,
            "More advanced:\n"
            "--remove-input\tremove the input file as soon as possible\n"
            "--ram-limit <nnn>[kmgKMG]\tfix maximum memory usage\n"
            "--keep-temps\tkeep all temporary files\n"
            "--subdir <d>\tchdir to <d> beforehand (mkdir if not found)\n"
            "--legacy\tproduce only one jumbo matrix file\n"
           );
    exit(1);
}

struct row {
    uint32_t    i;
    int         w;
};

uint32_t nr;
struct row * row_table /* = NULL */;

uint32_t nc;
struct row * col_table /* = NULL */;

uint32_t nmax;

/* These two really only have to do with optimizing I/O for the matrix.
 * In case we're able to skip parsing entirely, we save line sizes for
 * reuse.
 */
size_t * line_bytes;
size_t header_bytes;    // used for seeking in the input file
char * header /* = NULL */;

struct slice * row_slices;
struct slice * col_slices;

struct file_info_s {
    uint32_t coeffs;
    char * s;
};
typedef struct file_info_s file_info[1];
file_info * final_file_info;

/* {{{ Several data types */
/* {{{ Datatype containing description for slices of the matrix */
struct slice {
    uint32_t * r;
    uint32_t nrows;
    uint32_t coeffs;
    uint32_t i0;
};

int slice_finding_helper(const uint32_t * key, const struct slice * elem)
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
}


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

    uint32_t min_size = water / n;

    for(i = 0 ; i < n ; i++) {
        res[i].nrows = min_size + (i < (water % n));
        res[i].r = (uint32_t *) malloc(res[i].nrows * sizeof(uint32_t));
    }

    return res;
}
/* }}} */
/* {{{ copy buffers -- useful for copying large chunks of data when there
 * is no parsing involved at all
 */
struct copybuf_s {
    char * buf;
    size_t siz;
};

typedef struct copybuf_s copybuf[1];

void cb_init(copybuf cb) {
    memset(cb, 0, sizeof(copybuf));
}
void cb_clear(copybuf cb) {
    free(cb->buf);
}
void cb_ensure(copybuf cb, size_t amount__)
{
    if (cb->siz >= (amount__)) return;
    if (cb->siz == 0) cb->siz=1024;
    assert(amount__ < (1 << 30));
    for( ; cb->siz < (amount__) ; cb->siz <<= 1);
    cb->buf = (char *) realloc(cb->buf, cb->siz);
}
/* }}} */
/* {{{ A data type for managing sets of files */
enum fileset_status { INPUT, TEMP, OUTPUT };
struct fileset_s {
    char ** names;
    unsigned int n;
    enum fileset_status status;
};

typedef struct fileset_s fileset[1];
typedef struct fileset_s * fileset_ptr;

void fileset_init(fileset x, unsigned int n)
{
    x->names = (char **) malloc(n * sizeof(char *));
    memset(x->names, 0, n * sizeof(char *));
    x->n = n;
}

/* This possibly unlinks one of the files within an file set, but only
 * if the flags keep_temps, remove_input are in accordance with the
 * status of this particular input set.
 */
void fileset_dispose(fileset x, unsigned int i)
{
    int yes=0;
    switch(x->status) {
        case INPUT:
            assert(x->n == 1 && i == 0);
            yes = remove_input;
            break;
        case TEMP:
            yes = !keep_temps;
            break;
        case OUTPUT:
        default:
            break;
    }
    if (yes) {
        unlink(x->names[i]);
    }
}

void fileset_clear(fileset x)
{
    unsigned int i;
    for(i = 0 ; i < x->n ; i++) {
        free(x->names[i]);
    }
    free(x->names);
    memset(x, 0, sizeof(fileset));
}

/* Fixup things such as tmp-/some/dir/matrix.txt into
 * /some/dir/tmp-matrix.txt */
void prefix_fixup(char * s, const char * pfx)
{
    size_t s1,s2;
    char * last_slash;

    s1 = strlen(pfx);

    if (s1 == 0)
        return;
    assert(strncmp(s, pfx, s1) == 0);
    last_slash = strrchr(s, '/');
    if (last_slash == NULL)
        return;
    last_slash++;
    s2 = last_slash - s;

    memmove(s, s+s1, s2-s1);
    strncpy(s+s2-s1, pfx, s1);
}

/* POSIX dirname sucks */
char * my_dirname(const char * s)
{
    const char * last_slash = strrchr(s, '/');
    if (last_slash == NULL)
        return NULL;
    else
        return strndup(s, last_slash - s);
}

const char * my_basename(const char * s)
{
    const char * last_slash = strrchr(s, '/');
    if (last_slash == NULL)
        return s;
    else
        return last_slash+1;
}

#if 0
int basename_has_prefix(const char * s, const char * pfx)
{
    size_t s1;
    const char * last_slash;
    s1 = strlen(pfx);
    if (s1 == 0)
        return 1;
    last_slash = strrchr(s, '/');
    if (last_slash == NULL)
        last_slash = s;
    else
        last_slash++;
    return strncmp(last_slash, pfx, s1) == 0;
}
int has_suffix(const char * s, const char * x)
{
    size_t ls = strlen(s);
    size_t lx = strlen(x);
    size_t j;
    s += ls-1;
    x += lx-1;
    for(j = 0 ; j < lx && j < ls && *s-- == *x-- ; j++) ;
    return j == lx;
}
#endif


void remove_prefix(char * s, const char * pfx)
{
    size_t s1;
    char * last_slash;
    s1 = strlen(pfx);
    if (s1 == 0)
        return;
    last_slash = strrchr(s, '/');
    if (last_slash == NULL)
        last_slash = s;
    else
        last_slash++;
    ASSERT_ALWAYS(strncmp(last_slash, pfx, s1) == 0);
    memmove(last_slash, last_slash + s1, strlen(last_slash + s1) + 1);
}

void fileset_name(fileset x, const char * name, const char * key)
{
    unsigned int i;
    for(i = 0 ; i < x->n ; i++) {
        free(x->names[i]);
    }
    const char * prefix = x->status == TEMP ? "tmp-" : "";
    if (x->n == 1) {
        asprintf(&(x->names[0]), "%s%s", prefix, name);
        prefix_fixup(x->names[0], prefix);
    } else {
        for(i = 0 ; i < x->n ; i++) {
            asprintf(&(x->names[i]), "%s%s.%s%u", prefix, name, key, i);
            prefix_fixup(x->names[i], prefix);
        }
    }
}

void fileset_init_transfer(fileset y,
        unsigned int n, const char * key,
        fileset x,
        unsigned int i,
        enum fileset_status s)
{
    fileset_init(y,n);
    y->status = s;
    if (x->status == INPUT) {
        fileset_name(y, working_filename, key);
    } else {
        assert(x->status == TEMP);
        const char * v = x->names[i];
        char * vv = strdup(v);
        remove_prefix(vv, "tmp-");
        fileset_name(y, vv, key);
        free(vv);
    }
}

FILE ** fileset_open(fileset x, const char * mode)
{
    unsigned int i;
    FILE ** res = (FILE **) malloc(x->n * sizeof(FILE *));
    for(i = 0 ; i < x->n ; i++) {
        res[i] = fopen(x->names[i], mode);
        DIE_ERRNO_DIAG(res[i] == NULL, "fopen", x->names[i]);
    }
    return res;
}
void fileset_close(fileset x, FILE ** f)
{
    unsigned int i;
    for(i = 0 ; i < x->n ; i++) {
        fclose(f[i]);
    }
    free(f);
}
void fileset_swap(fileset x, fileset y)
{
    fileset z;
    memcpy(z,x,sizeof(fileset));
    memcpy(x,y,sizeof(fileset));
    memcpy(y,z,sizeof(fileset));
}
void fileset_init_isolate(fileset y, fileset x, unsigned int i)
{
    y->status = x->status;
    y->n = 1;
    y->names=(char **) malloc(sizeof(char*));
    y->names[0] = strdup(x->names[i]);
}
FILE * fileset_open_one(fileset x, unsigned int i, const char * mode)
{
    FILE * res = fopen(x->names[i], mode);
    DIE_ERRNO_DIAG(res == NULL, "fopen", x->names[i]);
    if (x->status == INPUT) {
        fseek(res, header_bytes, SEEK_SET);
    }
    return res;
}
void fileset_close_one(fileset x MAYBE_UNUSED, unsigned int i MAYBE_UNUSED, FILE * f)
{
    fclose(f);
}

/* }}} */
/* {{{ column permutation buffers */
struct cperm_buf_s {
    copybuf * bufs;
    unsigned int * pos;
    unsigned int * w;
    unsigned int n;
};
typedef struct cperm_buf_s cperm_buf[1];
/* This merely does a reset for later reuse. */
void cperm_reset(cperm_buf cp)
{
    unsigned int i;
    for(i = 0 ; i < cp->n ; i++) {
        cp->pos[i] = cp->w[i] = 0;
    }
}
void cperm_init(cperm_buf cp, unsigned int n)
{
    unsigned int i;
    cp->n=n;
    cp->w = (unsigned int *) malloc(n * sizeof(unsigned int));
    cp->pos = (unsigned int *) malloc(n * sizeof(unsigned int));
    cp->bufs = (copybuf *) malloc(n * sizeof(copybuf));
    for(i = 0 ; i < cp->n ; i++) {
        cb_init(cp->bufs[i]);
    }
    cperm_reset(cp);
}
void cperm_clear(cperm_buf cp)
{
    unsigned int i;
    for(i = 0 ; i < cp->n ; i++) {
        cb_clear(cp->bufs[i]);
    }
    free(cp->bufs);
    free(cp->w);
    free(cp->pos);
}
void cperm_append(cperm_buf cp, unsigned int i, const char * text, unsigned int sz)
{
    cb_ensure(cp->bufs[i], cp->pos[i] + sz);
    memcpy(cp->bufs[i]->buf + cp->pos[i], text, sz);
    cp->pos[i] += sz;
}
/* }}} */
/* }}} */


// {{{ give mean, std dev and so on
void give_stats(const char * text, const struct row * data, unsigned int n)
{
    uint32_t dmin = UINT_MAX;
    uint32_t dmax = 0;
    double sumcf = 0;
    double sumcf2 = 0;
    uint32_t k;
    unsigned int nz = 0;
#define STORE_ZEROES    10
    unsigned int zz[STORE_ZEROES];

    for(k = 0 ; k < n ; k++) {
        uint32_t cf = data[k].w;
        if (data[k].w == 0) {
            if (nz < STORE_ZEROES)
                zz[nz]=k;
            nz++;
        }
        double dcf = cf;
        if (cf < dmin) dmin = cf;
        if (cf > dmax) dmax = cf;
        sumcf  += dcf;
        sumcf2 += dcf * dcf;
    }

    double mean = sumcf / n;
    /* weighted mean is the average weight of the row/column of a
     * coefficient sits in */
    double wmean = sumcf2 / sumcf;
    double sdev = sqrt(sumcf2 / n - mean*mean);

    printf("%s weights: [%"PRIu32"..%"PRIu32
            "], mean=%.1f, wmean=%.1f, sdev=%.1f\n",
            text, dmin, dmax, mean, wmean, sdev);

    if (nz) {
        unsigned int i;
        printf("%u zero %s:", nz, text);
        for(i = 0 ; i < nz && i < STORE_ZEROES ; i++) {
            printf(" %u", zz[i]);
        }
        if (i < nz) {
            printf(" ...");
        }
        printf("\n");
    }
    /*
       for(k = 0 ; k < nc ; k++) {
       printf(" %"PRIu32, col_freq[k]);
       }
       printf(" \n");
       */
}
// }}}

/* {{{ read the matrix */
void read_matrix_trailer()
{
    unsigned int i;

    if (permute_rows) give_stats("row", row_table, nr);
    if (permute_cols) give_stats("column", col_table, nc);

    if (weight_sort_in_cells) {
        assert(permute_rows);
        for(i = 1 ; i < nr ; i++) {
            if (row_table[i].w > row_table[i-1].w) {
                rows_are_weight_sorted=0;
                break;
            }
        }
        if (rows_are_weight_sorted) {
            printf("rows are already sorted by weight\n");
            permute_rows = (nhslices > 1);
        }

        assert(permute_cols);
        for(i = 1 ; i < nc ; i++) {
            if (col_table[i].w > col_table[i-1].w) {
                cols_are_weight_sorted=0;
                break;
            }
        }
        if (cols_are_weight_sorted) {
            printf("cols are already sorted by weight\n");
            permute_cols = (nvslices > 1);
        }
    }
}

void update_header()
{
    /* nmax is the largest dimension. If we pad to square size, we will
     * work with a matrix of size nmax * nmax.
     */
    nmax = nc < nr ? nr : nc;

    fprintf(stderr, "Matrix has %"PRIu32" rows, %"PRIu32" columns\n",
            nr, nc);

    /* Make this only a warning now -- it's actually going to disappear,
     * as this program is meant to become completely generic. */
    if (nc < nr) {
        fprintf(stderr, "Matrix has more rows than columns\n"
                "Perhaps the matrix should have been transposed first\n");
    /*
        exit(1);
    */
    }

    if (pad_to_square) {
        fprintf(stderr, "Output matrix will be padded to"
                " %"PRIu32" rows,"
                " %"PRIu32" columns\n",
                nmax, nmax);
        /* update the header accordingly ; it's no longer correct to copy
         * the old header verbatim */
        header = realloc(header, header_bytes + 10);
        if (strncmp(header, "//", 2) == 0) {
            /* untested ! */
            int pos;
            sscanf(header, "// %*" SCNu32 " ROWS %*" SCNu32 "%n", &pos);
            char * trailer = strdup(header + pos);
            snprintf(header, header_bytes + 10,
                    "// %" PRIu32 " ROWS %" PRIu32 "%s", nmax, nmax, trailer);
            free(trailer);
        } else {
            snprintf(header, header_bytes + 10,
                    "%" PRIu32 " %" PRIu32 "\n", nmax, nmax);
        }
    }
}

int read_matrix()
{
    reading_buffer b;

    rb_open(b, pristine_filename);
    uint32_t i;
    int rc;

    /* Read matrix header */
    rb_read_line(b);
    if (strncmp(b->buf, "//", 2) == 0) {
        rc = sscanf(b->buf, "// %" SCNu32 " ROWS %" SCNu32, &nr, &nc);
        PARSE_CHECK(rc < 2, "header", b->buf);
    } else {
        rc = sscanf(b->buf, "%" SCNu32 " %" SCNu32, &nr, &nc);
        PARSE_CHECK(rc < 2, "header", b->buf);
    }
    rb_gobble_long_line(b);
    header_bytes = b->o;
    header = strndup(b->buf, b->o);

    update_header();

    if (permute_rows) {
        row_table = (struct row*) malloc(nmax * sizeof(struct row));
        if (row_table == NULL) abort();
    }
    if (permute_cols) {
        col_table = (struct row *) malloc(nmax * sizeof(struct row));
        if (col_table == NULL) abort();
    }

    if (permute_cols) {
        unsigned int k;
        memset(col_table, 0, nmax * sizeof(struct row));
        for(k = 0 ; k < nmax ; k++) {
            col_table[k].i=k;
        }
    }

    line_bytes = (size_t *) malloc(nmax * sizeof(size_t));
    if (line_bytes == NULL) abort();

    assert(permute_rows);

    fprintf(stderr, "Reading matrix (first pass, counting weights)\n");
    for(i = 0 ; i < nr ; i++) {
        int w;
        int pos;
        rb_read_line(b);
        char * ptr = b->buf;
        char * nptr;
        w = strtoul(ptr, &nptr, 10);
        if (ptr == nptr) {
            fprintf(stderr, "Parse error while reading line %"PRIu32" of %s\n",
                    i+1, b->filename);
            exit(1);
        }
        ptr = nptr;
        row_table[i].i = i;
        row_table[i].w = w;
        line_bytes[i] = -b->o;
        if (i) {
            line_bytes[i-1] += b->o;
        }
        if (permute_cols) {
            /* Then we need to fill out the column information. This
             * implies parsing rows of the matrix unfortunately.
             */
            int nread;
            for(nread = 0 ; ; nread++) {
                uint32_t j;
                /* A column index + whitespace + exponent is never going to take
                 * more than 20 bytes */
                pos = ptr - b->buf;
                rb_feed_buffer_again_if_lowwater(b,&pos,20);
                ptr = b->buf + pos;
                /* remove leading ws */
                for( ; *ptr && (isspace(*ptr) && *ptr != '\n') ; ptr++) ;
                if (*ptr == '\n')
                    break;

                j = strtoul(ptr, &nptr, 10);
                if (nptr == ptr) {
                    fprintf(stderr, "Parse error while reading line %"
                            PRIu32" of %s\n", i+1, b->filename);
                    exit(1);
                }
                ptr = nptr;

                assert(j < nc);
                col_table[j].w++;
                /* discard as well the possible exponent information */
                for( ; *ptr && !isspace(*ptr) ; ptr++) ;
            }
            ASSERT_ALWAYS(nread == w);
        }

        rb_gobble_long_line(b);
    }
    line_bytes[i-1] += b->o;
    for( ; i < nmax ; i++) {
        row_table[i].i = i;
        row_table[i].w = 0;
        line_bytes[i] = 0;
    }

    /*
    fprintf(stderr, "done\n");
    fflush(stderr);
    */

    rb_close(b);

    read_matrix_trailer();
    if (pad_to_square) nr = nc = nmax;

    return 0;
}


void read_shuffled_matrix(const char * bigfile, const char * base)
{
    FILE * info;
    FILE * big;
    char * infoname;
    char rbuf[BIGCHUNK];

    asprintf(&infoname, "%s.info", base);

    info = fopen(infoname, "r");
    DIE_ERRNO_DIAG(info == NULL, "fopen", infoname);

    big = fopen(bigfile, "w");
    DIE_ERRNO_DIAG(big == NULL, "fopen", bigfile);

    fgets(rbuf, sizeof(rbuf), info);
    if (strncmp(rbuf, "//", 2) == 0) {
        int rc = sscanf(rbuf, "// %" SCNu32 " ROWS %" SCNu32, &nr, &nc);
        PARSE_CHECK(rc < 2, "header", rbuf);
    } else {
        int rc = sscanf(rbuf, "%" SCNu32 " %" SCNu32, &nr, &nc);
        PARSE_CHECK(rc < 2, "header", rbuf);
    }
    header_bytes = ftello(info);
    header = strdup(rbuf);
    // Of course, the new header does _NOT_ get copied to the big file !

    update_header();

    unsigned int nhs, nvs;
    fscanf(info, "%u %u", &nhs, &nvs);
    fprintf(stderr, "Existing matrix is split %ux%u\n", nhs, nvs);

    struct onefile {
        unsigned int i0;
        unsigned int j0;
        unsigned int i1;
        unsigned int j1;
        unsigned int ncoeffs;
        char locfile[80];
    };

    struct onefile ** files;

    unsigned int si;
    unsigned int sj;

    files = malloc(nhs * sizeof(struct onefile*));
    for(si = 0 ; si < nhs ; si++) {
        files[si] = malloc(nvs * sizeof(struct onefile));
    }

    for(si = 0 ; si < nhs ; si++) {
        for(sj = 0 ; sj < nvs ; sj++) {
            int rc;
            rc = fscanf(info, "%u %u", &si, &sj);
            FATAL_ERROR_CHECK(rc != 2, "parse error in .info file");
            rc = fscanf(info, "%u %u %u %u %u %80s\n",
                    &files[si][sj].i0,
                    &files[si][sj].j0,
                    &files[si][sj].i1,
                    &files[si][sj].j1,
                    &files[si][sj].ncoeffs,
                    files[si][sj].locfile);
            FATAL_ERROR_CHECK(rc != 6, "parse error in .info file");
            files[si][sj].i1 += files[si][sj].i0;
            files[si][sj].j1 += files[si][sj].j0;
        }
    }
    fclose(info);

    if (permute_rows) {
        row_table = (struct row*) malloc(nmax * sizeof(struct row));
        if (row_table == NULL) abort();
    }
    if (permute_cols) {
        col_table = (struct row *) malloc(nmax * sizeof(struct row));
        if (col_table == NULL) abort();
    }

    if (permute_cols) {
        unsigned int k;
        memset(col_table, 0, nmax * sizeof(struct row));
        for(k = 0 ; k < nmax ; k++) {
            col_table[k].i=k;
        }
    }

    line_bytes = (size_t *) malloc(nmax * sizeof(size_t));
    if (line_bytes == NULL) abort();

    assert(permute_rows);

    fprintf(stderr, "Reading matrix (first pass, counting weights)\n");

    FILE ** fps;
    unsigned int * fws;
    fps = malloc(nvs * sizeof(FILE *));
    fws = malloc(nvs * sizeof(unsigned int));

    off_t pos = ftello(big);

    unsigned int i = 0;

    char * input_dirname = my_dirname(base);

    for(si = 0 ; si < nhs ; si++) {
        for(sj = 0 ; sj < nvs ; sj++) {
            const char * nm = files[si][sj].locfile;
            if (input_dirname) {
                snprintf(rbuf,sizeof(rbuf),"%s/%s", input_dirname, nm);
                nm = rbuf;
            }
            fps[sj] = fopen(nm, "r");
            DIE_ERRNO_DIAG(fps[sj] == NULL, "fopen", nm);
            /* discard the header. */
            fgets(rbuf, sizeof(rbuf), fps[sj]);
            // fos[sj] = ftello(fps[sj]);
        }

        ASSERT_ALWAYS(i == files[si][0].i0);

        for( ; i < files[si][0].i1 ; i++) {
            unsigned int w = 0;
            for(sj = 0 ; sj < nvs ; sj++) {
                int rc = fscanf(fps[sj], "%u", &fws[sj]);
                FATAL_ERROR_CHECK(rc != 1, "parsing split files");
                w += fws[sj];
            }

            fprintf(big, "%u", w);
            for(sj = 0 ; sj < nvs ; sj++) {
                for(unsigned int k = 0 ; k < fws[sj] ; k++) {
                    unsigned int index;
                    int rc = fscanf(fps[sj], "%u", &index);
                    FATAL_ERROR_CHECK(rc != 1, "parsing split files");
                    index += files[si][sj].j0;
                    fprintf(big, " %u", index);
                    if (permute_cols) {
                        col_table[index].w++;
                    }
                    for(char c ; !isspace(c = fgetc(fps[sj])) ; fputc(c, big));
                }
            }
            fputc('\n', big);
            row_table[i].i = i;
            row_table[i].w = w;
            off_t npos = ftello(big);
            line_bytes[i] = npos - pos;
            pos = npos;
        }

        for(sj = 0 ; sj < nvs ; sj++) {
            fclose(fps[sj]);
            if (remove_input) {
                const char * nm = files[si][sj].locfile;
                if (input_dirname) {
                    snprintf(rbuf,sizeof(rbuf),"%s/%s", input_dirname, nm);
                    nm = rbuf;
                }
                unlink(nm);
            }
        }
    }

    free(input_dirname);

    for( ; i < nmax ; i++) {
        row_table[i].i = i;
        row_table[i].w = 0;
        line_bytes[i] = 0;
    }

    free(fws);
    free(fps);

    for(si = 0 ; si < nhs ; si++) {
        free(files[si]);
    }
    free(files);

    fclose(big);

    read_matrix_trailer();
    if (pad_to_square) nr = nc = nmax;

    if (remove_input) {
        unlink(infoname);
    }
    free(infoname);
}


/* }}} */

/*{{{ utilities for shuffle_rtable */
int decr_weight_cmp(const struct row * a, const struct row * b)
{
    int dw = b->w - a->w;
    if (dw) return dw;
    if (a->i < b->i)
        return 1;
    else if (b->i < a->i)
        return -1;
    return 0;
}

#if 0
int row_compare_index(const struct row * a, const struct row * b)
{
    return a->i - b->i;
}
#endif


/*}}}*/

/* {{{ shuffle_rtable: This is the basic procedure which dispatches rows
 * (or columns) in buckets according to our preferred strategy
 */
struct slice * shuffle_rtable(
        const char * text,
        struct row * rt,
        uint32_t n,
        unsigned int ns)
{
    uint32_t i;
    uint32_t ni;
    struct bucket * heap;
    struct slice * slices;

    qsort(rt, n, sizeof(struct row), (sortfunc_t) &decr_weight_cmp);

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
    for(i = 0 ; i < n ; i = ni) {
        int chunk;
        /* If it turns out that we have a great many rows of the same
         * weight (this sort of thing happens), then we'll dispatch them
         * in large chunks and not in a completely alternating way */
        for(ni = i + 1 ; ni < n && rt[i].w == rt[ni].w ; ni++);
        chunk = (ni-i) / ns;
        if (chunk) {
            int oversize=0;
            unsigned int k;
            for(k = 0 ; k < ns ; k++) {
                int l;
                for(l = 0 ; l < chunk && heap[k].room ; l++) {
                    int j = heap[k].i;
                    int pos = slices[j].nrows-heap[k].room;
                    slices[j].r[pos]=rt[i].i;
                    heap[k].s += rt[i].w;
                    heap[k].room--;
                    i++;
                }
                /* If one bucket becomes full, then this will change the
                 * comparison order inevitably */
                if (heap[k].room == 0)
                    oversize=1;
            }
            if (oversize)
                make_heap(heap, heap + ns);
        }
        for( ; i < ni ; i++) {
            int j = heap[0].i;
            int pos = slices[j].nrows-heap[0].room;
            assert(heap[0].room);
            slices[j].r[pos] = rt[i].i;
            heap[0].s += rt[i].w;
            heap[0].room--;
            pop_heap(heap, heap + ns);
            push_heap(heap, heap + ns);
        }
    }

    qsort(heap, ns, sizeof(struct bucket), (sortfunc_t) heap_index_compare);

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

struct slice * replicate_permutation(
        const char * text MAYBE_UNUSED,
        struct row * rt,
        uint32_t n,
        unsigned int ns,
        const struct slice * reference,
        unsigned int ns_reference)
{
    /* pick the dispatching from the reference array, and output a
     * slicing which corresponds to the same permutation, but in the
     * other dimension.
     */
    ASSERT_ALWAYS(nr == nc);
    struct slice * slices;
    uint32_t i, j, jr, k, kr;
    slices = alloc_slices(n, ns);

    j = k = 0;
    for(kr = 0 ; kr < ns_reference ; kr++) {
        for(jr = 0 ; jr < reference[kr].nrows ; jr++) {
            i = reference[kr].r[jr];
            if (j == slices[k].nrows) {
                k++;
                j=0;
            }
            slices[k].r[j]=i;
            assert(rt[i].i == i);
            slices[k].coeffs+=rt[i].w;
            j++;
        }
    }

    for(i = 0 ; i < ns ; i++) {
        printf("%s slice %d, span=%" PRIu32", weight=%"PRIu32"\n",
                text,
                i, slices[i].nrows,
                slices[i].coeffs);
    }

    uint32_t i0 = 0;
    for(i = 0 ; i < ns ; i++) {
        slices[i].i0=i0;
        i0 += slices[i].nrows;
    }

    return slices;
}

/* }}} */


void compute_permutation()
{
    if (replicate_rows_perm_for_columns && permute_rows) {
        row_slices = shuffle_rtable("horizontal", row_table, nr, nhslices);
        col_slices = replicate_permutation("vertical", col_table, nc,
                nvslices, row_slices, nvslices);
    } else if (replicate_columns_perm_for_rows && permute_cols) {
        col_slices = shuffle_rtable("vertical", col_table, nc, nvslices);
        row_slices = replicate_permutation("horizontal", row_table, nr,
                nhslices, col_slices, nvslices);
    } else {
        if (permute_rows) {
            row_slices = shuffle_rtable("horizontal", row_table, nr, nhslices);
        }
        if (permute_cols) {
            col_slices = shuffle_rtable("vertical", col_table, nc, nvslices);
        }
    }
    if (nhslices == 1 && !permute_rows && permute_cols) {
        row_slices = alloc_slices(nr, nhslices);
        unsigned int i;
        for(i = 0 ; i < nr ; i++) {
            row_slices[0].r[i]=i;
        }
        row_slices[0].i0 = 0;
        row_slices[0].nrows = nr;
    }
    final_file_info = malloc(nhslices*nvslices*sizeof(file_info));
    memset(final_file_info, 0, nhslices*nvslices*sizeof(file_info));
}


#if 0 /* {{{ old hslices code */
void open_hslices_files(const char * oname, const char * wname,
        const char * key,
        char *** pnames, FILE *** g,
        int n, const char * mode)
{
    int ii;
    char ** names;
    names = (char **) malloc(n * sizeof(char *));
    if (pnames) *pnames = names;
    if (g) *g = (FILE **) malloc(n * sizeof(FILE *));
    if (n == 1) {
        ii=0;
        if (strcmp(mode, "w") == 0) {
            fprintf(stderr, "No point in writing one slice only\n");
            /* Not only does it seem pointless, but it would as well be
             * dangerous because of the risk of erasing the original
             * input file. */
            exit(1);
        }
        names[ii] = strdup(oname);
        if (g) {
            (*g)[ii] = fopen(names[ii], mode);
            DIE_ERRNO_DIAG((*g)[ii] == NULL, "fopen", names[ii]);
        }
        return;
    }
    for(ii = 0 ; ii < n ; ii++) {
        size_t sz = strlen(wname);
        names[ii] = (char *) malloc(sz + 20);
        snprintf(names[ii], sz + 20, "%s.%s%d", wname, key, ii);
        if (g) {
            (*g)[ii] = fopen(names[ii], mode);
            DIE_ERRNO_DIAG((*g)[ii] == NULL, "fopen", names[ii]);
        }
    }
    return;
}

void close_hslices_files(char *** names, FILE *** g, int n)
{
    int ii;
    if (g) {
        for(ii = 0 ; ii < n ; ii++) {
            fclose((*g)[ii]);
        }
        free(*g);
    }

    if (names) {
        for(ii = 0 ; ii < n ; ii++) {
            free((*names)[ii]);
        }
        free(*names);
    }
}
#endif /* }}} */

void dispatcher(
        fileset dst, fileset src,
        size_t * sizes,
        unsigned int *dispatch,
        unsigned int nr)
{
    FILE ** g;
    unsigned int ii;
    FILE * f;
    unsigned int i;

    assert(src->n == 1);

    copybuf cb;

    cb_init(cb);
    cb_ensure(cb, BIGCHUNK);

    f = fileset_open_one(src, 0, "r");
    g = fileset_open(dst, "w");

    for(i = 0 ; i < nr ; i++) {
        ii = dispatch[i];
        size_t sz = sizes[i];
        cb_ensure(cb, sz);
        if (sz == 0) {
            /* special provision for writing zero rows */
            fprintf(g[ii], "0\n");
            continue;
        }
        fread(cb->buf, 1, sz, f);
        assert(cb->buf[sz-1] == '\n');
        fwrite(cb->buf, 1, sz, g[ii]);
    }
    fclose(f);
    cb_clear(cb);

    /* We've now done the dispatching, so the input file may be removed
     * if it's a temporary */
    fileset_dispose(src, 0);

    fileset_close(dst, g);
}


/* We're only doing this if nhslices > 1 or if we are asked to size-sort
 * rows
 */
void dispatch_row_slices(fileset fs)
{
    uint32_t i;
    unsigned int * dispatch;

    dispatch = (unsigned int *) malloc(nr * sizeof(unsigned int));
    for(i = 0 ; i < nhslices ; i++) {
        uint32_t k;
        for(k = 0 ; k < row_slices[i].nrows ; k++) {
            dispatch[row_slices[i].r[k]]=i;
        }
    }

    fileset fstmp;
    fileset_init_transfer(fstmp, nhslices, "h", fs, 0, TEMP);

    dispatcher(fstmp, fs, line_bytes, dispatch, nr);

    fileset_swap(fs, fstmp);
    fileset_clear(fstmp);

    free(dispatch);
}

/*{{{ write_permutation */
void write_permutation(const char * filename,
        struct slice * sl, unsigned int ns /*, int already_sorted */)
{
    FILE * f;
    uint32_t i,ii;
    f = fopen(filename, "w");
    i = 0;

    for(ii = 0 ; ii < ns ; ii++) {
        const struct slice * r = &(sl[ii]);
        for(i = 0 ; i < r->nrows ; i++) {
            /* This assertion is bogus. What really matters is the
             * weight, not the index in the file. So here, the data we
             * access might show a drift, while actually there is none
             * because the two rows have the same weight.
            if (already_sorted && i) {
                ASSERT_ALWAYS(r->r[i] > r->r[i-1]);
            }
            */
            // fprintf(f, "%" PRIu32 "\n", r->r[i]);
        }
        fwrite(r->r, sizeof(uint32_t), r->nrows, f);
    }
    fclose(f);
}
/*}}}*/

/*{{{ read_permutation */
unsigned int * read_permutation(const char * name, unsigned int n)
{
    FILE * f;
    unsigned int * res;
    if ((f = fopen(name, "r")) == NULL)
        return NULL;

    res = malloc(n * sizeof(unsigned int));
    unsigned int rc = fread(res, sizeof(unsigned int), n, f);
    FATAL_ERROR_CHECK(rc != n, "short read");
    return res;
}
/* }}} */

/*{{{ inverse_permutation */
void inverse_permutation(unsigned int * s, unsigned int n)
{
    if (s == NULL)
        return;

    unsigned int * t;

    t = malloc(n * sizeof(unsigned int));
    unsigned int i;
    for(i = 0 ; i < n ; i++) {
        ASSERT(s[i] < n);
        t[s[i]]=i;
    }
    memcpy(s,t,n * sizeof(unsigned int));
    free(t);
}
/* }}} */

void reindex_slices(struct slice * sl, unsigned int ns, unsigned int * s)
{
    if (s == NULL)
        return;

    uint32_t i,ii;
    for(ii = 0 ; ii < ns ; ii++) {
        const struct slice * r = &(sl[ii]);
        for(i = 0 ; i < r->nrows ; i++) {
            r->r[i] = s[r->r[i]];
        }
    }
}



/* {{{ comparison functions for weight_sort_hslice */
struct idxpair {
    unsigned int x;
    unsigned int y;
};

int idxcmp2(const struct idxpair * a, const struct idxpair * b)
{
    if (a->y < b->y) {
        return -1;
    } else if (b->y < a->y) {
        return 1;
    }
    return 0;
}
/* }}} */

/* {{{ logic for deciding how the split into RAM-sized chunk goes for
 * sorting */
struct chunk_info_s {
    unsigned int i0;
    unsigned int i1;
    size_t sz;
};
typedef struct chunk_info_s chunk_info[1];
typedef struct chunk_info_s * chunk_info_ptr;
struct chunk_splitting_s {
    chunk_info * c;
    unsigned int n;
};
typedef struct chunk_splitting_s chunk_splitting[1];

void clear_splitting(chunk_splitting s)
{
    free(s->c);
}

/* Decide in how many chunks we would have to split a slice in order to
 * perform sorting. This depends on the RAM limit which has been set.
 */
void compute_chunk_splitting(chunk_splitting s, unsigned int ii)
{
    size_t tot = 0;
    unsigned int i;
    unsigned int nrs = row_slices[ii].nrows;

    /* Split the file into RAM-sized chunks */
    for(i = 0 ; i < nrs ; i++) {
        tot += line_bytes[row_slices[ii].r[i]];
    }

    unsigned int nsubchunks = ((tot+ram_limit-1)/ram_limit);
    unsigned int nsc_max = (2*nsubchunks+1); // allow round-off error.

    s->c = (chunk_info *) malloc(nsc_max*sizeof(chunk_info));
    s->n = nsubchunks;

    if (nsubchunks == 1) {
        printf("Sorting hslice %u directly into memory\n", ii);
        s->c[0]->i0 = 0;
        s->c[0]->sz = tot;
        s->c[0]->i1 = nrs;
        return;
    }

    /* The 1024 here should be sufficient to avoid round-off
     * misses. They would not be catastrophic anyway. It's
     * just a pain to have nsubchunks+1 cuts instead of
     * nsubchunks (especially when the last one really is so
     * tiny). */
    size_t my_ram_limit = 1024 + (tot+nsubchunks-1) / nsubchunks;
    unsigned int ram_chunk=0;
    s->c[0]->i0 = 0;
    tot=0;
    for(i = 0 ; i < nrs ; i++) {
        size_t d = line_bytes[row_slices[ii].r[i]];
        if (tot + d >= my_ram_limit) {
            s->c[ram_chunk]->i1 = i;
            s->c[ram_chunk]->sz = tot;
            ram_chunk++;
            ASSERT_ALWAYS(ram_chunk < nsc_max);
            s->c[ram_chunk]->i0 = i;
            tot = 0;
        }
        tot += d;
    }
    s->c[ram_chunk]->i1 = i;
    s->c[ram_chunk]->sz = tot;
    s->n = ram_chunk + 1;
    s->c = realloc(s->c, s->n * sizeof(chunk_info));

    printf("Splitting hslice %u into %u intermediary files for sorting\n",
            ii, s->n);
}
/* }}} */

unsigned int * compute_dispatch_from_splitting(
        chunk_splitting s,
        unsigned int * shuffle,
        unsigned int nrs)
{
    unsigned int * dispatch;
    unsigned int k;

    if (s->n == 0)
        return NULL;
    dispatch = (unsigned int *) malloc(nrs * sizeof(unsigned int));
    for(k = 0 ; k < s->n ; k++) {
        unsigned int i;
        for(i = s->c[k]->i0 ; i < s->c[k]->i1 ; i++) {
            dispatch[shuffle[i]]=k;
        }
    }

    return dispatch;
}


/* This function returns the exact permutation of rows which is presently
 * relected in the horizontal slices (there could be only one such slice,
 * in which case it's the identity permutation).
 */
uint32_t * rowperm_after_dispatch(unsigned int ii)
{
    uint32_t nrs = row_slices[ii].nrows;
    unsigned int i;

    /* At this moment, we know that rows indexed within
     * row_slices[ii].data[] are effectively in fs->names[ii].
     * However, row at position k in this file is just the
     * k-earliest row from the input file that actually goes to
     * this slice.
     *
     * We need to fix this. Therefore we need to re-establish a
     * proper mapping that will give us information about the
     * rows currently in the slice, and where they will go.
     */

    struct idxpair * wherefrom;
    wherefrom = (struct idxpair *) malloc(nrs*sizeof(struct idxpair));
    for(i = 0 ; i < nrs ; i++) {
        wherefrom[i].x=i;
        wherefrom[i].y=row_slices[ii].r[i];
    }
    qsort(wherefrom, nrs, sizeof(struct idxpair), (sortfunc_t) idxcmp2);
    /* now wherefrom[kk].x is the index in the current file of the row
     * that will go to position [kk] eventually.
     * wherefrom[kk].y is the index of that same row in the
     * original big file */

    uint32_t * shuffle = malloc(nrs * sizeof(uint32_t));
    for(i = 0 ; i < nrs ; i++) {
        shuffle[wherefrom[i].x]=i;
    }
    free(wherefrom); wherefrom = NULL;

    return shuffle;
}

/* {{{ sink -- the way we eventually process the final data that will go
 * to the output files */

struct sink_s {
    fileset out;
    unsigned int last;
    /* fileset for final repartition of output files */
    fileset_ptr fh;
    fileset fv;
    FILE ** curr;
    FILE * f;   /* only for legacy jumbo output */
    /* Which row slice is being processed ? */
    unsigned int hnum;
    /* Only if we do column manipulation */
    uint32_t * cdirect;
    cperm_buf cp;
    /* this boolean triggers configuration, to happen as late as possible
     */
    int configured;
};
typedef struct sink_s * sink_ptr;
typedef struct sink_s sink[1];

void sink_configure(sink s)
{
    unsigned int jj;
    if (s->configured)
        return;
    if (permute_cols) {
        cperm_init(s->cp, nvslices);
        s->cdirect = (uint32_t *) malloc(nc * sizeof(uint32_t));
        memset(s->cdirect, 0, nc * sizeof(uint32_t));
        for(jj = 0 ; jj < nvslices ; jj++) {
            const struct slice * r = &(col_slices[jj]);
            unsigned int i;
            for(i = 0 ; i < r->nrows ; i++) {
                s->cdirect[r->r[i]]=r->i0+i;
            }
        }
    }
    s->configured = 1;
}

void sink_init(sink s, fileset_ptr fh)
{
    memset(s, 0, sizeof(sink));
    s->fh = fh;
}

void sink_hook(sink s, unsigned int i0, unsigned int i1)
{
    ASSERT_ALWAYS(i0 == s->last);

    if (i0 == 0) {
        sink_configure(s);
    }
    if (s->hnum == nhslices) {
        ASSERT_ALWAYS(i0 == i1);
        ASSERT_ALWAYS(i0 == nr);
        return;
    }
    if (legacy) {
        if (i0 == 0) {
            s->f = fopen(working_filename, "w");
            fprintf(s->f, "%s", header);
        }
        if (i0 == nr) {
            fclose(s->f);
            s->f = NULL;
        }
    }
    if (i0 == row_slices[s->hnum].i0 + row_slices[s->hnum].nrows) {
        if (!legacy) {
            fileset_close(s->fv, s->curr);
            fileset_clear(s->fv);
        }
        s->hnum++;
    }
    if (s->hnum == nhslices) {
        ASSERT_ALWAYS(i0 == i1);
        return;
    }
    if (i0 == row_slices[s->hnum].i0 && !legacy) {
        unsigned int i;
        fileset_init_transfer(s->fv, nvslices, "v", s->fh, s->hnum, OUTPUT);
        s->curr = fileset_open(s->fv, "w");
        for(i = 0 ; i < nvslices ; i++) {
            /* Print a header which could help for transposing */
            fprintf(s->curr[i], "%u %u\n", row_slices[s->hnum].nrows,
                    col_slices ? col_slices[i].nrows : nc);
            final_file_info[nvslices*s->hnum + i]->s = strdup(s->fv->names[i]);
        }
    }
}

void sink_shunt(sink s)
{
    s->last = nr;
    s->hnum = nhslices;
}

void sink_clear(sink s)
{
    assert(s->configured);
    if (permute_cols) {
        cperm_clear(s->cp);
        free(s->cdirect);
    }
    ASSERT_ALWAYS(s->last == nr);
    sink_hook(s, nr, nr);
    ASSERT_ALWAYS(s->hnum == nhslices);
}


/* }}} */

/* This sink version is eventually used when column manipulation is
 * performed */
void sink_feed_row(sink s,
        const char * data0,
        size_t sz)
{
    unsigned int nj = 0;
    unsigned int jj;
    const char * data = data0;
    char * ndata;
    cperm_reset(s->cp);

    if (sz) {
        /* Don't use sscanf here. Its complexity is
         * linear in strlen() of its argument */
        nj=strtoul(data, &ndata, 10);
        PARSE_CHECK(data==ndata, "row", strndup(data0,sz));
        data=ndata;
        for(jj = 0 ; jj < nj ; jj++) {
            unsigned int j;
            assert(*data == ' ');
            data++;
            j=strtoul(data, &ndata, 10);
            PARSE_CHECK(data==ndata, "row", strndup(data0, sz));
            data=ndata;
            unsigned int w;
            char scrap[10];
            w = which_slice(col_slices, nvslices, s->cdirect[j]);
            snprintf(scrap, sizeof(scrap), " %u",
                    s->cdirect[j]-(legacy?0:col_slices[w].i0));
            cperm_append(s->cp, w, scrap, strlen(scrap));
            if (*data == ':') {
                int e;
                data++;
                e=strtol(data, &ndata, 10);
                PARSE_CHECK(data==ndata, "row",
                        strndup(data0,sz));
                cperm_append(s->cp, w, data-1, ndata-data+1);
                data = ndata;
            }
            s->cp->w[w]++;
            final_file_info[nvslices*s->hnum + w]->coeffs++;
        }
        if (*data == ' ') { data++; }
        assert(*data == '\n');
        data=NULL;
    }

    if (legacy) {
        fprintf(s->f, "%u", nj);
        for(jj = 0 ; jj < nvslices ; jj++) {
            fwrite(s->cp->bufs[jj]->buf, 1, s->cp->pos[jj], s->f);
        }
        fputc('\n', s->f);
    } else {
        for(jj = 0 ; jj < nvslices ; jj++) {
            fprintf(s->curr[jj], "%u", s->cp->w[jj]);
            fwrite(s->cp->bufs[jj]->buf, 1, s->cp->pos[jj], s->curr[jj]);
            fputc('\n', s->curr[jj]);
        }
    }
#if 0
    if (nvslices > 1) {
        for(jj = 0 ; jj < nvslices ; jj++) {
            fprintf(fv[jj], "%u", s->cp->w[jj]);
            fwrite(s->cp->bufs[jj]->buf, 1, s->cp->pos[jj], fv[jj]);
            fputc('\n', fv[jj]);
        }
    } else {
        fprintf(f, "%u", nj);
        fwrite(s->cp->bufs[0]->buf, 1, s->cp->pos[0], f);
        fputc('\n', f);
    }
#endif
    s->last++;
}

void sink_feed_memory(sink s,
        const char * inmem,
        unsigned int ii,
        chunk_info_ptr sc)
{
    sink_hook(s, row_slices[ii].i0 + sc->i0, row_slices[ii].i0 + sc->i1);
    ASSERT_ALWAYS(s->hnum == ii);
    if (permute_cols) {
        unsigned int i;
        for(i = sc->i0 ; i < sc->i1 ; i++) {
            size_t sz = line_bytes[row_slices[s->hnum].r[i]];
            sink_feed_row(s, inmem, sz);
            inmem += sz;
        }
    } else {
        assert(nvslices == 1);
        if (legacy) {
            fwrite(inmem, 1, sc->sz, s->f);
        } else {
            fwrite(inmem, 1, sc->sz, s->curr[0]);
        }
        s->last = row_slices[ii].i0 + sc->i1;
    }
}


void weight_sort_hslice(sink datasink, fileset fs, unsigned int ii)
{
    copybuf cb;
    cb_init(cb);
    cb_ensure(cb, BIGCHUNK);


    uint32_t nrs = row_slices[ii].nrows;
    uint32_t * shuffle = rowperm_after_dispatch(ii);
    unsigned int i;

    /* We now know that row shuffle[j] in the hslice number ii is in fact
     * the one whose real index is row_slices[ii].r[j]
     *
     * our purpose is to make it go to position j in that same file
     */

    chunk_splitting split;

    compute_chunk_splitting(split, ii);


    /* To use the dispatcher() routine, we must create a
     * dispatch[] array, and an nbytes[] one.
     * Note that while both are available, one could do without
     * if needed:
     * nbytes[] could be obtained on the fly
     * dispatch[] is really so trivial (mere index comparison)...
     */

    size_t * nbytes = malloc(nrs * sizeof(size_t));
    uint32_t tot = 0;
    for(i = 0 ; i < nrs ; i++) {
        nbytes[shuffle[i]]=line_bytes[row_slices[ii].r[i]];
    }

    unsigned int * dispatch;
    dispatch = compute_dispatch_from_splitting(split, shuffle, nrs);

    fileset fs_local;
    fileset_init_isolate(fs_local, fs, ii);

    if (split->n > 1) {
        fileset fstmp;
        fileset_init_transfer(fstmp, split->n, "sort", fs, ii, TEMP);
        dispatcher(fstmp, fs_local, nbytes, dispatch, nrs);
        fileset_swap(fs_local, fstmp);
        fileset_clear(fstmp);
    }

    /* We'll now need something else: an array indicating the
     * offsets at which the input rows in the sub-chunks must be
     * put. This array, at first, is indexed by destination row.
     */
    off_t * poking_in_ordered_file;
    poking_in_ordered_file = (off_t *) malloc(nrs * sizeof(off_t));
    tot=0;
    for(i = 0 ; i < nrs ; i++) {
        size_t d = nbytes[shuffle[i]];
        poking_in_ordered_file[i]=tot;
        tot += d;
    }

    off_t * poking_place = (off_t *) malloc(nrs * sizeof(off_t));
    int * tmp = (int *) malloc(split->n*sizeof(int));
    unsigned int s;
    for(s = 0 ; s < split->n ; s++) tmp[s] = split->c[s]->i0;
    tot=0;
    /* Hard to avoid for the computation to come */
    unsigned int * direct = (unsigned int*) malloc(nrs * sizeof(unsigned int));
    for(i = 0 ; i < nrs ; i++) {
        direct[shuffle[i]]=i;
    }
    /* We'll now change the use of the nbytes[] array. Instead of
     * being previously indexed by row numbers in the source
     * file, it's now indexed by row number in the concatenation
     * of the sub-chunks. This means that if the input file has
     * row r4,r2,r1,r3, and that the sub-chunks are (r4,r1) and
     * then (r2,r3), then the new nbytes will be s4,s1,s2,s3
     * instead of s4,s2,s1,s3 (s_i is the byte size of r_i).
     */
    /* To do so, we look at rows starting from the ordered
     * version. We'll maintain one ``current index'' for each
     * sub-chunk. This amounts to emulating the dispatching
     * again, but only in-memory. */

    for(s = 0 ; s < split->n ; s++) tmp[s] = split->c[s]->i0;
    for(i = 0 ; i < nrs ; i++) {
        unsigned int which = split->n > 1 ? dispatch[i] : 0;
        ASSERT_ALWAYS(which < split->n);
        /* This row is in sub chunk-[[which]] */
        /* The next unassigned byte-size from that sub-chunk is
         * nbytes[tmp[which]].
         */
        size_t sz = line_bytes[row_slices[ii].r[direct[i]]];
        /* Note that it is
         * important that we process this array in direct order.
         */
        nbytes[tmp[which]] = sz;
        poking_place[tmp[which]]=
            poking_in_ordered_file[direct[i]] -
            poking_in_ordered_file[split->c[which]->i0];
        tmp[which]++;
    }

    /* Free as much as we can before doing the large-mem stuff */
    free(tmp); tmp=NULL;
    free(direct); direct = NULL;
    free(poking_in_ordered_file); poking_in_ordered_file = NULL;
    free(shuffle); shuffle = NULL;
    free(dispatch); dispatch = NULL;

    for(s = 0 ; s < split->n ; s++) {
        FILE * g = fileset_open_one(fs_local, s, "r");
        chunk_info_ptr sc = split->c[s];
        char * inmem = malloc(sc->sz);
        for(i = sc->i0 ; i < sc->i1 ; i++) {
            off_t sz = nbytes[i];
            if (sz == 0)
                continue;
            /* We may have zero rows, in which case nothing is to be
             * written. We have sufficient information available for
             * knowing later on that zero rows should be written out
             */
            cb_ensure(cb, sz);
            fread(cb->buf, 1, sz, g);
            assert(cb->buf[sz-1] == '\n');
            assert(poking_place[i] >= 0);
            assert(poking_place[i] + sz <= (off_t) sc->sz);
            memcpy(inmem + poking_place[i], cb->buf, sz);
        }
        fileset_close_one(fs_local, s, g);
        /* At this point, we have a sorted part of an horizontal hslice
         * available. We are now ready for producing the final output
         * file, but this could mean several different things.
         *
         * - If a column splitting is expected, then all rows must go in
         *   order towards the final files.
         *
         * - If only a column sort is wanted, then there's hardly any
         *   difference. Only the naming of the output file eventually
         *   changes.
         *
         * - In case we're not going to do anything with columns, then
         *   the whole bunch of data can be treated in one go.
         */
        fileset_dispose(fs_local, s);
        /* It's important to call fileset_dispose _before_ the sink is
         * activated, because doing the converse could lead to late
         * deletion of a file which should actually not be deleted.
         */

        sink_feed_memory(datasink, inmem, ii, sc);
        free(inmem);
    }

    clear_splitting(split);
    free(poking_place); poking_place = NULL;
    free(nbytes);
}

void do_per_hslice_stuff(sink datasink, fileset fs)
{
    if (weight_sort_in_cells && !rows_are_weight_sorted) {
        unsigned int ii;
        for(ii = 0 ; ii < nhslices ; ii++) {
            weight_sort_hslice(datasink, fs, ii);
        }
    } else if (permute_cols) {
        unsigned int ii;
        // (nvslices > 1 || (weight_sort_in_cells && !cols_are_weight_sorted)
        for(ii = 0 ; ii < nhslices ; ii++) {
            /* Then we'll go linearly, but this will imply column work
             * (hence row by row).
             */
            uint32_t * shuffle = rowperm_after_dispatch(ii);
            unsigned int nrs = row_slices[ii].nrows;
            size_t * nbytes = malloc(nrs * sizeof(size_t));
            FILE * f = fileset_open_one(fs, ii, "r");
            unsigned int i;

            struct slice * sl = &(row_slices[ii]);

            copybuf cb;
            cb_init(cb);
            cb_ensure(cb, BIGCHUNK);

            for(i = 0 ; i < nrs ; i++) {
                nbytes[shuffle[i]]=line_bytes[sl->r[i]];
            }
            sink_hook(datasink, sl->i0, sl->i0 + sl->nrows);
            for(i = 0 ; i < nrs ; i++) {
                size_t sz = nbytes[i];
                cb_ensure(cb, sz);
                if (sz == 0) {
                    const char zrow[] = "0\n";
                    sz = strlen(zrow);
                    cb_ensure(cb, sz);
                    strcpy(cb->buf, zrow);
                } else {
                    fread(cb->buf, 1, sz, f);
                }
                assert(cb->buf[sz-1] == '\n');
                sink_feed_row(datasink, cb->buf, sz);
            }
            fileset_close_one(fs, ii, f);
            fileset_dispose(fs, ii);
            cb_clear(cb);
            free(shuffle);
            free(nbytes);
        }
    } else {
        /* Well, easy !  */
        sink_shunt(datasink);
    }

    /* Just for fun, here's the (zsh) shell loop that I used to check
     * that this program is not completely nuts (on an older version
     * which did not subtract an offset in column indices within
     * vslices).
cat /tmp/mat.h0.debug* | (n=0; while read x ; do echo "row $n" >&2 ; echo $x | (read n a ; for x in ${=a} ; do grep -n '^'$x'$' /tmp/col_perm.txt | cut -d: -f1 | (read z ; expr $z - 1); done) | sort >! /tmp/orig-$n ; for i in {5..8} ; do read a b <&$i ; echo $b | xargs -n 1 echo ; done | sort >! /tmp/sorted-$n ; diff -u /tmp/orig-$n /tmp/sorted-$n ;  let n+=1; done 5</tmp/mat.h0.v0 6</tmp/mat.h0.v1 7</tmp/mat.h0.v2 8</tmp/mat.h0.v3)
     * (well, actually, _who_ is nuts ?)
     */

}

void write_info_file(int argc, char * argv[])
{
    char * info;
    FILE * f;
    unsigned int i,j;
    asprintf(&info, "%s.info", working_filename);
    f = fopen(info, "w");
    DIE_ERRNO_DIAG(f == NULL, "fopen", info);
    fprintf(f, "%s", header);
    fprintf(f, "%u %u\n", nhslices, nvslices);
    for(i = 0 ; i < nhslices ; i++) {
        for(j = 0 ; j < nvslices ; j++) {
            fprintf(f, "%u %u"  /* total i, j */
                    " %"PRIu32" %"PRIu32        /* i0 j0 for this block */
                    " %"PRIu32" %"PRIu32        /* nr nc for this block */
                    " %"PRIu32  /* ncoeffs */
                    " %s\n",    /* filename */
                    i,j,
                    row_slices ? row_slices[i].i0 : 0,
                    col_slices ? col_slices[j].i0 : 0,
                    row_slices ? row_slices[i].nrows : nr,
                    col_slices ? col_slices[j].nrows : nc,
                    final_file_info[i*nvslices+j]->coeffs,
                    my_basename(final_file_info[i*nvslices+j]->s));
        }
    }
    time_t t = time(NULL);
    fprintf(f, "# %s", ctime(&t));
#ifdef  REV
    fprintf(f, "# revision " REV "\n");
#else
    fprintf(f, "# revision (unspecified)\n");
#endif
    fprintf(f, "#");
    for(i = 0 ; i < (unsigned int) argc ; i++) {
        fprintf(f, " %s", argv[i]);
    }
    fprintf(f, "\n");
    fclose(f);
    free(info);
}
void cleanup()
{
    unsigned int i;
    if (row_slices)     free_slices(row_slices, nhslices);
    if (col_slices)     free_slices(col_slices, nvslices);
    free(header);
    free(line_bytes);
    for(i = 0 ; i < nhslices * nvslices ; i++) {
        free(final_file_info[i]->s);
    }
    free(final_file_info);
}

void writeout_both_permutations()
{
    if (permute_rows) {
        char * rp;
        asprintf(&rp, "%s.row_perm", working_filename);
        write_permutation(rp, row_slices, nhslices);
        // write_permutation(rp, row_slices, nhslices, rows_are_weight_sorted);
        free(rp);
    }
    if (permute_cols) {
        char * cp;
        asprintf(&cp, "%s.col_perm", working_filename);
        write_permutation(cp, col_slices, nvslices);
        // write_permutation(cp, col_slices, nvslices, cols_are_weight_sorted);
        free(cp);
    }
}

int main(int argc, char * argv[])
{
    int argc0 = argc;
    char ** argv0 = argv;
    param_list pl;

    param_list_init(pl);

    argv++,argc--;
    param_list_configure_knob(pl, "--legacy", &legacy);
    param_list_configure_knob(pl, "--remove-input", &remove_input);
    param_list_configure_knob(pl, "--pad", &pad_to_square);
    param_list_configure_alias(pl, "--pad", "--square");
    param_list_configure_alias(pl, "out", "output-name");
    param_list_configure_alias(pl, "in", "matrix");
    param_list_configure_alias(pl, "ram-limit", "ramlimit");
    param_list_configure_alias(pl, "nbuckets", "nslices");
    param_list_configure_alias(pl, "nbuckets", "slices");
    param_list_configure_knob(pl, "--keep-temps", &keep_temps);
    param_list_configure_knob(pl, "--permute-rows-like-columns",
            &replicate_columns_perm_for_rows);
    param_list_configure_knob(pl, "--permute-columns-like-rows",
            &replicate_rows_perm_for_columns);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage();
    }

    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "subdir")) != NULL) {
        fprintf(stderr, "--subdir is no longer supported."
                " prepend --out with a path instead.\n");
        exit(1);
#if 0
        /* chdir'ing to the destination directory when input files may
         * reside elsewhere is sort of buggy */
        int rc = chdir(tmp);
        if (rc < 0 && errno == ENOENT) {
            rc = mkdir(tmp, 0777);
            DIE_ERRNO_DIAG(rc < 0, "mkdir", tmp);
            rc = chdir(tmp);
            DIE_ERRNO_DIAG(rc < 0, "chdir", tmp);
        }
#endif
    }

    if ((tmp = param_list_lookup_string(pl, "ram-limit")) != NULL) {
        char * eptr;
        ram_limit = strtoul(tmp, &eptr, 10);
        if (ram_limit == 0) usage();
        switch (*eptr) {
            case 'g': case 'G': ram_limit <<= 10;
            case 'm': case 'M': ram_limit <<= 10;
            case 'k': case 'K': ram_limit <<= 10;
                                eptr++;
            default:
                                if (*eptr != '\0') usage();
        }
    }

    working_filename = strdup(param_list_lookup_string(pl, "out"));
    pristine_filename = strdup(param_list_lookup_string(pl, "in"));
    if (working_filename == NULL) {
        fprintf(stderr, "Required argument --out is missing\n");
        exit(1);
    }
    if (pristine_filename == NULL) {
        fprintf(stderr, "Required argument --out is missing\n");
        exit(1);
    }

    {
        char * output_dirname = my_dirname(working_filename);
        if (output_dirname) {
            if (access(output_dirname,X_OK) < 0) {
                int rc = mkdir(output_dirname, 0777);
                DIE_ERRNO_DIAG(rc < 0, "mkdir", output_dirname);
            }
            free(output_dirname);
        }
    }

    /* We accept a sole argument to mean n hslices and one vslice, or an
     * argument of the form MMMxNNN
     */
    int split[2] = {1,1};

    tmp = param_list_lookup_string(pl, "nbuckets");
    if (tmp == NULL) {
        fprintf(stderr, "Required argument --nslices is missing\n");
        exit(1);
    }

    if (strchr(tmp,'x')) {
        param_list_parse_intxint(pl, "nbuckets", split);
    } else {
        param_list_parse_int(pl, "nbuckets", &split[0]);
    }

    param_list_parse_int(pl, "hslices", &split[0]);
    param_list_parse_int(pl, "vslices", &split[1]);

    if (param_list_warn_unused(pl)) {
        usage();
    }
    param_list_clear(pl);

    nhslices = split[0];
    nvslices = split[1];


    /* 10 is overestimated here. I think 4 suffices: 0,1,2, and f.
     */
    if (sysconf(_SC_OPEN_MAX) <= (long int) (nhslices + 10)) {
        fprintf(stderr, "Not enough open files allowed (%ld)\n",
                sysconf(_SC_OPEN_MAX));
        exit(1);
    }

#if 0
    /* That's a possible cause of disaster, so we'll avoid this case */
    if (access("matrix.txt",F_OK) == 0) {
        fprintf(stderr, "matrix.txt already exists in destination directory\n");
        exit(1);
    }
#endif

    ASSERT_ALWAYS(replicate_columns_perm_for_rows + replicate_rows_perm_for_columns <= 1);

    weight_sort_in_cells = 1;

    permute_rows = weight_sort_in_cells || nhslices > 1;
    permute_cols = weight_sort_in_cells || nvslices > 1;

    if (replicate_columns_perm_for_rows) {
        permute_rows = permute_cols;
    }
    if (replicate_rows_perm_for_columns) {
        permute_cols = permute_rows;
    }

    if (!(permute_rows || permute_cols || weight_sort_in_cells)) {
        /* Then we've got absolutely nothing to do ! */
        fprintf(stderr, "Nothing to do ? Really ? Why not use cp then ?\n");
        exit(1);
    }

    int recycle_shuffled_matrix = 0;
    if (access(pristine_filename, F_OK) < 0) {
        char * tmp;
        asprintf(&tmp, "%s.info", pristine_filename);
        if (access(tmp, F_OK) < 0) {
            fprintf(stderr, "%s: file not found\n", pristine_filename);
            exit(1);
        }
        free(tmp);
        recycle_shuffled_matrix = 1;
    }

    fileset work;
    fileset_init(work, 1);

    if (recycle_shuffled_matrix) {
        /* special case. We rebuild a full matrix from the existing
         * version. Note that it is somewhat disk-intensive. we're not
         * doing it the totally stupid way, which would incur re-doing
         * balance in the reverse order. However, we're not doing it the
         * I/O-smart way, in that we create a temporary file which could
         * be avoided. This would be a nightmare, however.
         */
        asprintf(&(work->names[0]), "tmp-%s", working_filename);
        prefix_fixup(work->names[0],"tmp-");

        read_shuffled_matrix(work->names[0], pristine_filename);

        work->status = TEMP;

        char * oldrp = NULL;
        asprintf(&oldrp, "%s.row_perm", pristine_filename);
        unsigned int * old_rowperm = read_permutation(oldrp, nr);
        free(oldrp);

        char * oldcp = NULL;
        asprintf(&oldcp, "%s.col_perm", pristine_filename);
        unsigned int * old_colperm = read_permutation(oldcp, nc);
        free(oldcp);

        unsigned int i;

        /* Apply to row_table and col_table */
        if (old_rowperm) {
            for(i = 0 ; i < nr ; i++)
                row_table[i].i = old_rowperm[i];
        }
        if (old_colperm) {
            for(i = 0 ; i < nc ; i++)
                col_table[i].i = old_colperm[i];
        }

        /* Haven't dared to check whether the current code exploits these
         * flags correctly */
        rows_are_weight_sorted = 0;
        cols_are_weight_sorted = 0;

        fprintf(stderr, "Computing permutation\n");
        compute_permutation();
        if (row_table) free(row_table);
        if (col_table) free(col_table);

        writeout_both_permutations();

        /* Now permute back row_slices and col_slices */
        inverse_permutation(old_rowperm, nr);
        inverse_permutation(old_colperm, nc);

        reindex_slices(row_slices, nhslices, old_rowperm);
        reindex_slices(col_slices, nvslices, old_colperm);

        free(old_rowperm);
        free(old_colperm);
    } else {
        read_matrix();
        work->names[0] = strdup(pristine_filename);
        work->status = INPUT;

        fprintf(stderr, "Computing permutation\n");
        compute_permutation();
        if (row_table) free(row_table);
        if (col_table) free(col_table);
    }

    writeout_both_permutations();

    if (nhslices > 1) {
        fprintf(stderr, "Dispatching matrix in row slices\n");
        dispatch_row_slices(work);
    }

    sink datasink;

    sink_init(datasink, work);

    /* Now we're going to handle all the horizontal slices in turn. We'll
     * eventually take this opportunity to do actual sorting by weight if
     * required. But prior to that, we'll also split slices vertically
     * first. */
    do_per_hslice_stuff(datasink, work);

    sink_clear(datasink);
    fileset_clear(work);

    if (!legacy) {
        write_info_file(argc0, argv0);
    }

    free(pristine_filename);
    free(working_filename);

    fprintf(stderr, "Done\n");
    cleanup();
    return 0;
}
