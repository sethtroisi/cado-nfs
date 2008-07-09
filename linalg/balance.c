#define __STDC_FORMAT_MACROS    /* For C++ to have PRId64 and friends */
#define _GNU_SOURCE    /* for strndup only -- really easy to make a workalike */

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
#include <math.h>       /* XXX could probably remove it */

/* platform headers */
#include <sys/stat.h>   /* mkdir */
#include <sys/types.h>  /* mkdir */
#include <unistd.h>     /* chdir */


/*{{{ macros */
#ifdef  __cplusplus
using namespace std;
#endif

/* I'd prefer the program to stay completely standalone */
#ifndef ASSERT_ALWAYS
#define ASSERT_ALWAYS(x) do { if (!(x)) { fprintf(stderr, "Assertion " #x " failed\n"); abort(); } } while (0)
#endif
#if defined(__GNUC__)
#ifndef	MAYBE_UNUSED
#define MAYBE_UNUSED __attribute__ ((unused))
#endif
#else
#ifndef	MAYBE_UNUSED
#define MAYBE_UNUSED
#endif
#endif

#define DIE_ERRNO_DIAG(tst, func, arg) do {				\
    if ((tst)) {					        	\
        fprintf(stderr, func "(%s): %s\n", arg, strerror(errno));       \
        exit(1);					        	\
    }							        	\
} while (0)
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
unsigned int nhslices = 1;
unsigned int nvslices = 1;

int keep_temps = 0;
int remove_input = 0;
int pad_to_square = 0;
int legacy = 0;

int permute_rows = 0;
int permute_cols = 0;

size_t ram_limit = 1 << 28;     /* 256 MB */

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


/* full path to the input matrix */
const char * pristine_filename;

/* basename of the resulting matrix filename */
const char * working_filename = "matrix.txt";
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

/* These two really only have to do with optimizing I/O for the matrix.
 * In case we're able to skip parsing entirely, we save line sizes for
 * reuse.
 */
size_t * line_bytes;
size_t header_bytes;
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
/* {{{ reading buffer stuff */
struct reading_buffer_s {
    char buf[1024];
    FILE * f;
    int siz;
    off_t o;
    char * filename;
};

typedef struct reading_buffer_s reading_buffer[1];


static void rb_open(reading_buffer b, const char * filename)
{
    b->filename = strdup(filename);
    b->f = fopen(filename, "r");
    DIE_ERRNO_DIAG(b->f == NULL, "fopen", b->filename);
}
static void rb_close(reading_buffer b)
{
    fclose(b->f);
    free(b->filename);
}

static void rb_read_line(reading_buffer b)
{
    char * ptr;
    b->o = ftello(b->f);                                                  
    ptr = fgets(b->buf, sizeof(b->buf), b->f);
    if (ptr == NULL) {
        fprintf(stderr, "Unexpected %s at position %ld in %s\n",    
                feof(b->f) ? "EOF" : "error", b->o, b->filename);            
        exit(1);
    }
    b->siz = strlen(ptr);
}

static void rb_gobble_long_line(reading_buffer b)
{
    char * ptr;
    while (!(b->siz < (int) (sizeof(b->buf)-1) ||
                b->buf[sizeof(b->buf)-2] == '\n'))
    {
        ptr = fgets(b->buf, sizeof(b->buf), b->f);
        if (ptr == NULL) {
            fprintf(stderr, "Unexpected %s at position %ld in %s\n",    
                    feof(b->f) ? "EOF" : "error",
                    ftello(b->f), b->filename);    
            exit(1);
        }
        b->siz = strlen(ptr);
    }                                                                   
    b->o = ftello(b->f);
}

void rb_feed_buffer_again_if_lowwater(reading_buffer b,
        int * ppos, int low)
{
    char * ptr;
    if (*ppos+low <= (int) sizeof(b->buf) ||
            (b->siz < (int) (sizeof(b->buf)-1) ||
                b->buf[sizeof(b->buf)-2] == '\n')) {
        return;
    }

    memcpy(b->buf, b->buf+*ppos, b->siz-*ppos);
    b->siz = b->siz-*ppos;
    *ppos=0;
    ptr = fgets(b->buf+b->siz,sizeof(b->buf)-b->siz,b->f);
    if (ptr == NULL) {
        fprintf(stderr, "Unexpected %s at position %ld in %s\n",    
                feof(b->f) ? "EOF" : "error", ftello(b->f), b->filename);    
        exit(1);
    }
    b->siz += strlen(b->buf + b->siz);
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

#if 0
int has_prefix(const char * s, const char * pfx)
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
    memcpy(last_slash, last_slash + s1, strlen(last_slash + s1) + 1);
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
    for(k = 0 ; k < n ; k++) {
        uint32_t cf = data[k].w;
        if (data[k].w == 0) nz++;
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
        printf("%u zero %s\n", nz, text);
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

    /* nmax is the largest dimension. If we pad to square size, we will
     * work with a matrix of size nmax * nmax.
     */
    uint32_t nmax = nc < nr ? nr : nc;

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
        if (sscanf(b->buf, "%d%n", &w, &pos) < 1) {
            fprintf(stderr, "Parse error while reading line %"PRIu32" of %s\n",
                    i+1, b->filename);
            exit(1);
        }
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
                int dpos;
                int rc;
                /* A column index is never going to take more than that */
                rb_feed_buffer_again_if_lowwater(b,&pos,10);
                /* We know that leading whit space will be removed */
                rc = sscanf(b->buf+pos, "%" SCNu32 "%n", &j, &dpos);
                if (rc == EOF) {
                    break;
                }
                if (rc < 1) {
                    fprintf(stderr, "Parse error while reading line %"
                            PRIu32" of %s\n", i+1, b->filename);
                    exit(1);
                }

                assert(j < nc);
                col_table[j].w++;
                pos += dpos;
                /* discard as well the possible exponent information */
                for(;b->buf[pos]&&!isspace(b->buf[pos]);pos++);
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


    if (pad_to_square) {
        nr = nmax;
        nc = nmax;
        fprintf(stderr, "Output matrix will be padded to"
                " %"PRIu32" rows,"
                " %"PRIu32" columns\n",
                nr, nc);
    }

    return 0;
}
/* }}} */

/*{{{ utilities for shuffle_rtable */
int decr_weight_cmp(const struct row * a, const struct row * b)
{
    return b->w - a->w;
}

#if 0
int row_compare_index(const struct row * a, const struct row * b)
{
    return a->i - b->i;
}
#endif

/* {{{ Data type for buckets used in the priority queue stuff (w/ tests) */
/* We'll arrange our slices in a priority heap, so that we'll always
 * easily access the lightest one.
 *
 * The number of rows in a bucket does not really matter, except that we
 * absolutely want to make sure that it never exceeds the bucket size.
 * Hence if a bucket becomes full, then it's always considered heavier
 * than anything.
 * */

struct bucket {
    unsigned long s;    /* total weight */
    unsigned long room; /* nb of rows that can still be stored */
    int i;              /* bucket index */
};

int heap_index_compare(const struct bucket * a, const struct bucket * b)
{
    return (int32_t) a->i - (int32_t) b->i;
}

int heap_compare(const struct bucket * a, const struct bucket * b)
{
    /* Returns whether b is a better candidate to be on top of the heap
     * IOW: Whether b is lighter than a
     * */
    if (b->room == 0) return 0;
    if (a->room == 0) return 1;
    return b->s < a->s;
}
/* }}} */

/* The code below is the glibc heap implementation, as copied from
 * /usr/include/c++/4.1.2/bits/stl_heap.h */
/* {{{ */
void
push_heap0(struct bucket * __first,
        ptrdiff_t hole, ptrdiff_t __topIndex,
        const struct bucket * __value)
{
    ptrdiff_t __parent = (hole - 1) / 2;
    while (hole > __topIndex && heap_compare(__first + __parent, __value))
    {
        *(__first + hole) = *(__first + __parent);
        hole = __parent;
        __parent = (hole - 1) / 2;
    }
    *(__first + hole) = *__value;
}

static inline void
push_heap(struct bucket * __first, struct bucket * __last)
{
    struct bucket ref = __last[-1];
    push_heap0(__first, (__last - __first) - 1, 0, & ref);
}

void
adjust_heap(struct bucket * __first, ptrdiff_t hole,
        ptrdiff_t __len, const struct bucket * __value)
{
    const ptrdiff_t __topIndex = hole;
    ptrdiff_t right = 2 * hole + 2;
    while (right < __len)
    {
        if (heap_compare(__first + right,
                    __first + (right - 1)))
            right--;
        *(__first + hole) = *(__first + right);
        hole = right;
        right = 2 * (right + 1);
    }
    if (right == __len)
    {
        *(__first + hole) = *(__first + (right - 1));
        hole = right - 1;
    }
    push_heap0(__first, hole, __topIndex, __value);
}

static inline void
pop_heap(struct bucket * __first, struct bucket * __last)
{
    struct bucket ref = __last[-1];
    __last[-1] = __first[0];
    adjust_heap(__first, 0, __last - __first, & ref);
}

static inline void
make_heap(struct bucket * __first, struct bucket * __last)
{
    if (__last - __first < 2)
        return;

    const ptrdiff_t __len = __last - __first;
    ptrdiff_t __parent = (__len - 2) / 2;
    for(;;) {
        struct bucket ref = __first[__parent];
        adjust_heap(__first, __parent, __len, &ref);
        if (__parent == 0)
            return;
        __parent--;
    }
}
/* end of heap code */
/* }}} */

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
                oversize = l < chunk;
            }
            if (oversize)
                make_heap(heap, heap + ns);
        }
        for( ; i < ni ; i++) {
            int j = heap[0].i;
            int pos = slices[j].nrows-heap[0].room;
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
/* }}} */


void compute_permutation()
{
    if (permute_rows) {
        row_slices = shuffle_rtable("horizontal", row_table, nr, nhslices);
    }
    if (permute_cols) {
        col_slices = shuffle_rtable("vertical", col_table, nc, nvslices);
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
#define BIGCHUNK        16384
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
        struct slice * slice, unsigned int ns, int already_sorted)
{
    FILE * f;
    uint32_t i,ii;
    f = fopen(filename, "w");
    i = 0;

    for(ii = 0 ; ii < ns ; ii++) {
        const struct slice * r = &(slice[ii]);
        for(i = 0 ; i < r->nrows ; i++) {
            if (already_sorted && i) {
                ASSERT_ALWAYS(r->r[i] > r->r[i-1]);
            }
            fprintf(f, "%" PRIu32 "\n", r->r[i]);
        }
    }
    fclose(f);
}
/*}}}*/

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

void write_info_file()
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
                    final_file_info[i*nvslices+j]->s);
        }
    }
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

int main(int argc, char * argv[])
{
    int rc;
    argv++,argc--;
    for(;argc;argv++,argc--) {
        if (strcmp(argv[0], "--subdir") == 0) {
            argv++,argc--;
            if (!argc) usage();
            rc = chdir(argv[0]);
            if (rc < 0 && errno == ENOENT) {
                rc = mkdir(argv[0], 0777);
                DIE_ERRNO_DIAG(rc < 0, "mkdir", argv[0]);
                rc = chdir(argv[0]);
                DIE_ERRNO_DIAG(rc < 0, "chdir", argv[0]);
            }
        } else if (strcmp(argv[0], "--legacy") == 0) {
            legacy = 1;
        } else if (strcmp(argv[0], "--remove-input") == 0) {
            remove_input = 1;
        } else if (strcmp(argv[0], "--pad") == 0
                || strcmp(argv[0], "--square") == 0)
        {
            pad_to_square = 1;
        } else if (strcmp(argv[0], "--keep-temps") == 0) {
            keep_temps = 1;
        } else if (strcmp(argv[0], "--ram-limit") == 0) {
            argv++,argc--;
            if (!argc) usage();
            char * eptr;
            ram_limit = strtoul(argv[0], &eptr, 10);
            if (ram_limit == 0) usage();
            switch (*eptr) {
                case 'g': case 'G': ram_limit <<= 10;
                case 'm': case 'M': ram_limit <<= 10;
                case 'k': case 'K': ram_limit <<= 10;
                eptr++;
                default:
                if (*eptr != '\0') usage();
            }
        } else if (strcmp(argv[0], "--output-name") == 0
                || strcmp(argv[0], "--out") == 0)
        {
            argv++,argc--;
            if (!argc) usage();
            working_filename = argv[0];
        } else if (strcmp(argv[0], "--matrix") == 0
                || strcmp(argv[0], "--in") == 0)
        {
            argv++,argc--;
            if (!argc) usage();
            pristine_filename = argv[0];
        } else if (strcmp(argv[0], "--nbuckets") == 0
                || strcmp(argv[0], "--nslices") == 0
                || strcmp(argv[0], "--slices") == 0)
        {
            /* We accept a sole argument to mean n hslices and one
             * vslice, or an argument of the form MMMxNNN
             */
            argv++,argc--;
            if (!argc) usage();
            {
                char * ptr;
                nhslices = strtoul(argv[0], &ptr, 10);
                if (nhslices == 0)
                    usage();
                if (*ptr == '\0') {
                    nvslices = 1;
                } else if (*ptr != 'x') {
                    usage();
                } else {
                    ptr++;
                    nvslices = strtoul(ptr, &ptr, 10);
                    if (nvslices == 0)
                        usage();
                    if (*ptr != '\0')
                        usage();
                }
            }
        } else if (strcmp(argv[0], "--hslices") == 0) {
            argv++,argc--;
            if (!argc) usage();
            nhslices = atoi(argv[0]);
        } else if (strcmp(argv[0], "--vslices") == 0) {
            argv++,argc--;
            if (!argc) usage();
            nvslices = atoi(argv[0]);
        } else {
            usage();
        }
    }


    if (pristine_filename == NULL) {
        fprintf(stderr, "which matrix whall I read ?\n");
        exit(1);
    }

    /* 10 is overestimated here. I think 4 suffices: 0,1,2, and f.
     */
    if (sysconf(_SC_OPEN_MAX) <= nhslices + 10) {
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

    weight_sort_in_cells = 1;

    permute_rows = weight_sort_in_cells || nhslices > 1;
    permute_cols = weight_sort_in_cells || nvslices > 1;

    if (!(permute_rows || permute_cols || weight_sort_in_cells)) {
        /* Then we've got absolutely nothing to do ! */
        fprintf(stderr, "Nothing to do ? Really ? Why not use cp then ?\n");
        exit(1);
    }

    read_matrix(pristine_filename);

    fprintf(stderr, "Computing permutation\n");
    compute_permutation();

    if (row_table) free(row_table);
    if (col_table) free(col_table);

    if (permute_rows) {
        write_permutation("row_perm.txt",
                row_slices, nhslices, rows_are_weight_sorted);
    }
    if (permute_cols) {
        write_permutation("col_perm.txt",
                col_slices, nvslices, cols_are_weight_sorted);
    }

    fileset work;
    fileset_init(work, 1);
    work->names[0] = strdup(pristine_filename);
    work->status = INPUT;

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
        write_info_file();
    }

    fprintf(stderr, "Done\n");

    cleanup();
    return 0;
}
