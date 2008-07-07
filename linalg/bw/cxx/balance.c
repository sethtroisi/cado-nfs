#define __STDC_FORMAT_MACROS    /* For C++ to have PRId64 and friends */
#define _GNU_SOURCE    /* for strndup only -- really easy to make a workalike */

#include <stddef.h>     /* ptrdiff_t */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <inttypes.h>   /* SCNu32 etc */
#include <assert.h>
#include <limits.h>
#include <ctype.h>


#include <math.h>       /* XXX probably best to remove it */

#ifdef  __cplusplus
using namespace std;
#endif

/* I'd prefer the program to stay completely standalone */
#ifndef ASSERT_ALWAYS
#define ASSERT_ALWAYS(x) do { if (!(x)) abort(); } while (0)
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


unsigned int nhslices = 1;
unsigned int nvslices = 1;

int keep_temps = 0;
int remove_input = 0;

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
const char * working_filename = "mat";

void usage()
{
    fprintf(stderr, "Usage: ./bw-balance <options>\n"
            "Options:\n"
            "--hslices <n>\toptimize for <n> horizontal strides\n"
            "--subdir <d>\tfetch data in directory <d>\n"
            "--matrix <m>\tmatrix filename\n"
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

// {{{ give mean, std dev and so on
void give_stats(const char * text, const struct row * data, unsigned int n)
{
    uint32_t dmin = UINT_MAX;
    uint32_t dmax = 0;
    double sumcf = 0;
    double sumcf2 = 0;
    uint32_t k;
    for(k = 0 ; k < n ; k++) {
        uint32_t cf = data[k].w;
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
    /*
       for(k = 0 ; k < nc ; k++) {
       printf(" %"PRIu32, col_freq[k]);
       }
       printf(" \n");
       */
}
// }}}

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

    /* nmax is the largest dimension. We're going to work with a matrix of
     * size nmax * nmax. If nr is specified as being smaller than this, then 
     * we'll pad with zeroes
     */
    uint32_t nmax = nc < nr ? nr : nc;

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
        int k;
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


    nr = nmax;
    nc = nmax;

    return 0;
}

int decr_weight_cmp(const struct row * a, const struct row * b)
{
    return b->w - a->w;
}

int row_compare_index(const struct row * a, const struct row * b)
{
    return a->i - b->i;
}

typedef int (*sortfunc_t)(const void *, const void *);



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

inline void
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

inline void
pop_heap(struct bucket * __first, struct bucket * __last)
{
    struct bucket ref = __last[-1];
    __last[-1] = __first[0];
    adjust_heap(__first, 0, __last - __first, & ref);
}

inline void
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

/* {{{ Datatype containing description for slices of the matrix */
struct slice {
    struct row * data;
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
        free(slices[i].data);
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
        res[i].data = (struct row *) malloc((res[i].nrows + 1)
                    * sizeof(struct row));
        /* This is a safeguard for what we're doing in shipout() */
        // XXX still used ?
        res[i].data[res[i].nrows].i = UINT_MAX;
    }

    return res;
}
/* }}} */

/* This is the basic procedure which dispatches rows (or columns) in
 * buckets according to our preferred strategy
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
        int k;
        /* If it turns out that we have a great many rows of the same
         * weight (this sort of thing happens), then we'll dispatch them
         * in large chunks and not in a completely alternating way */
        for(ni = i + 1 ; ni < n && rt[i].w == rt[ni].w ; ni++);
        chunk = (ni-i) / ns;
        if (chunk) {
            int oversize=0;
            for(k = 0 ; k < ns ; k++) {
                int l;
                for(l = 0 ; l < chunk && heap[k].room ; l++) {
                    int j = heap[k].i;
                    int pos = slices[j].nrows-heap[k].room;
                    slices[j].data[pos]=rt[i];
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
            slices[j].data[pos] = rt[i];
            heap[0].s += rt[i].w;
            heap[0].room--;
            pop_heap(heap, heap + ns);
            push_heap(heap, heap + ns);
        }
    }

    qsort(heap, ns, sizeof(struct bucket), (sortfunc_t) heap_index_compare);

    for(i = 0 ; i < ns ; i++) {
        int j = heap[i].i;
        assert(heap[i].i == i);
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

struct slice * row_slices;
struct slice * col_slices;

void do_sorting()
{
    if (permute_rows) {
        row_slices = shuffle_rtable("horizontal", row_table, nr, nhslices);
    }
    if (permute_cols) {
        col_slices = shuffle_rtable("vertical", col_table, nc, nvslices);
    }
}

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

void dispatcher(
        const char * dst, const char * key, unsigned int n,
        const char * src,
        size_t hdr, size_t * sizes,
        unsigned int *dispatch,
        unsigned int nr)
{
    char ** chunks;
    FILE ** g;
    unsigned int ii;
    FILE * f;
    int i;

    copybuf cb;

    cb_init(cb);
#define BIGCHUNK        16384
    cb_ensure(cb, BIGCHUNK);

    f = fopen(src, "r");
    DIE_ERRNO_DIAG(f == NULL, "fopen", src);

    open_hslices_files(src, dst, key, &chunks, &g, n, "w");

    if (hdr) {
        fseek(f, hdr, SEEK_SET);
    }

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


    close_hslices_files(&chunks, &g, n);
}


/* We're only doing this if nhslices > 1 or if we are asked to size-sort
 * rows
 */
void dispatch_row_slices()
{
    uint32_t i;
    unsigned int * dispatch;

    dispatch = (unsigned int *) malloc(nr * sizeof(unsigned int));
    for(i = 0 ; i < nhslices ; i++) {
        uint32_t k;
        for(k = 0 ; k < row_slices[i].nrows ; k++) {
            dispatch[row_slices[i].data[k].i]=i;
        }
    }

    dispatcher(working_filename, "h", nhslices,
            pristine_filename, header_bytes, line_bytes,
            dispatch, nr);

    free(dispatch);
}

void write_permutation(const char * filename,
        struct slice * slice, int already_sorted)
{
    FILE * f;
    uint32_t i,ii;
    f = fopen(filename, "w");
    i = 0;

    for(ii = 0 ; ii < nhslices ; ii++) {
        const struct slice * r = &(slice[ii]);
        for(i = 0 ; i < r->nrows ; i++) {
            if (already_sorted && i) {
                ASSERT_ALWAYS(r->data[i].i > r->data[i-1].i);
            }
            fprintf(f, "%" PRIu32 "\n", r->data[i].i);
        }
    }
    fclose(f);
}

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


void do_per_hslice_stuff()
{
    char ** chunks;
    uint32_t i, ii;
    copybuf cb;
    uint32_t * cdirect;
    cperm_buf cp;

    if (permute_cols) {
        cperm_init(cp, nvslices);
        cdirect = (uint32_t *) malloc(nc * sizeof(uint32_t));
        memset(cdirect, 0, nc * sizeof(uint32_t));
        for(ii = 0 ; ii < nvslices ; ii++) {
            const struct slice * r = &(col_slices[ii]);
            for(i = 0 ; i < r->nrows ; i++) {
                cdirect[r->data[i].i]=r->i0+i;
            }
        }
    }

    cb_init(cb);

    open_hslices_files(pristine_filename, working_filename, "h", 
            &chunks, NULL, nhslices, "r");

    for(ii = 0 ; ii < nhslices ; ii++) {
        if (weight_sort_in_cells && !rows_are_weight_sorted) {
            const char * writable_name;
            if (nhslices == 1) {
                writable_name=working_filename;
            } else {
                writable_name=chunks[ii];
            }
            cb_ensure(cb, BIGCHUNK);
            /* We have to do the sorting */
            uint32_t nrs = row_slices[ii].nrows;
            /* At this moment, we know that rows indexed within
             * row_slices[ii].data[] are effectively in chunks[ii].
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
                wherefrom[i].y=row_slices[ii].data[i].i;
            }
            qsort(wherefrom, nrs, sizeof(struct idxpair), (sortfunc_t) idxcmp2);
            /* now wherefrom[kk].x is the index in the current file of the row
             * that will go to position [kk] eventually.
             * wherefrom[kk].y is the index of that same row in the
             * original big file */

            /* So we need to decide who goes where. Of course we could
             * invert again the wherefrom array, but we'll pay too much
             * for this, and anyway this array will become cumbersome. So
             * let's rebuild the data in a more tidy way.
             *
             * We store the inverse permutation in shuffle[]
             */
            uint32_t * shuffle = malloc(nrs * sizeof(uint32_t));
            for(i = 0 ; i < nrs ; i++) {
                shuffle[wherefrom[i].x]=i;
            }
            free(wherefrom); wherefrom = NULL;
            /* We now know that row shuffle[j] in the file chunks[ii] is in fact
             * the one whose real index is
             * row_slices[ii].data[j].i
             *
             * our purpose is to make it go to position j in
             * that same file
             */


            /* To use the dispatcher() routine, we must create a
             * dispatch[] array, and an nbytes[] one.
             * Note that while both are available, one could do without
             * if needed:
             * nbytes[] could be obtained on the fly
             * dispatch[] is really so trivial (mere index comparison)...
             */

            /* nbytes[] gives the size in bytes of the rows found in the
             * chunk
             */
            size_t * nbytes = malloc(nrs * sizeof(size_t));

            uint32_t tot = 0;
            for(i = 0 ; i < nrs ; i++) {
                nbytes[shuffle[i]]=line_bytes[row_slices[ii].data[i].i];
            }

            /* Split the file into RAM-sized chunks */
            for(i = 0 ; i < nrs ; i++) {
                tot += nbytes[i];
            }
            /* Decide exactly which will be the chunk size ; otherwise
             * we'll do some ridiculous cuts */
            unsigned int nsubchunks = ((tot+ram_limit-1)/ram_limit);
            unsigned int nsc_max = (2*nsubchunks+1); // allow round-off err.
            struct chunk_info_s {
                unsigned int i0;
                unsigned int i1;
                size_t sz;
            };
            typedef struct chunk_info_s chunk_info[1];
            chunk_info * sc_inf = (chunk_info *) malloc(nsc_max*sizeof(chunk_info));

            unsigned int * dispatch = malloc(nrs * sizeof(unsigned int));
            /* We're not going to split over many files if there's only
             * one such file */
            if (nsubchunks > 1) {
                /* The 1024 here should be sufficient to avoid round-off
                 * misses. They would not be catastrophic anyway. It's
                 * just a pain to have nsubchunks+1 cuts instead of
                 * nsubchunks (especially when the last one really is so
                 * tiny). */
                size_t my_ram_limit = 1024 + (tot+nsubchunks-1) / nsubchunks;
                unsigned int ram_chunk=0;
                sc_inf[0]->i0 = 0;
                tot=0;
                for(i = 0 ; i < nrs ; i++) {
                    size_t d = nbytes[shuffle[i]];
                    if (tot + d >= my_ram_limit) {
                        sc_inf[ram_chunk]->i1 = i;
                        sc_inf[ram_chunk]->sz = tot;
                        ram_chunk++;
                        ASSERT_ALWAYS(ram_chunk < nsc_max);
                        sc_inf[ram_chunk]->i0 = i;
                        tot = 0;
                    }
                    tot += d;
                    dispatch[shuffle[i]] = ram_chunk;
                }
                sc_inf[ram_chunk]->i1 = i;
                sc_inf[ram_chunk]->sz = tot;
                nsubchunks = ram_chunk + 1;
                printf("Splitting %s into %d intermediary files for sorting\n",
                        chunks[ii], nsubchunks);
                dispatcher(writable_name, "sort", nsubchunks,
                        chunks[ii], nhslices == 1 ? header_bytes : 0,
                        nbytes,
                        dispatch, nrs);

                if (remove_input && nhslices == 1) {
                    unlink(pristine_filename);
                }
            } else {
                printf("Sorting %s directly into memory\n", chunks[ii]);
                sc_inf[0]->i0 = 0;
                sc_inf[0]->sz = tot;
                sc_inf[0]->i1 = nrs;
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

            off_t * poking_place;
            poking_place = (off_t *) malloc(nrs * sizeof(off_t));

            int * tmp;
            tmp = (int *) malloc(nsubchunks*sizeof(int));
            unsigned int s;
            for(s = 0 ; s < nsubchunks ; s++) tmp[s] = sc_inf[s]->i0;
            tot=0;
            /* Hard to avoid for the computation to come */
            int * direct = (int*) malloc(nrs * sizeof(int));
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

#if 0
            for(i = 0 ; i < 100 ; i++) {
                printf("%u(%lu) %u\n", direct[i], nbytes[i], dispatch[i]);
            }
#endif
            for(s = 0 ; s < nsubchunks ; s++) tmp[s] = sc_inf[s]->i0;
            for(i = 0 ; i < nrs ; i++) {
                int which = nsubchunks > 1 ? dispatch[i] : 0;
                /* This row is in sub chunk-[[which]] */
                /* The next unassigned byte-size from that sub-chunk is
                 * nbytes[tmp[which]].
                 */
                size_t sz = line_bytes[row_slices[ii].data[direct[i]].i];
                /* Note that it is
                 * important that we process this array in direct order.
                 */
                nbytes[tmp[which]] = sz;
                poking_place[tmp[which]]=
                    poking_in_ordered_file[direct[i]] -
                    poking_in_ordered_file[sc_inf[which]->i0];
                tmp[which]++;
            }
            char ** subchunks;
            open_hslices_files(chunks[ii], writable_name, "sort",
                    &subchunks, NULL, nsubchunks, "r");

            /* Free as much as we can before doing the large-mem stuff */
            free(dispatch); dispatch=NULL;
            free(tmp); tmp=NULL;
            free(direct); direct = NULL;
            free(poking_in_ordered_file); poking_in_ordered_file = NULL;
            free(shuffle); shuffle = NULL;

            FILE * f;
            FILE ** fv;
            char ** vchunks;
            
            if (nvslices > 1) {
                open_hslices_files(chunks[ii], writable_name, "v",
                        &vchunks, &fv, nvslices, "w");
                f = NULL;
            } else if (nhslices > 1) {
                /* Only one vslice, but several hslices. Each hslice ends
                 * up in a different file. */
                /* writable_name is chunks[ii] in this case */
                assert(writable_name == chunks[ii]);
                f = fopen(writable_name, "w");
                DIE_ERRNO_DIAG(f == NULL, "fopen", writable_name);
            } else {
                /* Well, we're only doing sorting here. So we'll create a
                 * new copy of the matrix.  Notice that chunks[ii] is
                 * actually the source matrix. We should not overwerite
                 * it ! */
                /* writable_name is working_filename in this case */
                assert(writable_name == working_filename);
                assert(strcmp(chunks[ii], pristine_filename) == 0);
                assert(permute_rows || permute_cols);
                f = fopen(writable_name, "w");
                DIE_ERRNO_DIAG(f == NULL, "fopen", writable_name);
                /* We need to recover the header !!! */
                fwrite(header, 1, header_bytes, f);
            }


            for(s = 0 ; s < nsubchunks ; s++) {
                FILE * g = fopen(subchunks[s],"r");
                DIE_ERRNO_DIAG(g == NULL, "fopen", subchunks[s]);
                struct chunk_info_s * sc = sc_inf[s];
                char * inmem = malloc(sc->sz);
                if (nhslices == 1 && nsubchunks == 1) {
                    fseek(g, header_bytes, SEEK_SET);
                }
                for(i = sc->i0 ; i < sc->i1 ; i++) {
                    off_t sz = nbytes[i];
                    ASSERT_ALWAYS(sz != 0);
                    cb_ensure(cb, sz);
                    fread(cb->buf, 1, sz, g);
                    assert(cb->buf[sz-1] == '\n');
                    assert(poking_place[i] >= 0);
                    assert(poking_place[i] + sz <= sc->sz);
                    memcpy(inmem + poking_place[i], cb->buf, sz);
                }
                fclose(g);
                if (remove_input && nhslices == 1 && nsubchunks == 1) {
                    /* last corner case in which we might have actually
                     * needed the source file until here. */
                    assert(strcmp(subchunks[s], pristine_filename) == 0);
                    unlink(pristine_filename);
                }
                /* Here we go ; copy the file in correct order now */
                if (permute_cols) {
#if 0
                    {
                        char * x;
                        asprintf(&x, "%s.debug%d", chunks[ii], s);
                        f =fopen(x,"w");
                        fwrite(inmem, 1, sc->sz, f);
                        fclose(f);
                    }

    /* Just for fun, here's the (zsh) shell loop that I used to check
     * that this program is not completely nuts
cat /tmp/mat.h0.debug* | (n=0; while read x ; do echo "row $n" >&2 ; echo $x | (read n a ; for x in ${=a} ; do grep -n '^'$x'$' /tmp/col_perm.txt | cut -d: -f1 | (read z ; expr $z - 1); done) | sort >! /tmp/orig-$n ; for i in {5..8} ; do read a b <&$i ; echo $b | xargs -n 1 echo ; done | sort >! /tmp/sorted-$n ; diff -u /tmp/orig-$n /tmp/sorted-$n ;  let n+=1; done 5</tmp/mat.h0.v0 6</tmp/mat.h0.v1 7</tmp/mat.h0.v2 8</tmp/mat.h0.v3)
     * (well, actually, _who_ is nuts ?)
     */

#endif

                    /* We have to parse now :-( */
                    off_t pos = 0;
                    for(i = sc->i0 ; i < sc->i1 ; i++) {
                        /* FIXME */
                        size_t sz = line_bytes[row_slices[ii].data[i].i];
                        unsigned int nj;
                        unsigned int jj;
                        const char * data = inmem+pos;
                        const char * data0 = data;
                        char * ndata;
                        cperm_reset(cp);
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
                            w = which_slice(col_slices, nvslices, cdirect[j]);
                            snprintf(scrap, sizeof(scrap), " %u",
                                    cdirect[j]-col_slices[w].i0);
                            cperm_append(cp, w, scrap, strlen(scrap));
                            if (*data == ':') {
                                int e;
                                data++;
                                e=strtol(data, &ndata, 10);
                                PARSE_CHECK(data==ndata, "row",
                                        strndup(data0,sz));
                                cperm_append(cp, w, data-1, ndata-data+1);
                                data = ndata;
                                /* OK, we're at the innermost block, 8
                                 * levels below earth. Shame on me. */
                            }
                            cp->w[w]++;
                        }
                        assert(*data == '\n');
                        data=NULL;
                        pos+=sz;
                        if (nvslices > 1) {
                            for(jj = 0 ; jj < nvslices ; jj++) {
                                fprintf(fv[jj], "%u", cp->w[jj]);
                                fwrite(cp->bufs[jj]->buf, 1, cp->pos[jj], fv[jj]);
                                fputc('\n', fv[jj]);
                            }
                        } else {
                            fprintf(f, "%u", nj);
                            fwrite(cp->bufs[0]->buf, 1, cp->pos[0], f);
                            fputc('\n', f);
                        }
                    }
                } else {
                    fwrite(inmem, 1, sc->sz, f);
                }
                free(inmem);
                if (!keep_temps && nsubchunks > 1) unlink(subchunks[s]);
            }
            if (nvslices > 1) {
                close_hslices_files(&vchunks, &fv, nvslices);
                if (!keep_temps) unlink(writable_name);
            } else {
                fclose(f);
            }

            free(sc_inf); sc_inf = NULL;
            free(poking_place); poking_place = NULL;
            close_hslices_files(&subchunks, NULL, nsubchunks);
            
            free(nbytes);
        }
    }

    close_hslices_files(&chunks, NULL, nhslices);

    cb_clear(cb);

    if (permute_cols) {
        cperm_clear(cp);
        free(cdirect);
    }
}

void cleanup()
{
    /* free() only based on pointer, because we might have put off the
     * flag midway in case we noticed that the input was sorted already
     */
    if (row_slices)     free_slices(row_slices, nhslices);
    if (row_table)      free(row_table);
    if (col_slices)     free_slices(col_slices, nvslices);
    if (col_table)      free(col_table);

    free(line_bytes);
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
            DIE_ERRNO_DIAG(rc < 0, "chdir", argv[0]);
        } else if (strcmp(argv[0], "--remove-input") == 0) {
            argv++,argc--;
            remove_input = 1;
        } else if (strcmp(argv[0], "--keep-temps") == 0) {
            argv++,argc--;
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
        } else if (strcmp(argv[0], "--output-name") == 0) {
            argv++,argc--;
            if (!argc) usage();
            working_filename = argv[0];
        } else if (strcmp(argv[0], "--matrix") == 0) {
            argv++,argc--;
            if (!argc) usage();
            pristine_filename = argv[0];
        } else if (strcmp(argv[0], "--nbuckets") == 0 ||
                strcmp(argv[0], "--nslices") == 0 ||
                strcmp(argv[0], "--slices") == 0) {
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
    do_sorting();

    if (permute_rows) {
        write_permutation("row_perm.txt", row_slices, rows_are_weight_sorted);
    }
    if (permute_cols) {
        write_permutation("col_perm.txt", col_slices, cols_are_weight_sorted);
    }

    if (nhslices > 1) {
        fprintf(stderr, "Dispatching matrix in row slices\n");
        dispatch_row_slices();
        if (remove_input) {
            unlink(pristine_filename);
        }
    }

    /* Now we're going to handle all the horizontal slices in turn. We'll
     * eventually take this opportunity to do actual sorting by weight if
     * required. But prior to that, we'll also split slices vertically
     * first. */
    do_per_hslice_stuff();

    fprintf(stderr, "Done\n");

    cleanup();
    return 0;
}
#if 0
    if (header) {
        char * newfile = malloc(strlen(dst) + 20);
        snprintf(newfile, strlen(dst) + 20, "%s.%sslices", dst);

        f = fopen(newfile, "w");
        DIE_ERRNO_DIAG(f == NULL, "fopen", newfile);
        fwrite(header, 1, header_bytes, f);
        fprintf(f, "%d\n", nhslices);

        i = 0;
        for(ii = 0 ; ii < nhslices ; ii++) {
            /* Want I J i0 ncoeffs filename */
            fprintf(f, "%u %u %"PRIu32" %"PRIu32 " %" PRIu32 " %s\n",
                    ii, 0,
                    i,
                    i+row_slices[ii].nrows, row_slices[ii].coeffs, chunks[ii]);
            i += row_slices[ii].nrows;
        }
        fclose(f);
        free(newfile);
        free(header);
    }
#endif
