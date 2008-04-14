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

int nslices = 1;

void usage()
{
    fprintf(stderr, "Usage: ./bw-balance <options>\n"
            "Options:\n"
            "--nslices <n>\toptimize for <n> horizontal strides\n"
            "--subdir <d>\tfetcch data in directory <d>\n"
            "--matrix <m>\tmatrix filename"
                " (defaults to matrix.txt when --subdir is set)\n"
           );
    exit(1);
}

struct row {
    uint32_t    i;
    int         w;
    off_t       o;
};

struct row * row_table;
uint32_t nr;

char * legacy_header_modulus;


int read_matrix(const char * filename)
{
    char buf[1024];

    FILE * f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        exit(1);
    }

    char * ptr;
    int siz;
    off_t o;

#define READ_LINE() do {						\
        o = ftello(f);                                                  \
        ptr = fgets(buf, sizeof(buf), f);				\
        if (ptr == NULL) {						\
            fprintf(stderr, "Unexpected %s at position %ld in %s\n",    \
                    feof(f) ? "EOF" : "error", o, filename);            \
            exit(1);							\
        }								\
        siz = strlen(ptr);						\
    } while (0)

#define GOBBLE_LONG_LINE() do {                                         \
    while (!(siz < (sizeof(buf)-1) || buf[sizeof(buf)-2] == '\n')) {	\
        ptr = fgets(buf, sizeof(buf), f);				\
        if (ptr == NULL) {						\
            fprintf(stderr, "Unexpected %s at position %ld in %s\n",    \
                    feof(f) ? "EOF" : "error", ftello(f), filename);    \
            exit(1);							\
        }								\
        siz = strlen(ptr);						\
    }                                                                   \
    } while (0)								\

    
    uint32_t nc;
    uint32_t i;

    /* Read matrix header */
    READ_LINE();
    legacy_header_modulus = NULL;
    if (strncmp(buf, "//", 2) == 0) {
        int pos;
        if (sscanf(buf, "// %" SCNu32 " ROWS %" SCNu32 "%n",
                    &nr, &nc, &pos) < 2)
        {
            fprintf(stderr, "Parse error while reading header: %s\n", buf);
            exit(1);
        }
        legacy_header_modulus = strdup(buf + pos);
    } else if (sscanf(buf, "%" SCNu32 " %" SCNu32, &nr, &nc) != 2) {
        fprintf(stderr, "Parse error while reading header: %s\n", buf);
        exit(1);
    }
    GOBBLE_LONG_LINE();

    if (nc < nr) {
        fprintf(stderr, "Uh ? Too many rows.\n"
                "Perhaps the matrix should have been transposed first\n");
        exit(1);
    }

    /* nc is the largest dimension. We're going to work with a matrix of
     * size nc * nc. If nr is specified as being smaller than this, then 
     * we'll pad with zeroes
     */

    row_table = malloc((nc+1) * sizeof(struct row));
    if (row_table == NULL) abort();

    fprintf(stderr, "Row table: %"PRIu32" entries, %lu MB\n",
            nc, (nc * sizeof(struct row)) >> 20);
    fprintf(stderr, "Reading row table...");
    fflush(stderr);

    for(i = 0 ; i < nr ; i++) {
        int w;
        READ_LINE();
        if (sscanf(buf, "%d", &w) != 1) {
            fprintf(stderr, "Parse error while reading line %"PRIu32" of %s\n",
                    i+1, filename);
            exit(1);
        }
        row_table[i].i = i;
        row_table[i].w = w;
        row_table[i].o = o;
        GOBBLE_LONG_LINE();
    }
    for( ; i < nc ; i++) {
        row_table[i].i = i;
        row_table[i].w = 0;
        row_table[i].o = o;
    }
    row_table[nc].o = ftello(f);

    fprintf(stderr, "done\n");
    fflush(stderr);

    fclose(f);

    nr = nc;

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


struct row ** slices;
uint32_t * slice_sizes;


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

struct bucket * heap;

int heap_index_compare(const struct bucket * a, const struct bucket * b)
{
    return (int32_t) a->i - (int32_t) b->i;
}

int heap_compare(const struct bucket * a, const struct bucket * b)
{
    /* Returns whether b is a batter candidate to be on top of the heap
     * IOW: Whether b is lighter than a
     * */
    if (b->room == 0) return 0;
    if (a->room == 0) return 1;
    return b->s < a->s;
}


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


void do_sorting()
{
    fprintf(stderr, "sorting...");
    fflush(stderr);

    qsort(row_table, nr, sizeof(struct row), (sortfunc_t) &decr_weight_cmp);

    fprintf(stderr, "done\n");
    fflush(stderr);

    slices = malloc(nslices * sizeof(struct row *));
    int i;

    uint32_t slice_min_size = nr / nslices;

    slice_sizes = malloc(nslices * sizeof(uint32_t));
    heap = malloc(nslices * sizeof(struct bucket));

    for(i = 0 ; i < nslices ; i++) {
        slice_sizes[i] = slice_min_size + (i < (nr % nslices));
        slices[i] = malloc((slice_sizes[i] + 1)* sizeof(struct row));
        /* This is a safeguard for what we're doing in shipout() */
        slices[i][slice_sizes[i]].i = UINT_MAX;
        heap[i].s = 0;
        heap[i].i = i;

        heap[i].room = slice_sizes[i];
    }

    make_heap(heap, heap + nslices);

    for(i = 0 ; i < nr ; i++) {
        int j = heap[0].i;
        struct row * tail = slices[j] + (slice_sizes[j] - heap[0].room);
        *tail = row_table[i];
        heap[0].room--;
        heap[0].s += row_table[i].w;
        pop_heap(heap, heap + nslices);
        push_heap(heap, heap + nslices);
    }

    qsort(heap, nslices, sizeof(struct bucket),
            (sortfunc_t) heap_index_compare);

    for(i = 0 ; i < nslices ; i++) {
        int j = heap[i].i;
        printf("Slice %d, %ld rows, weight %ld\n",
                i, slice_sizes[j] - heap[i].room,
                heap[i].s);
        if (heap[i].room != 0) {
            abort();
        }

        fprintf(stderr, "Sorting slice %d...", j);fflush(stderr);
        qsort(slices[j], slice_sizes[j], sizeof(struct row),
                (sortfunc_t) row_compare_index);
        fprintf(stderr, "done\n");fflush(stderr);
    }
}



void shipout(const char * filename)
{
    uint32_t i;
    int j;
    struct row ** ptr;
    char ** chunks;
    FILE ** g;
    FILE * f;
    char * tbuf = NULL;
    int tbufsiz = 0;

#define ENSURE(amount__)						\
        do {								\
            if (tbufsiz < (amount__)) {					\
                if (tbufsiz == 0) tbufsiz=1024;				\
                for( ; tbufsiz < (amount__) ; tbufsiz <<= 1);		\
                tbuf = realloc(tbuf, tbufsiz);				\
            }								\
        } while (0)

#define BIGCHUNK        16384
    ENSURE(BIGCHUNK);

    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        exit(1);
    }

    ptr = malloc(nslices * sizeof(struct row *));
    chunks = malloc(nslices * sizeof(char *));
    g = malloc(nslices * sizeof(FILE *));
    for(j = 0 ; j < nslices ; j++) {
        ptr[j] = slices[j];
        chunks[j] = malloc(strlen(filename) + 20);
        snprintf(chunks[j], strlen(filename) + 20, "%s.slice%d", filename, j);
        g[j] = fopen(chunks[j], "w");
        if (g[j] == NULL) {
            fprintf(stderr, "%s: %s\n", chunks[j], strerror(errno));
            exit(1);
        }
    }
    j = 0;
    off_t previous = 0;
    off_t next;
    for(i = 0 ; i < nr ; i++) {
        assert(i == 0 || ftello(f) == previous);
        /* pain. If I weren't tired enough of C++-to-C, I'd do some more
         * heap stuff here */
        for(j = 0 ; j < nslices && ptr[j]->i != i ; j++);
        if (j == nslices) {
            abort();
        }
        next = ptr[j]->o;
        ptr[j]++;
        if (i == 0) {
            previous = next;
            fseeko(f, previous, SEEK_SET);
            continue;
        }
        ENSURE(next - previous);
        if (next == previous) {
            /* special provision for writing zero rows */
            fprintf(g[j], "0\n");
            continue;
        }
        fread(tbuf, next - previous, 1, f);
        fwrite(tbuf, next - previous, 1, g[j]);
        previous = next;
    }

    /* Do not forget the last line */
    next = row_table[nr].o;
    ENSURE(next - previous);
    if (next == previous) { fprintf(g[j], "0\n"); }
    fread(tbuf, next - previous, 1, f);
    fwrite(tbuf, next - previous, 1, g[j]);

    for(j = 0 ; j < nslices ; j++) {
        fclose(g[j]);
    }
    fclose(f);

    /* Now, the final stuff. We rebuild a complete matrix from the
     * slices. It should be avoided, since we have already done something
     * which will be useful afterwards.
     */
    char * renamed = malloc(strlen(filename) + 20);
    snprintf(renamed, strlen(filename) + 20, "%s.old", filename);
    char * newfile = malloc(strlen(filename) + 20);
    snprintf(newfile, strlen(filename) + 20, "%s.new", filename);


    f = fopen(newfile, "w");
    if (f == NULL) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        exit(1);
    }
    if (legacy_header_modulus) {
        fprintf(f, "// %" PRIu32 " ROWS %" PRIu32"%s",
                nr, nr, legacy_header_modulus);
    } else {
        fprintf(f, "%" PRIu32 " %" PRIu32 "\n", nr, nr);
    }

    for(j = 0 ; j < nslices ; j++) {
        g[j] = fopen(chunks[j], "r");
        if (g[j] == NULL) {
            fprintf(stderr, "%s: %s\n", chunks[j], strerror(errno));
            exit(1);
        }
        for(;;) {
            size_t nread;
            size_t nwritten;
            nread = fread(tbuf, 1, BIGCHUNK, g[j]);
            if (nread == 0) break;
            nwritten = fwrite(tbuf, 1, nread, f);
            if (nread != nwritten) {
                fprintf(stderr, "Copying %s into %s: %s\n",
                        chunks[j], newfile, strerror(errno));
                exit(1);
            }
        }
        fclose(g[j]);
    }
    fclose(f);

    if (rename(filename, renamed) < 0) {
        fprintf(stderr, "rename(%s,%s) : %s\n", filename, renamed,
                strerror(errno));
    }
    if (rename(newfile, filename) < 0) {
        fprintf(stderr, "rename(%s,%s) : %s\n", newfile, filename,
                strerror(errno));
    }

    /* FIXME: Really a pity to do that. */
    for(j = 0 ; j < nslices ; j++) {
        unlink(chunks[j]);
    }
    /* XXX We still keep the old file, because pristine source are
     * valuable ; it's up to the calling script to remove it if desired */



    free(renamed);
    free(newfile);
    for(j = 0 ; j < nslices ; j++) {
        free(chunks[j]);
        /* ptr[j] is really only a pointer */
    }
    free(g);
    free(ptr);
    free(chunks);
    free(tbuf);
}

void cleanup()
{
    int i;
    free(slice_sizes);
    for(i = 0 ; i < nslices ; i++) {
        free(slices[i]);
    }
    free(slices);
    free(heap);
    free(row_table);
    free(legacy_header_modulus);
}

int main(int argc, char * argv[])
{
    const char * matrix_filename = NULL;
    argv++,argc--;
    for(;argc;argv++,argc--) {
        if (strcmp(argv[0], "--subdir") == 0) {
            argv++,argc--;
            if (!argc) usage();
            if (chdir(argv[0]) < 0) {
                fprintf(stderr, "chdir(%s): %s\n", argv[0], strerror(errno));
                exit(1);
            }
            if (matrix_filename == NULL) {
                /* provide this default only when --subdir is given */
                matrix_filename = "matrix.txt";
            }
        } else if (strcmp(argv[0], "--matrix") == 0) {
            argv++,argc--;
            if (!argc) usage();
            matrix_filename = argv[0];
        } else if (strcmp(argv[0], "--nbuckets") == 0 ||
                strcmp(argv[0], "--nslices")) {
            argv++,argc--;
            if (!argc) usage();
            nslices = atoi(argv[0]);
        } else {
            usage();
        }
    }

    if (matrix_filename == NULL) {
        fprintf(stderr, "which matrix whall I read ?\n");
        exit(1);
    }

    read_matrix(matrix_filename);
    do_sorting();
    shipout(matrix_filename);
    cleanup();
    return 0;
}
