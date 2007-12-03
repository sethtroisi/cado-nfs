#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>

#include "readmat.h"
#include "manu.h"

/* Number of unsigned ints that fit in RAM */
#define	MEM_LIMIT	(1UL << 20)	

/* Rows processed simultaneously when outputting the result */
#define NROWS_CORE      20

unsigned int nfiles = 0;

const char * tmpdir;
const char * filename_in;
const char * filename_out;

struct vec_uint_struct {
    unsigned int * data;
    unsigned int size;
    unsigned int alloc;
};

typedef struct vec_uint_struct vec_uint_t[1];

sparse_mat_t mat;
vec_uint_t * transposed;

unsigned int vec_next_size(unsigned int x)
{
    return MAX(x * 2, x + 16);
}

void vec_uint_init(vec_uint_t x)
{
    memset(x, 0, sizeof(vec_uint_t));
}
void vec_uint_clear(vec_uint_t x)
{
    free(x->data);
    memset(x, 0, sizeof(vec_uint_t));
}
void vec_uint_init_array(vec_uint_t ** x, unsigned int n)
{
    *x = (vec_uint_t *) malloc(n * sizeof(vec_uint_t));
    memset(*x, 0, n * sizeof(vec_uint_t));
}
void vec_uint_clear_array(vec_uint_t ** x, unsigned int n)
{
    unsigned int j;
    for(j = 0 ; j < n ; j++) {
        free((*x)[j]->data);
    }
    memset((*x), 0, n * sizeof(vec_uint_t));
}
void vec_uint_push_back(vec_uint_t ptr, unsigned int x)
{
    if (ptr->size >= ptr->alloc) {
        ptr->alloc = vec_next_size(ptr->alloc);
        ptr->data = realloc(ptr->data, ptr->alloc * sizeof(unsigned int));
        BUG_ON(ptr->data == NULL);
    }
    ptr->data[ptr->size++] = x;
}

void transpose_and_save(unsigned int i0, unsigned int i1)
{
    char dstfile[80];

    const unsigned int * src = mat->data;
    unsigned int i, j, c;

    for(i = i0 ; i < i1 ; i++) {
        unsigned int nc = * src++;
        for(j = 0 ; j < nc ; j++) {
            c = *src++;
            ASSERT(c < mat->ncols);
            vec_uint_push_back(transposed[c], i);
            ASSERT(transposed[c]->size < mat->nrows);
        }
    }

    snprintf(dstfile, sizeof(dstfile), "%s/temp.%u.%04d",
            tmpdir, getpid(), nfiles);

    FILE * out;
    out = fopen(dstfile, "w"); BUG_ON_MSG(!out, "cannot open temp file");
    for(j = 0 ; j < mat->ncols ; j++) {
        fprintf(out, "%d", transposed[j]->size);
        for(c = 0 ; c < transposed[j]->size ; c++) {
            fprintf(out, " %d", transposed[j]->data[c]);
        }
        fprintf(out, "\n");
        // reset the size pointer for next chunk
        transposed[j]->size = 0;
    }
    fclose(out);

    /*
    fprintf(stderr, "Written %s (rows %u to %u, transposed)\n",
            dstfile, i0, i1);
            */
}


void rebuild()
{
    FILE * fout;
    FILE ** fin_vchunks;
    unsigned int i, j;
    int ret;

    fout = fopen(filename_out, "w");
    BUG_ON_MSG(fout == NULL, "Cannot open output matrix");

    fprintf(fout, "%u %u\n", mat->ncols, mat->nrows);

    fin_vchunks = malloc(nfiles * sizeof(FILE *));
    for(i = 0 ; i < nfiles ; i++) {
        char dstfile[80];
        snprintf(dstfile, sizeof(dstfile), "%s/temp.%u.%04d",
                tmpdir, getpid(), i);
        fin_vchunks[i] = fopen(dstfile, "r");
        BUG_ON_MSG(fin_vchunks[i] == NULL, "Cannot reopen temp file");
    }

    vec_uint_t * rows;
    vec_uint_init_array(&rows, NROWS_CORE);

    for(j = 0 ; j < mat->ncols ; j+= NROWS_CORE) {
        unsigned int k;
        unsigned int kmax = MIN(NROWS_CORE, mat->ncols - j);
        for(k = 0 ; k < NROWS_CORE ; k++) {
            rows[k]->size = 0;
        }
        for(i = 0 ; i < nfiles ; i++) {
            unsigned int nc;
            unsigned int c;
            for(k = 0 ; k < kmax ; k++) {
                ret = fscanf(fin_vchunks[i], "%u", &nc);
                BUG_ON(ret != 1);
                ret = 0;
                for(c = 0 ; c < nc ; c++) {
                    unsigned int idx;
                    ret += fscanf(fin_vchunks[i], "%u", &idx);
                    vec_uint_push_back(rows[k], idx);
                }
                BUG_ON(ret != nc);
            }
        }
        for(k = 0 ; k < kmax ; k++) {
            unsigned int c;
            fprintf(fout, "%d", rows[k]->size);
            for(c = 0 ; c < rows[k]->size ; c++) {
                fprintf(fout, " %d", rows[k]->data[c]);
            }
            fprintf(fout, "\n");
        }
    }
    vec_uint_clear_array(&rows, NROWS_CORE);
    for(i = 0 ; i < nfiles ; i++) {
        fclose(fin_vchunks[i]);
    }
    free(fin_vchunks);
}


int main(int argc, char * argv[])
{

    FILE * fin;

    if (argc != 3 && argc != 4) {
        fprintf(stderr, "usage: ./transpose <input matrix> <output matrix> [ <tmp dir> ]\n");
    }

    if (argc >= 3) {
        filename_in = argv[1];
        filename_out = argv[2];
    }
    if (argc == 4) {
        tmpdir = argv[3];
    } else {
        tmpdir = "/tmp";
    }

    fin = fopen(filename_in, "r");
    BUG_ON_MSG(fin == NULL, "Cannot open input matrix");

    sparse_mat_init(mat);

    read_matrix_header(fin, mat);

    unsigned int i;
    unsigned int i0;
    unsigned int * dst = mat->data;
    unsigned int wt0 = 0;

    vec_uint_init_array(&transposed, mat->ncols);

    for(i = i0 = 0 ; i < mat->nrows ; ) {
        dst = read_matrix_row(fin, mat, dst, 1);
        i++;

        if (mat->wt >  wt0 + MEM_LIMIT) {
            fprintf(stderr, "Read %u rows and %u coeffs into core,"
                    " saving to tempfile\n", i, mat->wt - wt0);
            transpose_and_save(i0, i);
            i0 = i;
            nfiles++;
            dst = mat->data;
            wt0 = mat->wt;
        }
    }
    if (i != i0) {
        transpose_and_save(i0, i);
        i0 = i;
        nfiles++;
        dst = mat->data;
    }
    vec_uint_clear_array(&transposed, mat->ncols);

    rebuild();

    return 0;
}
/* vim: set sw=4 sta et: */
