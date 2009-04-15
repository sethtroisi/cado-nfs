#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>

#include "readmat.h"
#include "readmat-easy.h"

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

uint32_t * read_easy(const char * filename,
        unsigned int * p_nr, unsigned int * p_nc)
{
    sparse_mat_t smat;
    FILE * f;

    sparse_mat_init(smat);

    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", filename, strerror(errno));
        exit(1);
    }
    read_matrix_header(f, smat);

    uint32_t * data;
    size_t alloc = 0;
    size_t size = 0;

    data = NULL;

    if (p_nr) *p_nr = smat->nrows;
    if (p_nc) *p_nc = smat->ncols;

    /* do it the bovine way. read just about everything in memory. */
    for(unsigned int i = 0 ; i < smat->nrows ; i++) {
        read_matrix_row(f,smat,smat->data,1);
        CHECK_REALLOC(data, alloc, size, 1 + smat->data[0]);
        data[size++] = smat->data[0];
        qsort(smat->data + 1, smat->data[0],
                sizeof(smat->data[0]), (sortfunc_t) uint_cmp);
        for(unsigned int j = 0 ; j < smat->data[0] ; j++) {
            data[size++] = smat->data[1+j];
        }
    }

    fclose(f);

    return data;
}

uint32_t * read_easy_transposed(const char * filename,
        unsigned int * p_nr, unsigned int * p_nc)
{
    sparse_mat_t smat;
    FILE * f;

    sparse_mat_init(smat);

    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", filename, strerror(errno));
        exit(1);
    }
    read_matrix_header(f, smat);

    uint32_t * data;
    size_t alloc = 0;
    size_t size = 0;

    uint32_t * colweights;
    uint32_t ** colptrs;
    colweights = malloc(smat->ncols * sizeof(uint32_t));
    memset(colweights, 0, smat->ncols * sizeof(uint32_t));
    colptrs = malloc(smat->ncols * sizeof(uint32_t *));

    data = NULL;

    if (p_nr) *p_nr = smat->ncols;
    if (p_nc) *p_nc = smat->nrows;

    /* do it the bovine way. read just about everything in memory. */
    for(unsigned int i = 0 ; i < smat->nrows ; i++) {
        read_matrix_row(f,smat,smat->data,1);
        CHECK_REALLOC(data, alloc, size, 1 + smat->data[0]);
        data[size++] = smat->data[0];
        qsort(smat->data + 1, smat->data[0],
                sizeof(smat->data[0]), (sortfunc_t) uint_cmp);
        for(unsigned int j = 0 ; j < smat->data[0] ; j++) {
            data[size++] = smat->data[1+j];
            colweights[smat->data[1+j]]++;
        }
    }

    fclose(f);

    uint32_t * res = malloc(size * sizeof(uint32_t));
    uint32_t * ptr = res;
    for(unsigned int j = 0 ; j < smat->ncols ; j++) {
        uint32_t w = colweights[j];
        *ptr++ = w;
        colptrs[j] = ptr;
        ptr += w;
    }

    uint32_t * v = data;
    for(unsigned int i = 0 ; i < smat->nrows ; i++) {
        uint32_t rlen = *v++;
        for(unsigned int j = 0 ; j < rlen ; j++) {
            uint32_t c = *v++;
            *(colptrs[c]++) = i;
        }
    }

    free(data);
    free(colweights);
    free(colptrs);

    return res;
}

/* vim: set sw=4: */
