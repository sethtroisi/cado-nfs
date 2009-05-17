#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>
#include <stddef.h>

#include "macros.h"
#include "readmat.h"
#include "readmat-easy.h"

#define CHECK_REALLOC(ptr__,alloc__,size__,n__) do {    		\
        if (size__ + n__ >= alloc__) {					\
            alloc__ = size__ + n__;                                     \
            alloc__ += alloc__ / 2;                                     \
            ptr__ = realloc(ptr__, alloc__ * sizeof(uint32_t));		\
        }								\
    } while (0)

void read_easy(const char * filename,
        uint32_t ** p_direct, uint32_t ** p_transposed,
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
    colweights = malloc(smat->ncols * sizeof(uint32_t));
    memset(colweights, 0, smat->ncols * sizeof(uint32_t));

    data = NULL;

    if (p_nr) *p_nr = smat->nrows;
    if (p_nc) *p_nc = smat->ncols;

    /* do it the bovine way. read just about everything in memory. */
    for(unsigned int i = 0 ; i < smat->nrows ; i++) {
        read_matrix_row(f,smat,smat->data,1);
        CHECK_REALLOC(data, alloc, size, 1 + smat->data[0]);
        data[size++] = smat->data[0];
        /* read_matrix_row now returns a sorted row */
        for(unsigned int j = 0 ; j < smat->data[0] ; j++) {
            data[size++] = smat->data[1+j];
            colweights[smat->data[1+j]]++;
        }
    }
    fclose(f);

    if (p_direct) {
        *p_direct = data;
    }

    if (p_transposed == NULL) {
        free(colweights);
        return;
    }

    uint32_t ** colptrs;
    colptrs = malloc(smat->ncols * sizeof(uint32_t *));

    if (smat->ncols > smat->nrows) {
        size += smat->ncols - smat->nrows;
    }

    uint32_t * res = malloc(size * sizeof(uint32_t));
    uint32_t * ptr = res;
    for(unsigned int j = 0 ; j < smat->ncols ; j++) {
        uint32_t w = colweights[j];
        ASSERT((size_t) (ptr - res) < size);
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

    if (p_direct == NULL) {
        free(data);
    }

    ASSERT(p_transposed != NULL);
    *p_transposed = res;

    free(colweights);
    free(colptrs);
}

/* vim: set sw=4: */
