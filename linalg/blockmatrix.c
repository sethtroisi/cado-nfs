#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blockmatrix.h"
#include "bit_matrices.h"
#include "macros.h"

blockmatrix blockmatrix_alloc(unsigned int nrows, unsigned int ncols)
{
    blockmatrix res = malloc(sizeof(struct blockmatrix_s));
    res->nrows = nrows;
    res->ncols = ncols;
    res->nrblocks = iceildiv(nrows, 64);
    res->ncblocks = iceildiv(ncols, 64);
    res->stride = res->nrblocks;
    if (nrows == 0 || ncols == 0) {
        res->mb = NULL;
    } else {
        res->mb = malloc(res->nrblocks * res->ncblocks * sizeof(mat64));
    }
    res->sub = 0;
    return res;
}

void blockmatrix_free(blockmatrix b)
{
    ASSERT_ALWAYS(b->sub == 0);
    free(b->mb); b->mb = NULL;
    b->nrblocks = 0;
    b->ncblocks = 0;
    b->nrows = 0;
    b->ncols = 0;
    b->stride = 0;
    free(b);
}

void blockmatrix_print(blockmatrix b, const char *vname)
{
    FILE * f = fopen("/tmp/debug", "a");
    if (f == NULL) {
        fprintf(stderr, "Cannot append to /tmp/debug");
        return;
    }
    fprintf(f,"%s:=Matrix(GF(2),%u,%u,[\n", vname, b->nrows, b->ncols);
    for (unsigned int i = 0; i < b->nrows; i += 64) {
	for (unsigned int ii = 0; ii < 64 && i + ii < b->nrows; ii++) {
	    if (i + ii)
		fprintf(f,",\n");
	    for (unsigned int j = 0; j < b->ncols; j += 64) {
		for (unsigned int jj = 0; jj < 64 && j + jj < b->ncols; jj++) {
		    mat64_ptr bl = b->mb[i / 64 + (j / 64) * b->stride];
		    if (j + jj)
			fprintf(f,",");
		    fprintf(f,"%lu", (bl[ii] >> jj) & 1UL);
		}
	    }
	}
    }
    fprintf(f,"]);\n");
    fclose(f);
}

void blockmatrix_zero(blockmatrix b)
{
    if (!b->sub) {
        memset(b->mb, 0, b->nrblocks * b->ncblocks * sizeof(mat64));
    } else {
        for(unsigned int j = 0 ; j < b->ncblocks ; j++)  
            memset(b->mb + j * b->stride, 0, b->nrblocks * sizeof(mat64));
    }
}

/* Computes transpose(a) * b */
void blockmatrix_mul_Ta_b(blockmatrix c,
        const blockmatrix a,
        const blockmatrix b)
{
    ASSERT_ALWAYS(c->nrows == a->ncols);
    ASSERT_ALWAYS(c->ncols == b->ncols);
    ASSERT_ALWAYS(a->nrows == b->nrows);

    ASSERT_ALWAYS(a != b);
    ASSERT_ALWAYS(a != c);

    blockmatrix_zero(c);

    for(unsigned int i = 0 ; i < c->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < c->ncblocks ; j++) {
            addmul_TN64_N64(c->mb[i + j * c->stride],
                    (uint64_t *) (a->mb + i * a->stride),
                    (uint64_t *) (b->mb + j * b->stride),
                    a->nrows);
        }
    }
}

void blockmatrix_mul_smallb(blockmatrix c,
        const blockmatrix a,
        const blockmatrix b)
{
    ASSERT_ALWAYS(c->nrows == a->nrows);
    ASSERT_ALWAYS(c->ncols == b->ncols);
    ASSERT_ALWAYS(a->ncols == b->nrows);

    ASSERT_ALWAYS(a != b);
    ASSERT_ALWAYS(a != c);

    blockmatrix_zero(c);

    for(unsigned int j = 0 ; j < b->ncblocks ; j++) {
        for(unsigned int i = 0 ; i < b->nrblocks ; i++) {
            addmul_N64_6464(
                    (uint64_t *) (c->mb + j * c->stride),
                    (uint64_t *) (a->mb + i * a->stride),
                    b->mb[i + j * b->stride], a->nrows);
        }
    }
}

/* Takes the block matrix m, and copy its data at position (i0, j0) inside
 * the tiny matrix.
 */
void blockmatrix_copy_to_flat(uint64_t * tiny, unsigned int stride,
        int i0, int j0, blockmatrix m)
{
    for(unsigned int i = 0 ; i < m->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < m->ncblocks ; j++) {
            mat64_ptr tm = m->mb[i + j * m->stride];
            /* if needed, easy to transpose tm here and swap (i,j) */
            uint64_t * tp = tiny + (i0+i*64) * stride + j0/64 + j;
            for(unsigned int k = 0 ; k < 64 ; k++)
                tp[k*stride] = tm[k];
        }
    }
}

void blockmatrix_copy_transpose_to_flat(uint64_t * tiny, unsigned int stride,
        int i0, int j0, blockmatrix m)
{
    for(unsigned int i = 0 ; i < m->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < m->ncblocks ; j++) {
            mat64 tm;
            transp_6464(tm, m->mb[i + j * m->stride]);
            uint64_t * tp = tiny + (i0+j*64) * stride + j0/64 + i;
            /* Note that the tiny matrix must have space allocated for rows and
             * cols multiples of 64, otherwise disaster will occur */
            for(unsigned int k = 0 ; k < 64 ; k++)
                tp[k*stride] = tm[k];
        }
    }
}

void blockmatrix_copy_transpose_from_flat(blockmatrix m, uint64_t * tiny, unsigned int stride, int i0, int j0)
{
    for(unsigned int i = 0 ; i < m->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < m->ncblocks ; j++) {
            mat64 tm;
            uint64_t * tp = tiny + (i0+j*64) * stride + j0/64 + i;
            for(unsigned int k = 0 ; k < 64 ; k++)
                tm[k] = tp[k*stride];
            transp_6464(m->mb[i + j * m->stride], tm);
        }
    }
}

void blockmatrix_copy_from_flat(blockmatrix m, uint64_t * tiny, unsigned int stride, int i0, int j0)
{
    for(unsigned int i = 0 ; i < m->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < m->ncblocks ; j++) {
            mat64_ptr tm = m->mb[i + j * m->stride];
            uint64_t * tp = tiny + (i0+i*64) * stride + j0/64 + j;
            for(unsigned int k = 0 ; k < 64 ; k++)
                tm[k] = tp[k*stride];
        }
    }
}

void blockmatrix_read_from_flat_file(blockmatrix k, int i0, int j0, const char * name, unsigned int fnrows, unsigned int fncols)
{
    FILE * f = fopen(name, "r");
    ASSERT_ALWAYS(f);
    ASSERT_ALWAYS(i0 + fnrows <= k->nrows);
    ASSERT_ALWAYS(j0 + fncols <= k->ncols);
    ASSERT_ALWAYS(j0 % 64 == 0);
    ASSERT_ALWAYS(fncols % 64 == 0);
    for(unsigned int r = 0 ; r < fnrows ; r++) {
        for(unsigned int g = 0 ; g < fncols ; g+=64) {
            uint64_t v;
            int rc = fread(&v, sizeof(uint64_t), 1, f);
            ASSERT_ALWAYS(rc == 1);
            k->mb[((i0+r)/64) + ((j0+g)/64) * k->stride][(i0+r)%64] = v;
        }
    }
    fclose(f);
}

void blockmatrix_write_to_flat_file(const char * name, blockmatrix k, int i0, int j0, unsigned int fnrows, unsigned int fncols)
{
    FILE * f = fopen(name, "w");
    ASSERT_ALWAYS(f);
    ASSERT_ALWAYS(i0 + fnrows <= k->nrows);
    ASSERT_ALWAYS(j0 + fncols <= k->ncols);
    ASSERT_ALWAYS(j0 % 64 == 0);
    // for our convenience, we allow larger files.
    // ASSERT_ALWAYS(fncols % 64 == 0);
    for(unsigned int r = 0 ; r < fnrows ; r++) {
        for(unsigned int g = 0 ; g < fncols ; g+=64) {
            uint64_t v;
            v = k->mb[((i0+r)/64) + ((j0+g)/64) * k->stride][(i0+r)%64];
            int rc = fwrite(&v, sizeof(uint64_t), 1, f);
            ASSERT_ALWAYS(rc == 1);
        }
    }
    fclose(f);
}

void blockmatrix_transpose(blockmatrix b, blockmatrix a)
{
    ASSERT_ALWAYS(b->nrows == a->ncols);
    ASSERT_ALWAYS(a->nrows == b->ncols);
    for(unsigned int i = 0 ; i < a->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < a->ncblocks ; j++) {
            transp_6464(b->mb[j + i * b->stride], a->mb[i + j * a->stride]);
        }
    }
}

blockmatrix blockmatrix_submatrix(blockmatrix k, int i0, int j0, unsigned int nrows, unsigned int ncols)
{
    ASSERT_ALWAYS(i0 % 64 == 0);
    ASSERT_ALWAYS(j0 % 64 == 0);
    ASSERT_ALWAYS(i0 + nrows <= k->nrows);
    ASSERT_ALWAYS(j0 + ncols <= k->ncols);
    blockmatrix res = malloc(sizeof(struct blockmatrix_s));
    res->nrows = nrows;
    res->ncols = ncols;
    res->stride = k->nrblocks;
    res->nrblocks = iceildiv(nrows, 64);
    res->ncblocks = iceildiv(ncols, 64);
    /* Do ***NOT*** call blockmatrix_free on a submatrix !!! */
    res->mb = k->mb + (i0/64) + (j0/64) * k->stride;
    if (nrows == 0 || ncols == 0) res->mb = NULL;
    res->sub = 1;
    return res;
}
