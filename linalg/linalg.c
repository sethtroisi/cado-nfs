#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>

int kernel(mp_limb_t* mat, mp_limb_t** ker, int nrows, int ncols,
    int limbs_per_row, int limbs_per_col);

// data is stored in a huge table:
//   [ nb_coeff c0 c1 ... ck nb_coeff c0 c1 ..... ]
typedef struct {
  unsigned int nrows;
  unsigned int ncols;
  unsigned int *data;
  double avg_wt;  // average nb of coeff per row.
} sparse_mat_t;


typedef struct {
  unsigned int nrows;
  unsigned int ncols;
  mp_limb_t *data;
  unsigned int limbs_per_row;
  unsigned int limbs_per_col;
} dense_mat_t;

void
readmat(FILE *file, sparse_mat_t *mat) {
  int ret;
  int alloced;
  int used;
  int i, j;
  double wt = 0;
  ret = fscanf(file, "%d %d", &(mat->nrows), &(mat->ncols));
  assert (ret == 2);
  alloced = mat->nrows*10;
  mat->data = (unsigned int *)malloc(alloced*sizeof(unsigned int));
  for (i = 0; i < mat->nrows; ++i) {
    long a;
    unsigned long b;
    int nc;
    ret = fscanf(file, "%d %d", &a, &b);
    assert (ret == 2);
    ret = fscanf(file, "%d", &nc);
    assert (ret == 1);
    wt += nc;
    if (used + nc + 1 >= alloced) {
      alloced += 1000;
      mat->data = (unsigned int *)realloc(mat->data, alloced*sizeof(unsigned int));
    }
    mat->data[used++] = nc;
    for (j = 0; j < nc; ++j) {
      int x;
      ret = fscanf(file, "%d", &x);
      assert (ret == 1);
      mat->data[used++] = x;
    }
  }
  mat->avg_wt = wt / mat->nrows;
}


void setcoeff(dense_mat_t *Dmat, unsigned int i, unsigned int j) {
  int index;
  int j0, j1;
  mp_limb_t mask;
  index = Dmat->limbs_per_row*i;
  j0 = j/GMP_NUMB_BITS;
  index += j0;
  j1 = j - (j0*GMP_NUMB_BITS);
  mask = 1UL<<j1;
  Dmat->data[index] |= mask;
}

void
sparse2densetranspose(dense_mat_t *Dmat, sparse_mat_t *Smat) {
  int i, j;
  unsigned int *ptr;
  Dmat->nrows = Smat->ncols;
  Dmat->ncols = Smat->nrows;
  Dmat->limbs_per_row = (Dmat->ncols / GMP_NUMB_BITS) + 1;
  Dmat->limbs_per_col = (Dmat->nrows / GMP_NUMB_BITS) + 1;
  
  Dmat->data = (mp_limb_t *)malloc(Dmat->limbs_per_row*Dmat->nrows*sizeof(mp_limb_t));
  for (i = 0; i < Dmat->limbs_per_row*Dmat->nrows; ++i)
    Dmat->data[i] = 0UL;

  ptr = Smat->data;
  for (i = 0; i < Smat->nrows; ++i) {
    unsigned int nc = *ptr++;
    for (j = 0; j < nc; ++j) {
      unsigned int coeff = *ptr++;
      setcoeff(Dmat, coeff, i);
    }
  }
}

int main(int argc, char **argv) {
  FILE *file;
  sparse_mat_t mat;
  dense_mat_t dmat;
  mp_limb_t **ker;
  int dim;

  if (argc == 1) {
    file = stdin;
  } else if (argc == 2) {
    file = fopen(argv[1], "r");
    assert (file != NULL);
  } else {
    fprintf(stderr, "usage: %s [filename]\n", argv[0]);
    exit(1);
  }
  
  readmat(file, &mat);
  fprintf(stderr, "have read matrix: nrows = %d, ncols = %d\n", mat.nrows, mat.ncols);
  fprintf(stderr, "average wt of rows is %f\n", mat.avg_wt);

  sparse2densetranspose(&dmat, &mat);

  dim = kernel(dmat.data, NULL, dmat.nrows, dmat.ncols,
      dmat.limbs_per_row, dmat.limbs_per_col);

  fprintf(stderr, "dim of kernel is %d\n", dim);

  free(mat.data);
  fclose(file);
}
