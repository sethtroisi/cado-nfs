/* This program sorts a matrix produced by procrels into another matrix where
   the columns are sorted by increasing weight, in order to speed up the
   linear algebra (in particular when using Gaussian elimination).

   Example of use:
   ./procrels rels > mat
   ./matsort mat > mat2
   ./linalg mat2 > ker
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>

/* FIXME: should be shared with linalg.c.
   Data is stored in a huge table:
   [ a b nb_coeff c0 c1 ... ck a b nb_coeff c0 c1 ..... ] */
typedef struct {
  unsigned int nrows;
  unsigned int ncols;
  int *data;
  double avg_wt;  // average nb of coeff per row.
} sparse_mat_t;

/* the following is adapted from linalg.c */
void
readmat (FILE *file, sparse_mat_t *mat)
{
  int ret;
  int alloced;
  int used = 0;
  int i, j;
  double wt = 0;

  ret = fscanf (file, "%d %d", &(mat->nrows), &(mat->ncols));
  assert (ret == 2);
  alloced = mat->nrows * 10;
  mat->data = (unsigned int *) malloc (alloced * sizeof (unsigned int));
  for (i = 0; i < mat->nrows; ++i) {
    long a;
    unsigned long b;
    int nc;
    ret = fscanf (file, "%d %d", &a, &b);
    assert (ret == 2);
    mat->data[used++] = a;
    mat->data[used++] = b;
    ret = fscanf (file, "%d", &nc);
    assert (ret == 1);
    wt += nc;
    if (used + nc + 3 >= alloced)
      {
        alloced += 1000;
        mat->data = (unsigned int *) realloc (mat->data,
                                              alloced * sizeof (unsigned int));
    }
    mat->data[used++] = nc;
    for (j = 0; j < nc; ++j)
      {
        int x;
        ret = fscanf(file, "%d", &x);
        assert (ret == 1);
        mat->data[used++] = x;
        assert (0 <= x && x < mat->ncols);
      }
  }
  mat->avg_wt = wt / mat->nrows;
  /* sort columns by increasing weight */
}

typedef struct
{
  unsigned int index;
  int weight;
} column_t;

int
cmpcol (const void *c1, const void *c2)
{
  return ((column_t*) c1)->weight - ((column_t*) c2)->weight;
}

void
sort_columns (sparse_mat_t *mat, column_t *C)
{
  int i, j, k, l, nc;

  for (j = 0; j < mat->ncols; j++)
    {
      C[j].index = j;
      C[j].weight = 0;
    }
  for (i = k = 0; i < mat->nrows; i++)
    {
      k += 2; /* skip a and b */
      nc = mat->data[k++];
      for (j = 0; j < nc; j++)
        C[mat->data[k++]].weight ++;
    }
  for (j = 0; j < mat->ncols; j++)
    {
      if (C[j].weight == 0)
        fprintf (stderr, "index %d does not appear\n", j);
      else if (C[j].weight == 1)
        fprintf (stderr, "index %d appears only once\n", j);
    }
  fprintf (stderr, "weight of first column is %d\n", C[0].weight);
  fprintf (stderr, "weight of last column is %d\n", C[mat->ncols - 1].weight);
  qsort (C, mat->ncols, sizeof (column_t), cmpcol);
  /* now invert the permutation: indeed, if column j in now in C[k],
     then C[k].index = j, but we need to replace every j by k in the
     output file, i.e., we want C[j].index = k */
  for (k = 0; k < mat->ncols; k++)
    if (C[k].weight > 0)
      {
        i = k;
        j = C[i].index;
        while (j != k)
          {
            l = C[j].index;
            C[j].index = i;
            C[j].weight = 0;
            i = j;
            j = l;
          }
        C[k].index = i;
        C[k].weight = 0;
      }
}

void
print_matrix (sparse_mat_t *mat, column_t *C)
{
  int i, j, k, nc;

  printf("%d %d\n", mat->nrows, mat->ncols);
  for (i = k = 0; i < mat->nrows; i++)
    {
      nc = mat->data[k + 2];
      printf ("%d %d %d", mat->data[k], mat->data[k + 1], nc); /* a b nc */
      k += 3;
      for (j = 0; j < nc; j++)
        printf (" %d", C[mat->data[k++]].index);
      printf ("\n");
    }
}

int
main (int argc, char *argv[])
{
  FILE *file;
  sparse_mat_t mat;
  column_t *C; /* C[j].weight is the number of occurrences of index j
                  in the relations */

  fprintf (stderr, "%s rev. %s\n", argv[0], REV);

  if (argc == 1)
    file = stdin;
  else if (argc == 2)
    {
      file = fopen (argv[1], "r");
      assert (file != NULL);
    }
  else
    {
      fprintf (stderr, "usage: %s [filename]\n", argv[0]);
      exit (1);
  }
  
  readmat (file, &mat);
  fprintf (stderr, "have read matrix: nrows = %d, ncols = %d\n",
           mat.nrows, mat.ncols);
  fprintf (stderr, "average wt of rows is %f\n", mat.avg_wt);

  C = (column_t*) malloc (mat.ncols * sizeof (column_t));

  sort_columns (&mat, C);

  print_matrix (&mat, C);

  free (mat.data);
  free (C);
  fclose (file);

  return 0;
}
