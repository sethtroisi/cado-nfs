/* Program to reduce the excess of relations.
   Usage: matprune -keep <value> mat > new_mat
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>
#include <string.h>

typedef struct {
  long a;
  unsigned long b;
  int len;    /* number of non-zero entries */
  int *val;   /* array of size len containing column indices of non-zero
                 entries */
} rel_t;

typedef struct {
  int ind;        /* row index of starting of component */
  int dif;        /* edges[i] - nodes[i] >= -1 */
  int edg;        /* edges[i] */
} component_t;

typedef struct {
  int nrows;
  int ncols;
  int rem_nrows;     /* number of remaining rows */
  int rem_ncols;     /* number of remaining columns */
  rel_t *data;
  unsigned long* wt; /* weight of prime j, <= 1 for a deleted prime */
  unsigned long* ad; /* ad[j] contains the sum of indices of all rows
                        containing prime j; when wt[j]=2, and we know
                        that row i contains j, then the other row is
                        simply k = ad[j] - i. Overflow is not a problem,
                        since the value is always exact mod 2^32 or 2^64. */
  int* nodes;        /* nodes[i] is the number of rows in the connected
                        component of row i, in the graph where two rows
                        are connected if they share one prime appearing only
                        twice; nodes[i] = -1 for a deleted relation */
  int* edges;         /* edges[i] is the number of edges in the connected
                         component of row i */
  int alloc;          /* allocated size for connected components */
  int ncomps;         /* number of connected components */
  component_t* comps; /* list of connected components */
} sparse_mat_t;

/* Reads a matrix file, and puts in mat->wt[j], 0 <= j < ncols, the
   weight of column j (adapted from matsort.c), and in mat->ad[j],
   0 <= j < ncols, the sum of row indices where j appears. */
void
readmat (FILE *file, sparse_mat_t *mat)
{
  int ret;
  int i, j, l2;
  long a;
  unsigned long b;
  int nc, x;

  ret = fscanf (file, "%d %d", &(mat->nrows), &(mat->ncols));
  assert (ret == 2);

  fprintf (stderr, "Reading matrix of %d rows and %d columns: excess is %d\n",
           mat->nrows, mat->ncols, mat->nrows - mat->ncols);
  mat->rem_nrows = mat->nrows;

  mat->wt = (unsigned long*) malloc (mat->ncols * sizeof (unsigned long));
  mat->ad = (unsigned long*) malloc (mat->ncols * sizeof (unsigned long));
  for (j = 0; j < mat->ncols; j++)
    mat->wt[j] = mat->ad[j] = 0;

  mat->data = (rel_t*) malloc (mat->nrows * sizeof (rel_t));

  for (i = 0; i < mat->nrows; i++)
    {
      ret = fscanf (file, "%ld %lu", &a, &b);
      assert (ret == 2);
      mat->data[i].a = a;
      mat->data[i].b = b;
      ret = fscanf (file, "%d", &nc);
      assert (ret == 1);
      mat->data[i].len = nc;
      mat->data[i].val = (int*) malloc (nc * sizeof (int));
      for (j = 0; j < nc; j++)
        {
          ret = fscanf(file, "%d", &x);
          assert (ret == 1);
          assert (0 <= x && x < mat->ncols);
          mat->data[i].val[j] = x;
          mat->wt[x] ++;
          mat->ad[x] += i; /* computed mod 2^32 or 2^64 */
        }
    }

  /* count the number of primes appearing at least twice */
  for (j = l2 = 0; j < mat->ncols; j++)
    {
      if (mat->wt[j] >= 2)
        l2 ++;
      else
        fprintf (stderr, "Warning: prime %d appears only %lu time(s)\n", j,
                 mat->wt[j]);
    }
  fprintf (stderr, "Found %d primes appearing at least twice\n", l2);
  mat->rem_ncols = l2;
}

/* visit connected component of relation i, and compute nodes[i] and
   edges[i], which is twice the numbers of edges of that component */
void
visit (int i, sparse_mat_t *mat)
{
  int m, j, k;

  mat->nodes[i] = 1; /* also mark as visited */
  mat->edges[i] = 0;
  for (m = 0; m < mat->data[i].len; m++)
    {
      j = mat->data[i].val[m];
      if (mat->wt[j] == 2)
        {
          k = mat->ad[j] - i;
          mat->edges[i] ++;
          if (mat->nodes[k] == 0)
            {
              visit (k, mat);
              mat->nodes[i] += mat->nodes[k];
              mat->edges[i] += mat->edges[k];
            }
        }
    }
}

/* Delete relations and primes in component of i.
   FIXME: if some wt[j] becomes 2, we should join the connected components
   of the two remaining relations containing j.
*/
void
delete (int i, sparse_mat_t *mat)
{
  int m, j, k;

  mat->nodes[i] = -1; /* mark as deleted and visited */
  mat->rem_nrows --;
  for (m = 0; m < mat->data[i].len; m++)
    {
      j = mat->data[i].val[m];
      mat->ad[j] -= i;
      mat->wt[j] --;
      if (mat->wt[j] == 1)
        {
          mat->rem_ncols --;
          k = mat->ad[j];
          if (mat->nodes[k] >= 0)
            delete (k, mat);
        }
    }
}

/* Compares two connected components: the best one is the one with the
   largest edge-node difference (i.e., more primes than relations are
   deleted), and in case of equality, with the largest value of edge.
   With qsort, the best component will be last.
*/
int
compare (const void *v1, const void *v2)
{
  component_t *c1 = (component_t*) v1;
  component_t *c2 = (component_t*) v2;

  if (c1->dif != c2->dif)
    return c1->dif - c2->dif;
  else
    return c1->edg - c2->edg;
}

/* Compute connected components. Assume mat->nodes and mat->edges have
   already been allocated, with mat->nrows entries each. */
void
compute_components (sparse_mat_t *mat)
{
  int i;

  for (i = 0; i < mat->nrows; i++)
    if (mat->nodes[i] >= 0)
      mat->nodes[i] = 0;
  /* we use nodes[i] for visited[i], since nodes[i] >= 1 */
  mat->ncomps = 0;
  for (i = 0; i < mat->nrows; i++)
    if (mat->nodes[i] == 0)
      {
        visit (i, mat);
        assert ((mat->edges[i] & 1) == 0);
        mat->edges[i] >>= 1;
        if (mat->ncomps >= mat->alloc)
          {
            mat->alloc += mat->alloc / 10 + 1;
            mat->comps = (component_t*) realloc (mat->comps,
                                            mat->alloc * sizeof (component_t));
          }
        mat->comps[mat->ncomps].ind = i;
        mat->comps[mat->ncomps].dif = mat->edges[i] - mat->nodes[i];
        mat->comps[mat->ncomps].edg = mat->edges[i];
        mat->ncomps ++;
      }
  qsort (mat->comps, mat->ncomps, sizeof (component_t), compare);
#ifdef DEBUG
  fprintf (stderr, "Found %d connected components\n", mat->ncomps);
#endif
}

/* We consider the graph whose vertices are rows of the matrix
   (i.e., relations), where there is an edge between vertex i and vertex k
   if there is a common prime of index j in relations i and k, which appear
   nowhere else. (Note this is a multigraph, since several primes can appear
   only in i and k.)

   For each vertex i, we compute the number of vertices nod[i] in the
   connected component of i (nod[i] = 1 if no prime j appears only in i
   and in another relation k), and the number of edges edg[i] in that
   connected component (edg[i] = 0 if no prime j appears only in i and
   in another relation k).
   The algorithm used to compute nod[i] and edg[i] is as follows:

     Input: wt[j] is the weight of column j
            ad[j] is the sum i+k of relations where j appears (wt[j]=2)
            visited[i] is initialized to 0
     for i from 0 to rows-1 do
        if visited[i]=0 then visit(i) end if
     end for

     Visit := function(i) # return number of nodes of connected component of i
                          # and twice the number of edges
        visited[i]=1
        nodes = 1 # number of nodes of connected component of i
        edges = 0 # number of half-edges of connected component of i
        for j in rel[i] do
           if wt[j] = 2 then
              k = ad[j] - i # other relation where j appears
              edges ++ # count the edge i->k
              if visited[k]=0 then
                 (nodes,edges) += visit(k)
              end if
           end if
        return (nodes,edges)
     end function
*/
void
prune (sparse_mat_t *mat, int keep)
{
  int i, t, excess, old_nrows, old_ncols, mid;
  component_t *c;

  mat->nodes = (int*) malloc (mat->nrows * sizeof (int));
  mat->edges = (int*) malloc (mat->nrows * sizeof (int));

  for (i = 0; i < mat->nrows; i++)
    mat->nodes[i] = 0;

  mat->alloc = 0;
  mat->comps = NULL;

  compute_components (mat);

  excess = mat->nrows - mat->ncols;
  mid = (excess + keep) / 2;
  while (mat->ncomps > 0 && excess + mat->comps[mat->ncomps - 1].dif >= keep)
    {
      t = mat->ncomps - 1;
      c = mat->comps + t;
      i = c->ind;
      /* Warning: this component might already have been deleted */
      if (mat->nodes[i] >= 0) /* it was not deleted */
        {
          old_nrows = mat->rem_nrows;
          old_ncols = mat->rem_ncols;
          delete (i, mat);
#ifdef DEBUG
          fprintf (stderr, "Removed component %d with %d nodes and %d edges\n",
                   i, old_nrows - mat->rem_nrows, old_ncols - mat->rem_ncols);
          fprintf (stderr, "Remains %d relations and %d primes (excess %d)\n",
              mat->rem_nrows, mat->rem_ncols, mat->rem_nrows - mat->rem_ncols);
#endif
          excess += c->dif;
          if (excess <= mid)
            {
              fprintf (stderr, "Remains %d relations and %d primes (excess %d)\n",
              mat->rem_nrows, mat->rem_ncols, mat->rem_nrows - mat->rem_ncols);
              compute_components (mat);
              mid = (excess + keep) / 2;
            }
        }
      mat->ncomps = t;
    }
}

/* Renumber columns to 0..rem_ncols-1. The permutation is stored in
   mat->wt. */
void
renumber_cols (sparse_mat_t *mat)
{
  unsigned long *perm; /* perm[j] = k if old column j coes to new column k */
  int j, s;

  perm = mat->wt;
  for (j = s = 0; j < mat->ncols; j++)
    if (mat->wt[j] >= 2)
      perm[j] = s++;
  assert (s == mat->rem_ncols);
}

void
print_matrix (sparse_mat_t *mat)
{
  int i, j;
  unsigned long *perm = mat->wt;
  rel_t r;

  printf("%d %d\n", mat->rem_nrows, mat->rem_ncols);
  for (i = 0; i < mat->nrows; i++)
    if (mat->nodes[i] >= 0)
      {
        r = mat->data[i];
        printf ("%ld %lu %d", r.a, r.b, r.len);
        for (j = 0; j < r.len; j++)
          printf (" %lu", perm[r.val[j]]);
        printf ("\n");
      }
}

int
count_remaining_cols (sparse_mat_t *mat)
{
  int s, j;
  for (j = s = 0; j < mat->ncols; j++)
    if (mat->wt[j] >= 2)
      s++;
  return s;
}

int
main (int argc, char *argv[])
{
  FILE *file;
  sparse_mat_t mat;
  int keep = 128; /* number of excess to keep */
  char **argv0 = argv;

  fprintf (stderr, "%s rev. %s\n", argv[0], REV);

  if (argc >= 3 && strcmp (argv[1], "-keep") == 0)
    {
      keep = atoi (argv[2]);
      argc -= 2;
      argv += 2;
    }

  if (argc == 1)
    file = stdin;
  else if (argc == 2)
    {
      file = fopen (argv[1], "r");
      assert (file != NULL);
    }
  else
    {
      fprintf (stderr, "usage: %s [-keep nnn] [filename]\n", argv0[0]);
      exit (1);
  }
  
  readmat (file, &mat);

  prune (&mat, keep);

  renumber_cols (&mat);

  print_matrix (&mat);

  free (mat.data);
  free (mat.wt);
  free (mat.ad);
  free (mat.nodes);
  free (mat.edges);
  free (mat.comps);
  fclose (file);

  return 0;
}
