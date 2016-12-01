#include "cado.h"
#include <stdio.h>

#include "portability.h"
#include "filter_config.h"
#include "utils.h"
#include "merge_replay_matrix.h"
#include "sparse.h"
#include "mst.h"

#ifdef TIMINGS
double tfill[MERGE_LEVEL_MAX] = {0,}, tmst[MERGE_LEVEL_MAX] = {0,};
double nfill[MERGE_LEVEL_MAX] = {0,};
#endif

//////////////////////////////////////////////////////////////////////
// heap data structure to store the remaining edges to be considered
// in Prim's algorithm

// #define DEBUG 1

/* The queue is a heap of 32-bit entries, with the upper 16 bits being the
   weight of the edge, the next 8 bits being the end vertex, and the last
   8 bits the start vertex. In such a way by comparing two values we compare
   the weight. The heap starts at Q[1], Q[0] being the number of elements of
   the heap, thus the heap values are Q[1]...Q[Q[0]]. */

#define isQueueEmpty(Q) (Q[0] == 0)
#define W(x) ((x) >> 16)
#define S(x) ((x) & 255)
#define T(x) (((x) >> 8) & 255)

#define MAX_QUEUE_SIZE (MERGE_LEVEL_MAX*MERGE_LEVEL_MAX/2)

#if DEBUG >= 1
static void
printQueue (uint32_t *Q)
{
  for (unsigned int i = 1; i <= Q[0]; i++)
    fprintf (stderr, "w=%d s=%d t=%d\n", W(Q[i]), S(Q[i]), T(Q[i]));
}
#endif

static void
upQueue (uint32_t *Q, int k)
{
  uint32_t x = Q[k];

  while ((k > 1) && (x < Q[k/2]))
    {
      // we are at level > 0 and the father is > son
      // the father replaces the son
      Q[k] = Q[k/2];
      k /= 2;
    }
  // we found the place of x
  Q[k] = x;
}

static void
addEdge (uint32_t *Q, int s, int t, int Auv)
{
  Q[0] ++;
  ASSERT(Q[0] < MAX_QUEUE_SIZE);
  Q[Q[0]] = (Auv << 16) | (t << 8) | s;
  upQueue (Q, Q[0]);
}

// Move Q[k] down.
static void
downQueue (uint32_t *Q, unsigned int k)
{
  uint32_t x = Q[k];
  unsigned int j;

  /* the left son of k is 2k, the right son is 2k+1 */
  while (2*k <= Q[0])
    {
      // k has at least a left son
      j = 2*k;
      if (j < Q[0])
        // k has a right son
        if (Q[j] > Q[j+1])
          j++;
      // at this point, Q[j] is the smallest son
      if (x <= Q[j])
        break; /* x is smaller than smallest son */
      else
        {
          // the father takes the place of the son
          Q[k] = Q[j];
          k = j;
        }
    }
  // we found the place of x
  Q[k] = x;
}

static void
popQueue(int *s, int *t, uint32_t *Q)
{
  *s = S(Q[1]);
  *t = T(Q[1]);
  Q[1] = Q[Q[0]];
  Q[0]--;
  downQueue (Q, 1);
}

// Add all edges (u,v) with v in V
static void
addAllEdges (uint32_t *Q, int u, int *V, int nV,
             int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
  for (int i = 0; i < nV; i++)
    addEdge (Q, u, V[i], A[u][V[i]]);
}

/* naive implementation of Prim's algorithm:
   we put in start[i] and end[i] the values s and t of each edge (s, t) */
static int
minimalSpanningTreePrimNaive (int *start, int *end,
                              int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
  int n, k, i, j, w = 0, imin, jmin, wmin;
  int S[MERGE_LEVEL_MAX], T[MERGE_LEVEL_MAX];

  /* S is the set of vertices in the current tree, T is the set of remaining
     vertices */
  S[0] = 0; /* S = {0} */
  for (i = 1; i < m; i++)
    T[i-1] = i; /* T = {1, 2, ..., m-1} */
  n = 1;     /* number of vertices in S */
  k = m - 1; /* number of vertices in T */
  while (k)
    {
      int s, t;
      /* find the edge with minimal weight from S to T */
      wmin = INT_MAX;
      imin = jmin = -1;
      for (i = 0; i < n; i++)
        for (j = 0; j < k; j++)
          if (A[S[i]][T[j]] < wmin)
            {
              imin = i;
              jmin = j;
              wmin = A[S[i]][T[j]];
            }
      ASSERT(imin != -1 && jmin != -1);
      s = S[imin];
      t = T[jmin];
      w += wmin;
      S[n] = t;
      T[jmin] = T[k - 1];
      start[n-1] = s;
      end[n-1] = t;
      n++;
      k--;
    }
  return w;
}

/* Return the weight of the minimal spanning tree.
   For each (s,t) which is part of the tree, we have
   start[i] = s and end[i] = t for 0 <= i < m-1. */
static int
minimalSpanningTreeWithPrim(int *start, int *end,
			    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
  /* the queue contains at most m-1 edges when nU=1 (those connected to the
     root node), then at most 2*(m-2) when nU=2 (we remove one node and add
     m-2), ... More generally whehn nU=k it is at most k*m - k*(k-1)/2 - (k-1).
     The maximum is for k=m-2 or m-2 when it equals m^2/2-3/2*m+2 < m^2/2. */
    uint32_t Q[MAX_QUEUE_SIZE];
    int u, s, t, i, nU, nV, w;
    int V[MERGE_LEVEL_MAX]; /* edges remaining to be dealt with */
    int index[MERGE_LEVEL_MAX];

    // nodes of T
    for(i = 0; i < m; i++){
        V[i] = i;
        index[i] = i; /* V[index[i]] = i if i is in V, -1 otherwise */
    }
    u = 0;
    index[u] = -1;
    index[V[m-1]] = u;
    V[u] = V[m-1];
    nU = 1;   /* number of edges already in the MST */
    nV = m-1; /* number of edges remaining to be dealt with */
    Q[0] = 0; /* heap is empty initially */
    addAllEdges (Q, u, V, nV, A);
#if DEBUG >= 1
    printQueue(Q);
#endif
    w = 0;
    while (!isQueueEmpty(Q)) /* while queue is non empty */
      {
        popQueue (&s, &t, Q); // pop queue
#if DEBUG >= 1
	fprintf(stderr, "Popping a = (%d, %d) of weight %d\n", s, t, A[s][t]);
#endif
	if (index[t] != -1){
	    // t does not close a cycle
	    // T[u] <- T[u] union (s, t)
#if DEBUG >= 1
	    fprintf(stderr, "new edge: (%d, %d) of weight %d\n",s,t,A[s][t]);
#endif
	    w += A[s][t];
            start[nU - 1] = s;
            end[nU - 1] = t;
            ASSERT(V[index[t]] == t);
            index[V[nV-1]] = index[t];
            V[index[t]] = V[nV-1];
            index[t] = -1;
	    nV--;
            nU++;
	    if (nV == 0)
		break;
            addAllEdges (Q, t, V, nV, A);
#if DEBUG >= 1
	    printQueue(Q);
#endif
	}
    }
    return w;
}

#if !defined(FOR_DL) && 0
/* special fillRowAddMatrix() code for factorization */
/* return a list l of all (unique and sorted) ideals appearing in relations
   ind[0] to ind[j-1], and put in *n the length of l */
static index_t*
get_ideals (int *n, index_t **rows, int32_t *ind, int j)
{
  if (j == 1)
    {
      *n = rowLength (rows, ind[0]);
      return rows[ind[0]] + 1;
    }
  int k = j / 2, n1, n2, k1, k2;
  index_t *l1 = get_ideals (&n1, rows, ind, k);
  index_t *l2 = get_ideals (&n2, rows, ind + k, j - k);
  /* merge the two lists */
  index_t *l;
  l = malloc ((n1 + n2) * sizeof (index_t));
  for (k1 = 0, k2 = 0, k = 0; k1 < n1 && k2 < n2;)
    {
      if (l1[k1] <= l2[k2])
        {
          k2 += (l1[k1] == l2[k2]);
          l[k++] = l1[k1++];
        }
      else /* l2[k2] < l1[k1] */
        l[k++] = l2[k2++];
    }
  while (k1 < n1)
    l[k++] = l1[k1++];
  while (k2 < n2)
    l[k++] = l2[k2++];
  *n = k;
  if (j >= 4) /* then floor(j/2) >= 2 */
    free (l1);
  if (j >= 3) /* then ceil(j/2) >= 2 */
    free (l2);
  return l; /* not need to reallocate since l will be freed right away */
}

/* for each ideal j in r[0] to r[k-1], find the unique u such that l[u] = j,
   and set bit u of t to 1 */
static void
mpz_set_bits (mpz_t t, index_t *r, int k, index_t *l, int n)
{
  int i, u, v;
  index_t j;

  for (i = 0; i < k; i++)
    {
      j = r[i];
      u = 0;
      v = n;
      /* we have l[u] <= j < l[v] */
      while (u + 1 < v)
        {
          int w = (u + v) / 2;
          if (j < l[w])
            v = w;
          else
            u = w;
        }
      /* now l[u] = j */
      ASSERT(l[u] == j);
      mpz_setbit (t, u);
    }
}

void
fillRowAddMatrix(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], filter_matrix_t *mat,
                 int m, int32_t *ind, int32_t ideal MAYBE_UNUSED)
{
  int i, j, n;
  index_t *l;
  mpz_t *t, u;

#ifdef TIMINGS
  tfill[m] -= seconds ();
  nfill[m] += 1.0;
#endif
  l = get_ideals (&n, mat->rows, ind, m);
  mpz_init (u);
  t = malloc (m * sizeof (mpz_t));
  for (i = 0; i < m; i++)
    {
      mpz_init (t[i]); /* set to 0 */
      mpz_set_bits (t[i], mat->rows[ind[i]] + 1, mat->rows[ind[i]][0], l, n);
    }
  for (i = 0; i < m; i++)
    {
      A[i][i] = 0;
      for (j = i + 1; j < m; j++)
        {
          mpz_xor (u, t[i], t[j]);
          A[i][j] = mpz_popcount (u);
          A[j][i] = A[i][j];
        }
    }
  for (i = 0; i < m; i++)
    mpz_clear (t[i]);
  free (t);
  free (l);
  mpz_clear (u);
#ifdef TIMINGS
  tfill[m] += seconds ();
#endif
}
#else
/* given an ideal of weight m, fills the m x m matrix A so that
   A[i][j] is the weight of the sum of the i-th and j-th rows
   containing the ideal, for 0 <= i, j < m */
void
fillRowAddMatrix(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], filter_matrix_t *mat,
                 int m, int32_t *ind, int32_t ideal)
{
    int i, j;

#ifdef TIMINGS
    tfill[m] -= seconds ();
    nfill[m] += 1.0;
#endif
    /* A[i][i] is not used, thus we don't initialize it. */
    for(i = 0; i < m; i++)
	for(j = i+1; j < m; j++){
	    A[i][j] = weightSum (mat, ind[i], ind[j], ideal);
	    A[j][i] = A[i][j];
	}
#ifdef TIMINGS
    tfill[m] += seconds ();
#endif
}
#endif

int
minimalSpanningTree(int *start, int *end,
		    int m, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
  int ret;

#ifdef TIMINGS
  tmst[m] -= seconds ();
#endif
  if (m <= 16)
    ret = minimalSpanningTreePrimNaive (start, end, A, m);
  else
    ret = minimalSpanningTreeWithPrim (start, end, A, m);
#ifdef TIMINGS
  tmst[m] += seconds ();
#endif
  return ret;
}

int
minCostUsingMST (filter_matrix_t *mat, int m, int32_t *ind, int32_t j)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], w;
    int sons[MERGE_LEVEL_MAX];
    int father[MERGE_LEVEL_MAX];

    fillRowAddMatrix (A, mat, m, ind, j);
    w = minimalSpanningTree (father, sons, m, A);
    /* w is the cost of all merges, we should subtract the cost of the
       initial relations */
    for (int i = 0; i < m; i++)
      w -= matLengthRow(mat, ind[i]);
    return w;
}

void
printMST (int *father, int *sons, int m, int32_t *ind)
{
  for (int i = 0; i < m; i++)
    printf ("father=%d(%d) son=%d(%d)\n", father[i], ind[father[i]],
            sons[i], ind[sons[i]]);
}

#if 0
/* for each active ideal m, and each pair of active rows i and j containing m,
   if weight (row[i] + row[j]) < weight (row[i]), replaces row[i] by
   row[i] + row[j] */
void
finalOptimize (filter_matrix_t *mat)
{
  uint32_t h, n;
  index_t i, j, k;
  unsigned int w;
  uint32_t *ind;
  typerow_t **rows = mat->rows;
  static unsigned int diff = 0;

  printf ("enter finalOptimize, weight %lu\n", mat->weight);
  for (h = 0; h < mat->ncols; h++)
    {
      if (mat->R[h] != NULL) /* ideal is still active */
        {
          int ok = 1;
          n = mat->R[h][0];
          ind = mat->R[h] + 1;
          for (i = 0; i < n && ok; i++)
            for (j = i + 1; j < n && ok; j++)
              {
                w = weightSum (mat, ind[i], ind[j], h);
                if (rowLength (rows, ind[i]) >= rowLength (rows, ind[j]))
                  k = i;
                else
                  k = j;
                /* k is the index of the heaviest row */
                if (w < rowLength (rows, ind[k]))
                  {
                    if (rowLength (rows, ind[k]) - w > diff)
                      {
                        diff = rowLength (rows, ind[k]) - w;
                        printf ("reduce weight by %u\n", diff);
                      }
                    addRowsAndUpdate (mat, ind[k], ind[(i+j)-k], h);
                    ok = 0; /* since addRowsAndUpdate modifies the R structure,
                               we stop working on that ideal */
                  }
              }
        }
    }
}
#endif
