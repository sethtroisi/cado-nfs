#include "cado.h"
#include <stdio.h>

#include "portability.h"
#include "filter_config.h"
#include "utils.h"
#include "merge_replay_matrix.h"
#include "sparse.h"
#include "mst.h"

//////////////////////////////////////////////////////////////////////
// Prim

#define QUEUE_TYPE 1 // 0 for naive, 1 for heap

#if QUEUE_TYPE == 0

// we put fQ in Q[0][0], lQ in Q[0][1].

#define isQueueEmpty(Q) (Q[0][0] == Q[0][1])

static void
popQueue(int *s, int *t, int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3])
{
    int fQ = Q[0][0];

    *s = Q[fQ][0];
    *t = Q[fQ][1];
    Q[0][0]++;
}

#if DEBUG >= 1
static void
printQueue(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3])
{
    int i;

    for(i = Q[0][0]; i < Q[0][1]; i++)
	fprintf(stderr, "Q[%d] = [%d, %d, %d]\n", i,Q[i][0], Q[i][1], Q[i][2]);
}
#endif

static void
addEdge(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3], int u, int v, int Auv)
{
    int i, j, fQ = Q[0][0], lQ = Q[0][1];

    for(i = fQ; i < lQ; i++)
	if(Auv < Q[i][2])
	    break;
    // shift everybody: Auv >= Q[i][2]
    for(j = lQ; j > i; j--){
	Q[j][0] = Q[j-1][0];
	Q[j][1] = Q[j-1][1];
	Q[j][2] = Q[j-1][2];
    }
    Q[i][0] = u;
    Q[i][1] = v;
    Q[i][2] = Auv;
    Q[0][1]++;
}

#else /* heap queue type */

// Q[0][0] contains the number of items in Q[], so that useful part of Q
// is Q[1..Q[0][0]].

#define isQueueEmpty(Q) (Q[0][0] == 0)

#if DEBUG >= 1
static void
printQueue(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3])
{
    int level = 0, imax = 1, i;

    fprintf(stderr, "L0:");
    for(i = 1; i <= Q[0][0]; i++){
	fprintf(stderr, " %d", Q[i][2]);
	if(i == imax){
	    imax = (imax<<1)+1;
	    fprintf(stderr, "\nL%d:", ++level);
	}
    }
    fprintf(stderr, "\n");
}
#endif

static void
upQueue(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3], int k)
{
    int x = Q[k][0], y = Q[k][1], v = Q[k][2];

    while((k > 1) && (Q[k/2][2] > v)){
	// we are at level > 0 and the father is > son
	// the father replaces the son
	Q[k][0] = Q[k/2][0];
	Q[k][1] = Q[k/2][1];
	Q[k][2] = Q[k/2][2];
	k /= 2;
    }
    // we found the place of (x, y, v)
    Q[k][0] = x;
    Q[k][1] = y;
    Q[k][2] = v;
}

static void
addEdge(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3], int u, int v, int Auv)
{
    Q[0][0]++;
    Q[Q[0][0]][0] = u;
    Q[Q[0][0]][1] = v;
    Q[Q[0][0]][2] = Auv;
    upQueue(Q, Q[0][0]);
#if DEBUG >= 1
    if(Q[0][0] >= MERGE_LEVEL_MAX*MERGE_LEVEL_MAX-2)
	fprintf(stderr, "size(Q) -> %d\n", Q[0][0]);
#endif
}

// Move Q[1] down.
static void
downQueue(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3], int k)
{
    int x = Q[k][0], y = Q[k][1], v = Q[k][2], j;

    while(2*k <= Q[0][0]){
	// k has at least a left son
	j = 2*k;
	if(j < Q[0][0])
	    // k has a right son
	    if(Q[j][2] > Q[j+1][2])
		j++;
	// at this point, Q[j] is the smallest son
	if(v <= Q[j][2])
	    break;
	else{
	    // the father takes the place of the son
	    Q[k][0] = Q[j][0];
	    Q[k][1] = Q[j][1];
	    Q[k][2] = Q[j][2];
	    k = j;
	}
    }
    // we found the place of v
    Q[k][0] = x;
    Q[k][1] = y;
    Q[k][2] = v;
}

// smallest edge (s, t) is always in Q[1]
static void
popQueue(int *s, int *t, int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3])
{
    *s = Q[1][0];
    *t = Q[1][1];
    Q[1][0] = Q[Q[0][0]][0];
    Q[1][1] = Q[Q[0][0]][1];
    Q[1][2] = Q[Q[0][0]][2];
    Q[0][0]--;
    downQueue(Q, 1);
}
#endif /* QUEUE_TYPE == 0 */

// Add all neighbors of u, which is already in T; all edges with v in T
// are not useful anymore.
static void
addAllEdges(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3],
	    int u, int *father, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
    int v;

    for(v = 0; v < m; v++)
	if((v != u) && (father[v] < 0))
	    addEdge(Q, u, v, A[u][v]);
}

/* we simply put in father[i] and sons[i][0]
   the values s and t of each edge (s, t) */

/* naive implementation of Prim's algorithm */
static int
minimalSpanningTreePrimNaive (int *start, int *end,
                              int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
  int n, k, i, j, w = 0, imin, jmin, wmin;
  static int S[MERGE_LEVEL_MAX], T[MERGE_LEVEL_MAX];

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
    static int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3]; // over-conservative
    int u, s, t, i, nV, w;
    static int father[MERGE_LEVEL_MAX];
    static int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];

    // nodes of T
    for(i = 0; i < m; i++){
	father[i] = -1;
	sons[i][0] = 0; // last occupied position for a son of i
    }
    u = 0;
    father[u] = u; // for the root, isn't it?
    nV = m-1;
#if QUEUE_TYPE == 0
    Q[0][0] = 1;
    Q[0][1] = 1;
    addAllEdges(Q, u, father, A, m);
    // active part of Q is Q[fQ..lQ[
    ASSERT(Q[0][0] == 1);
    ASSERT(Q[0][1] == m);
#else /* heap */
    Q[0][0] = 0;
    addAllEdges(Q, u, father, A, m);
#endif
#if DEBUG >= 1
    printQueue(Q);
#endif
    w = 0;
    while(! isQueueEmpty(Q)){
	// while queue is non empty
	// pop queue
        popQueue(&s, &t, Q);
#if DEBUG >= 1
	fprintf(stderr, "Popping a = (%d, %d) of weight %d\n", s, t, A[s][t]);
#endif
	if(father[t] == -1){
	    // t does not close a cycle
	    // T[u] <- T[u] union (s, t)
#if DEBUG >= 1
	    fprintf(stderr, "new edge: (%d, %d) of weight %d\n",s,t,A[s][t]);
#endif
	    w += A[s][t];
	    father[t] = s;
            start[m - 1 - nV] = s;
            end[m - 1 - nV] = t;
	    // store new son for s
	    sons[s][0]++;
	    sons[s][sons[s][0]] = t;
	    nV--;
	    if(nV == 0)
		break;
	    addAllEdges(Q, t, father, A, m);
#if DEBUG >= 1
	    printQueue(Q);
#endif
	}
    }
    return w;
}

#if 0
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
mpz_set_bits (mpz_t t, index_t *r, int k, index_t *l, int n, int32_t ideal)
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
      if (l[u] != j)
        {
          printf ("ideal=%d u=%d v=%d l[u]=%u j=%u\n", ideal, u, v, l[u], j);
        }
      ASSERT_ALWAYS(l[u] == j);
      mpz_setbit (t, u);
    }
}

/* special code for factorization */
void
fillRowAddMatrix2(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], filter_matrix_t *mat,
                 int m, int32_t *ind, int32_t ideal MAYBE_UNUSED)
{
  int i, j, n;
  index_t *l = get_ideals (&n, mat->rows, ind, m);
  mpz_t *t, u;

  mpz_init (u);
  t = malloc (m * sizeof (mpz_t));
  for (i = 0; i < m; i++)
    {
      mpz_init(t[i]); /* set to 0 */
      mpz_set_bits (t[i], mat->rows[ind[i]] + 1, mat->rows[ind[i]][0], l, n,
                    ideal);
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
}
#endif

/* given an ideal of weight m, fills the m x m matrix A so that
   A[i][j] is the weight of the sum of the i-th and j-th rows
   containing the ideal, for 0 <= i, j < m */
void
fillRowAddMatrix(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], filter_matrix_t *mat,
                 int m, int32_t *ind, int32_t ideal)
{
    int i, j;

    /* A[i][i] is not used, thus we don't initialize it. */
    for(i = 0; i < m; i++)
	for(j = i+1; j < m; j++){
	    A[i][j] = weightSum (mat, ind[i], ind[j], ideal);
	    A[j][i] = A[i][j];
	}
}

int
minimalSpanningTree(int *start, int *end,
		    int m, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
  if (m <= 16)
    return minimalSpanningTreePrimNaive (start, end, A, m);
  else
    return minimalSpanningTreeWithPrim (start, end, A, m);
}

int
minCostUsingMST (filter_matrix_t *mat, int m, int32_t *ind, int32_t j)
{
    static int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], w;
    static int sons[MERGE_LEVEL_MAX];
    static int father[MERGE_LEVEL_MAX];

    fillRowAddMatrix (A, mat, m, ind, j);
    w = minimalSpanningTree (father, sons, m, A);
    /* w is the cost of all merges, we should subtract the cost of the
       initial relations */
    for (int i = 0; i < m; i++)
      w -= matLengthRow(mat, ind[i]);
    return w;
}

void
printMST (int father[MERGE_LEVEL_MAX],
          int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1], int m, int32_t *ind)
{
  for (int i = 0; i < m; i++)
    printf ("father=%d(%d) son=%d(%d)\n", father[i], ind[father[i]],
            sons[i][0], ind[sons[i][0]]);
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
