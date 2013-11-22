#include "cado.h"
#include <stdio.h>

#include "portability.h"
#include "utils.h"
#include "filter_utils.h"
#include "filter_matrix.h"
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

#else /* QUEUE_TYPE == 0 */

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

    while((k > 1) && (Q[k/2][2] >= v)){
	// we are at level > 0 and the father is >= son
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
    if(Q[0][0] >= MERGE_LEVEL_MAX*MERGE_LEVEL_MAX-2)
	fprintf(stderr, "size(Q) -> %d\n", Q[0][0]);
}

// Move Q[1] down.
static void
downQueue(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3], int k)
{
    int x = Q[k][0], y = Q[k][1], v = Q[k][2], j;

    while(k <= Q[0][0]/2){
	// k has at least a left son
	j = 2*k;
	if(j < Q[0][0])
	    // k has a right son
	    if(Q[j][2] > Q[j+1][2])
		j++;
	// at this point, Q[j] is the largest son
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

// height[0] = 0, etc. Returns hmax; *w will contain the minimal sum.
static int
minimalSpanningTreeWithPrim(int *w, int *father, int *height, 
			    int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
			    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
    int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3]; // over-conservative
    int u, s, t, i, nV, hmax = 0;

    // nodes of T
    for(i = 0; i < m; i++){
	father[i] = -1;
	height[i] = 0;
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
#else
    Q[0][0] = 0;
    addAllEdges(Q, u, father, A, m);
#endif
#if DEBUG >= 1
    printQueue(Q);
#endif
    *w = 0;
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
	    *w += A[s][t];
	    father[t] = s;
	    // store new son for s
	    sons[s][0]++;
	    sons[s][sons[s][0]] = t;
	    height[t] = height[s]+1;
	    if(height[t] > hmax)
		hmax = height[t];
	    nV--;
	    if(nV == 0)
		break;
	    addAllEdges(Q, t, father, A, m);
#if DEBUG >= 1
	    printQueue(Q);
#endif
	}
    }
    return hmax;
}

/* given an ideal of weight m, fills the m x m matrix A so that
   A[i,j] is the weight of the sum of the i-th and j-th rows
   containing the ideal, for 0 <= i, j < m */
void
fillRowAddMatrix(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], filter_matrix_t *mat,
                 int m, int32_t *ind, int32_t ideal)
{
    int i, j;

    for(i = 0; i < m; i++)
	A[i][i] = 0;
    for(i = 0; i < m; i++)
	for(j = i+1; j < m; j++){
	    A[i][j] = weightSum(mat, ind[i], ind[j], ideal);
	    A[j][i] = A[i][j];
	}
}

int
minimalSpanningTree(int *w, int *father, int *height, 
		    int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		    int m, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
    return minimalSpanningTreeWithPrim(w, father, height, sons, A, m);
    //    minimalSpanningTreeWithKruskal(A, m);
}

int
minCostUsingMST(filter_matrix_t *mat, int m, int32_t *ind, int32_t j, 
                double *tfill, double *tMST)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX];
    int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];
    int father[MERGE_LEVEL_MAX], height[MERGE_LEVEL_MAX], w;
    // int hmax;

    *tfill = seconds();
    fillRowAddMatrix(A, mat, m, ind, j);
    *tfill = seconds()-*tfill;
    *tMST = seconds();
    /* hmax = */ minimalSpanningTree(&w, father, height, sons, m, A);
    *tMST = seconds()-*tMST;
    return w;
}

void
printMST(int father[MERGE_LEVEL_MAX], int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1], int m)
{
    int i, k;

    for(i = 0; i < m; i++)
	fprintf(stderr, "father[%d] = %d\n", i, father[i]);
    for(i = 0; i < m; i++){
	fprintf(stderr, "Sons of %d:", i);
	for(k = 1; k <= sons[i][0]; k++)
	    fprintf(stderr, " %d", sons[i][k]);
	fprintf(stderr, "\n");
    }
}

