/* 
 * Program: pruning from matprune.c converted to the use of merge.c by FM
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>

#define WANT_ASSERT

#define DEBUG 0

#include "utils/utils.h"
#include "sparse.h"
#include "merge.h"
#include "prune.h"

typedef struct {
  int ind;        /* row index of starting of component */
#if 0
  int dif;        /* edges[i] - nodes[i] >= -1 */
#endif
  int edg;        /* edges[i] */
} component_t;

// We know that M[i][j] exists and j has weight 2. We need the other row
// containing j.
int
getOtherRow(sparse_mat_t *mat, int i, int j)
{
    int k;

    for(k = 1; k <= mat->R[j][0]; k++)
	if((mat->R[j][k] != -1) && (mat->R[j][k] != i))
	    break;
    return mat->R[j][k];
}

/* visit connected component of relation i, and compute nodes[i] and
   edges[i], which is twice the numbers of edges of that component */
// TODO: force i to have some j of weight 2...
void
visit(int i, int *nodes, int *edges, sparse_mat_t *mat)
{
    int m, j, k;
    
    nodes[i] = 1; /* also mark as visited */
    edges[i] = 0;
#if USE_TAB == 0
    for(m = 0; m < mat->data[i].len; m++){
#else
    for(m = 1; m <= lengthRow(mat, i); m++){
#endif
	j = cell(mat, i, m);
	if(mat->wt[j] == 2){
	    edges[i]++;
	    k = getOtherRow(mat, i, j);
	    if(nodes[k] == 0){
		visit(k, nodes, edges, mat);
		nodes[i] += nodes[k];
		edges[i] += edges[k];
            }
        }
    }
}

/* Delete relations and primes in component of i.
   FIXME: if some wt[j] becomes 2, we should join the connected components
   of the two remaining relations containing j.
*/
void
delete(int i, int *nodes, sparse_mat_t *mat)
{
    int j, k, v;
    
    nodes[i] = -1; /* mark as deleted and visited */
    // inspired from removeRowSWAR(mat, i) but for the recursive call
    // to delete...!
    mat->weight -= lengthRow(mat, i);
#if USE_TAB == 0
    for(k = 0; k < mat->data[i].len; k++){
#else
    for(k = 1; k <= lengthRow(mat, i); k++){
#endif
	j = cell(mat, i, k);
	removeCellSWAR(mat, i, j);
	if(mat->wt[j] == 1){
	    v = getOtherRow(mat, i, j);
	    if(nodes[v] >= 0)
		delete(v, nodes, mat);
	}
    }
    destroyRow(mat, i);
    report1(i);
    mat->rem_nrows--;
}

/* Compares two connected components. Since most components are trees,
   i.e., should have dif := #edges - #nodes = -1, we simply compare
   the number of edges. With qsort, the best component will be last.
*/
int
compare(const void *v1, const void *v2)
{
  component_t *c1 = (component_t*) v1;
  component_t *c2 = (component_t*) v2;

  return c1->edg - c2->edg;
}

/* Compute connected components. Assume mat->nodes and mat->edges have
   already been allocated, with mat->nrows entries each. */
void
compute_components(component_t **comps, int *ncomps, int *alloc, int *nodes, int *edges, sparse_mat_t *mat)
{
    component_t *lcomps = *comps;
    int i, lncomps = 0, lalloc = *alloc;
    
    for(i = 0; i < mat->nrows; i++)
	if(nodes[i] >= 0)
	    nodes[i] = 0;
    /* we use nodes[i] for visited[i], since nodes[i] >= 1 */
    for(i = 0; i < mat->nrows; i++)
	if(nodes[i] == 0){
	    visit(i, nodes, edges, mat);
	    // i is interesting if nodes[i] > 1
	    if(nodes[i] == 1)
		continue;
#if DEBUG >= 1
	    fprintf(stderr, "(nodes, edges)[%d] = (%d, %d)\n",
		    i, nodes[i], edges[i]);
#endif
	    ASSERT ((edges[i] & 1) == 0);
	    edges[i] >>= 1;
	    if(lncomps >= lalloc){
		lalloc += lalloc / 10 + 1;
		lcomps = (component_t*) realloc(lcomps,
						lalloc * sizeof(component_t));
	    }
	    lcomps[lncomps].ind = i;
#if 0
	    lcomps[lncomps].dif = edges[i] - nodes[i];
#endif
	    lcomps[lncomps].edg = edges[i];
	    lncomps++;
	}
    /* FIXME: since we only keep the largest components in prune, we do not
       need to sort, just to put the largest components at the end, which
       can be done in O(n) in a single pass:
       i = 0 # current component
       j = n-1 # last component, we assume size(c[j])=...=size(c[n-1])
       while (i < j)
          if size(c[i])<size(c[j])
             j = n - 1
             swap(c[i], c[j])
          elif size(c[i])=size(c[j])
             j = j - 1
             swap(c[i], c[j])
          i = i + 1
    */
    qsort(lcomps, lncomps, sizeof(component_t), compare);
#if DEBUG >= 0
    fprintf (stderr, "****** Found %d connected components\n", lncomps);
#endif
    *alloc = lalloc;
    *ncomps = lncomps;
    *comps = lcomps;
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
prune(sparse_mat_t *mat, int keep)
{
    int i, t, excess, old_nrows, old_ncols, mid;
    int *nodes, *edges, alloc, ncomps;
    component_t *c, *comps;
    int cur_nodes, old_nodes = 0;
    int count_recompute = 1;
    int iter = 0;
    
    nodes = (int*) malloc(mat->nrows * sizeof (int));
    edges = (int*) malloc(mat->nrows * sizeof (int));
    
    for(i = 0; i < mat->nrows; i++)
	nodes[i] = 0;
    
    alloc = 0;
    comps = NULL;
    ncomps = 0;
    
    compute_components(&comps, &ncomps, &alloc, nodes, edges, mat);
    excess = mat->nrows - mat->ncols;

    mid = (excess + keep) / 2;
    /* in theory, we should write excess + comps[ncomps - 1].dif >= keep,
       but in most cases dif is -1, thus we simply write excess > keep */
    while (ncomps > 0 && excess > keep)
      {
	t = ncomps - 1;
	c = comps + t;
	i = c->ind;
	/* Warning: this component might already have been deleted */
	if(nodes[i] >= 0){ /* it was not deleted */
	    old_nrows = mat->rem_nrows;
	    old_ncols = mat->rem_ncols;
	    delete(i, nodes, mat);
	    // humf!
            deleteAllColsFromStack(mat, 0);
            cur_nodes = old_nrows - mat->rem_nrows;
	    excess = mat->rem_nrows - mat->rem_ncols;
#if DEBUG >= 0
            if (iter++ < 10 || (iter & 2047) == 0)
              {
                fprintf(stderr,"Removed %d rels/%d primes, ",
                        cur_nodes, old_ncols - mat->rem_ncols);
                fprintf(stderr, "remains %d rels/%d primes (excess %d)\n",
                        mat->rem_nrows, mat->rem_ncols, excess);
              }
#endif
            if (old_nodes == 0)
              old_nodes = cur_nodes;
            /* New strategy: instead of recomputing the components when the
               excess is divided by two, recompute whenever the size of the
               removed components decreases. This should avoid we remove very
               small components at the first iteration. */
	    else if (cur_nodes < old_nodes || excess < mid)
              {
		compute_components(&comps, &ncomps, &alloc, nodes, edges, mat);
                count_recompute ++;
                old_nodes = 0;
                iter = 0;
                mid = (excess + keep) / 2;
            }
        }
	ncomps = t;
    }
    fprintf (stderr, "Number of connected components computations: %d\n",
             count_recompute);
}

