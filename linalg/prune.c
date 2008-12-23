/* 
 * Program: pruning from matprune.c converted to the use of merge.c by FM
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <limits.h> /* for INT_MAX */

#define DEBUG 0

#include "utils/utils.h"

#include "merge_opts.h"
#include "sparse.h"
#include "dclist.h"
#include "sparse_mat.h"
#include "report.h"
#include "swar.h"
#include "merge_mono.h"
#include "prune.h"

#if USE_TAB == 0
#error "This file assumes USE_TAB is non-zero"
#endif

double removeCellSWAR_time = 0.0;

typedef struct {
  int ind;        /* row index of starting of component */
#if 0
  int dif;        /* edges[i] - nodes[i] >= -1 */
#endif
  int edg;        /* edges[i] */
} component_t;

/* visit connected component of relation i, and compute:
   nodes[i] : the number of nodes (i.e., rows or relations) of the connected
              component of i
   edges[i] : twice the numbers of edges of that component */
void
visit (int i, int *nodes, int *edges, sparse_mat_t *mat, unsigned int *rsum)
{
    int m, j, k;
    
    nodes[i] = 1; /* also mark as visited */
    edges[i] = 0;
    for(m = 1; m <= lengthRow(mat, i); m++){
	j = cell(mat, i, m);
	if(mat->wt[GETJ(mat, j)] == 2){
	    edges[i]++;
            k = rsum[j] - i;
	    if(nodes[k] == 0){
                visit(k, nodes, edges, mat, rsum);
		nodes[i] += nodes[k];
		edges[i] += edges[k];
            }
        }
    }
}

/* Delete relations and primes in component of relation/row i.
   FIXME: if some wt[j] becomes 2, for some prime/column j, we should join
   the connected components of the two remaining relations containing j.
   RETURN VALUE: number of deletions performed, hence the number of j
   removed (i.e. which have become of weight 0).
*/
int
delete (report_t *rep, int i, int *nodes, sparse_mat_t *mat, unsigned int *rsum)
{
    INT j;
    int k, v, nd = 0;
    int l = lengthRow(mat, i); /* number of non-zero entries in row i */
    
    nodes[i] = -1; /* mark row i as deleted and visited */
    // inspired from removeRowSWAR(mat, i) but for the recursive call
    // to delete...!
    mat->weight -= l; /* decrease the total weight of the matrix */
    for(k = 1; k <= l; k++){
        j = cell(mat, i, k); /* column index of the kth non-zero coefficient */
        rsum[j] -= i;        /* update rsum */
        removeCellSWAR_time -= seconds ();
	removeCellAndUpdate(mat, i, j); /* defined in swar.c */
        removeCellSWAR_time += seconds ();
	if(mat->wt[GETJ(mat, j)] == 0)
	    nd++;
	if(mat->wt[GETJ(mat, j)] == 1){
            v = rsum[j]; /* other row containing j */
	    if (nodes[v] >= 0)
		nd += delete (rep, v, nodes, mat, rsum);
	}
    }
    destroyRow(mat, i);
    report1(rep, i);
    mat->rem_nrows--;
    return nd;
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
compute_components (component_t **comps, int *ncomps, int *alloc, int *nodes,
                    int *edges, sparse_mat_t *mat, unsigned int *rsum)
{
    component_t *lcomps = *comps;
    int i, lncomps = 0, lalloc = *alloc;
    
    for(i = 0; i < mat->nrows; i++)
      /* if nodes[i]=-1, row i was deleted, don't consider it any more */
	if(nodes[i] >= 0)
	    nodes[i] = 0;
    /* we use nodes[i] for visited[i], since nodes[i] >= 1 */
    for(i = 0; i < mat->nrows; i++)
      if (nodes[i] == 0){
            visit (i, nodes, edges, mat, rsum);
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
#if DEBUG >= 1
    fprintf (stderr, "****** Found %d connected components\n", lncomps);
#endif
    *alloc = lalloc;
    *ncomps = lncomps;
    *comps = lcomps;
}

/* We consider the graph whose vertices are rows of the matrix
   (i.e., relations), where there is an edge between vertex i and vertex k
   if there is a common prime of index j in relations i and k, which appears
   nowhere else, i.e., the prime/column j is of weight 2.
   (This is a multigraph, since several primes can appear only in i and k.)

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
prune (report_t *rep, sparse_mat_t *mat, int keep)
{
    int i, j, k, t, excess, old_nrows, old_ncols, mid;
    int *nodes, *edges, alloc, ncomps;
    component_t *c, *comps;
    int count_recompute = 1;
    double prune_time, component_time = 0.0, rsum_time;
    unsigned int *rsum; /* rsum[j] is the sum (modulo 2^32 or 2^64) of rows
                           where prime of index j appears. When only two
                           rows remain, and one is known, the other one can
                           be recovered by rsum[j] - i, as long as all row
                           indices are < 2^32 or 2^64. */
    int cur_nodes; /* nb of nodes of the current (deleted) component */
    int deleted_component = 0;
    int min_size = INT_MAX, max_size = 0;

    prune_time = seconds ();
    
    fprintf(stderr, "Entering prune with keep=%d\n", keep);
    /* nodes[i] is the number of nodes (i.e., rows or relations) in the
       connected component of row i in the above-defined graph */
    nodes = (int*) malloc(mat->nrows * sizeof (int));
    /* edges[i] is twice the number of edges in the connected component of i */
    edges = (int*) malloc(mat->nrows * sizeof (int));

    rsum_time = seconds ();
    rsum = (uint32_t*) malloc(mat->ncols * sizeof (uint32_t));
    for (j = 0; j < mat->ncols; j++)
      {
        rsum[j] = 0;
        if (mat->R[GETJ(mat, j)] != NULL)
          for (k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
            /* no R[j][k] should be -1 here, since we just read the matrix */
            rsum[j] += mat->R[GETJ(mat, j)][k];
      }
    rsum_time = seconds () - rsum_time;
    fprintf (stderr, "Computing rsum took %1.2e seconds\n", rsum_time);
    
    for(i = 0; i < mat->nrows; i++)
      nodes[i] = 0; /* mark all rows as non-visited */
    
    alloc = 0;
    comps = NULL;
    ncomps = 0;

    component_time -= seconds ();
    compute_components (&comps, &ncomps, &alloc, nodes, edges, mat, rsum);
    component_time += seconds ();
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
	if (nodes[i] >= 0){ /* it was not deleted */
	    old_nrows = mat->rem_nrows;
	    old_ncols = mat->rem_ncols;
	    mat->rem_ncols -= delete (rep, i, nodes, mat, rsum);
            deleted_component ++;
            cur_nodes = old_nrows - mat->rem_nrows;
            if (cur_nodes < min_size)
              min_size = cur_nodes;
            if (cur_nodes > max_size)
              max_size = cur_nodes;
	    excess = mat->rem_nrows - mat->rem_ncols;
#if DEBUG >= 1
            fprintf(stderr,"Removed %d rels/%d primes, ",
                    cur_nodes, old_ncols - mat->rem_ncols);
            fprintf(stderr, "remains %d rels/%d primes (excess %d)\n",
                    mat->rem_nrows, mat->rem_ncols, excess);
#endif
            if (excess < mid)
              {
                fprintf (stderr, "%d-%d:%d ", max_size, min_size,
                         deleted_component);
                min_size = INT_MAX;
                max_size = 0;
                deleted_component = 0;
                component_time -= seconds ();
                compute_components (&comps, &ncomps, &alloc, nodes, edges,
                                    mat, rsum);
                component_time += seconds ();
                count_recompute ++;
                mid = (excess + keep) / 2;
              }
        }
	else
          ncomps = t;
      }
    fprintf (stderr, "\n");
#if 1 // old code replace by this...
    // humf!
    t = mat->rem_ncols;
    deleteEmptyColumns(mat);
    mat->rem_ncols = t;
#endif
    fprintf (stderr, "Number of connected components computations: %d\n",
             count_recompute);
    free(nodes);
    free(edges);
    free (rsum);

    prune_time = seconds () - prune_time;
    fprintf (stderr, "Prune time=%1.2e (components %1.2e, removeCellSWAR %1.2e)\n", prune_time, component_time, removeCellSWAR_time);
}

