/* 
 * Program: history
 * Author : F. Morain
 * Purpose: managing history of merges
 * 
 * Algorithm:
 *
 */


#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "portability.h"
#include "utils.h"
#include "macros.h"
#include "filter_utils.h"
#include "filter_matrix.h"
#include "sparse.h"

#define DEBUG 0

void
fprintRow(FILE *file, typerow_t *row)
{
  index_t i;
#ifdef FOR_DL
    fprintf(file, "[%" PRid "]", row[0].id);
    for(i = 1; i <= row[0].id; i++)
        fprintf(file, " %" PRid "(%" PRId32 ")", row[i].id, row[i].e);
#else
    fprintf(file, "[%" PRid "]", row[0]);
    for(i = 1; i <= row[0]; i++)
        fprintf(file, " %" PRid "", row[i]);
#endif
}

// i1 += i2
// A row is row[0..max] where row[0] = max and the real components are
// row[1..max].
// If len != -1, then it is the real length of row[i1]+row[i2].
//
// If j is given, it is the index of the column that is used for
// pivoting in the case of DL. Then, the operation is
//   i1 = e2*i1 + e1*i2
// where e1 and e2 are adjusted so that the j-th column is zero in i1.
//
// Also update the data for the index file, if needed (i.e. if the given
// pointer is not NULL).
void
addRowsUpdateIndex(typerow_t **rows, index_data_t index_data, int i1, int i2,
        MAYBE_UNUSED int32_t j)
{
    uint32_t k1, k2, k, len;
    typerow_t *tmp;

    ASSERT(rows[i1] != NULL);
    ASSERT(rows[i2] != NULL);
#if DEBUG >= 1
    fprintf(stderr, "R[%d] =", i1); fprintRow(stderr, rows[i1]); 
    fprintf(stderr, "\n");
    fprintf(stderr, "R[%d] =", i2); fprintRow(stderr, rows[i2]);
    fprintf(stderr, "\n");
#endif
    len = 1 + rowLength(rows, i1) + rowLength(rows, i2);
    tmp = (typerow_t *)malloc(len * sizeof(typerow_t));
    k = k1 = k2 = 1;

    int e1 = 1, e2 = 1;  // default value for non-DL

#ifdef FOR_DL /* look for the exponents of j in i1 i2*/
    e1 = 0, e2 = 0;
    int d;
    unsigned int l;
    for (l = 1 ; l <= rowLength(rows, i1) ; l++)
        if ((int) rowCell(rows, i1, l) == j)
            e1 = rows[i1][l].e;
    for (l = 1 ; l <= rowLength(rows, i2) ; l++)
        if ((int) rowCell(rows, i2, l) == j)
            e2 = rows[i2][l].e;

    ASSERT (e1 != 0 && e2 != 0);

    d  = (int) gcd_int64 ((int64_t) e1, (int64_t) e2);
    e1 /= -d;
    e2 /= d;

    for (l = 1 ; l <= rowLength(rows, i2) ; l++)
        rows[i2][l].e *= e1;
    for (l = 1 ; l <= rowLength(rows, i1) ; l++)
        rows[i1][l].e *= e2;

#if DEBUG >= 1
    fprintf (stdout, "Computing %d*rows[%d] + %d*rows[%d] for j=%d\n", 
                      e2, i1, e1, i2, j);
#endif

#endif


    // loop while everybody is here
    while((k1 <= rowLength(rows, i1)) && (k2 <= rowLength(rows, i2))){
        if(rowCell(rows, i1, k1) < rowCell(rows, i2, k2))
            tmp[k++] = rows[i1][k1++];
        else if(rowCell(rows, i1, k1) > rowCell(rows, i2, k2))
            tmp[k++] = rows[i2][k2++];
        else{
#ifdef FOR_DL
            if (rows[i1][k1].e + rows[i2][k2].e != 0)
            {
                tmp[k].id = rows[i1][k1].id;
                tmp[k++].e = rows[i1][k1].e + rows[i2][k2].e;
            }
#else
#if DEBUG >= 1
            fprintf(stderr, "WARNING: j1=j2=%d in addRows\n", k1);
#endif
#endif
            k1++;
            k2++;
        }
    }
    // finish with k1
    for( ; k1 <= rowLength(rows, i1); k1++)
	tmp[k++] = rows[i1][k1];
    // finish with k2
    for( ; k2 <= rowLength(rows, i2); k2++)
	tmp[k++] = rows[i2][k2];
    ASSERT(k <= len);

    // copy back
    free(rows[i1]);
#ifdef FOR_DL
        tmp[0].id = k-1;
#else
        tmp[0] = k-1;
#endif
	if(k == len)
	    rows[i1] = tmp;
	else
	    rows[i1] = realloc(tmp, k * sizeof(typerow_t));
#ifdef FOR_DL
    /* restore old coeff for row i2 */
    for (l = 1 ; l <= rowLength(rows, i2) ; l++)
        rows[i2][l].e /= e1;
#endif


    // Now, deal with the index_data.
    if (index_data != NULL) {
        k = k1 = k2 = 0;   // in index_data_t, we count from 0...

        relset_t r1 = index_data[i1];
        relset_t r2 = index_data[i2];
        relset_t tmp;
        tmp.n = 0;
        tmp.rels = (multirel_t *) malloc((r1.n+r2.n)*sizeof(multirel_t));
        while ((k1 < r1.n) && (k2 < r2.n)) {
            if (r1.rels[k1].ind_row < r2.rels[k2].ind_row) {
                tmp.rels[k].ind_row = r1.rels[k1].ind_row;
                tmp.rels[k++].e = e2*r1.rels[k1++].e;
            } else if (r1.rels[k1].ind_row > r2.rels[k2].ind_row) { 
                tmp.rels[k].ind_row = r2.rels[k2].ind_row;
                tmp.rels[k++].e = e1*r2.rels[k2++].e;
            } else {
#ifdef FOR_DL
                int32_t e = e2*r1.rels[k1].e + e1*r2.rels[k2].e;
                if (e != 0) {
                    tmp.rels[k].ind_row = r1.rels[k1].ind_row;
                    tmp.rels[k++].e = e;
                }
#endif
            k1++;
            k2++;
            }
        }
        // finish with k1 and k2
        for( ; k1 < r1.n; k1++) {
            tmp.rels[k].ind_row = r1.rels[k1].ind_row;
            tmp.rels[k++].e = e2*r1.rels[k1].e;
        }
        for( ; k2 < r2.n; k2++) {
            tmp.rels[k].ind_row = r2.rels[k2].ind_row;
            tmp.rels[k++].e = e1*r2.rels[k2].e;
        }
        ASSERT (k <= r1.n + r2.n);
        tmp.n = k;

        // copy back to i1
        free (index_data[i1].rels);
        index_data[i1] = tmp;
    }

#if DEBUG >= 1
    fprintf(stderr, "row[%d]+row[%d] =", i1, i2);
    fprintRow(stderr, rows[i1]); fprintf(stderr, "\n");
#endif
}

int
hasCol(int32_t **rows, int i, int32_t j)
{
    int32_t k;

    for(k = 1; k <= rows[i][0]; k++)
	if(rows[i][k] == j)
	    return 1;
    return 0;
}

// A line is "[-]i i1 ... ik [#j]"
int parse_hisfile_line (int32_t *ind, char *t, int32_t *j)
{
  int ni = 0, sg = 1;

  if(*t == '-')
    {
      sg = -1;
	    t++;
    }

  ind[0] = 0;
  while(1)
    {
      if((*t == '\n') || (*t == '#'))
          break;
      else if (*t == ' ')
        {
          ni++;
          ind[ni] = 0;
        }
	    else
          ind[ni] = 10 * ind[ni] + (*t - '0');

	    t++;
    }

  ind[0] = sg * ind[0];
  *j=0;

  if (*t == '#')
      for (t++ ; (*t != '\n') ; t++)
          *j = 10 * (*j) + (*t - '0');
  else
    {
      ni++;
      *j = -1;
    }

  return ni;
}
