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
#include "macros.h"
#include "filter_matrix.h"
#include "sparse.h"
#include "utils.h"

#define DEBUG 0

void
fprintRow(FILE *file, typerow_t *row)
{
    int i;
#ifdef FOR_FFS
    fprintf(file, "[%d]", row[0].id);
    for(i = 1; i <= row[0].id; i++)
        fprintf(file, " %d(%d)", row[i].id, row[i].e);
#else
    fprintf(file, "[%d]", row[0]);
    for(i = 1; i <= row[0]; i++)
        fprintf(file, " %d", row[i]);
#endif
}

// row[0..row[0]] is of lenfth row[0]+1.
int32_t*
copyRow(int32_t *row)
{
    int32_t *tmp = (int32_t *)malloc((1+row[0]) * sizeof(int32_t));
    
    memcpy(tmp, row, (1+row[0]) * sizeof(int32_t));
    return tmp;
}

// i1 += i2
// A row is row[0..max] where row[0] = max and the real components are
// row[1..max].
// If len != -1, then it is the real length of row[i1]+row[i2].
void
addRows(typerow_t **rows, int i1, int i2, MAYBE_UNUSED int32_t j)
{
    int32_t k1, k2, k, len;
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

#ifdef FOR_FFS /* look for the exponents of j in i1 i2*/
    int e1 = 0, e2 = 0;
    int d;
    int l;
    for (l = 1 ; l <= rowLength(rows, i1) ; l++)
        if (rowCell(rows, i1, l) == j)
            e1 = rows[i1][l].e;
    for (l = 1 ; l <= rowLength(rows, i2) ; l++)
        if (rowCell(rows, i2, l) == j)
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
#ifdef FOR_FFS
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
        /* FIXME: why not use realloc here instead? Since k <= len,
           it suffices to shrink the tmp[] array to k entries.
           Also, we might detect the special case k = len. */
#if 0
	int *tmp2 = (int32_t *)malloc(k * sizeof(int32_t));
	memcpy(tmp2, tmp, k * sizeof(int32_t));
	tmp2[0] = k-1;
	rows[i1] = tmp2;
	free(tmp);
#else
#ifdef FOR_FFS
        tmp[0].id = k-1;
#else
        tmp[0] = k-1;
#endif
	if(k == len)
	    rows[i1] = tmp;
	else
	    rows[i1] = realloc(tmp, k * sizeof(typerow_t));
#ifdef FOR_FFS
    /* restore old coeff for row i2 */
    for (l = 1 ; l <= rowLength(rows, i2) ; l++)
        rows[i2][l].e /= e1;
#endif
#endif
#if DEBUG >= 1
    fprintf(stderr, "row[%d]+row[%d] =", i1, i2);
    fprintRow(stderr, rows[i1]); fprintf(stderr, "\n");
#endif
}

void
removeWeight(int32_t **rows, int *wt, int i)
{
    int32_t k;

    for(k = 1; k <= rows[i][0]; k++)
	wt[rows[i][k]]--;
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

int 
cmp (const void *p, const void *q)
{
  int x = *((int *)p);
  int y = *((int *)q);
  return (x <= y ? -1 : 1);
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
