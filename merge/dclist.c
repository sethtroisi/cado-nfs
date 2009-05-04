/* dclist.c --- routines dealing with doubly-chained lists

Copyright 2008-2009 Francois Morain.
Reviewed by Paul Zimmermann, February 2009.

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* Doubly-chained lists are used to store all the columns of same weight j;
   typically S[w] points to the doubly-chained list for weight w,
   and A[j] points to the unique cell in some S[w] containing j.
   We use two pointers to be able to remove a cell in the middle of the list.

   For technical reasons, the first element of the list is always -1.
   As a consequence, a doubly-chained list is never empty (NULL).
*/

#include "utils.h" /* for ASSERT */
#include "dclist.h"

/* initialize a doubly-chained list with element j */
dclist
dclistCreate (int32_t j)
{
    dclist dcl = (dclist) malloc (sizeof (struct dclist));

    dcl->j = j;
    dcl->prev = NULL;
    dcl->next = NULL;
    return dcl;
}

/* clear a doubly-chained list */
void
dclistClear (dclist dcl)
{
  ASSERT (dcl != NULL);
  if (dcl->next != NULL)
    dclistClear (dcl->next);
  free (dcl);
}

/* output a doubly-chained list */
void
dclistPrint (FILE *file, dclist dcl)
{
  ASSERT (dcl != NULL);
  fprintf (file, "%ld", (long int) dcl->j);
  dcl = dcl->next;
  while (dcl != NULL)
    {
      fprintf (file, " -> %ld", (long int) dcl->j);
      dcl = dcl->next;
    }
}

/* return the number of elements of a doubly-chained list, excluding the
   initial -1 */
int
dclistLength (dclist dcl)
{
    int l = -1;

    ASSERT (dcl != NULL);
    do
      {
	l++;
	dcl = dcl->next;
      }
    while (dcl != NULL);
    return l;
}

/* insert j in doubly-chained list dcl (between cell of dcl and dcl->next),
   and return pointer at cell containing j */
dclist
dclistInsert (dclist dcl, int32_t j)
{
    dclist newdcl;

    ASSERT (dcl != NULL);

    newdcl = dclistCreate (j);
    newdcl->next = dcl->next;
    newdcl->prev = dcl;
    if (dcl->next != NULL)
      dcl->next->prev = newdcl;
    dcl->next = newdcl;
    return newdcl;
}

/* remove current cell from doubly-chained list dcl, reconnect prev
   and next, and free dcl; the new list is not empty, since it remains
   at least the -1 entry
*/
void
dclistRemove (dclist dcl)
{
  dclist p, n;

  ASSERT (dcl != NULL);
  p = dcl->prev;
  if (p != NULL)
    p->next = dcl->next;
  n = dcl->next;
  if (n != NULL)
    n->prev = dcl->prev;
  free (dcl);
}

/* same as dclistRemove, but does not reconstruct prev pointers,
   and returns the deleted value of j */
int
dclistRemoveNext(dclist dcl)
{
  dclist foo;
  int j;

  ASSERT (dcl != NULL);
  foo = dcl->prev;
  j = dcl->j;
  if (foo != NULL)
    foo->next = dcl->next;
  free (dcl);
  return j;
}

/* return the first cell of dcl (avoiding the initial cell -1) */
dclist
dclistFirst (dclist dcl)
{
  ASSERT (dcl != NULL);
  return (dcl->j == -1) ? dcl->next : dcl;
}
