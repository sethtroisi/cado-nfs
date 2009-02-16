#include "utils/utils.h"
#include "sparse.h"
#include "dclist.h"

dclist
dclistCreate(int32_t j)
{
    dclist dcl = (dclist) malloc(sizeof(struct dclist));

    dcl->j = j;
    dcl->prev = NULL;
    dcl->next = NULL;
    return dcl;
}

void
dclistClear (dclist dcl)
{
  if (dcl != NULL)
    {
      dclistClear (dcl->next);
      free (dcl);
    }
}

void
dclistPrint(FILE *file, dclist dcl)
{
    while(dcl != NULL){
      fprintf(file, " -> %ld", (long int) dcl->j);
	dcl = dcl->next;
    }
}

int
dclistLength(dclist dcl)
{
    int l = 0;

    while(dcl != NULL){
	l++;
	dcl = dcl->next;
    }
    return l;
}

/* insert j in doubly-chained list dcl (between cell of dcl and dcl->next),
   and return pointer at cell containing j */
dclist
dclistInsert(dclist dcl, int32_t j)
{
    dclist newdcl = dclistCreate(j);

    newdcl->next = dcl->next;
    newdcl->prev = dcl;
    if(newdcl->next != NULL)
	newdcl->next->prev = newdcl;
    dcl->next = newdcl;
    return newdcl;
}

// connect dcl to S
void
dclistConnect(dclist S, dclist dcl)
{
    dcl->next = S->next;
    dcl->prev = S;
    S->next = dcl;
}

