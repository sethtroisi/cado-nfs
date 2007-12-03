#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"

#include <assert.h>

#include "macros.h"
#include "types.h"
#include "endian.h"
#include "auxfuncs.h"
#include "relation.h"

unsigned char relation_degLP(type32 * d)
{
    type32 a,m;
    int i,s;

    if (d[1])
    {
	a=d[1];
	s=32;
    }
    else
    {
	a=d[0];
	s=0;
    }

    m=UC32(0xFF000000U);
    for(i=3;i>=0 && (a&m)==UC32(0);i--) m>>=8;
    s+=i*8;

    m=UC32(1<<((s+7)&31));
    for(i=7;i>=0 && (a&m)==UC32(0);i--) m>>=1;
    return s+i;
}

void _relation_alloc (relation * p_r, type32 lnginfo)
{
    int data_length=0;

    (*p_r)	 = my_malloc(sizeof(relation_s));
    (*p_r)->info = lnginfo;
    if (lnginfo==0)
	(*p_r)->data = NULL;
    else
    {
	data_length  = RELATION_DATASIZE(*p_r);
	(*p_r)->data = my_malloc(data_length*sizeof(type32));
    }
}

void _ext_relation_alloc (ext_relation * p_r, type32 lnginfo)
{
    int hascs;

    (*p_r)	 = my_malloc(sizeof(ext_relation_s));
    (*p_r)->info = lnginfo;
    (*p_r)->nSPs = RELATION_NSPS(*p_r);
    (*p_r)->nLPs = RELATION_NLPS(*p_r);
    hascs        = RELATION_HAS_CHECKSUM(*p_r);
    (*p_r)->data_length = (*p_r)->nSPs*2+(*p_r)->nLPs*3+hascs*2;

    if (lnginfo==0)
    {
	(*p_r)->data = NULL;
	(*p_r)->LP   = NULL;
	(*p_r)->SP   = NULL;
	(*p_r)->checksum   = NULL;
    }
    else
    {
	(*p_r)->data = my_malloc((*p_r)->data_length*sizeof(type32));
	(*p_r)->LP	 = RELATION_LP(*p_r);
	(*p_r)->SP	 = RELATION_SP(*p_r);
	(*p_r)->checksum = RELATION_CS(*p_r);
    }
    (*p_r)->degLP[0] = 0;
    (*p_r)->degLP[1] = 0;
}


void _relation_free (relation *p_r)
{
    if ((*p_r)->info)
	free((*p_r)->data);
    free(*p_r);
    *p_r=NULL;
}

void _ext_relation_free (ext_relation *p_r)
{
    if ((*p_r)->info)
	free((*p_r)->data);
    free(*p_r);
    *p_r=NULL;
}

int relation_read (relation r, FILE *f)
{
    if (r->info)
	free(r->data);
    fread(&(r->info),sizeof(type32),1,f);
    DO_BIG_ENDIAN(mswap32(r->info));

    if (r->info==0)
    {
	r->data     = NULL;
	return 1;
	/* return 0; */
    }

#if (SPMLINE_MAXIMUM_LENGTH_ON_INPUT!=0)
    if (RELATION_NSPS(r)>SPMLINE_MAXIMUM_LENGTH_ON_INPUT)
    {
	fprintf(stderr,"Ligne aberrante (trop longue), ignorée\n");
	r->data     = NULL;
	return 0;
    }
#endif

    r->data=my_malloc(RELATION_DATASIZE(r)*sizeof(type32));
    fread(r->data,sizeof(type32),RELATION_DATASIZE(r),f);
    DO_BIG_ENDIAN(
	unsigned int i;
	for(i=0;i<RELATION_DATASIZE(r);i++)
	    mswap32(r->data[i])
    );
    return 1;
}

int ext_relation_read(ext_relation r, FILE *f)
{
    int hascs,n;

    if (r->info)
	free(r->data);
    fread(&(r->info),sizeof(type32),1,f);
    DO_BIG_ENDIAN(mswap32(r->info));

    r->nLPs        = RELATION_NLPS(r);
    r->nSPs        = RELATION_NSPS(r);
    hascs          = RELATION_HAS_CHECKSUM(r);
    r->data_length = r->nSPs*2+r->nLPs*3+hascs*2;

#if (SPMLINE_MAXIMUM_LENGTH_ON_INPUT!=0)
    if (r->nSPs>SPMLINE_MAXIMUM_LENGTH_ON_INPUT)
    {
	fprintf(stderr,"Ligne aberrante (trop longue), ignorée\n");
	goto return_null;
    }
#endif
    
    if (r->info==0)
    {
    r->data     = NULL;
    r->LP       = NULL;
    r->SP       = NULL;
    r->checksum = NULL;
	return 1;
	/* return 0; */
    }

    r->data=my_malloc(r->data_length*sizeof(type32));
    n=fread(r->data,sizeof(type32),r->data_length,f);
    if (((unsigned int) n) < r->data_length)
    {
	fprintf(stderr,"Ligne aberrante (trop longue), ignorée\n");
	free(r->data);
	goto return_null;
    }

    DO_BIG_ENDIAN(
	unsigned int i;
	for(i=0;i<RELATION_DATASIZE(r);i++)
	    mswap32(r->data[i])
    );

    r->LP	= RELATION_LP(r);
    r->SP	= RELATION_SP(r);
    r->checksum = RELATION_CS(r);
    r->degLP[0] = (r->nLPs>0)?relation_degLP(r->LP[0].poly):0;
    r->degLP[1] = (r->nLPs>1)?relation_degLP(r->LP[1].poly):0;

    return 1;

return_null:
    r->data     = NULL;
    r->LP       = NULL;
    r->SP       = NULL;
    r->checksum = NULL;

    return 0;
}

int relation_write (FILE *f, relation r)
{
    unsigned int n=RELATION_DATASIZE(r);

    DO_BIG_ENDIAN(mswap32(r->info));
    fwrite(&(r->info),sizeof(type32),1,f);
    DO_BIG_ENDIAN(mswap32(r->info));

    if (r->data==NULL) return 0;

    DO_BIG_ENDIAN(unsigned int i; for(i=0;i<n;i++) mswap32(r->data[i]));
    fwrite(r->data,sizeof(type32),n,f);
    DO_BIG_ENDIAN(unsigned int i; for(i=0;i<n;i++) mswap32(r->data[i]));

    return 0;
}

int ext_relation_write (FILE *f, ext_relation r)
{
    DO_BIG_ENDIAN(mswap32(r->info));
    fwrite(&(r->info),sizeof(type32),1,f);
    DO_BIG_ENDIAN(mswap32(r->info));

    if (r->data==NULL) return 0;

    DO_BIG_ENDIAN(unsigned int i; for(i=0;i<r->data_length;i++) mswap32(r->data[i]));
    fwrite(r->data,sizeof(type32),r->data_length,f);
    DO_BIG_ENDIAN(unsigned int i; for(i=0;i<r->data_length;i++) mswap32(r->data[i]));

    return 0;
}

#ifdef CHECKSUM_ENABLED
void relation_checksum (type32 * dest, const type32 * desc, relation r)
{
    compute_checksum(dest,desc,(void*)
	    (r->data+(RELATION_HAS_CHECKSUM(r)?2:0)),
	    ((RELATION_NLPS(r)*3)+(RELATION_NSPS(r)<<1)),sizeof(type32));
}

void ext_relation_checksum (type32 * dest, const type32 * desc, ext_relation r)
{
    compute_checksum(dest,desc,(void*)(r->data+(r->checksum?2:0)),
	    (r->data_length-(r->checksum?2:0)),sizeof(type32));
}
#endif

void _relation_extend (ext_relation * p_s, relation * p_r)
{
    ext_relation n;

    n=my_malloc(sizeof(ext_relation_s));
    memcpy(n,*p_r,sizeof(relation_s));

    n->nSPs	   = RELATION_NSPS(*p_r);
    n->nLPs	   = RELATION_NLPS(*p_r);
    n->SP	   = RELATION_SP(*p_r);
    n->LP	   = RELATION_LP(*p_r);
    n->checksum	   = RELATION_CS(*p_r);
    n->data_length = RELATION_DATASIZE(*p_r);
    n->degLP[0]	   = (n->nLPs>0)?relation_degLP(n->LP[0].poly):0;
    n->degLP[1]	   = (n->nLPs>1)?relation_degLP(n->LP[1].poly):0;

    free(*p_r);
    *p_r=NULL;
    *p_s=n;
}

void _ext_relation_dup(ext_relation * p_s, ext_relation r)
{
    *p_s	 = my_malloc(sizeof(ext_relation_s));
    memcpy(*p_s,r,sizeof(ext_relation_s));
    (*p_s)->data = my_malloc(r->data_length*sizeof(type32));
    memcpy((*p_s)->data,r->data,r->data_length*sizeof(type32));
    (*p_s)->SP	   = RELATION_SP(*p_s);
    (*p_s)->LP	   = RELATION_LP(*p_s);
    (*p_s)->checksum = RELATION_CS(*p_s);
}

void ext_relation_reap_checksum(ext_relation s)
{
    type32 * alt;

    if (s->checksum==NULL)
	return;

    s->data_length -= 2;

    alt	= my_malloc(s->data_length*sizeof(type32));
    memcpy(alt,s->data+2,s->data_length*sizeof(type32));
    free(s->data);
    s->data	 = alt;
    s->checksum	 = NULL;
    if (s->LP!=NULL)
	s->LP	-= 2;
    s->SP	-= 2;

    s->info	 = s->nSPs | (s->checksum?RELATION_W_CHECKSUM:0) |
	(s->nLPs==1?RELATION_PF_FP:(s->nLPs==2?RELATION_PP:RELATION_FF));
}


void ext_relation_adjust(ext_relation r)
{
    /* In the case where r->nSPs might have been modified */
    unsigned int data_length;

    data_length	= r->nSPs*2+r->nLPs*3+(r->checksum?2:0);
    if (data_length == r->data_length)
	return;

    r->data        = realloc(r->data,data_length*sizeof(type32));
    r->data_length = data_length;

    r->info	 = r->nSPs | (r->checksum?RELATION_W_CHECKSUM:0) |
	(r->nLPs==1?RELATION_PF_FP:(r->nLPs==2?RELATION_PP:RELATION_FF));

    r->SP	= RELATION_SP(r);
    r->LP	= RELATION_LP(r);
    r->checksum	= RELATION_CS(r);
}

int sp_compare(sp_data a, sp_data b)
{
    return a->index-b->index;
}
