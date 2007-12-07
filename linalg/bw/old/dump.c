#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "params.h"
#include <assert.h>
#include "types.h"
#include "macros.h"
#include "old-endian.h"
#include "auxfuncs.h"

const char * indexes_fmt="run/output/%d/indexes";
const char * values_fmt="run/output/%d/values";

#define FILENAME_LENGTH 80

struct line {
	unsigned int weight;
	size_t	idx_off;
	size_t	val_off;
	size_t	idx_size;
	size_t	val_size;
};

type32 * indexes;
stype32 * values;
type32 nlines;

int read_line(struct line * l, stype32 * p_idx, stype32 * p_val)
{
	stype32 idx,val;
	idx=*p_idx;
	val=*p_val;
	l->weight=0;
	l->idx_off=idx;
	l->val_off=val;
	l->idx_size=0;
	l->val_size=0;
	if (indexes[idx]==0 && values[val]==0) {
		val++;
		l->val_size++;
	} else {
		do {
			idx++;
			l->idx_size++;
			val++;
			l->val_size++;
			l->weight++;
		} while (indexes[idx]!=0);
	}
	idx++;
	l->idx_size++;
	*p_idx=idx;
	*p_val=val;
	return l->weight;
}

void display_line(struct line * l)
{
	int i,k=0,i0;
	for(i=i0=0;i<nlines && k<l->weight;i++) {
		if (i-i0==indexes[l->idx_off+k]) {
			printf("%d, ", values[l->val_off+k]);
			k++;
			i0=i;
		} else {
			printf("0, ");
		}
	}
	for(;i<nlines;i++) {
		printf("0, ");
	}
	printf("\n");
}

void display_line_sparse(struct line * l)
{
	int k,i0=0;
	printf("%d", l->weight);
	for(k=0;k<l->weight;k++) {
		i0+=indexes[l->idx_off+k];
		printf(" %d:%d", i0, values[l->val_off+k]);
	}
	printf("\n");
}

struct line * linetab = NULL;
type32 alloc_lines=0;

struct line * gimme_lineptr(type32 i)
{
	if (i >= alloc_lines) {
		alloc_lines=(alloc_lines!=0)?(2*alloc_lines):(alloc_lines+1);
		linetab=realloc(linetab,alloc_lines*sizeof(struct line));
	}
	return linetab + i;
}

int main(int argc, char *argv[])
{
	char fn_indexes[FILENAME_LENGTH];
	char fn_values[FILENAME_LENGTH];
	int fd;
	struct stat sbuf;
	type32 i;
	type32 nc;
	int r;
	stype32 idx,val;
	size_t isize;
	size_t vsize;
	int sparse=0;
	int quiet=0;

	for(;argc>2;) {
		if (strcmp(argv[1],"--sparse")==0) {
			argc--; argv++; sparse=1; continue;
		}
		if (strcmp(argv[1],"-q")==0) {
			argc--; argv++; quiet=1; continue;
		}
		break;
	}

	if (argc!=2)
		die("Please give a bank number\n",1);

	sprintf(fn_indexes,indexes_fmt,atoi(argv[1]));
	sprintf(fn_values,values_fmt,atoi(argv[1]));

	fd=open(fn_indexes,O_RDONLY);
	if (fd<0) {perror("open(indexes)"); exit(errno);}
	if (fstat(fd,&sbuf)<0) {perror("fstat(indexes)"); exit(errno);}
	isize=sbuf.st_size;
	if (!quiet) printf("indexes : %lu bytes\n",(unsigned long)isize);
	indexes=(void*)mmap(NULL,isize,PROT_READ,MAP_SHARED,fd,0);
	if (indexes==NULL) {perror("mmap(indexes)"); exit(errno);}
	if (close(fd)<0) {perror("close(indexes)"); exit(errno);}

	fd=open(fn_values,O_RDONLY);
	if (fd<0) {perror("open(values)"); exit(errno);}
	if (fstat(fd,&sbuf)<0) {perror("fstat(values)"); exit(errno);}
	vsize=sbuf.st_size;
	if (!quiet) printf("values : %lu bytes\n",(unsigned long)vsize);
		values=(void*)mmap(NULL,vsize,PROT_READ,MAP_SHARED,fd,0);
	if (values==NULL) {perror("mmap(values)"); exit(errno);}
	if (close(fd)<0) {perror("close(values)"); exit(errno);}

	nc=0;
	idx=val=0;
	for(i=0;;i++) {
		r=read_line(gimme_lineptr(i),&idx,&val);
		nc+=r;
		if (	idx * sizeof(type32) == isize &&
			val * sizeof(type32) == vsize)
		{
			i++;
			break;
		}
		if (r==0) {
			if (!quiet)
				printf("Blank line before end of file\n");
			continue;
		}
	}
	nlines=i;

	if (!quiet)
		printf("Read %lu lines, totalizing %lu coeffs\n",
			(unsigned long)nlines,
			(unsigned long)nc);

	if (idx * sizeof(type32) != isize || val * sizeof(type32) != vsize) {
		printf("File size not reached.\n");
		printf("Stopped at offset %lu(indexes):%lu(values)\n",
				(unsigned long)idx * sizeof(type32),
				(unsigned long)val * sizeof(type32));
	}

	if (sparse) {
		for(i=0;i<nlines;i++) display_line_sparse(&(linetab[i]));
	} else {
		for(i=0;i<nlines;i++) display_line(&(linetab[i]));
	}

	return 0;
}



