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
#include "endian.h"
#include "auxfuncs.h"
#include "variables.h"
#include "filenames.h"
#include "options.h"
#include "tagfile.h"

const char * indexes_fmt="run/output/%d/indexes";
const char * values_fmt="run/output/%d/values";
const char * mvec_fmt="run/output/%d/M-vector";
const char * x0_fmt="run/output/%d/X0-vector";
const char * lambda_fmt="run/output/%d/eigenvalue";

#define FILENAME_LENGTH 80

#define DO_CORRECTION
#define MAX_CORRECTION	255


stype32 * eigenvec;
stype32 * mvec;


type32 * indexes;
type32 * values;
size_t idx,val;

int careful_add(stype32 * dst, stype32 v)
{
	stype64 t;
	t=(stype64)*dst;
	t+=(stype64)v;
	if (ABS(t) & ~((UC64(1)<<31)-UC64(1))) {
		fprintf(stderr,"These two numbers are too big\n");
		exit(99);
	}
	*dst=(stype32)t;
	return (t != 0);
}

int register_coeff(type32 j, stype32 v)
{
	return careful_add(&mvec[j],v);
}

/* works only in /32/ mode */
int read_line_register(stype32 u)
{
	type32 j=0;
	int reset=0;
	type32	t;

	if (indexes[idx]==0 && values[val]==0) {
		val++;
	} else {
		do {
			t=indexes[idx];
			DO_BIG_ENDIAN(mswap32(t));
			j+=t;
			t=values[val];
			DO_BIG_ENDIAN(mswap32(t));
			reset|=(register_coeff(j, u*t) == 0);
			idx++;
			val++;
		} while (indexes[idx]!=0);
	}
	idx++;
	return reset;
}

stype32 read_line_nocancel(int i, stype32 u)
{
#ifdef DO_CORRECTION
	int reset;
	int try;
	size_t i_base;
	size_t v_base;

	i_base=idx;
	v_base=val;
	for(try=0;;try++) {
		idx=i_base;
		val=v_base;
		if (read_line_register(u)==0)
			break;
		if (try==0) {
			printf("LINE %d : %d", i, u);
		}

		/* Step back !  */
		u=-u;
		idx=i_base;
		val=v_base;
		reset=read_line_register(u);

		/* now choose some other value */
		u=((rand()&1)?1:-1)*(1+(rand()&MAX_CORRECTION));
		printf(" -> %d",u);
	}
	if (try)
		printf(" -> ok\n");
#else
	read_line_register(u);
#endif
	return u;
}

int main(int argc, char *argv[])
{
	int fd;
	struct stat sbuf;
	type32 i;
	size_t isize;
	size_t vsize;
	type32 ii;
	FILE *f;
	struct opt_desc * opts = NULL;
	int n_opts = 0;

	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	new_option(&opts,&n_opts,
		OPT_INTPRM(NULL,&bank_num,OPT_INDEP|OPT_REQUIRED));

	check_endianness(stdout);

	process_options(argc,argv,n_opts,opts);

	set_all_filenames(bank_num);

	read_tag_file();


	printf("Reading indexes\n");
	fd=open(w_indexes_filename,O_RDONLY);
	if (fd<0) {perror("open(indexes)"); exit(errno);}
	if (fstat(fd,&sbuf)<0) {perror("fstat(indexes)"); exit(errno);}
	isize=sbuf.st_size;
	printf("indexes : %lu bytes\n",(unsigned long)isize);
	indexes=(void*)mmap(NULL,isize,PROT_READ,MAP_SHARED,fd,0);
	if (indexes==NULL) {perror("mmap(indexes)"); exit(errno);}
	if (close(fd)<0) {perror("close(indexes)"); exit(errno);}

	printf("Reading values\n");
	fd=open(w_values_filename,O_RDONLY);
	if (fd<0) {perror("open(values)"); exit(errno);}
	if (fstat(fd,&sbuf)<0) {perror("fstat(values)"); exit(errno);}
	vsize=sbuf.st_size;
	printf("values : %lu bytes\n",(unsigned long)vsize);
	values=(void*)mmap(NULL,vsize,PROT_READ,MAP_SHARED,fd,0);
	if (values==NULL) {perror("mmap(values)"); exit(errno);}
	if (close(fd)<0) {perror("close(values)"); exit(errno);}

	/* A short note about the compressed matrix format (this note
	 * should go in c/prep/prep.c as well).
	 *
	 * a line is made up of several modes. There can be up to 3
	 * different modes. The modes govern the size of the index shift
	 * field. The value field is always 32bit-wide. At any given
	 * point, during the scan of the line, we hava access to the next
	 * I (index shift) value, as well as the V (coeff) value.
	 *
	 * If (I==0), in general this indicates a mode change. However,
	 * as long as no non-zero coefficient has been met so far in the
	 * course of building the matrix, we cannot guarantee that this
	 * doesn't simply mean that a coefficient is there at index 0.
	 * Hence, a mode change on a line which is still empty is
	 * indicated by (I==0) && (V==0). Changing mode, when there is
	 * only one mode, or when we are already in the last processing
	 * mode, means ending the line. This implies that an empty line
	 * is indicated by (I==0) && (V==0) as many times as we have
	 * modes.
	 *
	 * A line with k non-zero coeffs, featuring z mode changes at
	 * points where the line is still empty, for a total of m
	 * available modes, occupies k+m different indices, and k+z
	 * different values. Assuming z is zero most of the time, the
	 * difference vsize-isize is N*m, where N is the number of rows.
	 */

	printf("we have %d lines\n", nrows);
	
	/*
	 * For now, we use (1,1,1,...) as the check vector. Intermixing
	 * \pm 1's could be a good idea.
	 */
	printf("Building a random vector with %d coefficients\n", nrows);
	eigenvec = my_malloc(nrows*sizeof(stype32));
	mvec	 = my_malloc(nrows*sizeof(stype32));
	for(i=0;i<nrows;i++)
		eigenvec[i]=1;
	memset(mvec,0,nrows*sizeof(stype32));

	printf("Multiplying by the matrix\n");
	idx=val=0;
	for(i=0;i<nrows;i++) {
		eigenvec[i]=read_line_nocancel(i, eigenvec[i]);
	}

	printf("Counting zero coefficients\n");
	ii=0;
	for(i=0;i<nrows;i++) {
		if (mvec[i]) ii++;
		else {
			printf("m_{%d}==0\n",i);
		}
	}
	
	printf("Computed vector m has %d non-zero coeffs, dim=%d\n",ii,nrows);

	f=fopen(mvec_filename,"w");
	fwrite(mvec,sizeof(stype32),nrows,f);
	fclose(f);

	f=fopen(x0_filename,"w");
	fwrite(eigenvec,sizeof(stype32),nrows,f);
	fclose(f);

	munmap((void*)indexes,isize);
	munmap((void*)values,vsize);

	return 0;
}



