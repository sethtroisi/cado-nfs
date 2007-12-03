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

const char * indexes_fmt="run/output/%d/indexes";
const char * values_fmt="run/output/%d/values";

#define FILENAME_LENGTH 80

type32 * indexes;
type32 * values;
size_t idx,val;

struct vec_coeff {
	type32 j;
	stype32 v;
};

struct vec_coeff *vec_coeff_tab=NULL;
type32 alloc_coeffs=0;
type32 present_coeffs=0;
void register_coeff(type32 j, stype32 v)
{
	if (present_coeffs >= alloc_coeffs) {
		alloc_coeffs=(alloc_coeffs!=0)?
				(2*alloc_coeffs):
				(alloc_coeffs+1);
		vec_coeff_tab=realloc(vec_coeff_tab,
				alloc_coeffs*sizeof(struct vec_coeff));
	}
	vec_coeff_tab[present_coeffs].j=j;
	vec_coeff_tab[present_coeffs].v=v;
	present_coeffs++;
}

int cmp_vec_coeffs(const void *a, const void *b)
{
	return  (int)(((struct vec_coeff*)a)->j)-
		(int)(((struct vec_coeff*)b)->j);
}

void read_line_register(void)
{
	type32 j=0;

	if (indexes[idx]==0 && values[val]==0) {
		val++;
	} else {
		do {
			j+=indexes[idx];
			register_coeff(j, values[val]);
			idx++;
			val++;
		} while (indexes[idx]!=0);
	}
	idx++;
}

void print_line(int nl)
{
	type32 j=(type32)-1,jj;
	if (indexes[idx]==0 && values[val]==0) {
		val++;
	} else {
		do {
			for(jj=j+1;jj<j+indexes[idx]+(j==(type32)-1);jj++)
				printf(" 0, ");
			if (j==(type32)-1)
				j++;
			j+=indexes[idx];
			printf(" %d, ", values[val]);
			idx++;
			val++;
		} while (indexes[idx]!=0);
	}
	for(jj=j+1;jj<nl;jj++)
		printf(" 0, ");
	idx++;
}

void sprint_line(int nl)
{
	type32 j=0;
	if (indexes[idx]==0 && values[val]==0) {
		val++;
	} else {
		do {
			j+=indexes[idx];
			printf("%d:%d ", j,values[val]);
			idx++;
			val++;
		} while (indexes[idx]!=0);
	}
	idx++;
}


int main(int argc, char *argv[])
{
	char fn_indexes[FILENAME_LENGTH];
	char fn_values[FILENAME_LENGTH];
	int fd;
	struct stat sbuf;
	type32 i,k,last_j;
	type32 nc, nl;
	size_t isize;
	size_t vsize;
	type32 ii;
	int print=0;
	int sprint=0;

	if (argc>=3 && strcmp(argv[1],"print")==0) {
		argv++;
		argc--;
		print=1;
	}

	if (argc>=3 && strcmp(argv[1],"sprint")==0) {
		argv++;
		argc--;
		sprint=1;
	}

	if (argc!=2)
		die("Please give a bank number\n",1);

	sprintf(fn_indexes,indexes_fmt,atoi(argv[1]));
	sprintf(fn_values,values_fmt,atoi(argv[1]));

	printf("Reading indexes\n");
	fd=open(fn_indexes,O_RDONLY);
	if (fd<0) {perror("open(indexes)"); exit(errno);}
	if (fstat(fd,&sbuf)<0) {perror("fstat(indexes)"); exit(errno);}
	isize=sbuf.st_size;
	printf("indexes : %lu bytes\n",(unsigned long)isize);
	indexes=(void*)mmap(NULL,isize,PROT_READ,MAP_SHARED,fd,0);
	if (indexes==NULL) {perror("mmap(indexes)"); exit(errno);}
	if (close(fd)<0) {perror("close(indexes)"); exit(errno);}

	printf("Reading values\n");
	fd=open(fn_values,O_RDONLY);
	if (fd<0) {perror("open(values)"); exit(errno);}
	if (fstat(fd,&sbuf)<0) {perror("fstat(values)"); exit(errno);}
	vsize=sbuf.st_size;
	printf("values : %lu bytes\n",(unsigned long)vsize);
	values=(void*)mmap(NULL,vsize,PROT_READ,MAP_SHARED,fd,0);
	if (values==NULL) {perror("mmap(values)"); exit(errno);}
	if (close(fd)<0) {perror("close(values)"); exit(errno);}

	nc = vsize / sizeof(type32) - 1;
	nl = isize / sizeof(type32) - nc;

	printf("Assuming the system was extended from an inhomogeneous one,\n"
		"we have %lu coefficients in %lu lines\n",
		(unsigned long) nc,
		(unsigned long) nl);
	
	if (print) {
		printf("Printing the matrix\n[ ");
		idx=val=0;
		for(i=0;i<nl;i++) {
			print_line(nl);
			printf("\n  ");
		}
		printf("]\n");
	}

	if (sprint) {
		printf("Printing the matrix (sparse form)\n");
		idx=val=0;
		for(i=0;i<nl;i++) {
			sprint_line(nl);
			printf("\n");
		}
	}

	printf("Analysing the matrix\n");
	idx=val=0;
	for(i=0;i<nl;i++) {
		read_line_register();
	}

	printf("Matrix has %d terms. sorting\n",present_coeffs);

	qsort(vec_coeff_tab,present_coeffs,
			sizeof(struct vec_coeff),&cmp_vec_coeffs);

	last_j=(type32)-1;
	ii=(type32)-1;

	assert(present_coeffs);

	for(k=0;k<present_coeffs;k++) {
		if (vec_coeff_tab[k].j!=last_j) {
			if (last_j!=(type32)-1) {
				printf("Column %d :", last_j);
				printf(" <%d coeffs>",k-ii);
				if (k-ii < 32) {
					for(i=ii;i<k;i++) {
						printf(" %d",
							vec_coeff_tab[i].v);
					}
				}
				printf("\n");
			}
			ii=k;
			last_j=vec_coeff_tab[k].j;
		}
	}
	printf("Column %d :", last_j);
	printf(" <%d coeffs>",k-ii);
	if (k-ii < 32) {
		for(i=ii;i<k;i++) {
			printf(" %d", vec_coeff_tab[i].v);
		}
	}
	printf("\n");
	
	munmap((void*)indexes,isize);
	munmap((void*)values,vsize);

	return 0;
}



