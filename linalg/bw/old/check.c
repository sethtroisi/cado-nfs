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

const char * c_fmt="run/output/%d/C<%02d>.%03d";
const char * v_fmt="run/output/%d/V<%02d>.%03d";
const char * x0_fmt="run/output/%d/X0-vector";
const char * lambda_fmt="run/output/%d/eigenvalue";

#define FILENAME_LENGTH 80

struct vec_coeff {
	type32 j;
	stype32 v;
};

struct vec_coeff *vec_coeff_tab=NULL;
type32 alloc_coeffs=0;
type32 present_coeffs=0;

int main(int argc, char *argv[])
{
	char fn_c[FILENAME_LENGTH];
	char fn_v[FILENAME_LENGTH];
	char fn_x0[FILENAME_LENGTH];
	char fn_lambda[FILENAME_LENGTH];
	int fd;
	struct stat sbuf;
	type32 i;
	type32 nc, nl;
	type32 ii,k,last_j;
	type32 eigenvec[SPARSE_VECTOR_WEIGHT];
	stype32 eigenvalue;
	FILE *f;
	int banknum,colnum;

	if (argc!=3)
		die("Usage : ./check <bank> <column>\n",1);
	banknum=atoi(argv[1]);
	colnum=atoi(argv[2]);

	printf("Checking column %d in bank %d\n",banknum, colnum);

	sprintf(fn_x0,x0_fmt,banknum);
	f=fopen(fn_x0,"w");
	for(i=0;i<SPARSE_VECTOR_WEIGHT;i++) {
		fscanf(f,"%lu",(unsigned long*) &(eigenvec[i]));
	}
	fclose(f);

	sprintf(fn_lambda,lambda_fmt,banknum);
	f=fopen(fn_lambda,"w");
	fscanf(f,"%ld",(long*) &eigenvalue);
	fclose(f);

	printf("Read eigenvector : ");
	for(i=0;i<SPARSE_VECTOR_WEIGHT;i++) {
		printf("e%d + ",eigenvec[i]);
	}
	printf("\n");
	printf("Read eigenvalue is %ld\n",(long)eigenvalue);
	
	return 0;
}



