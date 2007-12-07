#ifdef	KNOWN_LITTLE_ENDIAN
#define USE_MMAP
#endif

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "slave_params.h"
#include "params.h"
#include <assert.h>
#include "types.h"
#include "macros.h"
#include "old-endian.h"
#include "auxfuncs.h"
#include "filenames.h"
#include "bw_scalar.h"
#include "bw_lvblock.h"
#include "slave.h"
#include "bw_lvblock_steps.h"
#include "matrix.h"
#include "threaded.h"
#include "addmul.h"

size_t indexes_size;
size_t values_size;

#ifndef KNOWN_LITTLE_ENDIAN
#define reorder_indexes() swaporder_raw_32bit(matrix.indexes,indexes_size)
#endif

#ifndef KNOWN_LITTLE_ENDIAN
void swaporder_raw_32bit(type32 *x, size_t n)
{
	size_t p;
	for(p = 0 ; p * sizeof(type32) < n ; p++)
		mswap32(x[p]);
}
#endif

void multiply(bw_vector_block w, bw_vector_block v)
{
	coord_t i;
	type32 *now;
	stype32 *val;
	bw_vector_block  tw, tv;
	type32 step;

	now=matrix.indexes;
	val=matrix.values;

	tw=w;
	
	now+=tsv()->idx_off;
	val+=tsv()->val_off;
	bw_lvblock_step_n00(tw,tsv()->i0);

	for(i=tsv()->i0;i<tsv()->i1;i++,bw_lvblock_step_p00(tw),now++) {
		tv=v;
		step=*now;
		if (EXPECT_FALSE(step==0)) {
			if (*val==0) {
				val++;
				continue;
			}
		}
		do {
			now++;
			bw_lvblock_step_n00(tv,step);
			addmultiply(tw,tv,(stype32)*val);
			val++;
			step=*now;
		} while (step!=0);
	}
}

int compute_thread_offsets(void)
{
	type32 *now;    
	stype32 *val;    
	coord_t i;
	size_t ioff;
	size_t voff;
	int	nc=0;

	now=matrix.indexes;
	val=matrix.values;

	ioff=0;
	voff=0;
	for(i=0;i<tsv()->i0;i++,now++,ioff++) {
		if (*now==0 && *val==0) {
			val++;
			voff++;
			/* We fall here when there is a blank line. This
			 * is indeed the case when the system has been
			 * extended to account for it being inhomogoneous
			 */
			continue;
		}
		do { now++; ioff++; val++; voff++; } while (*now!=0);
	}
	fprintf(stderr,"Thread %d : idx_off=%08lx, val_off=%08lx\n",
			tseqid(),(unsigned long)ioff,(unsigned long)voff);
	tsv()->idx_off=ioff;
	tsv()->val_off=voff;
	for(;i<tsv()->i1;i++,now++) {
		if (*now==0 && *val==0) {
			val++;
			/* We fall here when there is a blank line. This
			 * is indeed the case when the system has been
			 * extended to account for it being inhomogoneous
			 */
			continue;
		}
		do { now++; val++; nc++; } while (*now!=0);
	}
	return nc;
}

void load_matrix(void)
{
#ifdef USE_MMAP
	int fd;
#else
	FILE *f;
#endif
	struct stat sbuf;

	/* values (easier) */

	if (stat(w_values_filename,&sbuf)<0)
		die("stat(%s) : %s\n",1,w_values_filename,strerror(errno));

	values_size=sbuf.st_size;
	
#ifdef USE_MMAP
	DO_BIG_ENDIAN(die("Please convert the data to big-endian first\n",1););
	fd=open(w_values_filename,O_RDONLY);
	if (fd<0) {perror("open(values)"); exit(errno);}
	matrix.values=(void*)mmap(NULL,values_size,PROT_READ,MAP_SHARED,fd,0);
	if (matrix.values==NULL) {perror("mmap(values)"); exit(errno);}
	if (close(fd)<0) {perror("close(values)"); exit(errno);}
#else
	f=fopen(w_values_filename,"r");
	if (f==NULL)
		die("fopen(%s) : %s\n",1,w_values_filename,strerror(errno));
	
	matrix.values=mymalloc(values_size);

	fprintf(stderr,"Loading the matrix into memory (values)...");
	fflush(stderr);
	fread(matrix.values,values_size,1,f);
	fprintf(stderr,"done\n");
	fclose(f);

	DO_BIG_ENDIAN(
		fprintf(stderr,"Adapting to host byte order (big endian)...");
		swaporder_raw_32bit((type32*)matrix.values,values_size);
		fprintf(stderr,"done\n")
	);
#endif

	/* indexes */
	
	if (stat(w_indexes_filename,&sbuf)<0)
		die("stat(%s) : %s\n",1,w_indexes_filename,strerror(errno));

	indexes_size=sbuf.st_size;

#ifdef USE_MMAP
	DO_BIG_ENDIAN(die("Please convert the data to big-endian first\n",1););
	fd=open(w_indexes_filename,O_RDONLY);
	if (fd<0) {perror("open(indexes)"); exit(errno);}
	matrix.indexes=(void*)mmap(NULL,indexes_size,PROT_READ,MAP_SHARED,fd,0);
	if (matrix.indexes==NULL) {perror("mmap(indexes)"); exit(errno);}
	if (close(fd)<0) {perror("close(indexes)"); exit(errno);}
#else
	f=fopen(w_indexes_filename,"r");
	if (f==NULL)
		die("fopen(%s) : %s\n",1,w_indexes_filename,strerror(errno));
	
	matrix.indexes=mymalloc(indexes_size);

	fprintf(stderr,"Loading the matrix into memory (indexes)...");
	fflush(stderr);
	fread(matrix.indexes,indexes_size,1,f);
	fprintf(stderr,"done\n");
	fclose(f);

	DO_BIG_ENDIAN(
		fprintf(stderr,"Adapting to host byte order (big endian)...");
		reorder_indexes();
		fprintf(stderr,"done\n");
	);
#endif
}

/* vim:set sw=8: */
