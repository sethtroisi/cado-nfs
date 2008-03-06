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
#include "types.h"
#include "macros.h"
#include "old-endian.h"
#include "auxfuncs.h"
#include "variables.h"
#include "tagfile.h"
#include "filenames.h"

const char * indexes_fmt="run/output/%d/indexes";
const char * permutation_fmt="run/output/%d/balancing";
const char * values_fmt="run/output/%d/values";

#define FILENAME_LENGTH 80

type32 * indexes;
type32 * values;
size_t isize;
size_t vsize;

#define IMBRICATED_OPTS	6

const int give_details[]={2,3,4,6,8,12,24,1024};
const int opt_chain[IMBRICATED_OPTS]={1,3072,24,8,4,2};

struct single_line {
	size_t	idx_off;
	size_t	val_off;
	size_t	idx_size;
	size_t	val_size;
	type32	w;
};

struct line_packet {
	type32 n;
	type32 w;
	type32 * lines;
	type32 alloc;
};

struct packet_family {
	type32 np;
	struct line_packet * p;
};


struct single_line * line_pool;
struct packet_family family[IMBRICATED_OPTS];

type32 nlines;
type32 ncoeffs;
type32 pool_alloc;

int * flattened_lines;
int nflattened;


int read_line(struct single_line * l, stype32 * p_idx, stype32 * p_val)
{
	stype32 idx,val;

	idx=*p_idx;
	val=*p_val;
	l->w=0;
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
			l->w++;
		} while (indexes[idx]!=0);
	}
	idx++;
	l->idx_size++;
	*p_idx=idx;
	*p_val=val;
	return l->w;
}

void ext_ptr(void ** pp, type32 i, type32 * plen, size_t s)
{
	if (i >= *plen) {
		*plen=(*plen!=0)?(2* *plen):(1+*plen);
		*pp=realloc(*pp,*plen*s);
	}
}

void give_info_before(void)
{
	int i,j,k;
	type32 w;
	int g;

	printf("Repartition scheme, depending on the number of CPUs\n");
	for(i=0;i<sizeof(give_details)/sizeof(give_details[0]);i++) {
		g=give_details[i];
		printf("N=%-2d | ",g);
		for(j=0;j<g;j++) {
			w=0;
			for(k=j*nlines/g;k<(j+1)*nlines/g;k++)
				w+=line_pool[k].w;
			printf("%2.2f%% ",100.0*(double)w/ncoeffs);
		}
		printf("\n");
	}
}

int compute_ordering(int * order, int flow, int fnum, int pnum)
{
	int n=0,i,k;
	if (fnum==flow) {
		order[n++]=pnum;
	} else for(i=0;i<family[fnum].p[pnum].n;i++) {
		k=family[fnum].p[pnum].lines[i];
		n+=compute_ordering(order+n,flow,fnum-1,k);
	}
	return n;
}

void give_info_after(void)
{
	int i,j,k,g,n,w,nl,flow;
	int * order;

	printf("Repartition scheme, depending on the number of CPUs\n");

	for(i=0;i<sizeof(give_details)/sizeof(give_details[0]);i++) {
		g=give_details[i];
		for(flow=IMBRICATED_OPTS-1;flow>0;flow--) {
			if (family[flow].np % g == 0)
				break;
		}
		/* We don't want to work on the family [0]. Therefore, we
		 * can give_details only for g=divisor of the
		 * opt_chain[1] */
		ASSERT(flow>0);
		
		order=malloc(family[flow].np * sizeof(int));
		for(k=n=0;k<family[IMBRICATED_OPTS-1].np;k++) {
			n+=compute_ordering(order+n,flow,IMBRICATED_OPTS-1,k);
		}

		printf("N=%-2d | ",g);
		for(j=0;j<g;j++) {
			w=0;
			for(k=j*family[flow].np/g;k<(j+1)*family[flow].np/g;k++)
				w+=family[flow].p[order[k]].w;
			printf("	%2.2f%%",100.0*(double)w/ncoeffs);
		}
		printf("\n");
		printf("N=%-2d |+",g);
		for(j=0;j<g;j++) {
			w=0;
			for(k=j*family[flow].np/g;k<(j+1)*family[flow].np/g;k++)
				w+=family[flow].p[order[k]].w;
			printf("	%d",w);
		}
		printf("\n");
		printf("N=%-2d ||",g);
		for(j=0;j<g;j++) {
			nl=0;
			for(k=j*family[flow].np/g;k<(j+1)*family[flow].np/g;k++)
				nl+=family[flow].p[order[k]].n;
			printf("	%d",nl);
		}
		printf("\n");
		free(order);
	}
}

void map_files(int num)
{
	char fn_indexes[FILENAME_LENGTH];
	char fn_values[FILENAME_LENGTH];
	int fd;
	struct stat sbuf;

	sprintf(fn_indexes,indexes_fmt,num);
	sprintf(fn_values,values_fmt,num);

	fd=open(fn_indexes,O_RDONLY);
	if (fd<0) {perror("open(indexes)"); exit(errno);}
	if (fstat(fd,&sbuf)<0) {perror("fstat(indexes)"); exit(errno);}
	isize=sbuf.st_size;
	printf("indexes : %lu bytes\n",(unsigned long)isize);
	indexes=(void*)mmap(NULL,isize,PROT_READ,MAP_SHARED,fd,0);
	if (indexes==NULL) {perror("mmap(indexes)"); exit(errno);}
	if (close(fd)<0) {perror("close(indexes)"); exit(errno);}

	fd=open(fn_values,O_RDONLY);
	if (fd<0) {perror("open(values)"); exit(errno);}
	if (fstat(fd,&sbuf)<0) {perror("fstat(values)"); exit(errno);}
	vsize=sbuf.st_size;
	printf("values : %lu bytes\n",(unsigned long)vsize);
	values=(void*)mmap(NULL,vsize,PROT_READ,MAP_SHARED,fd,0);
	if (values==NULL) {perror("mmap(values)"); exit(errno);}
	if (close(fd)<0) {perror("close(values)"); exit(errno);}
}

void unmap_files(void)
{
	munmap((void*)indexes,isize);
	munmap((void*)values,vsize);
}

/********************************************************************/

void read_single_lines(void)
{
	stype32 idx,val;
	int i,r;

	ncoeffs=nlines=0;
	pool_alloc=0;
	line_pool=NULL;

	idx=val=0;
	for(i=0;;i++) {
		void * foo = line_pool;
		ext_ptr(&foo,i,&pool_alloc,sizeof(struct single_line));
		line_pool = foo;
		r=read_line(&line_pool[i],&idx,&val);
		ncoeffs+=r;
		if (	idx * sizeof(stype32) == isize &&
			val * sizeof(stype32) == vsize)
		{
			i++;
			break;
		}
		ASSERT(idx * sizeof(stype32) < isize);
		ASSERT(val * sizeof(stype32) < vsize);

		if (r==0) {
			printf("Blank line before end of file\n");
			continue;
		}
	}

	nlines=i;

	printf("Read %lu lines, totalizing %lu coeffs\n",
			(unsigned long)nlines,
			(unsigned long)ncoeffs);
	if (idx * sizeof(type32) != isize || val * sizeof(type32) != vsize) {
		printf("File size not reached.\n");
		printf("Stopped at offset %lu(indexes):%lu(values)\n",
				(unsigned long)idx * sizeof(type32),
				(unsigned long)val * sizeof(type32));
	}

	/* consolidate pointer */

	line_pool = realloc(line_pool,nlines*sizeof(struct single_line));
}

int cmp_packets(const void *a, const void *b)
{
	return  (int)(((struct line_packet*)b)->w)-
		(int)(((struct line_packet*)a)->w);
}

void balance_packets(struct packet_family * dst, struct packet_family * src)
{
	type32 i,j,k,*order;

	/* Sort the packets */
	qsort(src->p,src->np,sizeof(struct line_packet),&cmp_packets);
	
	order=malloc(dst->np*sizeof(type32));
	for(i=0;i<dst->np;i++) {
		order[i]=i;
	}

	for(i=0;i<src->np;i++) {
		j=order[0];
		/* Put the heaviest available packet in the lightest
		 * bucket */
		ext_ptr((void**)&(dst->p[j].lines),dst->p[j].n,
				&(dst->p[j].alloc),sizeof(type32));
		dst->p[j].w+=src->p[i].w;
		dst->p[j].lines[dst->p[j].n++]=i;
		for(k=1;k<dst->np && dst->p[order[k]].w<dst->p[j].w;k++);
		memmove(order,order+1,(k-1)*sizeof(int));
		order[k-1]=j;
	}

	free(order);
}

void optimization_chain(void)
{
	type32 i,j;

	/* Build the original packet structure */
	family[0].np = nlines;
	family[0].p  = malloc(nlines*sizeof(struct line_packet));

	for(j=0;j<nlines;j++) {
		family[0].p[j].n=1;
		family[0].p[j].w=line_pool[j].w;
		family[0].p[j].lines=malloc(sizeof(type32));
		family[0].p[j].lines[0]=j;
		family[0].p[j].alloc=1;
	}

	/***** Now start the chain of optimizations ******/

	for(i=1;i<IMBRICATED_OPTS;i++) {
		family[i].np=opt_chain[i];
		family[i].p=malloc(opt_chain[i]*sizeof(struct line_packet));
		for(j=0;j<family[i].np;j++) {
			family[i].p[j].n=0;
			family[i].p[j].w=0;
			family[i].p[j].lines=NULL;
			family[i].p[j].alloc=0;
		}

		balance_packets(&family[i],&family[i-1]);
	}
}


void flatten_lines(int fnum, int pnum)
{
	int i;
	for(i=0;i<family[fnum].p[pnum].n;i++) {
		int k=family[fnum].p[pnum].lines[i];
		if (fnum==0) {
			flattened_lines[nflattened++]=k;
		} else {
			flatten_lines(fnum-1,k);
		}
	}
}

void write_permutation(FILE *f, int fnum, int pnum)
{
	type32 i,k;
	for(i=0;i<family[fnum].p[pnum].n;i++) {
		k=family[fnum].p[pnum].lines[i];
		if (fnum==0) fprintf(f, " %u", k);
		else write_permutation(f,fnum-1,k);
	}
}

void write_packet_indexes(FILE *f, int fnum, int pnum)
{
	type32 i,k;

	for(i=0;i<family[fnum].p[pnum].n;i++) {
		k=family[fnum].p[pnum].lines[i];
		if (fnum==0) {
			fwrite(indexes+line_pool[k].idx_off,
					sizeof(type32),
					line_pool[k].idx_size,
					f);
		} else {
			write_packet_indexes(f,fnum-1,k);
		}
	}
}

void write_packet_values(FILE *f, int fnum, int pnum)
{
	type32 i,k;

	for(i=0;i<family[fnum].p[pnum].n;i++) {
		k=family[fnum].p[pnum].lines[i];
		if (fnum==0) {
			fwrite(values+line_pool[k].val_off,
					sizeof(type32),
					line_pool[k].val_size,
					f);
		} else {
			write_packet_values(f,fnum-1,k);
		}
	}
}

void final_flush(int num)
{
	char fn_indexes_old[FILENAME_LENGTH];
	char fn_values_old[FILENAME_LENGTH];
	char fn_indexes[FILENAME_LENGTH];
	char fn_values[FILENAME_LENGTH];
	char fn_permutation[FILENAME_LENGTH];

	FILE * f_i;
	FILE * f_v;
	FILE * f_p;
	int j;

	/* Change file names */
	sprintf(fn_indexes,indexes_fmt,num);
	sprintf(fn_permutation,permutation_fmt,num);
	sprintf(fn_values,values_fmt,num);
	strcpy(fn_indexes_old,fn_indexes);
	strcpy(fn_values_old,fn_values);
	strcpy(fn_indexes_old + strlen(fn_indexes), ".old");
	strcpy(fn_values_old  + strlen(fn_values),  ".old");

	if (rename(fn_indexes,fn_indexes_old)<0) {
		fprintf(stderr, "Error while renaming %s : %s\n",
				fn_indexes,strerror(errno));
		exit(errno);
	}
	if (rename(fn_values,fn_values_old)<0) {
		fprintf(stderr, "Error while renaming %s : %s\n",
				fn_values,strerror(errno));
		exit(errno);
	}

	printf("Writing permutation to %s\n",fn_permutation);
	f_p=fopen(fn_permutation,"w");
	fprintf(f_p,"[");
	for(j=0;j<family[IMBRICATED_OPTS-1].np;j++) {
		write_permutation(f_p,IMBRICATED_OPTS-1,j);
	}
	fprintf(f_p," ]\n");
	fclose(f_p);

	printf("Writing new indexes to %s\n",fn_indexes);
	f_i=fopen(fn_indexes,"w");
	for(j=0;j<family[IMBRICATED_OPTS-1].np;j++) {
		write_packet_indexes(f_i,IMBRICATED_OPTS-1,j);
	}
	fclose(f_i);

	printf("Writing new values to %s\n",fn_values);
	f_v=fopen(fn_values,"w");
	for(j=0;j<family[IMBRICATED_OPTS-1].np;j++) {
		write_packet_values(f_v,IMBRICATED_OPTS-1,j);
	}
	fclose(f_v);
}

void fixup_x_files()
{
	coord_t	* tab;
	int i,d;

	tab = malloc(m_param * sizeof(coord_t));
	load_x_vectors(tab, 0, m_param);

	d=0;
	for(i=0;i<m_param;i++) {
		for(;line_pool[flattened_lines[tab[i]+d]].w==0;d++) {
			fprintf(stderr,"WARNING: zero line (%d) on X-coordinate X%d. Avoiding.\n", tab[i]+d, i);
		}
		tab[i]+=d;
	}

	if (d!=0) {
		FILE * f;
		char filename[FILENAME_LENGTH];
		fprintf(stderr,"Fixing up X coordinates\n");
		for(i=0;i<m_param;i++) {
			snprintf(filename,FILENAME_LENGTH,x_meta_filename,i);
			f=fopen(filename,"w");
			if (f==NULL) die("%s : %s",1,filename,strerror(errno));
			fprintf(f,"e%03d\n",tab[i]);
			fclose(f);
		}
	}
}


int main(int argc, char *argv[])
{
	int num;
	int j;

	if (argc!=2) {
		fprintf(stderr,"give a bank number\n");
		exit(1);
	}
	num=atoi(argv[1]);

	set_all_filenames(num);

	read_tag_file();
	map_files(num);
	read_single_lines();
	give_info_before();
	optimization_chain();
	give_info_after();
	flattened_lines=malloc(nlines*sizeof(int));
	nflattened=0;
	for(j=0;j<family[IMBRICATED_OPTS-1].np;j++) {
		flatten_lines(IMBRICATED_OPTS-1,j);
	}
	ASSERT(nflattened==nlines);
	fixup_x_files();
	final_flush(num);
	unmap_files();

	return 0;
}
