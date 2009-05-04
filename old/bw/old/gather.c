#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "slave_params.h"
#include "params.h"
#include "types.h"
#include "old-endian.h"
#include "macros.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "bw_lvblock.h"
#include "bw_lvblock_steps.h"
#include "variables.h"
#include "options.h"
#include "slave.h"
#include "filenames.h"
#include "matrix.h"
#include "tagfile.h"
#include "timer.h"
#include "threaded.h"
#include "addmul.h"
/* #include "version.h" */
#include "gmp-hacks.h"

int computed_nbxs=1;
int computed_nbys=1;

static int differential=1;	/* on by default */

coord_t *bw_x;		/* These are just elements of the canonical basis */

void showuse(void)
{
	die("Usage : bw-gather <bank#>\n",1);
}


int write_final_result(bw_vector_block t)
{
	FILE *f;

	f=fopen(w_filename,"w");
	if (f==NULL)
		die("fopen(%s) : %s\n",1,w_filename,strerror(errno));
	bw_lvblock_write(f,t);
	fclose(f);

	return 0;
}

int also_write_de_extended_thing(bw_vector_block t)
{
	mpz_t	c0;
	mpz_t	c;
	bw_vector_block	xv;
	FILE	* f;
	int	j;
	char	filename[FILENAME_LENGTH];

	xv = t;
	
	mpz_init(c0);
	mpz_init(c);

	MPZ_SET_MPN(c0,xv,bw_allocsize+2);

	mpz_neg(c0,c0);
	mpz_invert(c0,c0,modulus);

	for(j=0;j<ncols;j++) {
		MPZ_SET_MPN(c,xv,bw_allocsize+2);
		mpz_mul(c,c,c0);
		mpz_mod(c,c,modulus);
		MPN_SET_MPZ(xv,bw_allocsize,c);
		bw_lvblock_step_p00(xv);
	}

	mpz_clear(c);
	mpz_clear(c0);

	sprintf(filename,"%s.unext",w_filename);
	f=fopen(filename,"w");
	if (f==NULL)
		die("fopen(%s) : %s\n",1,filename,strerror(errno));
	xv=t;
	bw_lvblock_step_p00(xv);
	ncols--;
	bw_lvblock_write(f,xv);
	ncols++;
	fclose(f);

	return 0;
}

int add_file(bw_vector_block dst, const char * name)
{
	bw_vector_block tmp, xw, xv;
	mp_limb_t c;
	FILE *f;
	int	j;

	bw_lvblock_alloc(tmp);
	bw_lvblock_set_zero_separated(tmp);

	f=fopen(name,"r");
	if (f==NULL)
		die("fopen(%s) : %s\n",1,name,strerror(errno));

	printf("Reading %s\n",name);

	bw_lvblock_read(tmp,f);

	xw=tmp;
	xv=dst;

	for(j=0;j<ncols;j++) {
		c=mpn_add_n(xv,xv,xw,bw_allocsize+2);
		if (c) {
			fprintf(stderr,"Uh-uh, overflow.\n");
			eternal_sleep();
		}
		bw_lvblock_step_p00(xw);
		bw_lvblock_step_p00(xv);
	}
	fclose(f);
	bw_lvblock_free(tmp);

	return 0;
}

void add_differential_files(bw_vector_block dst)
{
	DIR * dirp;
	struct dirent * currdir;
	char pattern[FILENAME_LENGTH];
	char filename[FILENAME_LENGTH];
	char sfilename[FILENAME_LENGTH];
	char *s;
	int v,i,r;

	dirp=opendir(wdir_filename);
	s=strrchr(ds_meta_filename,'/');
	if (s!=NULL) {
		s++;
	} else {
		s=ds_meta_filename;
	}
	mkfname(pattern,"%s",s);
	
	for(;(currdir=readdir(dirp))!=NULL;) {
		if (sscanf(currdir->d_name,pattern,&v,&i,&r)==3) {
			mkfname(filename,pattern,v,i,r);
			if (strcmp(filename,currdir->d_name)!=0)
				continue;
			mkfname(filename,ds_meta_filename,v,i,r);
			mkfname(sfilename,"%s",filename);
			add_file(dst,sfilename);
		}
	}
	closedir(dirp);
}


int main(int argc, char *argv[])
{
	char name[FILENAME_LENGTH];
	int i=-1;
	bw_vector_block bw_v, bw_w, sol;
	int l;
	struct opt_desc * opts = NULL;
	int n_opts=0;

	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	/* puts(version_string); */
	check_endianness(stdout);

	new_option(&opts, &n_opts,
		OPT_FLAG("--differential",&differential));
	new_option(&opts, &n_opts,
		OPT_INTPRM(NULL,&bank_num,OPT_INDEP|OPT_REQUIRED));

	process_options(argc,argv,n_opts,opts);

	set_all_filenames(bank_num);
	read_tag_file();

#ifdef HARDCODE_PARAMS
	consistency_check("nbys",computed_nbys,nbys);
	consistency_check("n_param",computed_n_param,n_param);
	consistency_check("bw_allocsize",computed_bw_allocsize,bw_allocsize);
	consistency_check("bw_longsize",computed_bw_longsize,bw_longsize);
#endif

	load_matrix();

	configure_threads(1,ncols);
	compute_thread_offsets();

	bw_x=malloc(nbxs*sizeof(coord_t));
	load_x_vectors(bw_x,0,nbxs);

	bw_lvblock_alloc(bw_v);
	bw_lvblock_set_zero_separated(bw_v);

	if (differential) {
		add_differential_files(bw_v);
	} else {
		for(l=0;l<n_param;l++) {
			sprintf(name,h_meta_filename,l);
			add_file(bw_v,name);
		}
	}

	bw_lvblock_reduce_separated(bw_v);

	bw_lvblock_alloc(bw_w);
	bw_lvblock_set_zero_separated(bw_w);

	sol = NULL;

	for(i=0;i<100;) {
		printf("Trying B^%d * w\n",i);
		multiply(bw_w,bw_v);
		bw_lvblock_reduce_separated(bw_w);
		if (bw_lvblock_is_zero_separated(bw_w)) {
			printf("B^%d * w is a solution.\n",i);
			sol=bw_v;
			break;
		} else {
			int j;
			bw_vector_block l=bw_w;
			for(j=0;j<nbxs;j++) {
				bw_lvblock_step_n00(l,bw_x[j]);
				if (!bw_scalar_is_zero(l)) {
					/* If this one is consistently
					 * non-null, then there's a bug
					 * */
					printf("(B^%d*w)[%d] is not 0! means BUG if consistently repeated :-((\n",
					i,bw_x[j]);
				}
			}
		}

		bw_lvblock_set_zero_separated(bw_v);

		i++;

		printf("Trying B^%d * w\n",i);
		multiply(bw_v,bw_w);
		bw_lvblock_reduce_separated(bw_v);
		if (bw_lvblock_is_zero_separated(bw_v)) {
			printf("B^%d * w is a solution.\n",i);
			sol=bw_w;
			break;
		} else {
			int j;
			bw_vector_block l=bw_v;
			for(j=0;j<nbxs;j++) {
				bw_lvblock_step_n00(l,bw_x[j]);
				if (!bw_scalar_is_zero(l)) {
					/* If this one is consistently
					 * non-null, then there's a bug
					 * */
					printf("(B^%d*w)[%d] is not 0! means BUG if consistently repeated :-((\n",
					i,bw_x[j]);
				}
			}
		}
		bw_lvblock_set_zero_separated(bw_w);

		i++;
	}
	
	if (sol == NULL) {
		fprintf(stderr,
			"Mmmh, no solution before B^100 * w, I give up.\n");
		exit(0);
	}

	if (bw_lvblock_is_zero_separated(sol)) {
		fprintf(stderr,"Trivial solution encountered !\n");
		exit(0);
	}
	write_final_result(sol);
	
	also_write_de_extended_thing(sol);
		
	free(bw_x);

	return 0;
}

/* vim:set sw=8: */
