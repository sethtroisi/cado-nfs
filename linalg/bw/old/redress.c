#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "slave_params.h"
#include "params.h"
#include <assert.h>
#include "types.h"
#include "endian.h"
#include "macros.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "bw_lvblock.h"
#include "bw_lvblock_steps.h"
#include "variables.h"
#include "options.h"
#include "filenames.h"
#include "slave.h"
#include "matrix.h"
#include "tagfile.h"
#include "timer.h"
#include "threaded.h"
#include "addmul.h"
/* #include "version.h" */

/* This is a quick hack to redress the dS files that have been badly
 * modified...
 */
int computed_nbxs=4;
int computed_nbys=4;

coord_t *bw_x;

struct fn_cell {
	char		* fn;
	struct fn_cell  * next;
};

struct fn_cell * filename_ll = NULL;

void	read_vec(bw_vector_block dst, const char * fn)
{
	FILE *f;
	f=fopen(fn,"r");
	if (f==NULL) {
		perror(fn);
		exit(errno);
	}
	bw_lvblock_read(dst,f);
	fclose(f);
}

void	write_vec(const char * fn, bw_vector_block dst)
{
	FILE *f;
	f=fopen(fn,"w");
	if (f==NULL) {
		perror(fn);
		exit(errno);
	}
	bw_lvblock_write(f,dst);
	fclose(f);
}


int 
append_ll_filename(char **argv,
		struct extended_option_desc * conf UNUSED_VARIABLE,
		void * p_val)
{
	struct fn_cell ** pl=p_val;

	/* Go to the end... */
	for(;*pl!=NULL;pl=&((*pl)->next));

	*pl=malloc(sizeof(struct fn_cell));
	(*pl)->next=NULL;
	(*pl)->fn=strdup(argv[0]);
	return 1;
}

void	nds_name(char * dest, const char * nm)
{
	const char * s;
	int	bl;

	s=strrchr(nm,'/');
	if (s==NULL)
		s=nm;
	else
		s++;

	bl=(s-nm);
	memcpy(dest,nm,(s-nm));

	dest[bl++]='n';

	strcpy(dest+bl,s);
	return;
}

void showuse(void)
{
	die("Usage : redress <bank> /32/ <ncpus> [<file> ...]\n",1);
}

void compute_nds(bw_vector_block nds,
		bw_vector_block d0,
		bw_vector_block d1,
		int kval)
{
	int	i;
	bw_vector_block	nds_orig=nds;
#ifndef NDEBUG
	int	l;
#endif

	for(i=0;i<kval;i++) {
		memcpy(nds,d1,bw_allocsize*sizeof(mp_limb_t));
#ifndef NDEBUG
		for(l=bw_allocsize;l<bw_longsize;l++)
			assert(d1[l]==0);
#endif
		bw_lvblock_step_p00(d1);
		bw_lvblock_step_p00(d0);
		bw_lvblock_step_p00(nds);
	}

	for(;i<ncols;i++) {
		memcpy(nds,d1,bw_allocsize*sizeof(mp_limb_t));
		memcpy(nds+bw_allocsize+2,d0,bw_allocsize*sizeof(mp_limb_t));
#ifndef NDEBUG
		for(l=bw_allocsize;l<bw_longsize;l++) {
			assert(d0[l]==0);
			assert(d1[l]==0);
			/* Of course ndl[l] is not zero yet. It will
			 * become zero after the reduction. */
		}
#endif
		bw_lvblock_step_p00(d1);
		bw_lvblock_step_p00(d0);
		bw_lvblock_step_p00(nds);
	}

	bw_lvblock_reduce_separated(nds_orig);
}

int main(int argc, char *argv[])
{
	bw_vector_block d0, d1, nds;
	int ncpus;
	struct opt_desc * opts = NULL;
	int n_opts=0;
	struct fn_cell * lp;
	char	nname[FILENAME_LENGTH];
	int	kval;

	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	/* puts(version_string); */

	new_option(&opts, &n_opts,
		OPT_INTPRM(NULL,&bank_num,OPT_INDEP|OPT_REQUIRED));
	new_option(&opts, &n_opts,
		OPT_INTPRM(NULL,&ncpus,OPT_INDEP|OPT_REQUIRED));
	new_option(&opts,&n_opts,"--inputfile -i",
			1,	PARAM_TEXTUAL(NULL),
			&append_ll_filename,
			&filename_ll,
			OPT_MATCH_NULL|OPT_REQUIRED,
			"sets input file");

	process_options(argc,argv,n_opts,opts);

	if (ncpus==0) {
		showuse();
	}

	set_all_filenames(bank_num);
	read_tag_file();

#ifdef HARDCODE_PARAMS
	consistency_check("bw_allocsize",computed_bw_allocsize,bw_allocsize);
	consistency_check("bw_longsize",computed_bw_longsize,bw_longsize);
#endif

	kval=ncols/ncpus;
	configure_threads(1,ncols);

	printf("The %d first coefficients are straight\n", kval);

	/* no matrix stuff neeeded here */
	/* no bw_x stuff either */


	/* We need three vectors.
	 *
	 * 	d0	the previous dS file
	 * 	d1	the current  dS file
	 * 	nds	the resulting file.
	 *
	 * nds is computed as [ d1[0]...d1[k-1] (d1-d0)[k]...(d1-d0)[N-1]];
	 *
	 * where k = N / ncpus;
	 *
	 * The program proceeds as such
	 *
	 * allocate
	 *
	 * set everybody to zero
	 *
	 * LOOP:	read d1
	 * 		compute nds
	 * 		write nds
	 * 		compute d1 onto d0
	 * 		loop
	 */

	bw_lvblock_alloc(d0);
	bw_lvblock_alloc(d1);
	bw_lvblock_alloc(nds);
	bw_lvblock_set_zero_full(d0);

	for(lp=filename_ll;lp;lp=lp->next) {
		bw_lvblock_set_zero_full(d1);
		bw_lvblock_set_zero_full(nds);
		printf("Reading %s\n",lp->fn);
		read_vec(d1,lp->fn);
		compute_nds(nds,d0,d1,kval);
		nds_name(nname,lp->fn);
		printf("Writing %s\n",nname);
		write_vec(nname,nds);
		memcpy(d0,d1,ncols*bw_longsize*sizeof(mp_limb_t));
	}

	bw_lvblock_free(d0);
	bw_lvblock_free(d1);
	bw_lvblock_free(nds);

	return 0;
}

/* vim:set sw=8: */
