#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "master_params.h"
#include "params.h"
#include <assert.h>
#include "types.h"
#include "macros.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "filenames.h"
#include "options.h"
#include "tagfile.h"
#include "timer.h"
#include "variables.h"
#include "structure.h"
#include "gmp-hacks.h"
#include "modulus_hacks.h"
#include "right_action.h"
#include "e_polynomial.h"
#include "twisting_polynomials.h"
#include "fft_on_matrices.h"
#include "field_def.h"
#include "field_prime.h"
#include "field_quad.h"
#include "field_usage.h"
/* #include "version.h" */

int * degnom;
int * valuation;
int * degree;
int dmax;
unsigned int * clist;
bw_nbpoly f_rev_poly;
int t_counter;

extern void p_nbmat(bw_nbmat);

/* reload f polynomials from iteration inum */
static int
load_f_matrix(void)
{
	int i,j,k;
	struct stat sbuf;
	char filename[FILENAME_LENGTH];
	int d;
	FILE * f;

	/* build the degnom table */

	dmax=-1;
	for(j=0;j<bigdim;j++) {
		degnom[j]=-1;
		d=-1;
		for(i=0;i<n_param;i++) {
			sprintf(filename,f_meta_filename,i,j,t_counter);
			if (stat(filename,&sbuf)<0) {
				perror(filename);
				exit(errno);
			}
			d=(sbuf.st_size/(bw_filesize*sizeof(mp_limb_t)))-1;
			if (degnom[j]!=-1 && degnom[j]!=d) {
				fprintf(stderr,
					"The degree of column %d is unclear\n",
					j);
				exit(1);
			}
			degnom[j]=d;
		}
		if (d>dmax) 
			dmax=d;
	}

	nbpoly_alloc(f_rev_poly,dmax);
	nbpoly_zero(f_rev_poly,dmax);

	for(j=0;j<bigdim;j++) 
	for(i=0;i<n_param;i++) {
		sprintf(filename,f_meta_filename,i,j,t_counter);
		f = fopen(filename,"r");
		if (f == NULL) {
			perror(filename);
			exit(errno);
		}
		for(k=0;k<=degnom[j];k++) {
			bw_scalar_read(nbmat_scal(nbpoly_coeff(
							f_rev_poly,degnom[j]-k),
						i,j),f);
		}
		fclose(f);
	}

	for(j=0;j<bigdim;j++) {
		int val,deg;
		for(val=0;val<=degnom[j];val++) {
			if (!ncol_is_zero(nbmat_col(nbpoly_coeff(f_rev_poly,val),j)))
				break;
		}
		for(deg=degnom[j];deg>=0;deg--) {
			if (!ncol_is_zero(nbmat_col(nbpoly_coeff(f_rev_poly,deg),j)))
				break;
		}
		if (deg<0 || val>degnom[j]) {	/* f_j identically zero */
			printf("Column f_%d is zero, should not happen\n",j);
			eternal_sleep();
		}
		assert(deg>=val);
		printf("Column f_rev_%d (=X^%df_%d(1/X)) : val=%d, deg=%d\n",
				j,degnom[j],j,val,deg);
		valuation[j]=val;
		degree[j]=deg;
	}

	for(j=0;j<bigdim;j++) {
		if (valuation[j]!=0) {
			printf("Use column %d\n",j);
			exit(0);
		}
	}

	return 0;
}

static int stupid_sort(const void * a, const void * b)
{
	return (degnom[*(int*)a]-degnom[*(int*)b]);
}

static int
locate_solution_columns(void)
{
	int t_0;
	int j,k;

	clist=my_malloc(n_param*sizeof(unsigned int));
	t_0=(m_param+n_param-1)/n_param;

	k=0;
	for(j=0;j<bigdim;j++) {
		if (degnom[j]*bigdim<(n_param*t_0+m_param*t_counter)) {
			/* Column is below its expected degree, therefore
			 * it's probably a solution */
			if (k>=n_param) {
				printf("There are more than %d candidates\n",
						n_param);
				exit(1);
			}
			printf("Using column %d\n",j);
			clist[k++]=j;
		}
	}
	if (k<n_param) {
		printf("There are less than %d candidates\n", n_param);
		exit(1);
	}
	qsort(clist,n_param,sizeof(unsigned int),&stupid_sort);
	/* Therefore, the degnoms won't be changed */
	return 0;
}

static int do_gauss_on_f(void)
{
	int i,j,k,l,s;
	bw_scalar piv;
	int rank;
	mp_limb_t * inv;
	mp_limb_t * lambda;
	bw_nbmat f_0;
	unsigned int * zero;
	int nz=0;

	inv	= my_malloc(k_size * sizeof(mp_limb_t));
	lambda	= my_malloc(k_size * sizeof(mp_limb_t));
	zero	= my_malloc(n_param* sizeof(unsigned int));

	f_0=nbpoly_coeff(f_rev_poly,0);

#ifndef NDEBUG
	p_nbmat(f_0);
#endif

	/* Pay attention here, this is a gaussian elimination on
	 * *columns* */

	rank = 0 ;
	for(j = 0 ; j < n_param ; j++) {
		/* Find the pivot inside the column. */
		for(i = 0 ; i < n_param ; i++) {
			piv = nbmat_scal(f_0,i,clist[j]);
			if (!k_is_zero(piv))
				break;
		}
		if (i == n_param) {
			zero[nz++]=j;
			continue;
		}
		rank++;
		k_inv(inv,nbmat_scal(f_0,i,clist[j]));
		/* Cancel this coeff in all other columns. */
		for(k = j + 1 ; k < n_param ; k++) {
			k_mul(lambda,nbmat_scal(f_0,i,clist[k]),inv);
			k_neg(lambda,lambda);
			assert(degnom[j]<=degnom[k]);
			k_set_zero(nbmat_scal(f_0,i,clist[k]));
			for(l=i+1;l<m_param;l++) {
				addmul( nbmat_scal(f_0,l,clist[k]),
					nbmat_scal(f_0,l,clist[j]),
					lambda);
			}
			for(s=1;s<=dmax;s++) for(l=0;l<n_param;l++) {
				addmul( nbmat_scal(nbpoly_coeff(f_rev_poly,s),
							l,clist[k]),
					nbmat_scal(nbpoly_coeff(f_rev_poly,s),
							l,clist[j]),
					lambda);
			}
		}
	}

	free(inv);
	free(lambda);
	
	printf("Rank : %d\n",rank);

	for(j=0;j<nz;j++) {
		char filename[FILENAME_LENGTH];
		int jr;
		FILE *f;

		jr=clist[zero[j]];

		printf("Column %d has valuation increased\n",jr);
		printf("Writing column as column number %d\n",bigdim+j);

		/* Rewrite them in the specified order */
		for(i=0;i<n_param;i++) {
			sprintf(filename,f_meta_filename,i,bigdim+j,t_counter);
			f = fopen(filename,"w");
			if (f == NULL) {
				perror("writing f");
				eternal_sleep();
			}
			printf("Writing %s\n",filename);
			for(k=0;k<=degnom[jr];k++) {
				bw_scalar_write(f,nbmat_scal(
						nbpoly_coeff(f_rev_poly,
								degnom[jr]-k),
							i,jr));
			}
			fclose(f);
		}
		for(k=0;k<=degnom[jr];k++) {
			if (!ncol_is_zero(nbmat_col(nbpoly_coeff(f_rev_poly,
							k), jr)))
				break;
		}
		
		/* k is now the valuation of the (reversed) column. */
	
		sprintf(filename,valu_meta_filename,bigdim+j,t_counter);
		f = fopen(filename,"w");
		if (f == NULL) {
			perror("Writing valuation file");
			eternal_sleep();
		}
		printf("Writing %s\n",filename);
		fprintf(f,"COLUMN %d VALUATION %d DEGNOM %d\n",
				bigdim+j ,k, degnom[jr]);
		fclose(f);
	}
	return 0;
}



int main(int argc, char * argv[])
{
	struct opt_desc * opts = NULL;
	int n_opts=0;
	int bank_num;
	
	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	/* puts(version_string); */

	coredump_limit(0);

	new_option(&opts,&n_opts,
			OPT_INTPRM(NULL,&bank_num,OPT_INDEP|OPT_REQUIRED));
	new_option(&opts,&n_opts,
			OPT_INTPRM(NULL,&t_counter,OPT_INDEP|OPT_REQUIRED));

	process_options(argc,argv,n_opts,opts);

	set_all_filenames(bank_num);

	read_tag_file();

#ifdef HARDCODE_PARAMS
	consistency_check("m_param",computed_m_param,m_param);
	consistency_check("n_param",computed_n_param,n_param);
	consistency_check("bw_allocsize",computed_bw_allocsize,bw_allocsize);
	consistency_check("bw_longsize",computed_bw_longsize,bw_longsize);
#endif

	degnom=my_malloc(bigdim*sizeof(int));
	degree=my_malloc(bigdim*sizeof(int));
	valuation=my_malloc(bigdim*sizeof(int));
	
	field_k=new_prime_field(modulus_plain,bw_allocsize);
	
	load_f_matrix();

	locate_solution_columns();

	do_gauss_on_f();

	return 0;
}
