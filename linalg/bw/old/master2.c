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
#include "timer.h"
#include "variables.h"
#include "structure.h"
#include "gmp-hacks.h"
#include "modulus_hacks.h"
#include "e_polynomial.h"
#include "twisting_polynomials.h"
#include "fft_on_matrices.h"
#include "field_def.h"
#include "field_prime.h"
#include "field_quad.h"
#include "field_usage.h"
/* #include "version.h" */
#include "master_common.h"

bw_nbpoly	f_poly;
int		t_counter;
int	      * global_delta;
int		global_sum_delta;
unsigned int  * chance_list;
static int	recursion_level;	/* static, hence 0 on init */
static int	rec_threshold=1;
static int	print_min=10;
static int	preferred_quadratic_algorithm;	/* Defaults to 0 (old) */
static int	t_init;
static int	no_check_input;
static int	enable_cplx;
static int	fermat_prime;

const char f_base_filename[] = "F_INIT";

#define SAVE_LEVEL_THRESHOLD	4

void reclevel_prolog(void)
{
	int i;
	printf("%2d  ",recursion_level);
	for(i=0;i<recursion_level;i++) printf("  ");
}

int sum_delta(int * delta)
{
	int i,res=0;
	for(i=0;i<bigdim;i++) res += delta[i];
	return res;
}

int max_delta(int * delta)
{
	int i,res=-1;
	for(i=0;i<bigdim;i++) if (delta[i]>res) res=delta[i];
	return res;
}


static int
save_pi(struct t_poly * pi, int t_start, int t_middle, int t_end)
{
	char filename[FILENAME_LENGTH];
	FILE *f;
	int res;

	sprintf(filename, pi_meta_filename, t_start, t_end);
	f=fopen(filename,"w");
	if (f==NULL)
		return -1;
	res=tp_write(f,pi);
	if (fclose(f)<0)
		res=-1;

	if (res==0) {
		printf("Saved file %s\n",filename);
	} else {
		fprintf(stderr,"Failure to save %s: %s\n",
				filename,strerror(errno));
		return -1;
	}

	if (t_middle==-1)
		return res;

#if 0
	/* Unlinking not done, for safety */
	sprintf(filename, pi_meta_filename, t_start, t_middle);

	if (unlink(filename)<0) {
		fprintf(stderr,"Cannot unlink %s: %s\n",
				filename,strerror(errno));
	} else {
		printf("Unlinked %s\n",filename);
	}

	sprintf(filename, pi_meta_filename, t_middle, t_end);

	if (unlink(filename)<0) {
		fprintf(stderr,"Cannot unlink %s: %s\n",
				filename,strerror(errno));
	} else {
		printf("Unlinked %s\n",filename);
	}
#endif



	return res;
}

static void
core_if_null(const void * p, const char * v)
{
	if (p==NULL) {
		printf("Could not allocate space for %s, dumping core\n",v);
		abort();
	}
}


static int
retrieve_pi_files(struct t_poly ** p_pi, int t_start)
{
	DIR * pi_dir;
	struct dirent * curr;
	int n_pi_files;
	struct couple {int s; int e;} * pi_files;
	const char *pattern;
	int i;
	struct t_poly * left = NULL, * right = NULL;
	int order;
	struct dft_bb * dft_left, * dft_right, * dft_prod;
	double tt;

	*p_pi=NULL;

	if ((pi_dir=opendir("."))==NULL) {
		perror(".");
		return t_start;
	}

	printf("Scanning directory %s for pi files\n", ".");

	pattern = strrchr(pi_meta_filename,'/');
	if (pattern == NULL) {
		pattern = pi_meta_filename;
	} else {
		pattern++;
	}

	for(n_pi_files=0;(curr=readdir(pi_dir))!=NULL;) {
		int s,e;
		if (sscanf(curr->d_name,pattern,&s,&e)==2) {
			printf("Found %s\n", curr->d_name);
			if (s>e) {
				printf("but that's a stupid one\n");
				continue;
			}
			n_pi_files++;
		}
	}

	if (n_pi_files==0) {
		printf("Found no pi files\n");
		return t_start;
	}

	pi_files=malloc(n_pi_files*sizeof(struct couple));
	
	rewinddir(pi_dir);

	for(i=0;i<n_pi_files && (curr=readdir(pi_dir))!=NULL;) {
		if (sscanf(curr->d_name,
					pattern,
					&(pi_files[i].s),
					&(pi_files[i].e))==2)
		{
			if (pi_files[i].s>pi_files[i].e)
				continue;
			i++;
		}
	}
	n_pi_files=i;
	closedir(pi_dir);

	/* The rule is: only look at the best candidate. It's not worth
	 * bothering about more subtle cases */
	for(;;) {
		int t_max;
		int best=-1;
		FILE *f;
		char filename[FILENAME_LENGTH];

		printf("Scanning for data starting at t=%d\n",t_start);
		t_max=-1;
		for(i=0;i<n_pi_files;i++) {
			if (pi_files[i].s==t_start && pi_files[i].e != -1) {
				printf("candidate : ");
				printf(pattern,pi_files[i].s,pi_files[i].e);
				printf("\n");
				if (pi_files[i].e>t_max) {
					t_max=pi_files[i].e;
					best=i;
				}
			}
		}
		if (t_max==-1) {
			printf("Could not find such data\n");
			break;
		}

		sprintf(filename,pi_meta_filename,t_start,t_max);
		printf("trying %s\n", filename);
		f=fopen(filename,"r");
		if (f==NULL) {
			perror(filename);
			pi_files[best].e=-1;
			continue;
		}
		/* Which degree can we expect for t_start..t_max ?
		 */

		unsigned int pideg;
		pideg = iceildiv(m_param * (t_max - t_start), bigdim);
		pideg += 10;
		if (t_max > total_work) {
			pideg += t_max - total_work;
		}

		right=tp_read(f, pideg);
		fclose(f);

		if (right==NULL) {
			printf("%s : bad or nonexistent data\n",filename);
			pi_files[best].e=-1;
			continue;
		}

		if (left==NULL) {
			left=right;
			right=NULL;
			t_start=t_max;
			continue;
		}

		printf("Beginning multiplication\n");
		*p_pi=tp_comp_alloc(left,right);
		core_if_null(*p_pi,"*p_pi");

		order=ceil_log2((*p_pi)->degree+1);
		
		dft_left=fft_tp_dft(left,order,&tt);
		printf("DFT(pi_left,%d) : %.2fs\n",order,tt);
		core_if_null(dft_left,"dft_left");

		dft_right=fft_tp_dft(right,order,&tt);
		printf("DFT(pi_right,%d) : %.2fs\n",order,tt);
		core_if_null(dft_right,"dft_right");

		dft_prod=fft_bbb_conv(dft_left,dft_right,&tt);
		printf("CONV(pi_left,pi_right,%d) : %.2fs\n",order,tt);
		core_if_null(dft_prod,"dft_prod");

		fft_tp_invdft(*p_pi,dft_prod,&tt);
		printf("IDFT(pi,%d) : %.2fs\n",order,tt);

		tp_free(left);
		tp_free(right);
		dft_bb_free(dft_left);
		dft_bb_free(dft_right);
		dft_bb_free(dft_prod);
		left=*p_pi;
		right=NULL;
		t_start=t_max;
	}
	free(pi_files);

	*p_pi=left;
	return t_start;
}

/*
 * bw_init
 *
 * fill in the structures with the data available on disk. a(X) is
 * fetched this way. f(X) is chosen at random, in order to make [X^0]e(X)
 * nonsingular. Once this condition is satisfied, e(X) is computed. All
 * further computations will be done with e(X).
 *
 * a(X) can be thrown away as soon as e(X) is computed.
 *
 * At a given point of the algorithm, one can bet that most of the data
 * for e(X) is useless. Of this, the biggest part is likely to be cleanly
 * swapped out.
 *
 */
static struct e_coeff * bw_init(void)
{
	int		  reached_degree;
	int		  i,j,k,s,t,t0,test;
	mpz_t		  blah;
	bw_mnpoly	  a_poly;
	bw_mbpoly	  e_poly;
	gmp_randstate_t	  randstate;
	unsigned long int	  seed;
	unsigned int		* clist;
	struct e_coeff		* res;
	FILE		* f;

#ifdef	CORRECT_STUPID_WOE
	printf("Using A(X) div X in order to consider Y as starting point\n");
	total_work--;
#endif

	printf("Reading scalar data in polynomial ``a''\n");
	mnpoly_alloc(a_poly,total_work);
	mnpoly_zero(a_poly,total_work);
	
	reached_degree = read_data_for_series(a_poly);

	t_counter = t0 = t_init = iceildiv(m_param,n_param);

	if (reached_degree < t_counter) {
		printf( "This amount of data is insufficient. Synchronious\n"
			"revision up to stamp %d is mandatory\n",t_counter);
		exit(1);
	}
	if (reached_degree != total_work) {
		printf( "Since we do not have the full information yet, we\n"
			"merely can tell if the system looks like doable\n");
	}
			
	/* Data read stage completed. */

	printf("Preparing FFT engine\n");
	if (fermat_prime) {
#ifdef	HAS_NATIVE_FFT
		enable_cplx=-1;
#else
		die("To use the --fermat-prime option, "
			"recompile with -DHAS_NATIVE_FFT\n",1);
#endif
	}


	prepare_fft_engine(ceil_log2((total_work<<1)+2),enable_cplx);

	printf("Setting up primary polynomial f (t0=%d)\n",t0);

	f=fopen(f_base_filename,"r");
	test=(f!=NULL);

	if (test) {
		printf("Found data on disk. Restoring\n");
		f_poly=nbpoly_read(f, t0);
		fclose(f);
		if (STRICTTYPE_VAL(f_poly)==NULL) {
			perror(f_base_filename);
			test=0;
		}
	}

	if (!test) {
		nbpoly_alloc(f_poly,t0);
		nbpoly_zero(f_poly,t0);
		gmp_randinit(randstate, GMP_RAND_ALG_DEFAULT, 128);
		/*
		seed=time(NULL);
		*/
		seed=0x3b7b99de;
		gmp_randseed_ui(randstate,seed);
		mpz_init(blah);
		printf("Building f from random data (seed=0x%08lx)\n",seed);
		for(k=0; k < t0  ;k++)
		for(j=0; j < m_param;j++)
		for(i=0; i < n_param;i++) {
			mpz_urandomb(blah,
					randstate,
					bw_allocsize * mp_bits_per_limb);
			/* mpz_random(blah, bw_allocsize); */
			mpz_fdiv_r(blah,blah,modulus);
			MPN_SET_MPZ(nbmat_scal(nbpoly_coeff(f_poly,k),i,j),
					bw_allocsize,blah);
		}
		gmp_randclear(randstate);
		mpz_clear(blah);

		for(i = 0 ; i < subdim ; i++) {
			bw_scalar_set_one(nbmat_scal(nbpoly_coeff(f_poly, t0),
						i, m_param+ i));
		}

		f=fopen(f_base_filename,"w");
		if (f==NULL) {
			perror(f_base_filename);
		} else {
			test=nbpoly_write(f,f_poly,t0);
			fclose(f);
			if (test==0) {
				printf("Written f to %s\n",f_base_filename);
			}
		}
	}

	clist		= malloc(bigdim * sizeof(unsigned int));
	global_delta	= malloc(bigdim * sizeof(int));
	chance_list	= malloc(bigdim * sizeof(unsigned int));
	for(j=0;j<bigdim;j++) {
		clist[j]	= j;
		global_delta[j]	= t_counter;
		chance_list[j]	= 0;
	}
	
	/* Here, the basic setup is done. Unless we are very unlucky, and
	 * ctaf is singular. We will check against that later, since this
	 * is not very likely. (``later'' means : in the course of the
	 * gaussian elimination work).
	 *
	 * In the case of GF(2), the odds are bigger of course. Choice of
	 * another ``f'' polynomial is likely to trim the probabilities
	 * down to a minimum. If this still doesn't work, increasing t0
	 * might help. Otherwise, we're in trouble. */

#ifdef ACTIVATE_HACKS
	printf("Doing some precomputations on the modulus\n");
	do_modulus_precomps();
#else	
	choke me;
#endif

	printf("Computing value of e(X)=a(X)f(X) (degree %d)\n",total_work);
	mbpoly_alloc(e_poly,total_work);
	mbpoly_zero(e_poly,total_work);

	for(t = 0 ; t <= total_work; t++) {
	for(s = 0 ; s <= t0;         s++) {
		if (s+t < t0 || s+t==0 || s+t > total_work)
			continue;
		for(i = 0     ; i <  m_param;    i++) {
		for(j = 0     ; j <  bigdim ;    j++) {
		for(k = 0     ; k <  n_param;    k++) {
			addmul( mbmat_scal(mbpoly_coeff(e_poly,s+t),i,j),
				mnmat_scal(mnpoly_coeff(a_poly,t),i,k),
				nbmat_scal(nbpoly_coeff(f_poly,s),k,j));
		}}}
	}}

	printf("Throwing out a(X)\n");
	mnpoly_free(a_poly);

	res =  ec_encapsulate(mbpoly_subpoly(e_poly,1),clist,total_work-1);
	
	return res;
}

/* It might seem merely cosmetic and useless to sort w.r.t both the
 * global and local nominal degrees. In fact, it is crucial for the
 * corectness of the computations. (Imagine a 2-step increase, starting
 * with uneven global deltas, and hitting an even situation in the
 * middle. One has to sort out the local deltas to prevent trashing the
 * whole picture).
 *
 * The positional sort, however, *is* cosmetic (makes debugging easier).
 */

struct xcol_id {
	unsigned int pos;
	int deg;
	int sdeg;
};

static int
col_cmp(const struct xcol_id * x, const struct xcol_id * y)
{
	int diff;
	int sdiff;
	diff = x->deg - y->deg;
	sdiff = x->sdeg - y->sdeg;
	return diff?diff:(sdiff?sdiff:(x->pos - y ->pos));
}

/* Sort the degree table delta, with local ordering obeying the columns
 * of pi. The permutation is applied to delta, but not to pi. It is
 * returned, as a malloc'ed int*
 */
static void
column_order(unsigned int * perm, int * delta, struct t_poly * pi)
{
	struct xcol_id * tmp_delta;
	int i;

	tmp_delta   = malloc(bigdim  * sizeof(struct xcol_id));

	for(i=0;i<bigdim;i++) {
		tmp_delta[i].pos=i;
		tmp_delta[i].deg=delta[i];
		tmp_delta[i].sdeg=pi->degnom[pi->clist[i]];
	}

	qsort(tmp_delta,bigdim,sizeof(struct xcol_id),(sortfunc_t)&col_cmp);

	for(i=0;i<bigdim;i++) {
		perm[i]=tmp_delta[i].pos;
		delta[i]=tmp_delta[i].deg;
	}

	free(tmp_delta);
}


static void e_transvec(bw_mbmat e,
		unsigned int * clist,
		int o_j1,
		int o_j2,
		int ikill,
		bw_scalar lambda)
{
	int i;
	int j1,j2;

	j1=clist[o_j1];
	j2=clist[o_j2];

	bw_scalar_set_zero(mbmat_scal(e,ikill,j2));
	
	for(i=ikill+1;i<m_param;i++) {
		addmul( mbmat_scal(e,i,j2),
			mbmat_scal(e,i,j1),
			lambda);
	}
}

/*
 * ec : A raw polynomial matrix
 * pi : The permutation applied to previous data to obtain ec.
 * e  : [X^t] ec * pi
 * t  : the unknown above
 * v  : v==0 : transport the transvections automatically on ec
 *      v==1 : don't and recompute extensively e for next time if applicable
 *
 *                    v==1 : [X^t] ec * pi == e (ec==ec_base)
 */
static void
bw_gauss_onestep(bw_mbmat e,
		struct e_coeff * ec,
		struct t_poly * pi,
		int * delta,
		unsigned int * pivlist,
		double * tt)
{
	int i,j,k,jr;
	bw_scalar piv;
	int rank;
	mp_limb_t * inv;
	mp_limb_t * lambda;
	unsigned int * pivots_list;
	struct timeval tv;

	timer_r(&tv,TIMER_SET);

	if (pivlist)
		pivots_list = pivlist;
	else
		pivots_list = malloc(m_param * sizeof(int));
	
	inv	= malloc(k_size * sizeof(mp_limb_t));
	lambda	= malloc(k_size * sizeof(mp_limb_t));

	/* Pay attention here, this is a gaussian elimination on
	 * *columns* */

	rank = 0 ;
	for(j = 0 ; j < bigdim ; j++) {
		jr=pi->clist[j];
		/* Find the pivot inside the column. */
		for(i = 0 ; i < m_param ; i++) {
			piv = mbmat_scal(e,i,jr);
			if (!k_is_zero(piv))
				break;
		}
		if (i == m_param)
			continue;
		assert(rank<m_param);
		pivots_list[rank++] = j;
		k_inv(inv,mbmat_scal(e,i,jr));
		/* Cancel this coeff in all other columns. */
		for(k = j + 1 ; k < bigdim ; k++) {
			k_mul(lambda,mbmat_scal(e,i,pi->clist[k]),inv);
			k_neg(lambda,lambda);
			assert(delta[j]<=delta[k]);
			if (!pivlist)
				ec_transvec(ec,j,k,i,lambda);
			else
				e_transvec(e,pi->clist,j,k,i,lambda);
			tp_transvec(pi,j,k,i,lambda);
		}

	}

	free(inv);
	free(lambda);
	
	if (rank!=m_param) {
		fprintf(stderr,"Duh, rank is not == m !\n");
		exit(1);
	}

	for(j = 0 ; j < m_param ; j++) {
		if (!pivlist)
			ec_x_multiply(ec,pivots_list[j]);
		tp_x_multiply(pi,pivots_list[j]);
		delta[pivots_list[j]]++;
		global_sum_delta++;
	}

	if (!pivlist)
		free(pivots_list);

	*tt=timer_r(&tv,TIMER_ASK);
}

static int
bw_check_chance(bw_mbmat e, unsigned int * clist)
{
	int i,j;
	int maxchance;

	maxchance=0;

	for(j=0;j<bigdim;j++) {
		for(i=0;i<m_param;i++) {
			bw_reduce_short_scalar(mbmat_scal(e,i,j));
		}
	}

	for(j=0;j<bigdim;j++) {
		if (mcol_is_zero(mbmat_col(e,clist?clist[j]:j)))
		{
			if (++chance_list[j] > maxchance)
				maxchance=chance_list[j];
			printf("Column %d happens to be zero ! (%d)\n",
					j,chance_list[j]);
		} else
			chance_list[j]=0;
	}

	return maxchance;
}

static void
banner_traditional(int t, int deg, double inner, double * last)
{
	*last+=inner;
	if ((t+1)%print_min==0 || t==deg) {
		reclevel_prolog();
		printf("t=%d	avg=%.1f	"
			"step:	%.2fs	last %d : %.2fs\n",
			t_counter,((double)global_sum_delta)/bigdim,inner,
			1+(t%print_min),*last);
		*last=0.0;
	}
}

static void
compute_ctaf(bw_mbmat e,
		struct e_coeff * ec,
		struct t_poly * pi,
		int t,
		unsigned int * known_cols,
		double * tt)
{
	int i,j,k,s;
	int * bounds;
	struct timeval tv;

	timer_r(&tv,TIMER_SET);

	bounds=malloc(bigdim*sizeof(int));
	memset(bounds,0,bigdim*sizeof(int));
	if (known_cols) for(j=0;j<m_param;j++) {
		bounds[pi->clist[j]]=-1;
	}
	for(j=0;j<bigdim;j++) {
		if (bounds[j]==-1)
			continue;
		bounds[j]=pi->degnom[j];
		mcol_zero(mbmat_col(e,j));
	}

	for(j = 0 ; j < bigdim ; j++)
	for(s = 0 ; s <= bounds[j] ; s++)
	for(i = 0 ; i <  m_param; i++)
	for(k = 0 ; k <  bigdim; k++)	{
		addmul( mbmat_scal(e,i,j),
			mbmat_scal(mbpoly_coeff(ec->p,t-s),i,k),
			bbmat_scal(bbpoly_coeff(pi->p,s),k,j));
	}
	free(bounds);

	*tt=timer_r(&tv,TIMER_SET);
}

static void
bw_traditional_algo_1(struct e_coeff * ec, int * delta,
		struct t_poly * pi, int check_chance)
{
	unsigned int * perm;
	int t;
	bw_mbmat e;
	unsigned int * pivlist;
	double inner_1,inner_2,last;
	int k;

	perm=malloc(bigdim*sizeof(unsigned int));
	pivlist=malloc(m_param*sizeof(unsigned int));

	mbmat_alloc(e);
	mbmat_zero(e);

	assert(!ec_is_twisted(ec));
	last=0.0;

	for(t=0;t<=ec->degree;t++) {
		compute_ctaf(e,ec,pi,t,t?pivlist:NULL, &inner_1);
		column_order(perm,delta,pi);
		tp_apply_perm(pi,perm);
		if (check_chance)
			bw_check_chance(e,pi->clist);
		bw_gauss_onestep(e,ec,pi,delta,pivlist,	&inner_2);
		t_counter++;
		banner_traditional(t, ec->degree, inner_1 + inner_2, &last);
	}
	printf("DELTA : ( ");
	for(k=0;k<bigdim;k++) printf("%d ",delta[k]);
	printf(")\n");

	ec_advance(ec,ec->degree+1);	/* cosmetic */
	mbmat_free(e);
	free(pivlist);
	free(perm);
}

static void
bw_traditional_algo_2(struct e_coeff * ec, int * delta,
		struct t_poly * pi, int check_chance)
{
	unsigned int * perm;
	int t;
	int deg;
	double inner,last;

	perm=malloc(bigdim*sizeof(unsigned int));

	assert(!ec_is_twisted(ec));

	deg=ec->degree;
	last=0.0;

	for(t=0;t<=deg;t++) {
		if (check_chance)
			bw_check_chance(mbpoly_coeff(ec->p,0),ec->clist);
		column_order(perm,delta,pi);
		tp_apply_perm(pi,perm);
		ec_apply_perm(ec,perm);
		bw_gauss_onestep(mbpoly_coeff(ec->p,0),ec,pi,delta,NULL,&inner);
		ec_advance(ec,1);
		t_counter++;
		banner_traditional(t, deg, inner, &last);
	}

	free(perm);
}

/*
 * Rule : in the following, ec can be trashed at will.
 *
 * Compute pi_left, of degree (ec->degree + 1)*(m/(m+n)), such that ec * pi is
 * divisible by X^(ec->degree + 1) (that is, all coefficients up to
 * degree ec->degree in the product are forced to 0).
 *
 */

/* Forward declaration, it's used by the recursive version */
static double bw_lingen(struct e_coeff *, int *, struct t_poly **);

static double
bw_traditional_algorithm(struct e_coeff * ec,
		int * delta,
		struct t_poly ** p_pi,
		int check_chance)
{
	struct timeval tv;
	double tt;
	int deg;

	timer_r(&tv,TIMER_SET);

	*p_pi=tp_alloc(1 + ec->degree);
	tp_set_ident(*p_pi);

	deg=1+ec->degree;

	if (preferred_quadratic_algorithm==0) {
		bw_traditional_algo_1(ec,delta,*p_pi,check_chance);
	} else {
		bw_traditional_algo_2(ec,delta,*p_pi,check_chance);
	}

	tt=timer_r(&tv,TIMER_ASK);
	printf("constants : c_su=%.4e				# %d\n",
			deg?(tt/(double)(deg*deg)):-1.0,deg);

	return tt;
}

static double
bw_recursive_algorithm(struct e_coeff * ec,
		int * delta,
		struct t_poly ** p_pi)
{
	struct t_poly * pi_left, * pi_right;
	int deg,ldeg,rdeg;
	struct dft_mb *dft_e_left,*dft_e_middle;
	struct dft_bb *dft_pi_left,*dft_pi_right;
	struct dft_bb *dft_pi;
	unsigned int sub_order;
	int expected_pi_deg;
	struct timeval tv;
	double	t_dft_e_l,  t_dft_pi_l, t_conv_e, t_idft_e,
		t_dft_pi_r, t_conv_pi,  t_idft_pi, t_ft, t_cv, t_sub;

	timer_r(&tv,TIMER_SET);

	deg=ec->degree;

	/* Repartition of the job:
	 *
	 *		left		right
	 * deg==0	1		0	(never recursive)
	 * deg==1	1		1
	 * deg==2	2		1
	 * deg==n	n/2 + 1		(n+1)/2
	 * 
	 * The figures are for the number of steps, each one corres-
	 * ponding to a m/(m+n) increase of the average degree of pi.
	 */

	ldeg=(deg   /2)+1;
	rdeg=(deg+1)/2;

	assert(ldeg && rdeg && ldeg + rdeg == deg + 1);
	
	/* We aim at computing ec * pi / X^ldeg. The degree of this
	 * product will be
	 *
	 * ec->degree + pi->degree - ldeg
	 *
	 * (We are actually only interested in the low (ec->degree-ldeg)
 	 * degree part of the product, but the whole thing is required)
	 *
	 * The expected value of pi->degree is 
	 * 	ceil(ldeg*m/(m+n))
	 *
	 * The probability that pi exceeds this expected degree
	 * depends on the base field, but is actually low.
	 * However, by the end of the computations, this does
	 * happen because the degrees increase unevenly.
	 *
	 * The DFTs of e and pi can be computed using only the
	 * number of points given above, *even if their actual
	 * degree is higher*. The FFT routines need to have
	 * provision for this.
	 *
	 * The number of points will then be the smallest power
	 * of 2 above deg+ceil(ldeg*m/(m+n))-ldeg+1
	 */
	
	expected_pi_deg = iceildiv(ldeg*m_param, bigdim);
	sub_order=ceil_log2(deg+expected_pi_deg-ldeg+1);

	dft_e_left	= fft_ec_dft(ec,sub_order,&t_dft_e_l);
	reclevel_prolog();
	printf("DFT(e,%d) : %.2fs\n", sub_order, t_dft_e_l);
	core_if_null(dft_e_left,"dft_e_left");

	ec->degree	= ldeg - 1;
	t_sub=bw_lingen(ec,delta,&pi_left);

	dft_pi_left	= fft_tp_dft(pi_left,sub_order,&t_dft_pi_l);
	reclevel_prolog();
	printf("DFT(pi_l,%d) : %.2fs\n", sub_order, t_dft_pi_l);
	core_if_null(dft_pi_left,"dft_pi_left");

	printf("deg(pi_l)=%d, bound is %d\n",pi_left->degree,expected_pi_deg);
	if ((1<<sub_order) < deg+pi_left->degree-ldeg + 1) {
		printf("Argl. pi grows above its expected degree...\n");
		printf("%d is the bound, while :\n"
				"deg=%d\n"
				"deg(pi_left)=%d\n"
				"ldeg-1=%d\n"
				"hence, %d is too big\n",
				(1<<sub_order),deg,pi_left->degree,ldeg - 1,
				deg+pi_left->degree-ldeg + 1);
		eternal_sleep();
	}

	dft_e_middle = fft_mbb_conv_sp(dft_e_left,dft_pi_left,ldeg,&t_conv_e);
	reclevel_prolog();
	printf("CONV(e*pi_l,%d) : %.2fs\n", sub_order, t_conv_e);
	core_if_null(dft_e_middle,"dft_e_middle");

	/* This is a special convolution in the sense that we
	 * compute f(w)*g(w) / w^k for k=ldeg, since we are
	 * interested in fg div X^k (we know fg mod X^k==0)
	 */
	
#if 0
	ec_park(ec);				/* moderately useful... */
#endif
	ec_untwist(ec);
	fft_mb_invdft(ec->p,dft_e_middle,deg-ldeg,&t_idft_e);
	reclevel_prolog();
	printf("IDFT(e,%d) : %.2fs\n", dft_e_middle->order, t_idft_e);

	ec->degree=deg-ldeg;

	dft_mb_free(dft_e_middle);
	dft_mb_free(dft_e_left);

	assert(ec->degree==rdeg-1);

	t_sub+=bw_lingen(ec,delta,&pi_right);
	printf("deg(pi_r)=%d, bound is %d\n",pi_right->degree,expected_pi_deg);

	*p_pi		= tp_comp_alloc(pi_left,pi_right);
	core_if_null(*p_pi,"*p_pi");
	
	printf("deg(pi_prod)=%d (max=%d)\n",(*p_pi)->degree, (1<<sub_order));
	assert((*p_pi)->degree < (1<<sub_order));
	
	dft_pi_right	= fft_tp_dft(pi_right, sub_order, &t_dft_pi_r);
	reclevel_prolog();
	printf("DFT(pi_r,%d) : %.2fs\n", sub_order, t_dft_pi_r);
	core_if_null(dft_pi_right,"dft_pi_right");

	dft_pi		= fft_bbb_conv(dft_pi_left,dft_pi_right, &t_conv_pi);
	reclevel_prolog();
	printf("CONV(pi_l*pi_r,%d) : %.2fs\n", sub_order, t_conv_pi);
	core_if_null(dft_pi,"dft_pi");


	fft_tp_invdft(*p_pi,dft_pi,&t_idft_pi);
	reclevel_prolog();
	printf("IDFT(pi,%d) : %.2fs\n",dft_pi->order,t_idft_pi);

	dft_bb_free(dft_pi);
	dft_bb_free(dft_pi_right);
	dft_bb_free(dft_pi_left);
	tp_free(pi_left);
	tp_free(pi_right);

	reclevel_prolog();

	t_ft=t_dft_e_l+t_dft_pi_l+t_idft_e+t_dft_pi_r+t_idft_pi;
	t_cv=t_conv_e+t_conv_pi;

	printf("proper : %.2fs (%.2fs FT + %.2fs CV), sub : %.2fs\n",
			t_ft+t_cv,t_ft,t_cv,t_sub);
	printf("constants : c_ft=%.4e c_cv=%.4e		# %d,%d\n",
			t_ft/(double)(sub_order<<sub_order),
			t_cv/(double)(1<<sub_order),
			deg,sub_order);
	printf("Different values for M1:");
	printf("   e_left: M1=%.3e\n",
			t_dft_e_l/ (sub_order<<sub_order)/(m_param*bigdim));
	printf("  pi_left: M1=%.3e\n",
			t_dft_pi_l/(sub_order<<sub_order)/(bigdim*bigdim));
	printf("    e_inv: M1=%.3e\n",
			t_idft_e/  (sub_order<<sub_order)/(m_param*bigdim));
	printf(" pi_right: M1=%.3e\n",
			t_dft_pi_r/(sub_order<<sub_order)/(bigdim*bigdim));
	printf("   pi_inv: M1=%.3e\n",
			t_idft_pi/ (sub_order<<sub_order)/(bigdim*bigdim));
	printf("   e_conv: M1=%.3e\n",
			t_conv_e/ (1<<sub_order)  /(m_param*bigdim*bigdim));
	printf("  pi_conv: M1=%.3e\n",
			t_conv_pi/ (1<<sub_order) / (bigdim*bigdim*bigdim));
	return timer_r(&tv,TIMER_ASK);
}

static double
bw_lingen(struct e_coeff * ec,
		int * delta,
		struct t_poly ** p_pi)
{
	int check_chance;
	double inner;
	int deg;
	int t_before;
	int did_rec=0;

	t_before=t_counter;
	deg=ec->degree;
	reclevel_prolog();
	printf("Degree %d\n",deg);

	recursion_level++;
	check_chance=(t_counter > total_work - 5);
	if (ec->degree < rec_threshold || check_chance) {
		check_chance=(t_counter > total_work - 5 - ec->degree);
		inner=bw_traditional_algorithm(ec,delta,p_pi,check_chance);
	} else {
		did_rec=1;
		inner=bw_recursive_algorithm(ec,delta,p_pi);
	}
	recursion_level--;

	reclevel_prolog();
	printf("Degree %d : took %.2fs\n",deg,inner);

	if (recursion_level <= SAVE_LEVEL_THRESHOLD) {
		if (!did_rec || recursion_level == SAVE_LEVEL_THRESHOLD) {
			save_pi(*p_pi,t_before,-1,t_counter);
		} else {
			int t_middle;
			t_middle=t_before+(deg/2)+1;
			save_pi(*p_pi,t_before,t_middle,t_counter);
		}
	}

	return inner;
}

void showuse(void)
{
	die("Usage : bw-master <bank#>\n",1);
}

void compute_f_final(struct t_poly * pi_prod)
{
	bw_nbpoly original_f;
	int i,j,k,s,t,n;

	printf("Computing value of f(X)=f0(X)pi(X) (degree %d)\n",total_work);

	original_f=f_poly;
	n=t_init+pi_prod->degree;

	nbpoly_alloc(f_poly,n);
	nbpoly_zero(f_poly,n);

	for(s = 0 ; s <= t_init     ; s++) {
	for(t = 0 ; t <= pi_prod->degree ; t++) {
		for(i = 0     ; i <  n_param;    i++) {
		for(j = 0     ; j <  bigdim ;    j++) {
		for(k = 0     ; k <  bigdim;    k++) {
		if (t<=pi_prod->degnom[pi_prod->clist[j]])
			addmul( nbmat_scal(nbpoly_coeff(f_poly,s+t),i,j),
				nbmat_scal(nbpoly_coeff(original_f,s),i,k),
				bbmat_scal(bbpoly_coeff(pi_prod->p,t),k,
					pi_prod->clist[j]));
		}}}
	}}
}

void
block_wiedemann(void)
{
	struct e_coeff * ec;
	struct t_poly * pi_left, * pi_right, *pi_prod;
	int j, t_start, new_t;
	double tt;

	ec=bw_init();
	pi_left=pi_right=0;

	t_start=t_counter;
	new_t=retrieve_pi_files(&pi_left,t_counter);

	printf("ec->degree=%d, new_t=%d, t_counter=%d\n",
			ec->degree,new_t,t_counter);

	if (pi_left!=NULL && !(no_check_input && ec->degree<new_t-t_counter)) {
		int dg_kill;
		struct dft_bb * dft_pi_left;
		struct dft_mb * dft_e_left, * dft_e_middle;
		int order;

		printf("Beginning multiplication (e_left*pi_left)\n");
		
		/* That's a pity to bring the DFT so far, but we have to
		 * check the input is correct */

		dg_kill=new_t-t_counter;

		if (no_check_input==0) {
			order=ceil_log2(1+pi_left->degree+ec->degree);
		} else {
			order=ceil_log2(1+pi_left->degree+ec->degree-dg_kill);
		}
		
		dft_e_left=fft_ec_dft(ec,order,&tt);
		printf("DFT(e_left,%d) : %.2fs\n",order,tt);
		core_if_null(dft_e_left,"dft_e_left");

		dft_pi_left=fft_tp_dft(pi_left,order,&tt);
		printf("DFT(pi_left,%d) : %.2fs\n",order,tt);
		core_if_null(dft_pi_left,"dft_pi_left");

		ec_untwist(ec);
		if (no_check_input==0) {
			dft_e_middle=fft_mbb_conv_sp(dft_e_left,dft_pi_left,
					0,&tt);
		} else {
			dft_e_middle=fft_mbb_conv_sp(dft_e_left,dft_pi_left,
					dg_kill,&tt);
			ec_advance(ec,dg_kill);
		}
		printf("CONV(e_left,pi_left,%d) : %.2fs\n",order,tt);
		core_if_null(dft_e_middle,"dft_e_middle");

		if (no_check_input==0) {
			fft_mb_invdft(ec->p,dft_e_middle,
					ec->degree,&tt);
			printf("IDFT(e,%d) : %.2fs\n",order,tt);
			printf("Verifying the product\n");
			for(;t_counter<new_t;t_counter++) {
				int res=1;
				for(j=0;res && j<bigdim;j++) {
					res = res && mcol_is_zero(
							mbmat_col(
							mbpoly_coeff(ec->p,0),
							j));
				}
				if (!res) {
					die("Incorrect input\n",1);
				}
				ec_advance(ec,1);
			}
			printf("Input OK\n");
		} else {
			fft_mb_invdft(ec->p,dft_e_middle,
					ec->degree-dg_kill,&tt);
			printf("IDFT(e,%d) : %.2fs\n",order,tt);
		}
		tp_act_on_delta(pi_left,global_delta);

		dft_bb_free(dft_pi_left);
		dft_mb_free(dft_e_left);
		dft_mb_free(dft_e_middle);
		t_counter=new_t;

		save_pi(pi_left,t_start,-1,t_counter);
	}
	if (pi_left!=NULL && no_check_input && ec->degree<new_t-t_counter) {
		printf("We are not interested in the computation of e(X)\n");
		ec_advance(ec,new_t-t_counter);
		tp_act_on_delta(pi_left,global_delta);
		t_counter=new_t;
		save_pi(pi_left,t_start,-1,t_counter);
	}

	global_sum_delta=sum_delta(global_delta);

	if (ec->degree>=0) {
		bw_lingen(ec,global_delta,&pi_right);
	} else {
		pi_right=pi_left;
		pi_left=NULL;
	}


	if (pi_left!=NULL) {
		struct dft_bb * dft_pi_left, * dft_pi_right, * dft_pi_prod;
		int order;

		printf("Beginning multiplication (pi_left*pi_right)\n");
		pi_prod=tp_comp_alloc(pi_left,pi_right);
		core_if_null(pi_prod,"pi_prod");

		order=ceil_log2(pi_prod->degree+1);
		
		dft_pi_left=fft_tp_dft(pi_left,order,&tt);
		printf("DFT(pi_left,%d) : %.2fs\n",order,tt);
		core_if_null(dft_pi_left,"dft_pi_left");

		dft_pi_right=fft_tp_dft(pi_right,order,&tt);
		printf("DFT(pi_right,%d) : %.2fs\n",order,tt);
		core_if_null(dft_pi_right,"dft_pi_right");

		dft_pi_prod=fft_bbb_conv(dft_pi_left,dft_pi_right,&tt);
		printf("CONV(pi_left,pi_right,%d) : %.2fs\n",order,tt);
		core_if_null(dft_pi_prod,"dft_pi_prod");

		fft_tp_invdft(pi_prod,dft_pi_prod,&tt);
		printf("IDFT(pi,%d) : %.2fs\n",order,tt);

		tp_free(pi_left);
		tp_free(pi_right);
		dft_bb_free(dft_pi_left);
		dft_bb_free(dft_pi_right);
		dft_bb_free(dft_pi_prod);
		pi_left=NULL;
		pi_right=NULL;
		save_pi(pi_prod,t_start,new_t,t_counter);
		/* new_t is the inner value that has to be discarded */
	} else {
		pi_prod=pi_right;
		/* Don't save here, since it's already been saved by
		 * lingen (or it comes from the disk anyway).
		 */
	}

	compute_f_final(pi_prod);

	bw_commit_f(f_poly, global_delta);

	print_chance_list(total_work, chance_list);

	/* Now I can clean up everything if I want... if I want... */
}

void usage()
{
	die("Usage: ./master2 -t <n> [options]\n",1);
}

int
main(int argc, char *argv[])
{
	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	argv++, argc--;

	int pop = 0;

	for( ; argc ; ) {
		if (strcmp(argv[0], "-t") == 0) {
			if (argc <= 1) usage();
			argv++, argc--;
			rec_threshold = atoi(argv[0]);
			argv++, argc--;
			continue;
		}
		if (strcmp(argv[0], "-p") == 0) {
			if (argc <= 1) usage();
			argv++, argc--;
			print_min = atoi(argv[0]);
			argv++, argc--;
			continue;
		}
		if (strcmp(argv[0], "--subdir") == 0) {
			if (argc <= 1) usage();
			int rc = chdir(argv[1]);
			if (rc < 0) {
				perror(argv[1]);
				exit(errno);
			}
			argv+=2;
			argc-=2;
			continue;
		}
		if (strcmp(argv[0], "--enable-complex-field") == 0) {
			enable_cplx = 1;
			argv++, argc--;
			continue;
		}
		if (strcmp(argv[0], "--fermat-prime") == 0) {
			fermat_prime = 1;
			argv++, argc--;
			continue;
		}
		if (strcmp(argv[0], "--no-check-input") == 0) {
			no_check_input = 1;
			argv++, argc--;
			continue;
		}
		if (strcmp(argv[0], "-q") == 0) {
			if (argc <= 1) usage();
			argv++, argc--;
			preferred_quadratic_algorithm = atoi(argv[0]);
			argv++, argc--;
			continue;
			/* changes the quadratic algorithm employed
		{"coppersmith old uneven",0},
		{"new even",1},
			 */
		}

		if (pop == 0) {
			read_mat_file_header(argv[0]);
			pop++;
			argv++, argc--;
			continue;
		} else if (pop == 1) {
			m_param = atoi(argv[0]);
			pop++;
			argv++, argc--;
			continue;
		} else if (pop == 2) {
			n_param = atoi(argv[0]);
			pop++;
			argv++, argc--;
			continue;
		}

		usage();
	}

	if (m_param == 0 || n_param == 0) {
		usage();
	}

	total_work = Lmacro(nrows, m_param, n_param);

	if (rec_threshold == 1)
		usage();

	coredump_limit(1);

#ifdef HARDCODE_PARAMS
	consistency_check("m_param",computed_m_param,m_param);
	consistency_check("n_param",computed_n_param,n_param);
	consistency_check("bw_allocsize",computed_bw_allocsize,bw_allocsize);
	consistency_check("bw_longsize",computed_bw_longsize,bw_longsize);
#endif

	/********************************************************************/

	block_wiedemann();

	return 0;
}
