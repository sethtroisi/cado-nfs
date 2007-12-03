#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef ENABLE_PTHREADS
#include <pthread.h>
#endif

#include "master_params.h"
#include "params.h"
#include <assert.h>
#include "macros.h"
#include "types.h"
#include "endian.h"
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
/* #include "version.h" */
#include "barrier.h"


#define	xxxVERY_VERBOSE
/* Warning ! ; the (pos) field here is only used from within the gauss
 * stuff, in order to have pivot columns virtually first. It is
 * meaningless to other parts of the program (the columns being
 * interchangeable with each other anyway).
 */

struct col_id {
	int degnom;
	int pos;
};

struct bw_context {
	bw_mnpoly	a;
	int		t0;
	int		cur_deg;	/* actually the number of coeffs
					   available ; that is, the degree
					   reached + 1 */
};

struct bw_iterator {
	bw_nbpoly	f;
	bw_mbmat	ctaf;
	int		t;
	struct col_id * clist;
	int	      * pivots_list;
	int	      * chance_list;
};

static struct bw_context	*cx;
static struct bw_iterator	*it;

#define ACTIVATE_HACKS
#define OPT_POWEROF2
#define SKIP_KNOWNCOLS
#define BOUND_F 

#if 0
static void
commit_nbpoly_col(bw_nbpoly poly, int j, int deg, int date)
{
	int i;
	int t;
	char filename[FILENAME_LENGTH];
	FILE *f;
	
	for(i=0;i<n_param;i++) {
		sprintf(filename,f_meta_filename,i,j,date);
		f = fopen(filename,"w");
		if (f == NULL) {
			perror("writing f");
			eternal_sleep();
		}
		printf("writing %s\n",filename);
		for(t=0;t<=deg;t++) {
			bw_scalar_write(f,nbmat_scal(
				nbpoly_coeff(poly,deg-t),i,j));
		}
		/* Pour etre cohérent avec ce que je fais dans master2,
		 * j'écris les coeffs en sens inverse. Comme ça bw-gather
		 * se porte mieux au final...
		 */
		fclose(f);
	}
}
#endif

static void
read_nbpoly_col(bw_nbpoly poly, int j, int deg, int date)
{
	int i;
	int t;
	char filename[FILENAME_LENGTH];
	FILE *f;
	for(i=0;i<n_param;i++) {
		sprintf(filename,f_meta_filename,i,j,date);
		f = fopen(filename,"r");
		if (f == NULL) {
			fprintf(stderr,"read(%s) : %s\n",
					filename,strerror(errno));
			eternal_sleep();
		}
		printf("reading %s\n",filename);
		for(t=0;t<=deg;t++) {
			bw_scalar_read(nbmat_scal(nbpoly_coeff(poly,t),i,j),f);
		}
		fclose(f);
	}
}
	
typedef int (*sortfunc_t)(const void*, const void*);

static void
bw_commit_f(struct bw_iterator * it)
{
	int i;
	int j;
	int t;
	char filename[FILENAME_LENGTH];
	FILE *f;
	int * delta;

	delta=my_malloc(bigdim*sizeof(int));
	for(j=0;j<bigdim;j++) {
		delta[it->clist[j].pos]=it->clist[j].degnom;
	}

	for(j=0;j<bigdim;j++) {
		int deg=delta[j];
		for(i=0;i<n_param;i++) {
			sprintf(filename,f_meta_filename,i,j,it->t);
			f = fopen(filename,"w");
			if (f == NULL) {
				perror("writing f");
				eternal_sleep();
			}
			printf("Writing %s\n",filename);
			for(t=0;t<=deg;t++) {
				bw_scalar_write(f,nbmat_scal(
					nbpoly_coeff(it->f,deg-t), i,j));
			}
			/* Je les écris dans le sens inverse. Pour une
			 * fois... En fait, c'est plus cohérent pour
			 * l'utilisation qui suit, vu que l'évaluation
			 * bête-et-conne n'est pas moins rapide que par
			 * Horner */
			fclose(f);
		}
		for(t=0;;t++) {
			if (!ncol_is_zero(nbmat_col(nbpoly_coeff(it->f,
							deg-t), j)))
				break;
		}
		
		/* t is now the valuation of the (reversed) column. */
	
		sprintf(filename,valu_meta_filename,j,it->t);
		f = fopen(filename,"w");
		if (f == NULL) {
			perror("Writing valuation file");
			eternal_sleep();
		}
		printf("Writing %s\n",filename);
		fprintf(f,"COLUMN %d VALUATION %d DEGNOM %d HAPPY %d TIMES\n",
				j,t, deg, it->chance_list[j]);
		fclose(f);
	}

	free(delta);

	/*
	printf("#RESULT T=%d\n", it->t);
	for(j=0;j<bigdim;j++) {
		if (it->chance_list[j]) printf("#RESULT J=%d\n", j);
	}
	*/
	/* 20060420 : in sync with newer code */
	printf("// step %d LOOK [", it->t);
	for(j=0;j<bigdim;j++) {
		if (it->chance_list[j]) printf(" %d", j);
	}
	printf(" ]\n");
}

#ifdef	VERY_VERBOSE
int bw_dump_f(FILE *f, struct bw_iterator * it)
{
	int i,j,k;
	mp_limb_t * s;
	bw_nbpoly p = it->f;
	char * str;
	int * delta;

	delta=my_malloc(bigdim*sizeof(int));
	for(j=0;j<bigdim;j++) {
		delta[it->clist[j].pos]=it->clist[j].degnom;
	}


	s=malloc(bw_allocsize*sizeof(mp_limb_t));
	str=malloc(bw_allocsize*sizeof(mp_limb_t)*3+2);

	for(i=0;i<n_param;i++) {
		for(j=0;j<bigdim;j++) {
			int deg=delta[j];
			for(k=0;k<=deg;k++) {
				int l,ii,n;
				bw_scalar c=nbmat_scal(nbpoly_coeff(p,k),i,j);
				memcpy(s,c,bw_allocsize*sizeof(mp_limb_t));
				n=mpn_get_str(str,10,s,bw_allocsize);
				if (n==0) str[n++]=0;
				str[n]=0;
				for(ii=0;ii<n-1 && str[ii]==0;ii++);
				for(l=ii;l<n;l++) {
					str[l]+='0';
				}
				if (k==0) {
					fprintf(f,"%s",str+ii);
				} else {
					if (k==1)
						fprintf(f,"+%s*X", str+ii);
					else
						fprintf(f,"+%s*X^%d",str+ii,k);
				}

			}
			fprintf(f,"%s",
				(j!=(bigdim-1))?", ":(
					(i==(n_param-1))?"\n":",\n"));
		}
	}
	free(delta);
	free(s);
	free(str);
	return 0;
}
#endif

/* reload f polynomials from iteration inum */
static int
bw_reload_f(struct bw_iterator * it, int inum)
{
	int j;

	for(j=0;j<bigdim;j++) {
		 char filename[FILENAME_LENGTH];
		 FILE * f;
		 sprintf(filename,valu_meta_filename,j,inum);
		 f = fopen(filename,"r");
		 if (f == NULL) {
			 fprintf(stderr,"read(%s) : %s\n",
					 filename,strerror(errno));
			 eternal_sleep();
		 }
		 printf("Reading %s\n",filename);
		 fscanf(f,"COLUMN %*d VALUATION %*d DEGNOM %d\n",
				 &(it->clist[j].degnom));
		 it->clist[j].pos=j;
		 fclose(f);
		 read_nbpoly_col(it->f,j,it->clist[j].degnom,inum);
	}
	return 0;
}

/* do_ctaf.
 *
 * Most of the time is spent in this functions. Prior to 20010815, there
 * used to be a simple-and stupid version of this function that switched
 * back and forth to mpz, and did not use the already known f's.
 *
 * We know replace this function by a slightly parallelized variant...
 *
 */

static int * ctaf_bounds;
static int want_threads=1;
#ifdef ENABLE_PTHREADS
static barrier_t ctaf_barrier;
static pthread_t * secondaries;
static int nthreads=1;

void spawn_threads(void *(*f)(void *))
{
	int i;
	int blah=PTHREAD_CANCEL_ENABLE;

	if (nthreads==want_threads) {
		(*f)(NULL);
	} else {
		printf("Spawning %d helper threads for ctaf\n",want_threads-1);
		nthreads=want_threads;
		barrier_init(&ctaf_barrier,want_threads,&blah);
		secondaries=my_malloc(want_threads*sizeof(pthread_t));
		for(i=1;i<want_threads;i++) {
			pthread_create(&secondaries[i],NULL,f,NULL);
		}
		secondaries[0]=pthread_self();
		(*f)(NULL);
	}
}

void stop_threads(void)
{
	int i;
	int err;

	for(i=1;i<want_threads;i++) {
		err=pthread_cancel(secondaries[i]);
		fprintf(stderr,"Canceling thread %d : %s\n",i,strerror(err));
	}
	/* The cleanup handlers inside the barrier ensure that the mutex
	 * is released once the threads takes the path out of the
	 * barrier. Therefore, the barrier_destroy should never return
	 * EBUSY, neither should it wait longer than reasonable.
	 */
	barrier_destroy(&ctaf_barrier);
	nthreads=1;
	free(secondaries);
	free(ctaf_bounds);
}
#else
void spawn_threads(void *(*f)(void *))
{
	(*f)(NULL);
}
void stop_threads(void)
{
}
#endif
			
static void
ctaf_column(int j, int bound)
{
	int i,k,t;
	bw_mcol c;

	mcol_alloc(c);
	mcol_zero(c);

	for(t = 0 ; t <= bound ; t++) {
	for(i = 0 ; i <  m_param; i++)	{
	for(k = 0 ; k <  n_param; k++)	{
		addmul( mcol_scal(c,i),
			mnmat_scal(mnpoly_coeff(cx->a,it->t-t),i,k),
			nbmat_scal(nbpoly_coeff(it->f,t),k,j));
	}}}

	mcol_copy(mbmat_col(it->ctaf,j),c);
	mcol_free(c);
}

static void *
bw_do_ctaf_inner(void * a)
{
	int co,j;
#ifdef ENABLE_THREADS
	double twait;
	static double maxwait=0.0;

secondary_thread:
	barrier_wait(&ctaf_barrier,&twait,NULL);
	if (twait>2*maxwait) {
		printf("Thread 0x%lx : waited %.2fs\n",
				(long int)pthread_self(),twait);
		maxwait=twait;
	}
#endif

	for(j = co = 0 ; j <  bigdim ; j++)	{
		if (ctaf_bounds[j]==-1)
			continue;
#ifdef ENABLE_THREADS
		if (pthread_self()==secondaries[co%nthreads]) {
#if 0
			printf("COL%d (L=%d) -> CPU%d\n",
					j,ctaf_bounds[j],co%nthreads);
#endif
			ctaf_column(j,ctaf_bounds[j]);
		}
#else
		ctaf_column(j,ctaf_bounds[j]);
#endif
		co++;
	}
	
#ifdef ENABLE_THREADS
	barrier_wait(&ctaf_barrier,&twait,NULL);
	if (twait>2*maxwait) {
		printf("Thread 0x%lx : waited %.2fs\n",
				(long int)pthread_self(),twait);
		maxwait=twait;
	}

	/* Now, unless I am the master thread, I step back, waiting to
	 * redo the same job.
	 */
	if (pthread_self()!=secondaries[0])
		goto secondary_thread;
#endif

	return NULL;
}

static void
bw_do_ctaf(void)
{
	int rj, j;
	static int not_first=0;


	if (not_first++) {
		memset(ctaf_bounds,0,bigdim*sizeof(int));
		for(rj = 0 ; rj < m_param ; rj++) {
			j=it->clist[it->pivots_list[rj]].pos;
			ctaf_bounds[j]=-1;
		}
	} else {
		ctaf_bounds=my_malloc(bigdim*sizeof(int));
		memset(ctaf_bounds,0,bigdim*sizeof(int));
	}

	for(rj = 0 ; rj <  bigdim ; rj++)	{
		j=it->clist[rj].pos;
		if (ctaf_bounds[j]==0) {
			ctaf_bounds[j]=it->clist[rj].degnom;
			mcol_zero(mbmat_col(it->ctaf,j));
		}
	}

	spawn_threads(&bw_do_ctaf_inner);
}

/* Semantics of current, asked, etc... degrees : number of coeffs,
 * actually. That is, when calling read_data_for_series(a,3,5,...), we
 * know [X^0]...[X^2], and we want [X^3],[X^4].
 */
static int
read_data_for_series(bw_mnpoly a, int cur_deg, int ask_deg, int blocking)
{
	int	i,j,k;
	FILE	*a_disk=NULL;
	char	filename[FILENAME_LENGTH];
	mpz_t	blah;
	int	avail_rev,minrev,askrev,maxrev;
	int	reached_deg;

	mpz_init(blah);

retry_read:

	/* Find out the highest revision commonly available */

	minrev=(MAX(cur_deg-1,0))/periodicity;
	askrev=(MAX(ask_deg-1,0))/periodicity;
	maxrev=(total_work-1)/periodicity;	/* total_work > 0 anyway */
	maxrev+=2;	/* There's a final revision increase due to the
			   last commit. It doesn't hurt to increase, anyway. */

	avail_rev=maxrev;
	reached_deg=total_work;

	for(i=0;i<m_param;i++) for(j=0;j<n_param;j++) {
		int rev;

		for(rev=maxrev;rev>=minrev;rev--) {
			struct stat sbuf;
			sprintf(filename,a_meta_filename,i,j,rev);
			if (stat(filename,&sbuf)<0 || access(filename,R_OK)<0)
				continue;
			reached_deg=MIN(reached_deg,
					sbuf.st_size/(bw_filesize*sizeof(mp_limb_t)));
			break;
		}
		if (rev < avail_rev) {
			avail_rev=rev;
		}
	}

	if (reached_deg==cur_deg && blocking) {
		/* Don't bother try to advance if the revision hasn't
		 * even increased. Don't clutter the display with
		 * messages either */
		sleep(60);
		goto retry_read;
	}

	if (avail_rev<minrev) {
		die("We don't even have revision %d available (availrev=%d, askrev=%d, maxrev=%d, curdeg=%d, askdeg=%d) !\n",63,minrev,avail_rev, askrev,maxrev,cur_deg, ask_deg);
	} else if (avail_rev<askrev) {
		printf("Highest revision available: %d (demanded %d)\n",
				avail_rev,askrev);
	}

	assert(reached_deg>=cur_deg);

	minrev=avail_rev;

	
	/* Now try to actually read the data, up to the highest available
	 * amount. */

	for(i=0;i<m_param;i++) for(j=0;j<n_param;j++) {
		int rev;

		for(rev=minrev;rev<=maxrev;rev++) {
			sprintf(filename,a_meta_filename,i,j,rev);
			a_disk=fopen(filename,"r");
			if (a_disk) break;
			if (errno == ENOENT) continue;
			die("fopen(%s) : %s",errno,filename,strerror(errno));
		}

		if (rev>maxrev)
			die("Data for (%d,%d) disappeared\n",1,i,j);

		printf("Reading file %s\n",filename);

#ifdef	CORRECT_STUPID_WOE
		fseek(a_disk,bw_filesize*sizeof(mp_limb_t),SEEK_SET);
#endif
		fseek(a_disk,cur_deg*bw_filesize*sizeof(mp_limb_t),SEEK_CUR);

		for(k=cur_deg;k<reached_deg && !feof(a_disk);k++) {
			bw_scalar_read(mnmat_scal(mnpoly_coeff(a,k),i,j),
					a_disk);
			MPZ_SET_MPN(blah,mnmat_scal(mnpoly_coeff(a,k),i,j),
					bw_allocsize);
			mpz_fdiv_r(blah,blah,modulus);
			MPN_SET_MPZ(mnmat_scal(mnpoly_coeff(a,k),i,j),
					bw_allocsize,blah);
		}
		if (k<reached_deg) {
			printf("File %s is shorter than expected\n",filename);
			printf("Degree truncated from %d to %d\n",
					reached_deg,k);
			reached_deg=k;
		}
		fclose(a_disk);
	}

	cur_deg=reached_deg;

	if (cur_deg < ask_deg) {
	       if (blocking) {
		       printf("We only have degree %d, "
				"while %d is demanded. Waiting\n",
				cur_deg,ask_deg);
		       sleep(60);	/* one minute */
		       goto	retry_read;
	       } else {
		       printf("Degree reached: %d (wanted %d)\n",
			       cur_deg,ask_deg);
	       }
	}

	mpz_clear(blah);

	return reached_deg;
}
	
#if 0
/* This is the randomized init. code
 */
static void
bw_init(struct bw_iterator * it, struct bw_context * cx, int start)
{
	int	i,j,k;
	mpz_t	blah;

#ifdef	CORRECT_STUPID_WOE
	printf("Using A(X) div X in order to consider Y as starting point\n");
	total_work--;
#endif

	printf("Reading scalar data in polynomial ``a''\n");
	mnpoly_alloc(cx->a,total_work);
	mnpoly_zero(cx->a,total_work);
	
	cx->cur_deg = read_data_for_series(cx->a,0,total_work,0);

	it->t = cx->t0 = CDIV(m_param,n_param);

	if (cx->cur_deg < it->t) {
		printf( "This amount of data is insufficient.\n"
				"Synchronious revision up to stamp %d is "
				"mandatory\n",it->t);
		exit(1);
	}
			
	/* Data read stage completed. */

	nbpoly_alloc(it->f,total_work);
	nbpoly_zero(it->f,total_work);
	
	it->clist = my_malloc(bigdim * sizeof(struct col_id));

	if (start==-1) {
		gmp_randstate_t		randstate;
		unsigned long int	seed;

		mpz_init(blah);

		gmp_randinit(randstate, GMP_RAND_ALG_DEFAULT, 128);
		/*
		seed=time(NULL);
		*/
		seed=0x3b7b99de;
		gmp_randseed_ui(randstate,seed);
		printf("Using random generator with seed=0x%08lx\n",seed);
		for(k=0; k < it->t  ;k++)
		for(j=0; j < m_param;j++)
		for(i=0; i < n_param;i++) {
			mpz_urandomb(blah,
					randstate,
					bw_allocsize * mp_bits_per_limb);
			mpz_fdiv_r(blah,blah,modulus);
			MPN_SET_MPZ(nbmat_scal(nbpoly_coeff(it->f,k),i,j),
					bw_allocsize,blah);
		}
		gmp_randclear(randstate);
		mpz_clear(blah);

		for(i=0;i<n_param;i++) {
			bw_scalar_set_one(nbmat_scal(nbpoly_coeff(it->f,it->t),
						i,m_param+i));
		}
		for(j=0;j<bigdim;j++) {
			it->clist[j].degnom = it->t;
			it->clist[j].pos    = j;
		}
	} else {
		printf("Trying to recover from iteration # %d\n",start);
		bw_reload_f(it,start);
		it->t=start;
	}


	it->chance_list = my_malloc(bigdim * sizeof(int));
	it->pivots_list = my_malloc(m_param * sizeof(int));
	memset(it->chance_list,0,bigdim * sizeof(int));
	
	mbmat_alloc(it->ctaf);
	
	/* Here, the basic setup is done. Unless we are very unlucky, and
	 * ctaf is singular. We will check against that later, since this
	 * is not very likely.
	 *
	 * In the case of GF(2), the odds are bigger of course. Choice of
	 * another ``f'' polynomial is likely to trim the probabilities
	 * down to a minimum. If this still doesn't work, increasing t0
	 * might help. Otherwise, we're in trouble. */

#ifdef ACTIVATE_HACKS
	printf("Doing some precomputations on the modulus\n");
	do_modulus_precomps();
#endif
}
#endif

static void
bw_init(struct bw_iterator * it, struct bw_context * cx, int start)
{
	int		i,j,k;
	bw_mbmat	reduced_rank;
	int		* ik;
	int		* rk;
	int		* pivots;
	int		r;
	bw_mcol		acol;
	bw_scalar	h_coeff;
	mpz_t h;


#ifdef	CORRECT_STUPID_WOE
	printf("Using A(X) div X in order to consider Y as starting point\n");
	total_work--;
#endif

	printf("Reading scalar data in polynomial ``a''\n");
	mnpoly_alloc(cx->a,total_work);
	mnpoly_zero(cx->a,total_work);
	
	cx->cur_deg = read_data_for_series(cx->a,0,total_work,0);

	/* Data read stage completed. */

	nbpoly_alloc(it->f,total_work);
	nbpoly_zero(it->f,total_work);

	it->clist = my_malloc(bigdim * sizeof(struct col_id));

	if (start==-1) {
		printf("Computing t0\n");
		mbmat_alloc(reduced_rank);
		mcol_alloc(acol);
		ik = malloc(m_param * sizeof(int));
		rk = malloc(m_param * sizeof(int));
		pivots = malloc(m_param * sizeof(int));
		bw_scalar_alloc(h_coeff,BW_SCALAR_SHORT);
		mpz_init(h);
		r = 0;

		for(k=0;r < m_param && k<cx->cur_deg;k++) {
			for(j=0;r < m_param && j<n_param;j++) {
				int u,v;
				for(i=0;i<m_param;i++) {
					memcpy(mcol_scal(acol,i),mnmat_scal(mnpoly_coeff(cx->a,k),i,j),bw_allocsize*sizeof(mp_limb_t));
				}
				rk[r]=j;
				/* kill as many coeffs as we can */
				for(v=0;v<r;v++) {
					bw_mcol pcol;
					u=pivots[v];
					/* use pcol[u] to cancel acol[u] */
					memcpy(h_coeff,mcol_scal(acol,u),bw_allocsize*sizeof(mp_limb_t));

					pcol=mbmat_col(reduced_rank,v);
					for(u=0;u<m_param;u++) {
						addmul(mcol_scal(acol,u),mcol_scal(pcol,u),h_coeff);
						bw_reduce_short_scalar(mcol_scal(acol,u));
					}
				}
				for(u=0;u<m_param;u++) {
					if (!bw_scalar_is_zero(mcol_scal(acol,u)))
						break;
				}

				if (u==m_param) {
					printf("[X^%d] col %d -> (rank == %d)\n"
							,k,j,r);
				} else {
					mpz_t w;
					mpz_init(w);
					pivots[r]=u;
					MPZ_SET_MPN(h,mcol_scal(acol,u),bw_allocsize);
					mpz_neg(h,h);
					mpz_invert(h,h,modulus);

					for(v=0;v<m_param;v++) {
						MPZ_SET_MPN(w,mcol_scal(acol,v),bw_allocsize);
						mpz_mul(w,w,h);
						mpz_mod(w,w,modulus);
						MPN_SET_MPZ(mcol_scal(acol,v),bw_allocsize,w);
					}
					mcol_copy(mbmat_col(reduced_rank,r),acol);
					printf("col. %d, [X^%d] -> rank++ (head row %d)\n",j,k,u);
					ik[r]=k;
					r++;
					mpz_clear(w);
				}
			}
		}
			
		if (r!=m_param) {
			printf("This amount of data is insufficient.\n");
			exit(1);
		}
		it->t = cx->t0 = ik[r-1]+1;
				
		/* Now build f */

		/* First n columns: identity matrix */
		for(i=0;i<n_param;i++) 
			bw_scalar_set_one(nbmat_scal(nbpoly_coeff(it->f,0),i,i));

		/* rest: X^(s-ik)rk's */
		for(j=0;j<m_param;j++) {
			bw_scalar_set_one(nbmat_scal(nbpoly_coeff(it->f,it->t-ik[j]),rk[j],n_param+j));
		}

		for(j=0;j<bigdim;j++) {
			it->clist[j].degnom = it->t;
			it->clist[j].pos    = j;
		}
		free(ik);
		free(rk);
		free(pivots);
		free(h_coeff);
		mbmat_free(reduced_rank);
		mcol_free(acol);
		mpz_clear(h);
	} else {
		printf("Trying to recover from iteration # %d\n",start);
		bw_reload_f(it,start);
		it->t=start;
	}


	it->chance_list = my_malloc(bigdim * sizeof(int));
	it->pivots_list = my_malloc(m_param * sizeof(int));
	memset(it->chance_list,0,bigdim * sizeof(int));
	
	mbmat_alloc(it->ctaf);
	
#ifdef ACTIVATE_HACKS
	printf("Doing some precomputations on the modulus\n");
	do_modulus_precomps();
#endif
}

/*
 * Normalize a column.
 * Give -1/mat(i0,j) (modulo ad hoc.)
 */
static void
col_normalizer(mpz_t dst, bw_mbmat mat, int i0, int j)
{
	mpz_t pz;

	MPZ_INIT_SET_MPN(pz,mbmat_scal(mat,i0,j),bw_allocsize);

	if (mpz_invert(dst,pz,modulus) == 0) {
		mpz_gcd(dst,pz,modulus);
		mpz_out_str(stdout,10,dst);
		printf(" is a factor of the modulus\n");
		exit(191);
	}
	mpz_clear(pz);
	mpz_neg(dst,dst);
	mpz_fdiv_r(dst,dst,modulus);
}

/* Cancel a column.
 * The reference column is jpiv. The column j is added the appropriate
 * multiple of column jpiv so that the coefficient at (i0,j) cancels out.
 *
 * ctaf version :
 *
 * On input, (i0,jpiv) must be the first nonzero coefficient on the jpiv
 * column. Thus, the coefficients for i<i0 in column j are left
 * untouched. The value in norm is -1/mat(i0,jpiv).
 * hnorm is set to -mat(i0,j)/mat(i0,jpiv)
 *
 * f version : 
 *
 * No such thing is required. All computations are carried out. norm is
 * not needed. hnorm is used.
 *
 */
static void
ctaf_col_cancel(bw_mbmat mat, int i0, int j, int jpiv,
		mpz_t norm, mpz_t hnorm)
{
	mpz_t cz,ez;
	int i;

	MPZ_SET_MPN(hnorm,mbmat_scal(mat,i0,j),bw_allocsize);
	mpz_init(cz);
	mpz_init(ez);
			
	mpz_mul(hnorm,hnorm,norm);
	mpz_fdiv_r(hnorm,hnorm,modulus);
	bw_scalar_set_zero(mbmat_scal(mat,i0,j));
	
	for(i=i0+1;i<m_param;i++) {
		MPZ_SET_MPN(cz,mbmat_scal(mat,i,j),bw_allocsize);
		MPZ_SET_MPN(ez,mbmat_scal(mat,i,jpiv),bw_allocsize);
		mpz_mul(ez,hnorm,ez);
		mpz_add(cz,cz,ez);
		mpz_fdiv_r(cz,cz,modulus);
		MPN_SET_MPZ(mbmat_scal(mat,i,j),bw_allocsize,cz);
	}
	
	mpz_clear(cz);
	mpz_clear(ez);
}

static void
f_col_cancel(bw_nbmat mat, int j, int jpiv, mpz_t hnorm)
{
	mpz_t cz,ez;
	int i;

	mpz_init(cz);
	mpz_init(ez);
			
	for(i=0;i<n_param;i++) {
		MPZ_SET_MPN(cz,nbmat_scal(mat,i,j),bw_allocsize);
		MPZ_SET_MPN(ez,nbmat_scal(mat,i,jpiv),bw_allocsize);
		mpz_mul(ez,hnorm,ez);
		mpz_add(cz,cz,ez);
		mpz_fdiv_r(cz,cz,modulus);
		MPN_SET_MPZ(nbmat_scal(mat,i,j),bw_allocsize,cz);
	}
	
	mpz_clear(cz);
	mpz_clear(ez);
}
	

/*
 * We much prefer a human-readable ordering.
 */
static int col_cmp(const struct col_id * x, const struct col_id * y)
{
	int diff;
	diff = x->degnom - y->degnom;
	return diff?diff:(x->pos - y ->pos);
}

static void
bw_gauss(struct bw_iterator * it)
{
	int i,j,jr,k,t;
	bw_scalar piv;
	int rank;
	mpz_t norm,hnorm;

	/* Sort the columns in ascending order of nominal degre. */
	qsort(it->clist,bigdim,sizeof(struct col_id),(sortfunc_t)&col_cmp);

	mpz_init(norm);
	mpz_init(hnorm);

	/* Pay attention here, this is a gaussian elimination on
	 * *columns* */
	rank = 0 ;
	for(j = 0 ; j < bigdim ; j++) {
		jr = it->clist[j].pos;
		/* Find the pivot inside the column. */
		for(i = 0 ; i < m_param ; i++) {
			piv = mbmat_scal(it->ctaf,i,jr);
			/*bw_reduce_short_scalar(piv); */
			/* It's done in check_chance */
			if (!bw_scalar_is_zero(piv))
				break;
		}
		if (i == m_param)
			continue;
		assert(rank<m_param);
		it->pivots_list[rank++] = j;
		col_normalizer(norm,it->ctaf,i,jr);
		/* Cancel this coeff in all other columns. */
		for(k = j + 1 ; k < bigdim ; k++) {
			ctaf_col_cancel(it->ctaf,i,it->clist[k].pos,jr,
					norm,hnorm);
			for(t = 0 ; t <= it->clist[k].degnom ; t++) {
				f_col_cancel(nbpoly_coeff(it->f,t),
						it->clist[k].pos,jr,
						hnorm);
			}
		}

	}

	mpz_clear(norm);
	mpz_clear(hnorm);
	
	if (rank!=m_param) {
		fprintf(stderr,"Duh, rank is not == m !\n");
		exit(1);
	}
}

/*
 * Multiply the specified column by ``X''
 */
static void
col_multiply_x(bw_nbpoly f, int j, int pdeg)
{
	int k;

	for ( k = pdeg ; k >= 0 ; k-- ) {
		ncol_copy(	nbmat_col(nbpoly_coeff(f,k+1),j),
				nbmat_col(nbpoly_coeff(f,k),j));
	}
	ncol_zero(nbmat_col(nbpoly_coeff(f,0),j));
}


/*
 * Multiply my the appropriate ``D'' matrix.
 *
 * This operation takes place on the columns, of course.
 */
static void
bw_iterate(struct bw_iterator * it)
{
	int j;

	for(j = 0 ; j < m_param ; j++) {
		col_multiply_x(it->f,
				it->clist[it->pivots_list[j]].pos,
				it->clist[it->pivots_list[j]].degnom++);
	}

	it->t++;
}

/*
 * Check if some columns, by chance, are already 0 in ctaf. This is an
 * indicator that we found a combination of vectors which is orthogonal
 * to the space spanned by the t(x)M^i (and thus, we hope, is =0).
 *
 * Two consecutive ``chances'' on the same column are considered as a
 * great incentive to try out this column.
 */
static int
bw_check_chance(struct bw_iterator * it)
{
	int i,j;
	int maxchance;

	maxchance=0;

	for(j=0;j<bigdim;j++) {
		for(i=0;i<m_param;i++) {
			bw_reduce_short_scalar(mbmat_scal(it->ctaf,i,j));
		}
	}

	for(j=0;j<bigdim;j++) {
		if (mcol_is_zero(mbmat_col(it->ctaf,j))) {
			if (++it->chance_list[j] > maxchance)
				maxchance=it->chance_list[j];
			printf("Column %d happens to be zero ! (%d)\n",
					j,it->chance_list[j]);
		} else
			it->chance_list[j]=0;
	}

	return maxchance;
}


/*
 * Free the internal structures before leaving.
 */
static void
bw_free(struct bw_iterator * it, struct bw_context * cx)
{
	mnpoly_free(cx->a);
	nbpoly_free(it->f);
	mbmat_free(it->ctaf);
	free(it->clist);
	free(it->pivots_list);
	free(it->chance_list);
	free(it);
	free(cx);
}

void showuse(void)
{
	die("Usage : bw-master <bank#>\n",1);
}

#define MASTER_NUM_OPTIONS 1

void
block_wiedemann(int start)
{
	int mavg;
	double totaltime;

	struct timeval start_time;

	it = my_malloc(sizeof(struct bw_iterator));
	cx = my_malloc(sizeof(struct bw_context));
	
	bw_init(it,cx,start);
	timer_r(&start_time, TIMER_SET);
	totaltime=0.0;


	for(	mavg=bigdim*cx->t0+(it->t-cx->t0)*m_param;
#if 0
		it->t * bigdim - mavg <= bigdim * CDIV(ncols,m_param); 
#else
		it->t < total_work;
#endif
		mavg += m_param)
	{
		double asympt;
		double tt;
		struct timeval step_time;

		if (it->t >= cx->cur_deg) {
			printf("Now fetching additional data\n");
			/* First read as much as we can */
			cx->cur_deg = read_data_for_series(cx->a,
					cx->cur_deg,
					total_work,
					0);
		}
		if (it->t >= cx->cur_deg) {
			/* If not enough, wait to have at least a decent
			 * amount of data available */
			cx->cur_deg = read_data_for_series(cx->a,
					cx->cur_deg,
					MIN(total_work,it->t+2*periodicity),
					1);
		}

#ifdef	VERY_VERBOSE
		bw_dump_f(stdout,it);
#endif
		timer_r(&step_time, TIMER_SET);
		bw_do_ctaf();
		if (bw_check_chance(it) > 2)
			break;
		bw_gauss(it);
		bw_iterate(it);
		tt=timer_r(&step_time, TIMER_ASK);
		totaltime=timer_r(&start_time, TIMER_ASK);

		asympt=(double)it->t;
		asympt=asympt*(asympt+1.0);
		asympt/=2.0;
		asympt*=(double)(m_param*n_param);
		asympt*=(double)(m_param*n_param);
		asympt/=(double)bigdim;

		/* given c, the predicted running time is:
		   c * (mn)^2 / (m+n) * \sum t
		   c * (mn)^2 / (m+n) * tmax^2/2
		 */
		printf("t=%d, average diff=%.2f ; took %.2fs ; total %.2fs ; c=%.4e\n",
				it->t,
				((double)mavg)/((double)(bigdim)),
				tt,totaltime,
				totaltime/asympt);
	}
	bw_commit_f(it);
	stop_threads();
	bw_free(it,cx);
}


int
main(int argc, char *argv[])
{
	struct opt_desc * opts = NULL;
	int n_opts=0;
	int start=-1;

	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	/* puts(version_string); */
	check_endianness(stdout);


	coredump_limit(1);

	new_option(&opts,&n_opts,
		OPT_INTPRM(NULL,&bank_num,OPT_INDEP|OPT_REQUIRED));
	new_option(&opts,&n_opts,
		OPT_INTPRM("--recover -r",&start,OPT_INDEP));
	new_option(&opts,&n_opts,
		OPT_INTPRM("--nthreads -n",&want_threads,OPT_ONCE));

	process_options(argc,argv,n_opts,opts);

	set_all_filenames(bank_num);

	read_tag_file();

#ifdef HARDCODE_PARAMS
	consistency_check("m_param",computed_m_param,m_param);
	consistency_check("n_param",computed_n_param,n_param);
	consistency_check("bw_allocsize",computed_bw_allocsize,bw_allocsize);
	consistency_check("bw_longsize",computed_bw_longsize,bw_longsize);
#endif

	/********************************************************************/

	block_wiedemann(start);

	return 0;
}
