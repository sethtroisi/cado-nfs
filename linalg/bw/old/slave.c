#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>
#ifdef ENABLE_PTHREADS
#include <pthread.h>
#endif

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
#include "corners.h"
#include "slave.h"
#include "filenames.h"
#include "matrix.h"
#include "tagfile.h"
#include "timer.h"
#include "threaded.h"
#include "barrier.h"
#include "state.h"
#include "addmul.h"
/* #include "version.h" */

/* Values from recovery.c */
extern int keep_periodicity;
extern barrier_t vs_checkpoint_barrier;



#define TIMER_MTH	TIMER_MTN(number_of_threads)

/* This variable has a default value, but the corresponding parameter is
 * required anyway, so this default value doesn't have much meaning */

enum task_type_t task_type=TASK_DOTPRODUCTS;

const struct textflag_desc task_flag[]={
		{"mksol polynomial poly p m",TASK_POLYNOMIAL},
		{"slave dotproduct dot d s",TASK_DOTPRODUCTS},
		{NULL,0}
};

bw_vector_block bw_w, bw_v;
coord_t *bw_x;		/* These are just elements of the canonical basis */
stype32 * bw_m=NULL;
stype32 * bw_x0=NULL;
bw_scalar bw_s;

mp_limb_t * scalar_block;
bw_vector_block bw_sum;

/* The 1 at the ends puts differential mode to be the default. */
struct mksol_info_block	mksd={-1,-1,-1,NULL,1};

int computed_nbys;
int computed_nbxs;

int recovered_iterations;

corner_s first,last;

void showuse(void)
{
	die("Usage : bw-slave <bank#> <i0>,<j0> <i1>,<j1>\n",1);
}

struct ttime {
	double v;
	char unit;
};

void get_ttime(struct ttime * p, time_t s)
{
	if (s<60) {
		p->v=s;
		p->unit='s';
	} else if (s<3600) {
		p->v=(double)s/60.0;
		p->unit='m';
	} else if (s<86400) {
		p->v=(double)s/3600.0;
		p->unit='h';
	} else if (s<604800) {
		p->v=(double)s/86400.0;
		p->unit='d';
	} else { 
		p->v=(double)s/604800.0;
		p->unit='w';
	}
}

int display_decision(int n, int p)
{
	int per[]={
		1,2,5,
		10,20,50,
		100,200,500,
		1000,2000,5000,
		10000,20000,50000,
		100000,200000,500000,
		1000000,2000000,5000000
	};
	int i;

	for(i=0;i<sizeof(per)/sizeof(per[0]);i++) {
		if (n>p+per[i])
			continue;
		if (n-p==per[i]-(p%per[i]))
			return 1;
	}
	return 0;
}

void ptime(char *s, struct timeval * tv, int offset)
{
	time_t t;
	struct tm toto;

	t=tv->tv_sec+offset;
	localtime_r(&t,&toto);
	strftime(s,80,"%b %d %H:%M", &toto);
	/* sprintf(s+strlen(s),".%04d",tv->tv_usec/10000); */
}



void display_stats(int n, int sync)
{
	int i;
	char s0[80];
	char s1[80];
	char s2[80];
	double avg;
	time_t twe;
	struct ttime ttwe;
	double wall;	/* wall-clock seconds elapsed */
	double cputime;
	double pcpu;

	if (n<=tsv()->go_mark) {
		printf("I can't give you any estimates right now. Wait\n");
		return;
	}

	if (sync) {
		if (!display_decision(n,tsv()->go_mark)&&!tsv()->reset_notice)
			return;
		switch(n) {
			case 1:
			case 100:
				tsv()->reset_notice=1;
		}
	}


	cputime	= timer_r(	&(tsv()->elapsed_time),	TIMER_ASK | TIMER_MTH);
	wall	= wall_clock_r(	&(tsv()->start_time),	TIMER_ASK | TIMER_MTH);

	/*
	printf("T%d CPU time : %.4f\n",tseqid(),cputime);
	printf("T%d wall clock time : %.4f\n",tseqid(),wall);
	*/

	/* cputime is the time spent *in the given thread/process* since
	 * the timer has been reset. This is obtained by a getrusage
	 * difference count from the value stored in tsv()->elapsed_time.
	 */

#define CLOCK_THRESHOLD	0.0001

	if (wall==0) {
		wall=CLOCK_THRESHOLD;
	}

	if (cputime==0) {
		cputime=wall;
	}

	pcpu = 100.0 * cputime / wall;
	avg=cputime/(n-tsv()->go_mark);

	if (!sync) {
		double avg_low;
		avg_low=cputime / (n+1-tsv()->go_mark);
		sprintf(s0,"%.2fs<av<%.2fs",avg_low,avg);
	} else {
		sprintf(s0,"av=%.3fs M0=%.2e",avg,avg/tsv()->ncoeffs);
	}

	twe = (time_t) (total_work * avg);
	get_ttime(&ttwe,twe);

	ptime(s1,&(tsv()->start_time),wall+(total_work-n)*avg);
	ptime(s2,&(tsv()->start_time),wall+(total_work-n)*avg/(pcpu/100.0));

	for(i=0;i<strlen(s1) && i<strlen(s2) && s1[i]==s2[i] ; i++);

	switch(strlen(s2+i)) {
		case 0:	break;
		case 1: i--;
		case 2: i--;
		case 3: break;
		case 4: i--;
		case 5: break;
		case 6: /* should not happen */ break;
		case 7: i--;
		case 8: break;
		default: i=0;
	}
	if (i<0) i=0;	/* should not happen */

	/* That's silly, but we have to be kind with the reentrant stdio */
	if (strlen(s1+i)) {
		printf("T%d N=%d %s %d%% tw=%.2f%c eta=<%s>..<%s>\n",
				tseqid(),
				n,
				s0,
				(int)floor(pcpu),
				ttwe.v, ttwe.unit,
				s1,s2+i);
	} else {
		printf("T%d N=%d %s %d%% tw=%.2f%c eta=<%s>\n",
				tseqid(),
				n,
				s0,
				(int)floor(pcpu),
				ttwe.v, ttwe.unit,
				s1);
	}

	if (sync && tsv()->reset_notice) {
		timer_r(	&(tsv()->elapsed_time),	TIMER_SET|TIMER_MTH);
		wall_clock_r(	&(tsv()->start_time),	TIMER_SET|TIMER_MTH);
		tsv()->go_mark=tsv()->iter;
		ptime(s2,&(tsv()->start_time),0);
		printf("T%d <%s> : synchronous timer reset.\n",
				tseqid(),s2);
		tsv()->reset_notice=0;
	}
}
			
void sighandler(int signum)
{
	if (!is_a_working_thread())
		return;

	display_stats(tsv()->iter,0);
	fflush(stdout);
	fflush(stderr);
	if (signum==SIGUSR1 || signum==SIGUSR2) {
		tsv()->reset_notice=1;
		printf("A synchronous timer reset has been scheduled\n");
	}
}

#ifndef __SVR4
/* the solaris headers are lame, and this code crashes gcc */
void fpe_sighandler(int signum, siginfo_t * info, void * addr)
{
	fflush(stdout);
	fflush(stderr);
	fprintf(stderr,"Caught SIGFPE !!\n");
#ifdef __GNUC__
	{
		void * p;
		p=__builtin_return_address(0);
		if (p==NULL) goto endframe;
		fprintf(stderr,"address %p\n",p);
		p=__builtin_return_address(1);
		if (p==NULL) goto endframe;
		fprintf(stderr,"address %p\n",p);
		p=__builtin_return_address(2);
		if (p==NULL) goto endframe;
		fprintf(stderr,"address %p\n",p);
		p=__builtin_return_address(3);
		if (p==NULL) goto endframe;
		fprintf(stderr,"address %p\n",p);
		p=__builtin_return_address(4);
		if (p==NULL) goto endframe;
		fprintf(stderr,"address %p\n",p);
		p=__builtin_return_address(5);
		if (p==NULL) goto endframe;
		fprintf(stderr,"address %p\n",p);
		p=__builtin_return_address(6);
		if (p==NULL) goto endframe;
		fprintf(stderr,"address %p\n",p);
	}
endframe:
#endif
	fprintf(stderr,"signo:	%d\n",info->si_signo);
	fprintf(stderr,"errno:	%d\n",info->si_errno);
	fprintf(stderr,"code:	%d\n",info->si_code);
	fprintf(stderr,"pid:	%d\n",info->si_pid);
	fprintf(stderr,"uid:	%d\n",info->si_uid);
	fprintf(stderr,"status:	%d\n",info->si_status);
	fprintf(stderr,"utime:	%ld\n",info->si_utime);
	fprintf(stderr,"stime:	%ld\n",info->si_stime);
	/* fprintf(stderr,"value:	%d\n",info->si_value); */
	fprintf(stderr,"int:	%d\n",info->si_int);
	fprintf(stderr,"ptr:	%p\n",info->si_ptr);
	fprintf(stderr,"addr:	%p\n",info->si_addr);
	fprintf(stderr,"band:	%ld\n",info->si_band);
	fprintf(stderr,"fd:	%d\n",info->si_fd);

	if (signum!=SIGFPE) {
		fprintf(stderr, "(SIGFPE in disguise)\n");
	}
}


void setup_sighandler(void)
{
	struct sigaction act;
	act.sa_handler=&sighandler;
	sigemptyset(&act.sa_mask);
	act.sa_flags=0;
	sigaction(SIGUSR1, &act, NULL);
	sigaction(SIGUSR2, &act, NULL);
	sigaction(SIGHUP,  &act, NULL);
	sigaction(SIGQUIT,  &act, NULL);
#ifdef NDEBUG
	sigaction(SIGINT,  &act, NULL);
#endif
	sigaction(SIGTERM, &act, NULL);
}

void setup_fpe_sighandler(void)
{
	struct sigaction act;
	sigemptyset(&act.sa_mask);
	act.sa_sigaction=&fpe_sighandler;
	act.sa_flags=SA_SIGINFO;
	sigaction(SIGFPE,  &act, NULL);
}
#endif

void raw_partial_dot_product(bw_scalar s, stype32 * x, bw_vector_block v)
{
	size_t i;
	bw_long_scalar_set_zero(s);

	bw_lvblock_step_n00(v,tsv()->i0);
	x+=tsv()->i0;
	for(i=tsv()->i0;i<tsv()->i1;i++,x++) {
		addmultiply(s,v,*x);
		/* if (*x) {
			addmultiply(s,v,*x);
		} else {
			addmultiply(s+bw_allocsize+2,v,-*x);
		} */
		bw_lvblock_step_n00(v,1);
	}
}

void pscal(FILE *f, bw_scalar x)
{
	int j;
	fprintf(f,"0x");
	for(j=0;j<bw_allocsize+2;j++)
#if defined(__i386__)
		fprintf(f,"%08lX",x[bw_allocsize+2-1-j]);
#elif defined(__ia64__) || defined(__alpha) || defined(__alpha__)
		fprintf(f,"%016lX",x[bw_allocsize+2-1-j]);
#elif defined(sparc) || defined(__sparc)
		fprintf(f,"%08lX",x[bw_allocsize+2-1-j]);
#elif defined(__x86_64__)
		fprintf(f,"%016lX",x[bw_allocsize+2-1-j]);
#else
#error "Which is your byte size ?"
#endif
}

void dot_products(void)
{
	mp_limb_t carry;
	bw_scalar s,t;
#define n bw_allocsize
#ifdef GROSBUG
	FILE * dbstr=tsv()->dbstr;
#endif

	s=tsv()->s;
	t=tsv()->t;

	raw_partial_dot_product(s,bw_m, tsv()->src);
	raw_partial_dot_product(t,bw_x0,tsv()->dst);
#ifdef GROSBUG
	fprintf(dbstr,"Thread %d iteration %d:\n",tseqid(),tsv()->iter);
	fprintf(dbstr,"S%d: ",tseqid());
	pscal(dbstr,s);
	fprintf(dbstr," - ");
	pscal(dbstr,s+n+2);
	fprintf(dbstr,"\n");
	fprintf(dbstr,"T%d: ",tseqid());
	pscal(dbstr,t);
	fprintf(dbstr," - ");
	pscal(dbstr,t+n+2);
	fprintf(dbstr,"\n");
#endif

	carry=mpn_add_n(s,     s,     t+n+2, n+2);
	assert(carry==0);
	carry=mpn_add_n(s+n+2, s+n+2, t,     n+2);
	assert(carry==0);

	/* s[n+2] is writable, OK for routine. */
	reduce_p_minus_q(s,s,s+n+2);
	/* Now done by reduce_p_minus_q */
	/* bw_reduce_short_scalar(s); */
	memset(s+n,0,(n+4)*sizeof(mp_limb_t));
#undef n
}

void dump_regs(bw_scalar s)
{
	int i;

	fflush(stdout);
	fflush(stderr);

	fprintf(stderr,"Here is the dump of the relevant registers\n");

	fflush(stderr);

	for(i=0;i<number_of_threads;i++) {
		printf("V%d: s=",i);
		pscal(stdout,tsv_table[i].s);
		printf("\n");
	}
	printf("total=");
	pscal(stdout,s);
	printf("\n");
	

	fflush(stdout);
	fflush(stderr);
}

int addup(int n)
{
	bw_scalar s;
	s=bw_s;
	if (n==0) {
		/*memset(s,0,(bw_allocsize+2)*sizeof(mp_limb_t)); */
		memcpy(s,tsv()->s,(bw_allocsize+2)*sizeof(mp_limb_t));
	} else {
		mpn_add_n(s, s, tsv()->s, bw_allocsize+2);
#if 0
		mpn_divrem(s+bw_allocsize,0,
				s,bw_allocsize+2,
				modulus_hbs,bw_allocsize);
#endif
	}
	if (n == number_of_threads - 1) {
		int i;
		bw_shorten_long_scalar(s);
		bw_reduce_short_scalar(s);
		for(i=0;i<bw_allocsize;i++) {
			if (s[i]!=0) {
				fprintf(stderr,
					"s[%d]==0x%lX, not 0 !!!\n",i,s[i]);
				dump_regs(s);
				fprintf(stderr,
					"Error after computation of "
					"iteration %d -> redoing this one.\n",
					tsv()->iter);
				bw_long_scalar_set_zero(s);
				return 0;
			}
		}
		bw_long_scalar_set_zero(s);
	}
	return 1;
}

void check_bw_s(void)
{
}

void show_timers(void)
{
	char s[80];
	ptime(s,&(tsv()->start_time),0);
	printf("Start time set to <%s>\n",s);
	/*
	ptime(s,&(tsv()->elapsed_time),0);
	printf("Elapsed time set to <%s>\n",s);
	*/
}

/*
 * NOTA : m and x0 are dense word-sized vectors. They are stored as a raw
 * binary array.
 */
static void read_m_vector(void)
{
	FILE *f;
	struct stat sbuf;

	f=fopen(mvec_filename,"r");
	if (f==NULL) {
		perror(mvec_filename);
		exit(errno);
	}
	if (fstat(fileno(f),&sbuf)<0) {
		fprintf(stderr, "Cannot stat %s : %s\n",
				mvec_filename,strerror(errno));
		exit(errno);
	}
	bw_m=malloc(sbuf.st_size);
	fread(bw_m,sizeof(stype32),sbuf.st_size/sizeof(stype32),f);
	DO_BIG_ENDIAN(int i;
		for(i=0;i<sbuf.st_size/sizeof(stype32);i++) {
			type32 blah = bw_m[i];
			mswap32(blah);
			bw_m[i]=blah;
		});
	fclose(f);
}

static void read_x0_vector(void)
{
	FILE *f;
	struct stat sbuf;

	f=fopen(x0_filename,"r");
	if (f==NULL) {
		perror(x0_filename);
		exit(errno);
	}
	if (fstat(fileno(f),&sbuf)<0) {
		fprintf(stderr, "Cannot stat %s : %s\n",
				x0_filename,strerror(errno));
		exit(errno);
	}
	bw_x0=malloc(sbuf.st_size);
	fread(bw_x0,sizeof(stype32),sbuf.st_size/sizeof(stype32),f);
	DO_BIG_ENDIAN(int i;
		for(i=0;i<sbuf.st_size/sizeof(stype32);i++) {
			type32 blah = bw_x0[i];
			mswap32(blah);
			bw_x0[i]=blah;
		});
	fclose(f);
}

/* From a *completely stabilized* vector v, do either:
 *
 * - the computation of the nbxs dot products, and append them to the
 * checkpointing engine (WHEN task==TASK_DOTPRODUCTS)
 *
 * - pick a coefficient in the polynomial pool, and add the result to the
 * destination vector (WHEN task==TASK_POLYNOMIAL)
 *
 * in case v is a block of vectors, do the same for each (the result has
 * to be a block of vectors, for parallelization concerns, and also
 * because summing all the components up at the end is the job of
 * bw-gather. Last but not least, it's not bad to be able to resume
 * computing with a different setting than the original one).
 *
 */
void use_vector(bw_vector_block v)
{
	mp_limb_t     * s;
	int		res;
	int		l;

	switch(task_type) {
		case TASK_DOTPRODUCTS:
			commit_more_values(v);
			break;
		case TASK_POLYNOMIAL:
			s=scalar_block;
			res=0;
			for(l=0;l<nbys;l++) {
#ifndef NDEBUG
			if (ftell(tsv()->aux_files[l])!=
				(mksd.valuation+tsv()->iter)*
				bw_filesize*sizeof(mp_limb_t))
			{
				fprintf(stderr,"Bad file position\n");
			}
#endif
				res += bw_scalar_read(s,tsv()->aux_files[l]);
				/*eofs_met += (feof(mksd.poly_files[l]) !=
				 * 0); */
				s+=bw_allocsize;
			}
			if (res != nbys * bw_filesize) {
				fprintf(stderr,
					"Incoherent number of bytes read: %d\n",
					res);
				eternal_sleep();
			}
			bw_multiproduct_separated(bw_sum,v,scalar_block);
			break;
	}
}

void open_polynomials(int n)
{
	int	l;
	char	filename[FILENAME_LENGTH];

	tsv()->aux_files = malloc(nbys * sizeof(FILE*));
	n+=mksd.valuation;
	for(l=0;l<nbys;l++) {
		FILE *f;
		sprintf(filename,f_meta_filename,
				first.j+l,mksd.solution_col,mksd.t_value);
		f = fopen(filename,"r");
		if (f == NULL) {
			perror(filename);
			exit(1);
		}
		fseek(f, n*bw_filesize*sizeof(mp_limb_t),SEEK_SET);
		if (ftell(f)!=n*bw_filesize*sizeof(mp_limb_t)) {
			die("Bad seek in %s\n",1,filename);
		}
		tsv()->aux_files[l]=f;
	}
}

			

/* Warning. Positioning of the barriers is subject to extreme caution
 * here. bw_lvblock_set_zero_separated, multiply, bw_lvblock_reduce all *write* to
 * their first argument, and can work unsynchronized.
 *
 * However, use_vector *reads* from its first argument, in such a
 * manner that a synchronization step is required beforehand.
 *
 * Last but not least, the ++iter needs a mutex, since it is read inside
 * the loop. So the writing and reading are spanned on both sides of the
 * synchronization point.
 */

barrier_t	main_loop_barrier;

void * thread(void * blah)
{
	double * t;
	bw_vector_block src_vec,dst_vec;
	double twait;
	double maxwait=0.0;
	int ok;

	if (task_type==TASK_POLYNOMIAL) {
		if (recovered_iterations) {
			open_polynomials(recovered_iterations+1);
		} else {
			open_polynomials(0);
		}
	}

	bw_scalar_alloc(tsv()->s,BW_SCALAR_LONG);
	bw_scalar_alloc(tsv()->t,BW_SCALAR_LONG);
	bw_long_scalar_set_zero(tsv()->s);
	bw_long_scalar_set_zero(tsv()->t);

	/* NO ! bw_long_scalar_set_zero(bw_s); */

	timer_r(	&(tsv()->elapsed_time),	TIMER_SET | TIMER_MTH);
	wall_clock_r(	&(tsv()->start_time),	TIMER_SET | TIMER_MTH);

	show_timers();

	tsv()->go_mark=tsv()->iter=recovered_iterations;

#ifndef __SVR4	/* see above */
	setup_sighandler();
#endif
	tsv()->ncoeffs=compute_thread_offsets();

	dst_vec=bw_v;

	if (task_type != TASK_POLYNOMIAL || !recovered_iterations)
		use_vector(bw_v);	/* First coeff ! */	/* <v+ */

	for(;tsv()->iter<total_work;)
	{
		/*
		 * Invariant: bw_v = B^iter * y, already commited.
		 *
		 * iter+1 vectors commited
		 */
iteration_1:
		/* set up src and dst */
		tsv()->src=src_vec=bw_v;
		tsv()->dst=dst_vec=bw_w;
		
		/* zero out dst */
		bw_lvblock_set_zero_separated(dst_vec);	/* >wl */

		/* multiply */
		multiply(dst_vec,src_vec);	/* >wl, <v+ */

		/* reduce */
		bw_lvblock_reduce_separated(dst_vec);	/* >wl */

		dot_products();			/* <vl, <wl */

		ok=barrier_wait(&main_loop_barrier,&twait,&addup);
		if (twait>maxwait) {
			printf("Thread %d : waited %.2fs\n",tseqid(),twait);
			maxwait=twait;
		}
		if (!ok)
			goto iteration_1;
			
		tsv()->iter++;

		/* Append to the ``A'' files (buffered) */
		use_vector(dst_vec);	/* <w+ */

		/* Do we need to write the ``V'' file again ? */
		if (tsv()->iter%periodicity==0) state_checkpoint(dst_vec,0);

		/************ first iteration completed ********/
		display_stats(tsv()->iter,1);

		if (tsv()->iter==total_work)
			break;
		/*
		 * Invariant: bw_w = B^iter * y, already commited.
		 *
		 * iter+1 vectors commited
		 */
iteration_2:
		/* set up src and dst */
		tsv()->src=src_vec=bw_w;
		tsv()->dst=dst_vec=bw_v;

		/* zero out dst */
		bw_lvblock_set_zero_separated(dst_vec);	/* >vl */

		/* multiply */
		multiply(dst_vec,src_vec);	/* >vl, <w+ */
		
		/* reduce */
		bw_lvblock_reduce_separated(dst_vec);	/* >vl */
		
		dot_products();			/* <wl, <vl */

		ok=barrier_wait(&main_loop_barrier,&twait,&addup);
		if (twait>maxwait) {
			printf("Thread %d : waited %.2fs\n",tseqid(),twait);
			maxwait=twait;
		}
		if (!ok)
			goto iteration_2;
			
		tsv()->iter++;

		/* Append to the ``A'' files (buffered) */
		use_vector(dst_vec);	/* <v+ */

		/* Do we need to write the ``V'' file again ? */
		if (tsv()->iter%periodicity==0) state_checkpoint(dst_vec,0);
		
		/************ BOTH iterations completed ********/
		display_stats(tsv()->iter,1);
	}

		/*
		 * Invariant: bw_v = B^total_work * y, already commited.
		 *
		 * total_work+1 vectors commited
		 */

	/* Unless we have just commited a checkpoint, issue a final one
	if (tsv()->iter%periodicity!=0)
	 */
	state_checkpoint(dst_vec,1);

	t=malloc(sizeof(double));
	*t=timer_r(&(tsv()->elapsed_time),TIMER_ASK | TIMER_MTH);


	return (void*)t;
}

/* This part is only relevant to the mksol part */
static void mksol_setup(void)
{
	FILE * f;
	char filename[FILENAME_LENGTH];
	int check_col, degnom;

	if (mksd.t_value==-1) {
		die("t value unspecified (option --tv)\n",1);
	}

	if (mksd.solution_col==-1) {
		die("solution col unspecified (option --sc) \n",1);
	}

	bw_lvblock_alloc(bw_sum);
	bw_lvblock_set_zero_separated(bw_sum);

	scalar_block=malloc(nbys*bw_allocsize*sizeof(mp_limb_t));
	memset(scalar_block,0,nbys*bw_allocsize*sizeof(mp_limb_t));

	sprintf(filename,valu_meta_filename,mksd.solution_col,mksd.t_value);
	f = fopen(filename,"r");
	if (f == NULL) {
		perror(filename);
		exit(1);
	}
	if (fscanf(f,"COLUMN %d VALUATION %d DEGNOM %d\n",
				&check_col, &mksd.valuation, &degnom) != 3)
	{
		die("Parsing of valuation file failed.\n",1);
	}
	if (check_col != mksd.solution_col) {
		die("The file %s contains irrelevant data\n",1,filename);
	}
	fclose(f);

	mksd.sum=bw_sum;
	total_work=degnom-mksd.valuation;
	/* Biggest power computed (see invariants in thread()) :
	 *
	 * B^(degnom-valuation)*y;
	 *
	 */
}

int loose_on_limits=0;

int main(int argc, char *argv[])
{
	struct opt_desc * opts = NULL;
	int n_opts = 0;
	int tmp;
	double t;
	int force_work=0;
	int spawned_threads=NTHREADS;

	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	/* puts(version_string); */
	check_endianness(stdout);

	new_option(&opts,&n_opts,
		OPT_INTPRM("--total-work -t",&force_work,OPT_ONCE));
	new_option(&opts,&n_opts,
		OPT_INTPRM(NULL,&bank_num,OPT_INDEP|OPT_REQUIRED));
	new_option(&opts,&n_opts,
		NULL,
		1, PARAM_TEXTUAL(NULL),
		&process_corner,
		&first,
		OPT_INDEP|OPT_REQUIRED);
	new_option(&opts,&n_opts,
		NULL,
		1, PARAM_TEXTUAL(NULL),
		&process_corner,
		&last,
		OPT_INDEP|OPT_REQUIRED);
	new_option(&opts, &n_opts,
		OPT_FLAG("--differential -d",&mksd.differential));
	new_option(&opts, &n_opts,
		OPT_INTPRM("--sc",&mksd.solution_col,OPT_INDEP));
	new_option(&opts, &n_opts,
		OPT_INTPRM("--tv",&mksd.t_value,OPT_INDEP));
	new_option(&opts,&n_opts,
		OPT_FLAG("--loose  -l",&loose_on_limits));
	new_option(&opts,&n_opts,
		OPT_INTPRM("--keep  -k",&keep_periodicity,OPT_ONCE));
	new_option(&opts,&n_opts,
		OPT_INTPRM("--nthreads -n",&spawned_threads,OPT_ONCE));
	new_option(&opts,&n_opts,"--task",
			1,PARAM_TEXTFLAG("slave",task_flag),
			builtin_process_textflag,
			&task_type,
			OPT_ONCE|OPT_REQUIRED,	/* for the sake of clarity */
			"changes the behavior of the program");


	process_options(argc,argv,n_opts,opts);

	coredump_limit(1);
	if (!loose_on_limits) {
		check_limit_requirements();
	} else {
		fprintf(stderr,"Warning, limit checking disabled...\n");
	}


	if (keep_periodicity==0) {
		keep_periodicity=DEFAULT_KEEP;
	}

	if (first.i>last.i) {
		tmp=first.i; first.i=last.i; last.i=tmp;
	}

	if (first.j>last.j) {
		tmp=first.j; first.j=last.j; last.j=tmp;
	}

	computed_nbys=last.j-first.j+1;
	computed_nbxs=last.i-first.i+1;

	/*
	printf("Understood : Bank %d, going from (%d,%d) to (%d,%d)\n",
			bank_num,first.i,first.j,last.i,last.j);
			*/

	set_all_filenames(bank_num);

	read_tag_file();

	if (force_work) {
		printf("Forcing %d iterations instead of %d\n",
				force_work,total_work);
		total_work=force_work;
	}

#ifdef HARDCODE_PARAMS
	consistency_check("nbys",computed_nbys,nbys);
	consistency_check("nbxs",computed_nbxs,nbxs);
	consistency_check("m_param",computed_m_param,m_param);
	consistency_check("n_param",computed_n_param,n_param);
	consistency_check("bw_allocsize",computed_bw_allocsize,bw_allocsize);
	consistency_check("bw_longsize",computed_bw_longsize,bw_longsize);
#endif

	if (first.i<0 || last.i>=m_param || first.j<0 || last.j>=n_param)
		die("The specified area is out of range\n",1);

	load_matrix();

	configure_threads(1,ncols);

	bw_x=malloc(nbxs*sizeof(coord_t));
	bw_lvblock_alloc(bw_v);
	bw_lvblock_alloc(bw_w);
	bw_lvblock_set_zero_separated(bw_w);
	bw_lvblock_set_zero_separated(bw_v);
	bw_scalar_alloc(bw_s,BW_SCALAR_LONG);


	switch(task_type) {
		case TASK_DOTPRODUCTS:
		    	mksd.differential=0;	/* pointless anyway */
			load_x_vectors(bw_x,first.i,nbxs);
			break;
		case TASK_POLYNOMIAL:
			mksol_setup();
			break;
	}

	recovered_iterations=recover_state(bw_v, &mksd);

	read_m_vector();
	read_x0_vector();

	printf("Beginning the job\n");
	printf("%d iterations are needed\n", total_work-recovered_iterations);

	if (recovered_iterations>total_work) {
		printf("There is nothing to do\n");
		return 0;
	}

	fflush(stdout);

#ifdef ENABLE_PTHREADS
	printf("Going multithread : spawning %d threads\n",spawned_threads);
	configure_threads(spawned_threads,ncols);
#else
	configure_threads(1,ncols);
#endif

#ifndef __SVR4	/* see above */
	setup_fpe_sighandler();
#endif

	barrier_init(&main_loop_barrier,number_of_threads,NULL);
	barrier_init(&vs_checkpoint_barrier,number_of_threads,NULL);
	t=start_threads((threadfunc_t) &thread, total_work);
	barrier_destroy(&main_loop_barrier);
	barrier_destroy(&vs_checkpoint_barrier);

	printf("Whole work : avg %.3fs cumulative user time\n", t/total_work);
		
	bw_lvblock_free(bw_w);
	bw_lvblock_free(bw_v);
	free(bw_x);

	return 0;
}

/* vim:set sw=8: */
