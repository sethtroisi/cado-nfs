#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>

#include "master_params.h"
#include "params.h"
#include <assert.h>
#include "macros.h"
#include "types.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "timer.h"
#include "variables.h"
#include "structure.h"
#include "gmp-hacks.h"
#include "modulus_hacks.h"
/* #include "version.h" */
#include "master_common.h"

#define	xxxVERY_VERBOSE
/* Warning ! ; the (pos) field here is only used from within the gauss
 * stuff, in order to have pivot columns virtually first. It is
 * meaningless to other parts of the program (the columns being
 * interchangeable with each other anyway).
 */

#define ACTIVATE_HACKS
#define OPT_POWEROF2
#define SKIP_KNOWNCOLS
#define BOUND_F 

static struct bw_context	*cx;
static struct bw_iterator	*it;

#ifdef	VERY_VERBOSE
int bw_dump_f(FILE *f, struct bw_iterator * it)
{
	int i,j,k;
	mp_limb_t * s;
	bw_nbpoly p = it->f;
	char * str;
	int * delta;

	delta=malloc(bigdim*sizeof(int));
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
bw_do_ctaf_inner()
{
	int co,j;

	for(j = co = 0 ; j <  bigdim ; j++)	{
		if (ctaf_bounds[j]==-1)
			continue;
		ctaf_column(j,ctaf_bounds[j]);
		co++;
	}
	

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
		ctaf_bounds=malloc(bigdim*sizeof(int));
		memset(ctaf_bounds,0,bigdim*sizeof(int));
	}

	for(rj = 0 ; rj <  bigdim ; rj++)	{
		j=it->clist[rj].pos;
		if (ctaf_bounds[j]==0) {
			ctaf_bounds[j]=it->clist[rj].degnom;
			mcol_zero(mbmat_col(it->ctaf,j));
		}
	}

	bw_do_ctaf_inner();
}

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
	
	cx->cur_deg = read_data_for_series(cx->a);

	/* Data read stage completed. */

	nbpoly_alloc(it->f,total_work);
	nbpoly_zero(it->f,total_work);

	it->clist = malloc(bigdim * sizeof(struct col_id));

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


	it->chance_list = malloc(bigdim * sizeof(int));
	it->pivots_list = malloc(m_param * sizeof(int));
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

double muops(int t)
{
	double a;
	a=(double)t;
	a=a*(a+1.0);
	a/=2.0;
	a*=(double)(m_param*n_param);
	a*=(double)(m_param*n_param);
	a/=(double)bigdim;
	return a;
}

void
block_wiedemann(int start)
{
	int mavg;
	double totaltime;
	double lasttime;

	struct timeval start_time;

	it = malloc(sizeof(struct bw_iterator));
	cx = malloc(sizeof(struct bw_context));
	
	bw_init(it,cx,start);
	timer_r(&start_time, TIMER_SET);
	totaltime=0.0;
	lasttime=0.0;

	double total_muops = muops(total_work);

	for(	mavg=bigdim*cx->t0+(it->t-cx->t0)*m_param;
#if 0
		it->t * bigdim - mavg <= bigdim * CDIV(nrows,m_param); 
#else
		it->t < total_work;
#endif
		mavg += m_param)
	{
		double asympt;
		double tt;
		struct timeval step_time;
		int chance;

#ifdef	VERY_VERBOSE
		bw_dump_f(stdout,it);
#endif
		timer_r(&step_time, TIMER_SET);
		bw_do_ctaf();
		chance = bw_check_chance(it);
		if (chance > 2)
			break;
		bw_gauss(it);
		bw_iterate(it);
		tt=timer_r(&step_time, TIMER_ASK);
		totaltime=timer_r(&start_time, TIMER_ASK);

		asympt = muops(it->t);

		/* given c, the predicted running time is:
		   c * (mn)^2 / (m+n) * \sum t
		   c * (mn)^2 / (m+n) * tmax^2/2
		 */
		if (totaltime > lasttime + 1.0 || chance || !(it->t % 100)) {
			printf("t=%d dbar=%.2f step=%.2fs total=%.2fs c=%.4e %d%% of %ds\n",
					it->t,
					((double)mavg)/((double)(bigdim)),
					tt,totaltime,
					totaltime/asympt,
					(int) (100.0 * asympt / total_muops),
			      (int) (totaltime * total_muops / asympt));
			lasttime = totaltime;
		}
	}

	{
		int j;
		int * delta = malloc(bigdim * sizeof(int));
		for (j = 0; j < bigdim; j++) {
			delta[it->clist[j].pos] =it-> clist[j].degnom;
		}
		bw_commit_f(it->f, delta);
		free(delta);
	}


	print_chance_list(it->t, it->chance_list);
	bw_free(it,cx);
}

void block_signals()
{
	struct sigaction act;
	act.sa_handler=SIG_IGN;
	sigemptyset(&act.sa_mask);
	act.sa_flags=0;
	sigaction(SIGUSR1, &act, NULL);
	sigaction(SIGUSR2, &act, NULL);
	sigaction(SIGHUP,  &act, NULL);
	sigaction(SIGQUIT, &act, NULL);
	sigaction(SIGINT,  &act, NULL);
	sigaction(SIGTERM, &act, NULL);
}

int
main(int argc, char *argv[])
{
	int start=-1;

	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	// coredump_limit(1);
	
	if (argc > 4 && strcmp(argv[1], "--subdir") == 0) {
		int rc = chdir(argv[2]);
		if (rc < 0) {
			perror(argv[2]);
			exit(errno);
		}
		argv+=2;
		argc-=2;
	}

	if (argc != 4) {
		die("Usage: mater-old <mat-file> <m> <n>\n", 1);
	}

	read_mat_file_header(argv[1]);
	m_param = atoi(argv[2]);
	n_param = atoi(argv[3]);
#ifndef HARDCODE_PARAMS
        bigdim = m_param + n_param;
#endif

	total_work = Lmacro(nrows, m_param, n_param);

#ifdef HARDCODE_PARAMS
	consistency_check("m_param",computed_m_param,m_param);
	consistency_check("n_param",computed_n_param,n_param);
	consistency_check("bw_allocsize",computed_bw_allocsize,bw_allocsize);
	consistency_check("bw_longsize",computed_bw_longsize,bw_longsize);
#endif

	block_signals();
	block_wiedemann(start);

	return 0;
}
