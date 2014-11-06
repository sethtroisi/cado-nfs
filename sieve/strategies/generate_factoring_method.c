#include <float.h>
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

#include "cado.h"
#include "portability.h"
#include "utils.h"
#include "facul.h"
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "getprime.h"

#include "generate_factoring_method.h"
#include "convex_hull.h"

/*
  BOUND_SIGMA is used when you generate a random value of sigma.
*/
const int BOUND_SIGMA = 100;

/* 
   The test fails with probabily equals to 1/4^15 (approx 0) (MILLER RABIN).
*/
const int NB_PRIMALITY_TEST = 15;

const int BOUND_FILTER = 20;	//todo: unused!!! 

const double EPSILON_DBL = LDBL_EPSILON;

/*
  If the probability to find a prime number of p bits is less than 3%, 
  we can suppose that the probability to find p+i bits (i>0) is equal to 0.
*/
const double BENCH_MIN_PROBA = 0.03;

/************************************************************************/
/*                      To generate numbers n=pq                        */
/************************************************************************/

/*
  this function allows to compute an array with the probability of 
  scoring of each size of prime number in [min, ..., max].
 */
double *distribution_prime_number(int min, int max)
{
    int len = max - min +1;
    ASSERT(len >= 0);

    double *tab = malloc(len * sizeof(double));
    ASSERT(tab != NULL);

    unsigned long elem = pow(2, min);
    int log_elem = min;
    for (int i = 0; i < len; i++) {
	tab[i] = (elem - 1) / log_elem;
	elem *= 2;
	log_elem++;
    }

    double proba = 1;
    double sum = 0;
    for (int i = 0; i < len; i++) {
	tab[i] *= proba;
	proba /= 2;
	sum += tab[i];
    }

    ASSERT(sum != 0);
    for (int i = 0; i < len; i++)
	tab[i] /= sum;

    return tab;
}

/*
  This function generates in 'res' a random prime number of 'lenFact' bits.
*/

void generate_prime_factor(mpz_t res, gmp_randstate_t state, int lenFact)
{
    mpz_t min;
    mpz_t max;
    mpz_t diff;
    mpz_init(min);
    mpz_init(max);
    mpz_init(diff);
    mpz_ui_pow_ui(min, 2, lenFact - 1);
    do {
	mpz_urandomm(res, state, min);
	mpz_add(res, res, min);

    }
    while (mpz_probab_prime_p(res, NB_PRIMALITY_TEST) == 0);

    //clear
    mpz_clear(min);
}

/*
  This function generates in 'res' a random composite number of 'lenFactall' 
  bits with a prime factor of 'lenFact1' bits.
*/
void
generate_composite_integer(mpz_t res, gmp_randstate_t state,
			   int lenFact1, int lenFactall)
{
    mpz_t bound_max, bound_min;
    mpz_init(bound_max);
    mpz_init(bound_min);
    mpz_ui_pow_ui(bound_min, 2, lenFactall - 1);	//2^(lenfactall)
    mpz_mul_ui(bound_max, bound_min, 2);
    mpz_sub_ui(bound_max, bound_max, 1);	//2^(lenfactall+1)-1

    mpz_t p;
    mpz_init(p);

    mpz_t q;
    mpz_init(q);
    do {
	/*p */
	generate_prime_factor(p, state, lenFact1);
	/*q */
	generate_prime_factor(q, state, lenFactall - lenFact1);

	/* compute n */
	mpz_mul(res, q, p);

    } while (mpz_cmp(bound_min, res) > 0 || mpz_cmp(bound_max, res) < 0);

    // clear
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(bound_max);
    mpz_clear(bound_min);
}


int
select_random_index_according_dist(double *dist, int len)
{
    int precision = 100000;
    //100000 to consider approximation of distribution
    int alea = rand() % precision;
    int i = 0;
    int bound = (int)(dist[0] * precision);
    while (i < (len-1) && alea > bound) {
	i++;
	bound += (int)(dist[i] * precision);
    }
    return i;
}

/*
  This function generates in 'res' a random composite number of 'lenFactall' 
  bits with a prime factor where this size depends to the distrubution 
  of prime numbers (in 'dist')!
*/
int
generate_composite_integer_interval(mpz_t res, gmp_randstate_t state,
				    double *dist, int lenFact1_min,
				    int lenFact1_max, int lenFactall)
{
    int len = lenFact1_max -lenFact1_min +1;
    int index = select_random_index_according_dist(dist, len);
    generate_composite_integer(res, state, lenFact1_min + index, lenFactall);
    return lenFact1_min + index;
}

/************************************************************************/
/*                To model our factoring methods                        */
/************************************************************************/

/* 
   This function allows to create a type facul_strategy_t which contains all
   informations for our factoring method.
   This is necessary, when you would use the function facul () to realize 
   our bench.
*/

facul_strategy_t *generate_fm(int method, int curve, unsigned long B1,
			      unsigned long B2)
{
    facul_strategy_t *strategy;
    strategy = malloc(sizeof(facul_strategy_t));
    strategy->methods = malloc(2 * sizeof(facul_method_t));
    /*
       without this second method, the function
       facul_clear_strategy (strategy) failed!
     */
    strategy->methods[1].method = 0;

    strategy->lpb = ULONG_MAX;
    strategy->assume_prime_thresh = 0;
    strategy->BBB = 0;
    strategy->methods[0].method = method;

    if (method == PM1_METHOD) {
	strategy->methods[0].plan = malloc(sizeof(pm1_plan_t));
	ASSERT(strategy->methods[0].plan != NULL);
	pm1_make_plan(strategy->methods[0].plan, B1, B2, 0);
    } else if (method == PP1_27_METHOD || method == PP1_65_METHOD) {
	strategy->methods[0].plan = malloc(sizeof(pp1_plan_t));
	ASSERT(strategy->methods[0].plan != NULL);
	pp1_make_plan(strategy->methods[0].plan, B1, B2, 0);
    } else if (method == EC_METHOD) {
	long sigma;
	if (curve == 3) {
	    curve = rand() % 3;
	}
	if (curve == MONTY16) {
	    sigma = 1;
	} else
	    sigma = 2 + rand();

	strategy->methods[0].plan = malloc(sizeof(ecm_plan_t));
	ASSERT(strategy->methods[0].plan != NULL);
	ecm_make_plan(strategy->methods[0].plan, B1, B2, curve,
		      labs(sigma), 0, 0);
    } else {
	exit(EXIT_FAILURE);
    }
    return strategy;
}

/************************************************************************/
/*            BENCH THE PROBABILITY AND TIME OF OUR METHODS             */
/************************************************************************/

/*
  This function collects the probability of success to find a prime number 
  (of 'len_p' bits) in a composite number (of 'len_n' bits) 
  with the strategy.
*/

static double
bench_proba_fm(facul_strategy_t * strategy,
	       gmp_randstate_t state, unsigned long len_p, unsigned long len_n)
{
    int nb_succes_max = 1000, nb_succes = 0;
    int nb_test = 0;
    int nb_test_max = 10 * nb_succes_max;
    mpz_t N;
    mpz_init(N);

    while (nb_succes < nb_succes_max && (nb_test < nb_test_max)) {
	/* generation of the integer N (which will be factoring) */
	generate_composite_integer(N, state, len_p, len_n);

	if (strategy->methods[0].method == EC_METHOD) {
	    ecm_plan_t *plan = strategy->methods[0].plan;
	    //RANDOM_SIGMA
	    if (plan->parameterization != MONTY16) {
		plan->sigma = 2 + rand() % BOUND_SIGMA;
		strategy->methods[0].plan = plan;
	    }
	}
	/* 
	   f will contain the prime factor of N that the strategy
	   found.  Note that N is composed by two prime factors by the
	   previous function.
	 */
	unsigned long f[2];
	f[0] = 0;

	facul(f, N, strategy);

	if (f[0] != 0)
	    nb_succes++;
	nb_test++;
    }

    //clear
    mpz_clear(N);

    return nb_succes / ((double)nb_test);
}

static double
sub_routine_bench_time(facul_strategy_t * method, gmp_randstate_t state,
		       int len_n)
{
    double tps = 0;
    int nb_test = 10000;
    mpz_t N;
    mpz_init(N);
    unsigned long f[2];

    for (int i = 0; i < nb_test; i++) {
	generate_prime_factor(N, state, len_n);
	//compute the the time of execution
	uint64_t starttime, endtime;
	starttime = microseconds();

	facul(f, N, method);

	endtime = microseconds();

	tps += endtime - starttime;
    }
    //clear
    mpz_clear(N);

    return tps / (nb_test);

}

/*
  This function collects the times to run the strategy st.
  These time are computed for integers of bits size MODREDCUL_MAXBITS, 
  MODREDC15UL_MAXBITS, MODREDC2UL2_MAXBITS and 3*MODREDCUL_MAXBITS.
 */
static double *bench_time_fm(facul_strategy_t * st, gmp_randstate_t state)
{
    double *res = calloc(4, sizeof(double));
    ASSERT(res != NULL);
    res[0] = sub_routine_bench_time(st, state, MODREDCUL_MAXBITS - 1);
    res[1] = sub_routine_bench_time(st, state, MODREDC15UL_MAXBITS - 1);
    res[2] = sub_routine_bench_time(st, state, MODREDC2UL2_MAXBITS - 1);
    res[3] = sub_routine_bench_time(st, state, MODREDC15UL_MAXBITS * 2);
    return res;
}

/*
  This function allows to collect for each factoring methods, 
  the probability to find a prime number of a certain size, 
  from 'len_p_min' until the probability become null. 
*/
void bench_proba(gmp_randstate_t state, tabular_fm_t * fm, int len_p_min)
{
    int len = fm->index;	//number of methods!
    int p_max = 100;
    double *proba = malloc(p_max * sizeof(double));
    ASSERT(proba != NULL);

    unsigned long *param;
    fm_t *elem;
    for (int i = 0; i < len; i++) {
	elem = tabular_fm_get_fm(fm, i);
	param = fm_get_method(elem);
	unsigned long method = param[0];
	unsigned long curve = param[1];
	unsigned long B1 = param[2];
	unsigned long B2 = param[3];

	facul_strategy_t *st = generate_fm(method, curve, B1, B2);

	int ind_proba = 0;
	do {
	    if (B1 == 0 && B2 == 0)
		proba[ind_proba] = 0;
	    else {
		int len_p = len_p_min + ind_proba;
		/*
		   TODO: check if it's sufficient to don't be disturbed by
		   the second factor of n.
		 */
		int len_n = 60 + len_p;
		proba[ind_proba] = bench_proba_fm(st, state, len_p, len_n);
	    }
	    ind_proba++;
	} while ((proba[ind_proba - 1] - BENCH_MIN_PROBA) > EPSILON_DBL
		 && ind_proba < p_max);

	fm_set_proba(elem, proba, ind_proba - 1, len_p_min);
	//free
	facul_clear_strategy(st);
    }

    free(proba);
}

/*
  This function allows to collect for each factoring methods, 
  the differents time to find a prime number in a interger
  of differentbits size.
*/
void bench_time(gmp_randstate_t state, tabular_fm_t * fm)
{
    unsigned long *param;
    fm_t *elem;
    int len = fm->index;	//number of methods!

    for (int i = 0; i < len; i++) {
	elem = tabular_fm_get_fm(fm, i);
	param = fm_get_method(elem);
	unsigned long method = param[0];
	unsigned long curve = param[1];
	unsigned long B1 = param[2];
	unsigned long B2 = param[3];
	if (B1 != 0 || B2 != 0) {
	    facul_strategy_t *st = generate_fm(method, curve, B1, B2);
	    double *time = bench_time_fm(st, state);
	    fm_set_time(elem, time, 4);
	    //free
	    free(time);
	    facul_clear_strategy(st);
	} else {
	    double time[4] = { 0, 0, 0, 0 };
	    fm_set_time(elem, time, 4);
	}
    }
}

/************************************************************************/
/*       BENCH according to an interval of size of prime numbers        */
/************************************************************************/

/*
  This function allows to compute the probability and the time of a strategy
  to find a prime number in an interval [2**len_p_min, 2**len_p_max].
*/

static double *sub_routine_bench_proba_cost_interval(facul_strategy_t *
						     strategy,
						     gmp_randstate_t state,
						     int len_p_min,
						     int len_p_max, int len_n)
{
    double *res = malloc(2 * sizeof(double));
    ASSERT(res != NULL);
    int nb_succes_lim = 1000, nb_succes = 0;
    int nb_test = 0;
    int nb_test_max = 10 * nb_succes_lim;

    mpz_t N;
    mpz_init(N);
    double tps = 0;
    double *disp = distribution_prime_number(len_p_min, len_p_max);

    while ((nb_succes < nb_succes_lim) && (nb_test < nb_test_max)) {
	/* generation of the integer N (which will be factoring) */
	generate_composite_integer_interval(N, state, disp,
					    len_p_min, len_p_max, len_n);

	if (strategy->methods[0].method == EC_METHOD) {
	    ecm_plan_t *plan = strategy->methods[0].plan;
	    //RANDOM_SIGMA
	    if (plan->parameterization != MONTY16) {
		plan->sigma = 2 + rand() % BOUND_SIGMA;
		strategy->methods[0].plan = plan;
	    }
	}
	//compute the the time of execution
	double starttime, endtime;
	starttime = microseconds();
	/* 
	   f will contain the prime factor of N that the strategy
	   found.  Note that N is composed by two prime factors by the
	   previous function.
	 */
	unsigned long f[2];
	f[0] = 0;

	facul(f, N, strategy);

	if (f[0] != 0)
	    nb_succes++;

	endtime = microseconds();
	tps += endtime - starttime;

	nb_test++;
    }

    //clear
    mpz_clear(N);
    free(disp);

    res[0] = nb_succes / ((double)nb_test);
    res[1] = tps / ((double)nb_test);

    return res;
}

/* 
   This function returns default parameters for the collection of 
   factoring methods and stores their in 'res' such that: 
   res[0]=b1_min 
   res[1]=b1_max
   res[2]=b1_step 
   res[3]=c_min 
   res[4]=c_max 
   res[5]=c_step 
   TODO: these parameters could be certainly increased :)
*/
double *choice_parameters(int method, int len_p_min)
{
    double b1_min, b1_max, b1_step, c_min, c_max, c_step;
    if (method == PM1_METHOD ||
	method == PP1_27_METHOD || method == PP1_65_METHOD) {
	if (len_p_min < 20) {
	    b1_min = 20;
	    b1_max = 1000;
	    b1_step = 1.2;
	    c_min = 10;
	    c_max = 50;
	    c_step = 5;
	} else if (len_p_min <= 25) {
	    b1_min = 20;
	    b1_max = 4500;
	    b1_step = 1.4;
	    c_min = 40;
	    c_max = 80;
	    c_step = 5;
	} else if (len_p_min < 30) {
	    b1_min = 100;
	    b1_max = 10000;
	    b1_step = 1.4;
	    c_min = 20;
	    c_max = 200;
	    c_step = 10;
	} else {
	    b1_min = 100;
	    b1_max = 30000;
	    b1_step = 1.4;
	    c_min = 20;
	    c_max = 500;
	    c_step = 10;
	}

    } else if (method == EC_METHOD) {
	if (len_p_min < 20) {
	    b1_min = 20;
	    b1_max = 5000;
	    b1_step = 1.2;
	    c_min = 10;
	    c_max = 100;
	    c_step = 10;
	} else if (len_p_min <= 25) {
	    b1_min = 100;
	    b1_max = 10000;
	    b1_step = 1.4;
	    c_min = 10;
	    c_max = 200;
	    c_step = 40;
	} else if (len_p_min < 30) {
	    b1_min = 100;
	    b1_max = 13000;
	    b1_step = 1.4;
	    c_min = 10;
	    c_max = 500;
	    c_step = 20;
	} else {
	    b1_min = 200;
	    b1_max = 30000;
	    b1_step = 1.4;
	    c_min = 10;
	    c_max = 1000;
	    c_step = 20;
	}

    } else {
	return NULL;
    }
    double *result = malloc(6 * sizeof(double));
    ASSERT(result != NULL);
    result[0] = b1_min;
    result[1] = b1_max;
    result[2] = b1_step;
    result[3] = c_min;
    result[4] = c_max;
    result[5] = c_step;
    return result;
}

/*
  This function returns an array which contains, for each method, 
  the probability and the cost to find a prime number for each size between 
  'len_p_min' and 'len_p_max'. The type of these methods is specified by 
  two parameters : 'method' and 'curve', and the parameters B1 and B2 are 
  given by 'param_region' or, if it's equal to NULL, use the default sieve.
*/

static tabular_fm_t *bench_proba_cost_interval(int method, int curve,
					       gmp_randstate_t state,
					       int len_p_min, int len_p_max,
					       int len_n, double *param_region)
{
    //define the sieve region
    int c_min, c_max;
    int b1_min, b1_max;
    double c_pas, b1_pas;
    if (param_region != NULL) {
	b1_min = param_region[0];
	b1_max = param_region[1];
	b1_pas = param_region[2];

	c_min = param_region[3];
	c_max = param_region[4];
	c_pas = param_region[5];
    } else {
	//default parameters for the sieve region.
	double *param = choice_parameters(method, len_p_min);
	b1_min = param[0];
	b1_max = param[1];
	b1_pas = param[2];

	c_min = param[3];
	c_max = param[4];
	c_pas = param[5];
	free(param);
    }

    double *result;
    double max_proba = 0.9;

    tabular_fm_t *tab_fusion = tabular_fm_create();

    //add zero method
    unsigned long tmp_method[4] = { method, 0, 0, 0 };
    double zero = 0;
    tabular_fm_add(tab_fusion, tmp_method, 4, &zero, 1, &zero, 1, len_p_min);

    for (int c = c_min; c <= c_max; c += c_pas) {
	int B1;
	int B2;

	B1 = b1_min;
	B2 = B1 * c;

	tabular_fm_t *tab = tabular_fm_create();
	unsigned long elem[4];
	double proba = 0;
	double tps = 0;
	while (B1 < b1_max && proba < max_proba) {
	    facul_strategy_t *fm = generate_fm(method, curve, B1, B2);

	    result = sub_routine_bench_proba_cost_interval
		(fm, state, len_p_min, len_p_max, len_n);

	    proba = result[0];
	    tps = result[1];

	    elem[0] = method;
	    elem[1] = curve;
	    elem[2] = B1;
	    elem[3] = B2;

	    tabular_fm_add(tab, elem, 4, &proba, 1, &tps, 1, len_p_min);

	    facul_clear_strategy(fm);

	    free(result);

	    if (b1_pas <= 2)	//that is a percent.
		B1 = (B1 * b1_pas);
	    else
		B1 = B1 + b1_pas;
	    B2 = B1 * c;

	}
	//merge arrays
	tabular_fm_concat(tab_fusion, tab);
	tabular_fm_free(tab);
    }
    return tab_fusion;
}

/************************************************************************/
/*                     GENERATE FACTORING METHODS                       */
/************************************************************************/

/*
  This function returns for each sub-interval of prime number, the probability
  and the time to find a prime number in N (an interger of len_n bits).
  With this function, we can specify the method, the curve and the sieve 
  region of parameters for B1 and B2. Note that an option 'opt_ch' 
  allows to apply the selection by convex hull.
*/
tabular_fm_t *generate_factoring_methods_mc(gmp_randstate_t state,
					    int len_p_min, int len_p_max,
					    int len_n, int method, int curve,
					    int opt_ch, double *param_sieve)
{
    ASSERT(len_p_min <= len_p_max);

    tabular_fm_t *gfm;

    tabular_fm_t *collect = bench_proba_cost_interval
	(method, curve, state, len_p_min, len_p_max, len_n, param_sieve);

    if (opt_ch) {
	//apply the convex hull
	gfm = convex_hull_fm(collect);

	tabular_fm_free(collect);
    } else {
	double *param;
	if (param_sieve == NULL)
	    param = choice_parameters(method, len_p_min);
	else
	    param = param_sieve;

	gfm = collect;
	/*
	   Warning: the function choice_parameters() allocate dynamic
	   memory. So, don't forget to free it!
	 */
	if (param_sieve == NULL)
	    free(param);
    }

    return gfm;
}

/*
  This function calls the previous function (generate_factoring_methods_mc ())
  for each one of our functions.
*/

tabular_fm_t *generate_factoring_methods(gmp_randstate_t state, int len_p_min,
					 int len_p_max, int len_n, int opt_ch,
					 double *param_sieve)
{

    tabular_fm_t *gfm = tabular_fm_create();

    //we begin by the first method:
    int method = PM1_METHOD;
    int curve = 0;

    while (method <= NB_METHOD && curve < NB_CURVE) {
	printf("method = %d, curve = %d\n", method, curve);
	tabular_fm_t *res = generate_factoring_methods_mc
	    (state, len_p_min, len_p_max, len_n, method, curve, opt_ch,
	     param_sieve);

	tabular_fm_concat(gfm, res);

	//free
	tabular_fm_free(res);

	//index
	if (method != EC_METHOD)
	    method++;
	else
	    curve++;
    }

    return gfm;
}

/************************************************************************/
/*                    MERGE and compute the convex hull                 */
/************************************************************************/
/*
  This function collects factoring methods from the file 'file_in' 
  and make a selection by convex hull and prints them in file_out.
 */
tabular_fm_t *convex_hull_from_file(FILE * file_in, FILE * file_out)
{
    tabular_fm_t *all_st = tabular_fm_fscan(file_in);
    if (all_st == NULL)
	return NULL;

    tabular_fm_t *res = convex_hull_fm(all_st);

    tabular_fm_free(all_st);

    int err = tabular_fm_fprint(file_out, res);
    if (err < 0)
	return NULL;

    return res;
}

/************************************************************************/
/*                      FILTERING                                       */
/************************************************************************/

//this function filter methods in fm!
//nb_prime_nb is the max of len_proba!
// TODO: modify his content!!! unused in this moment
tabular_fm_t *filtering(tabular_fm_t * fm, int len_p_min, int nb_prime_nb)
{
    int len = fm->index;
    //create compromis
    double **compromis = malloc((len) * sizeof(double *));
    ASSERT(compromis != NULL);
    for (int i = 0; i < len; i++) {
	compromis[i] = malloc((nb_prime_nb) * sizeof(double));
	ASSERT(compromis[i] != NULL);

	fm_t *el = tabular_fm_get_fm(fm, i);
	if (fm_is_zero(el))
	    for (int j = 0; j < nb_prime_nb; j++)
		compromis[i][j] = 0;
	else {
	    for (int j = 0; j < nb_prime_nb; j++) {
		if (el->proba[j] == 0 || j > el->len_proba) {
		    compromis[i][j] = INFINITY;
		} else {
		    compromis[i][j] = el->time[j] / el->proba[j];
		}
	    }
	}
    }

    double **dist = malloc(len * sizeof(double *));
    ASSERT(dist != NULL);
    for (int i = 0; i < len; i++) {
	dist[i] = malloc((len) * sizeof(double));
	ASSERT(dist[i] != NULL);
    }
    //compute the matrix
    for (int i = 0; i < len; i++) {
	dist[i][i] = 0;
	for (int j = i + 1; j < len; j++) {
	    double moy_dist = 0;
	    int nb_el = 0;
	    for (int p = 0; p < nb_prime_nb; p++) {
		double tmp = 0;
		if (compromis[i][p] != INFINITY && compromis[j][p] != INFINITY) {
		    tmp = compromis[i][p] - compromis[j][p];	//etablir
								//un
								//pourcentage!!!tmp/compromis[i][p]
		    nb_el++;
		}
		if (tmp < 0)
		    tmp *= -1;
		moy_dist += tmp;
	    }
	    if (nb_el == 0)
		moy_dist = INFINITY;
	    moy_dist /= nb_el;
	    dist[i][j] = moy_dist;
	    dist[j][i] = moy_dist;
	}
    }
    double seuil = BOUND_FILTER;
    double min;
    int curve_near[2];
    curve_near[0] = 0;
    curve_near[1] = 0;
    double *disp =
	distribution_prime_number(len_p_min, len_p_min + nb_prime_nb);

    do {
	min = INFINITY;
	for (int i = 0; i < len; i++) {
	    for (int j = i + 1; j < len; j++) {
		if (dist[i][j] < min) {
		    curve_near[0] = i;
		    curve_near[1] = j;
		    min = dist[i][j];
		}
	    }
	}
	if (min <= seuil) {
	    /* //keep the curve wich have the best probability */
	    printf("les deux courbes proches sont: %d, %d\n valeur de min %f\n",
		   curve_near[0], curve_near[1], min);

	    double moy_proba_0 = 0;
	    double moy_proba_1 = 0;
	    fm_t *el0 = tabular_fm_get_fm(fm, curve_near[0]);
	    fm_t *el1 = tabular_fm_get_fm(fm, curve_near[1]);
	    double *proba0 = fm_get_proba(el0);
	    double *proba1 = fm_get_proba(el1);
	    int nb_min =
		(el0->len_proba <
		 el0->len_proba) ? el0->len_proba : el1->len_proba;
	    for (int p = 0; p < nb_min; p++) {
		moy_proba_0 += proba0[p] * disp[p];
		moy_proba_1 += proba1[p] * disp[p];
	    }
	    printf("moy_proba = %f\t, %f\n", moy_proba_0, moy_proba_1);
	    if (moy_proba_0 < moy_proba_1) {
		//cancel the method curve_near[0]
		tabular_fm_put_zero(fm, curve_near[0]);
		for (int i = 0; i < len; i++) {
		    dist[i][curve_near[0]] = INFINITY;
		    dist[curve_near[0]][i] = INFINITY;
		}
	    } else {
		//cancel the method curve_near[1]
		tabular_fm_put_zero(fm, curve_near[1]);
		for (int i = 0; i < len; i++) {
		    dist[i][curve_near[1]] = INFINITY;
		    dist[curve_near[1]][i] = INFINITY;
		}
	    }
	}
    } while (min <= seuil);

    //filter ECM METHODS:: take advantage that ECM can be repeat over
    //several curves!!!!
    for (int i = 0; i < len; i++) {
	fm_t *elem = tabular_fm_get_fm(fm, i);
	unsigned long *method_i = fm_get_method(elem);
	double *proba_i = fm_get_proba(elem);
	double *temps_i = fm_get_time(elem);
	int len_i = fm_get_len_proba(elem);
	if (method_i[0] == EC_METHOD && method_i[1] != MONTY16) {
	    for (int j = i + 1; j < len; j++) {
		elem = tabular_fm_get_fm(fm, j);
		unsigned long *method_j = fm_get_method(elem);
		double *proba_j = fm_get_proba(elem);
		double *temps_j = fm_get_time(elem);
		int len_j = fm_get_len_proba(elem);
		int nb_min = (len_i < len_j) ? len_i : len_j;
		if (method_j[0] == EC_METHOD) {
		    int cpt_b = 0;
		    //becarefull: the data t isn't sorted!!! so we
		    //need to know which is better than the other.
		    int min = (temps_i[0] < temps_j[0]) ? i : j;
		    int max = (min == i) ? j : i;
		    double comp_min = 0;
		    double comp_max = 0;
		    double p1, p2, t1, t2;
		    for (int p = 0; p < nb_min; p++) {	//gestion des min/max
			if (min == i) {
			    //printf ("1---> %d\n", i);
			    p1 = proba_i[p];
			    t1 = temps_i[p];
			    p2 = proba_j[p];
			    t2 = temps_j[p];
			} else {
			    //printf ("1---> %d\n", j);
			    p1 = proba_j[p];
			    t1 = temps_j[p];
			    p2 = proba_i[p];
			    t2 = temps_i[p];
			}
			int ratio = round(t2 / t1);
			/* printf ("numero des courbes: %d, %d\n", i,j); */
			/* printf ("t2 = %f, t1 = %f\n", t2, t1); */
			/* printf("ratio = %d, %f, %f\n",ratio,
			   pow(1-p1, ratio), 1-p2); */
			/* printf ("p1 = %f, t1= %f, p2 = %f, t2 =
			   %f\n", p1, t1, p2, t2); */
			/* getchar(); */
			if (pow(1 - p1, ratio) < 1 - p2)
			    cpt_b++;
			comp_min += ratio * t1 / (1 - pow(1 - p1, ratio));
			comp_max += t2 / p2;
		    }
		    //delete the method max!
		    //printf("compromis = %f, %f\n", comp_min, comp_max);
		    if (cpt_b > floor(nb_min / 2) &&
			comp_min <= comp_max + seuil) {
			tabular_fm_put_zero(fm, max);
			//printf ("supprimer: %d\n", max);
		    }
		    /* printf ("numero des courbes: %d, %d\n", i,j); */
		    /* printf ("succés: %d/%d\n", cpt_b, nb_prime_nb); */
		    //getchar();
		}
	    }
	}
    }
    // collect the best methods!!!
    //tabular_fm_print (t);

    printf("cleaning up!!!!\n");

    tabular_fm_t *res = tabular_fm_create();
    int nb_zero = 0;

    char res_file[] = "filtering_iii";
    int len_str = strlen(res_file);
    FILE *select_fm;
    fm_t *elem;

    for (int i = 0; i < len; i++) {
	if (!tabular_fm_is_zero(fm, i)) {
	    printf("i = %d\n", i);
	    elem = tabular_fm_get_fm(fm, i);
	    int nb_prime = fm_get_len_proba(elem);
	    tabular_fm_add_fm(res, elem);
	    int c100 = (i) / 100;
	    int c10 = (i - 100 * c100) / 10;
	    int c1 = i - 10 * c10 - 100 * c100;
	    res_file[len_str - 3] = c100 + 48;
	    res_file[len_str - 2] = c10 + 48;
	    res_file[len_str - 1] = c1 + 48;
	    select_fm = fopen(res_file, "w");
	    fprintf(select_fm,
		    "#filtering fm[%d]: [%d, ...]\n # len_p\t proba\t tps \t compromis\n",
		    i, len_p_min);

	    for (int k = 0; k < nb_prime; k++)
		fprintf(select_fm, "%d \t %f\t %f\t %f\n", len_p_min + k,
			elem->proba[k], elem->time[k], compromis[i][k]);
	    /* fclose (select_fm); */
	} else if (nb_zero == 0) {
	    tabular_fm_add_fm(res, tabular_fm_get_fm(fm, i));
	    nb_zero++;
	}
    }

    //free
    free(disp);
    for (int i = 0; i < len; i++) {
	free(dist[i]);
	free(compromis[i]);
    }
    free(dist);
    free(compromis);
    return res;
}

/************************************************************************/
/*                      CONVEX_HULL_FM                                  */
/************************************************************************/

/*
  These differents functions allow to use the module convex_hull, 
  to compute the convex hull of a set of factoring methods.
 */

tabular_point_t *convert_tab_point_to_tab_fm(tabular_fm_t * t)
{
    tabular_point_t *res = tabular_point_create();
    int len = tabular_fm_get_index(t);
    fm_t *elem;
    for (int i = 0; i < len; i++) {
	elem = tabular_fm_get_fm(t, i);
	tabular_point_add(res, i, elem->proba[0], elem->time[0]);
    }
    return res;
}

tabular_fm_t *convert_tab_fm_to_tab_point(tabular_point_t * t,
					  tabular_fm_t * init)
{
    tabular_fm_t *res = tabular_fm_create();
    int len = tabular_point_get_index(t);
    for (int i = 0; i < len; i++) {
	int index = point_get_number(tabular_point_get_point(t, i));
	tabular_fm_add_fm(res, init->tab[index]);
    }
    return res;
}

tabular_fm_t *convex_hull_fm(tabular_fm_t * t)
{
    tabular_point_t *tmp = convert_tab_point_to_tab_fm(t);
    tabular_point_t *res = convex_hull(tmp);
    tabular_fm_t *res_fm = convert_tab_fm_to_tab_point(res, t);
    tabular_point_free(tmp);
    tabular_point_free(res);
    return res_fm;
}
