#include "cado.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>



#include "portability.h"
#include "utils.h"
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "facul.h"
#include "finding_good_strategy.h"
#include "generate_factoring_method.h"
#include "generate_strategies.h"
#include "tab_strategy.h"
#include "tab_fm.h"
#include "tab_decomp.h"

int CONST_TEST_R = 55;

/*
This binary allows to test our procedure choosing optimal
strategies. In fact, using the strategy of CADO, we can use it to
approximate the "theorical" number of relations found per second found
with certains parameters. Then, comparing the "theorical" value with the
real value, we could detect if a problem exists in our procedure.
 */

tabular_fm_t *generate_methods_cado(const int lpb)
{
    int n = nb_curves(lpb);
    tabular_fm_t *res = tabular_fm_create();
    fm_t *fm = fm_create();
    unsigned long method[4];

    /* run one P-1 curve with B1=315 and B2=2205 */
    method[0] = PM1_METHOD;	//method
    method[1] = 0;		//curve
    method[2] = 315;		//B1
    method[3] = 2205;		//B2
    fm_set_method(fm, method, 4);
    tabular_fm_add_fm(res, fm);

    /* run one P+1 curve with B1=525 and B2=3255 */
    method[0] = PP1_27_METHOD;	//method
    method[1] = 0;		//curve
    method[2] = 525;		//B1
    method[3] = 3255;		//B2
    fm_set_method(fm, method, 4);
    tabular_fm_add_fm(res, fm);

    /* run one ECM curve with Montgomery parametrization, B1=105, B2=3255 */
    method[0] = EC_METHOD;	//method
    method[1] = MONTY12;	//curve
    method[2] = 105;		//B1
    method[3] = 3255;		//B2
    fm_set_method(fm, method, 4);
    tabular_fm_add_fm(res, fm);

    if (n > 0) {
	method[0] = EC_METHOD;	//method
	method[1] = BRENT12;	//curve
	method[2] = 315;	//B1
	method[3] = 5355;	//B2
	fm_set_method(fm, method, 4);
	tabular_fm_add_fm(res, fm);
    }

    /* heuristic strategy where B1 is increased by sqrt(B1) at each curve */
    double B1 = 105.0;
    for (int i = 4; i < n + 3; i++) {
	double B2;
	unsigned int k;
	B1 += sqrt(B1);
	B2 = 17.0 * B1;
	/* we round B2 to (2k+1)*105, thus k is the integer nearest to
	   B2/210-0.5 */
	k = B2 / 210.0;
	method[0] = EC_METHOD;	//method
	method[1] = MONTY12;	//curve
	method[2] = B1;
	method[3] = (2 * k + 1) * 105;	//B2
	fm_set_method(fm, method, 4);
	tabular_fm_add_fm(res, fm);
    }
    fm_free(fm);
    return res;
}

/*
This function generates the strategy of cado and computes the
  probability and the time to find a prime divisor in a cofactor of
  'r' bits with the bound fbb and lpb.  This strategy is the
  concatenation of all methods in 'methods'.
*/
tabular_strategy_t *generate_strategy_cado(tabular_fm_t * methods,
					   tabular_decomp_t * tab_dec, int fbb,
					   int lpb, int r)
{
    tabular_strategy_t *tab_strat = tabular_strategy_create();
    strategy_t *strat = strategy_create();

    int lim = 2 * fbb - 1;
    if (r < lim) {
	fm_t *zero = fm_create();
	unsigned long method[4] = { PM1_METHOD, 0, 0, 0 };
	fm_set_method(zero, method, 4);
	strategy_add_fm(strat, zero);
	strategy_set_time(strat, 0.0);

	if (r != 1 && (r < fbb || r > lpb))
	    strategy_set_proba(strat, 0.0);
	else			// r==1 or fbb<= r0 <= lpb
	    strategy_set_proba(strat, 1.0);
	fm_free(zero);
    } else {
      int len = 3 + nb_curves(lpb);
      //printf ("len  = %d\n", len);
      ASSERT(len <= methods->index);
      /* strategy_add_fm(strat, methods->tab[0]); */
      //strategy_add_fm(strat, methods->tab[2]);
      for (int i = 0; i < len; i++)
	strategy_add_fm(strat, methods->tab[i]);

	//eval
	double p = compute_proba_strategy(tab_dec, strat, fbb, lpb);

	double t = compute_time_strategy(tab_dec, strat, r);
	/* if (r == CONST_TEST_R) */
	/*   { */
	/*     printf ("p =  %lf, t = %lf\n", p, t); */
	/*     getchar(); */
	/*   } */
	strategy_set_proba(strat, p);
	strategy_set_time(strat, t);
    }

    tabular_strategy_add_strategy(tab_strat, strat);
    strategy_free(strat);

    return tab_strat;
}

/*
  generate the matrix with the strategy of CADO.
*/

tabular_strategy_t ***generate_matrix_cado(const char *name_directory_decomp,
					   tabular_fm_t * methods,
					   int fbb0, int lpb0, int mfb0,
					   int fbb1, int lpb1, int mfb1)
{
    /*
       allocates the matrix which contains all optimal strategies for
       each couple (r0, r1): r0 is the lenght of the cofactor in
       the first side, and r1 for the second side.
     */
    tabular_strategy_t ***matrix = malloc(sizeof(*matrix) * (mfb0 + 1));
    ASSERT(matrix != NULL);

    for (int r0 = 0; r0 <= mfb0; r0++) {
	matrix[r0] = malloc(sizeof(*matrix[r0]) * (mfb1 + 1));
	ASSERT(matrix[r0] != NULL);
    }

    /*
       For each r0 in [fbb0..mfb0], we precompute and strore the data
       for all strategies.  Whereas, for each r1, the values will be
       computed sequentially. 
     */
    fm_t *zero = fm_create();
    unsigned long method_zero[4] = { 0, 0, 0, 0 };
    fm_set_method(zero, method_zero, 4);

    tabular_strategy_t **data_rat = malloc(sizeof(*data_rat) * (mfb0 + 1));
    ASSERT(data_rat);

    int lim0 = 2 * fbb0 - 1;
    for (int r0 = 0; r0 <= mfb0; r0++) {
	tabular_decomp_t *tab_decomp = NULL;
	if (r0 >= lim0) {
	    char name_file[200];
	    sprintf(name_file,
		    "%s/decomp_%d_%d", name_directory_decomp, fbb0, r0);
	    FILE *file = fopen(name_file, "r");

	    tab_decomp = tabular_decomp_fscan(file);

	    if (tab_decomp == NULL) {
		fprintf(stderr, "impossible to read '%s'\n", name_file);
		exit(EXIT_FAILURE);
	    }
	    fclose(file);
	}
	data_rat[r0] = generate_strategy_cado(methods, tab_decomp,
					      fbb0, lpb0, r0);
	tabular_decomp_free(tab_decomp);
    }

    /*
       read good elements for r_2 in the array data_r1 and compute the
       data for each r_1. So :
     */
    int lim1 = 2 * fbb1 - 1;
    for (int r1 = 0; r1 <= mfb1; r1++) {
	tabular_decomp_t *tab_decomp = NULL;
	if (r1 >= lim1) {
	    char name_file[200];
	    sprintf(name_file,
		    "%s/decomp_%d_%d", name_directory_decomp, fbb1, r1);
	    FILE *file = fopen(name_file, "r");

	    tab_decomp = tabular_decomp_fscan(file);

	    if (tab_decomp == NULL) {
		fprintf(stderr, "impossible to read '%s'\n", name_file);
		exit(EXIT_FAILURE);
	    }
	    fclose(file);
	}

	tabular_strategy_t *strat_r1 =
	    generate_strategy_cado(methods, tab_decomp, fbb1, lpb1, r1);

	tabular_decomp_free(tab_decomp);

	for (int r0 = 0; r0 <= mfb0; r0++) {
	    tabular_strategy_t *res =
		generate_strategy_r0_r1(data_rat[r0], strat_r1);
	    matrix[r0][r1] = res;
	}
	tabular_strategy_free(strat_r1);
    }

    //free
    for (int r0 = 0; r0 <= mfb0; r0++)
	tabular_strategy_free(data_rat[r0]);
    free(data_rat);

    fm_free(zero);
    return matrix;
}


/************************************************************************/
/*                  To bench our strategies                             */
/************************************************************************/


int
select_random_index_dec(double sum_nb_elem, tabular_decomp_t* t)
{
    //100000 to consider approximation of distribution
    double alea = rand() % 10000;
    int i = 0;
    double bound = (t->tab[0]->nb_elem/sum_nb_elem) * 10000;
    int len = t->index;
    while (i < (len-1) && (alea - bound) >= 1) {
	i++;
	bound += (t->tab[i]->nb_elem/sum_nb_elem) * 10000;
    }
    return i;
}


/*
  This function compute directly the probabity and the time to find a
  factor in a good decomposition in a cofactor of length r.
 */
double*
bench_proba_time_st(gmp_randstate_t state, facul_strategy_t* strategy,
		    tabular_decomp_t* init_tab, int r, MAYBE_UNUSED int lpb)
{
    double* res = malloc (2*sizeof (double));
    int nb_test = 0, nb_success = 0;
    int nb_test_max = 100000;
    double time = 0;
    
  
    double sum_dec = 0;
    for (int i = 0; i < init_tab->index; i++)
	sum_dec += init_tab->tab[i]->nb_elem;

    mpz_t N[nb_test_max];
    for (int i =0; i < nb_test_max; i++)
      {
	mpz_init (N[i]);
	int index = select_random_index_dec(sum_dec, init_tab);
	int len_p = init_tab->tab[index]->tab[1];
	//generation of the integer N (which will be factoring)
	generate_composite_integer(N[i],state, len_p, r);
	//printf ("%d, %d\n", len_p, r);
	//generate_composite_integer(N[i],state, 18,r);
	//generate_prime_factor(N[i], state, 50);
      }
    //}}
    double starttime = microseconds();
    while (nb_test < nb_test_max)
	{
	    /*
	       f will contain the prime factor of N that the strategy
	       found.  Note that N is composed by two prime factors by the
	       previous function.
	    */
	    unsigned long f[2];
	    f[0] = 0;
	    facul(f, N[nb_test], strategy);
	    if (f[0] != 0)
		{
		  nb_success++;
		}
	    nb_test++;
	}
    double endtime = microseconds();
    time = endtime - starttime;
    
    res[0] = nb_success / ((double)nb_test);
    res[1] = time / ((double)nb_test);
    
    for (int i = 0; i < nb_test_max; i++)
      mpz_clear(N[i]);

    return res;
}

/************************************************************************/
/*                            USAGE                                     */
/************************************************************************/

static void declare_usage(param_list pl)
{
    param_list_usage_header(pl,
			    "This binary allows to test the strategy of cado,"
			    "and especially compute the theorical number of "
			    "relations found per second by this strategy.\n");

    param_list_decl_usage(pl, "fbb0",
			  "set rationnal factor base bound to 2^fbb0\n");
    param_list_decl_usage(pl, "fbb1",
			  "set algebraic factor base bound to 2^fbb1\n");
    param_list_decl_usage(pl, "lpb0",
			  "set rational large prime bound to 2^lpb0");
    param_list_decl_usage(pl, "lpb1",
			  "set algebraic large prime bound to 2^lpb1");
    param_list_decl_usage(pl, "mfb0", "set the first cofactor bound to 2^mfb0");
    param_list_decl_usage(pl, "mfb1",
			  "set the second cofactor bound to 2^mfb1");
    param_list_decl_usage(pl, "decomp",
			  "to locate the file or the directory , according to\n"
			  "\t \t if you need one or several files,\n"
			  "\t \t which contain(s) the file(s) of cofactors decompositions.");
    param_list_decl_usage(pl, "dist",
			  "the pathname of our file which contains the distribution\n"
			  "\t\t of our pairs of cofactors.");
    param_list_decl_usage(pl, "t",
			  "specify the time (seconds) to optain cofactors in the file\n"
			  "\t\t given by the option 'dist'.");
    param_list_decl_usage(pl, "out",
			  "the output file which contain our strategies\n");
}


/************************************************************************/
/*     MAIN                                                             */
/************************************************************************/

int main(int argc, char *argv[])
{

    param_list pl;
    param_list_init(pl);
    declare_usage(pl);

    if (argc <= 1) {
	param_list_print_usage(pl, argv[0], stderr);
	exit(EXIT_FAILURE);
    }

    argv++, argc--;
    for (; argc;) {
	if (param_list_update_cmdline(pl, &argc, &argv)) {
	    continue;
	}
	/* Could also be a file */
	FILE *f;
	if ((f = fopen(argv[0], "r")) != NULL) {
	    param_list_read_stream(pl, f);
	    fclose(f);
	    argv++, argc--;
	    continue;
	}
	fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
	param_list_print_usage(pl, argv[0], stderr);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }

    //option parser
    int fbb0 = -1;
    int lpb0 = -1;
    int fbb1 = -1;
    int lpb1 = -1;
    int mfb0 = -1;
    int mfb1 = -1;
    double C0 = -1;

    param_list_parse_int(pl, "lpb0", &lpb0);
    param_list_parse_int(pl, "lpb1", &lpb1);
    param_list_parse_int(pl, "fbb1", &fbb1);
    param_list_parse_int(pl, "mfb1", &mfb1);
    param_list_parse_int(pl, "fbb0", &fbb0);
    param_list_parse_int(pl, "mfb0", &mfb0);
    param_list_parse_double(pl, "t", &C0);

    if (fbb0 == -1 || lpb0 == -1 || mfb0 == -1 ||
	fbb0 == -1 || lpb0 == -1 || mfb0 == -1 || C0 == -1) {
	fputs("ALL parameters are mandatory!\n", stderr);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }

    //int fbb = (fbb0 < fbb1) ? fbb0 : fbb1;
    //    int lpb = (lpb0 > lpb1) ? lpb0 : lpb1;
    //convert the time in micro-s. because all previous binaries
    //compute their times in micro-s.
    C0 *= 1000000;		//s-->micro-s 

    //option: tab decomp
    const char *name_directory_decomp;
    //  "/localdisk/trichard/results/decomp_cofactor/decomp_tmp";
    if ((name_directory_decomp =
	 param_list_lookup_string(pl, "decomp")) == NULL) {
	fputs("Parser error: Please re-run with the option "
	      "-decomp and a valid directory name.\n", stderr);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    //option: distribution cofactors 
    const char *name_file_cofactor;
    //"/localdisk/trichard/cado768/cofactors";
    if ((name_file_cofactor = param_list_lookup_string(pl, "dist")) == NULL) {
    	fputs("Parser error: Please re-run with the option -dist"
    	      "followed by the pathname of the file which stores the "
    	      "distribution of our cofactors.\n", stderr);
    	param_list_clear(pl);
    	exit(EXIT_FAILURE);
    }

    FILE *file_C = fopen(name_file_cofactor, "r");
    unsigned long **distrib_C = extract_matrix_C(file_C, mfb0 + 1, mfb1 + 1);
    if (distrib_C == NULL) {
    	fprintf(stderr, "Error while reading file %s\n", name_file_cofactor);
    	param_list_clear(pl);
    	exit(EXIT_FAILURE);
    }
    fclose(file_C);

    gmp_randstate_t state;
    gmp_randinit_default(state);

    //select our methods
    //tabular_fm_t *methods = generate_methods_cado(lpb);
    //benchmark
    /* bench_proba(state, methods, fbb); */
    /* bench_time(state, methods); */
    /* tabular_fm_print (methods); */
    FILE* file_in = fopen ("/localdisk/trichard/cadoRSA155/data_fm_25", "r");
    tabular_fm_t *methods = tabular_fm_fscan (file_in);
    printf ("len  = %d, (%d)\n", methods->index, 3+nb_curves(lpb0));
    //methods->index = 3+nb_curves(lpb0);
    //change it!!!pb for free and the nb is diff for each side
    
    //generate our strategies
    //remark: for each pair (r0, r1), we have only one strategy!!
    tabular_strategy_t ***matrix =
	generate_matrix_cado(name_directory_decomp, methods,
			     fbb0, lpb0, mfb0,
			     fbb1, lpb1, mfb1);

    //eval our strategy!
    double Y = 0, T = C0;
    for (int r0 = 0; r0 <= mfb0; r0++) {
    	for (int r1 = 0; r1 <= mfb1; r1++) {
    	    Y += distrib_C[r0][r1] * matrix[r0][r1]->tab[0]->proba;
    	    T += distrib_C[r0][r1] * matrix[r0][r1]->tab[0]->time;
    	}
    }
    //print the result!
    printf(" Y = %lf relations, T = %lf s., yt = %1.10lf rel/s\n", Y,
    	   T / 1000000, Y / T * 1000000);

    const char *pathname_output;
    pathname_output = param_list_lookup_string(pl, "out");

    if (pathname_output != NULL) {
	FILE *file_output = fopen(pathname_output, "w");
	for (int r0 = 0; r0 <= mfb0; r0++)
	    for (int r1 = 0; r1 <= mfb1; r1++)
		if (matrix[r0][r1]->tab[0] != NULL) {
		    fprintf(file_output,
			    "[r0=%d, r1=%d] : (p = %lf, t = %lf)\n",
			    r0, r1, matrix[r0][r1]->tab[0]->proba,
			    matrix[r0][r1]->tab[0]->time);
		    strategy_fprint_design(file_output, matrix[r0][r1]->tab[0]);
		}
	fclose(file_output);
    }

    //{{test bench strategy!
    /* int r0 = CONST_TEST_R; */
    /* int r1 = CONST_TEST_R; */
    /* printf ("r0 = %d, r1 = %d, lpb0 = %d, lpb1 = %d\n", r0, r1, lpb0, lpb1); */
    /* //printf ("nb call = %lu\n", distrib_C[r0][r1]); */
    /* printf ("proba found: %lf, %lf\n", matrix[r0][r1]->tab[0]->proba, matrix[r0][r1]->tab[0]->time); */
    /* char name_file[200]; */
    /* sprintf(name_file, */
    /* 	    "%s/decomp_%d_%d", name_directory_decomp, fbb0, r0); */
    /* FILE *file = fopen(name_file, "r"); */
    
    /* tabular_decomp_t* tab_decomp = tabular_decomp_fscan(file); */
    
    /* tabular_decomp_print (tab_decomp); */
    /* fclose (file); */
    
    /* facul_strategy_t* facul_st = facul_make_strategy (fbb0, lpb0, 0, 0); */

    /* double* res = bench_proba_time_st(state, facul_st, tab_decomp, r0, lpb0); */
    /* double p0 = res[0], t0 = res[1]; */
    /* printf ("side = 0, proba = %lf, time = %lf\n", res[0], res[1]); */
    /* tabular_decomp_free (tab_decomp); */
    /* facul_clear_strategy(facul_st); */
    /* //r1 */
    /* sprintf(name_file, */
    /* 	    "%s/decomp_%d_%d", name_directory_decomp, fbb1, r1); */
    /* file = fopen(name_file, "r"); */
    
    /* tab_decomp = tabular_decomp_fscan(file); */
    /* fclose (file); */
    
    /* facul_st = facul_make_strategy (fbb1, lpb1, 0, 0); */
    /* res = bench_proba_time_st(state, facul_st, tab_decomp, r1, lpb1); */
    /* printf ("side = 1, proba = %lf, time = %lf\n", res[0], res[1]); */
    /* printf ("two sides, proba = %lf, time = %lf\n", res[0]*p0, t0+p0*res[1]); */
    /* tabular_decomp_free (tab_decomp); */
    /* //Free facil_st */
    /* facul_clear_strategy(facul_st); */
    //}}
    
    //free
    for (int r0 = 0; r0 <= mfb0; r0++) {
	for (int r1 = 0; r1 <= mfb1; r1++)
	    tabular_strategy_free(matrix[r0][r1]);
	free(distrib_C[r0]);
	free(matrix[r0]);
    }
    free(distrib_C);
    free(matrix);
    gmp_randclear(state);
    param_list_clear(pl);
    return EXIT_SUCCESS;
}
