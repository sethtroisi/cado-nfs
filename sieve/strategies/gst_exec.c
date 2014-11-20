#include "cado.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>

#include "portability.h"
#include "utils.h"
#include "facul.h"
#include "ecm.h"
#include "tab_strategy.h"
#include "generate_strategies.h"

/************************************************************************/
/*                            USAGE                                     */
/************************************************************************/

static void declare_usage(param_list pl)
{
    param_list_usage_header(pl,
			    "This binary allows to build the best strategies for each couple (r0, r1)\n"
			    "where (r0,r1) are the bits size for our couple of cofactors.\n");

    param_list_decl_usage(pl, "gdc",
    "(switch)  to precompute all decompositions of cofactors in the \n "
"\t \t sieve region. So, you must specify the sieve region with these options:\n"
			  "\t \t -fbb1, -mfb1, -fbb0, -mfb0\n");
    param_list_decl_usage(pl, "gst_r",
			  "(switch)  to precompute the best strategies \n"
			  "\t \t for one bit size cofactor.\n "
			  "\t \t You must specify these options:\n"
			  "\t \t -fbb0, -lpb0, -r -decomp\n");
    param_list_decl_usage(pl, "gst",
	"(switch)  to merge two (or all) precomputing did by the option 'gst_r',\n "
	"\t \t and thus find the best strategie(s) for one (or each) couple (r0,r1).\n"
	"\t \t So, you must specify these options:\n"
	"\t \t -mfb0, -mfb1, -fbb0, -fbb1 ,-in, (-r0 -r1)\n");

    param_list_decl_usage(pl, "r0",
			  "with r1 specify the bit sizes of the cofactor couple"
			  " that you want to find the best strategies.\n");
    param_list_decl_usage(pl, "r1",
			  "with r0 specify the bit sizes of the cofactor couple"
			  " that you want to find the best strategies.\n");
    param_list_decl_usage(pl, "r",
			  "specify the bit size of the studied cofactor.\n");
    param_list_decl_usage(pl, "r",
			  "specify the bit size of the studied cofactor.\n");
    param_list_decl_usage(pl, "fbb1",
			  "set algebraic factor base bound to 2^fbb1\n");
    param_list_decl_usage(pl, "fbb0",
			  "set rationnal factor base bound to 2^fbb0\n");
    param_list_decl_usage(pl, "lpb0",
			  "set rational large prime bound to 2^lpb0");
    param_list_decl_usage(pl, "lpb1",
			  "set algebraic large prime bound to 2^lpb1");
    param_list_decl_usage(pl, "mfb0", "set the first cofactor bound to 2^mfb0");
    param_list_decl_usage(pl, "mfb1", "set the second cofactor bound to 2^mfb1");
    param_list_decl_usage(pl, "in", "to locate the file which contains\n "
	  "\t \t our factoring methods, or locate the directory \n"
	  "\t \t where the precomputed files for option 'gst' are stored");
    param_list_decl_usage(pl, "out", "to locate the directory where the "
			  "file(s) will be stored.");
    param_list_decl_usage(pl, "decomp",
	  "to locate the file or the directory , according to\n"
	  "\t \t if you need one or several files,\n"
	  "\t \t which contain(s) the file(s) of cofactors decompositions.");
}

/************************************************************************/
/*     MAIN                                                             */
/************************************************************************/

int main(int argc, char *argv[])
{
    param_list pl;
    param_list_init(pl);
    declare_usage(pl);
    /* 
       Passing NULL is allowed here. Find value with
       param_list_parse_switch later on 
     */
    param_list_configure_switch(pl, "gdc", NULL);
    param_list_configure_switch(pl, "gst", NULL);
    param_list_configure_switch(pl, "gst_r", NULL);

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

    /*default values */
    int fbb1 = -1;
    int mfb1 = -1;
    int fbb0 = -1;
    int mfb0 = -1;
    int lpb0 = -1;
    int lpb1 = -1;
    int r = -1;

    param_list_parse_int(pl, "lpb0", &lpb0);
    param_list_parse_int(pl, "lpb1", &lpb1);
    param_list_parse_int(pl, "fbb1", &fbb1);
    param_list_parse_int(pl, "mfb1", &mfb1);
    param_list_parse_int(pl, "fbb0", &fbb0);
    param_list_parse_int(pl, "mfb0", &mfb0);
    param_list_parse_int(pl, "r", &r);

    //store data
    const char *directory_out;
    if ((directory_out = param_list_lookup_string(pl, "out")) == NULL) {
	directory_out = "./";
    }

    int gdc = param_list_parse_switch(pl, "-gdc");
    int gst_r = param_list_parse_switch(pl, "-gst_r");
    int gst = param_list_parse_switch(pl, "-gst");

    //todo: to use this option, it's necessary to include the library gsl!!!!
    //and precompute just for one side!
    if (gdc)			//precompute all decompositions!
    {
	fprintf(stdout, "This feature isn't avaliable. Wait few days!!!\n");
	exit(EXIT_SUCCESS);
	//todo: do this for only one side!!!
	/*   //check parameters */
	/*   if (fbb1 == -1 || mfb1 == -1 || */
	/*          fbb0 == -1 || mfb0 == -1)  */
	/*        { */
	/*          fprintf (stderr,"for this option '-gdc', you must also " */
	/*                   "specify these options:\t -fbb1, -mfb1, -fbb0, -mfb0\n"); */
	/*          exit (EXIT_FAILURE); */
	/*        } */
	/*   if (fbb1 !=fbb0) */
	/*        { */
	/*          precompute_tabular_decomp_files (2*fbb1-1, mfb1, fbb1); */
	/*          precompute_tabular_decomp_files (2*fbb0-1, mfb0, fbb0); */
	/*        } */
	/*   else */
	/*        { */
	/*          int max_cof = (mfb0 > mfb1)? mfb0 : mfb1; */
	/*          precompute_tabular_decomp_files (2*fbb0-1, max_cof, fbb0); */
	/*        } */

    } else if (gst) {
	//check parameters!

	const char *directory_in;
	if ((directory_in = param_list_lookup_string(pl, "in")) == NULL) {
	    fputs("Parser error: Please re-run with the option -in"
		  "with a directory which contains the precomputed files.\n",
		  stderr);
	    exit(EXIT_FAILURE);
	}


	int r0 = -1, r1 = -1;
	param_list_parse_int(pl, "r0", &r0);
	param_list_parse_int(pl, "r1", &r1);

	if (r0 != -1 && r1 != -1)
	    {
		ASSERT_ALWAYS(fbb0 != -1 && fbb1 != -1 );
		//just one couple will be studied!
		char name_file_in[strlen(directory_in) + 20];
		FILE * file_in;
		//get back the best strategies for r0!
		sprintf(name_file_in, "%s/strategies(%d)_%d",
			directory_in, fbb0, r0);
		file_in = fopen(name_file_in, "r");
		tabular_strategy_t *strat_r0 = tabular_strategy_fscan(file_in);
		if (strat_r0 == NULL) {
		    fprintf(stderr,
			    "Parser error: can't read the file '%s'\n",
			    name_file_in);
		    exit(EXIT_FAILURE);
		}
		fclose(file_in);

		//get back the best strategies for r1!
		sprintf(name_file_in, "%s/strategies(%d)_%d",
			directory_in, fbb1, r1);
		file_in = fopen(name_file_in, "r");
		tabular_strategy_t *strat_r1 = tabular_strategy_fscan(file_in);
		if (strat_r1 == NULL) {
		    fprintf(stderr,
			    "Parser error: can't read the file '%s'\n",
			    name_file_in);
		    exit(EXIT_FAILURE);
		}

		fclose(file_in);

		//compute the best strategies for (r0, r1)!
		char res_file[200];
		tabular_strategy_t *res =
		    generate_strategy_r0_r1(strat_r0, strat_r1);
		//print
		sprintf(res_file, 
			"%s/strategies_%d_%d", directory_out, r0, r1);
		FILE *file = fopen(res_file, "w");
		tabular_strategy_fprint(file, res);
		fclose(file);
		//free
		tabular_strategy_free (strat_r0);
		tabular_strategy_free (strat_r1);
		tabular_strategy_free (res);

	    } else {
	    ASSERT_ALWAYS(fbb0 != -1 && fbb1 != -1 &&
			  mfb0 != -1 && mfb1 != -1 &&
			  fbb0 < mfb0 && fbb1 < mfb1);
	    /* For the both side and each size of cofactor the best
	       strategies had been choosed. Now, it's time to merge these
	       datas to compute the best strategies for each couple (r0, r1)!
	    */

	    tabular_strategy_t **data_rat = malloc(sizeof(*data_rat) * (mfb0 + 1));
	    ASSERT_ALWAYS(data_rat);
	    char name_file_in[strlen(directory_in) + 20];
	    for (int r0 = 0; r0 <= mfb0; r0++) {
		sprintf(name_file_in, "%s/strategies(%d)_%d",
			directory_in, fbb0, r0);
		FILE *file_in = fopen(name_file_in, "r");
		data_rat[r0] = tabular_strategy_fscan(file_in);
		if (data_rat[r0] == NULL) {
		    fprintf(stderr,
			    "Parser error: impossible to read '%s'\n",
			    name_file_in);
		    exit(EXIT_FAILURE);
		}
		fclose(file_in);
	    }

	    for (int r1 = 0; r1 <= mfb1; r1++) {
		sprintf(name_file_in, "%s/strategies(%d)_%d",
			directory_in, fbb1, r1);
		FILE *file_in = fopen(name_file_in, "r");
		tabular_strategy_t *strat_r1 = tabular_strategy_fscan(file_in);
		if (strat_r1 == NULL) {
		    fprintf(stderr,
			    "Parser error: impossible read the file '%s'\n",
			    name_file_in);
		    exit(EXIT_FAILURE);
		}

		fclose(file_in);

		for (int r0 = 0; r0 <= mfb0; r0++) {
		    char res_file[200];
		    tabular_strategy_t *res =
			generate_strategy_r0_r1(data_rat[r0], strat_r1);
		    //print
		    sprintf(res_file, "%s/strategies_%d_%d", directory_out, r0, r1);
		    FILE *file = fopen(res_file, "w");
		    tabular_strategy_fprint(file, res);
		    fclose(file);

		    tabular_strategy_free(res);
		}
		tabular_strategy_free(strat_r1);
	    }
	    for (int r0 = 0; r0 <= mfb0; r0++)
		tabular_strategy_free(data_rat[r0]);
	    free(data_rat);
	}
    } else {

	const char *name_file_in;
	if ((name_file_in = param_list_lookup_string(pl, "in")) == NULL) {
	    fputs("Parser error: Please re-run with the option "
		  "-in and a valid file name.\n", stderr);
	    exit(EXIT_FAILURE);
	}


	FILE *file_in = fopen(name_file_in, "r");
	tabular_fm_t *c = tabular_fm_fscan(file_in);
	if (c == NULL) {
	    fprintf(stderr, "impossible to read %s\n", name_file_in);
	    param_list_clear(pl);
	    exit(EXIT_FAILURE);
	}
	fclose(file_in);

	tabular_fm_t *data_pp1_27 = extract_fm_method(c, PP1_27_METHOD, 0);
	tabular_fm_t *data_pp1_65 = extract_fm_method(c, PP1_65_METHOD, 0);

	tabular_fm_t *data_pm1 = extract_fm_method(c, PM1_METHOD, 0);

	tabular_fm_t *data_ecm_m16 = extract_fm_method(c, EC_METHOD, MONTY16);

	tabular_fm_t *data_ecm_m12 = extract_fm_method(c, EC_METHOD, MONTY12);

	tabular_fm_t *data_ecm_b12 = extract_fm_method(c, EC_METHOD, BRENT12);

	tabular_fm_t *data_ecm_rc = tabular_fm_create();
	tabular_fm_concat(data_ecm_rc, data_ecm_b12);
	tabular_fm_concat(data_ecm_rc, data_ecm_m12);

	tabular_fm_t *data_pp1 = tabular_fm_create();
	tabular_fm_concat(data_pp1, data_pp1_27);
	tabular_fm_concat(data_pp1, data_pp1_65);

	tabular_fm_sort(data_pm1);
	tabular_fm_sort(data_pp1);
	tabular_fm_sort(data_ecm_m16);
	tabular_fm_sort(data_ecm_rc);

	if (gst_r) {
	    ASSERT_ALWAYS(fbb0 < lpb0 && r != -1 && fbb0 != -1 && lpb0 != -1);

	    //precompute the convex hull for one bit size of cofactor
	    int lim = 2 * fbb0 - 1;
	    tabular_decomp_t *tab_decomp = NULL;
	    if (r >= lim) {
		const char *name_file_decomp;
		if ((name_file_decomp =
		     param_list_lookup_string(pl, "decomp")) == NULL) {
		    fputs("Parser error: Please re-run with the option "
			  "-decomp and a valid file name.\n", stderr);
		    exit(EXIT_FAILURE);
		}
		FILE *file_decomp = fopen(name_file_decomp, "r");

		tab_decomp = tabular_decomp_fscan(file_decomp);

		if (tab_decomp == NULL) {
		    fprintf(stderr, "impossible to read '%s'\n",
			    name_file_decomp);
		    exit(EXIT_FAILURE);
		}
		fclose(file_decomp);
	    }
	    fm_t *zero = fm_create();
	    unsigned long method_zero[4] = { 0, 0, 0, 0 };
	    fm_set_method(zero, method_zero, 4);

	    tabular_strategy_t *res =
		generate_strategies_oneside(tab_decomp, zero, data_pm1,
					    data_pp1, data_ecm_m16,
					    data_ecm_rc, fbb0, lpb0, r);

	    tabular_decomp_free(tab_decomp);

	    char name_file[strlen (directory_out) + 20];
	    sprintf(name_file, "%s/strategies(%d)_%d", directory_out, fbb0, r);
	    FILE *file_out = fopen(name_file, "w");
	    tabular_strategy_fprint(file_out, res);
	    fclose(file_out);
	    tabular_strategy_free(res);
	    fm_free(zero);

	} else {

	    ASSERT_ALWAYS(fbb1 < lpb1 && lpb1 < mfb1 &&
			  fbb0 < lpb0 && lpb0 < mfb0);
	    /*
	      computes the matrix of strategies where for each couple
	      (r0, r1) we will optain the best strategies to find a
	      relation.
	    */
	    const char *name_directory_decomp;
	    if ((name_directory_decomp =
		 param_list_lookup_string(pl, "decomp")) == NULL) {
		fputs("Parser error: Please re-run with the option "
		      "-decomp and a valid directory name.\n", stderr);
		exit(EXIT_FAILURE);
	    }

	    tabular_strategy_t ***matrix =
		generate_matrix(name_directory_decomp,
				data_pm1, data_pp1,
				data_ecm_m16, data_ecm_rc,
				fbb0, lpb0, mfb0,
				fbb1, lpb1,mfb1);

	    char res_file[strlen(directory_out) + 20];
	    for (int r0 = 0; r0 <= mfb0; r0++)
		for (int r1 = 0; r1 <= mfb1; r1++) {
		    sprintf(res_file, "%s/strategies_%d_%d", directory_out,
			    r0, r1);
		    FILE *file = fopen(res_file, "w");
		    tabular_strategy_fprint(file, matrix[r0][r1]);
		    fclose(file);
		}

	    //free
	    for (int r0 = 0; r0 <= mfb0; r0++) {
		for (int r1 = 0; r1 <= mfb1; r1++)
		    tabular_strategy_free(matrix[r0][r1]);
		free(matrix[r0]);
	    }
	    free(matrix);
	}
	//free
	tabular_fm_free(data_pp1_65);
	tabular_fm_free(data_pp1_27);
	tabular_fm_free(data_pp1);
	tabular_fm_free(data_pm1);
	tabular_fm_free(data_ecm_m16);
	tabular_fm_free(data_ecm_m12);
	tabular_fm_free(data_ecm_b12);
	tabular_fm_free(data_ecm_rc);
	tabular_fm_free(c);
    }
    param_list_clear(pl);

    return EXIT_SUCCESS;
}
