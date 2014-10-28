#include <math.h>
#include <float.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>

#include "facul.h"
#include "ecm.h"
#include "cado.h"
#include "portability.h"
#include "utils.h"
#include "tab_strategy.h"
#include "generate_strategies.h"

/************************************************************************/
/*             USAGE   + CHECK_PARAMETERS                              */
/************************************************************************/

static void declare_usage(param_list pl)
{
    param_list_usage_header(pl,
			    "This binary allows to build the best strategies for each couple (r0, r1)\n"
			    "where (r0,r1) are the bits size for our couple of cofactors.\n");

    param_list_decl_usage(pl, "afbb",
			  "set algebraic factor base bound to 2^afbb\n");
    param_list_decl_usage(pl, "rfbb",
			  "set rationnal factor base bound to 2^rfbb\n");
    param_list_decl_usage(pl, "rlpb",
			  "set rational large prime bound to 2^rlpb");
    param_list_decl_usage(pl, "alpb",
			  "set algebraic large prime bound to 2^alpb");
    param_list_decl_usage(pl, "rmfb", "set rational cofactor bound to 2^rmfb");
    param_list_decl_usage(pl, "amfb", "set algebraic cofactor bound to 2^amfb");
    param_list_decl_usage(pl, "in", "to locate the file which contains "
			  "our factoring methods.");
    param_list_decl_usage(pl, "out", "to locate the directory where the "
			  "file(s) will be stored.");

    param_list_decl_usage(pl, "gdc",
			  "(switch)  to precompute all decompositions of cofactors in the \n "
			  "sieve region. So, you must specify the sieve region with these options:\n"
			  "\t-afb, -amfb, -rfb, -rmfb\n");
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
    int afbb = -1;
    int amfb = -1;
    int rfbb = -1;
    int rmfb = -1;
    param_list_parse_int(pl, "afbb", &afbb);
    param_list_parse_int(pl, "amfb", &amfb);
    param_list_parse_int(pl, "rfbb", &rfbb);
    param_list_parse_int(pl, "rmfb", &rmfb);

    int gdc = param_list_parse_switch(pl, "-gdc");

    //todo: to use this option, it's necessary to include the library gsl!!!!
    if (gdc)			//precompute all decompositions!
    {
	fprintf(stdout, "This feature isn't avalaible. Wait few days!!!\n");
	exit(EXIT_SUCCESS);
	/*   //check parameters */
	/*   if (afbb == -1 || amfb == -1 || */
	/*          rfbb == -1 || rmfb == -1)  */
	/*        { */
	/*          fprintf (stderr,"for this option '-gdc', you must also " */
	/*                   "specify these options:\t -afbb, -amfb, -rfbb, -rmfb\n"); */
	/*          exit (EXIT_FAILURE); */
	/*        } */
	/*   if (afbb !=rfbb) */
	/*        { */
	/*          precompute_tabular_decomp_files (2*afbb-1, amfb, afbb); */
	/*          precompute_tabular_decomp_files (2*rfbb-1, rmfb, rfbb); */
	/*        } */
	/*   else */
	/*        { */
	/*          int max_cof = (rmfb > amfb)? rmfb : amfb; */
	/*          precompute_tabular_decomp_files (2*rfbb-1, max_cof, rfbb); */
	/*        } */
    } else {
	int rlpb = -1;
	int alpb = -1;
	param_list_parse_int(pl, "rlpb", &rlpb);
	param_list_parse_int(pl, "alpb", &alpb);

	ASSERT_ALWAYS(afbb < alpb && alpb < amfb && rfbb < rlpb && rlpb < rmfb);

	const char *name_file_in;
	if ((name_file_in = param_list_lookup_string(pl, "in")) == NULL) {
	    if (name_file_in != NULL)
		printf("%s\n", name_file_in);
	    fputs("Parser error: Please re-run with the option "
		  "-in with a valid file name.\n", stderr);
	    exit(EXIT_FAILURE);
	}

	tabular_fm_t *c = tabular_fm_create();

	FILE *file_in = fopen(name_file_in, "r");
	int err = fm_fscan(file_in, c);
	if (err < 0) {
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

	printf("all: (%d)\n", c->index);
	printf("pp1: (%d)\n", data_pp1->index);
	printf("pm1: (%d)\n", data_pm1->index);
	printf("ecm_m16: (%d)\n", data_ecm_m16->index);
	printf("ecm_rc: (%d)\n", data_ecm_rc->index);

	//todo: subdivid in several sequences!!!!
	tabular_strategy_t ***matrix = generate_matrix(data_pm1, data_pp1,
						       data_ecm_m16,
						       data_ecm_rc,
						       rfbb, rlpb,
						       rmfb, afbb, alpb, amfb);

	//store data
	const char *directory_out;
	if ((directory_out = param_list_lookup_string(pl, "out")) == NULL) {
	    directory_out = "./";
	}

	char res_file[200];
	for (int r1 = 0; r1 <= rmfb; r1++)	//todo: kill all
	    for (int r2 = 0; r2 <= amfb; r2++) {
		sprintf(res_file, "%s/strategies_%d_%d", directory_out, r1, r2);
		FILE *file = fopen(res_file, "w");
		tabular_strategy_print_file(matrix[r1][r2], file);
		fclose(file);
	    }

	//free
	for (int r1 = 0; r1 <= rmfb; r1++) {
	    for (int r2 = 0; r2 <= amfb; r2++)
		tabular_strategy_free(matrix[r1][r2]);
	    free(matrix[r1]);
	}
	free(matrix);

	//free
	printf("free!!\n");
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
