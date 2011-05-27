/*
  Main to call the rootsieve5.c, root sieve for degree 5 and degree 6 polynomials.
*/
#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "auxiliary.h"
#include "murphyE.h"
#include "rho.h"
#include "rootsieve5.h"

void rootsieve_stdin ( param_t param);
void rootsieve_file ( FILE *file, param_t param);
void rootsieve_file_msieve ( FILE *file, param_t param);
void param_init (param_t param);
void param_clear (param_t param);

/*
  Usage
*/
static void
usage (char **argv)
{
	 fprintf (stderr, "Error: Unexpected argument: %s\n", argv[1]);
	 fprintf (stderr, "\nUsage: %s -f F -w W -l L [OPTION]\n", argv[0]);
	 fprintf (stderr, "Usage (STDIN): %s -w W -l L [OPTION]\n", argv[0]);

	 fprintf (stderr, "\nRequired options\n");
	 fprintf (stderr, " -f F                   F is the poly file in CADO format\n");
	 fprintf (stderr, " -fm F                  F is the poly file in msieve format, must also have parameters -n and -d\n");
	 fprintf (stderr, " -w W                   W is the lower bound of the qudratic rotation\n");
	 fprintf (stderr, " -l L                   L is the length of the qudratic rotation\n");
	 fprintf (stderr, "\nOther options:\n");
	 fprintf (stderr, " -e N E_1 E_2 ... E_N   N is the number of small primes in finding sublattices\n");
	 fprintf (stderr, "                        E_1 to E_N are the prime powers in finding sublattices\n");
	 fprintf (stderr, " -umax U                U and -U are the bounds of the sieving region for linear rotation\n");
	 fprintf (stderr, " -vmax V                V and -V are the bounds of the sieving region for constant rotation\n");
	 fprintf (stderr, " -norm M                M is the (estimated) lognorm upper bound in root sieve\n");
	 fprintf (stderr, " -n N                   N is integer (only in -fm mode)\n");
	 fprintf (stderr, " -d D                   D is the degree of polynomial (only in -fm mode)\n");
	 fprintf (stderr, " --s2                   Sieve only mode, usage is explained in Example 3\n");

	 fprintf (stderr, "\nExample 1: %s -f rsa768.poly -w -512 -l 1024  -norm 70\n", argv[0]);

	 fprintf (stderr, "Rootsieve with qudratic rotation between -512 and 512 and within\nthe sieving region bounded by norm 70.\n");

	 fprintf (stderr, "\nExample 2: %s -f rsa768.poly -w -512 -l 1024 -e 5 7 4 3 2 2 -umax 16 -vmax 10000000\n", argv[0]);

	 fprintf (stderr, "As above, but we request five prime factors (2, 3, 5, 7, 11) in\nthe sublattice with powers (7, 4, 3, 2, 2).\n");

	 fprintf (stderr, "\nExample 3: %s -f rsa768.poly --s2 -w 12 -u 345 -v 6789 -mod 1814400 -umax 16 -vmax 10000000\n", argv[0]);

	 fprintf (stderr, "Special use. Assume that we know the polynomial has good root \nproperty at rotation (12*x^2 + 345*x + 6789), we want to search\n12*x^2 + (345 + 1814400*i)*x + (6789 + 1814400*j)*x for i, j\nbounded by -umax and -vmax.\n");

	 fprintf (stderr, "\nExample 4: %s -fm rsa768.poly -w -512 -l 1024 -e 5 7 4 3 2 2 -umax 16 -vmax 10000000 -n N -d 6\n", argv[0]);

	 fprintf (stderr, "Read msieve format where each line is (ad l m). The same \nparameter as in Example 2, but need to have -n integer \nand -d degree parameters.\n");

	 fprintf (stderr, "\nUse Example 2 as a starting point for a raw polynomial. If the \nnumber is larger, use larger integers in \"-e N E_1 E_2 ... E_N\". \n");
	 exit(1);
}


/*
  Main
*/
int
main (int argc, char **argv)
{
	 int i;
	 param_t param;

	 /* print command-line arguments */
	 fprintf (stderr, "# %s.r%s", *argv, CADO_REV);
	 for (i = 1; i < argc; i++)
		  fprintf (stderr, " %s", *(argv+i));
	 fprintf (stderr, "\n");
	 param_init (param);

	 /* read polynomials in cado format */
	 if (argc >= 5 && strcmp(argv[1], "-f") == 0) {

		  FILE *file = NULL;
		  char *filename = NULL;
		  filename = argv[2];
		  argv += 2;
		  argc -= 2;

		  /* manually s2 parameters */
		  if (strcmp (argv[1], "--s2") == 0) {
			   param->flag = 2;
			   argv += 1;
			   argc -= 1;
			   while (argc >= 3 && argv[1][0] == '-') {

					if (argc >= 3 && strcmp (argv[1], "-w") == 0)
					{
						 param->s2_w = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-u") == 0)
					{
						 mpz_set_str (param->s2_u, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-v") == 0)
					{
						 mpz_set_str (param->s2_v, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-mod") == 0)
					{
						 mpz_set_str (param->s2_mod, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-umax") == 0)
					{
						 param->s2_Amax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-vmax") == 0)
					{
						 param->s2_Bmax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-n") == 0)
					{
						 mpz_set_str (param->n, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-norm") == 0)
					{
						 param->lognorm_bound = atof (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else {
						 usage (argv);
					}
			   }
		  }
		  /* auto mode and/or use manually s1 parameters */
		  else {
			   param->flag = 0;
			   while (argc >= 3 && argv[1][0] == '-') {

					if (argc >= 3 && strcmp (argv[1], "-w") == 0)
					{
						 param->w_left_bound = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-l") == 0)
					{
						 param->w_length = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-umax") == 0)
					{
						 param->s2_Amax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-vmax") == 0)
					{
						 param->s2_Bmax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-n") == 0)
					{
						 mpz_set_str (param->n, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-norm") == 0)
					{
						 param->lognorm_bound = atof (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-e") == 0)
					{
						 param->flag = 1;
						 param->s1_num_e_sl = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;

						 int j = 0;
						 for (j = 0; j < param->s1_num_e_sl; j ++) {
							  param->s1_e_sl[j] = (unsigned short) atoi (argv[1]);
							  argv += 1;
							  argc -= 1;
						 }
					}
					else {
						 usage (argv);
					}
			   }
		  }

		  file = fopen(filename, "r");
		  if (file == NULL) {
			   fprintf(stderr, "# Error in reading file\n");
			   exit (1);
		  }

		  /* root sieve raw polys in the file */
		  rootsieve_file (file, param);
		  fclose (file);
		  param_clear (param);

		  return 0;
	 }

	 /* read polynomials in msieve format */
	 if (argc >= 5 && strcmp(argv[1], "-fm") == 0) {

		  FILE *file = NULL;
		  char *filename = NULL;
		  filename = argv[2];
		  argv += 2;
		  argc -= 2;

		  /* manually s2 parameters */
		  if (strcmp (argv[1], "--s2") == 0) {
			   param->flag = 2;
			   argv += 1;
			   argc -= 1;
			   while (argc >= 3 && argv[1][0] == '-') {

					if (argc >= 3 && strcmp (argv[1], "-w") == 0)
					{
						 param->s2_w = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-u") == 0)
					{
						 mpz_set_str (param->s2_u, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-v") == 0)
					{
						 mpz_set_str (param->s2_v, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-mod") == 0)
					{
						 mpz_set_str (param->s2_mod, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-umax") == 0)
					{
						 param->s2_Amax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-vmax") == 0)
					{
						 param->s2_Bmax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-n") == 0)
					{
						 mpz_set_str (param->n, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-norm") == 0)
					{
						 param->lognorm_bound = atof (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-d") == 0)
					{
						 param->d = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else {
						 usage (argv);
					}
			   }
		  }
		  /* auto mode and/or use manually s1 parameters */
		  else {
			   param->flag = 0;
			   while (argc >= 3 && argv[1][0] == '-') {

					if (argc >= 3 && strcmp (argv[1], "-w") == 0)
					{
						 param->w_left_bound = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-l") == 0)
					{
						 param->w_length = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-umax") == 0)
					{
						 param->s2_Amax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-vmax") == 0)
					{
						 param->s2_Bmax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-n") == 0)
					{
						 mpz_set_str (param->n, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-norm") == 0)
					{
						 param->lognorm_bound = atof (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-d") == 0)
					{
						 param->d = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-e") == 0)
					{
						 param->flag = 1;
						 param->s1_num_e_sl = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;

						 int j = 0;
						 for (j = 0; j < param->s1_num_e_sl; j ++) {
							  param->s1_e_sl[j] = (unsigned short) atoi (argv[1]);
							  argv += 1;
							  argc -= 1;
						 }
					}
					else {
						 usage (argv);
					}
			   }
		  }

		  file = fopen(filename, "r");
		  if (file == NULL) {
			   fprintf(stderr, "# Error in reading file\n");
			   exit (1);
		  }

		  if ( mpz_cmp_ui(param->n, 0) == 0 ) {
			   fprintf(stderr, "# Error: please input parameter \"-n number\"\n");
			   exit (1);
		  }
		  if ( param->d == 0 ) {
			   fprintf(stderr, "# Error: please input parameter \"-n degree\"\n");
			   exit (1);
		  }


		  /* root sieve raw polys in the file */
		  rootsieve_file_msieve (file, param);
		  fclose (file);
		  param_clear (param);

		  return 0;
	 }

	 /* read polynomials from stdin */
	 else if (argc >= 5)
	 {

		  /* manually s2 parameters */
		  if (strcmp (argv[1], "--s2") == 0) {
			   argv += 1;
			   argc -= 1;
			   param->flag = 2;
			   while (argc >= 3 && argv[1][0] == '-') {

					if (argc >= 3 && strcmp (argv[1], "-w") == 0)
					{
						 param->s2_w = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-u") == 0)
					{
						 mpz_set_str (param->s2_u, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-v") == 0)
					{
						 mpz_set_str (param->s2_v, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-mod") == 0)
					{
						 mpz_set_str (param->s2_mod, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-umax") == 0)
					{
						 param->s2_Amax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-vmax") == 0)
					{
						 param->s2_Bmax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-n") == 0)
					{
						 mpz_set_str (param->n, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-norm") == 0)
					{
						 param->lognorm_bound = atof (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else {
						 usage (argv);
					}
			   }
		  }
		  /* auto mode and/or use manually s1 parameters */
		  else {
			   param->flag = 0;
			   while (argc >= 3 && argv[1][0] == '-') {

					if (argc >= 3 && strcmp (argv[1], "-w") == 0)
					{
						 param->w_left_bound = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-l") == 0)
					{
						 param->w_length = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-umax") == 0)
					{
						 param->s2_Amax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-vmax") == 0)
					{
						 param->s2_Bmax = atol (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-n") == 0)
					{
						 mpz_set_str (param->n, argv[2], 10);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-norm") == 0)
					{
						 param->lognorm_bound = atof (argv[2]);
						 argv += 2;
						 argc -= 2;
					}
					else if (argc >= 3 && strcmp (argv[1], "-e") == 0)
					{
						 param->flag = 1;
						 param->s1_num_e_sl = atoi (argv[2]);
						 argv += 2;
						 argc -= 2;

						 int j = 0;
						 for (j = 0; j < param->s1_num_e_sl; j ++) {
							  param->s1_e_sl[j] = (unsigned short) atoi (argv[1]);
							  argv += 1;
							  argc -= 1;
						 }
					}
					else {
						 usage (argv);
					}
			   }
		  }

		  /* root sieve raw polys in the file */
		  rootsieve_stdin (param);
		  param_clear (param);
		  return 0;
	 }
	 else
		  usage (argv);

}
