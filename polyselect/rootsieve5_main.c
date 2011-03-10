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
void param_init (param_t param);
void param_clear (param_t param);

/*
  Usage
*/
static void
usage (char **argv)
{
	 fprintf (stderr, "# Error: Unexpected argument: %s\n", argv[1]);
	 fprintf (stderr, "# Usage 1 (a): %s -f FILE -w WLEFT -l WLENGTH\n", argv[0]);
	 fprintf (stderr, "# Usage 1 (b): %s -f FILE -w WLEFT -l WLENGTH -E N E_1 E_2 ... E_N \n", argv[0]);
	 fprintf (stderr, "# Parameters: \n#\t\"WLEFT\" is the left-bound of a qudratic rotation; \n#\t\"WLENGTH\" is the amount for qudratic rotation; \n#\t\"N\" and \"e_i\" are the number of sublattice primes and exponents. \n");
	 fprintf (stderr, "# Usage 2: %s -f FILE --s2 -w W -u U -v V -mod MOD -umax MAX2 -vmax MAX1\n", argv[0]);
	 fprintf (stderr, "# Parameters: \n#\t\"W, U, V, MOD\" defines a sublattice; \n#\t\"MAX2, MAX1\" are the rotation bounds;\n");
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

	 /* read polynomials from "-f file" */
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
