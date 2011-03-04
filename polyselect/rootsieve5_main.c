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

void rootsieve_stdin ( int, int );
void rootsieve_file ( FILE *file, int, int );

/*
  Usage
*/
static void
usage (char **argv)
{
	 fprintf (stderr, "# Error: Unexpected argument: %s\n", argv[1]);
	 fprintf (stderr, "# Usage: \n# %s < STDIN -w WLEFT -l WLENGTH\n", argv[0]);
	 fprintf (stderr, "# %s -f FILE -w WLEFT -l WLENGTH\n", argv[0]);
	 fprintf (stderr, "# Note: WLEFT is the left-bound of a qudratic rotation and WLENGTH is the amount for qudratic rotation.\n");
	 exit(1);
}


/*
  Main
*/
int
main (int argc, char **argv)
{

	 int i, w_left_bound = 0, w_length = 1;
	 /* print command-line arguments */
	 fprintf (stderr, "# %s.r%s", *argv, CADO_REV);
	 for (i = 1; i < argc; i++)
		  fprintf (stderr, " %s", *(argv+i));
	 fprintf (stderr, "\n");

	 /* read polynomials from "-f file" */
	 if (argc == 7) {

		  FILE *file = NULL;
		  char *filename = NULL;

		  while (argc >= 3 && argv[1][0] == '-') {
			   if (argc >= 3 && strcmp(argv[1], "-f") == 0) {
					filename = argv[2];
					argv += 2;
					argc -= 2;
			   }
			   else if (argc >= 3 && strcmp (argv[1], "-w") == 0)
			   {
					w_left_bound = atoi (argv[2]);
					argv += 2;
					argc -= 2;
			   }
			   else if (argc >= 3 && strcmp (argv[1], "-l") == 0)
			   {
					w_length = atoi (argv[2]);
					argv += 2;
					argc -= 2;
			   }
			   else {
					usage (argv);
			   }
		  }

		  // read
		  file = fopen(filename, "r");
		  if (file == NULL) {
			   fprintf(stderr, "# Error in reading file\n");
			   exit (1);
		  }

		  // optimize the raw polys in the file
		  rootsieve_file (file, w_left_bound, w_length);
		  fclose (file);
		  return 0;
	 }

	 /* read polynomials from stdin*/
	 else if (argc == 5)
	 {
		  while (argc >= 3 && argv[1][0] == '-') {
			   if (argc >= 3 && strcmp (argv[1], "-w") == 0)
			   {
					w_left_bound = atoi (argv[2]);
					argv += 2;
					argc -= 2;
			   }
			   else if (argc >= 3 && strcmp (argv[1], "-l") == 0)
			   {
					w_length = atoi (argv[2]);
					argv += 2;
					argc -= 2;
			   }
			   else {
					usage (argv);
			   }
		  }

		  rootsieve_stdin (w_left_bound, w_length);
		  return 0;
	 }
	 else {
		  usage (argv);
	 }
}
