/**
   Root optimization for polynomials in the number field sieve.

   [1. Run]

   Input polynomials are in "CADO" or "Msieve" format. They are read 
   from STDIN/FILE. Output polynomials are in "CADO" format. Examples:

   ropt < POLY -w 0 -l 8 > OUTPUT 2> ERRINFO
   ropt -f POLY -w 0 -l 8 > OUTPUT 2> ERRINFO

   If deg=5, quadratic rotations "-w" and "-l" will be ignored.
   If deg=6, "-w" is the leftmost of a quadratic rotation and 
   "-l" is the step from the "-w".

   [2. Algorithm]

   There are two steps in the root sieve.

   Given a polynomial, the code first tries to find good
   (u, v) (mod (p1^e1 * p2^e2* \cdots *pn^en)) for small prime powers.
   This step is done by hensel-lift-like procedure in a p^2-ary tree data
   structure (for each such p) and discard any bad/impossible (u, v)
   during the tree-building. The multiple roots are updated in a way
   following Emmanuel Thome's idea. Finally, Finally we use CRT to
   find some good (u, v) pairs, those with small u. Note if there are
   too many pi^ei, the CRTs will domimated the running time.

   In the second step, we do the actually root sieve for all primes
   up to some bound. We only consider the (u, v) points lying on the
   sublattice within the sieving range (U, V).

   The method is described in "Root optimization of polynomials in 
   the number field sieve." on http://wwwmaths.anu.edu.au/~bai.

   [3. Degree 6]

   For degree 6 polynomial, the following processes are called in order,

   (a) - For each quadratic rotation w;
   (b) -- Tune parameters; ("WANT_TUNE" in "ropt_param.h".)
   (c) -- Find good sublattices (w, u, v)
   (d) - Compare (priority queue) all good sublattices and pick up top ones.
   (e) - For each such sublattice (w, u, v)
   (f) -- Do the root sieve.

   Details:
   (a) In the following command,

   ropt < POLY -w 0 -l 8

   "-w" defines the leftmost point of quadratic rotation.
   f'(x) = f(x) + w*x^2*g(x)

   "-l" defines the steps for quadratic rotation, (w+l-1)

   (b, c)
   The code starts by looking each quadratic rotation, say
   [w, \cdots, w+l-1]. For each rotated polynomial, we will
   have to find a suitable set of parameters (p1^e1 * p2^e2*
   \cdots *pn^en) to produce the sublattice. Note that, if it 
   is larger, then we could find polynomials with better alpha,
   but probably worse size; reversly, a small parameter gives 
   worse alpha, but probably better size. It is hard to prebuilt
   a universel parameter. Hence we may tune the parameters (p1^e1 
   * p2^e2* \cdots *pn^en) using a trial sieving. The tune function
   is disabled by default. Toggle "WANT_TUNE" in "ropt_param.h".

   (d)
   After good sublattices (w, u, v) are found for all w in the
   permitted range (as you set by "-w" and "-l"), we will compare
   the alpha values between all these sublattices. At the moment,
   we only pick up the top "rsparam->nbest_sl = 128" sublattices
   for sieving. You may change this in "ropt_param.h".

   (e, f)
   For each survived sublattice (w, u, v), we do the root sieve.
   The permitted bound for "v" could be huge which runs forever.
   Therefore, function "ropt_s2param_setup_range()" limits the sieving
   range for "v". You may change this in "ropt_param.h".

   Even the range "v" is limited, it may take much memory.
   Therefore, we divide it into "SIEVEARRAY_SIZE " in "ropt_param.h".
   Setting it larger take more memory and may reduce (setup) time.

   For each such "SIEVEARRAY_SIZE", we actually sieving in blocks
   of "L1_SIZE 12288". You may change this in "ropt_param.h".

   Another important parameter:  "NUM_TOPALPHA_SIEVEARRAY 16".
   For each sieve array of "SIEVEARRAY_SIZE", compute the MurphyE of
   16 polynomials whose alpha values are the best among this array.

   These poynomials will be then filtered into another global priority
   queue which records "NUM_TOPE_SUBLATTICE=8" polynomials with top
   MurphyE.

   Note that, for polynomials of small skewness. Size can be more
   important, hence you may want to set "NUM_TOPALPHA_SIEVEARRAY 8"
   larger. However, this may reduce the performance since MurphyE
   computation is slow.

   [4. Tuning]

   The "ropt_param.c" contains functions: "ropt_s2param_setup_range()"
   and  "ropt_param_setup ()". They can be tuned. Please also see its
   header file for other parameters.
  
   [5. Log]

   -- (Dec, 2010) block sieving.
   -- (Dec, 2010) reduced memory usage in return_combined_sublattice().
   -- (Jan, 2011) addded priority queue, changed precision in return_combined_sublattice().
   -- (Feb, 2011) some tunnings, correct bugs in rootsieve_v().
   -- (Mar, 2011) readme.
   -- (Apr, 2011) changed 1d sieving to 2d.
   -- (Apr, 2012) spitting into separate c files.
   -- (May, 2012) rewrite.

   [6. Bugs]

   Please report bugs to Shi Bai (shih.bai AT gmail.com).
*/


#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "rho.h"
#include "auxiliary.h"
#include "murphyE.h"
#include "ropt.h"
#include "portability.h"
#include "area.h"

/**
 * Usage
 */
void
usage (char **argv)
{
  fprintf (stderr, "Error: Unhandled parameter: %s\n", argv[1]);
  fprintf (stderr, "\n");
  fprintf (stderr, "Usage: %s -f fname [options]\n", argv[0]);
  fprintf (stderr, "       %s -f fname --s2 [options]\n", argv[0]);
  fprintf (stderr, "       %s -fm fname [options]\n", argv[0]);
  fprintf (stderr, "       %s -fm fname --s2 [options]\n", argv[0]);
  fprintf (stderr, "       %s [options]\n", argv[0]);
  fprintf (stderr, "       %s --s2 [options]\n", argv[0]);
  fprintf (stderr, "\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, " -f fname      read polynomials in CADO format.\n");
  fprintf (stderr, " -fm fname     read polynomials in MSIEVE format (use together\n");
  fprintf (stderr, "                with options -n and -d).\n");
  fprintf (stderr, " -n N          (only in -fm mode) the integer to be factored.\n");
  fprintf (stderr, " -d D          (only in -fm mode) the degree of the polynomial.\n");
  fprintf (stderr, " -amin L       L is the lower bound for quadratic rotation.\n");
  fprintf (stderr, " -amax R       R is the upper bound for quadratic rotation.\n");
  fprintf (stderr, " -bmax U       U and -U are the sieving length for linear rotation.\n");
  fprintf (stderr, " -cmax V       V and -V are the sieving length for constant rotation.\n");
  fprintf (stderr, " -e N E_1 E_2 ... E_N\n");
  fprintf (stderr, "               N is the number of small primes in finding sublattices.\n");
  fprintf (stderr, "               E_1 to E_N are the prime powers in finding sublattices.\n");
  fprintf (stderr, " -norm M       M is the (estimated) lognorm upper bound in ropt.\n");
  fprintf (stderr, " -v            number of occurrences defines verbose level {0, 1, 2, 3} (default 0).\n");
  fprintf (stderr, " --s2          (switch) sieve-only mode (use together with the following\n");
  fprintf (stderr, "               four options: -a, -b, -c and -mod.\n");
  fprintf (stderr, " -a A          Fix the quadratic rotation of the sublattice by A.\n");
  fprintf (stderr, " -b B          Fix the linear rotation of the sublattice by B.\n");
  fprintf (stderr, " -c C          Fix the constant rotation of the sublattice by C.\n");
  fprintf (stderr, " -mod M        M is the sublattice modulus.\n");
  fprintf (stderr, " -mod M        M is the sublattice modulus.\n");
  fprintf (stderr, " -Bf F         algebraic smoothness bound (default %.2e).\n", BOUND_F);
  fprintf (stderr, " -Bg G         rational smoothness bound (default %.2e).\n", BOUND_G);
  fprintf (stderr, " -area A       sieving area (default %.2e).\n", AREA);
  fprintf (stderr, " -rseffort M   sieving effort ranging from 1 to 5 (default %d).\n", DEFAULT_RSEFFORT);

  fprintf (stderr, "\nExample 1: %s -f fname\n", argv[0]);
  fprintf (stderr, "Root optimization for all CADO-formatted polynomials in 'fname'.\n");

  fprintf (stderr, "\nExample 2: %s -f fname -amin -512 -amax 512  -norm 70\n", argv[0]);
  fprintf (stderr, "As above, but restricts the quadratic rotation between -512\n"
           "and 512 and the sieving region by norm 70.\n");

  fprintf (stderr, "\nExample 3: %s -f fname -amin -512 -amax 512 -e 5 7 4 3 2 2 -bmax 16 -cmax 10000000\n", argv[0]);
  fprintf (stderr, "As above, but uses five prime factors (2, 3, 5, 7, 11) in the\n"
           "sublattice with powers (7, 4, 3, 2, 2). It also tells ropt to\n"
           "root sieve a region of 16 by 10000000.\n");

  fprintf (stderr, "\nExample 4: %s -fm rsa768.poly -amin -512 -amax 512 -e 5 7 4 3 2 2 -bmax 16 -cmax 10000000 -n $N -d 6\n", argv[0]);
  fprintf (stderr, "As above, but reads msieve format where each line contains\n"
           "'c_d Y1 Y0'. The parameters -n and -d are compulsory.\n");

  fprintf (stderr, "\nExample 5: %s -f fname --s2 -a 12 -b 345 -c 6789 -mod 1814400 -bmax 16 -cmax 10000000\n", argv[0]);
  fprintf (stderr, "Sieve-only mode. Assume that we know the polynomial has good\n"
           "root property at rotation (12*x^2 + 345*x + 6789), we want to\n"
           "search 12*x^2 + (345 + 1814400*i)*x + (6789 + 1814400*j) where\n"
           "i, j are bounded by -bmax and -cmax.\n");
  exit(1);
}


/**
 * parse manually input parameters to param.
 */
static void
ropt_parse_param ( int argc,
                   char **argv,
                   ropt_param_t param )
{

  /* filename only */
  if (argc > 1) {
    /* stage 2 (root sieve) parameters only */
    if (strcmp (argv[1], "--s2") == 0) {

      param->stage_flag = 2;
      argv += 1;
      argc -= 1;
      while (argc >= 2 && argv[1][0] == '-') {
        if (strcmp (argv[1], "-v") == 0)
        {
          param->verbose ++;
          argv += 1;
          argc -= 1;
        }
        else if (argc >= 3 && strcmp (argv[1], "-a") == 0)
        {
          param->s2_w = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-b") == 0)
        {
          mpz_set_str (param->s2_u, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-c") == 0)
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
        else if (argc >= 3 && strcmp (argv[1], "-bmax") == 0)
        {
          param->s2_Amax = atol (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-cmax") == 0)
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
          param->bound_lognorm = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-d") == 0)
        {
          param->d = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-Bf") == 0)
        {
          bound_f = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-Bg") == 0)
        {
          bound_g = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-area") == 0)
        {
          area = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-rseffort") == 0)
        {
          param->effort = atoi (argv[2]);
          argv += 2;
          argc -= 2;
          if (param->effort < 1 || param->effort > 5) {
            fprintf (stderr, "Error: -rseffort not in range (1-5).\n");
            exit(1);
          }       
        }
        else {
          usage (argv);
        }
      }
    }
    /* stage 1 and stage 2 parameters */
    else {
      while (argc >= 2 && argv[1][0] == '-') {
        if (strcmp (argv[1], "-v") == 0)
        {
          param->verbose ++;
          argv += 1;
          argc -= 1;
        }
        else if (argc >= 3 && strcmp (argv[1], "-amin") == 0)
        {
          param->w_left_bound = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-amax") == 0)
        {
          param->w_length = atoi (argv[2]) - param->w_left_bound + 1;
          if (param->w_length < 0) {
            fprintf (stderr, "Error in options -amin and/or -amax.\n");
            exit(1);
          }
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-bmax") == 0)
        {
          param->s2_Amax = atol (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-cmax") == 0)
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
          param->bound_lognorm = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-d") == 0)
        {
          param->d = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-Bf") == 0)
        {
          bound_f = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-Bg") == 0)
        {
          bound_g = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-area") == 0)
        {
          area = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-rseffort") == 0)
        {
          param->effort = atoi (argv[2]);
          argv += 2;
          argc -= 2;
          if (param->effort < 1 || param->effort > 5) {
            fprintf (stderr, "Error: -rseffort not in range (1-5).\n");
            exit(1);
          }
        }
        else if (argc >= 3 && strcmp (argv[1], "-e") == 0)
        {
          param->stage_flag = 1;
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
  }
}


/**
 * Main function: call ropt_on_cadopoly() or ropt_on_cadopoly()
 * or ropt_on_stdin().
 */
int
main (int argc, char **argv)
{
  int i;

  /* L1 cache size */
  ropt_L1_cachesize ();

  ropt_param_t param;
  ropt_param_init (param);

  /* print command-line arguments */
  fprintf (stderr, "# %s.r%s", *argv, CADO_REV);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", *(argv+i));
  fprintf (stderr, "\n");

  /* read polynomials in cado format */
  if (argc >= 3 && strcmp(argv[1], "-f") == 0) {

    FILE *file = NULL;
    char *filename = NULL;
    filename = argv[2];
    argv += 2;
    argc -= 2;

    file = fopen(filename, "r");
    if (file == NULL) {
      fprintf(stderr, "# Error in reading file\n");
      exit (1);
    }

    /* parse parameters */
    ropt_parse_param (argc, argv, param);

    /* call ropt_on_cadopoly() */
    ropt_on_cadopoly (file, param);

    fclose (file);
    ropt_param_free (param);

    return 0;
  }
  /* read polynomials in msieve format */
  else if (argc >= 3 && strcmp(argv[1], "-fm") == 0) {

    FILE *file = NULL;
    char *filename = NULL;
    filename = argv[2];
    argv += 2;
    argc -= 2;

    file = fopen(filename, "r");
    if (file == NULL) {
      fprintf(stderr, "# Error in reading file\n");
      exit (1);
    }

    /* parse parameters */
    ropt_parse_param (argc, argv, param);

    if (mpz_cmp_ui(param->n, 0) == 0) {
      fprintf(stderr, "# Error: please input parameter \"-n number\"\n");
      exit (1);
    }
    if ( param->d == 0 ) {
      fprintf(stderr, "# Error: please input parameter \"-n degree\"\n");
      exit (1);
    }

    /* call ropt_on_msievepoly() */
    ropt_on_msievepoly (file, param);

    fclose (file);
    ropt_param_free (param);

    return 0;
  }
  /* read polynomials from stdin */
  else if (argc >= 3)
  {

    /* parse parameters */
    ropt_parse_param (argc, argv, param);

    /* call ropt_on_stdin() */
    ropt_on_stdin (param);

    ropt_param_free (param);

    return 0;
  }
  else
    usage (argv);
}
