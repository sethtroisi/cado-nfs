/*
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
  Therefore, function "rsbound_setup_AB_bound()" limits the sieving
  range for "v". You may change this in "ropt_param.h".

  Even the range "v" is limited, it may take much memory.
  Therefore, we divide it into "SIEVEARRAY_SIZE " in "ropt_param.h".
  Setting it larger take more memory and may reduce (setup) time.

  For each such "SIEVEARRAY_SIZE", we actually sieving in blocks
  of "L1_SIZE 12288". You may change this in "ropt_param.h".

  Another important parameter:  "TOPALPHA_EACH_SIEVEARRAY 16".
  For each sieve array of "SIEVEARRAY_SIZE", compute the MurphyE of
  16 polynomials whose alpha values are the best among this array.

  These poynomials will be then filtered into another global priority
  queue which records "TOPE_EACH_SUBLATTICE=8" polynomials with top
  MurphyE.

  Note that, for polynomials of small skewness. Size can be more
  important, hence you may want to set "TOPALPHA_EACH_SIEVEARRAY 8"
  larger. However, this may reduce the performance since MurphyE
  computation is slow.

  [4. Tuning]

  The "ropt_param.c" contains functions: "rsbound_setup_AB_bound()"
  and  "rsparam_setup ()". They can be tuned. Please also see its
  header file for other parameters.
  
  [5. Log]

  -- (Dec, 2010) block sieving.
  -- (Dec, 2010) reduced memory usage in return_all_sublattices().
  -- (Jan, 2011) addded priority queue, changed precion in return_all_sublattices().
  -- (Feb, 2011) some tunnings, correct bugs in rootsieve_v().
  -- (Mar, 2011) readme.
  -- (Apr, 2011) changed 1d sieving to 2d.
  -- (Apr, 2012) spitting into separate c files.

  [6. Bugs]

  Please report bugs to Shi Bai (shi.bai AT anu.edu.au).
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
#include "ropt_param.h"
#include "ropt_arith.h"
#include "ropt_io.h"
#include "ropt.h"

/*
  Do the root sieve on all polynomials in the file.
*/
void
ropt_file_cado ( FILE *file,
                 param_t param )
{
  unsigned flag = 0UL, count = 0;
  char str[MAX_LINE_LENGTH];

  /* rootsieve_struct */
  ropt_poly_t rs;
  rsstr_init (rs);

  /* for each polynomial, do the root sieve. */
  while (1) {

    /* read poly */
    if (fgets(str, MAX_LINE_LENGTH, file) == NULL)
      break;

    if ( str[0] == 'Y' ) {
      if ( str[1] == '1' ) {
        gmp_sscanf (str, "Y1: %Zd\n", rs->g[1]);
        (flag) ^= (1<<8);
      }
      else if ( str[1] == '0' ) {
        gmp_sscanf (str, "Y0: %Zd\n", rs->g[0]);
        (flag) ^= (1<<7);
      }
      else {
        fprintf (stderr, "Error in parsing line %s.\n", str);
        exit (1);
      }
    }
    else if ( str[0] == 'c') {
      if ( str[1] == '0' ) {
        gmp_sscanf (str, "c0: %Zd\n", rs->f[0]);
        (flag) ^= 1;
      }
      else if ( str[1] == '1' ) {
        gmp_sscanf (str, "c1: %Zd\n", rs->f[1]);
        (flag) ^= (1<<1);
      }
      else if ( str[1] == '2' ) {
        gmp_sscanf (str, "c2: %Zd\n", rs->f[2]);
        (flag) ^= (1<<2);
      }
      else if ( str[1] == '3' ) {
        gmp_sscanf (str, "c3: %Zd\n", rs->f[3]);
        (flag) ^= (1<<3);
      }
      else if ( str[1] == '4' ) {
        gmp_sscanf (str, "c4: %Zd\n", rs->f[4]);
        (flag) ^= (1<<4);
      }
      else if ( str[1] == '5' ) {
        gmp_sscanf (str, "c5: %Zd\n", rs->f[5]);
        (flag) ^= (1<<5);
      }
      else if ( str[1] == '6' ) {
        gmp_sscanf (str, "c6: %Zd\n", rs->f[6]);
        (flag) ^= (1<<6);
      }
      else
      {
        fprintf (stderr, "Error in parsing line %s.\n", str);
        exit (1);
      }
    }
    else if ( str[0] == 'n') {
      gmp_sscanf (str, "n: %Zd\n", rs->n);
      (flag) ^= (1<<9);
    }

    else
      continue;

    if (flag == 1023UL || flag == 959UL) {

      /* pre-compute and setup rs */
      rsstr_setup (rs);

      bestpoly_t bestpoly;
      bestpoly_init (bestpoly, rs->d);

      /* set best poly */
      int i;
      for (i = 0; i <= rs->d; i++)
      {
        mpz_set (bestpoly->f[i], rs->f[i]);
      }
      mpz_set (bestpoly->g[0], rs->g[0]);
      mpz_set (bestpoly->g[1], rs->g[1]);

      fprintf (stderr, "\n# Polynomial (# %5d).\n", count);
      print_poly_fg (rs->f, rs->g, rs->d, rs->n, 2);

      /* start main rootsieve function */
      ropt (rs, bestpoly, param, 2);
      bestpoly_free (bestpoly, rs->d);

      count += 1;
      flag = 0UL;
    }
  }

  /* free */
  rsstr_free (rs);
}

/*
  Do the root sieve on all polynomials in the file in msieve format.
*/
void
ropt_file_msieve ( FILE *file,
                   param_t param )
{
  unsigned count = 0;
  char str[MAX_LINE_LENGTH];

  /* rootsieve_struct */
  ropt_poly_t rs;
  rsstr_init (rs);

  mpz_t ad, l, m;
  mpz_init (ad);
  mpz_init (l);
  mpz_init (m);

  /* for each polynomial, do the root sieve. */
  while (1) {

    /* read poly */   
    if (fgets(str, MAX_LINE_LENGTH, file) == NULL)
      break;
    int ret = gmp_sscanf (str, "%Zd %Zd %Zd", ad, l, m);
    if (ret != 3)
      continue;

    /* set degree, n, ad, l, m */
    rs->d = param->d;
    mpz_set (rs->n, param->n);
    mpz_set (rs->f[rs->d], ad);
    Lemma21 (rs->f, rs->n, rs->d, l, m);
    mpz_set (rs->g[1], l);
    mpz_neg (rs->g[0], m);
    optimize (rs->f, rs->d, rs->g, 0, 1);
    fprintf (stderr, "# Polynomial (# %5d).\n", count);

#if SKIP_ROOTSIEVE_M
    print_poly_info_short (rs->f, rs->g, rs->d, rs->n);
#else
    print_poly_fg (rs->f, rs->g, rs->d, rs->n, 2);

    rsstr_setup (rs);
    bestpoly_t bestpoly;
    bestpoly_init (bestpoly, rs->d);
    int i;
    for (i = 0; i <= rs->d; i++)
    {
      mpz_set (bestpoly->f[i], rs->f[i]);
    }
    mpz_set (bestpoly->g[0], rs->g[0]);
    mpz_set (bestpoly->g[1], rs->g[1]);

    ropt (rs, bestpoly, param, 2);
    bestpoly_free (bestpoly, rs->d);
#endif
    count += 1;
  }

  /* free */
  mpz_clear (ad);
  mpz_clear (l);
  mpz_clear (m);
  rsstr_free (rs);
}

/*
  Do the root sieve for one polynomial from stdin.
  For debugging purpose.
*/
void
ropt_stdin ( param_t param )
{
  /* rootsieve_struct */
  ropt_poly_t rs;
  rsstr_init (rs);

  /* read poly to rs */
  read_ggnfs (rs->n, rs->f, rs->g, rs->m);

  /* pre-compute and setup rs */
  rsstr_setup (rs);

  bestpoly_t bestpoly;
  bestpoly_init (bestpoly, rs->d);

  /* set best poly */
  int i;
  for (i = 0; i <= rs->d; i++)
  {
    mpz_set (bestpoly->f[i], rs->f[i]);
  }
  mpz_set (bestpoly->g[0], rs->g[0]);
  mpz_set (bestpoly->g[1], rs->g[1]);

  fprintf (stderr, "\n# Polynomial (# 0).\n");
  print_poly_fg (rs->f, rs->g, rs->d, rs->n, 2);

  /* start main rootsieve function */
  ropt (rs, bestpoly, param, 2);

  rsstr_free (rs);
  bestpoly_free (bestpoly, rs->d);
}

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
  fprintf (stderr, " -w W                   W is the lower bound of the quadratic rotation\n");
  fprintf (stderr, " -l L                   L is the length of the quadratic rotation\n");
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

  fprintf (stderr, "Rootsieve with quadratic rotation between -512 and 512 and within\nthe sieving region bounded by norm 70.\n");

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
    ropt_file_cado (file, param);
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
    ropt_file_msieve (file, param);
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
    ropt_stdin (param);
    param_clear (param);
    return 0;
  }
  else
    usage (argv);
}
