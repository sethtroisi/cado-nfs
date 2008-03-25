/* Convert between different polynomial formats
   (CADO, Franke-Kleinjung, GGNFS, CWI).
   First version written by Paul Zimmermann, 20030804.
   Rewritten to allow different input/output formats, 20071220.

   Call with --help or any invalid option for usage information.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmp.h"

#define MAX_DEGREE 6

void
read_franke (mpz_t N, mpz_t *X, mpz_t Y1, mpz_t Y0, mpz_t M)
{
  int c;

  if (mpz_inp_str (N, stdin, 0) == 0)
    {
      fprintf (stderr, "Error while reading N\n");
      exit (1);
    }
  c = getchar ();
  if (c != '\n')
    {
      fprintf (stderr, "Error: end of line expected after N\n");
      exit (1);
    }
  while ((c = getchar ()) == '#')
    {
      while ((c = getchar ()) != '\n');
    }
  if (c != 'X' || getchar () != '5' || getchar () != ' ')
    {
      fprintf (stderr, "Error: X5 expected\n");
      exit (1);
    }
  mpz_inp_str (X[5], stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X5\n");
      exit (1);
    }
  if (getchar () != 'X' || getchar () != '4' || getchar () != ' ')
    {
      fprintf (stderr, "Error: X4 expected\n");
      exit (1);
    }
  mpz_inp_str (X[4], stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X4\n");
      exit (1);
    }
  if (getchar () != 'X' || getchar () != '3' || getchar () != ' ')
    {
      fprintf (stderr, "Error: X3 expected\n");
      exit (1);
    }
  mpz_inp_str (X[3], stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X3\n");
      exit (1);
    }
  if (getchar () != 'X' || getchar () != '2' || getchar () != ' ')
    {
      fprintf (stderr, "Error: X2 expected\n");
      exit (1);
    }
  mpz_inp_str (X[2], stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X2\n");
      exit (1);
    }
  if (getchar () != 'X' || getchar () != '1' || getchar () != ' ')
    {
      fprintf (stderr, "Error: X1 expected\n");
      exit (1);
    }
  mpz_inp_str (X[1], stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X1\n");
      exit (1);
    }
  if (getchar () != 'X' || getchar () != '0' || getchar () != ' ')
    {
      fprintf (stderr, "Error: X0 expected\n");
      exit (1);
    }
  mpz_inp_str (X[0], stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X0\n");
      exit (1);
    }
  if (getchar () != 'Y' || getchar () != '1' || getchar () != ' ')
    {
      fprintf (stderr, "Error: Y1 expected\n");
      exit (1);
    }
  mpz_inp_str (Y1, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after Y1\n");
      exit (1);
    }
  if (getchar () != 'Y' || getchar () != '0' || getchar () != ' ')
    {
      fprintf (stderr, "Error: Y0 expected\n");
      exit (1);
    }
  mpz_inp_str (Y0, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after Y0\n");
      exit (1);
    }
  if (getchar () != 'M' || getchar () != ' ')
    {
      fprintf (stderr, "Error: M expected\n");
      exit (1);
    }
  mpz_inp_str (M, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after M\n");
      exit (1);
    }
}

void
readline (void)
{
  int c;

  while ((c = getchar ()) != '\n');
}

void
readstring (char *s)
{
  int i, n;
  char *t;
  
  n = strlen (s);
  t = malloc (n * sizeof (char));
  for (i = 0; i < n; i++)
    t[i] = getchar ();
  if (strncmp (s, t, n) != 0)
    {
      fprintf (stderr, "Error, expected %s, got %s\n", s, t);
      exit (1);
    }
  free (t);
}

void
read_ggnfs (mpz_t N, mpz_t *X, mpz_t Y1, mpz_t Y0, mpz_t M)
{
  int i, ret;
  char s[100]; /* input buffer */

  while (feof (stdin) == 0)
    {
      ret = scanf ("%s ", s);
      if (ret != 1)
        break;
      if (strcmp (s, "n:") == 0) /* n: input number */
        {
          if (mpz_inp_str (N, stdin, 0) == 0)
            {
              fprintf (stderr, "Error while reading N: ");
              exit (1);
            }
          if (getchar () != '\n')
            {
              fprintf (stderr, "Error: end of line expected after N\n");
              exit (1);
            }
        }
      else if (sscanf (s, "c%d:", &i) == 1) /* ci: coeff of degree i */
        {
          if (i > MAX_DEGREE)
            {
              fprintf (stderr, "Error, too large degree %d\n", i);
              exit (1);
            }
          mpz_inp_str (X[i], stdin, 0);
          if (getchar () != '\n')
            {
              fprintf (stderr, "Error: end of line expected after c%d\n", i);
              exit (1);
            }
        }
      else if (strcmp (s, "Y1:") == 0)
        {
          mpz_inp_str (Y1, stdin, 0);
          if (getchar () != '\n')
            {
              fprintf (stderr, "Error: end of line expected after Y1\n");
              exit (1);
            }
        }
      else if (strcmp (s, "Y0:") == 0)
        {
          mpz_inp_str (Y0, stdin, 0);
          if (getchar () != '\n')
            {
              fprintf (stderr, "Error: end of line expected after Y0\n");
              exit (1);
            }
        }
      else if (strcmp (s, "#") == 0 && scanf ("%s ", s) == 1 &&
               (strcmp (s, "M") == 0 || strcmp (s, "m:") == 0))
        {
          mpz_inp_str (M, stdin, 0);
          if (getchar () != '\n')
            {
              fprintf (stderr, "Error: end of line expected after M\n");
              exit (1);
            }
        }
      else
        readline ();
    }
}

void
out_cwi (mpz_t N, mpz_t *X, int deg, mpz_t Y1, mpz_t Y0, mpz_t M)
{
  int i;

  mpz_out_str (stdout, 10, N);
  printf ("\n");
  mpz_out_str (stdout, 10, M);
  printf ("\n\n2\n\n1\n");
  mpz_out_str (stdout, 10, Y0);
  printf (" ");
  mpz_out_str (stdout, 10, Y1);
  printf ("\n\n%d\n", deg);
  for (i = 0; i <= deg; i++)
    {
      mpz_out_str (stdout, 10, X[i]);
      printf (" ");
    }
  printf ("\n\n10000000\n");
}

void
out_franke (mpz_t N, mpz_t *X, int deg, mpz_t Y1, mpz_t Y0, mpz_t M)
{
  int i;

  mpz_out_str (stdout, 10, N);  printf ("\n");
  for (i = deg; i >= 0; i--)
    {
      printf ("X%d ", i);
      mpz_out_str (stdout, 10, X[i]);
      printf ("\n");
    }
  printf ("Y1 "); mpz_out_str (stdout, 10, Y1); printf ("\n");
  printf ("Y0 "); mpz_out_str (stdout, 10, Y0); printf ("\n");
  printf ("M ");  mpz_out_str (stdout, 10, M);  printf ("\n");
  printf ("0 2000000 2.7 27 54\n");
  printf ("0 4000000 2.7 27 54\n");
}

void
usage (char *argv0)
{
  fprintf (stderr, "Usage: %s [options] < file1 > file2\n", argv0);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "        -if xxx - specify input format\n");
  fprintf (stderr, "        -of yyy - specify output format\n");
  fprintf (stderr, "Available formats:\n");
  fprintf (stderr, "        cado  - CADO-NFS format (default)\n");
  fprintf (stderr, "        fk    - Franke-Kleinjung format\n");
  fprintf (stderr, "        ggnfs - GGNFS format\n");
  fprintf (stderr, "        cwi - CWI format\n");
}

#define FORMAT_CADO  0
#define FORMAT_FK    1
#define FORMAT_GGNFS 2
#define FORMAT_CWI   3

int
main (int argc, char *argv[])
{
  mpz_t N, *X, Y1, Y0, M;
  int iformat = FORMAT_CADO;
  int oformat = FORMAT_CADO;
  int i, deg;
  
  while (argc > 1)
    {
      if (argc >= 3 && strcmp (argv[1], "-if") == 0)
        {
          if (strcmp (argv[2], "cado") == 0)
            iformat = FORMAT_CADO;
          else if (strcmp (argv[2], "fk") == 0)
            iformat = FORMAT_FK;
          else if (strcmp (argv[2], "ggnfs") == 0)
            iformat = FORMAT_GGNFS;
          else if (strcmp (argv[2], "cwi") == 0)
            iformat = FORMAT_CWI;
          else
            {
              fprintf (stderr, "Error, invalid format: %s\n", argv[2]);
              exit (1);
            }
          argc -= 2;
          argv += 2;
        }
      else if (argc >= 3 && strcmp (argv[1], "-of") == 0)
        {
          if (strcmp (argv[2], "cado") == 0)
            oformat = FORMAT_CADO;
          else if (strcmp (argv[2], "fk") == 0)
            oformat = FORMAT_FK;
          else if (strcmp (argv[2], "ggnfs") == 0)
            oformat = FORMAT_GGNFS;
          else if (strcmp (argv[2], "cwi") == 0)
            oformat = FORMAT_CWI;
          else
            {
              fprintf (stderr, "Error, invalid format: %s\n", argv[2]);
              exit (1);
            }
          argc -= 2;
          argv += 2;
        }
      else
        {
          usage (argv[0]);
          exit (1);
        }
    }

  mpz_init (N);
  X = (mpz_t*) malloc ((MAX_DEGREE + 1) * sizeof (mpz_t));
  for (i = 0; i <= MAX_DEGREE; i++)
    mpz_init (X[i]);
  mpz_init (Y1);
  mpz_init (Y0);
  mpz_init (M);

  if (iformat == FORMAT_FK)
    read_franke (N, X, Y1, Y0, M);
  else if (iformat == FORMAT_GGNFS || iformat == FORMAT_CADO)
    read_ggnfs (N, X, Y1, Y0, M);
  else
    {
      fprintf (stderr, "Input format not yet implemented\n");
      exit (1);
    }

  for (deg = MAX_DEGREE; deg > 0 && mpz_cmp_ui (X[deg], 0) == 0; deg --);

  if (oformat == FORMAT_CWI)
    out_cwi (N, X, deg, Y1, Y0, M);
  else if (oformat == FORMAT_FK)
    out_franke (N, X, deg, Y1, Y0, M);
  else
    {
      fprintf (stderr, "Output format not yet implemented\n");
      exit (1);
    }

  mpz_clear (M);
  mpz_clear (Y0);
  mpz_clear (Y1);
  for (i = 0; i <= MAX_DEGREE; i++)
    mpz_clear (X[i]);
  mpz_clear (N);
  free (X);

  return 0;
}
