/* Convert Franke-Kleinjung's or GGNFS polynomial file into CWI format.
   Written by Paul Zimmermann, 20030804
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmp.h"

void
read_franke (mpz_t N, mpz_t X5, mpz_t X4, mpz_t X3, mpz_t X2, mpz_t X1,
	     mpz_t X0, mpz_t Y1, mpz_t Y0, mpz_t M)
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
  mpz_inp_str (X5, stdin, 0);
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
  mpz_inp_str (X4, stdin, 0);
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
  mpz_inp_str (X3, stdin, 0);
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
  mpz_inp_str (X2, stdin, 0);
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
  mpz_inp_str (X1, stdin, 0);
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
  mpz_inp_str (X0, stdin, 0);
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
read_ggnfs (mpz_t N, mpz_t X5, mpz_t X4, mpz_t X3, mpz_t X2, mpz_t X1,
	    mpz_t X0, mpz_t Y1, mpz_t Y0, mpz_t M)
{
  int c;

  readline (); /* name: ... */
  readstring ("n: ");
  if (mpz_inp_str (N, stdin, 0) == 0)
    {
      fprintf (stderr, "Error while reading N\n");
      exit (1);
    }
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after N\n");
      exit (1);
    }
  readline (); /* skew: ... */
  readline (); /* # norm ... */
  readstring ("c5: ");
  mpz_inp_str (X5, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X5\n");
      exit (1);
    }
  readstring ("c4: ");
  mpz_inp_str (X4, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X4\n");
      exit (1);
    }
  readstring ("c3: ");
  mpz_inp_str (X3, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X3\n");
      exit (1);
    }
  readstring ("c2: ");
  mpz_inp_str (X2, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X2\n");
      exit (1);
    }
  readstring ("c1: ");
  mpz_inp_str (X1, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X1\n");
      exit (1);
    }
  readstring ("c0: ");
  mpz_inp_str (X0, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after X0\n");
      exit (1);
    }
  readline (); /* # alpha ... */
  readstring ("Y1: ");
  mpz_inp_str (Y1, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after Y1\n");
      exit (1);
    }
  readstring ("Y0: ");
  mpz_inp_str (Y0, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after Y0\n");
      exit (1);
    }
  readline (); /* # Murphy_E ... */
  readstring ("# M ");
  mpz_inp_str (M, stdin, 0);
  if (getchar () != '\n')
    {
      fprintf (stderr, "Error: end of line expected after M\n");
      exit (1);
    }
}

void
out_cwi (mpz_t N, mpz_t X5, mpz_t X4, mpz_t X3, mpz_t X2, mpz_t X1,
	 mpz_t X0, mpz_t Y1, mpz_t Y0, mpz_t M)
{
  mpz_out_str (stdout, 10, N);
  printf ("\n");
  mpz_out_str (stdout, 10, M);
  printf ("\n\n2\n\n1\n");
  mpz_out_str (stdout, 10, Y0);
  printf (" ");
  mpz_out_str (stdout, 10, Y1);
  printf ("\n\n5\n");
  mpz_out_str (stdout, 10, X0); printf (" ");
  mpz_out_str (stdout, 10, X1); printf (" ");
  mpz_out_str (stdout, 10, X2); printf (" ");
  mpz_out_str (stdout, 10, X3); printf (" ");
  mpz_out_str (stdout, 10, X4); printf (" ");
  mpz_out_str (stdout, 10, X5);
  printf ("\n\n10000000\n");
}

void
out_franke (mpz_t N, mpz_t X5, mpz_t X4, mpz_t X3, mpz_t X2, mpz_t X1,
	    mpz_t X0, mpz_t Y1, mpz_t Y0, mpz_t M)
{
  mpz_out_str (stdout, 10, N);  printf ("\n");
  printf ("X5 "); mpz_out_str (stdout, 10, X5); printf ("\n");
  printf ("X4 "); mpz_out_str (stdout, 10, X4); printf ("\n");
  printf ("X3 "); mpz_out_str (stdout, 10, X3); printf ("\n");
  printf ("X2 "); mpz_out_str (stdout, 10, X2); printf ("\n");
  printf ("X1 "); mpz_out_str (stdout, 10, X1); printf ("\n");
  printf ("X0 "); mpz_out_str (stdout, 10, X0); printf ("\n");
  printf ("Y1 "); mpz_out_str (stdout, 10, Y1); printf ("\n");
  printf ("Y0 "); mpz_out_str (stdout, 10, Y0); printf ("\n");
  printf ("M ");  mpz_out_str (stdout, 10, M);  printf ("\n");
  printf ("0 2000000 2.7 27 54\n");
  printf ("0 4000000 2.7 27 54\n");
}

int
main (int argc, char *argv[])
{
  mpz_t N, X5, X4, X3, X2, X1, X0, Y1, Y0, M;
  int ggnfs = 0;

  if (argc >= 2 && strcmp (argv[1], "-ggnfs") == 0)
    {
      ggnfs = 1;
      argc --;
      argv ++;
    }

  mpz_init (N);
  mpz_init (X5);
  mpz_init (X4);
  mpz_init (X3);
  mpz_init (X2);
  mpz_init (X1);
  mpz_init (X0);
  mpz_init (Y1);
  mpz_init (Y0);
  mpz_init (M);

  if (ggnfs == 0)
    read_franke (N, X5, X4, X3, X2, X1, X0, Y1, Y0, M);
  else
    read_ggnfs (N, X5, X4, X3, X2, X1, X0, Y1, Y0, M);

  if (ggnfs == 0)
    out_cwi (N, X5, X4, X3, X2, X1, X0, Y1, Y0, M);
  else
    out_franke (N, X5, X4, X3, X2, X1, X0, Y1, Y0, M);

  mpz_clear (M);
  mpz_clear (Y0);
  mpz_clear (Y1);
  mpz_clear (X0);
  mpz_clear (X1);
  mpz_clear (X2);
  mpz_clear (X3);
  mpz_clear (X4);
  mpz_clear (X5);
  mpz_clear (N);

  return 0;
}
