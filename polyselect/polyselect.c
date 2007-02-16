/* polyselect1 - polynomial selection (naive implementation)

Author: Paul Zimmermann

Copyright 2007 INRIA

This file is part of the CADO project.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

   Usage: polyselect < c80 > c80.poly

   References:

   [1] "On polynomial selection for the general number field sieve",
       Thorsten Kleinjung, Mathematics of Computation 75 (2006), p. 2037-2047.
   [2] "Integer factorization, part 4: polynomial selection", Dan Bernstein,
       invited lecture, Arizona Winter School, University of Arizona, Tucson,
       Arizona, 14 March 2006.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h> /* for ULONG_MAX */
#include <values.h> /* for DBL_MAX */
#include <math.h>   /* for log, fabs */
#include "cado.h"

static void
usage ()
{
  fprintf (stderr, "Usage: polyselect [-v] [-e nnn] < in > out\n\n");
  fprintf (stderr, "       -v     - verbose\n");
  fprintf (stderr, "       -e nnn - use effort nnn\n");
  fprintf (stderr, "       in     - input file (number to factor)\n");
  fprintf (stderr, "       out    - output file (polynomials)\n");
  exit (1);
}

static void
init (cado_input in)
{
  mpz_init (in->n);
  in->degree = 0;
}

static void
clear (cado_input in)
{
  mpz_clear (in->n);
}

struct sd {
  size_t s;
  unsigned long d;
};

static struct sd default_degrees[DEFAULT_DEGREES_LENGTH] = DEFAULT_DEGREES;

/* default degree, when not given by user */
static int
default_degree (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_degrees[i].s; i++);
  return default_degrees[i].d;
}

static struct sd default_rlim_table[DEFAULT_RLIM_LENGTH] = DEFAULT_RLIM;
static struct sd default_alim_table[DEFAULT_ALIM_LENGTH] = DEFAULT_ALIM;

/* default degree, when not given by user */
static unsigned long
default_rlim (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_rlim_table[i].s; i++);
  return default_rlim_table[i].d;
}

/* default degree, when not given by user */
static unsigned long
default_alim (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_alim_table[i].s; i++);
  return default_alim_table[i].d;
}

static struct sd default_lpbr_table[DEFAULT_LPBR_LENGTH] = DEFAULT_LPBR;
static struct sd default_lpba_table[DEFAULT_LPBA_LENGTH] = DEFAULT_LPBA;

/* default rational large prime bound, when not given by user */
static unsigned long
default_lpbr (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_lpbr_table[i].s; i++);
  return default_lpbr_table[i].d;
}

/* default algebraic large prime bound, when not given by user */
static unsigned long
default_lpba (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_lpba_table[i].s; i++);
  return default_lpba_table[i].d;
}

static struct sd default_mfbr_table[DEFAULT_MFBR_LENGTH] = DEFAULT_MFBR;
static struct sd default_mfba_table[DEFAULT_MFBA_LENGTH] = DEFAULT_MFBA;

/* default rational factor-residual bound, when not given by user */
static unsigned long
default_mfbr (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_mfbr_table[i].s; i++);
  return default_mfbr_table[i].d;
}

/* default algebraic factor-residual bound, when not given by user */
static unsigned long
default_mfba (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_mfba_table[i].s; i++);
  return default_mfba_table[i].d;
}

struct sf {
  size_t s;
  double f;
};

static struct sf default_rlambda_table[DEFAULT_RLAMBDA_LENGTH] = DEFAULT_RLAMBDA;
static struct sf default_alambda_table[DEFAULT_ALAMBDA_LENGTH] = DEFAULT_ALAMBDA;

/* default rational lambda, when not given by user */
static double
default_rlambda (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_rlambda_table[i].s; i++);
  return default_rlambda_table[i].f;
}

/* default algebraic lambda, when not given by user */
static double
default_alambda (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_alambda_table[i].s; i++);
  return default_alambda_table[i].f;
}

static struct sd default_qint_table[DEFAULT_QINT_LENGTH] = DEFAULT_QINT;

/* default rational factor-residual bound, when not given by user */
static unsigned long
default_qint (mpz_t n)
{
  size_t s = mpz_sizeinbase (n, 10);
  int i;

  for (i = 0; s > default_qint_table[i].s; i++);
  return default_qint_table[i].d;
}

static void
init_poly (cado_poly p, cado_input in)
{
  int i;

  strcpy (p->name, "");
  mpz_init_set (p->n, in->n);
  p->degree = in->degree;
  if (p->degree == 0)
    p->degree = default_degree (p->n);
  p->f = (mpz_t*) malloc ((p->degree + 1) * sizeof (mpz_t));
  for (i = 0; i <= p->degree; i++)
    mpz_init (p->f[i]);
  p->g = (mpz_t*) malloc (2 * sizeof (mpz_t));
  mpz_init (p->g[0]);
  mpz_init (p->g[1]);
  mpz_init (p->m);
  strcpy (p->type, "gnfs");
}

static void
clear_poly (cado_poly p)
{
  int i;

  mpz_clear (p->n);
  for (i = 0; i <= p->degree; i++)
    mpz_clear (p->f[i]);
  free (p->f);
  mpz_clear (p->g[0]);
  mpz_clear (p->g[1]);
  free (p->g);
  mpz_clear (p->m);
}

static void
parse_error (int c)
{
  if (c == 0)
    {
      fprintf (stderr, "Error, no field read\n");
      exit (1);
    }
  if (c == EOF)
    {
      fprintf (stderr, "Error, end of file\n");
      exit (1);
    }
}

/* read a string of characters [a-z]+, returns the number of characters read */
static int
read_string (char *s)
{
  int c, l = 0;
  
  while (1)
    {
      c = getchar ();
      if (c == EOF)
	return c;
      if (islower (c) == 0)
	{
	  ungetc (c, stdin);
	  s[l] = '\0';
	  return l;
	}
      s[l++] = c;
    }
}

static void
parse_input (cado_input in)
{
  int c;
  char s[256];

  while ((c = getchar ()) != EOF)
    {
      /* skip blank characters and lines */
      while (isspace (c) && c != EOF)
	c = getchar ();

      if (c == EOF)
	break;
      
      if (c != '#') /* comment line */
	{
	  ungetc (c, stdin);
	  c = read_string (s);
	  if (c == 0 || c == EOF)
	    parse_error (c);

	  /* read possible spaces or tabs */
	  while (isblank (c = getchar ()))
	    printf ("read '%c'\n", c);

	  /* read ':' */
	  if (c != ':')
	    {
	      fprintf (stderr, "Error, ':' expected\n");
	      exit (1);
	    }

	  if (strcmp (s, "n") == 0)
	    {
	      if (mpz_inp_str (in->n, stdin, 10) == 0)
		{
		  fprintf (stderr, "Error after n:\n");
		  exit (1);
		}
	    }
	  else if (strcmp (s, "deg") == 0)
	    {
	      c = scanf ("%d", &(in->degree));
	      if (c == 0 || c == EOF)
		parse_error (c);
	    }
	  else
	    {
	      fprintf (stderr, "Error, unrecognized field: %s\n", s);
	      exit (1);
	    }
	}

      /* read to end of line */
      while (!feof (stdin) && c != '\n')
	c = getchar ();
    }

  if (mpz_cmp_ui (in->n, 0) == 0)
    {
      fprintf (stderr, "Error, no input number read\n");
      exit (1);
    }
}

/* round to nearest, assume d > 0 */
static void
mpz_ndiv_qr (mpz_t q, mpz_t r, mpz_t n, mpz_t d)
{
  int s;
  mpz_fdiv_qr (q, r, n, d); /* round towards -inf */
  mpz_mul_2exp (r, r, 1);
  s = mpz_cmp (r, d);
  mpz_div_2exp (r, r, 1);
  if (s > 0)
    {
      mpz_add_ui (q, q, 1);
      mpz_sub (r, r, d);
    }
}

/* L-Inf skewness, using Definition 3.1 from reference [1]:
   S = min_{s > 0} max_{i} |a_i s^{i-d/2}|
*/
double
skewness (mpz_t *p, unsigned long degree)
{
  unsigned long i, j, k;
  double s, smin, S, *loga;

  loga = (double*) malloc ((degree + 1) * sizeof (double));
  for (i = 0; i <= degree; i++)
    loga[i] = log (fabs (mpz_get_d (p[i])));
  smin = DBL_MAX;
  for (i = 0; i < degree; i++)
    for (j = i + 1; j <= degree; j++)
      {
	/* compute line going through log(ai) and log(aj),
	   and check it goes over the other points */
	s = (loga[j] - loga[i]) / (double) (j - i);
	/* tentative convex line is log(ai) + s*(x-i) */
	for (k = 0; k <= degree; k++)
	  if (k != i && k != j && loga[k] > loga[i] + s * (double) (k - i))
	    break;
	if (k > degree && s < smin)
	  {
	    smin = s;
	    S = exp (loga[i] + s * (double) (i - (double) degree / 2.0));
	  }
      }
  free (loga);
  return exp (-smin);
}

/* mu(m,S), as defined in reference [2]:
   
   mu(m,S) = (m/S+S)*(|a_d S^d| + ... + |a_0 S^{-d}|)
*/
double
get_mu (mpz_t *p, unsigned long degree, mpz_t m, double S)
{
  double mu = 0.0;
  double S2 = S * S;
  unsigned long i;

  for (i = 0; i <= degree; i++)
    mu = mu * S2 + fabs (mpz_get_d (p[i]));
  for (i = 0; i <= degree; i++)
    mu /= S;
  return (mpz_get_d (m) / S + S) * mu;
}

void
generate_base_m (cado_poly p, mpz_t m)
{
  unsigned long i;

  mpz_ndiv_qr (p->f[1], p->f[0], p->n, m);
  for (i = 1; i < p->degree; i++)
    mpz_ndiv_qr (p->f[i+1], p->f[i], p->f[i], m);
  mpz_neg (p->g[0], m);
  mpz_set_ui (p->g[1], 1);
  mpz_set (p->m, m);
  p->skew = skewness (p->f, p->degree);
}

void
generate_sieving_parameters (cado_poly out)
{
  out->rlim = default_rlim (out->n);
  out->alim = default_alim (out->n);
  out->lpbr = default_lpbr (out->n);
  out->lpba = default_lpba (out->n);
  out->mfbr = default_mfbr (out->n);
  out->mfba = default_mfba (out->n);
  out->rlambda = default_rlambda (out->n);
  out->alambda = default_alambda (out->n);
  out->qintsize = default_qint (out->n);
}

/* very naive method, just to check the interface with other modules */
void
generate_poly_0 (cado_poly out)
{
  unsigned long d = out->degree;

  /* in base m, coefficients are bounded by m/2, thus we want the leading
     coefficient to be at most m/2 too */
  mpz_mul_2exp (out->m, out->n, 1);
  mpz_root (out->m, out->m, d + 1);
  generate_base_m (out, out->m);
  generate_sieving_parameters (out);
}

/* first method described by Bernstein in reference [2], slide 6:
   Take m \approx B^{1/(d-1)} n^{1/(d+1)}
   Then f_d \approx B^{-d/(d-1)} n^{1/(d+1)}
   and f_0...f_{d-1} are \approx B^{1/(d-1)} n^{1/(d+1)} by default.
   We can hope f_0...f_{d-1} \approx B^{-d/(d-1)} n^{1/(d+1)}
   after B^{d(d+1)/(d-1)} tries.
   Then mu is decreased by a factor B.

   T is the number of tries, thus B = T^{(d-1)/d/(d+1)}.
*/
void
generate_poly (cado_poly out, unsigned long T, int verbose)
{
  unsigned long d = out->degree;
  double B, mu, best_mu = DBL_MAX;
  mpz_t best_m;

  mpz_init (best_m);
  B = pow ((double) T, 1.0 / (double) d / (double) (d + 1));
  mpz_root (out->m, out->n, d + 1);
  B *= mpz_get_d (out->m);
  mpz_set_d (out->m, B);
  while (T-- > 0)
    {
      generate_base_m (out, out->m);
      mu = get_mu (out->f, out->degree, out->m, out->skew);
      if (mu < best_mu)
	{
	  best_mu = mu;
	  if (verbose)
	    {
	      fprintf (stderr, "m=");
	      mpz_out_str (stderr, 10, out->m);
	      fprintf (stderr, " mu=%1.3e\n", mu);
	    }
	  mpz_set (best_m, out->m);
	}
      mpz_add_ui (out->m, out->m, 1);
    }
  generate_base_m (out, best_m);
  generate_sieving_parameters (out);
  mpz_clear (best_m);
}

void
print_poly (FILE *fp, cado_poly p)
{
  int i;
  double S;

  if (strlen (p->name) != 0)
    fprintf (fp, "name: %s\n", p->name);
  fprintf (fp, "n: ");
  mpz_out_str (fp, 10, p->n);
  fprintf (fp, "\n");
  fprintf (fp, "skew: %1.3f\n", S = p->skew);
  fprintf (fp, "# mu(m,S): %1.3e\n", get_mu (p->f, p->degree, p->m, S));
  for (i = p->degree; i >= 0; i--)
    {
      fprintf (fp, "c%d: ", i);
      mpz_out_str (fp, 10, p->f[i]);
      fprintf (fp, "\n");
    }
  fprintf (fp, "Y1: ");
  mpz_out_str (fp, 10, p->g[1]);
  fprintf (fp, "\n");
  fprintf (fp, "Y0: ");
  mpz_out_str (fp, 10, p->g[0]);
  fprintf (fp, "\n");
  fprintf (fp, "# m: ");
  mpz_out_str (fp, 10, p->m);
  fprintf (fp, "\n");
  fprintf (fp, "type: %s\n", p->type);
  fprintf (fp, "rlim: %lu\n", p->rlim);
  fprintf (fp, "alim: %lu\n", p->alim);
  fprintf (fp, "lpbr: %d\n", p->lpbr);
  fprintf (fp, "lpba: %d\n", p->lpba);
  fprintf (fp, "mfbr: %d\n", p->mfbr);
  fprintf (fp, "mfba: %d\n", p->mfba);
  fprintf (fp, "rlambda: %1.1f\n", p->rlambda);
  fprintf (fp, "alambda: %1.1f\n", p->alambda);
  fprintf (fp, "qintsize: %d\n", p->qintsize);
  fprintf (fp, "q0: %d\n", p->alim + 1);
}

/* return L_{1/3}(n,c) = exp(c log(n)^{1/3} log(log(n))^{2/3})
   with c = (64/9)^{1/3} ~ 1.922999
*/
double
L (mpz_t n)
{
  double logn, e;
  double c = 1.9229994270765444764;
  
  logn = log (mpz_get_d (n));
  e = log (logn);
  e = logn * (e * e);
  e = c * pow (e, 0.33333333333333333333);
  return exp (e);
}

int
main (int argc, char *argv[])
{
  cado_input in;
  cado_poly out;
  int verbose = 0;
  unsigned long effort = 1000000;

  /* parse options */
  while (argc > 1 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-v") == 0)
	{
	  verbose = 1;
	  argv ++;
	  argc --;
	}
      else if (argc > 2 && strcmp (argv[1], "-e") == 0)
	{
	  effort = atoi (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else
	{
	  fprintf (stderr, "Unknown option: %s\n", argv[1]);
	  usage ();
	}
    }

  if (argc != 1)
    usage ();

  init (in);
  parse_input (in);
  if (verbose)
    fprintf (stderr, "L[1/3,(64/9)^(1/3)](n)=%1.3e\n", L (in->n));
  init_poly (out, in);
  clear (in);

  /* generate_poly_0 (out); */
  generate_poly (out, effort, verbose);
  print_poly (stdout, out);
  clear_poly (out);

  return 0;
}
