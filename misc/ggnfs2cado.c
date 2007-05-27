/* convert relations from GGNFS rels.bin.* files to CADO or Franke/Kleinjung's
   formats */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tiff.h>
#include "gmp.h"

#define MAX_PRIMES 255 /* maximal number of factor base primes */
#define MAX_LPRIMES 3  /* maximal number of large primes */
#define DEGF_MAX 5
#define REPS 10

#define OUTPUT_CADO 0
#define OUTPUT_FK   1 /* Franke-Kleinjung's format */

uint32
get_uint32 (FILE *fp)
{
  uint32 w;
  fread (&w, sizeof (uint32), 1, fp);
  return w;
}

int32
get_int32 (FILE *fp)
{
  int32 w;
  fread (&w, sizeof (int32), 1, fp);
  return w;
}

typedef struct
{
  int32 a, b;
  int rfb_entries;
  int afb_entries;
  int sp_entries;
  int num_lrp;
  int num_lap;
  unsigned long rprimes[MAX_PRIMES]; /* index i of prime p(i) */
  unsigned long rexp[MAX_PRIMES];
  unsigned long aprimes[MAX_PRIMES];
  unsigned long aexp[MAX_PRIMES];
  unsigned long sprimes[MAX_PRIMES];
  int32 sexp[MAX_PRIMES];
  uint32 large_rprimes[MAX_LPRIMES]; /* large rational primes are 32-bit */
  uint32 large_aprimes[MAX_LPRIMES]; /* large algebraic primes are 32-bit
					  too, but the file stores also the
					  corresponding root, which we don't
					  need here. */
  int32 qcb1, qcb2; /* quadratic characters */
} relation_t;

/* prints a large prime, given its low and high 32-bit parts */
void
print_large_prime (uint32 l, uint32 h)
{
  mpz_t P;

  mpz_init (P);
  mpz_set_ui (P, h);
  mpz_mul_2exp (P, P, 32);
  mpz_add_ui (P, P, l);
  mpz_out_str (stdout, 16, P);
  mpz_clear (P);
}

/* output a relation in CADO format:
   a,b:p1,p2,...,pn:q1,q2,...,qn
   where a and p are printed in decimal, p1,...,qn are printed in hexadecimal,
   p1,p2,...,pn are the primes on the rational side.
*/
void
print_relation_cado (relation_t *rel, int32 *rfb, int32 *afb)
{
  int i,j;

  printf ("%d,%d", rel->a, rel->b);
  /* rational side */
  for (i = 0; i < rel->rfb_entries; i++)
    for (j = 0; j < rel->rexp[i]; j++)
      {
	putchar ((i + j == 0) ? ':' : ',');
	printf ("%x", rfb[rel->rprimes[i]]);
      }
  for (i = 0; i < rel->num_lrp; i++)
    {
      putchar ((rel->rfb_entries + i == 0) ? ':' : ',');
      printf ("%x", rel->large_rprimes[i]);
    }
  /* algebraic side */
  for (i = 0; i < rel->afb_entries; i++)
    for (j = 0; j < rel->aexp[i]; j++)
      {
	putchar ((i + j == 0) ? ':' : ',');
	printf ("%x", afb[rel->aprimes[i]]);
      }
  for (i = 0; i < rel->sp_entries; i++)
    {
      if (rel->sexp[i] > 0)
	for (j = 0; j < rel->sexp[i]; j++)
	  {
	    putchar ((rel->afb_entries + i + j == 0) ? ':' : ',');
	    printf ("%x", rel->sprimes[i]);
	  }
    }
  for (i = 0; i < rel->num_lap; i++)
    {
      putchar ((rel->afb_entries + rel->sp_entries == 0) ? ':' : ',');
      printf ("%x", rel->large_aprimes[i]);
    }
  printf ("\n");
}

/* output a relation in Franke/Kleinjung format:
   - the first line is "F 0 X 5 1" where the two last integers are the degrees
   - then it comes by groups of 3 consecutive lines (W gives a and b in
     hexadecimal form, X gives the algebraic factors, Y the linear factors.
   W 2335e65b 1d9
   X 96f8323 1413403 1369 16B5 1A0C77 347 55D 3 5 65 71 9D F42439
   Y 444B A6CF C827 2C8627 B 25 209 527 B 2 2 2
*/
void
print_relation_fk (relation_t *rel, int32 *rfb, int32 *afb)
{
  int i,j;

  if (rel->a > 0)
    printf ("W %x %x\n", rel->a, rel->b);
  else
    printf ("W -%x %x\n", -rel->a, rel->b);
  /* algebraic side */
  printf ("X");
  for (i = 0; i < rel->afb_entries; i++)
    for (j = 0; j < rel->aexp[i]; j++)
      printf (" %x", afb[rel->aprimes[i]]);
  for (i = 0; i < rel->sp_entries; i++)
    for (j = 0; j < rel->sexp[i]; j++)
      printf (" %x", rel->sprimes[i]);
  for (i = 0; i < rel->num_lap; i++)
    printf (" %x", rel->large_aprimes[i]);
  printf ("\n");
  /* rational side */
  printf ("Y");
  for (i = 0; i < rel->rfb_entries; i++)
    for (j = 0; j < rel->rexp[i]; j++)
      printf (" %x", rfb[rel->rprimes[i]]);
  for (i = 0; i < rel->num_lrp; i++)
    printf (" %x", rel->large_rprimes[i]);
  printf ("\n");
}

/* norm <- |f(a,b)| = |f[d]*a^d + ... + f[0]*b^d| */
void
eval_algebraic (mpz_t norm, mpz_t *f, int d, long a, unsigned long b)
{
  mpz_t B; /* powers of b */
  mpz_init_set_ui (B, 1);
  mpz_set (norm, f[d]);
  while (d-- > 0)
    {
      mpz_mul_si (norm, norm, a);
      mpz_mul_ui (B, B, b);
      mpz_addmul (norm, f[d], B);
    }
  mpz_clear (B);
  mpz_abs (norm, norm);
}

/* read relations from file fp, with rational factor base into rfb,
   algebraic factor base into afb, and coefficients from algebraic
   polynomial into f[0]...f[degf].
   Then prints relations with given format.
*/
unsigned long
read_relations (FILE *fp, int32 *rfb, int32 *afb, mpz_t *f, int degf,
		int format, int verbose)
{
  uint32 size_field; /* 32 bits */
  relation_t rel[1];
  int32 relsInFile;
  int i, j, k;
  mpz_t norm;
  int32 p;

  mpz_init (norm);

  /* the first 4 bytes give the number of relations in the file */
  rewind (fp);
  fread (&relsInFile, sizeof (int32), 1, fp);
  if (verbose)
    fprintf (stderr, "Header says %d relations\n", relsInFile);

  /* see function dataConvertToRel() from GGNFS, file rels.c */
  for (j = 0; j < relsInFile; j++)
    {
      size_field = get_uint32 (fp);
      rel->a = (int32) get_int32 (fp);
      rel->b = (int32) get_uint32 (fp);
      eval_algebraic (norm, f, degf, rel->a, rel->b);
      rel->rfb_entries = size_field >> 24;
      for (i = 0; i < rel->rfb_entries; i++)
	{
	  rel->rprimes[i] = get_uint32 (fp);
	  rel->rexp[i] = get_uint32 (fp);
	}
      rel->afb_entries = (size_field >> 16) & 255;
      for (i = 0; i < rel->afb_entries; i++)
	{
	  rel->aprimes[i] = get_uint32 (fp);
	  p = afb[rel->aprimes[i]];
	  rel->aexp[i] = get_uint32 (fp);
	  for (k = 0; k < rel->aexp[i]; k++)
	    {
	      if (!mpz_divisible_ui_p (norm, p))
		{
		  fprintf (stderr, "Error, f(%d,%d) not divisible by factor base prime %d\n",
			   rel->a, rel->b, p);
		  exit (1);
		}
	      mpz_divexact_ui (norm, norm, p);
	    }
	}
      rel->sp_entries = (size_field >> 8) & 255;
      for (i = 0; i < rel->sp_entries; i++)
	{
	  /* since GGNFS format gives the *index* of the special primes
	     and not their value, we have no easy mean to recover them.
	     The workaround is to factor the norm residue after all factor
	     base and large primes have been divided out. We thus put the
	     exponents to 0 for now. */
	  rel->sprimes[i] = get_uint32 (fp);
	  rel->sexp[i] = get_int32 (fp);
	  rel->sexp[i] = 0;
	}
      rel->qcb1 = get_uint32 (fp);
      rel->qcb2 = get_uint32 (fp);
      rel->num_lrp = (size_field >> 6) & 3;
      for (i = 0; i < rel->num_lrp; i++)
	rel->large_rprimes[i] = get_uint32 (fp);
      rel->num_lap = (size_field >> 4) & 3;
      for (i = 0; i < rel->num_lap; i++)
	{
	  rel->large_aprimes[i] = get_uint32 (fp);
	  get_uint32 (fp); /* corresponding root, not used here */
	  p = rel->large_aprimes[i];
	  if (!mpz_divisible_ui_p (norm, p))
	    {
	      fprintf (stderr, "Error, f(%d,%d) not divisible by large prime %d\n",
			   rel->a, rel->b, p);
	      exit (1);
	    }
	  mpz_divexact_ui (norm, norm, p);
	}

      i = 0; /* number of special primes */
      k = 2; /* smallest prime factor to try */
      while (mpz_cmp_ui (norm, 1) != 0)
	{
	  if (mpz_probab_prime_p (norm, REPS))
	    { /* assume it is a special prime */
	      if (!mpz_fits_ulong_p (norm))
		{
		  gmp_fprintf (stderr, "Error, special prime of f(%d,%d) does not fit in an unsigned long: %Zd\n", rel->a, rel->b, norm);
		  exit (1);
		}
	      rel->sprimes[i] = mpz_get_ui (norm);
	      rel->sexp[i] = 1;
	      i ++;
	      mpz_set_ui (norm, 1);
	    }
	  else /* search for a prime factor of norm */
	    {
	      while (mpz_divisible_ui_p (norm, k) == 0)
		k = k + 2 - (k == 2);
	      rel->sprimes[i] = k;
	      rel->sexp[i] = 0;
	      do
		{
		  rel->sexp[i] ++;
		  mpz_divexact_ui (norm, norm, k);
		}
	      while (mpz_divisible_ui_p (norm, k));
	      i ++;
	    }
	}

      switch (format)
	{
	case OUTPUT_CADO:
	  print_relation_cado (rel, rfb, afb);
	  break;
	case OUTPUT_FK:
	  print_relation_fk (rel, rfb, afb);
          break;
	default:
	  fprintf (stderr, "Error, unknown format %d\n", format);
	  exit (1);
	}
      fflush (stdout);
    }

  rel->sp_entries = i;

  mpz_clear (norm);

  return relsInFile;
}

#define BF_DELIMITER 0x0a
#define END_POLY "END_POLY"
#define END_HEADER "END_HEADER"

/* read fp until the string s is found,
   the next char read by getc (fp) is the one following s. */
void
get_string (char *s, FILE *fp, int verbose)
{
  int l, i, c;

  l = strlen (s);
  i = 0;
  while (i < l)
    {
      c = getc (fp);
      if (verbose)
	putc (c, stderr);
      if (c == s[i])
	i ++;
      else
	i = 0;
    }
}

/* read factor base from file fp, puts rational factor base into *rfb,
   (size into *rfb_size), algebraic factor base into *afb (size into
   *afb_size), and put coefficients from algebraic polynomial in f[0]..f[d].
   Returns degree d.
*/
int
read_fb (FILE *fp, int32 **rfb, int32 *rfb_size, int32 **afb, int32 *afb_size,
	 int verbose, mpz_t *f)
{
  int c, i, degf;
  long RFBsize, AFBsize;

  /* read polynomial */
  get_string ("Y1: ", fp, verbose);
  mpz_inp_str (f[0], fp, 10);
  if (verbose)
    mpz_out_str (stderr, 10, f[0]);
  c = getc (fp); /* newline */
  fscanf (fp, "c%d: ", &degf);
  if (degf > DEGF_MAX)
    {
      fprintf (stderr, "Error, too large degree\n");
      exit (1);
    }
  if (verbose)
    fprintf (stderr, "%cc%d: ", c, degf);
  mpz_inp_str (f[degf], fp, 10);
  if (verbose)
    mpz_out_str (stderr, 10, f[degf]);
  for (i = degf - 1; i >= 0; i--)
    {
      getc (fp); /* newline */
      fscanf (fp, "c%d: ", &c);
      if (c != i)
	{
	  fprintf (stderr, "Error, missing coefficient of degree %d\n", i);
	  exit (1);
	}
      if (verbose)
	fprintf (stderr, "\nc%d: ", i);
      mpz_inp_str (f[i], fp, 10);
      if (verbose)
	mpz_out_str (stderr, 10, f[i]);
    }
  get_string (END_POLY, fp, verbose);
  c = getc (fp);
  if (c != BF_DELIMITER)
    {
      fprintf (stderr, "Error, delimited expected\n");
      exit (1);
    }
  if (verbose)
    putc ('\n', stderr);

  /* read header */
  get_string ("RFBsize: ", fp, verbose);
  fscanf (fp, "%ld", &RFBsize);
  if (verbose)
    fprintf (stderr, "%ld", RFBsize);
  get_string ("AFBsize: ", fp, verbose);
  fscanf (fp, "%ld", &AFBsize);
  if (verbose)
    fprintf (stderr, "%ld", AFBsize);
  get_string (END_HEADER, fp, verbose);
  c = getc (fp);
  if (c != BF_DELIMITER)
    {
      fprintf (stderr, "Error, delimited expected\n");
      exit (1);
    }
  if (verbose)
    fprintf (stderr, "\nRFBsize=%ld AFBsize=%ld\n", RFBsize, AFBsize);
  *rfb_size = RFBsize;
  *rfb = (int32*) malloc (RFBsize * sizeof(int32));
  *afb_size = AFBsize;
  *afb = (int32*) malloc (AFBsize * sizeof(int32));
  for (i = 0; i < RFBsize; i++)
    {
      (*rfb)[i] = get_uint32 (fp);
      get_uint32 (fp); /* don't need the root mod g */
    }
  for (i = 0; i < AFBsize; i++)
    {
      (*afb)[i] = get_uint32 (fp);
      get_uint32 (fp); /* don't need the root mod f */
    }
  return degf;
}

int
main (int argc, char *argv[])
{
  FILE *fp;
  unsigned long num_rels;
  int32 *rfb = NULL, *afb = NULL;
  int32 rfb_size, afb_size;
  int verbose = 0, i, degf;
  int format = OUTPUT_CADO; /* default output format */
  mpz_t f[DEGF_MAX + 1];
  char *FB = NULL; /* name of factor base file */
  char *RELS; /* name of binary file to read */
  int num_files = 0;

  while (argc > 1 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-fk") == 0)
	{
	  format = OUTPUT_FK;
	  argv ++;
	  argc --;
	}
      else if (strcmp (argv[1], "-v") == 0)
	{
	  verbose = 1;
	  argv ++;
	  argc --;
	}
      else if (argc > 2 && strcmp (argv[1], "-fb") == 0)
	{
	  FB = argv[2];
	  argv += 2;
	  argc -= 2;
	}
      else
	{
	  fprintf (stderr, "Unknown option: %s\n", argv[1]);
	  exit (1);
	}
    }

  if (FB == NULL || argc <= 1)
    {
      fprintf (stderr, "Usage: %s [-fk] [-v] -fb <file> rels.bin.0 rels.bin.1 ...\n", argv[0]);
      exit (1);
    }

  for (i = 0; i <= DEGF_MAX; i++)
    mpz_init (f[i]);

  fp = fopen (FB, "rb");
  if (fp == NULL)
    {
      fprintf (stderr, "Error, unable to open factor base file %s\n", FB);
      exit (1);
    }
  degf = read_fb (fp, &rfb, &rfb_size, &afb, &afb_size, verbose, f);
  fclose (fp);
  
  while (argc > 1)
    {
      RELS = argv[1];
      argv ++;
      argc --;
      fp = fopen (RELS, "rb");
      if (fp == NULL)
	{
	  fprintf (stderr, "Error, unable to open relation file %s\n", RELS);
	  exit (1);
	}
      if (format == OUTPUT_FK && num_files == 0)
	printf ("F 0 X %d 1\n", degf);
      num_files ++;
      num_rels = read_relations (fp, rfb, afb, f, degf, format, verbose);
      fclose (fp);
    }

  fprintf (stderr, "Output %lu relations from %d file(s)\n", num_rels,
	   num_files);

  free (rfb);
  free (afb);
  for (i = 0; i <= DEGF_MAX; i++)
    mpz_clear (f[i]);

  return 0;
}
