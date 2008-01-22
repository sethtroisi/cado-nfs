#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h> /* for isdigit */
#include "cado.h"


#define TYPE_STRING 0
#define TYPE_MPZ 1
#define TYPE_INT 2
#define TYPE_ULONG 3
#define TYPE_DOUBLE 4
#define PARSE_MATCH -1
#define PARSE_ERROR 0
#define PARSE_NOMATCH 1

/* Return value: 1: tag does not match, 0: error occurred, -1 tag did match */

static int
parse_line (void *target, char *line, const char *tag, int *have, 
	    const int type)
{
  char *lineptr = line;

  if (strncmp (lineptr, tag, strlen (tag)) != 0)
    return PARSE_NOMATCH;

  if (have != NULL && *have != 0)
    {
      fprintf (stderr, "parse_line: %s appears twice\n", tag);
      return PARSE_ERROR;
    }
  if (have != NULL)
    *have = 1;
  
  lineptr += strlen (tag);
  if (type == TYPE_STRING) /* character string of length up to 256 */
    {
      strncpy ((char *)target, lineptr, 256);
      ((char *)target)[255] = '0';
    }
  else if (type == TYPE_MPZ)
    {
      mpz_init (*(mpz_t *)target);
      mpz_set_str (*(mpz_t *)target, lineptr, 0);
    }
  else if (type == TYPE_INT)
    *(int *)target = atoi (lineptr);
  else if (type == TYPE_ULONG)
    *(unsigned long *)target = strtoul (lineptr, NULL, 10);
  else if (type == TYPE_DOUBLE)
    *(double *)target = atof (lineptr);
  else return PARSE_ERROR;
  
  return PARSE_MATCH;
}

// return 0 on failure, 1 if success.
int  
read_polynomial (cado_poly poly, char *filename)
{
  FILE *file;
  const int linelen = 512;
  char line[linelen];
  int have_name = 0, have_n = 0, have_Y0 = 0, have_Y1 = 0;
  int i, ok = PARSE_ERROR;

  file = fopen (filename, "r");
  if (file == NULL)
    {
      fprintf (stderr, "read_polynomial: could not open %s\n", filename);
      return 0;
    }

  poly->f = (mpz_t *) malloc (MAXDEGREE * sizeof (mpz_t));
  poly->g = (mpz_t *) malloc (2 * sizeof (mpz_t));

  poly->name[0] = '\0';
  poly->degree = -1;
  poly->type[0] = '\0';

  while (!feof (file))
    {
      ok = 1;
      if (fgets (line, linelen, file) == NULL)
	break;
      if (line[0] == '#')
	continue;

      ok *= parse_line (&(poly->name), line, "name: ", &have_name, TYPE_STRING);
      ok *= parse_line (&(poly->n), line, "n: ", &have_n, TYPE_MPZ);
      ok *= parse_line (&(poly->skew), line, "skew: ", NULL, TYPE_DOUBLE);
      ok *= parse_line (&(poly->g[0]), line, "Y0: ", &have_Y0, TYPE_MPZ);
      ok *= parse_line (&(poly->g[1]), line, "Y1: ", &have_Y1, TYPE_MPZ);
      ok *= parse_line (&(poly->type), line, "type: ", NULL, TYPE_STRING);
      ok *= parse_line (&(poly->rlim), line, "rlim: ", NULL, TYPE_ULONG);
      ok *= parse_line (&(poly->alim), line, "alim: ", NULL, TYPE_ULONG);
      ok *= parse_line (&(poly->lpbr), line, "lpbr: ", NULL, TYPE_INT);
      ok *= parse_line (&(poly->lpba), line, "lpba: ", NULL, TYPE_INT);
      ok *= parse_line (&(poly->mfbr), line, "mfbr: ", NULL, TYPE_INT);
      ok *= parse_line (&(poly->mfba), line, "mfba: ", NULL, TYPE_INT);
      ok *= parse_line (&(poly->rlambda), line, "rlambda: ", NULL, TYPE_DOUBLE);
      ok *= parse_line (&(poly->alambda), line, "alambda: ", NULL, TYPE_DOUBLE);
      ok *= parse_line (&(poly->qintsize), line, "qintsize: ", NULL, TYPE_INT);


      if (ok == 1 && line[0] == 'c' && isdigit (line[1]) && line[2] == ':' &&
	       line[3] == ' ')
	{
	  int index = line[1] - '0', i;
	  for (i = poly->degree + 1; i <= index; i++)
	    mpz_init (poly->f[i]);
	  if (index > poly->degree)
	    poly->degree = index;
	  mpz_set_str (poly->f[index], line + 4, 0);
	  ok = -1;
	}

      if (ok == PARSE_NOMATCH)
	{
	  fprintf (stderr, 
		   "read_polynomial: Cannot parse line %s\nIgnoring.\n", 
		   line);
	  continue;
	}
      
      if (ok == PARSE_ERROR)
	break;
    }
  
  if (ok != PARSE_ERROR)
    {
      if (have_n == 0)
	{
	  fprintf (stderr, "n ");
	  ok = PARSE_ERROR;
	}
      if (have_Y0 == 0)
	{
	  fprintf (stderr, "Y0 ");
	  ok = PARSE_ERROR;
	}
      if (have_Y1 == 0)
	{
	  fprintf (stderr, "Y1 ");
	  ok = PARSE_ERROR;
	}
      if (ok == PARSE_ERROR)
	fprintf (stderr, "are missing in polynomial file\n");
    }

  if (ok == PARSE_ERROR)
    {
      for (i = 0; i <= poly->degree; i++)
	mpz_clear (poly->f[i]);
      if (have_n)
	mpz_clear (poly->n);
      if (have_Y0)
	mpz_clear (poly->g[0]);
      if (have_Y1)
	mpz_clear (poly->g[1]);
      free (poly->f);
      free (poly->g);
      return 0;
    }
  
  // compute m, the common root of f and g mod n
  mpz_init(poly->m);
  mpz_invert(poly->m, poly->g[1], poly->n);
  mpz_mul(poly->m, poly->m, poly->g[0]);
  mpz_mod(poly->m, poly->m, poly->n);
  mpz_sub(poly->m, poly->n, poly->m);
  mpz_mod(poly->m, poly->m, poly->n);
  return 1;
}

void
fprint_polynomial (FILE *fp, mpz_t *f, const int d)
{
  int i, s, first = 1;
  mpz_t c;

  mpz_init (c);
  for (i = d; i >= 0; i--)
    {
      s = mpz_cmp_ui (f[i], 0);
      if (s != 0)
        {
          if (s > 0)
            {
              if (first == 0)
                gmp_fprintf (fp, " + ");
              gmp_fprintf (fp, "%Zd", f[i]);
            }
          else if (s < 0)
            {
              mpz_abs (c, f[i]);
              gmp_fprintf (fp, " - %Zd", c);
            }
          first = 0;
          if (i >= 2)
            gmp_fprintf (fp, "*x^%d", i);
          else if (i == 1)
            gmp_fprintf (fp, "*x");
        }
    }
  mpz_clear (c);
  fprintf (fp, "\n");
}

#undef TYPE_STRING 
#undef TYPE_MPZ 
#undef TYPE_INT 
#undef TYPE_ULONG 
#undef TYPE_DOUBLE 
#undef PARSE_MATCH 
#undef PARSE_ERROR 
#undef PARSE_NOMATCH 


