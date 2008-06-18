/***************************** input/output **********************************/

#include <ctype.h>  /* for islower() */
#include <string.h> /* for strcmp() */

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

/************************ arrays of mpz_t ************************************/

/* allocate an array of d coefficients, and initialize it */
mpz_t*
alloc_mpz_array (int d)
{
  mpz_t *f;
  int i;

  f = (mpz_t*) malloc (d * sizeof (mpz_t));
  for (i = 0; i < d; i++)
    mpz_init (f[i]);
  return f;
}

/* reallocate an array having d0 coefficients to d > d0 coefficients */
mpz_t*
realloc_mpz_array (mpz_t *f, int d0, int d)
{
  int i;

  f = (mpz_t*) realloc (f, d * sizeof (mpz_t));
  for (i = d0; i < d; i++)
    mpz_init (f[i]);
  return f;
}

/* free an array of d coefficients */
void
clear_mpz_array (mpz_t *f, int d)
{
  int i;
  
  for (i = 0; i < d; i++)
    mpz_clear (f[i]);
  free (f);
}


