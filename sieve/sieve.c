/* Apply the sieve updates to sievearray */
/* Output: sievearray. Must have enough memory allocated */
/* Inputs: amin, amax, b the $a$ range and the line $b$ to sieve */
/* factorbase, polynomial */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#define MAXDEGREE 10
#ifdef WANT_ASSERT
  #include <assert.h>
  #define ASSERT(x) assert(x)
#else
  #define ASSERT(x)
#endif


typedef uint32_t fbprime_t;
typedef fbprime_t fbroot_t;

typedef struct {
  fbprime_t p;            /* A prime or a prime power */
  unsigned char plog;     /* logarithm (to some suitable base) of this prime */
  unsigned char nr_roots; /* how many roots there are for this prime */
  unsigned char dummy[2]; /* For dword aligning the roots */
  fbroot_t roots[0];      /* the actual length of this array is determined
                             by nr_roots */
} factorbase32_t;


/* Hack to get around C's automatic multiplying constants to be added to 
   pointers by the pointer's base data type's size */
static inline factorbase32_t *
fb_skip (factorbase32_t *fb, size_t s)
{
  return (factorbase32_t *)((void *)fb + s);
}

static size_t
fb_entrysize (const factorbase32_t *fb)
{
  return (sizeof (factorbase32_t) + fb->nr_roots * sizeof (fbroot_t));
}


unsigned long
mulmod (const unsigned long a, const unsigned long b, const unsigned long m)
{
  unsigned long _r, _a = a;

  __asm__ ( "mulq %2\n\t"
            "divq %3"
            : "=&d" (_r), "+a" (_a)
            : "rm" (b), "rm" (m)
            : "cc");

  return _r;
}


void 
sieve (unsigned char *sievearray, factorbase32_t *fb, 
       long amin, long amax, unsigned long b)
{
  uint32_t i, amin_p, a, p, d;
  unsigned char plog;
  const unsigned long l = amax - amin;

  ASSERT (amax > amin);

  /* Init the array. */
  
  memset (sievearray, 0, l);

  /* Do the sieving */

  while (fb->p > 0)
    {
      p = fb->p;
      plog = fb->plog;

      /* This modular reduction should be simplified somehow. Do it once
         and store it in fb? */
      if (amin < 0)
        amin_p = p - ((unsigned long)(-amin) % p);
      else
        amin_p = (unsigned long) amin % p;

      for (i = 0; i < fb->nr_roots; i++)
        {
          /* Find first index in sievearray where p on this root divides.
             So we want a/b == r (mod p) <=> a == br (mod p). Then we want
             d so that amin + d == a (mod p) <=> d == a - amin (mod p)
          */

          a = mulmod (b, fb->roots[i], p); /* We have r_i < p, hence b*r_i/p 
	  its in register and there is no division overflow. If we keep r_i
          in Montgomery representation, a single mul/REDC will compute the
          product, reduce it mod p and return it as an integer. TBD. */

          d = (a >= amin_p) ? a - amin_p : a + p - amin_p; 
	  /* Now d is the first index into sievearray where p divides. */

          /* Now update the sieve array. There is no partitioning, blocking, 
             bucket sieving or anything atm. */
          for (; d < l; d += p)
            sievearray[d] += plog;
        }
      
      /* Move on to the next factor base prime */
      fb = fb_skip (fb, fb_entrysize (fb));
    }
}


void 
print_fb_entry (factorbase32_t *fb)
{
  int i;
  printf ("%lu: ", (unsigned long) fb->p);
  for (i = 0; i + 1 < fb->nr_roots; i++)
    printf ("%lu,", (unsigned long) fb->roots[i]);
  printf ("%lu\n", (unsigned long) fb->roots[i]);
}


/* Add fb_add to (void *)fb + fb_size. If a realloc failed, returns NULL. */

void *
add_to_fb (factorbase32_t *fb, size_t *fbsize, size_t *fballoc,
	   const size_t allocblocksize, factorbase32_t *fb_add)
{
  const size_t fb_addsize  = fb_entrysize (fb_add); 
  factorbase32_t *newfb = fb;

  ASSERT(fb_addsize <= allocblocksize); /* Otherwise we still might not have
					   enough mem after the realloc */

  /* Do we need more memory for fb? */
  if (*fballoc < *fbsize + fb_addsize)
    {
      *fballoc += allocblocksize;
      newfb = realloc (fb, *fballoc);
      if (newfb == NULL)
	{
	  fprintf (stderr, 
		   "Could not reallocate factor base to %lu bytes\n",
		   *fballoc);
	  return NULL;
	}
    }
  memcpy ((void *)newfb + *fbsize, fb_add, fb_addsize);
  *fbsize += fb_addsize;

  return newfb;
}


factorbase32_t *
read_fb (const char *filename, const double log_scale)
{
  factorbase32_t *fb = NULL, *fb_cur, *fb_new;
  FILE *fbfile;
  size_t fbsize = 0, fballoc = 0;
  const size_t linesize = 255;
  size_t linelen;
  char line[linesize];
  char *lineptr;
  unsigned int i;
  int ok;
  unsigned long linenr = 0;
  const size_t allocblocksize = 1<<20; /* Allocate in MB chunks */

  fbfile = fopen (filename, "r");
  if (fbfile == NULL)
    {
      fprintf (stderr, "Could not open file %s for reading\n", filename);
      return NULL;
    }

  fb_cur = malloc (sizeof (factorbase32_t) + MAXDEGREE * sizeof(fbroot_t));
  if (fb_cur == NULL)
    {
      fclose (fbfile);
      return NULL;
    }

  while (!feof(fbfile))
    {
      if (fgets (line, linesize, fbfile) == NULL)
	break;
      linenr++;
      linelen = strlen (line);
      if (linelen > 0 && line[linelen - 1] == '\n')
	linelen--; /* Remove newline */
      for (i = 0; i < linelen; i++) /* Skip comments */
	if (line[i] == '#')
	  break;
      linelen = i;
      if (linelen == 0) /* Skip empty/comment lines */
	continue;
      line[linelen] = '\0';

      /* Parse the line */
      lineptr = line;
      fb_cur->p = strtoul (lineptr, &lineptr, 10);
      ok = 0;
      if (fb_cur->p == 0)
	fprintf (stderr, "read_fb: prime is not an integer on line %lu\n", 
		 linenr);
      else if (*lineptr != ':')
	fprintf (stderr, "read_fb: prime is not followed by colon on line %lu",
		 linenr);
      else
	ok = 1;
      if (!ok)
	continue;

      /* FIXME: If this is a prime power, compute log(p^k)-log(p^(k-1)) */
      fb_cur->plog = round (log ((double) fb_cur->p) * log_scale);

      lineptr++; /* Skip colon */
      ok = 1;
      fb_cur->nr_roots = 0;
      /* Read roots */
      while (ok && *lineptr != '\0' && fb_cur->nr_roots < MAXDEGREE)
	{
	  fb_cur->roots[fb_cur->nr_roots++] = strtoul (lineptr, &lineptr, 10);
	  if (*lineptr != '\0' && *lineptr != ',')
	    ok = 0;
	  if (*lineptr == ',')
	    lineptr++;
	}
      
      fb_new = add_to_fb (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
      fb = fb_new;
    }      

  fb_cur->p = 0;
  fb->nr_roots = 0;
  
  fb_new = add_to_fb (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
  if (fb_new == NULL)
    free (fb);
  fb = fb_new;
  
  fclose (fbfile);
  free (fb_cur);
  return fb;
}

int
main (int argc, char **argv)
{
  factorbase32_t *fb;
  size_t s;
  long amin, amax, i;
  unsigned long b;
  unsigned char *sievearray;
  const double log_scale = 1. / log (2.); /* Lets use log_2() for a start */
  int threshold = 10;

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 2 && strcmp (argv[1], "-thres") == 0)
	{
	  threshold = atoi (argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else 
	break;
    }

  if (argc < 2)
    exit (EXIT_FAILURE);

  fb = read_fb (argv[1], log_scale);

  if (fb == NULL)
    {
      fprintf (stderr, "Could not read factor base\n");
      exit (EXIT_FAILURE);
    }

  for (s = 0; fb_skip(fb, s)->p != 0; s += fb_entrysize (fb_skip (fb, s)))
    print_fb_entry (fb_skip(fb, s));

  amin = -1000;
  amax = 1000;
  b = 100;

  sievearray = malloc ((amax - amin + 1) * sizeof (char));
  
  sieve (sievearray, fb, amin, amax, b);

  for (i = amin; i <= amax; i++)
    if ((int) sievearray[i - amin] > threshold)
      printf ("%ld: %d  ", i, (int) sievearray[i - amin]);



  free (fb);
  free (sievearray);

  return 0;
}
