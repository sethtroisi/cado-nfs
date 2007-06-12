#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>

#define FAST /* uses product tree instead of naive accumulation */

#define THRESHOLD 2 /* must be >= 2 */

/* accumulate up to THRESHOLD products in prd[0], 2^i*THRESHOLD in prd[i].
   nprd is the number of already accumulated values: if nprd = n0 + 
   n1 * THRESHOLD + n2 * THRESHOLD^2 + ..., then prd[0] has n0 entries,
   prd[1] has n1*THRESHOLD entries, and so on.
*/
mpz_t*
accumulate_fast (mpz_t *prd, mpz_t a, unsigned long *lprd, unsigned long nprd)
{
  unsigned long i;

  mpz_mul (prd[0], prd[0], a);
  nprd ++;

  for (i = 0; nprd % THRESHOLD == 0; i++, nprd /= THRESHOLD)
    {
      /* need to access prd[i + 1], thus i+2 entries */
      if (i + 2 > *lprd)
        {
          lprd[0] ++;
          prd = (mpz_t*) realloc (prd, *lprd * sizeof (mpz_t));
          mpz_init_set_ui (prd[i + 1], 1);
        }
      mpz_mul (prd[i + 1], prd[i + 1], prd[i]);
      mpz_set_ui (prd[i], 1);
    }

  return prd;
}

/* prd[0] <- prd[0] * prd[1] * ... * prd[lprd-1] */
void
accumulate_fast_end (mpz_t *prd, unsigned long lprd)
{
  unsigned long i;

  for (i = 1; i < lprd; i++)
    mpz_mul (prd[0], prd[0], prd[i]);
}

int main(int argc, char **argv) {
  FILE * matfile, *kerfile;
  mpz_t m, a, b, *prd;
  unsigned long lprd; /* number of elements in prd[] */
  unsigned long nprd; /* number of accumulated products in prd[] */
  int ret;
  unsigned long w;
  int i, j, nlimbs;
  char str[1024];
  int depnum = 0;

  if (argc > 2 && strcmp (argv[1], "-depnum") == 0)
    {
      depnum = atoi (argv[2]);
      argc -= 2;
      argv += 2;
    }

  if (argc != 4) {
    fprintf(stderr, "usage: %s [-depnum nnn] matfile kerfile m\n", argv[0]);
    exit(1);
  }

  matfile = fopen(argv[1], "r");
  assert (matfile != NULL);
  kerfile = fopen(argv[2], "r");
  assert (kerfile != NULL);

  mpz_init(m);
  ret = gmp_sscanf(argv[3], "%Zd", m);
  assert (ret == 1);

  lprd = 1;
  nprd = 0;
  prd = (mpz_t*) malloc (lprd * sizeof (mpz_t));
  mpz_init_set_ui(prd[0], 1);
  mpz_init(a);
  mpz_init(b);

  {
    int nrows, ncols;
    ret = fscanf(matfile, "%d %d", &nrows, &ncols);
    assert (ret == 2);
    fgets(str, 1024, matfile); // read end of first line
    nlimbs = (nrows / GMP_NUMB_BITS) + 1;
  }

  /* go to dependency depnum */
  while (depnum > 0)
    {
      int c;
      /* read one line */
      while ((c = fgetc (kerfile)) != '\n')
        if (c == EOF)
          break;
      depnum --;
    }

  if (depnum > 0)
    {
      fprintf (stderr, "Error, not enough dependencies\n");
      exit (1);
    }

  for (i = 0; i < nlimbs; ++i) {
    ret = fscanf(kerfile, "%lx", &w);
    assert (ret == 1);
    for (j = 0; j < GMP_NUMB_BITS; ++j) {
      if (fgets(str, 1024, matfile)) {
	if (w & 1UL) {
	  ret = gmp_sscanf(str, "%Zd %Zd", a, b);
	  assert (ret == 2);
	  mpz_mul(b, b, m);
	  mpz_sub(a, a, b);
#ifndef FAST
          mpz_mul (prd[0], prd[0], a);
#else
          prd = accumulate_fast (prd, a, &lprd, nprd++);
#endif
	}
      }
      w >>= 1;
    }
  }

#ifdef FAST
  accumulate_fast_end (prd, lprd);
#endif

  printf("size of prd = %d bits\n", mpz_sizeinbase(prd[0], 2));

  if (mpz_sgn (prd[0]) < 0)
    {
      fprintf (stderr, "Error, product is negative: try another dependency\n");
      exit (1);
    }

  mpz_sqrtrem(prd[0], a, prd[0]);
  gmp_printf("remainder = %Zd\n", a);

  mpz_mod(prd[0], prd[0], m);
  gmp_printf("rational square root is %Zd\n", prd[0]);

  for (i = 0; i < lprd; i++)
    mpz_clear (prd[i]);
  mpz_clear (a);
  mpz_clear (b);
  free (prd);

  return 0;
}
