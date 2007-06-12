#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>

int main(int argc, char **argv) {
  FILE * matfile, *kerfile;
  mpz_t m, a, b, prd;
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

  mpz_init_set_ui(prd, 1);
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
	  mpz_mul(prd, prd, a);
	}
      }
      w >>= 1;
    }
  }

  printf("size of prd = %d bits\n", mpz_sizeinbase(prd, 2));

  if (mpz_sgn (prd) < 0)
    {
      fprintf (stderr, "Error, product is negative: try another dependency\n");
      exit (1);
    }

  mpz_sqrtrem(prd, a, prd);
  gmp_printf("remainder = %Zd\n", a);

  mpz_mod(prd, prd, m);
  gmp_printf("rational square root is %Zd\n", prd);


  return 0;
}


