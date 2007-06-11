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

  if (argc != 4) {
    fprintf(stderr, "usage: %s matfile kerfile m\n", argv[0]);
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

  mpz_sqrtrem(prd, a, prd);
  gmp_printf("remainder = %Zd\n", a);

  mpz_mod(prd, prd, m);
  gmp_printf("rational square root is %Zd\n", prd);


  return 0;
}


