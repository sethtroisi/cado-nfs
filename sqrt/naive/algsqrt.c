#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>
#include "cado.h"
#include "utils/utils.h"
#include "poly.h"


void polymodF_from_ab(polymodF_t tmp, long a, unsigned long b) {
  tmp->v = 0;
  tmp->p->deg = 1;
  mpz_set_ui(tmp->p->coeff[1], b);
  mpz_neg(tmp->p->coeff[1], tmp->p->coeff[1]);
  mpz_set_si(tmp->p->coeff[0], a);
}
#if 0
void polymodF_reduce_mod_ui(poly_t z, polymodF_t x, unsigned long p) {

}

// z := x^n mod f (mod p)
void poly_power_mod_ui(poly_t z, poly_t x, poly_t f, mpz_t n, unsigned long p) {

}

void poly_reduce_makemonic_mod(poly_t f, poly_t F, mpz_t pk) {

}

void poly_reduce_mod(poly_t f, poly_t F, mpz_t pk) {

}

void polymodF_sqrt(polymodF_t res, polymodF_t x, poly_t F, unsigned long p) {
  poly_t tmp, tmp2, f;
  const int d=F->deg;
  mpz_t q, pk, aux;

  mpz_init(q);
  mpz_init(pk);
  mpz_init(aux);
  poly_alloc(tmp, d);
  poly_alloc(tmp2, d);
  poly_alloc(f, d);
  
  mpz_ui_pow_ui(q, p, (unsigned long)d);
  mpz_set_ui(pk, p);

  poly_reduce_makemonic_mod(f, F, pk);

  polymodF_reduce_mod_ui(tmp, x, p);	// tmp := x mod p
  mpz_add_ui(aux, q, 1);
  mpz_sub_ui(q, q, 2);
  mpz_mul(aux, aux, q);
  mpz_add_ui(q, q, 2);
  mpz_divexact_ui(aux, aux, 4);		// aux := (q-2)(q+1)/4
  poly_power_mod_ui(tmp, tmp, f, aux, p);	// tmp := tmp^((q-2)(q+1)/4) mod p

  // now tmp contains the inverse of the square root of x modulo p
  // Let's lift it.

  do {
    mpz_sqr(pk, pk);	// double the current precision
    poly_reduce_makemonic_mod(f, F, pk);  

    poly_mul(tmp2, tmp, tmp);
    poly_mod(tmp2, tmp2, f);
    poly_reduce_mod(tmp2, tmp2, pk);



    







  mpz_clear(q);
  mpz_clear(qk);
  mpz_clear(aux);
  poly_free(tmp);
  poly_free(tmp2);
  poly_free(f);
}
#endif


int main(int argc, char **argv) {
  FILE * matfile, *kerfile;
  cado_poly pol;
  long a;
  unsigned long b;
  poly_t F;
  polymodF_t prd, tmp;
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
    fprintf(stderr, "usage: %s [-depnum nnn] matfile kerfile polyfile\n", argv[0]);
    exit(1);
  }

  matfile = fopen(argv[1], "r");
  assert (matfile != NULL);
  kerfile = fopen(argv[2], "r");
  assert (kerfile != NULL);

  ret = read_polynomial(pol, argv[3]);
  assert (ret == 1);

  // Init F to be the algebraic polynomial
  poly_alloc(F, pol->degree);
  for (i = pol->degree; i >= 0; --i)
    poly_setcoeff(F, i, pol->f[i]);

  // Init prd to 1.
  poly_alloc(prd->p, pol->degree);
  mpz_set_ui(prd->p->coeff[0], 1);
  prd->p->deg = 0;
  prd->v = 0;

  // Allocate tmp
  poly_alloc(tmp->p, 1);

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
	  ret = sscanf(str, "%ld %lu", &a, &b);
	  assert (ret == 2);

	  polymodF_from_ab(tmp, a, b);
	  polymodF_mul(prd, prd, tmp, F);
	}
      }
      w >>= 1;
    }
  }

  fprintf(stderr, "v = %d\n", prd->v);
  fprintf(stderr, "sizeinbit of cst term = %d\n", mpz_sizeinbase(prd->p->coeff[0], 2));
  fprintf(stderr, "sizeinbit of leading term = %d\n", mpz_sizeinbase(prd->p->coeff[pol->degree-1], 2));

  for (i = 0; i < pol->degree; ++i) 
    gmp_printf("%Zx\n", prd->p->coeff[i]);
  printf("%d\n", prd->v);

  {
    FILE * formagma;
    mpz_t m;

    formagma = fopen("/tmp/formagma", "w");
    if (formagma == NULL) {
      fprintf(stderr, "can not create /tmp/formagma\n");
      return 1;
    }
    fprintf(formagma, "f := ");
    for (i = 0; i <= pol->degree; ++i)
      gmp_fprintf(formagma, "%+Zd*x^%d", pol->f[i], i);
    fprintf(formagma, ";\nm := ");
    mpz_init(m);
    mpz_neg(m, pol->g[0]);
    gmp_fprintf(formagma, "%Zd;\nN := %Zd;\n", m, pol->n);
    mpz_clear(m);
    fclose(formagma);
  }

  poly_free(F);
  poly_free(prd->p);
  poly_free(tmp->p);

  return 0;
}
