#include "cado.h"
#include "utils.h"

#define _GNU_SOURCE // for pthread_tryjoin_np()
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h> 
#include <sys/resource.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <gmp.h>

#include "smooth_detect.h"

/*
ECM=...
gcc -Wall -g -std=c99 -I${ECM}/include descent_init_Fp.c smooth_detect.c -o descent_init_Fp -L${ECM}/lib -lecm -lgmp -lm -lpthread
*/


void HalfGcd(mpz_t a, mpz_t b, mpz_t u) {
  mpz_t v, w, x, q, r;
  mpz_set_ui(u, 1);
  mpz_init_set_ui(w, 0);
  mpz_init_set_ui(v, 0);
  mpz_init_set_ui(x, 1);
  mpz_init(q);
  mpz_init(r);
  /* invariant: a = u*a0 + v*b0 */
  while (mpz_cmpabs(a, u) > 0) {
    mpz_tdiv_qr(q, r, a, b);
    mpz_swap(a, b);
    mpz_swap(b, r);
    mpz_submul(u, q, w);
    mpz_swap(u, w);
    mpz_submul(v, q, x);
    mpz_swap(v, x);
  }
  mpz_clear(v);
  mpz_clear(w);
  mpz_clear(x);
  mpz_clear(q);
  mpz_clear(r);
}

// Function and parameters to be passed to smooth_detect.
typedef struct {
  mpz_t p;
  mpz_t z;
  mpz_poly f;   // used only in JL mode
  mpz_t m;      // root of f mod p (the one common with g)
} Fp_param;

void next_cand_Fp_hgcd(cand_t cand, void *params) {
  Fp_param * param = (Fp_param *)params;
  unsigned long e = random();
  mpz_t u0, v0;
  mpz_init(u0);
  mpz_init(v0);
  mpz_powm_ui(u0, param->z, e, param->p);
  mpz_t tmp;
  mpz_init(tmp);
  mpz_set(tmp, param->p);
  HalfGcd(u0, tmp, v0);
  mpz_abs(u0, u0);
  mpz_abs(v0, v0);
  cand_set_original_values(cand, u0, v0, e);
  mpz_clear(tmp);
  mpz_clear(u0);
  mpz_clear(v0);
}

// JL version
// Returns a boolean meaning "failure, try again".
int get_JL_candiate_from_e(unsigned long e, mpz_t u, mpz_t v, Fp_param * param,
    int print)
{
  mpz_t ze;
  mpz_init(ze);
  mpz_powm_ui(ze, param->z, e, param->p);

  //**** Create matrix (see Barbulescu's PhD thesis, section 8.4.2)
  // Warning: the LLL code count indices starting with 1 and not 0.
  // (this is a bug, imho)
  int d = param->f->deg;
  mat_Z M;
  LLL_init(&M, 2*d, 2*d);  // allocate and set to 0.
  // Topleft d*d block
  mpz_set(M.coeff[1][1], param->p);
  for (int i = 1; i < d; ++i) {
    mpz_set_ui(M.coeff[i+1][i+1], 1);
    mpz_neg(M.coeff[i+1][i], param->m);
  }
  // Bottomleft d*d block
  for (int i = 0; i < d; ++i) {
    mpz_set(M.coeff[d+i+1][i+1], ze);
  }
  // Bottomright d*d block
  for (int i = 0; i < d; ++i) {
    mpz_set_ui(M.coeff[d+i+1][d+i+1], 1);
  }

  //**** Apply LLL.
  { 
    mpz_t det, a, b;
    mpz_init(det);
    mpz_init_set_ui (a, 1);
    mpz_init_set_ui (b, 1);
    LLL(det, M, NULL, a, b);
    mpz_clear(det);
    mpz_clear(a);
    mpz_clear(b);
  }

  //**** Recover rational reconstruction
  // z^e = U(alpha)/V(alpha) mod (p, x-m)
  mpz_poly U, V;
  mpz_poly_init(U, d-1);
  mpz_poly_init(V, d-1);
  for (int i = 0; i < d; ++i) {
    mpz_poly_setcoeff(U, i, M.coeff[1][i+1]);
    mpz_poly_setcoeff(V, i, M.coeff[1][d+i+1]);
  }
  mpz_poly_cleandeg(U, d-1);
  mpz_poly_cleandeg(V, d-1);

  //**** Compute norms
  mpz_poly_resultant(u, U, param->f);
  mpz_poly_resultant(v, V, param->f);
  mpz_abs(u, u);
  mpz_abs(v, v);
  if (print) {
    printf("U = ");
    mpz_poly_fprintf_coeffs(stdout, U, ',');
    printf("V = ");
    mpz_poly_fprintf_coeffs(stdout, V, ',');
    gmp_printf("u = %Zd\nv = %Zd\n", u, v);
  }
  // TODO: we don't want to deal with ideals of degree > 1.
  // The simplest way is to forbid squares at all. This is not optimal,
  // and furthermore, we can do it only for small factors...
  int fail = 0;
  unsigned long x;
  x = mpz_gcd_ui(NULL, u, 4);   fail |= (x == 4);
  x = mpz_gcd_ui(NULL, u, 9);   fail |= (x == 9);
  x = mpz_gcd_ui(NULL, u, 25);  fail |= (x == 25);
  x = mpz_gcd_ui(NULL, u, 49);  fail |= (x == 49);
  x = mpz_gcd_ui(NULL, u, 121); fail |= (x == 121);
  x = mpz_gcd_ui(NULL, u, 169); fail |= (x == 169);
  x = mpz_gcd_ui(NULL, u, 289); fail |= (x == 289);
  x = mpz_gcd_ui(NULL, u, 361); fail |= (x == 361);
  x = mpz_gcd_ui(NULL, u, 529); fail |= (x == 529);
  x = mpz_gcd_ui(NULL, v, 4);   fail |= (x == 4);
  x = mpz_gcd_ui(NULL, v, 9);   fail |= (x == 9);
  x = mpz_gcd_ui(NULL, v, 25);  fail |= (x == 25);
  x = mpz_gcd_ui(NULL, v, 49);  fail |= (x == 49);
  x = mpz_gcd_ui(NULL, v, 121); fail |= (x == 121);
  x = mpz_gcd_ui(NULL, v, 169); fail |= (x == 169);
  x = mpz_gcd_ui(NULL, v, 289); fail |= (x == 289);
  x = mpz_gcd_ui(NULL, v, 361); fail |= (x == 361);
  x = mpz_gcd_ui(NULL, v, 529); fail |= (x == 529);

  LLL_clear(&M);
  mpz_poly_clear(U);
  mpz_poly_clear(V);
  mpz_clear(ze);
  return fail;
}

void next_cand_Fp_jl(cand_t cand, void *params) {
  Fp_param * param = (Fp_param *)params;
  mpz_t u, v;
  mpz_init(u);
  mpz_init(v);
  unsigned long e;
  int fail;
  do {
    e = random();
    fail = get_JL_candiate_from_e(e, u, v, param, 0);
  } while (fail); 

  cand_set_original_values(cand, u, v, e);
  mpz_clear(u);
  mpz_clear(v);
}


void usage(char *argv0) {
  fprintf(stderr, "./%s [-poly polfile] [-sqside xxx] [-jl] [-mt n] [-mineff e] [-maxeff E] [-seed s] [-target t] [-v] p z\n", argv0);
  abort();
}

typedef struct {
  Fp_param * params;
  const smooth_detect_param_s * smooth_param;
  unsigned long target;
  int jl;
} thparam_s;
typedef thparam_s* thparam_t;

void * process_one_thread(void *thparams) {
  thparam_t thparam = (thparam_s*) thparams;
  cand_t C;
  cand_init(C);
  if (!thparam->jl) {
    smooth_detect(C, next_cand_Fp_hgcd, thparam->params, thparam->target,
        thparam->smooth_param);
  } else {
    smooth_detect(C, next_cand_Fp_jl, thparam->params, thparam->target,
        thparam->smooth_param);
  }
  cand_s *winner = (cand_s*) malloc(sizeof(cand_s));
  cand_init(winner);
  cand_set(winner, C);
  return (void *)winner;
}

int main(int argc, char **argv) {
  char *argv0 = argv[0];
  unsigned long seed = 0;
  unsigned long target = 0;
  unsigned long nthread = 1;
  double mineff = 2000.0;
  double maxeff = 1e20;
  double minB1 = 100.0;
  int verbose = 0;
  int jl = 0;
  int sqside = 1;
  cado_poly cpoly;
  int gotpoly = 0;
  clock_t tm = clock();
  
  cado_poly_init(cpoly);

  while (argc > 2 && argv[1][0] == '-') {
    if (strcmp(argv[1], "-seed") == 0) {
      seed = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-mt") == 0) {
      nthread = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-minB1") == 0) {
      minB1 = strtof(argv[2], NULL);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-mineff") == 0) {
      mineff = strtof(argv[2], NULL);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-sqside") == 0) {
      sqside = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-maxeff") == 0) {
      maxeff = strtof(argv[2], NULL);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-v") == 0) {
      verbose = 1;
      argc -= 1;
      argv += 1;
    } else if (strcmp(argv[1], "-jl") == 0) {
      jl = 1;
      argc -= 1;
      argv += 1;
    } else if (strcmp(argv[1], "-target") == 0) {
      target = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-poly") == 0) {
      if (!cado_poly_read(cpoly, argv[2])) {
        fprintf(stderr, "Error reading polynomial file %s\n", argv[2]);
        usage(argv0);
      }
      gotpoly = 1;
      argc -= 2;
      argv += 2;
    } else {
      fprintf(stderr, "Unknown option %s\n", argv[1]);
      usage(argv0);
    }
  }
  if (argc != 3)
    usage(argv0);

  if (jl && !gotpoly) {
    fprintf(stderr, "Error, must provide -poly when using -jl option\n");
    usage(argv0);
  }

  if (seed == 0)
    seed = getpid() + (time(NULL)<<16);
  srandom(seed);

  Fp_param params;
  mpz_init_set_str(params.p, argv[1], 10);
  mpz_init_set_str(params.z, argv[2], 10);
  if (jl) {
    mpz_poly_init(params.f, cpoly->pols[sqside]->deg);
    mpz_poly_set(params.f, cpoly->pols[sqside]);
    mpz_init(params.m);
    int ret = cado_poly_getm(params.m, cpoly, params.p);
    ASSERT_ALWAYS(ret);
  }
  
  const smooth_detect_param_s smooth_param = {mineff, maxeff, 10, verbose, minB1};
  thparam_s thparam[1];
  thparam->params = &params;
  thparam->smooth_param = &smooth_param;
  thparam->target = target;
  thparam->jl = jl;

  pthread_t *thid = malloc(nthread*sizeof(pthread_t));
  for (unsigned int i = 0; i < nthread; ++i) {
    int ret = pthread_create(&thid[i], NULL, process_one_thread, thparam);
    assert(ret == 0);
  }
  
  while (1) {
    for (unsigned int i = 0; i < nthread; ++i) {
      cand_s *winner;
      // Non-blocking join: non-portable, gnu-specific.
      int ret = pthread_tryjoin_np(thid[i], (void *)(&winner));
      if (ret == 0) {
        printf("Youpi: e = %lu is a winner\n", winner->id);
        if (jl) {
          mpz_t u, v;
          mpz_init(u); mpz_init(v);
          // do again the LLL thing and print result.
          int fail = get_JL_candiate_from_e(winner->id, u, v, &params, 1);
          ASSERT_ALWAYS(!fail);
          mpz_clear(u); mpz_clear(v);
        }
        goto end;
      }
    }
    sleep(1);
  }

end:
  printf("Total CPU time: %.1f s\n", ((double)(clock() - tm)) / CLOCKS_PER_SEC);
  free(thid);
  cado_poly_clear(cpoly);
  mpz_clear(params.p);
  mpz_clear(params.z);
  return EXIT_SUCCESS;
}
