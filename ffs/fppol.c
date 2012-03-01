#define _POSIX_C_SOURCE 199309L
#include <stdlib.h>
#include <stdarg.h>

#include "fppol.h"
#include "fppol_internal.h"
#include "macros.h"



/* Multiprecision polynomials.
 *****************************************************************************/

// Initialize polynomial.
void fppol_init(fppol_ptr r)
{
  r->alloc = 0;
  r->limbs = NULL;
}


// Initialize a NULL-terminated list of polynomials.
void fppol_inits(fppol_ptr r, ...)
{
  va_list ap;
  for (va_start(ap, r); r != NULL; r = va_arg(ap, fppol_ptr))
    fppol_init(r);
  va_end(ap);
}


// Initialize polynomial with space for n terms, ie degree n-1.
void fppol_init2(fppol_ptr r, unsigned n)
{
  r->alloc = (n+63)>>6;
  r->limbs = malloc(r->alloc * sizeof(fppol64_t));
  ASSERT_ALWAYS(!n || r->limbs != NULL);
}


// Free polynomial.
void fppol_clear(fppol_ptr r)
{
  free(r->limbs);
}


// Free a NULL-terminated list of polynomials.
void fppol_clears(fppol_ptr r, ...)
{
  va_list ap;
  for (va_start(ap, r); r != NULL; r = va_arg(ap, fppol_ptr))
    fppol_clear(r);
  va_end(ap);
}


// Reallocate polynomial so as to free unused limbs.
void fppol_trim(fppol_ptr r)
{
  __fppol_realloc(r, r->deg+1);
}









#ifdef FPPOL_MAIN

#define __USE_SVID
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>

#if   defined(USE_F2)
#  define GET(z) z[0] = GET_WORD;
#  define PUT(z) printf("0x%016lx", z[0]);
#elif defined(USE_F3)
#  define GET(z) z[0] = GET_WORD; z[1] = GET_WORD; \
                 { uint64_t m = z[0] & z[1]; z[0] ^= m; z[1] ^= m; }
#  define PUT(z) printf("0x%016lx 0x%016lx", z[0], z[1]);
#endif

#if 0

#define N 100000000

#define GET_WORD mrand48()<<32 ^ mrand48()

int main()
{
  fppol64_t x = {0}, y = {0}, r, s, *p, *q;

  p = malloc(0x400 * sizeof(fppol64_t));
  q = malloc(0x400 * sizeof(fppol64_t));

  srand48(time(NULL));
  for (unsigned i = 0; i < 0x400; ++i) {
    GET(p[i]);
    GET(q[i]);
  }

  struct timespec user1, user2;

  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &user1);
  for (unsigned i = 0; i < N; ++i) {
    __FP_MUL_128_64x64(r, s, p[i&0x3ff], q[i&0x3ff]);
    fppol64_add(x, x, r);
    fppol64_add(y, y, s);
  }
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &user2);

  double duser = (double)(user2.tv_sec  - user1.tv_sec) +
                 (double)(user2.tv_nsec - user1.tv_nsec) / 1.0e9;
  fprintf(stderr, "Time: %0.3fs.\n", duser);

  PUT(x); printf("\n");
  PUT(y); printf("\n");

  return EXIT_SUCCESS;
}

#else

#define GET_WORD strtoul(*++argv, NULL, 0)

int main(int argc, char **argv)
{
  if (argc <= 2)
    return EXIT_FAILURE;

  fppol64_t r, s, p, q;
  int ret;
  ret = fppol64_set_str(p, *++argv);
  if (!ret) {
      fprintf(stderr, "p is not valid\n");
      exit(EXIT_FAILURE);
  }
  ret = fppol64_set_str(q, *++argv);
  if (!ret) {
      fprintf(stderr, "q is not valid\n");
      exit(EXIT_FAILURE);
  }

  printf("p      = "); fppol64_out(stdout, p); printf("\n");
  printf("q      = "); fppol64_out(stdout, q); printf("\n");

  fppol64_set_zero(r);
  printf("0      = "); fppol64_out(stdout, r); printf("\n");

  fppol64_add(r, p, q);
  printf("p + q  = "); fppol64_out(stdout, r); printf("\n");

  fppol64_sub(r, p, q);
  printf("p - q  = "); fppol64_out(stdout, r); printf("\n");

  if (fppol64_deg(p) < fppol64_deg(q)) {
    fppol64_shl1mod(r, p, q);
    printf("p*t mod q  = "); 
    fppol64_out(stdout, r); printf("\n");
    fppol64_shl1mod(r, r, q);
    printf("p*t^2 mod q  = "); 
    fppol64_out(stdout, r); printf("\n");

    ret = fppol64_invmod(r, p, q);
    printf("1/p mod q  = "); 
    if (!ret)
      printf("N/A\n");
    else {
      fppol64_out(stdout, r); printf("\n");
      fppol64_mulmod(s, p, r, q);
      printf("mulmod o invmod = "); fppol64_out(stdout, s); printf("\n");
    }
  }

  __FP_MUL_128_64x64(r, s, p, q, );
  printf("p * q  = ");
  fppol64_out(stdout, r); printf(" ");
  fppol64_out(stdout, s); printf("\n");

  if (!fppol64_is_zero(q)) {
    // compute quotient and remainder
    fppol64_divrem(r, s, p, q);
    printf("p / q = "); fppol64_out(stdout, r);
    printf(" rem = ");  fppol64_out(stdout, s);
    printf("\n");

    // check the result
    fppol64_mul(r, r, q);
    fppol64_add(r, r, s);
    fppol64_sub(r, r, p);
    printf("p - (q*r+s) = "); fppol64_out(stdout, r); printf("\n");
  }

  fppol64_gcd(r, p, q);
  printf("gcd(p, q)  = "); fppol64_out(stdout, r); printf("\n");


  fppol_t mp, mq, mr, ms, ma, mb;
  fppol_init(mp);
  fppol_init(mq);
  fppol_init(mr);
  fppol_init(ms);
  fppol_init(ma);
  fppol_init(mb);

  ret = fppol_set_str(mp, *++argv);
  if (!ret) {
      fprintf(stderr, "p is not valid\n");
      exit(EXIT_FAILURE);
  }
  ret = fppol_set_str(mq, *++argv);
  if (!ret) {
      fprintf(stderr, "q is not valid\n");
      exit(EXIT_FAILURE);
  }

  printf("deg(p) = %d\n", fppol_deg(mp));
  printf("deg(q) = %d\n", fppol_deg(mq));
  fppol_add(mr, mp, mq);
  printf("p + q  = "); fppol_out(stdout, mr); printf(" [%d]\n", fppol_deg(mr));
  fppol_sub(mr, mp, mq);
  printf("p - q  = "); fppol_out(stdout, mr); printf(" [%d]\n", fppol_deg(mr));
  fppol_mul(ms, mp, mq);
  printf("p * q  = "); fppol_out(stdout, ms); printf(" [%d]\n", fppol_deg(ms));
  fppol_mul(ms, ms, mp);
  printf("p^2*q  = "); fppol_out(stdout, ms); printf(" [%d]\n", fppol_deg(ms));
  fppol_mul(ms, mq, ms);
  printf("p^2*q^2= "); fppol_out(stdout, ms); printf(" [%d]\n", fppol_deg(ms));

  fppol_divrem(ma, mb, ms, mr);
  printf("s / r  = "); fppol_out(stdout, ma); printf(" [%d]\n", fppol_deg(ma));
  printf("s %% r  = ");fppol_out(stdout, mb); printf(" [%d]\n", fppol_deg(mb));
  fppol_mul(ma, ma, mr);
  fppol_add(ma, ma, mb);
  fppol_sub(ma, ms, ma);
  printf("check  = "); fppol_out(stdout, ma); printf(" [%d]\n", fppol_deg(ma));
  
  fppol_gcd(mr, mp, mq);
  printf("gcd(p, q)  = "); fppol_out(stdout, mr);
  printf(" [%d]\n", fppol_deg(mr));

  ret = fppol_invmod(ma, mr, ms);
  printf("1/r %%s = ");
  if (!ret) printf("N/A\n");
  else      fppol_out(stdout, ma), printf(" [%d]\n", fppol_deg(ma));

  fppol_shl(mr, mp, 1);
  printf("p << 1 = "); fppol_out(stdout, mr); printf(" [%d]\n", fppol_deg(mr));
  fppol_shl(mr, mp, 2);
  printf("p << 2 = "); fppol_out(stdout, mr); printf(" [%d]\n", fppol_deg(mr));
  fppol_shl(mr, mp, 127);
  printf("p <<127= "); fppol_out(stdout, mr); printf(" [%d]\n", fppol_deg(mr));
  fppol_shr(mr, ms, 1);
  printf("s >> 1 = "); fppol_out(stdout, mr); printf(" [%d]\n", fppol_deg(mr));
  fppol_shr(mr, ms, 2);
  printf("s >> 2 = "); fppol_out(stdout, mr); printf(" [%d]\n", fppol_deg(mr));
  fppol_shr(mr, ms, 127);
  printf("s >>127= "); fppol_out(stdout, mr); printf(" [%d]\n", fppol_deg(mr));

  fppol_clear(mp);
  fppol_clear(mq);
  fppol_clear(mr);
  fppol_clear(ms);
  fppol_clear(ma);
  fppol_clear(mb);

  unsigned n = 0;
  for (uint64_t i = 0, j; i < 243<<8; ++i) {
    if (!fppol64_set_ui(r, i, 6))
      continue;
    if ((j = fppol64_get_ui(r, 6)) != i) {
      ++n;
      printf("%lu: ", i);
      fppol64_out(stdout, r);
      printf(" != %lu\n", j);
    }
  }
  for (uint64_t i = 0, j; i < 122<<8; ++i) {
    if (!fppol64_monic_set_ui(r, i, 6))
      continue;
    if ((j = fppol64_monic_get_ui(r, 6)) != i) {
      ++n;
      printf("%lu: ", i);
      fppol64_out(stdout, r);
      printf(" != %lu\n", j);
    }
  }
  printf("%u conversion errors\n", n);

  return EXIT_SUCCESS;
}

#endif
#endif
