#ifndef CADO_UTILS_GMP_AUX_H_
#define CADO_UTILS_GMP_AUX_H_

#include <gmp.h>
#include <stdint.h>

/* the following function is missing in GMP */
#ifndef mpz_addmul_si
#define mpz_addmul_si(a, b, c)                  \
  do {                                          \
    if (c >= 0)                                 \
      mpz_addmul_ui (a, b, c);                  \
    else                                        \
      mpz_submul_ui (a, b, -(c));               \
  }                                             \
  while (0)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* gmp_aux */
extern void mpz_set_uint64 (mpz_t, uint64_t);
extern void mpz_set_int64 (mpz_t, int64_t);
extern uint64_t mpz_get_uint64 (mpz_srcptr);
extern int64_t mpz_get_int64 (mpz_srcptr);
extern void mpz_mul_uint64 (mpz_t a, mpz_srcptr b, uint64_t c);
extern void mpz_mul_int64 (mpz_t a, mpz_srcptr b, int64_t c);
extern void mpz_addmul_int64 (mpz_t a, mpz_srcptr b, int64_t c);
extern int mpz_fits_int64_p(mpz_srcptr);
extern unsigned long ulong_nextprime (unsigned long);
extern int ulong_isprime (unsigned long);
extern void mpz_ndiv_q (mpz_t q, mpz_t n, mpz_t d);

/* return the number of bits of p, counting from the least significant end */
extern int nbits (uintmax_t p);


#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_GMP_AUX_H_ */
