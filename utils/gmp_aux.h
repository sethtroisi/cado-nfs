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
extern uint64_t mpz_get_uint64 (mpz_t);
extern uint64_t uint64_nextprime (uint64_t);
extern unsigned long ulong_nextprime (unsigned long);
extern int isprime (unsigned long);


/* return the number of bits of p, counting from the least significant end */
extern int nbits (uintmax_t p);


#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_GMP_AUX_H_ */
