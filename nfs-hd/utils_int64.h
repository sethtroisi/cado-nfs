#ifndef UTILSINT64_T_H
#define UTILSINT64_T_H

#include <stdint.h>

/*
  Compute d^e.
  IN:
    d: int64_t, a number.
    e: int64_t, a number.
  OUT:
    an int64_t.
*/
int64_t pow_int64_t(int64_t d, int64_t e);

/*
  Compute d^e.
  IN:
    d: uint64_t, a number.
    e: int64_t, a number.
  OUT:
    an uint64_t.
*/
uint64_t pow_uint64_t(uint64_t d, int64_t e);

/*
  Compute the factorial naively.
  IN:
    res: int64_t * , the result.
    f: int64_t, a number.
*/
void factorial(int64_t * res, int64_t f);

/*
  Compute the absolute value of a and set in res.
  IN:
    res: uint64_t * , the absolute value of a.
    a: int64_t, a number.

 */
void int64_abs(uint64_t * res, int64_t a);

#endif /* UTILSINT64_T_H */
