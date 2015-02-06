#include "utils_int64.h"

/* int64_t pow_int64_t(int64_t d, int64_t e) */
/* { */
/*   int i = 0; */
/*   int64_t tmp = 1; */
/*   for ( ; i < e; i++) { */
/*     tmp = tmp * d; */
/*   } */
/*   return tmp; */
/* } */

void factorial(int64_t * res, int64_t f)
{
  * res = 1;
  for (int64_t i = 2; i <= f; i++) {
    *res = *res * i;
  }
}

void int64_abs(uint64_t * res, int64_t a)
{
  if (a < 0) {
    * res = (uint64_t) -a;
    return;
  }
  * res = (uint64_t) a;
}

int64_t pow_int64_t(int64_t d, int64_t e)
{
  if (e == 1) {
    return d;
  } else if ((e & 1) == 1) {
    return d * pow_int64_t(d * d, (e - 1) / 2);
  } else {
    return pow_int64_t(d * d, e / 2);
  }
}

uint64_t pow_uint64_t(uint64_t d, int64_t e)
{
  if (e == 1) {
    return d;
  } else if ((e & 1) == 1) {
    return d * pow_int64_t(d * d, (e - 1) / 2);
  } else {
    return pow_int64_t(d * d, e / 2);
  }
}

uint64_t invmod_uint64(uint64_t xx, uint64_t mm)
{
  uint64_t yy;
  modulus_t m;
  residue_t x, y;

  ASSERT_ALWAYS(LONG_BIT == 64);

  modul_initmod_ul(m, mm);        //we work mod mm
  modul_init_noset0(x, m);
  modul_init_noset0(y, m);
  modul_set_ul_reduced(x, xx, m);  // we know that xx < m

  modul_inv(y, x, m);

  yy = modul_get_ul(y, m);
  modul_clear(x, m);
  modul_clear(y, m);
  modul_clearmod(m);
  return yy;
}
