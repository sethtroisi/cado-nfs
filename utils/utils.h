#ifndef	CADO_UTILS_H_
#define	CADO_UTILS_H_

#include "cado.h"
#include "mod_ul.h"
#include <stdint.h>

#include <limits.h>

/* It's awful, but I confess that this ULONG_BITS is not ``portable'',
 * norm-wise.  GMP_LIMB_BITS is at hand, but could differ. An #ifdef
 * switch depending on macros like __x86_64 is considerably more fragile.
 */
#define	ULONG_BITS	((int) (sizeof(unsigned long) * CHAR_BIT))

#include "fpoly.h"
// #include "plain_poly.h"
// #include "modul_poly.h"
#include "getprime.h"
#include "timing.h"
#include "relation.h"
#include "gmp_aux.h"
#include "cado_poly.h"
#include "rootfinder.h"
#include "params.h"
#include "gcd_int64.h"
#include "gcd_uint64.h"
#include "discriminant.h"
#include "random_generation.h"
#include "mpz_array.h"
#include "gzip.h"
#include "hashpair.h"
#include "misc.h"

#endif	/* CADO_UTILS_H_ */
