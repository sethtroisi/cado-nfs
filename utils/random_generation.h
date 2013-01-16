#ifndef RANDOM_H_
#define RANDOM_H_

#include <gmp.h>

/* Define USE_GMP_RANDOM to enable GMP's random number generator as a
 * backend for the randomness source
 */

#define xxxUSE_GMP_RANDOM


#ifdef __cplusplus
extern "C" {
#endif

extern mp_limb_t myrand ();
extern void myrand_area(void*, size_t);
extern void setup_seeding(unsigned int);

#ifdef __cplusplus
}
#endif

#endif	/* RANDOM_H_ */
