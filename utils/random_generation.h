#ifndef RANDOM_H_
#define RANDOM_H_

#include <gmp.h>


/* Cado has its own seed-able random generator, which can either use the
 * libc generator, or GMP's random state. The USE_GMP_RANDOM flags
 * controls which is used.
 */

#define xxxUSE_GMP_RANDOM

struct cado_random_state_s {
#ifdef  USE_GMP_RANDOM
    gmp_randstate_t g[1];
#endif
};
typedef struct cado_random_state_s cado_random_state[1];
typedef struct cado_random_state_s * cado_random_state_ptr;
typedef const struct cado_random_state_s * cado_random_state_srcptr;


#ifdef __cplusplus
extern "C" {
#endif

extern mp_limb_t cado_random (cado_random_state_ptr);
extern void cado_random_area(cado_random_state_ptr, void*, size_t);
extern void cado_random_init(cado_random_state_ptr, unsigned int);
extern void cado_random_clear(cado_random_state_ptr);

#ifdef __cplusplus
}
#endif

#endif	/* RANDOM_H_ */
