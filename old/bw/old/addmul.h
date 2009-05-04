#ifndef ADDMUL_H_
#define ADDMUL_H_

#if defined(__GNUC__) && defined(REALLY_WANT_ADDMUL_INLINES)
#define INLINING_PREFIX		extern try_inline
#define INCLUDING_ADDMUL_INLINES_H_
#include "addmul_inlines.h"
#else

#ifdef	__cplusplus
extern "C" {
#endif

extern void reduce_p_minus_q(mp_limb_t *, mp_limb_t *, mp_limb_t *);
extern void addmultiply(bw_vector_block, bw_vector_block, stype32);

#ifdef	__cplusplus
}
#endif

#endif

#endif	/* ADDMUL_H_ */
