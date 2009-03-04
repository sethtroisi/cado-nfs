#ifndef XDOTPROD_H_
#define XDOTPROD_H_

#include "matmul_top.h"
#include "abase.h"

/* This interface is relevant to both krylov and mksol, since it's used
 * both for checking and computing A files (obviously only krylov is
 * concerned by the latter aspect). */

#ifdef __cplusplus
extern "C" {
#endif

void x_dotprod(matmul_top_data_ptr mmt, uint32_t * xv, abt * v);

#ifdef __cplusplus
}
#endif

#endif	/* XDOTPROD_H_ */
