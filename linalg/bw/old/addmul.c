#include <stdio.h>
#include <string.h>
#include "slave_params.h"
#include "params.h"
#include <assert.h>
#include "types.h"
#include "macros.h"
#include "auxfuncs.h"
#include "variables.h"
#include "bw_scalar.h"
#include "bw_lvblock.h"
#include "slave.h"
#include "bw_lvblock_steps.h"

#if defined(__GNUC__) && defined(REALLY_WANT_ADDMUL_INLINES)
/* GCC is able to define addmul as a quasi-macro
 * Therefore, a standalone version is useless if we are optimizing.
 * Still, when not optimizing, extern inlines are not implemented as
 * inlines, and we need a separate version.
 */
#define INLINING_PREFIX try_inline
#define INCLUDING_ADDMUL_INLINES_H_
#include "addmul_inlines.h"
#else
#define INLINING_PREFIX
#define INCLUDING_ADDMUL_INLINES_H_
#include "addmul_inlines.h"
#endif
