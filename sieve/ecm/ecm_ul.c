#include "cado.h"
#include "modredc_ul_default.h"

#define ecm ecm_ul

#define ell_pointorder ell_pointorder_ul
#define ellM_curveorder_jacobi ellM_curveorder_jacobi_ul
#define ell_curveorder ell_curveorder_ul

/* The next ones are static so there's no need to rename them, but it's 
   nice to have these functions distinguishable e.g. in profiler output */
#define ellM_double ellM_double_ul
#define ellM_add ellM_add_ul
#define ellM_interpret_bytecode ellM_interpret_bytecode_ul
#define ecm_stage2 ecm_stage2_ul

#define ecmE ecmE_ul
#define ellE_double ellE_double_ul
#define ellE_add ellE_add_ul
#define ecmE_stage2 ecmE_stage2_ul

#include "ecm.c"
