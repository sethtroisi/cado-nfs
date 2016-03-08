#include "cado.h"
#include "modredc_15ul_default.h"

#define ecm ecm_15ul
#define ell_pointorder ell_pointorder_15ul
#define ellM_curveorder_jacobi ellM_curveorder_jacobi_15ul
#define ell_curveorder ell_curveorder_15ul
#define ellM_double ellM_double_15ul
#define ellM_add ellM_add_15ul
#define ellM_interpret_bytecode ellM_interpret_bytecode_15ul
#define ecm_stage2 ecm_stage2_15ul

#define ecmE ecmE_15ul
#define ellE_double ellE_double_15ul
#define ellE_add ellE_add_15ul
#define ecmE_stage2 ecmE_stage2_15ul

#include "ecm.c"
