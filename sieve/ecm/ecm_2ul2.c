#include "cado.h"
#include "modredc_2ul2_default.h"

#define ecm ecm_2ul2
#define ell_pointorder ell_pointorder_2ul2
#define ellM_curveorder_jacobi ellM_curveorder_jacobi_2ul2
#define ell_curveorder ell_curveorder_2ul2
#define ellM_double ellM_double_2ul2
#define ellM_add ellM_add_2ul2
#define ellM_interpret_bytecode ellM_interpret_bytecode_2ul2
#define ecm_stage2 ecm_stage2_2ul2

#define ecmE ecmE_2ul2
#define ellE_double ellE_double_2ul2
#define ellE_add ellE_add_2ul2
#define ecmE_stage2 ecmE_stage2_2ul2

#include "ecm.c"
