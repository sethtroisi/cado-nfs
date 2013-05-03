#include "cado.h"
#include "mod_mpz_default.h"
#define ecm ecm_mpz
#define ell_pointorder ell_pointorder_mpz
#define ellM_curveorder_jacobi ellM_curveorder_jacobi_mpz
#define ell_curveorder ell_curveorder_mpz
#define ellM_double ellM_double_mpz
#define ellM_add ellM_add_mpz
#define ellM_interpret_bytecode ellM_interpret_bytecode_mpz
#define ecm_stage2 ecm_stage2_mpz

#include "ecm.c"
