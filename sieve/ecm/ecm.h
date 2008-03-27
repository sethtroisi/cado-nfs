#include "mod_ul.h"

#define BRENT12 0
#define MONTY12 1

int ecm_stage1 (residue_t, const int, const residue_t, const int, 
		const modulus_t);
unsigned long ell_pointorder (const residue_t, const int, const modulus_t, 
			      const int);
unsigned long ellM_curveorder_jacobi (residue_t, residue_t, modulus_t);
