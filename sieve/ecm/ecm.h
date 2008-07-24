#include "modredc_ul.h"
#include "stage2.h"

#define BRENT12 0
#define MONTY12 1

typedef struct {
  char *bc;             /* Bytecode for the Lucas chain for stage 1 */
  unsigned int bc_len;  /* Number of bytes in bytecode */
  unsigned int exp2;    /* Exponent of 2 in stage 1 primes */
  unsigned int B1;
  int parameterization;
  unsigned long sigma;  /* Sigma parameter for Brent curves, or
			   multiplier for Montgomery torsion-12 curves */
  stage2_plan_t stage2;
} ecm_plan_t;


void ecm_make_plan (ecm_plan_t *, const unsigned int, const unsigned int, 
		    const int, const unsigned long, const int);
void ecm_clear_plan (ecm_plan_t *);
unsigned long ecm (residue_t, const modulus_t, const ecm_plan_t *);

unsigned long ell_pointorder (const residue_t, const int, const modulus_t, 
			      const int);
unsigned long ellM_curveorder_jacobi (residue_t, residue_t, modulus_t);
