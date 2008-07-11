#include "modredc_ul.h"
#include "stage2.h"

typedef struct {
  char *bc;             /* Bytecode for the Lucas chain for stage 1 */
  unsigned int bc_len;  /* Number of bytes in bytecode */
  unsigned int B1;
  stage2_plan_t stage2;
} pp1_plan_t;


unsigned long pp1 (residue_t x, const modulus_t, const pp1_plan_t *);
void pp1_stage1 (residue_t, const char *, const unsigned int, 
		 const residue_t, const modulus_t);
void pp1_make_plan (pp1_plan_t *, const unsigned int, const unsigned int, int);
void pp1_clear_plan (pp1_plan_t *);
