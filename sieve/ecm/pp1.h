#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "stage2.h"

typedef struct {
  char *bc;             /* Bytecode for the Lucas chain for stage 1 */
  unsigned int bc_len;  /* Number of bytes in bytecode */
  unsigned int B1;
  stage2_plan_t stage2;
} pp1_plan_t;


int pp1_ul (modintredcul_t, const modulusredcul_t, const pp1_plan_t *);
void pp1_stage2_ul (residueredcul_t, const residueredcul_t, 
                 const stage2_plan_t *, const residueredcul_t, 
                 const modulusredcul_t);
int pp1_15ul (modintredc15ul_t, const modulusredc15ul_t, const pp1_plan_t *);
void pp1_stage2_15ul (residueredc15ul_t, const residueredc15ul_t, 
                    const stage2_plan_t *, const residueredc15ul_t, 
                    const modulusredc15ul_t);
void pp1_make_plan (pp1_plan_t *, const unsigned int, const unsigned int, int);
void pp1_clear_plan (pp1_plan_t *);
