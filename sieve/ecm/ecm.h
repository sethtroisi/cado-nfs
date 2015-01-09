#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "mod_mpz.h"
#include "stage2.h"

#define BRENT12   1
#define MONTY12   2
#define MONTY16   4
#define FULLMONTY 7
#define TWED12    8
#define TWED16   16


/* Twisted Edwards curve with a = -1 */
typedef struct {
  unsigned long g_numer, g_denom;
  unsigned long d_numer, d_denom;        /* d parameter for twisted Edwards curves */
  unsigned long x_numer, x_denom, 
                y_numer, y_denom;        /* non-torsion point on the curve */
} Edwards_curve_t;

static const Edwards_curve_t Ecurve14 = {1, 4, -50625, 4096, 104329, 16384, 1630827, 262144};

/* This is not ok for 32-bit builds. Apparently it's not used for the
 * moment anyway, so that's not much of a problem. But what should we do?
 * Use int64_t coefficients (implies things such as mul_si64 down the
 * line), or stick to ulongs and disable this curve for 32-bits ?
static const Edwards_curve_t Ecurve45 = {4, 5, -6561, 2560000, 106564329, 501760000, -3715030917, 280985600000};
 */

typedef struct {
  char *bc;             /* Bytecode for the Lucas chain for stage 1 */
  unsigned int bc_len;  /* Number of bytes in bytecode */
  unsigned int exp2;    /* Exponent of 2 in stage 1 primes */
  unsigned int B1;
  int parameterization; /* BRENT12 or MONTY12 */
  unsigned long sigma;  /* Sigma parameter for Brent curves, or
			   multiplier for Montgomery torsion-12 curves */

  stage2_plan_t stage2;
} ecm_plan_t;


typedef struct {
  unsigned int exp2;
  unsigned int B1;
  int parameterization; /* TWED12 or TWED16 */
  const Edwards_curve_t *E;   /* Parameters for Edwards curve */
} ecmE_plan_t;


int ecm_ul (modintredcul_t, const modulusredcul_t, const ecm_plan_t *);
int ecm_15ul (modintredc15ul_t, const modulusredc15ul_t, const ecm_plan_t *);
int ecm_2ul2 (modintredc2ul2_t, const modulusredc2ul2_t, const ecm_plan_t *);
int ecm_mpz (modintmpz_t, const modulusmpz_t, const ecm_plan_t *);

unsigned long ell_pointorder_ul (const residueredcul_t, const int, \
                                 const unsigned long, const unsigned long, \
                                 const modulusredcul_t, const int);

unsigned long ellM_curveorder_jacobi_ul (residueredcul_t, residueredcul_t, \
                                         modulusredcul_t);

unsigned long ell_pointorder_15ul (const residueredc15ul_t, const int, \
                                   const unsigned long, const unsigned long, \
                                   const modulusredc15ul_t, const int);

unsigned long ellM_curveorder_jacobi_15ul (residueredc15ul_t, residueredc15ul_t, \
                                           modulusredc15ul_t);

unsigned long ell_pointorder_2ul2 (const residueredc2ul2_t, const int, \
                                   const unsigned long, const unsigned long,\
                                   const modulusredc2ul2_t, const int);

unsigned long ellM_curveorder_jacobi_2ul2 (residueredc2ul2_t, residueredc2ul2_t, 
                                           modulusredc2ul2_t);

unsigned long ell_pointorder_mpz (const residuempz_t, const int, 
                                   const unsigned long, const unsigned long,
                                   const modulusmpz_t, const int);

unsigned long ellM_curveorder_jacobi_mpz (residuempz_t, residuempz_t, 
                                           modulusmpz_t);

void ecm_make_plan (ecm_plan_t *, const unsigned int, const unsigned int, 
		    const int, const unsigned long, const int, const int);

void ecm_clear_plan (ecm_plan_t *);


void ecmE_make_plan (ecmE_plan_t *, const unsigned int, const int); 

void ecmE_clear_plan (ecmE_plan_t *);
