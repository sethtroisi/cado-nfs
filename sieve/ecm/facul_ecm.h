#ifndef ECM_H_
#define ECM_H_

#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "mod_mpz.h"
#include "bytecode.h"
#include "stage2.h"

typedef enum {
  BRENT12 = 1,
  MONTY12 = 2,
  MONTY16 = 4,
  MONTYTWED12 = 8,
  MONTYTWED16 = 16
} ec_parameterization_t;

#define FULLMONTY (BRENT12 | MONTY12 | MONTY16)
#define FULLMONTYTWED (MONTYTWED12 | MONTYTWED16)
#define ECM_TORSION16 (MONTY16 | MONTYTWED16)
#define ECM_TORSION12 (BRENT12 | MONTY12 | MONTYTWED12)


/* Twisted Edwards curve with a = -1 */
typedef struct {
  long g_numer, g_denom;
  long d_numer, d_denom;        /* d parameter for twisted Edwards curves */
  long x_numer, x_denom, 
       y_numer, y_denom;        /* non-torsion point on the curve */
} Edwards_curve_t;

static const Edwards_curve_t Ecurve14 = {
  1, 4, 
  -50625, 4096, 
  104329, 16384, 
  1630827, 262144
};

/* This is not ok for 32-bit builds. Apparently it's not used for the
 * moment anyway, so that's not much of a problem. But what should we do?
 * Use int64_t coefficients (implies things such as mul_si64 down the
 * line), or stick to ulongs and disable this curve for 32-bits ?
 *
 * [Laurent:] Well, most of the ECM friendly curves have even larger
 * coefficients. We'll have to use mpz_t anyways. Ecurve14 is just a toy example.
 */

 /* static const Edwards_curve_t Ecurve45 = { */
 /*   4, 5, */
 /*   -6561, 2560000, */
 /*   106564329, 501760000, */
 /*   -3715030917, 280985600000 */
 /* }; */


typedef struct {
  bytecode bc;          /* Bytecode for stage 1 */
  unsigned int exp2;    /* Exponent of 2 in stage 1 primes */
  unsigned int B1;
  ec_parameterization_t parameterization;
  unsigned long parameter;  /* Used to compute curve coefficients with the
                             * above parameterization. */
  stage2_plan_t stage2;
} ecm_plan_t;


#ifdef __cplusplus
extern "C" {
#endif

int ecm_ul (modintredcul_t, const modulusredcul_t, const ecm_plan_t *);
int ecm_15ul (modintredc15ul_t, const modulusredc15ul_t, const ecm_plan_t *);
int ecm_2ul2 (modintredc2ul2_t, const modulusredc2ul2_t, const ecm_plan_t *);
int ecm_mpz (modintmpz_t, const modulusmpz_t, const ecm_plan_t *);

unsigned long ell_pointorder_ul (const residueredcul_t,
                                 const ec_parameterization_t,
                                 const unsigned long, const unsigned long,
                                 const modulusredcul_t, const int);

unsigned long ellM_curveorder_jacobi_ul (residueredcul_t, residueredcul_t, \
                                         modulusredcul_t);

unsigned long ell_pointorder_15ul (const residueredc15ul_t,
                                   const ec_parameterization_t,
                                   const unsigned long, const unsigned long,
                                   const modulusredc15ul_t, const int);

unsigned long ellM_curveorder_jacobi_15ul (residueredc15ul_t, residueredc15ul_t, \
                                           modulusredc15ul_t);

unsigned long ell_pointorder_2ul2 (const residueredc2ul2_t,
                                   const ec_parameterization_t,
                                   const unsigned long, const unsigned long,
                                   const modulusredc2ul2_t, const int);

unsigned long ellM_curveorder_jacobi_2ul2 (residueredc2ul2_t, residueredc2ul2_t, 
                                           modulusredc2ul2_t);

unsigned long ell_pointorder_mpz (const residuempz_t,
                                  const ec_parameterization_t,
                                  const unsigned long, const unsigned long,
                                  const modulusmpz_t, const int);

unsigned long ellM_curveorder_jacobi_mpz (residuempz_t, residuempz_t, 
                                           modulusmpz_t);


void ecm_make_plan (ecm_plan_t *, const unsigned int, const unsigned int, 
		    const ec_parameterization_t, const unsigned long, const int, const int);

void ecm_clear_plan (ecm_plan_t *);

int ec_parameter_is_valid_ul (ec_parameterization_t, const unsigned long);
int ec_parameter_is_valid_15ul (ec_parameterization_t, const unsigned long);
int ec_parameter_is_valid_2ul2 (ec_parameterization_t, const unsigned long);
int ec_parameter_is_valid_mpz (ec_parameterization_t, const unsigned long);
/* In fact the ec_parameter_is_valid_* functions do not use modredc arithmetic,
 * we can use any of the four.
 */
#ifndef ec_parameter_is_valid
#define ec_parameter_is_valid ec_parameter_is_valid_ul
#endif

#ifdef __cplusplus
}
#endif

#endif	/* ECM_H_ */
