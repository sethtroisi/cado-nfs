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

unsigned long ell_pointorder_ul (const unsigned long,
                                 const ec_parameterization_t,
                                 const unsigned long, const unsigned long,
                                 const modulusredcul_t, const int);

unsigned long ellM_curveorder_jacobi_ul (residueredcul_t, residueredcul_t, \
                                         modulusredcul_t);

#if 0 /* FIXME Currently, only the _ul version of this function works */
unsigned long ell_pointorder_15ul (const unsigned long,
                                   const ec_parameterization_t,
                                   const unsigned long, const unsigned long,
                                   const modulusredc15ul_t, const int);
#endif

unsigned long ellM_curveorder_jacobi_15ul (residueredc15ul_t, residueredc15ul_t, \
                                           modulusredc15ul_t);

#if 0 /* FIXME Currently, only the _ul version of this function works */
unsigned long ell_pointorder_2ul2 (const unsigned long,
                                   const ec_parameterization_t,
                                   const unsigned long, const unsigned long,
                                   const modulusredc2ul2_t, const int);
#endif

unsigned long ellM_curveorder_jacobi_2ul2 (residueredc2ul2_t, residueredc2ul2_t, 
                                           modulusredc2ul2_t);

#if 0 /* FIXME Currently, only the _ul version of this function works */
unsigned long ell_pointorder_mpz (const unsigned long,
                                  const ec_parameterization_t,
                                  const unsigned long, const unsigned long,
                                  const modulusmpz_t, const int);
#endif

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
