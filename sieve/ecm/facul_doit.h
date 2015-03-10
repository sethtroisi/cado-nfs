#ifndef FACUL_DOIT_H
#define FACUL_DOIT_H

#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "mod_mpz.h"
#include "facul.h"

int primetest_ul (const modulusredcul_t m);
int primetest_15ul (const modulusredc15ul_t m);
int primetest_2ul2 (const modulusredc2ul2_t m);
int facul_doit_ul (mpz_t *, const modulusredcul_t, 
		   const facul_strategy_t *, const int);
int facul_doit_15ul (mpz_t *, const modulusredc15ul_t, 
		     const facul_strategy_t *, const int);
int facul_doit_2ul2 (mpz_t *, const modulusredc2ul2_t, 
		     const facul_strategy_t *, const int);
int facul_doit_mpz (mpz_t *, const modulusmpz_t, 
		    const facul_strategy_t *, const int);

int
facul_doit_onefm_ul (mpz_t*, const modulusredcul_t,
		     const facul_method_t, modset_t*,
		     modset_t*, unsigned long, double, double);
int
facul_doit_onefm_15ul (mpz_t*, const modulusredc15ul_t,
		       const facul_method_t, modset_t*,
		       modset_t*, unsigned long, double, double);
int
facul_doit_onefm_2ul2 (mpz_t*, const modulusredc2ul2_t,
		       const facul_method_t, modset_t*,
		       modset_t*, unsigned long, double, double);
int
facul_doit_onefm_mpz (mpz_t*, const modulusmpz_t,
		      const facul_method_t, modset_t*,
		      modset_t*, unsigned long, double, double);

/* int* */
/* facul_both (unsigned long**, mpz_t* , */
/* 	    const facul_strategies_t *); */


#endif /* FACUL_DOIT_H */
