#ifndef ROPT_PARAM_H
#define ROPT_PARAM_H


#include "ropt_str.h"


/* ------------------
   possible to change
   ------------------ */

/* Stage 2 memory cutoff: at most, use 2^28 of int16_t
   which is about 256M memory. Large sieving region will
   be segmented into smaller chunks. */
#define SIZE_SIEVEARRAY 268435456

/* Max default sieve length of v */
#define SIZE_SIEVEARRAY_V_MAX 4194304

/* Min default sieve length of v */
#define SIZE_SIEVEARRAY_V_MIN 1000000

/* Max default sieve length of u */
#define SIZE_SIEVEARRAY_U 64


/* -------------------------
   perhaps no need to change
   ------------------------- */

#define PI 3.14159265358979

#define SUP_ALPHA 4.843

#define BOUND_LOGNORM_RATIO 1.05

#define MAX_LINE_LENGTH 4096

#define MAX_DEGREE 6

/* number of default parameters for various sublattices */
#define NUM_DEFAULT_SUBLATTICE 26

#define DEBUG 0

/* Skip ropt, only used for debugging purpose. */
#define SKIP_ROPT 0

/* Don't change this. */
#define NUM_SUBLATTICE_PRIMES 9

/* Top 16 (alpha) in each "SIZE_SIEVEARRAY" */
#define NUM_TOPALPHA_SIEVEARRAY 16

/* Top 16 (E) for each sublattice, may contain several
   SIZE_SIEVEARRAY */
#define NUM_TOPE_SUBLATTICE 16

/* Either rank stage 1 sublattices by E or by alpha, the former
   seems to be more accurate */
#define RANK_SUBLATTICE_BY_E 1


/* --- declarations --- */


extern unsigned int L1_cachesize;

extern unsigned int size_tune_sievearray;

extern const unsigned int primes[];

extern const unsigned char next_prime_idx[];

extern const double exp_alpha[];

extern const unsigned int default_sublattice_pe[NUM_DEFAULT_SUBLATTICE][NUM_SUBLATTICE_PRIMES];

extern const unsigned long default_sublattice_prod[NUM_DEFAULT_SUBLATTICE];

extern const unsigned int size_each_sublattice[NUM_SUBLATTICE_PRIMES][NUM_SUBLATTICE_PRIMES];

extern const unsigned int size_each_sublattice_tune[NUM_SUBLATTICE_PRIMES];

extern const unsigned int size_total_sublattices[8][2];

#endif
