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

/* maximum lognorm+exp_E increment for each rotation */
#define BOUND_LOGNORM_INCR_MAX 1.005

#define BOUND_LOGNORM_INCR_MAX_TUNESTEP 0.005

#define TUNE_LOGNORM_INCR 1

#define MAX_LINE_LENGTH 4096

/* number of default parameters for various sublattices */
#define NUM_DEFAULT_SUBLATTICE 26

#define DEBUG 0

/* Don't change this. */
#define NUM_SUBLATTICE_PRIMES 9

/* Top 32 (alpha) in each "SIZE_SIEVEARRAY": this does not
   affect the running-time since it is not dominating */
#define NUM_TOPALPHA_SIEVEARRAY 16

/* Top 32 (E) for each sublattice, may contain several
   SIZE_SIEVEARRAY */
#define NUM_TOPE_SUBLATTICE 16

/* Either rank stage 1 sublattices by E or by alpha, the former
   seems to be more accurate */
#define RANK_SUBLATTICE_BY_E 1


/* -------------------------
   Parameters for test sieving
   ------------------------- */

/* This is the maxmum of the following values */
#define TUNE_NUM_SUBLATTICE 8

/* allow more sublattice in Stage 1 (ranking based on alpha values) */
#define TUNE_NUM_SUBLATTICE_STAGE1 8

#define TUNE_NUM_SUBLATTICE_STAGE2 2

/* Similar to above, but in tune mode */
#define TUNE_NUM_TOPALPHA_SIEVEARRAY 8

/* Similar to above, but in tune mode */
#define TUNE_NUM_TOPE_SUBLATTICE 8

/* Tuning parameter */
#define TUNE_EARLY_ABORT 1

#define TUNE_EARLY_ABORT_THR 2

#define TUNE_BOUND_ON_UV_TRIALS 64

#define TUNE_BOUND_ON_MOD_TRIALS 64


/* mu and sigma in function exp_alpha() */
#define MU 0.0
#define SIGMA 0.824


/* --- declarations --- */


extern unsigned int L1_cachesize;

/* set once in ropt_io.c and unchanged then */
extern unsigned int size_tune_sievearray;

extern const unsigned int primes[];

extern const unsigned char next_prime_idx[];

extern const unsigned int default_sublattice_pe[NUM_DEFAULT_SUBLATTICE][NUM_SUBLATTICE_PRIMES];

extern const unsigned long default_sublattice_prod[NUM_DEFAULT_SUBLATTICE];

extern const unsigned int s1_size_each_sublattice[NUM_SUBLATTICE_PRIMES][NUM_SUBLATTICE_PRIMES];

extern const unsigned int s1_size_each_sublattice_tune[NUM_SUBLATTICE_PRIMES];

extern const unsigned int size_total_sublattices[9][3];

double exp_alpha (double logK);

#endif
