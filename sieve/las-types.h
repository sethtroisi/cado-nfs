#ifndef LAS_TYPES_H_
#define LAS_TYPES_H_

#include <stdint.h>
#include "fb.h"
#include "trialdiv.h"
#include "las-config.h"
#include "cado_poly.h"
#include "ecm/facul.h"

struct sieve_info_s;
typedef struct sieve_info_s * sieve_info_ptr;
typedef const struct sieve_info_s * sieve_info_srcptr;

#include "las-unsieve.h"

struct sieve_side_info_s {
    unsigned char Bound[256]; /* zero for good lognorms, 127 otherwise */
    fbprime_t *trialdiv_primes;
    trialdiv_divisor_t *trialdiv_data;
    unsigned char lognorm_table[1 << NORM_BITS];
    factorbase_degn_t * fb;

    /* These fields are used for the norm initialization essentially.
     * Only the scale is also relevant to part of the rest, since it
     * determines the logp contributions for factor base primes */
    double scale;      /* norm scale used on the algebraic side */
    double logmax;     /* norms on the alg-> side are < 2^alg->logmax */

    mpz_t *fij;       /* coefficients of F(a0*i+a1*j, b0*i+b1*j)  */
    double *fijd;     /* coefficients of F_q (divided by q on the special q side) */
};

typedef struct sieve_side_info_s * sieve_side_info_ptr;
typedef const struct sieve_side_info_s * sieve_side_info_srcptr;
typedef struct sieve_side_info_s sieve_side_info[1];

// General information about the siever
struct sieve_info_s {
    cado_poly cpoly;

    // general operational flags
    int nb_threads;
    FILE *output;
    const char * outputname; /* keep track of whether it's gzipped or not */
    int verbose;
    int bench;
    int ratq;   // 0 means special q on alg side, otherwise on rat side

    // sieving area
    uint32_t I;
    uint32_t J;
    int logI; // such that I = 1<<logI

    // description of the q-lattice
    uint64_t q;
    uint64_t rho;
    int32_t a0, b0, a1, b1;

    // parameters for bucket sieving
    int bucket_thresh;    // bucket sieve primes >= bucket_thresh
    int nb_buckets;
    int bucket_limit;   // maximal number of bucket_reports allowed in one bucket.
    // unsigned int degree;   /* polynomial degree */
    sieve_side_info sides[2];
    double B;         /* bound for the norm computation */

    unsieve_aux_data us;

    facul_strategy_t *strategy;

    unsigned int degree;
};

typedef struct sieve_info_s sieve_info[1];

#endif	/* LAS_TYPES_H_ */
