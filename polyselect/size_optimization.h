#ifndef POLYSELECT_SIZE_OPTIMIZATION_H_
#define POLYSELECT_SIZE_OPTIMIZATION_H_

/* Default maximal number of steps in size_optimize_local_descent */
#define SOPT_DEFAULT_MAX_STEPS 300

/* Default value for sopt_effort */
#define SOPT_DEFAULT_EFFORT 0

/* values of q greater than 1e10 in absolute value do not help */
#define SOPT_MAX_VALUE_FOR_Q_ROOTS 1e10

/* Number of best rational approximation to keep for each q-root */
#define SOPT_NB_RAT_APPROX_OF_Q_ROOTS 8

/* Maximun value for the denominator of rational approximation of q-roots */
#define SOPT_MAX_DEN_IN_RAT_APPROX_OF_Q_ROOTS 100.0

/* Call LLL for skew in [ skew0^(e/(2*NSKEW)) for e in [NSKEW..3*NSKEW] ]*/
#define SOPT_NSKEW 2 /* 2 seems to be quite good */
#define SOPT_NB_OF_SKEWNESS_VALUES (2*SOPT_NSKEW + 1)

/* Only the SOPT_MAX_LLL_POLY_PROCESS best polynomials produced by LLL are
 * given to the local optimization algorithm (should be >= 1) */
#define SOPT_MAX_LLL_POLY_PROCESS 2

/* Maximum degree possible of the rotations in the local descent algorithm. */
#define SOPT_MAX_DEGREE_ROTATION 7

#define SOPT_LOCAL_DESCENT_GUARD 0.001


double sopt_local_descent (mpz_poly_ptr, mpz_poly_ptr, mpz_poly_srcptr,
                           mpz_poly_srcptr, int, int, unsigned int, int);
double size_optimization (mpz_poly_ptr, mpz_poly_ptr, mpz_poly_srcptr,
                          mpz_poly_srcptr, const unsigned int, const int);

#endif	/* POLYSELECT_SIZE_OPTIMIZATION_H_ */


