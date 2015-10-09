/**
 * @file ropt_param.c
 * Holds extern vars. The customizable parameters 
 * are in ropt_param.h.
 */


#include "cado.h"
#include "ropt_param.h"


/**
 * Default values for L1 cache size and tune sieve length. The code
 * will also try to auto-detect them.
 */
unsigned int L1_cachesize = 12288;
unsigned int size_tune_sievearray = 6144;


/* ------------------
   possible to change
   ------------------ */


/**
 * Total number of sublattices in the tuning and sieving steps.
 * If with the default parametes (SIZE_SIEVEARRAY_V_MAX 4194304),
 * each root sieve takes about 2-4 seconds. The first column is 
 * ranked by the digits of integers to be factored and the right
 * column is the number (actually number+1).
 *
 * Usually, there is one or more tuning steps before the final root
 * sieve. In that case, more sublattices are checked (e.g. double/quad
 * the following valeus).
 * 
 * The number is linearly scaled by param->effort, but a larger value
 * is not always necessary since the sieving is done in order with
 * the best sublattices first.
 */
const unsigned int size_total_sublattices[8][2] = {
  /* {digits, num_of_sublattices} */
  {80,  4},  /* for up to 79 digits */
  {100, 8},  /* up to 99 digits */
  {140, 16},  /* up to 139 digits */
  {170, 32}, /* up to 169 digits */
  {180, 64}, /* up to 179 digits */
  {220, 128}, /* up to 219 digits */
  {260, 256}, /* up to 259 digits */
  {300, 512} /* up to 299 digits */
};


/**
 * Number of top sublattice for individual primes[i] in stage 1,
 * where i < NUM_SUBLATTICE_PRIMES. The constrcution should depends
 * on the total num of primes in s1param->e_sl[]. 
 * The main purpose is to prevent too much crt computations in stage 1.
 * They will be passed to s1param->individual_nbest_sl[] later.
 */
const unsigned int
s1_size_each_sublattice[NUM_SUBLATTICE_PRIMES][NUM_SUBLATTICE_PRIMES] = {
  { 64,  0,  0,  0,  0,  0,  0,  0,  0 }, // tlen_e_sl = 1
  { 64, 64,  0,  0,  0,  0,  0,  0,  0 }, // tlen_e_sl = 2
  { 64, 64, 64,  0,  0,  0,  0,  0,  0 }, // tlen_e_sl = 3
  { 64, 32, 32, 32,  0,  0,  0,  0,  0 }, // tlen_e_sl = 4
  { 32, 32, 16,  8,  8,  0,  0,  0,  0 }, // tlen_e_sl = 5
  { 32, 16, 16,  8,  8,  4,  0,  0,  0 }, // tlen_e_sl = 6 (current bound)
  { 32, 16,  8,  8,  4,  4,  4,  0,  0 }, // tlen_e_sl = 7
  { 32, 16,  8,  8,  4,  4,  4,  2,  0 }, // tlen_e_sl = 8
  { 32, 16,  8,  8,  4,  4,  4,  2,  2 }, // tlen_e_sl = 9
};


/**
 * As above, but it's only used for tuning good w in 'ropt_quadratic.c'.
 * Therefore, the values are much smaller than above.
 */
const unsigned int
s1_size_each_sublattice_tune[NUM_SUBLATTICE_PRIMES] = {
  8,  8,  4,  4,  4,  2,  2,  2,  2
};


/* -------------------------
   perhaps no need to change
   ------------------------- */


/**
 * Default parameters for sublattice p_i^e^i.
 * Non-decreasing due to ropt_s1param_setup() in ropt_str.c
 */
const unsigned int 
default_sublattice_pe[NUM_DEFAULT_SUBLATTICE][NUM_SUBLATTICE_PRIMES] = {
  /* 2, 3, 5, 7, 11, 13, 17, 19, 23 */
  { 1, 0, 0, 0, 0, 0, 0, 0, 0 }, // 2
  { 2, 0, 0, 0, 0, 0, 0, 0, 0 }, // 4
  { 2, 1, 0, 0, 0, 0, 0, 0, 0 }, // 12
  { 3, 1, 0, 0, 0, 0, 0, 0, 0 }, // 24
  { 4, 1, 0, 0, 0, 0, 0, 0, 0 }, // 48
  { 3, 2, 0, 0, 0, 0, 0, 0, 0 }, // 72
  { 3, 1, 1, 0, 0, 0, 0, 0, 0 }, // 120
  { 4, 2, 0, 0, 0, 0, 0, 0, 0 }, // 144
  { 4, 1, 1, 0, 0, 0, 0, 0, 0 }, // 240
  { 5, 2, 0, 0, 0, 0, 0, 0, 0 }, // 288
  { 3, 2, 1, 0, 0, 0, 0, 0, 0 }, // 360
  { 5, 1, 1, 0, 0, 0, 0, 0, 0 }, // 480
  { 4, 2, 1, 0, 0, 0, 0, 0, 0 }, // 720
  { 5, 2, 1, 0, 0, 0, 0, 0, 0 }, // 1440
  { 5, 3, 1, 0, 0, 0, 0, 0, 0 }, // 4320
  { 5, 2, 1, 1, 0, 0, 0, 0, 0 }, // 10080
  { 6, 2, 1, 1, 0, 0, 0, 0, 0 }, // 20160
  { 5, 3, 2, 0, 0, 0, 0, 0, 0 }, // 21600
  { 5, 2, 2, 1, 0, 0, 0, 0, 0 }, // 50400
  { 6, 2, 2, 1, 0, 0, 0, 0, 0 }, // 100800
  { 6, 3, 2, 1, 0, 0, 0, 0, 0 }, // 302400
  { 6, 3, 2, 1, 1, 0, 0, 0, 0 }, // 3326400
  { 7, 3, 2, 1, 1, 0, 0, 0, 0 }, // 6652800
  { 8, 3, 2, 1, 1, 0, 0, 0, 0 }, // 13305600
  { 7, 4, 2, 1, 1, 0, 0, 0, 0 }, // 19958400
  { 6, 3, 2, 1, 1, 1, 0, 0, 0 }, // 43243200
};


/**
 * Hard-coded product for above default_sublattice_pe.
 */
const unsigned long
default_sublattice_prod[NUM_DEFAULT_SUBLATTICE] = {
  2, 4, 12, 24, 48, 72, 120, 144, 240, 288, 360, 480,
  720, 1440, 4320, 10080, 20160, 21600, 50400, 100800,
  302400, 3326400, 6652800, 13305600, 19958400, 43243200,
};


/**
 * Primes.
 */
const unsigned int primes[NP] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199
};


/**
 * next_prime_idx[3] = ind(5) = 2 in above table.
 */
const unsigned char next_prime_idx[] = {
  0, // 0
  0, 1, 2, 2, 3, 3, 4, 4, 4, 4, // 1 - 10
  5, 5, 6, 6, 6, 6, 7, 7, 8, 8,
  8, 8, 9, 9, 9, 9, 9, 9, 10, 10,
  11, 11, 11, 11, 11, 11, 12, 12, 12, 12,
  13, 13, 14, 14, 14, 14, 15, 15, 15, 15,
  15, 15, 16, 16, 16, 16, 16, 16, 17, 17,
  18, 18, 18, 18, 18, 18, 19, 19, 19, 19,
  20, 20, 21, 21, 21, 21, 21, 21, 22, 22,
  22, 22, 23, 23, 23, 23, 23, 23, 24, 24,
  24, 24, 24, 24, 24, 24, 25, 25, 25, 25,
  26, 26, 27, 27, 27, 27, 28, 28, 29, 29,
  29, 29, 30, 30, 30, 30, 30, 30, 30, 30,
  30, 30, 30, 30, 30, 30, 31, 31, 31, 31,
  32, 32, 32, 32, 32, 32, 33, 33, 34, 34,
  34, 34, 34, 34, 34, 34, 34, 34, 35, 35,
  36, 36, 36, 36, 36, 36, 37, 37, 37, 37,
  37, 37, 38, 38, 38, 38, 39, 39, 39, 39,
  39, 39, 40, 40, 40, 40, 40, 40, 41, 41,
  42, 42, 42, 42, 42, 42, 42, 42, 42, 42,
  43, 43, 44, 44, 44, 44, 45, 45 // 191 - 198
};


/**
 * Asymptotic estimate of minimum order statistics
 * for 2^K many rotations where 0 <= K <= 149.
 * function expected_alpha_est() in alpha.sage
 */
const double exp_alpha[] = {
  0.000,   -0.617,   -0.951,   -1.254,   -1.521,
  -1.759,   -1.976,   -2.176,   -2.362,   -2.536,
  -2.701,   -2.858,   -3.008,   -3.151,   -3.289,
  -3.422,   -3.550,   -3.674,   -3.794,   -3.911,
  -4.025,   -4.136,   -4.245,   -4.350,   -4.454,
  -4.555,   -4.654,   -4.751,   -4.847,   -4.940,
  -5.032,   -5.122,   -5.211,   -5.299,   -5.385,
  -5.470,   -5.553,   -5.636,   -5.717,   -5.797,
  -5.876,   -5.954,   -6.031,   -6.107,   -6.183,
  -6.257,   -6.330,   -6.403,   -6.475,   -6.546,
  -6.617,   -6.686,   -6.755,   -6.824,   -6.891,
  -6.958,   -7.025,   -7.091,   -7.156,   -7.220,
  -7.284,   -7.348,   -7.411,   -7.473,   -7.535,
  -7.597,   -7.658,   -7.718,   -7.778,   -7.838,
  -7.897,   -7.956,   -8.014,   -8.072,   -8.130,
  -8.187,   -8.244,   -8.300,   -8.356,   -8.411,
  -8.467,   -8.522,   -8.576,   -8.630,   -8.684,
  -8.738,   -8.791,   -8.844,   -8.896,   -8.949,
  -9.001,   -9.052,   -9.104,   -9.155,   -9.206,
  -9.256,   -9.307,   -9.357,   -9.407,   -9.456,
  -9.505,   -9.554,   -9.603,   -9.652,   -9.700,
  -9.748,   -9.796,   -9.843,   -9.891,   -9.938,
  -9.985,  -10.032,  -10.078,  -10.124,  -10.170,
  -10.216,  -10.262,  -10.307,  -10.353,  -10.398,
  -10.443,  -10.487,  -10.532,  -10.576,  -10.620,
  -10.664,  -10.708,  -10.752,  -10.795,  -10.838,
  -10.881,  -10.924,  -10.967,  -11.010,  -11.052,
  -11.094,  -11.136,  -11.178,  -11.220,  -11.262,
  -11.303,  -11.345,  -11.386,  -11.427,  -11.468,
  -11.509,  -11.549,  -11.590,  -11.630,  -11.670
};
