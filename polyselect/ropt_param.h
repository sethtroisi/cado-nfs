#ifndef ROPT_PARAM_H
#define ROPT_PARAM_H

#include "ropt_arith.h"

/* no need to change */
#define PI 3.14159265358979324
#define SUP_ALPHA 4.843
#define MAX_LINE_LENGTH 4096
#define MAX_DEGREE 6

/* possible to change */
#define L1_SIZE 12288 // ~ l1 cache.
#define SIEVEARRAY_SIZE 134217728 // 2^27 of int16_t uses about 270Mb memory.
#define MAX_SIEVEARRAY_SIZE SIEVEARRAY_SIZE
#define TUNE_SIEVEARRAY_SIZE L1_SIZE / 2
#define TOPALPHA_EACH_SIEVEARRAY 16 // for each "SIEVEARRAY_SIZE", record top N poly's alpha avalues.
#define TOPE_EACH_SUBLATTICE 8 // for sublattice, record top 8 poly's E avalues.
#define LEN_SUBLATTICE_PRIMES 10

/* for debug purpose */
#define WANT_TUNE 0
#define DEBUG_FIND_SUBLATTICE 0
#define SKIP_ROOTSIEVE_M 0
#define DEBUG 0

/* Some structs for parameters */
typedef struct {
  /* regarding find_sublattice() functions. */
  unsigned short len_e_sl;
  unsigned short tlen_e_sl;
  unsigned short *e_sl;
  unsigned short ncrts_sl;
  unsigned short len_p_rs;
  unsigned long nbest_sl;
  /* only those u < the following in sublattices().
     also choose e_sl and len_e_sl such that \prod p_i^{e_i}
     is larger than global_u_bound_rs, then we only need to
     do a line sieving 1x[-U, U]  */
  unsigned long global_w_bound_rs;
  unsigned long global_u_bound_rs;
  mpz_t global_v_bound_rs;
  /* regarding rootsieve_run() functions. */
  float sizebound_ratio_rs;
  float exp_min_alpha_rs;
  double init_lognorm;
  double lognorm_bound;
  mpz_t modulus;
} _rsparam_t;
typedef _rsparam_t rsparam_t[1];

/* Sieving data struct */
typedef struct {
  mpz_t *f;
  mpz_t *g;
} _bestpoly_t;
typedef _bestpoly_t bestpoly_t[1];

/* Parameters struct */
typedef struct {
  int w_left_bound;
  int w_length;
  unsigned short s1_num_e_sl;
  unsigned short *s1_e_sl;

  int s2_w;
  mpz_t s2_u;
  mpz_t s2_v;
  mpz_t s2_mod;
  long s2_Amax;
  long s2_Bmax;
  int s2_Vmax;

  double lognorm_bound;

  /* flag == 0, only use -w and -l
     flag == 1, use stage 1 params
     flag == 2, use stage 2 params */
  char flag;
  mpz_t n;
  int d;
} _param_t;
typedef _param_t param_t[1];

/* Sieving data struct */
typedef struct {
  mpz_t n;
  mpz_t m;
  mpz_t *f;
  mpz_t *g;
  mpz_t *fx;
  mpz_t *gx;
  mpz_t *numerator;
  double skew;
  double alpha_proj;
  int d;
} _rsstr_t;
typedef _rsstr_t rsstr_t[1];

/* Sieve bound for (A + MOD*i)*x + (B + MOD*j), mainly
   used for stage 2 */
typedef struct {
  int Wbound;
  mpz_t Umax;
  mpz_t Umin;
  mpz_t Vmax;
  mpz_t Vmin;
  long Amax;
  long Amin;
  long Bmax;
  long Bmin;
  long Blocksize;
  mpz_t A;
  mpz_t B;
  mpz_t MOD;
} _rsbound_t;
typedef _rsbound_t rsbound_t[1];


/* --- declarations --- */

extern const unsigned int primes[];

extern const unsigned char next_prime_idx[];

extern const double exp_alpha[];

void param_init ( param_t param );

void param_clear ( param_t param );

void rsparam_init ( rsparam_t rsparam,
                    rsstr_t rs,
                    param_t param );

void rsparam_setup ( rsparam_t rsparam,
                     rsstr_t rs,
                     param_t param,
                     int verbose );

void rsparam_reset_bounds ( rsparam_t rsparam,
                            rsstr_t rs,
                            param_t param,
                            int verbose );

void rsparam_free ( rsparam_t rsparam );

void bestpoly_init ( bestpoly_t bestpoly,
                     int d );
void bestpoly_free ( bestpoly_t bestpoly,
                     int d );

void rsstr_init ( rsstr_t rs );

void rsstr_setup ( rsstr_t rs );

void rsstr_free ( rsstr_t rs );

void rsbound_setup_AB_bound ( rsbound_t rsbound,
                              rsparam_t rsparam,
                              param_t param,
                              mpz_t mod,
                              int verbose );

void rsbound_setup_sublattice ( rsbound_t rsbound,
                                mpz_t sl_A,
                                mpz_t sl_B,
                                mpz_t mod );

void rsbound_init ( rsbound_t rsbound );

void rsbound_free ( rsbound_t rsbound );

void rsbound_print ( rsbound_t rsbound );

#endif
