/* common header file for the CADO project */

#include "gmp.h"

/* default degrees for polynomial selection, entries are in digits */
#define DEFAULT_DEGREES {{70, 3}, {90, 4}, {ULONG_MAX, 5}}
#define DEFAULT_DEGREES_LENGTH 3

/* default rational/algebraic factor base bounds */
#define DEFAULT_RLIM {{70, 200000}, {90, 400000}, {ULONG_MAX, 4294967295UL}}
#define DEFAULT_RLIM_LENGTH 4
#define DEFAULT_ALIM DEFAULT_RLIM
#define DEFAULT_ALIM_LENGTH DEFAULT_RLIM_LENGTH

/* default large prime bounds */
#define DEFAULT_LPBR {{70, 23}, {90, 24}, {ULONG_MAX, 25}}
#define DEFAULT_LPBR_LENGTH 3
#define DEFAULT_LPBA DEFAULT_LPBR
#define DEFAULT_LPBA_LENGTH DEFAULT_LPBR_LENGTH

/* default factor-residual bounds */
#define DEFAULT_MFBR {{70, 35}, {90, 37}, {ULONG_MAX, 39}}
#define DEFAULT_MFBR_LENGTH 3
#define DEFAULT_MFBA DEFAULT_MFBR
#define DEFAULT_MFBA_LENGTH DEFAULT_MFBR_LENGTH

/* default lambda values */
#define DEFAULT_RLAMBDA {{70, 1.5}, {90, 1.7}, {ULONG_MAX, 1.9}}
#define DEFAULT_RLAMBDA_LENGTH 3
#define DEFAULT_ALAMBDA DEFAULT_RLAMBDA
#define DEFAULT_ALAMBDA_LENGTH DEFAULT_RLAMBDA_LENGTH

/* default sieving block lengths */
#define DEFAULT_QINT {{70, 5000}, {90, 10000}, {ULONG_MAX, 20000}}
#define DEFAULT_QINT_LENGTH 3

typedef struct
{
  mpz_t n;    /* number to factor */
  int degree; /* (optional) wanted degree */
} __cado_input_struct;

typedef __cado_input_struct cado_input[1];

typedef struct
{
  char name[256]; /* name */
  mpz_t n;        /* number to factor */
  double skew;    /* skewness */
  int degree;     /* (algebraic) degree */
  mpz_t *f;       /* algebraic coefficients */
  mpz_t *g;       /* rational coefficients */
  mpz_t m;        /* common root of f and g mod n */
  char type[256]; /* type (gnfs or snfs) */
  unsigned long rlim; /* rational  factor base bound */
  unsigned long alim; /* algebraic factor base bound */
  int lpbr;           /* rational  large prime bound is 2^lpbr */
  int lpba;           /* algebraic large prime bound is 2^lpba */
  int mfbr;           /* bound for rational  residuals is 2^mfbr */
  int mfba;           /* bound for algebraic residuals is 2^mfba */
  double rlambda;     /* lambda sieve parameter on the rational  side */
  double alambda;     /* lambda sieve parameter on the algebraic side */
  int qintsize;       /* sieve block range */
} __cado_poly_struct;

typedef __cado_poly_struct cado_poly[1];

