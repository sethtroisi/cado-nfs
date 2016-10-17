#include "cado.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include <gmp.h>

#include "portability.h"
#include "macros.h"
#include "addchain_bc.h"
#include "getprime.h"
#include "gmp_aux.h" /* for nbits */

#define ADDCHAIN_DICT_NRENTRIES 1
static size_t addchain_dict_len[ADDCHAIN_DICT_NRENTRIES] = {2};
static literal_t *addchain_dict_entry[ADDCHAIN_DICT_NRENTRIES] =
                                           {ADDCHAIN_DBL_STR ADDCHAIN_DBL_STR };
static code_t addchain_dict_code[ADDCHAIN_DICT_NRENTRIES] = {ADDCHAIN_2DBL};
static bc_dict_t addchain_dict =
  {ADDCHAIN_DICT_NRENTRIES, addchain_dict_len, addchain_dict_entry, addchain_dict_code};

/* Only used with l <= 31 */
#define mpz_mod_ui_2exp(n,l) (((n)->_mp_size) ? ((n)->_mp_d[0]&((1<<l)-1)) : 0)

/* compute a <- (a - r)/2 */
#define mpz_sub_div2_si(a, r)                     \
  do {                                            \
    if (r >= 0)                                   \
      mpz_sub_ui (a, a, r);                       \
    else                                          \
      mpz_add_ui (a, a, -r);                      \
    mpz_fdiv_q_2exp (a, a, 1); /* divides by 2 */ \
  } while (0)

/* Given k and q compute the best value for r, i.e. the odd integer r in [-q..q]
 * with the biggest valuation in 2 for k-r.
 *  Let l = nbits (2*q)
 *  r = k % 2^l ( with remainder in [-2^(l-1)..2^(l-1)-1]
 *  if abs(r) <= q
 *    return r
 *  else
 *    recompute r with the same formula but with l-1
 *    we known that in this case abs(r) <= q
 *    return r
 */
int
addchain_find_best_r (mpz_srcptr k, const uint8_t q)
{
  int l = nbits(2*q);
  unsigned int r = mpz_mod_ui_2exp (k, l);
  if (r <= q)
    return (int) r;
  else if ((1 << l) - r <= q)
    return (int) (r - (1 << l));
  else
  {
    l--;
    r = mpz_mod_ui_2exp (k, l);
    if (r <= q)
      return (int) r;
    else /* we are sure that ((1 << l) - r <= q) */
      return (int) (r - (1 << l));
  }
}

#define ZEROPOINT ((uint8_t) 0x00)

/* If r = ZEROPOINT the starting point is kP. Else starting point is kP + rP
 * r must be != 2
 */
static inline double
_starting_point (uint8_t k, uint8_t r, const addchain_cost_t * opcost,
                 bc_state_t *state, int verbose)
{
  if (verbose > 1)
  {
    if (r != ZEROPOINT)
      printf ("  Q <- %uP + %uP [cost: %f]\n", k, r, opcost->add);
    else
      printf ("  Q <- %uP [cost: 0]\n", k);
  }
  ASSERT_ALWAYS (r != 2);
  if (state)
  {
    literal_t bc1 = (k == 2) ? 0x7f : (k >> 1);
    literal_t bc2 = (r == ZEROPOINT) ? 0x00 : (r >> 1);
    if (r != ZEROPOINT)
      bc1 |= 0x80; /* set the first bit to 1 */
    bytecoder (bc1, state);
    bytecoder (bc2, state);
  }
  return (r == ZEROPOINT) ? 0. : opcost->add;
}

/* r must be odd */
static inline double
_dbladd (int r, const addchain_cost_t * opcost, bc_state_t *state, int verbose)
{
  if (verbose > 1)
    printf ("  Q <- 2Q%+dP [cost: %f]\n", r, opcost->dbladd);

  uint8_t rabs = (r >= 0) ? (uint8_t) r : (uint8_t) -r;
  ASSERT_ALWAYS (rabs % 2 == 1);
  if (state)
  {
    literal_t bc = (rabs >> 1);
    if (r < 0)
      bc |= 0x80; /* set the first bit to 1 */
    bytecoder (bc, state);
  }

  return opcost->dbladd;
}

static inline double
_dbl (const addchain_cost_t * opcost, bc_state_t *state, int verbose)
{
  if (verbose > 1)
    printf ("  Q <- 2Q [cost: %f]\n", opcost->dbl);

  if (state)
    bytecoder (ADDCHAIN_DBL, state);

  return opcost->dbl;
}

static double
addchain_rec (mpz_t k, const uint8_t q, const addchain_cost_t * opcost,
              bc_state_t *state, int verbose)
{
  /* Case 0: k is <= q and (2 or odd) [ i.e., one of the precomputed points ] */
  if (mpz_cmp_ui (k, q) <= 0 && (mpz_cmp_ui (k, 2) == 0 || mpz_odd_p (k)))
  {
    /* starting point is kP */
    if (verbose > 1)
      gmp_printf ("# k = %Zd [base case 0]\n", k);
    return _starting_point (mpz_get_ui (k), ZEROPOINT, opcost, state, verbose);
  }
  /* Case 1: k == q + 2 */
  else if (mpz_cmp_ui (k, q+2) == 0)
  {
    /* starting point is qP + 2P */
    if (verbose > 1)
      gmp_printf ("# k = %Zd [base case 1]\n", k);
    return _starting_point (2, q, opcost, state, verbose);
  }
  /* Case 2:  q+4 <= k <= 3*q-2 and k % 6 == 1 */
  else if (0 <= mpz_cmp_ui (k, q+4) && mpz_cmp_ui (k, 3*q-2) <= 0
                                    && mpz_congruent_ui_p (k, 1, 6))
  {
    if (verbose > 1)
      gmp_printf ("# k = %Zd [base case 2]\n", k);
    unsigned int kui = mpz_get_ui (k);
    uint8_t s = (uint8_t) ((kui+2)/3);
    uint8_t r = (uint8_t) ((kui-4)/3);
    /* starting point is sP */
    double c = _starting_point (s, ZEROPOINT, opcost, state, verbose);
    /* Then we do add dbladd */
    return c + _dbladd (r, opcost, state, verbose);
  }
  /* Case 3:  q+4 <= k <= 3*q and k % 6 == 3 */
  else if (0 <= mpz_cmp_ui (k, q+4) && mpz_cmp_ui (k, 3*q) <= 0
                                    && mpz_congruent_ui_p (k, 3, 6))
  {
    if (verbose > 1)
      gmp_printf ("# k = %Zd [base case 3]\n", k);
    unsigned int kui = mpz_get_ui (k);
    uint8_t s = (uint8_t) (kui/3);
    /* starting point is sP */
    double c = _starting_point (s, ZEROPOINT, opcost, state, verbose);
    /* Then we do add dbladd */
    return c + _dbladd (s, opcost, state, verbose);
    /* XXX If one day, we have a tripling more efficient that add + dbladd, we
     * could use it here */
  }
  /* Case 4:  q+4 <= k <= 3*q-4 and k % 6 == 5
   * Note: for q = 1, 3*q-4 is < 0 so this cas never happens.
   */
  else if (0 <= mpz_cmp_ui (k, q+4) && q > 1 && mpz_cmp_ui (k, 3*q-4) <= 0
                                    && mpz_congruent_ui_p (k, 5, 6))
  {
    if (verbose > 1)
      gmp_printf ("# k = %Zd [base case 4]\n", k);
    unsigned int kui = mpz_get_ui (k);
    uint8_t s = (uint8_t) ((kui-2)/3);
    uint8_t r = (uint8_t) ((kui+4)/3);
    /* starting point is sP */
    double c = _starting_point (s, ZEROPOINT, opcost, state, verbose);
    /* Then we do add dbladd */
    return c + _dbladd (r, opcost, state, verbose);
  }
  /* Case 5: 4 <= k <= 2*q-2 and k % 4 == 0 */
  else if (0 <= mpz_cmp_ui (k, 4) && mpz_cmp_ui (k, 2*q-2) <= 0
                                  && mpz_mod_ui_2exp (k, 2) == 0)
  {
    if (verbose > 1)
      gmp_printf ("# k = %Zd [base case 5]\n", k);
    unsigned int kui = mpz_get_ui (k);
    uint8_t r1 = (uint8_t) ((kui-2)/2);
    uint8_t r2 = (uint8_t) ((kui+2)/2);
    /* starting point is r1P + r2P */
    return _starting_point (r1, r2, opcost, state, verbose);
  }
  /* Even case */
  else if (mpz_even_p (k))
  {
    if (verbose > 1)
      gmp_printf ("# k = %Zd [case even]\n", k);
    mpz_fdiv_q_2exp (k, k, 1); /* divides by 2 */
    double cost = addchain_rec (k, q, opcost, state, verbose);
    return cost + _dbl (opcost, state, verbose);
  }
  else
  {
    if (verbose > 1)
      gmp_printf ("# k = %Zd [odd case]\n", k);
    int r = addchain_find_best_r (k, q);
    mpz_sub_div2_si (k, r);
    double cost = addchain_rec (k, q, opcost, state, verbose);
    return cost + _dbladd (r, opcost, state, verbose);
  }
}

double
addchain (mpz_srcptr E, const uint8_t q, const addchain_cost_t * opcost,
          bc_state_t *state, int verbose)
{
  double cost;
  mpz_t k;

  FATAL_ERROR_CHECK (q > ADDCHAIN_Q_MAX,
                                "q must be smaller or equal to ADDCHAIN_Q_MAX");
  FATAL_ERROR_CHECK (q % 2 == 0, "q must be odd");

  mpz_init_set (k, E);

  /* We put q as the first byte in the bytecode */
  if (state)
    bytecoder ((literal_t) q, state);

  /* Cost of the precomputation: 1 DBL (to compute 2P) and (q-1)/2 ADD (to
   * compute (2k+1)P from (2k-1)P and 2P, for k in [1..(q-1)/2])
   */
  cost = opcost->dbl_precomp + (q >> 1) * opcost->add_precomp;
  if (verbose > 1)
    printf ("## cost of precomputation: %f\n", cost);

  /* Recursively compute the cost of the addition chain */
  cost += addchain_rec (k, q, opcost, state, verbose);

  mpz_clear (k);
  return cost;
}

/* Bytecode an addition chain for
 *        E = prod (p^floor(log(B1)/log(p)) for f odd prime <= B1)
 */
unsigned int
addchain_bytecode (char **bc, unsigned int B1, unsigned int pow2_nb,
                   unsigned int pow3_extra, const addchain_cost_t * opcost,
                   int compress, int verbose)
{
  mpz_t E;

  /* E = product of all primes between 3 and B1 */
  mpz_init (E);
  if (pow3_extra) /* we start with additional power of 3 */
    mpz_ui_pow_ui (E, 3, pow3_extra);
  else
    mpz_set_ui (E, 1);
  prime_info pi;
  prime_info_init (pi);
  unsigned long p = getprime_mt (pi);
  ASSERT (p == 3);

  for ( ; p <= B1; p = (unsigned int) getprime_mt (pi))
  {
    for (unsigned long q = 1; q <= B1 / p; q *= p)
      mpz_mul_ui (E, E, p);
  }
  prime_info_clear (pi);

  if (verbose > 1)
    gmp_printf ("### E = %Zd\n", E);
  if (verbose)
    gmp_printf ("### E has %zu bits\n", mpz_sizeinbase (E, 2));

  /* Computes the bytecode for E: try every odd q up to ADDCHAIN_M_MAX and keep
   * the one with the lowest cost.
   */
  double mincost = DBL_MAX;
  uint8_t best_q = 0;
  for (uint8_t q = 1 ; q <= ADDCHAIN_Q_MAX ; q += 2)
  {
    if (verbose > 1)
      printf ("### Compute addchain with q = %u\n", q);
    double cost = addchain (E, q, opcost, NULL, verbose);
    if (cost < mincost)
    {
      mincost = cost;
      best_q = q;
    }
    if (verbose)
      printf ("## Addchain: q = %u ; cost = %f\n", q, cost);
  }

  if (verbose)
    printf ("## Addchain: best_q = %u mincost = %f\n", best_q, mincost);

  /* Put the best addchain into bytecode */
  bc_state_t * bc_state = bytecoder_init ((compress) ? &addchain_dict : NULL);
  double chaincost = addchain (E, best_q, opcost, bc_state, verbose);
  bytecoder_flush (bc_state);
  unsigned int bc_len = bytecoder_size (bc_state);
  *bc = (char *) malloc (bc_len * sizeof (char));
  ASSERT_ALWAYS (*bc);
  bytecoder_read (*bc, bc_state);
  bytecoder_clear (bc_state);

  if (verbose)
  {
    /* The cost of the initial doublings */
    double power2cost = pow2_nb * opcost->dbl;
    /* Print the bytecode */
    printf ("Byte code for stage 1: ");
    addchain_bytecode_fprintf (stdout, *bc, bc_len);
    /* Check the bytecode */
    addchain_bytecode_check (*bc, bc_len, E, verbose);
    printf ("## Addchain: cost of power of 2: %f\n", power2cost);
    printf ("## Addchain: total cost: %f\n", chaincost + power2cost);
  }

  mpz_clear (E);
  return bc_len;
}

void
addchain_bytecode_fprintf (FILE *out, const char *bc, unsigned int len)
{
  ASSERT_ALWAYS (len >= 3);
  fprintf (out, "(len=%u) %x [m = %u]", len, bc[0], bc[0]);
  fprintf (out, ", %x, %x", bc[1], bc[2]);
  unsigned int k = ((bc[1] & 0x7f) == 0x7f) ? 2 :
                                            ((uint8_t) (bc[1] & 0x7f) << 1) + 1;
  if (bc[1] & 0x80)
  {
    unsigned int k2 = ((uint8_t) (bc[2] & 0x7f) << 1) + 1;
    fprintf (out, " [starting point is %uP + %uP ]", k, k2);
  }
  else
    fprintf (out, " [starting point is %uP]", k);
  for (unsigned int i = 3; i < len; i++)
  {
    char b = bc[i];
    int r;
    switch (b)
    {
      case ADDCHAIN_2DBL:
        fprintf (out, ", %x [2DBL]", (uint8_t) b);
        break;
      case ADDCHAIN_DBL:
        fprintf (out, ", %x [DBL]", (uint8_t) b);
        break;
      default:
        r = ((uint8_t) (b & 0x7f) << 1) + 1;
        if (b & 0x80)
          r = -r;
        fprintf (out, ", %x [DBLADD with r = %d]", (uint8_t) b, r);
        break;
    }
  }
  fprintf (out, "\n");
}

/* Return nonzero if the bytecode does not produce E. Return 0 otherwise. */
int
addchain_bytecode_check (const char *bc, unsigned int len, mpz_srcptr E,
                         int verbose)
{
  mpz_t Q;
  uint8_t q;

  if (len < 3)
  {
    printf ("## Addchain: bytecode_check failed: len must be at least 3\n");
    return 1;
  }

  q = (uint8_t) bc[0];
  if (q > ADDCHAIN_Q_MAX)
  {
    printf ("## Addchain: bytecode_check failed: q = %u is larger than "
            "ADDCHAIN_Q_MAX = %u\n", q, ADDCHAIN_Q_MAX);
    return 1;
  }

  unsigned int k1 = ((bc[1] & 0x7f) == 0x7f) ? 2 :
                                            ((uint8_t) (bc[1] & 0x7f) << 1) + 1;
  if (k1 > q && k1 != 2)
  {
    printf ("## Addchain: bytecode_check failed: bc[1] = %u is larger than "
            "q = %u\n", k1, q);
    return 1;
  }

  if (bc[1] & 0x80)
  {
    unsigned int k2 = ((uint8_t) (bc[2] & 0x7f) << 1) + 1;
    if (k2 > q)
    {
      printf ("## Addchain: bytecode_check failed: bc[2] = %u is larger than "
              "q = %u\n", k1, q);
      return 1;
    }
    mpz_init_set_ui (Q, k1+k2);
  }
  else
    mpz_init_set_ui (Q, k1);

  int ret = 0;
  for (unsigned int i = 3; i < len && ret == 0; i++)
  {
    char b = bc[i];
    unsigned int r;
    switch (b)
    {
      case ADDCHAIN_2DBL:
        mpz_mul_2exp (Q, Q, 2);
        break;
      case ADDCHAIN_DBL:
        mpz_mul_2exp (Q, Q, 1);
        break;
      default:
        mpz_mul_2exp (Q, Q, 1);
        r = ((uint8_t) (b & 0x7f) << 1) + 1;
        if (r > q)
        {
          printf ("## Addchain: bytecode_check failed: bc[%u] corresponds to a "
                  "dbladd with |r|=%u which is larger than q = %u\n", i, r, q);
          ret = 1;
        }
        if (b & 0x80)
          mpz_sub_ui (Q, Q, r);
        else
          mpz_add_ui (Q, Q, r);
        break;
    }
  }

  if (ret == 0)
  {
    if (mpz_cmp (Q, E) != 0)
    {
      gmp_printf ("## Addchain: bytecode_check failed:\n"
                  "  Expected %Zd\n  Got %Zd\n", E, Q);
      ret = 1;
    }
  }
  if (ret == 0 && verbose)
    printf ("## Addchain: bytecode_check: chain is ok\n");
  mpz_clear (Q);
  return ret;
}
