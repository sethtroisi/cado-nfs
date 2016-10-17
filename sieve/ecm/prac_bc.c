#include "cado.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include <pthread.h>

#include "portability.h"
#include "macros.h"
#include "getprime.h"
#include "prac_bc.h"

#if 0
#define PP1_DICT_NRENTRIES 6
static size_t pp1_dict_len[PP1_DICT_NRENTRIES] = {1, 1, 2, 2, 3, 4};
static literal_t *pp1_dict_entry[PP1_DICT_NRENTRIES] =
  {"\xB", "\xA", "\xB\xA", "\x3\x0", "\x3\xB\xA", "\x3\x0\x3\x0"};
static code_t pp1_dict_code[PP1_DICT_NRENTRIES] = {0, 0, 10, 11, 13, 14};

static bc_dict_t pp1_dict =
  {PP1_DICT_NRENTRIES, pp1_dict_len, pp1_dict_entry, pp1_dict_code};


#endif

/* This dict is used to compress bytecode:
 *  - replace "fi" (final add of subchain / subchain init) into one code (10)
 *  - replace "\x3s" (rule 3 followed by rule 's' (swap)) into one code (11)
 *  - replace "\x3fi" into one code (12) [ here we have to encode f as \x66,
 *    else the compiler interpret \x3f as 1 char ]
 *  - replace "\x3s\x3s" into one code (13)
 */
#define PRAC_DICT_NRENTRIES 4
static size_t prac_dict_len[PRAC_DICT_NRENTRIES] = {2, 2, 3, 4};
static literal_t *prac_dict_entry[PRAC_DICT_NRENTRIES] =
                                  { "fi", "\x3s", "\x3\x66i", "\x3s\x3s"};
static code_t prac_dict_code[PRAC_DICT_NRENTRIES] = {10, 11, 12, 13};

static bc_dict_t prac_dict =
          {PRAC_DICT_NRENTRIES, prac_dict_len, prac_dict_entry, prac_dict_code};

/* Table of multipliers for PRAC. prac_mul[i], 1<=<i<=9, has continued 
   fraction sequence of all ones but with a 2 in the (i+1)-st place, 
   prac_mul[i], 10<=<i<=17, has continued fraction sequence of all ones 
   but with a 2 in second and in the (i-7)-th place
   and prac_mul[0] is all ones, i.e. the golden ratio. */
#define PRAC_NR_MULTIPLIERS 10
static const double prac_mul[PRAC_NR_MULTIPLIERS] = 
  {1.61803398874989484820 /* 0 */, 1.38196601125010515179 /* 1 */, 
   1.72360679774997896964 /* 2 */, 1.58017872829546410471 /* 3 */, 
   1.63283980608870628543 /* 4 */, 1.61242994950949500192 /* 5 */,
   1.62018198080741576482 /* 6 */, 1.61721461653440386266 /* 7 */, 
   1.61834711965622805798 /* 8 */, 1.61791440652881789386 /* 9 */
#if 0
   ,
   1.41982127170453589529, 1.36716019391129371457,
   1.38757005049050499808, 1.37981801919258423518,
   1.38278538346559613734, 1.38165288034377194202,
   1.38208559347118210614, 1.38192033153010418805,
   1.70431400059211079746};
#else
   };
#endif

/************************* Cache mechanism ************************************/
/* Mutual exclusion lock for the global variables for the cache mechanism */
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
struct prac_cache_s
{
  unsigned int k;
  const prac_cost_t *opcost;
  double best_mul;
  double mincost;
};

typedef struct prac_cache_s prac_cache_t;

prac_cache_t *prac_cache = NULL;
unsigned int prac_cache_nb = 0;
unsigned int prac_cache_alloc = 0;

static int
prac_cache_is_in (unsigned int k, const prac_cost_t *opcost, double * best_mul,
                  double * mincost)
{
  pthread_mutex_lock (&lock);
  for (unsigned int i = 0; i < prac_cache_nb; i++)
    if (prac_cache[i].k == k && prac_cache[i].opcost == opcost)
    {
      *best_mul = prac_cache[i].best_mul;
      *mincost = prac_cache[i].mincost;
      pthread_mutex_unlock (&lock);
      return 1;
    }
  pthread_mutex_unlock (&lock);
  return 0;
}

static void
prac_cache_add_one (unsigned int k, const prac_cost_t *opcost, double best_mul,
                    double mincost)
{
  pthread_mutex_lock (&lock);
  if (prac_cache_nb >= prac_cache_alloc)
  {
    prac_cache_alloc += 16;
    prac_cache = (prac_cache_t *)
                 realloc (prac_cache, prac_cache_alloc * sizeof (prac_cache_t));
    ASSERT_ALWAYS (prac_cache != NULL);
  }

  prac_cache[prac_cache_nb++] = (prac_cache_t) { .k = k, .opcost= opcost,
                                     .best_mul = best_mul, .mincost = mincost };
  pthread_mutex_unlock (&lock);
}

void
prac_cache_free ()
{
  pthread_mutex_lock (&lock);
  if (prac_cache)
    free (prac_cache);
  prac_cache = NULL;
  prac_cache_alloc = 0;
  prac_cache_nb = 0;
  pthread_mutex_unlock (&lock);
}

/***********************************************************************
   Generating Lucas chains with Montgomery's PRAC algorithm. Code taken 
   from GMP-ECM, mostly written by Paul Zimmermann, and slightly 
   modified here
************************************************************************/


/* Produce a PRAC chain with initial multiplier v. 
   Returns its arithmetic cost or DBL_MAX if no chain could be computed.
   The cost of a differential addition is opcost->dadd, the cost a doubling is
   opcost->dbl. */

static double
prac_chain (const unsigned int n, const double v, const prac_cost_t *opcost,
            bc_state_t *state)
{
  unsigned int d, e, r;
  double cost;

  d = n;
  r = (unsigned int) ((double) d / v + 0.5);
  if (r >= n)
    return DBL_MAX;
  d = n - r;
  e = 2 * r - n;

  /* initial doubling (subchain init) */
  cost = opcost->dbl;
  if (state != NULL)
    bytecoder ((literal_t) 'i', state);

  while (d != e)
  {
    if (d < e)
    { /* swap d and e */
      r = d;
      d = e;
      e = r;
      if (state != NULL)
        bytecoder ((literal_t) 's', state);
    }
    if (4 * d <= 5 * e && ((d + e) % 3) == 0)
    { /* condition 1 */
      d = (2 * d - e) / 3;
      e = (e - d) / 2;
      cost += 3 * opcost->dadd; /* 3 additions */
      if (state != NULL)
        bytecoder ((literal_t) 1, state);
    }
    else if (4 * d <= 5 * e && (d - e) % 6 == 0)
    { /* condition 2 */
      d = (d - e) / 2;
      cost += opcost->dadd + opcost->dbl; /* one addition, one doubling */
      if (state != NULL)
        bytecoder ((literal_t) 2, state);
    }
    else if (d <= 4 * e)
    { /* condition 3 */
      d -= e;
      cost += opcost->dadd; /* one addition */
      if (state != NULL)
        bytecoder ((literal_t) 3, state);
    }
    else if ((d + e) % 2 == 0)
    { /* condition 4 */
      d = (d - e) / 2;
      cost += opcost->dadd + opcost->dbl; /* one addition, one doubling */
      if (state != NULL)
        bytecoder ((literal_t) 4, state);
    }
    /* now d+e is odd */
    else if (d % 2 == 0)
    { /* condition 5 */
      d /= 2;
      cost += opcost->dadd + opcost->dbl; /* one addition, one doubling */
      if (state != NULL)
        bytecoder ((literal_t) 5, state);
    }
    /* now d is odd and e even */
    else if (d % 3 == 0)
    { /* condition 6 */
      d = d / 3 - e;
      cost += 3 * opcost->dadd + opcost->dbl; /* three additions, one doubling */
      if (state != NULL)
        bytecoder ((literal_t) 6, state);
    }
    else if ((d + e) % 3 == 0)
    { /* condition 7 */
      d = (d - 2 * e) / 3;
      cost += 3 * opcost->dadd + opcost->dbl; /* three additions, one doubling */
      if (state != NULL)
        bytecoder ((literal_t) 7, state);
    }
    else if ((d - e) % 3 == 0)
    { /* condition 8 */
      d = (d - e) / 3;
      cost += 3 * opcost->dadd + opcost->dbl; /* three additions, one doubling */
      if (state != NULL)
        bytecoder ((literal_t) 8, state);
    }
    else /* necessarily e is even */
    { /* condition 9 */
      e /= 2;
      cost += opcost->dadd + opcost->dbl; /* one addition, one doubling */
      if (state != NULL)
        bytecoder ((literal_t) 9, state);
    }
  }

  /* final addition */
  cost += opcost->dadd;
  if (state != NULL)
    bytecoder ((literal_t) 'f', state);

  /* Here d = gcd(n, r). If d != 1, then the sequence cannot be
     reduced below d, i.e., the chain does not start at 1.
     It would be the end of a concatenated chain instead.
     Return DBL_MAX in this case. */
  return (d == 1) ? cost : DBL_MAX;
}

/* Write bytecode for an addition chain for odd k (repeated val times), and
 * return its cost.
 * k must be odd and prime (not exactly: prac algorithm always succeeds on
 * primes and sometimes fails on composite; it must exist at least one
 * multipliers for which the caller knonws that prac will succeeds; for k prime
 * it is always the case).
 */
static double
prac_bytecode_one (const unsigned int k, const unsigned int val,
                   const prac_cost_t *opcost, bc_state_t *state, int verbose)
{
  double mincost = DBL_MAX, best_mul = 0.;

  ASSERT_ALWAYS (k % 2 == 1);

  if (k == 3)
  {
    /* There is only one Lucas chain possible in this case so we do not care
     * what multipliers we use.
     */
     best_mul = prac_mul[0];
     mincost = opcost->dbl + opcost->dadd;

    if (verbose)
      printf ("## PRAC: k=3 best_mul=%f mincost=%f\n", best_mul, mincost);
  }
  else if (prac_cache_is_in (k, opcost, &best_mul, &mincost))
  {
    if (verbose)
      printf ("## PRAC: k=%u best_mul=%f mincost=%f [from cache]\n", k,
              best_mul, mincost);
  }
  else
  {
    /* Find the best multiplier for this k */
    for (unsigned int i = 0 ; i < PRAC_NR_MULTIPLIERS; i++)
    {
      double mul = prac_mul[i];
      if (verbose > 1)
        printf ("### Compute chain with PRAC for k=%u with mul=%f\n", k, mul);
      double cost = prac_chain (k, mul, opcost, NULL);

      if (verbose && cost == DBL_MAX)
        printf ("## PRAC: could not compute Lucas chain with mul=%f\n", mul);
      if (cost < mincost)
      {
        mincost = cost;
        best_mul = mul;
      }
      if (verbose > 1)
        printf ("## PRAC: k=%u mul=%f cost=%f\n", k, mul, cost);
    }

    if (verbose)
      printf ("## PRAC: k=%u best_mul=%f mincost=%f\n", k, best_mul, mincost);

    /* add it to the cache */
    prac_cache_add_one (k, opcost, best_mul, mincost);
  }

  /* If we have mincost == DBL_MAX, it means that this k is composite and prac
   * cannot make a valid chain for it. We could try to factor k and make a
   * composite chain from the prime factors. For now, we bail out - caller
   * shouldn't give a composite k that doesn't have a prac chain we can find.
   */
  ASSERT_ALWAYS (mincost != DBL_MAX);

  /* Write the best chain (repeated val times) and return the cost of the chain
   * (repeated val times)
   */
  double cost = 0.;
  for (unsigned int i = 0; i < val; i++)
    cost += prac_chain (k, best_mul, opcost, state);
  return cost;
}

unsigned int
prac_bytecode (char **bc, unsigned int B1, unsigned int pow2_nb,
               unsigned int pow3_extra, const prac_cost_t *opcost, int compress,
               int verbose)
{
  mpz_t E; /* only use when verbose > 0 for checking the chain */
  double prac_cost = 0.;
  bc_state_t *bc_state = bytecoder_init ((compress) ? &prac_dict : NULL);

  if (verbose)
  {
    mpz_init (E);
    if (pow3_extra) /* we start with additional power of 3 */
      mpz_ui_pow_ui (E, 3, pow3_extra);
    else
      mpz_set_ui (E, 1);
  }
  /* Do all the odd primes */
  prime_info pi;
  prime_info_init (pi);
  unsigned int p = (unsigned int) getprime_mt (pi);
  ASSERT_ALWAYS (p == 3);
  for ( ; p <= B1; p = (unsigned int) getprime_mt (pi))
  {
    unsigned int val = (p == 3) ? pow3_extra : 0;
    for (unsigned int q = 1; q <= B1 / p; q *= p)
    {
      val++;
      if (verbose)
        mpz_mul_ui (E, E, p);
    }

    prac_cost += prac_bytecode_one (p, val, opcost, bc_state, verbose);
  }
  prime_info_clear (pi);

  /* Write the bytecode in *bc */
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
    prac_bytecode_fprintf (stdout, *bc, bc_len);
    /* Check the bytecode */
    prac_bytecode_check (*bc, bc_len, E, verbose);
    printf ("## PRAC: cost of power of 2: %f\n", power2cost);
    printf ("## PRAC: total cost: %f\n", prac_cost + power2cost);
  }

  if (verbose)
    mpz_clear (E);
  return bc_len;
}

#define mpz_dadd(r, A, B, C) do { \
    mpz_add (r, A, B);            \
  } while (0)
void
prac_bytecode_fprintf (FILE *out, const char *bc, unsigned int len)
{
  printf (" (len=%u)", len);
  for (unsigned i = 0; i < len; i++)
  {
    char b = bc[i];
    fprintf (out, "%s %u", (i == 0) ? "" : ",", (unsigned int) b);
    switch (b)
    {
      case 's': /* Swap A, B */
        fprintf (out, " [swap]");
        break;
      case 'i':
        fprintf (out, " [i]");
        break;
      case 'f':
        fprintf (out, " [f]");
        break;
      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      case 7:
      case 8:
      case 9:
        break;
      case 10:
        fprintf (out, " [fi]");
        break;
      case 11:
        fprintf (out, " [3s]");
        break;
      case 12:
        fprintf (out, " [3fi]");
        break;
      case 13:
        fprintf (out, " [3s3s]");
        break;
      default:
        fprintf (out, " [unknown]");
    }
  }
  printf ("\n");
}

/* Return nonzero if the bytecode does not produce E. Return 0 otherwise. */
int
prac_bytecode_check (const char *bc, unsigned int len, mpz_srcptr E,
                     int verbose)
{
  if (bc[0] != 'i')
  {
    printf ("## PRAC: bytecode_check failed: must start with 'i'\n");
    return 1;
  }

  mpz_t A, B, C, t, t2;
  mpz_init_set_ui (A, 1);
  mpz_init (B);
  mpz_init (C);
  mpz_init (t);
  mpz_init (t2);

  int ret = 0;
  for (unsigned i = 0; i < len; i++)
  {
    char b = bc[i];
    switch (b)
    {
      case 's': /* Swap A, B */
        mpz_swap (A, B);
        break;
      case 'i': /* Start of a subchain */
        mpz_set (B, A);
        mpz_set (C, A);
        mpz_mul_2exp (A, A, 1);
        break;
      case 'f': /* End of a subchain */
        mpz_dadd (A, A, B, C);
        break;
      case 1:
        mpz_dadd (t, A, B, C);
        mpz_dadd (t2, t, A, B);
        mpz_dadd (B, B, t, A);
        mpz_set (A, t2);
        break;
      case 2:
        mpz_dadd (B, A, B, C);
        mpz_mul_2exp (A, A, 1);
        break;
      case 3:
        mpz_dadd (C, B, A, C);
        mpz_swap (B, C);
        break;
      case 4:
        mpz_dadd (B, B, A, C);
        mpz_mul_2exp (A, A, 1);
        break;
      case 5:
        mpz_dadd (C, C, A, B);
        mpz_mul_2exp (A, A, 1);
        break;
      case 6:
        mpz_mul_2exp (t, A, 1);
        mpz_dadd (t2, A, B, C);
        mpz_dadd (A, t, A, A);
        mpz_dadd (C, t, t2, C);
        mpz_swap (B, C);
        break;
      case 7:
        mpz_dadd (t, A, B, C);
        mpz_dadd (B, t, A, B);
        mpz_mul_2exp (t, A, 1);
        mpz_dadd (A, A, t, A);
        break;
      case 8:
        mpz_dadd (t, A, B, C);
        mpz_dadd (C, C, A, B);
        mpz_swap (B, t);
        mpz_mul_2exp (t, A, 1);
        mpz_dadd (A, A, t, A);
        break;
      case 9:
        mpz_dadd (C, C, B, A);
        mpz_mul_2exp (B, B, 1);
        break;
      case 10:
        /* Combined final add of old subchain and init of new subchain [=fi] */
        mpz_dadd (B, A, B, C);
        mpz_set (C, B);
        mpz_mul_2exp (A, B, 1);
        break;
      case 11:
        /* Combined rule 3 and rule 0 [=\x3s] */
        mpz_dadd (C, B, A, C);
        /* (B,C,A) := (A,B,C)  */
        mpz_swap (B, C);
        mpz_swap (A, B);
        break;
      case 12:
        /* Combined rule 3, then subchain end/start [=\x3fi] */
        mpz_dadd (t, B, A, C);
        mpz_dadd (C, A, t, B);
        mpz_set (B, C);
        mpz_mul_2exp (A, C, 1);
        break;
      case 13:
        /* Combined rule 3, swap, rule 3 and swap, merged a bit [=\x3s\x3s] */
        mpz_set (t, B);
        mpz_dadd (B, B, A, C);
        mpz_set (C, A);
        mpz_dadd (A, A, B, t);
        break;
      default:
        printf ("## PRAC: bytecode_check failed: Unknown bytecode %u\n",
                                                              (unsigned int) b);
        ret = 1;
    }
  }

  if (ret == 0)
  {
    if (mpz_cmp (A, E) != 0)
    {
      gmp_printf ("## PRAC: bytecode_check failed:\n"
                  "  Expected %Zd\n  Got %Zd\n", E, A);
      ret = 1;
    }
  }
  if (ret == 0 && verbose)
    printf ("## PRAC: bytecode_check: chain is ok\n");
  mpz_clear (A);
  mpz_clear (B);
  mpz_clear (C);
  mpz_clear (t);
  mpz_clear (t2);
  return ret;
}
