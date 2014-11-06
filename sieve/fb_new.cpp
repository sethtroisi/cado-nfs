
#ifndef FB_STRUCT_H
#define FB_STRUCT_H

#include "cado.h"
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <gmp.h>
#include "fb.h"
#include "mod_ul.h"
#include "verbose.h"
#include "getprime.h"

#define FB_MAX_PARTS 4
#define FB_MAX_ROOTS 8

static unsigned int fb_log_2 (fbprime_t);

/* "Normal" factor base entries. We imply q=p, k=1, oldexp=0, exp=1,
   projective=false for all roots. */
template <int nr_roots>
struct fb_entry_x_roots_s {
    fbprime_t p;
    fbroot_t r[nr_roots];
};

/* "Special" entries are anything that needs auxiliary information:
   Prime powers, ramified primes where exp != oldexp + 1, etc. */
struct fb_entry_special_s {
    fbprime_t q, p; /* q = p^k */
    fbroot_t r[FB_MAX_ROOTS];
    bool projective[FB_MAX_ROOTS];
    /* exp and oldexp are maximal such that:
       If not projective and a == br (mod p^k), then p^exp | F(a,b)
           -"-               a == br (mod p^(k-1)), then p^oldexp | F(a,b)
       If projective and ar == b  -"- */
    unsigned char k, nr_roots, exp, oldexp;
};

/* If q = p^k, with k maximal and k > 1, return q.
   Otherwise return 0. If final_k is not NULL, write k there. */
fbprime_t
fb_is_power (fbprime_t q, unsigned long *final_k)
{
  unsigned long maxk, k;
  uint32_t p;

  maxk = fb_log_2(q);
  for (k = maxk; k >= 2; k--)
    {
      double dp = pow ((double) q, 1.0 / (double) k);
      double rdp = trunc(dp + 0.5);
      if (fabs(dp - rdp) < 0.001) {
        p = (uint32_t) rdp ;
        if (q % p == 0) {
          // ASSERT (fb_pow (p, k) == q);
          if (final_k != NULL)
            *final_k = k;
          return p;
        }
      }
    }
  return 0;
}

template <int nr_roots>
class fb_vector: public std::vector<fb_entry_x_roots_s<nr_roots> > {
  public:
  using std::vector<fb_entry_x_roots_s<nr_roots> >::push_back; /* wtf C++ */

  /* Overloaded push_back() method that accepts a factorbase_degn_t, which is
     converted to this vector's fb_entry_x_roots_s type before adding it */
  void push_back(const factorbase_degn_t *fb_cur)
  {
    fb_entry_x_roots_s<nr_roots> f;
    // ASSERT(nr_roots == fb_cur->nr_roots);
    f.p = fb_cur->p;
    for (size_t i = 0; i < nr_roots; i++)
      f.r[i] = fb_cur->roots[i];
     this->push_back(f);
  }
};

class fb_special_vector: public std::vector<fb_entry_special_s> {
  public:
  using std::vector<fb_entry_special_s>::push_back;

  void push_back(const factorbase_degn_t *fb_cur, const bool *projective)
  {
    fb_entry_special_s f;
    f.p = fb_cur->p;
    f.q = fb_pow(fb_cur->p, fb_cur->k);
    for (size_t i = 0; i < f.nr_roots; i++) {
      f.r[i] = fb_cur->roots[i];
      f.projective[i] = projective[i];
    }
    f.k = fb_cur->k;
    f.nr_roots = fb_cur->nr_roots;
    f.exp = fb_cur->exp;
    f.oldexp = fb_cur->oldexp;
    this->push_back(f);
  }
};

template <int nr_roots>
struct fb_slice_x_roots_s {
  size_t nr_entries;
  fb_entry_x_roots_s<nr_roots> *entries;
};

/* A "part" is the set of all factor base primes that get sieved over a given
   bucket region size.
   E.g., when we have only 1 level of bucket sorting, then the factor base has
   2 parts: the line-sieved primes, and the bucket-sieved primes. */
struct fb_part_s {
  /* How do we identify the number of slices in each array?
     Should we have
     size_t nr_entries[9];
     here, or let each array end with a "magic value", such as a slice with
     0 entries, or a NULL pointer? */
  fb_slice_x_roots_s<0> *fb0;
  fb_slice_x_roots_s<1> *fb1;
  fb_slice_x_roots_s<2> *fb2;
  fb_slice_x_roots_s<3> *fb3;
  fb_slice_x_roots_s<4> *fb4;
  fb_slice_x_roots_s<5> *fb5;
  fb_slice_x_roots_s<6> *fb6;
  fb_slice_x_roots_s<7> *fb7;
  fb_slice_x_roots_s<8> *fb8;
  fb_entry_special_s *special;
  

  /* These vectors are filled when we read or generate the factor base.
     The slices point into the vectors' storage. */
  fb_vector<0> fb_0_vector;
  fb_vector<1> fb_1_vector;
  fb_vector<2> fb_2_vector;
  fb_vector<3> fb_3_vector;
  fb_vector<4> fb_4_vector;
  fb_vector<5> fb_5_vector;
  fb_vector<6> fb_6_vector;
  fb_vector<7> fb_7_vector;
  fb_vector<8> fb_8_vector;
  fb_special_vector special_vector;
};
typedef struct fb_part_s * fb_part_ptr;

/* Splits the factor base for a polynomial into disjoint parts which are
   sieved over different sieve region sizes.
   For example, 
   parts[0] will contain very small primes that are line-sieved,
   parts[1] contains bucket-sieved primes with 1 level of bucket sorting
   (i.e., hits get sorted into bucket regions of size 2^16)
*/
struct fb_factorbase_s {
  fb_part_s parts[FB_MAX_PARTS];
  fbprime_t thresholds[FB_MAX_PARTS];
};
typedef struct fb_factorbase_s fb_factorbase_t;
typedef struct fb_factorbase_s *fb_factorbase_ptr;


/* Return p^e. Trivial exponentiation for small e, no check for overflow */
fbprime_t
fb_pow (const fbprime_t p, const unsigned long e)
{
    fbprime_t r = 1;

    for (unsigned long i = 0; i < e; i++)
      r *= p;
    return r;
}

/* Returns floor(log_2(n)) for n > 0, and 0 for n == 0 */
static unsigned int
fb_log_2 (fbprime_t n)
{
  unsigned int k;
  for (k = 0; n > 1; n /= 2, k++);
  return k;
}

/* Append a factor base entry given in fb_cur to the
   correct vector, as determined by the number of roots.
   "Special" primes are added to the special-primes vector.
   This function is a de-multiplexer: Using "switch" to turn a run-time value
   (fb_cur->nr_roots) into a compile-time constant which can be used as a
   template specifier. */
static void
fb_part_append(fb_part_ptr part, const factorbase_degn_t *fb_cur, const bool *projective)
{
  if (fb_cur->k != 1 || projective != NULL || fb_cur->oldexp != 0 || fb_cur->exp != 1) {
    part->special_vector.push_back(fb_cur, projective);
    return;
  }

  switch (fb_cur->nr_roots) {
    case 0: part->fb_0_vector.push_back(fb_cur);
    case 1: part->fb_1_vector.push_back(fb_cur);
    case 2: part->fb_2_vector.push_back(fb_cur);
    case 3: part->fb_3_vector.push_back(fb_cur);
    case 4: part->fb_4_vector.push_back(fb_cur);
    case 5: part->fb_5_vector.push_back(fb_cur);
    case 6: part->fb_6_vector.push_back(fb_cur);
    case 7: part->fb_7_vector.push_back(fb_cur);
    case 8: part->fb_8_vector.push_back(fb_cur);
    default: abort();
  }
}

/* Append a factor base entry to the factor base.
   The new entry is inserted into the correct part, as determined by the
   size of the prime p, and within that part, into the correct slice, as
   determined by the number of roots. */
void
fb_append(fb_factorbase_ptr fb, const factorbase_degn_t *fb_cur, const bool *projective)
{
  for (int i = 0; i < FB_MAX_PARTS; i++) {
    if (fb_cur->p <= fb->thresholds[i]) {
      fb_part_append(&fb->parts[i], fb_cur, projective);
    }
  }
}

/* strtoul(), but with const char ** for second argument.
   Otherwise it's not possible to do, e.g., strtoul(p, &p, 10) when p is
   of type const char *
*/
static inline unsigned long int
strtoul_const(const char *nptr, const char **endptr, const int base)
{
  char *end;
  unsigned long r;
  r = strtoul(nptr, &end, base);
  *endptr = end;
  return r;
}

/* Make one factor base entry for a linear polynomial poly[1] * x + poly[0]
   and the prime (power) q. We assume that poly[0] and poly[1] are coprime.
   Non-projective roots a/b such that poly[1] * a + poly[0] * b == 0 (mod q)
   with gcd(poly[1], q) = 1 are stored as a/b mod q.
   If do_projective != 0, also stores projective roots with gcd(q, f_1) > 1,
   but stores the reciprocal root.
   Returns true if the roots was projective, and false otherwise. */

static bool
fb_make_linear1 (factorbase_degn_t *fb_entry, const mpz_t *poly,
		 const fbprime_t p, const unsigned char newexp,
		 const unsigned char oldexp, const unsigned char k)
{
  modulusul_t m;
  residueul_t r0, r1;
  bool is_projective;

  fb_entry->p = p;
  fb_entry->exp = newexp;
  fb_entry->oldexp = oldexp;
  fb_entry->k = k;
  fb_entry->nr_roots = 1;

  modul_initmod_ul (m, p);
  modul_init_noset0 (r0, m);
  modul_init_noset0 (r1, m);

  modul_set_ul_reduced (r0, mpz_fdiv_ui (poly[0], p), m);
  modul_set_ul_reduced (r1, mpz_fdiv_ui (poly[1], p), m);

  /* We want poly[1] * a + poly[0] * b == 0 <=>
     a/b == - poly[0] / poly[1] */
  is_projective = (modul_inv (r1, r1, m) == 0); /* r1 = 1 / poly[1] */

  if (is_projective)
    {
      ASSERT_ALWAYS(mpz_gcd_ui(NULL, poly[1], p) > 1);
      /* Set r1 = poly[0] % p, r0 = poly[1] (mod p) */
      modul_set (r1, r0, m);
      modul_set_ul_reduced (r0, mpz_fdiv_ui (poly[1], p), m);
      int rc = modul_inv (r1, r1, m);
      ASSERT_ALWAYS(rc != 0);
    }

  modul_mul (r1, r0, r1, m); /* r1 = poly[0] / poly[1] */
  modul_neg (r1, r1, m); /* r1 = - poly[0] / poly[1] */

  fb_entry->roots[0] = modul_get_ul (r1, m);
  if (p % 2 != 0) {
    ASSERT(sizeof(unsigned long) >= sizeof(redc_invp_t));
    fb_entry->invp = (redc_invp_t) (- ularith_invmod (modul_getmod_ul (m)));
  }

  modul_clear (r0, m);
  modul_clear (r1, m);
  modul_clearmod (m);

  return is_projective;
}

struct fb_power_t {
  fbprime_t p, q;
  unsigned char k;
};
typedef std::vector<fb_power_t> powers_vector_t;

static int
cmp_powers(fb_power_t a, fb_power_t b)
{
  return a.q < b.q;
}

/* Create a list of prime powers (with exponent >1) up to lim */
static powers_vector_t *
make_power_list (const fbprime_t lim)
{
  fbprime_t p;
  powers_vector_t *v = new std::vector<fb_power_t>;
  
  for (p = 2; p <= lim / p; p = getprime(1)) {
    fbprime_t q = p;
    unsigned char k = 1;
    do {
      q *= p;
      k++;
      fb_power_t new_entry = {p, q, k};
      v->push_back(new_entry);
    } while (p <= lim / q);
  }
  getprime(0);
  
  std::sort (v->begin(), v->end(), cmp_powers);
  
  return v;
}


/* Generate a factor base with primes <= bound and prime powers <= powbound
   for a linear polynomial. If projective != 0, adds projective roots
   (for primes that divide leading coefficient).
   Returns 1 on success, 0 on error. */

void
fb_make_linear (fb_factorbase_ptr fb, const mpz_t *poly,
                const fbprime_t bound, const fbprime_t powbound,
                const int do_projective)
{
  fbprime_t next_prime;
  factorbase_degn_t *fb_cur;
  unsigned char newexp, oldexp;
  unsigned char k;

  fb_cur = (factorbase_degn_t *) malloc (fb_entrysize_uc (1));
  ASSERT (fb_cur != NULL);

  powers_vector_t *powers = make_power_list(powbound);
  size_t next_pow = 0;

  verbose_output_vfprint(0, 1, gmp_vfprintf,
               "# Making factor base for polynomial g(x) = %Zd*x%s%Zd,\n"
               "# including primes up to %" FBPRIME_FORMAT
               " and prime powers up to %" FBPRIME_FORMAT ".\n",
               poly[1], (mpz_cmp_ui (poly[0], 0) >= 0) ? "+" : "",
               poly[0], bound, powbound);

  for (next_prime = 2; next_prime <= bound; ) {
    fbprime_t q;
    /* Handle any prime powers that are smaller than next_prime */
    if (next_pow < powers->size() && (*powers)[next_pow].p <= next_prime) {
      /* The list of powers must not include primes */
      ASSERT_ALWAYS((*powers)[next_pow].p < next_prime);
      q = (*powers)[next_pow].q;
      k = (*powers)[next_pow].k;
      next_pow++;
    } else {
      q = next_prime;
      k = 1;
      next_prime = getprime(1);
    }
    newexp = k;
    oldexp = k - 1U;

    bool is_projective = fb_make_linear1 (fb_cur, poly, q, newexp, oldexp, k);
    if (is_projective && !do_projective)
      continue; /* If root is projective and we don't want those,
                   skip to next prime */

    fb_append (fb, fb_cur, &is_projective);
  }

  getprime (0); /* free prime iterator */

  free (fb_cur);
  free (powers);
}

#endif
