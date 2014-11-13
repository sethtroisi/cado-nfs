#include "cado.h"
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cctype>
#include <gmp.h>
#include "fb.h"
#include "mod_ul.h"
#include "verbose.h"
#include "getprime.h"
#include "gmp_aux.h"
#include "gzip.h"

#define FB_MAX_PARTS 4
#define MAXDEGREE 8

static unsigned int fb_log_2 (fbprime_t);

/* "Normal" factor base entries. We imply q=p, k=1, oldexp=0, exp=1,
   and projective=false for all roots. */
template <int Nr_roots>
struct fb_entry_x_roots_s {
    fbprime_t p;
    fbroot_t roots[Nr_roots];
};

/* General entries are anything that needs auxiliary information:
   Prime powers, projective roots, ramified primes where exp != oldexp + 1,
   etc. */
struct fb_entry_general_s {
    fbprime_t q, p; /* q = p^k */
    redc_invp_t invq; /* invq = -1/q (mod 2^32), or (mod 2^64) */
    fbroot_t roots[MAXDEGREE];
    bool projective[MAXDEGREE];
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

template <int Nr_roots>
class fb_vector: public std::vector<fb_entry_x_roots_s<Nr_roots> > {
  public:
  using std::vector<fb_entry_x_roots_s<Nr_roots> >::push_back; /* wtf C++ */

  /* Overloaded push_back() method that accepts a fb_entry_general_s, which is
     converted to this vector's fb_entry_x_roots_s type before adding it */
  void push_back(const fb_entry_general_s &fb_cur)
  {
    fb_entry_x_roots_s<Nr_roots> f;
    // ASSERT(Nr_roots == fb_cur->nr_roots);
    f.p = fb_cur.p;
    for (size_t i = 0; i < Nr_roots; i++)
      f.roots[i] = fb_cur.roots[i];
     this->push_back(f);
  }
};

typedef std::vector<fb_entry_general_s> fb_general_vector;

template <int Nr_roots>
struct fb_slice_x_roots_s {
  size_t nr_entries;
  fb_entry_x_roots_s<Nr_roots> *entries;
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
  fb_entry_general_s *general;
  

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
  fb_general_vector general_vector;
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
   "Special" primes are added to the general-primes vector.
   This function is a de-multiplexer: Using "switch" to turn a run-time value
   (fb_cur->nr_roots) into a compile-time constant which can be used as a
   template specifier. */
static void
fb_part_append(fb_part_ptr part, const fb_entry_general_s &fb_cur)
{
  bool have_projective = false;
  for (unsigned char i = 0; i < fb_cur.nr_roots; i++) {
    have_projective |= (fb_cur.projective[i]);
  }
  if (fb_cur.k != 1 || have_projective || fb_cur.oldexp != 0 || fb_cur.exp != 1) {
    part->general_vector.push_back(fb_cur);
    return;
  }

  switch (fb_cur.nr_roots) {
    case 0: part->fb_0_vector.push_back(fb_cur); break;
    case 1: part->fb_1_vector.push_back(fb_cur); break;
    case 2: part->fb_2_vector.push_back(fb_cur); break;
    case 3: part->fb_3_vector.push_back(fb_cur); break;
    case 4: part->fb_4_vector.push_back(fb_cur); break;
    case 5: part->fb_5_vector.push_back(fb_cur); break;
    case 6: part->fb_6_vector.push_back(fb_cur); break;
    case 7: part->fb_7_vector.push_back(fb_cur); break;
    case 8: part->fb_8_vector.push_back(fb_cur); break;
    default: abort();
  }
}

/* Append a factor base entry to the factor base.
   The new entry is inserted into the correct part, as determined by the
   size of the prime p, and within that part, into the correct slice, as
   determined by the number of roots. */
void
fb_append(fb_factorbase_ptr fb, const fb_entry_general_s &fb_cur)
{
  int i;
  static bool printed_too_large_prime_warning = false;

  /* Find the smallest threshold t such that t >= q */
  for (i = 0; i < FB_MAX_PARTS && fb_cur.q > fb->thresholds[i]; i++);
  /* No prime > largest threshold should ever be added */
  if (i == FB_MAX_PARTS) {
    if (!printed_too_large_prime_warning) {
      verbose_output_print(1, 0, "Factor base entry %" FBPRIME_FORMAT 
                           " is above factor base bound, skipping it "
                           "(and all other too large entries)\n", fb_cur.q);
      printed_too_large_prime_warning = true;
    }
    return; /* silently skip this entry */
  }
  fb_part_append(&fb->parts[i], fb_cur);
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

static inline unsigned long long int
strtoull_const(const char *nptr, const char **endptr, const int base)
{
  char *end;
  unsigned long long r;
  r = strtoull(nptr, &end, base);
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
fb_linear_root (fbroot_t *root, const mpz_t *poly, fbprime_t q)
{
  modulusul_t m;
  residueul_t r0, r1;
  bool is_projective;


  modul_initmod_ul (m, q);
  modul_init_noset0 (r0, m);
  modul_init_noset0 (r1, m);

  modul_set_ul_reduced (r0, mpz_fdiv_ui (poly[0], q), m);
  modul_set_ul_reduced (r1, mpz_fdiv_ui (poly[1], q), m);

  /* We want poly[1] * a + poly[0] * b == 0 <=>
     a/b == - poly[0] / poly[1] */
  is_projective = (modul_inv (r1, r1, m) == 0); /* r1 = 1 / poly[1] */

  if (is_projective)
    {
      ASSERT_ALWAYS(mpz_gcd_ui(NULL, poly[1], q) > 1);
      /* Set r1 = poly[0] % q, r0 = poly[1] (mod q) */
      modul_set (r1, r0, m);
      modul_set_ul_reduced (r0, mpz_fdiv_ui (poly[1], q), m);
      int rc = modul_inv (r1, r1, m);
      ASSERT_ALWAYS(rc != 0);
    }

  modul_mul (r1, r0, r1, m); /* r1 = poly[0] / poly[1] */
  modul_neg (r1, r1, m); /* r1 = - poly[0] / poly[1] */

  *root = modul_get_ul (r1, m);

  modul_clear (r0, m);
  modul_clear (r1, m);
  modul_clearmod (m);

  return is_projective;
}

class fb_powers {
  struct fb_power_t {
    fbprime_t p, q;
    unsigned char k;
  };
  std::vector<fb_power_t> *powers;

  static int
  cmp_powers(fb_power_t a, fb_power_t b)
  {
    return a.q < b.q;
  }

  public:
  fb_powers(fbprime_t);
  ~fb_powers(){delete this->powers;};
  fb_power_t operator[] (const size_t i) {
    return (*this->powers)[i];
  }
  size_t size() {return this->powers->size();}
};

/* Create a list of prime powers (with exponent >1) up to lim */
fb_powers::fb_powers (const fbprime_t lim)
{
  fbprime_t p;
  this->powers = new std::vector<fb_power_t>;
  
  for (p = 2; p <= lim / p; p = getprime(1)) {
    fbprime_t q = p;
    unsigned char k = 1;
    do {
      q *= p;
      k++;
      fb_power_t new_entry = {p, q, k};
      this->powers->push_back(new_entry);
    } while (q <= lim / p);
  }
  getprime(0);
  
  std::sort (this->powers->begin(), this->powers->end(), cmp_powers);
}

/* Generate a factor base with primes <= bound and prime powers <= powbound
   for a linear polynomial. If projective != 0, adds projective roots
   (for primes that divide leading coefficient).
   Returns 1 on success, 0 on error. */

void
new_fb_make_linear (fb_factorbase_ptr fb, const mpz_t *poly,
                const fbprime_t bound, const fbprime_t powbound,
                const int do_projective)
{
  fbprime_t next_prime;
  fb_entry_general_s fb_cur;

  fb_powers *powers = new fb_powers(powbound);
  size_t next_pow = 0;

  verbose_output_vfprint(0, 1, gmp_vfprintf,
               "# Making factor base for polynomial g(x) = %Zd*x%s%Zd,\n"
               "# including primes up to %" FBPRIME_FORMAT
               " and prime powers up to %" FBPRIME_FORMAT ".\n",
               poly[1], (mpz_cmp_ui (poly[0], 0) >= 0) ? "+" : "",
               poly[0], bound, powbound);

  for (next_prime = 2; next_prime <= bound; ) {
    /* Handle any prime powers that are smaller than next_prime */
    if (next_pow < powers->size() && (*powers)[next_pow].p <= next_prime) {
      /* The list of powers must not include primes */
      ASSERT_ALWAYS(next_prime < (*powers)[next_pow].q);
      fb_cur.q = (*powers)[next_pow].q;
      fb_cur.p = (*powers)[next_pow].p;
      fb_cur.k = (*powers)[next_pow].k;
      next_pow++;
    } else {
      fb_cur.q = fb_cur.p = next_prime;
      fb_cur.k = 1;
      next_prime = getprime(1);
    }
    fb_cur.exp = fb_cur.k;
    fb_cur.oldexp = fb_cur.k - 1U;
    fb_cur.nr_roots = 1;

    fb_cur.projective[0] = fb_linear_root (&fb_cur.roots[0], poly, fb_cur.q);
    if (fb_cur.projective[0] && !do_projective)
      continue; /* If root is projective and we don't want those,
                   skip to next prime */
    if (fb_cur.q % 2 != 0) {
      ASSERT(sizeof(unsigned long) >= sizeof(redc_invp_t));
      fb_cur.invq = (redc_invp_t) (- ularith_invmod (fb_cur.q));
    }

    fb_append (fb, fb_cur);
  }

  getprime (0); /* free prime iterator */

  delete (powers);
}


/* READ */
/* Remove newline, comment, and trailing space from a line. Write a
   '\0' character to the line at the position where removed part began (i.e.,
   line gets truncated).
   Return length in characters or remaining line, without trailing '\0'
   character.
*/
static size_t
fb_read_strip_comment (char *const line)
{
    size_t linelen, i;

    linelen = strlen (line);
    if (linelen > 0 && line[linelen - 1] == '\n')
        linelen--; /* Remove newline */
    for (i = 0; i < linelen; i++) /* Skip comments */
        if (line[i] == '#') {
            linelen = i;
            break;
        }
    while (linelen > 0 && isspace((int)(unsigned char)line[linelen - 1]))
        linelen--; /* Skip whitespace at end of line */
    line[linelen] = '\0';

    return linelen;
}

/* Read roots from a factor base file line and store them in roots.
   line must point at the first character of the first root on the line.
   linenr and p are used only for printing error messages in case of parsing error.
   Returns the number of roots read.
*/
static size_t
fb_read_roots (fbroot_t * const roots, bool *projective, const char *lineptr,
               const unsigned long linenr, const fbprime_t p)
{
    size_t nr_roots = 0;
    unsigned long long last_t = 0;

    while (*lineptr != '\0')
    {
        if (nr_roots == MAXDEGREE) {
            fprintf (stderr,
                    "# Error, too many roots for prime %" FBPRIME_FORMAT
                    " in factor base line %lu\n", p, linenr);
            exit(EXIT_FAILURE);
        }
        /* Projective roots r, i.e., ar == b (mod p), are stored as q + r in
           the factor base file; since q can be a 32-bit value, we read the
           root as a 64-bit integer first and subtract q if necessary. */
        const unsigned long long t = strtoull_const (lineptr, &lineptr, 10);
        if (nr_roots > 0 && t <= last_t) {
            fprintf (stderr,
                "# Error, roots must be sorted in the fb file, line %lu\n",
                linenr);
            exit(EXIT_FAILURE);
        }

        projective[nr_roots] = (t >= p);
        roots[nr_roots++] = static_cast<fbroot_t>(t - ((t >= p) ? p : 0));
        if (*lineptr != '\0' && *lineptr != ',') {
            fprintf(stderr,
                    "# Incorrect format in factor base file line %lu\n",
                    linenr);
            exit(EXIT_FAILURE);
        }
        if (*lineptr == ',')
            lineptr++;
    }

    if (nr_roots == 0) {
        fprintf (stderr, "# Error, no root for prime %" FBPRIME_FORMAT
                " in factor base line %lu\n", p, linenr - 1);
        exit(EXIT_FAILURE);
    }
    return nr_roots;
}

/* Parse a factor base line.
   Return 1 if the line could be parsed and was a "short version", i.e.,
   without explicit old and new exponent.
   Return 2 if the line could be parsed and was a "long version".
   Otherwise return 0. */
static void
fb_parse_line (fb_entry_general_s *const fb_cur, const char * lineptr,
               const unsigned long linenr)
{
    fb_cur->q = strtoul_const (lineptr, &lineptr, 10);
    if (fb_cur->q == 0) {
        fprintf(stderr, "# fb_read: prime is not an integer on line %lu\n",
                linenr);
        exit (EXIT_FAILURE);
    } else if (*lineptr != ':') {
        fprintf(stderr,
                "# fb_read: prime is not followed by colon on line %lu",
                linenr);
        exit (EXIT_FAILURE);
    }

    lineptr++; /* Skip colon after q */
    const bool longversion = (strchr(lineptr, ':') != NULL);

    /* NB: a short version is not permitted for a prime power, so we
     * do the test for prime powers only for long version */
    fb_cur->p = fb_cur->q;
    fb_cur->k = fb_cur->exp = 1;
    fb_cur->oldexp = 0;
    if (longversion) {
        unsigned long k;
        const fbprime_t p = fb_is_power (fb_cur->q, &k);
        ASSERT(ulong_isprime(p != 0 ? p : fb_cur->q));
        /* If q is not a power, then p==0, and we use p = q */
        if (p != 0) {
            fb_cur->p = p;
            fb_cur->k = static_cast<unsigned char>(k);
        } 

        /* read the multiple of logp, if any */
        /* this must be of the form  q:nlogp,oldlogp: ... */
        /* if the information is not present, it means q:1,0: ... */
        fb_cur->exp = strtoul_const (lineptr, &lineptr, 10);
        /*
        if (fb_cur->exp == 0) {
            fprintf(stderr, "# Error in fb_read: could not parse the integer after the colon of prime %" FBPRIME_FORMAT "\n", q);
            exit (EXIT_FAILURE);
        }*/
        if (*lineptr != ',') {
            fprintf(stderr, "# fb_read: nlogp is not followed by comma on line %lu", linenr);
            exit (EXIT_FAILURE);
        }
        lineptr++; /* skip comma */
        fb_cur->oldexp = strtoul_const (lineptr, &lineptr, 10);
        /*
        if (fb_cur->oldexp == 0) {
            fprintf(stderr, "# Error in fb_read: could not parse the integer after the comma of prime %" FBPRIME_FORMAT "\n", q);
            exit (EXIT_FAILURE);
        }*/
        if (*lineptr != ':') {
            fprintf(stderr, "# fb_read: oldlogp is not followed by colon on line %lu", linenr);
            exit (EXIT_FAILURE);
        }
        ASSERT (fb_cur->exp > fb_cur->oldexp);
        lineptr++; /* skip colon */
    }

    /* Read roots */
    fb_cur->nr_roots = fb_read_roots(fb_cur->roots, fb_cur->projective, lineptr, linenr, fb_cur->q);
}

/* Read a factor base file, splitting it into pieces.
   
   Primes and prime powers up to smalllim go into fb_small. If smalllim is 0,
   all primes go into fb_small, and nothing is written to fb_pieces.
   
   If smalllim is not 0, then nr_pieces separate factor bases are made for
   primes/powers > smalllim; factor base entries from the file are written to 
   these pieces in round-robin manner.

   Pointers to the allocated memory of the factor bases are written to fb_small 
   and, if smalllim > 0, to fb_pieces[0, ..., nr_pieces-1].

   Returns 1 if everything worked, and 0 if not (i.e., if the file could not be 
   opened, or memory allocation failed)
*/

int 
new_fb_read (fb_factorbase_ptr fb, const char * const filename)
{
    fb_entry_general_s fb_cur;
    FILE *fbfile;
    // too small linesize led to a problem with rsa768;
    // it would probably be a good idea to get rid of fgets
    const size_t linesize = 1000;
    char line[linesize];
    unsigned long linenr = 0;
    fbprime_t maxprime = 0;
    unsigned long nr_primes = 0;
    int error = 0;

    fbfile = fopen_maybe_compressed (filename, "r");
    if (fbfile == NULL) {
        fprintf (stderr, "# Could not open file %s for reading\n", filename);
        return 0;
    }

    while (!feof(fbfile)) {
        /* Sadly, the size parameter of fgets() is of type int */
        if (fgets (line, static_cast<int>(linesize), fbfile) == NULL)
            break;
        linenr++;
        if (fb_read_strip_comment(line) == (size_t) 0) {
            /* Skip empty/comment lines */
            continue;
        }

        fb_parse_line (&fb_cur, line, linenr);

        /* Compute invp */
        if (fb_cur.q % 2 != 0) {
            ASSERT(sizeof(unsigned long) >= sizeof(redc_invp_t));
            fb_cur.invq = 
                (redc_invp_t) (- ularith_invmod ((unsigned long) fb_cur.q));
        }

        fb_append (fb, fb_cur);

        /* fb_fprint_entry (stdout, fb_cur); */
	if (fb_cur.p > maxprime)
	    maxprime = fb_cur.p;
        nr_primes++;
    }

    if (!error) {
        verbose_output_print (0, 2, "# Factor base successfully read, %lu primes, largest was %"
                FBPRIME_FORMAT "\n", nr_primes, maxprime);
    }

    fclose_maybe_compressed (fbfile, filename);

    return error ? 0 : 1;
}




int main(int argc, char **argv)
{
  fb_factorbase_t fb1, fb2;

  fb1.thresholds[0] = 200;
  fb1.thresholds[1] = 1000;
  fb1.thresholds[2] = 0;
  fb1.thresholds[3] = 0;

  fb2.thresholds[0] = 200;
  fb2.thresholds[1] = 1000;
  fb2.thresholds[2] = 0;
  fb2.thresholds[3] = 0;

  mpz_t poly[2];
  
  mpz_init(poly[0]);
  mpz_init(poly[1]);
  
  mpz_set_ui(poly[0], 727);
  mpz_set_ui(poly[1], 210); /* Bunch of projective primes */

  new_fb_make_linear(&fb1, poly, 1000, 100, 1);
  if (argc > 1)
    new_fb_read(&fb2, argv[1]);

  mpz_clear(poly[0]);
  mpz_clear(poly[1]);
  exit(EXIT_SUCCESS);
}
