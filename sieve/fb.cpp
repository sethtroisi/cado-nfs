#include "cado.h"
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cctype>
#include <gmp.h>
#include "fb_new.h"
#include "mod_ul.h"
#include "verbose.h"
#include "getprime.h"
#include "gmp_aux.h"
#include "gzip.h"

#define FB_MAX_PARTS 4
#define MAXDEGREE 8

static unsigned int fb_log_2 (fbprime_t);

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

static inline redc_invp_t
compute_invq(fbprime_t q)
{
  if (q % 2 != 0) {
    ASSERT(sizeof(unsigned long) >= sizeof(redc_invp_t));
    return (redc_invp_t) (- ularith_invmod (q));
  } else {
    return 0;
  }
}

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

/* General entries are anything that needs auxiliary information:
   Prime powers, projective roots, ramified primes where exp != oldexp + 1,
   etc. They could, of course, also store the simple cases, but for those we
   use the simple struct to conserve memory and to decide algorithms (batch
   inversion, etc.) statically. */
class fb_general_entry {
  void read_roots (const char *, unsigned long);
public:
  fbprime_t q, p; /* q = p^k */
  redc_invp_t invq; /* invq = -1/q (mod 2^32), or (mod 2^64), depending on
		       the size of redc_invp_t */
  fbroot_t roots[MAXDEGREE];
  unsigned char exp[MAXDEGREE], oldexp[MAXDEGREE];
  bool projective[MAXDEGREE];
  
  /* exp and oldexp are maximal such that:
     If not projective and a == br (mod p^k), then p^exp | F(a,b)
     -"-               a == br (mod p^(k-1)), then p^oldexp | F(a,b)
     If projective and ar == b  -"- */
  
  unsigned char k, nr_roots;
  void parse_line (const char *line, unsigned long linenr);
  void merge (const fb_general_entry &);
  void fprint(FILE *out);
};


/* Read roots from a factor base file line and store them in roots.
   line must point at the first character of the first root on the line.
   linenr is used only for printing error messages in case of parsing error.
   Returns the number of roots read. */
void
fb_general_entry::read_roots (const char *lineptr, const unsigned long linenr)
{
    size_t i_roots = 0;
    unsigned long long last_t = 0;

    while (*lineptr != '\0')
    {
        if (i_roots == MAXDEGREE) {
            fprintf (stderr,
                    "# Error, too many roots for prime %" FBPRIME_FORMAT
                    " in factor base line %lu\n", p, linenr);
            exit(EXIT_FAILURE);
        }
        /* Projective roots r, i.e., ar == b (mod p), are stored as q + r in
           the factor base file; since q can be a 32-bit value, we read the
           root as a 64-bit integer first and subtract q if necessary. */
        const unsigned long long t = strtoull_const (lineptr, &lineptr, 10);
        if (i_roots > 0 && t <= last_t) {
            fprintf (stderr,
                "# Error, roots must be sorted in the fb file, line %lu\n",
                linenr);
            exit(EXIT_FAILURE);
        }

        projective[i_roots] = (t >= p);
        roots[i_roots++] = static_cast<fbroot_t>(t - ((t >= p) ? p : 0));
        if (*lineptr != '\0' && *lineptr != ',') {
            fprintf(stderr,
                    "# Incorrect format in factor base file line %lu\n",
                    linenr);
            exit(EXIT_FAILURE);
        }
        if (*lineptr == ',')
            lineptr++;
    }

    if (i_roots == 0) {
        fprintf (stderr, "# Error, no root for prime %" FBPRIME_FORMAT
                " in factor base line %lu\n", p, linenr - 1);
        exit(EXIT_FAILURE);
    }
    nr_roots = i_roots;
}

/* Parse a factor base line.
   Return 1 if the line could be parsed and was a "short version", i.e.,
   without explicit old and new exponent.
   Return 2 if the line could be parsed and was a "long version".
   Otherwise return 0. */
void
fb_general_entry::parse_line (const char * lineptr, const unsigned long linenr)
{
    q = strtoul_const (lineptr, &lineptr, 10);
    if (q == 0) {
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
    p = q;
    exp[0] = 1;
    oldexp[0] = 0;
    k = 1;
    if (longversion) {
        unsigned long k_ul;
        const fbprime_t base = fb_is_power (q, &k_ul);
        ASSERT(ulong_isprime(base != 0 ? base : q));
        /* If q is not a power, then base==0, and we use p = q */
        if (base != 0) {
            p = base;
            k = static_cast<unsigned char>(k_ul);
        }

        /* read the multiple of logp, if any */
        /* this must be of the form  q:nlogp,oldlogp: ... */
        /* if the information is not present, it means q:1,0: ... */
        exp[0] = strtoul_const (lineptr, &lineptr, 10);

        if (exp[0] == 0) {
            fprintf(stderr, "# Error in fb_read: could not parse the integer "
		    "after the colon of prime %" FBPRIME_FORMAT "\n", q);
            exit (EXIT_FAILURE);
        }
        if (*lineptr != ',') {
            fprintf(stderr,
		    "# fb_read: exp is not followed by comma on line %lu",
		    linenr);
            exit (EXIT_FAILURE);
        }
        lineptr++; /* skip comma */
        oldexp[0] = strtoul_const (lineptr, &lineptr, 10);
        if (*lineptr != ':') {
            fprintf(stderr,
		    "# fb_read: oldlogp is not followed by colon on line %lu",
		    linenr);
            exit (EXIT_FAILURE);
        }
        ASSERT (exp[0] > oldexp[0]);
        lineptr++; /* skip colon */
    }

    read_roots(lineptr, linenr);

    /* exp and oldexp are a property of a root, not of a prime (power).
       The factor base file should specify them per root, but specifies
       them per prime instead - a bit of a design bug.
       For long version lines, we thus use the exp and oldexp values for all
       roots specified in that line. */
    for (unsigned char i = 1; i < nr_roots; i++) {
      exp[i] = exp[0];
      oldexp[i] = oldexp[0];
    }
}

void
fb_general_entry::fprint(FILE *out)
{
  fprintf(out, "%" FBPRIME_FORMAT ": ", q);
  for (unsigned char i_root = 0; i_root < nr_roots; i_root++) {
    unsigned long long t = roots[i_root];
    if (projective[i_root])
      t += q;
    fprintf(out, "%llu", t);
    if (oldexp[i_root] != 0 || exp[i_root] != 1)
      fprintf(out, ":%hhu:%hhu", oldexp[i_root], exp[i_root]);
    if (i_root + 1 < nr_roots)
      fprintf(out, ",");
  }
  fprintf(out, "\n");
}

void
fb_general_entry::merge (const fb_general_entry &other)
{
  ASSERT_ALWAYS(p == other.p && q == other.q && k == other.k);
  for (unsigned char i_root = 0; i_root < other.nr_roots; i_root++) {
    ASSERT_ALWAYS(nr_roots < MAXDEGREE);
    roots[nr_roots] = other.roots[i_root];
    exp[nr_roots] = other.exp[i_root];
    oldexp[nr_roots] = other.oldexp[i_root];
    projective[nr_roots] = other.projective[nr_roots];
    nr_roots++;
  }
}

/* "Simple" factor base entries. We imply q=p, k=1, oldexp=0, exp=1,
   and projective=false for all roots. */
template <int Nr_roots>
class fb_entry_x_roots_s {
public:
  fbprime_t p;
  fbroot_t roots[Nr_roots];
  fb_entry_x_roots_s(const fb_general_entry &);
  void fprint(FILE *);
};

template <int Nr_roots>
fb_entry_x_roots_s<Nr_roots>::fb_entry_x_roots_s(const fb_general_entry &fb_cur)
{
  // ASSERT(Nr_roots == fb_cur->nr_roots);
  p = fb_cur.p;
  for (size_t i = 0; i < Nr_roots; i++)
    roots[i] = fb_cur.roots[i];
}

template <int Nr_roots>
void
fb_entry_x_roots_s<Nr_roots>::fprint(FILE *out)
{
  fprintf(out, "%" FBPRIME_FORMAT ": ", p);
  for (size_t i = 0; i < Nr_roots; i++) {
    fprintf(out, "%" FBROOT_FORMAT "%s", roots[i],
	    (i + 1 < Nr_roots) ? "," : "");
  }
  fprintf(out, "\n");
}


template <int Nr_roots>
class fb_vector: public std::vector<fb_entry_x_roots_s<Nr_roots> > {
  public:
  using std::vector<fb_entry_x_roots_s<Nr_roots> >::push_back; /* wtf C++ */

  /* Overloaded push_back() method that accepts a fb_general_entry, which is
     converted to this vector's fb_entry_x_roots_s type before adding it */
  void push_back(const fb_general_entry &fb_cur)
  {
    fb_entry_x_roots_s<Nr_roots> f = fb_cur;
    push_back(f);
  }
  void fprint(FILE *);
};

template <int Nr_roots>
void
fb_vector<Nr_roots>::fprint(FILE *out)
{
  for (size_t i = 0; i < this->size(); i++)
    (*this)[i].fprint(out);
}

typedef std::vector<fb_general_entry> fb_general_vector;

/* http://stackoverflow.com/questions/24130093/gdb-could-not-find-operator */
template class std::vector<fb_general_entry>;


/* The fb_slices class has nr_slices vectors; when appending elements with
   append(), it appends to these vectors in a round-robin fashion. */
template <int Nr_roots>
class fb_slices {
  fb_vector<Nr_roots> *vectors;
  size_t nr_slices, next_slice;
public:
  fb_slices(size_t nr_slices);
  void append(const fb_general_entry &);
  fb_entry_x_roots_s<Nr_roots> *get_slice(size_t slice);
  void fprint(FILE *);
};

template <int Nr_roots>
fb_slices<Nr_roots>::fb_slices(const size_t nr_slices)
{
  this->nr_slices = nr_slices;
  this->vectors = new fb_vector<Nr_roots>[nr_slices];
  this->next_slice = 0;
}

template <int Nr_roots>
void
fb_slices<Nr_roots>::append(const fb_general_entry &fb_cur)
{
  vectors[next_slice++].push_back(fb_cur);
  if (next_slice == nr_slices)
    next_slice = 0;
}

template <int Nr_roots>
fb_entry_x_roots_s<Nr_roots> *
fb_slices<Nr_roots>::get_slice(size_t slice)
{
  ASSERT_ALWAYS(slice < nr_slices);
  return vectors[slice];
}

template <int Nr_roots>
void
fb_slices<Nr_roots>::fprint(FILE *out)
{
  for (size_t slice = 0; slice < nr_slices; slice++) {
    fprintf (out, "#    Slice %zu:\n", slice);
    vectors[slice].fprint(out);
  }
}

/* A "part" is the set of all factor base primes that get sieved over a given
   bucket region size.
   E.g., when we have only 1 level of bucket sorting, then the factor base has
   2 parts: the line-sieved primes, and the bucket-sieved primes. */
class fb_part {
  /* How do we identify the number of slices in each array?
     Should we have
     size_t nr_entries[9];
     here, or let each array end with a "magic value", such as a slice with
     0 entries, or a NULL pointer? */

  /* These vectors are filled when we read or generate the factor base.
     The slices point into the vectors' storage. */
  fb_slices<0> *fb0_slices; /* From 0 to MAXDEGREE */
  fb_slices<1> *fb1_slices;
  fb_slices<2> *fb2_slices;
  fb_slices<3> *fb3_slices;
  fb_slices<4> *fb4_slices;
  fb_slices<5> *fb5_slices;
  fb_slices<6> *fb6_slices;
  fb_slices<7> *fb7_slices;
  fb_slices<8> *fb8_slices;
  fb_general_vector general_vector;
  public:
  fb_part(size_t nr_slices);
  void append(const fb_general_entry &);
  void fprint(FILE *);
};

fb_part::fb_part(const size_t nr_slices) {
  fb0_slices = new fb_slices<0>(nr_slices);
  fb1_slices = new fb_slices<1>(nr_slices);
  fb2_slices = new fb_slices<2>(nr_slices);
  fb3_slices = new fb_slices<3>(nr_slices);
  fb4_slices = new fb_slices<4>(nr_slices);
  fb5_slices = new fb_slices<5>(nr_slices);
  fb6_slices = new fb_slices<6>(nr_slices);
  fb7_slices = new fb_slices<7>(nr_slices);
  fb8_slices = new fb_slices<8>(nr_slices);
}


/* Append a factor base entry given in fb_cur to the
   correct vector, as determined by the number of roots.
   "Special" primes are added to the general-primes vector.
   This function is a de-multiplexer: Using "switch" to turn a run-time value
   (fb_cur->nr_roots) into a compile-time constant which can be used as a
   template specifier. */
void
fb_part::append(const fb_general_entry &fb_cur)
{
  /* Is this a simple factor base prime? */
  bool is_simple = (fb_cur.k == 1);
  for (unsigned char i = 0; i < fb_cur.nr_roots; i++) {
    is_simple &= !fb_cur.projective[i];
    is_simple &= (fb_cur.oldexp[i] == 0);
    is_simple &= (fb_cur.exp[i] == 1);
  }

  /* Non-simple ones go in the general vector */
  if (!is_simple) {
    general_vector.push_back(fb_cur);
    return;
  }

  /* Simple ones go in the simple vector with the corresponding number of
     roots */
  switch (fb_cur.nr_roots) {
    case 0: fb0_slices->append(fb_cur); break;
    case 1: fb1_slices->append(fb_cur); break;
    case 2: fb2_slices->append(fb_cur); break;
    case 3: fb3_slices->append(fb_cur); break;
    case 4: fb4_slices->append(fb_cur); break;
    case 5: fb5_slices->append(fb_cur); break;
    case 6: fb6_slices->append(fb_cur); break;
    case 7: fb7_slices->append(fb_cur); break;
    case 8: fb8_slices->append(fb_cur); break;
    default: abort();
  }
}

void
fb_part::fprint(FILE *out)
{
  fprintf(out, "#   Entries with 0 roots:\n");
  fb0_slices->fprint(out);
  fprintf(out, "#   Entries with 1 root:\n");
  fb1_slices->fprint(out);
  fprintf(out, "#   Entries with 2 roots:\n");
  fb2_slices->fprint(out);
  fprintf(out, "#   Entries with 3 roots:\n");
  fb3_slices->fprint(out);
  fprintf(out, "#   Entries with 4 roots:\n");
  fb4_slices->fprint(out);
  fprintf(out, "#   Entries with 5 roots:\n");
  fb5_slices->fprint(out);
  fprintf(out, "#   Entries with 6 roots:\n");
  fb6_slices->fprint(out);
  fprintf(out, "#   Entries with 7 roots:\n");
  fb7_slices->fprint(out);
  fprintf(out, "#   Entries with 8 roots:\n");
  fb8_slices->fprint(out);

  fprintf(out, "#   General entries (powers, ramified primes or primes with projective roots):\n");
  for (size_t i = 0; i < general_vector.size(); i++) {
    general_vector[i].fprint(out);
  }
}

/* Splits the factor base for a polynomial into disjoint parts which are
   sieved over different sieve region sizes.
   For example, 
   parts[0] will contain very small primes that are line-sieved,
   parts[1] contains bucket-sieved primes with 1 level of bucket sorting
   (i.e., hits get sorted into bucket regions of size 2^16)
*/
class fb_factorbase {
  fb_part *parts[FB_MAX_PARTS];
  void append(const fb_general_entry &);
  fbprime_t thresholds[FB_MAX_PARTS];
  public:
  fb_factorbase(const fbprime_t *thresholds, const size_t *nr_slices);
  void read(const char * const filename);
  void make_linear (const mpz_t *poly, fbprime_t powbound, bool do_projective);
  void fprint(FILE *);
};

fb_factorbase::fb_factorbase(const fbprime_t *thresholds,
			     const size_t *nr_slices)
{
  for (size_t i = 0; i < FB_MAX_PARTS; i++) {
    this->thresholds[i] = thresholds[i];
    this->parts[i] = new fb_part(nr_slices[i]);
  }
}

/* Append a factor base entry to the factor base.
   The new entry is inserted into the correct part, as determined by the
   size of the prime p, and within that part, into the correct slice, as
   determined by the number of roots. */
void
fb_factorbase::append(const fb_general_entry &fb_cur)
{
  int i;
  static bool printed_too_large_prime_warning = false;

  /* Find the smallest threshold t such that t >= q */
  for (i = 0; i < FB_MAX_PARTS && fb_cur.q > thresholds[i]; i++);
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
  parts[i]->append(fb_cur);
}

/* Remove newline, comment, and trailing space from a line. Write a
   '\0' character to the line at the position where removed part began (i.e.,
   line gets truncated).
   Return length in characters or remaining line, without trailing '\0'
   character.
*/
size_t
read_strip_comment (char *const line)
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

void
fb_factorbase::read(const char * const filename)
{
  fb_general_entry fb_cur, fb_last;
  bool had_entry = false;
  FILE *fbfile;
  // too small linesize led to a problem with rsa768;
  // it would probably be a good idea to get rid of fgets
  const size_t linesize = 1000;
  char line[linesize];
  unsigned long linenr = 0;
  fbprime_t maxprime = 0;
  unsigned long nr_primes = 0;
  
  fbfile = fopen_maybe_compressed (filename, "r");
  if (fbfile == NULL) {
    fprintf (stderr, "# Could not open file %s for reading\n", filename);
    return;
  }
  
  while (!feof(fbfile)) {
    /* Sadly, the size parameter of fgets() is of type int */
    if (fgets (line, static_cast<int>(linesize), fbfile) == NULL)
      break;
    linenr++;
    if (read_strip_comment(line) == (size_t) 0) {
      /* Skip empty/comment lines */
      continue;
    }
    
    fb_cur.parse_line (line, linenr);
    fb_cur.invq = compute_invq(fb_cur.q);
    if (!had_entry) {
      fb_last = fb_cur;
      had_entry = true;
    } else if (fb_cur.q == fb_last.q) {
      fb_last.merge(fb_cur);
    } else {
      append(fb_last);
      fb_last = fb_cur;
    }
    
    /* fb_fprint_entry (stdout, fb_cur); */
    if (fb_cur.p > maxprime)
      maxprime = fb_cur.p;
    nr_primes++;
  }

  if (had_entry)
    append(fb_last);
  
  verbose_output_print (0, 2, "# Factor base successfully read, %lu primes, "
			"largest was %" FBPRIME_FORMAT "\n",
			nr_primes, maxprime);
  
  fclose_maybe_compressed (fbfile, filename);
  
  return;
}

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
  size_t size() {return powers->size();}
};

/* Create a list of prime powers (with exponent >1) up to lim */
fb_powers::fb_powers (const fbprime_t lim)
{
  fbprime_t p;
  powers = new std::vector<fb_power_t>;
  
  for (p = 2; p <= lim / p; p = getprime(1)) {
    fbprime_t q = p;
    unsigned char k = 1;
    do {
      q *= p;
      k++;
      fb_power_t new_entry = {p, q, k};
      powers->push_back(new_entry);
    } while (q <= lim / p);
  }
  getprime(0);
  
  std::sort (powers->begin(), powers->end(), cmp_powers);
}

/* Generate a factor base with primes <= bound and prime powers <= powbound
   for a linear polynomial. If projective != 0, adds projective roots
   (for primes that divide leading coefficient).
   Returns 1 on success, 0 on error. */

void
fb_factorbase::make_linear (const mpz_t *poly, const fbprime_t powbound,
			    const bool do_projective)
{
  fbprime_t next_prime;
  fb_general_entry fb_cur;

  fb_powers *powers = new fb_powers(powbound);
  size_t next_pow = 0;

  verbose_output_vfprint(0, 1, gmp_vfprintf,
               "# Making factor base for polynomial g(x) = %Zd*x%s%Zd,\n"
               "# including primes up to %" FBPRIME_FORMAT
               " and prime powers up to %" FBPRIME_FORMAT ".\n",
               poly[1], (mpz_cmp_ui (poly[0], 0) >= 0) ? "+" : "",
               poly[0], thresholds[FB_MAX_PARTS-1], powbound);

  for (next_prime = 2; next_prime <= thresholds[FB_MAX_PARTS-1]; ) {
    /* Handle any prime powers that are smaller than next_prime */
    if (next_pow < powers->size() && (*powers)[next_pow].q <= next_prime) {
      /* The list of powers must not include primes */
      ASSERT_ALWAYS((*powers)[next_pow].q < next_prime);
      fb_cur.q = (*powers)[next_pow].q;
      fb_cur.p = (*powers)[next_pow].p;
      fb_cur.k = (*powers)[next_pow].k;
      next_pow++;
    } else {
      fb_cur.q = fb_cur.p = next_prime;
      fb_cur.k = 1;
      next_prime = getprime(1);
    }
    fb_cur.nr_roots = 1;
    fb_cur.exp[0] = fb_cur.k;
    fb_cur.oldexp[0] = fb_cur.k - 1U;

    fb_cur.projective[0] = fb_linear_root (&fb_cur.roots[0], poly, fb_cur.q);
    if (fb_cur.projective[0] && !do_projective)
      continue; /* If root is projective and we don't want those,
                   skip to next prime */

    fb_cur.invq = compute_invq(fb_cur.q);
    append(fb_cur);
  }

  getprime (0); /* free prime iterator */

  delete (powers);
}

void
fb_factorbase::fprint(FILE *out)
{
  for (size_t part = 0; part < FB_MAX_PARTS; part++) {
    fprintf(out, "#  Factor base entries up to %" FBPRIME_FORMAT "\n", thresholds[part]);
    parts[part]->fprint(out);
  }
}



int main(int argc, char **argv)
{
  fbprime_t thresholds[4] = {200, 1000, 1000, 1000};
  size_t nr_slices[4] = {3, 5, 0, 0};
  fbprime_t powbound = 100;
  fb_factorbase *fb1 = new fb_factorbase(thresholds, nr_slices),
    *fb2 = new fb_factorbase(thresholds, nr_slices);

  mpz_t poly[2];

  mpz_init(poly[0]);
  mpz_init(poly[1]);

  mpz_set_ui(poly[0], 727);
  mpz_set_ui(poly[1], 210); /* Bunch of projective primes */

  fb1->make_linear(poly, powbound, true);
  fprintf (stdout, "# Linear factor base:\n");
  fb1->fprint(stdout);
  if (argc > 1) {
    fb2->read(argv[1]);
    fprintf (stdout, "# Factor base from file:\n");
    fb2->fprint(stdout);
  }

  mpz_clear(poly[0]);
  mpz_clear(poly[1]);
  exit(EXIT_SUCCESS);
}
