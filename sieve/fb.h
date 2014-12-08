/*****************************************************************
*                Functions for the factor base                  *
*****************************************************************/

#ifndef FB_H
#define FB_H

#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#ifdef HAVE_SYS_MMAN_H
#include <sys/mman.h>
#else
#define MAP_FAILED ((void *) -1)
#endif
#include <gmp.h>
#include <vector>
#include "las-config.h"
#include "cado_poly.h" // for MAXDEGREE

/* Data types */

typedef unsigned int fbprime_t; /* 32 bits should be enough for everyone */
#define FBPRIME_FORMAT "u"
#define FBPRIME_MAX UINT_MAX
#define FBPRIME_BITS 32
typedef fbprime_t fbroot_t;
#define FBROOT_FORMAT "u"
typedef unsigned long largeprime_t; /* On IA32 they'll only get 32 bit
                                       large primes */
#define LARGEPRIME_FORMAT "lu"

#define FB_MAX_PARTS 4

/* If SUPPORT_LARGE_Q is defined, 64-bit redc is used in the function that
   converts roots to the p-lattice, and the redc code needs a 64-bit
   precomputed inverse. If SUPPORT_LARGE_Q is not defined, we store only a
   32-bit inverse to conserve memory. */
#if defined(SUPPORT_LARGE_Q)
typedef uint64_t redc_invp_t;
#else
typedef uint32_t redc_invp_t;
#endif

#define LOG_SCALE 1.4426950408889634 /* 1/log(2) to 17 digits, rounded to
                                        nearest. This is enough to uniquely
                                        identify the corresponding IEEE 754
                                        double precision number */

/* Base class with private copy-constructor and assigmnet operator.
   Classes which are not copy-constructible can inherit this with:
   private NonCopyable */
class NonCopyable {
 protected:
   NonCopyable() {}
   ~NonCopyable() {}
 private:
   NonCopyable(const NonCopyable&);
   NonCopyable& operator=(const NonCopyable&);
};


/* Forward declaration so fb_general_entry can use it in constructors */
template <int Nr_roots>
class fb_entry_x_roots_s;

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
  fb_general_entry(){}
  template <int Nr_roots>
  fb_general_entry (const fb_entry_x_roots_s<Nr_roots> &e);
  void parse_line (const char *line, unsigned long linenr);
  void merge (const fb_general_entry &);
  void fprint(FILE *out) const;
  bool is_simple() const;
};


class fb_general_vector: public std::vector<fb_general_entry> {
public:
  void count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
};


/* "Simple" factor base entries. We imply q=p, k=1, oldexp=0, exp=1,
   and projective=false for all roots. */
template <int Nr_roots>
class fb_entry_x_roots_s {
public:
  fbprime_t p;
  fbroot_t roots[Nr_roots];
  fb_entry_x_roots_s(){};
  /* Allow assignment-construction from general entries */
  fb_entry_x_roots_s(const fb_general_entry &e) {
    ASSERT_ALWAYS(Nr_roots == e.nr_roots);
    p = e.p;
    for (size_t i = 0; i < Nr_roots; i++)
      roots[i] = e.roots[i];
  }
  void fprint(FILE *) const;
};


template <int Nr_roots>
class fb_vector: public std::vector<fb_entry_x_roots_s<Nr_roots> > {
  public:
  void fprint(FILE *) const;
  void count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
};


class fb_slices_interface {
public:
  virtual ~fb_slices_interface(){}
  virtual void append(const fb_general_entry &) = 0;
  virtual void fprint(FILE *) const = 0;
  virtual void count_entries(size_t *nprimes, size_t *nroots, double *weight)
    const = 0;
};

/* The fb_slices class has nr_slices vectors; when appending elements with
   append(), it appends to these vectors in a round-robin fashion. */
template <int Nr_roots>
class fb_slices : public fb_slices_interface, private NonCopyable {
  fb_vector<Nr_roots> *vectors;
  size_t nr_slices, next_slice;
 public:
  fb_slices(size_t nr_slices);
  ~fb_slices();
  fb_entry_x_roots_s<Nr_roots> *get_slice(size_t slice);
  virtual void append(const fb_general_entry &);
  virtual void fprint(FILE *) const;
  virtual void count_entries(size_t *nprimes, size_t *nroots, double *weight)
    const;
};


/* A "part" is the set of all factor base primes that get sieved over a given
   bucket region size.
   E.g., when we have only 1 level of bucket sorting, then the factor base has
   2 parts: the line-sieved primes, and the bucket-sieved primes. */
class fb_part: public fb_slices_interface, private NonCopyable {
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
  fb_slices<9> *fb9_slices;
  fb_slices<10> *fb10_slices;
  fb_general_vector general_vector;
  fb_slices_interface *choose(const unsigned int n) const {
    ASSERT_ALWAYS(n <= MAXDEGREE);
    switch (n) {
      case 0: return fb0_slices;
      case 1: return fb1_slices;
      case 2: return fb2_slices;
      case 3: return fb3_slices;
      case 4: return fb4_slices;
      case 5: return fb5_slices;
      case 6: return fb6_slices;
      case 7: return fb7_slices;
      case 8: return fb8_slices;
      case 9: return fb9_slices;
      case 10: return fb10_slices;
      default: abort();
    }
  }
public:
  fb_part(size_t nr_slices);
  ~fb_part();
  void append(const fb_general_entry &);
  void fprint(FILE *) const;
  void count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
};


/* Splits the factor base for a polynomial into disjoint parts which are
   sieved over different sieve region sizes.
   For example,
   parts[0] will contain very small primes that are line-sieved,
   parts[1] contains bucket-sieved primes with 1 level of bucket sorting
   (i.e., hits get sorted into bucket regions of size 2^16)
*/
class fb_factorbase: public fb_slices_interface, private NonCopyable {
  fb_part *parts[FB_MAX_PARTS];
  fbprime_t thresholds[FB_MAX_PARTS];
  public:
  fb_factorbase(const fbprime_t *thresholds, const size_t *nr_slices);
  ~fb_factorbase();
  void read(const char * const filename);
  void make_linear (const mpz_t *poly, fbprime_t powbound, bool do_projective);
  bool mmap_fbc(const char *) {return false;};
  void dump_fbc(const char *) {return;};
  unsigned long *extract_bycost(size_t &n, fbprime_t pmax, fbprime_t td_thresh);
  size_t size() {abort(); return 0;}
  void fprint(FILE *) const;
  void append(const fb_general_entry &);
  void count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
};


unsigned char	fb_log (double, double, double);
fbprime_t       fb_pow (fbprime_t, unsigned long);
fbprime_t       fb_is_power (fbprime_t, unsigned long *);
unsigned char   fb_make_steps(fbprime_t *, fbprime_t, double);

#endif
