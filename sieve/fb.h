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
#include "fb-types.h"
#include "las-qlattice.h"

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
class fb_entry_x_roots;

/* A root modulo a prime power q. q is specified externally */
struct fb_general_root {
  fbroot_t r;
  /* exp and oldexp are maximal such that:
     If not projective and a == br (mod p^k), then p^exp | F(a,b)
     -"-               a == br (mod p^(k-1)), then p^oldexp | F(a,b)
     If projective and ar == b  -"- */
  unsigned char exp, oldexp;
  bool proj;

  fb_general_root (){}
  /* Constructor for simple entries */
  fb_general_root (const fbroot_t r) : r(r), exp(1), oldexp(0), proj(false){} 
  fb_general_root (const fbroot_t r, const unsigned char nexp,
                   const unsigned char oldexp, const bool proj) :
                   r(r), exp(nexp), oldexp(oldexp), proj(proj) {}

  /* A root is simple if it not projective, and the exp goes from 0 to 1 */
  bool is_simple() const {return exp == 1 && oldexp == 0 && !proj;}
  /* Print one root. Projective roots are printed as r+q */
  void fprint(FILE *out, const fbprime_t q) const {
    fprintf(out, "%llu", (unsigned long long) r + (proj ? q : 0));
    if (oldexp != 0 || this->exp != 1)
      fprintf(out, ":%hhu:%hhu", oldexp, this->exp);
  }
  void transform(fb_general_root &result, const fbprime_t q,
                 const redc_invp_t invq,
                 qlattice_basis_srcptr basis) const {
    result.exp = this->exp;
    result.oldexp = oldexp;
    fbprime_t t = r + (proj ? q : 0);
    t = fb_root_in_qlattice(q, t, invq, basis);
    result.proj = (t >= q);
    result.r = t - ( (t >= q) ? q : 0);
  }
};

/* General entries are anything that needs auxiliary information:
   Prime powers, projective roots, ramified primes where exp != oldexp + 1,
   etc. They could, of course, also store the simple cases, but for those we
   use the simple struct to conserve memory and to decide algorithms (batch
   inversion, etc.) statically. */
class fb_general_entry {
  void read_roots (const char *, unsigned char, unsigned char, unsigned long);
public:
  fbprime_t q, p; /* q = p^k */
  redc_invp_t invq; /* invq = -1/q (mod 2^32), or (mod 2^64), depending on
		       the size of redc_invp_t */
  fb_general_root roots[MAXDEGREE];
  unsigned char k, nr_roots;

  fb_general_entry(){}
  template <int Nr_roots>
  fb_general_entry (const fb_entry_x_roots<Nr_roots> &e);
  void parse_line (const char *line, unsigned long linenr);
  void merge (const fb_general_entry &);
  void fprint(FILE *out) const;
  bool is_simple() const;
  void transform_roots(fb_general_entry &, qlattice_basis_srcptr) const;
};

class fb_general_vector;

class fb_interface {
public:
  virtual ~fb_interface(){}
  virtual void append(const fb_general_entry &) = 0;
  virtual void fprint(FILE *) const = 0;
  virtual void count_entries(size_t *nprimes, size_t *nroots, double *weight)
    const = 0;
};

/* The work queue of slices that need to be sieved must return a common base
   class of fb_general_vector and fb_vector<Nr_roots>; this base class is
   fb_vector_interface. It allows determining which kind of vector it is,
   so that the correct siever implementation can be called. */
class fb_vector_interface: public fb_interface {
  public:
  fb_vector_interface(){}
  virtual ~fb_vector_interface(){}
  virtual int get_nr_roots() = 0;
  virtual bool is_general() = 0;
  virtual fb_general_vector * transform_roots(qlattice_basis_srcptr) const = 0;
};

class fb_slices_interface: public fb_interface {
  public:
  virtual ~fb_slices_interface(){}
  virtual fb_vector_interface *get_vector(size_t) const = 0;
  virtual int get_nr_roots() = 0;
  virtual bool is_general() = 0;
};


class fb_general_vector: public std::vector<fb_general_entry>, public fb_vector_interface {
public:
  fb_general_vector(){}
  fb_general_vector(size_type n) : std::vector<fb_general_entry>(n){}
  void append(const fb_general_entry &e) {push_back(e);}
  void fprint(FILE *) const;
  void count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
  int get_nr_roots(){return 0;};
  bool is_general(){return true;};
  fb_general_vector * transform_roots(qlattice_basis_srcptr) const;
};


/* "Simple" factor base entries. We imply q=p, k=1, oldexp=0, exp=1,
   and projective=false for all roots. */
template <int Nr_roots>
class fb_entry_x_roots {
public:
  fbprime_t p;
  fbroot_t roots[Nr_roots];
  fb_entry_x_roots(){};
  /* Allow assignment-construction from general entries */
  fb_entry_x_roots(const fb_general_entry &e) {
    ASSERT_ALWAYS(Nr_roots == e.nr_roots);
    p = e.p;
    /* Should we assert is_simple() here for each root? */
    for (size_t i = 0; i < Nr_roots; i++)
      roots[i] = e.roots[i].r;
  }
  void fprint(FILE *) const;
  void transform_roots(fb_general_entry &, qlattice_basis_srcptr) const;
};


template <int Nr_roots>
class fb_vector: public std::vector<fb_entry_x_roots<Nr_roots> >, public fb_vector_interface {
  public:
  fb_vector(){}
  /* FIXME: using size_type here does not work, why? */
  fb_vector(size_t n) : std::vector<fb_entry_x_roots<Nr_roots> >(n){}
  void append(const fb_general_entry &e){this->push_back(e);};
  void fprint(FILE *) const;
  void count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
  int get_nr_roots(){return Nr_roots;};
  bool is_general(){return false;};
  fb_general_vector * transform_roots(qlattice_basis_srcptr) const;
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
  void append(const fb_general_entry &);
  void fprint(FILE *) const;
  void count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
  int get_nr_roots(){return Nr_roots;};
  bool is_general(){return false;};
  fb_vector<Nr_roots> *get_vector(const size_t n) const {
    return (n < nr_slices) ? &vectors[n] : NULL;
  }
};


/* A "part" is the set of all factor base primes that get sieved over a given
   bucket region size.
   E.g., when we have only 1 level of bucket sorting, then the factor base has
   2 parts: the line-sieved primes, and the bucket-sieved primes. */
class fb_part: public fb_interface, private NonCopyable {
  /* How do we identify the number of slices in each array?
     Should we have
     size_t nr_entries[9];
     here, or let each array end with a "magic value", such as a slice with
     0 entries, or a NULL pointer? */

  /* If true, all entries go in the general vector */
  bool only_general;

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
  fb_slices_interface *get_slices(const unsigned int n) const {
    ASSERT_ALWAYS(n <= MAXDEGREE);
    
    if (only_general)
     return NULL;
    
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
  fb_part(size_t nr_slices, bool only_general=false);
  ~fb_part();
  void append(const fb_general_entry &);
  void fprint(FILE *) const;
  void count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
  bool is_only_general() const {return only_general;}
  fb_general_vector::iterator begin() {
    ASSERT_ALWAYS(only_general);
    return general_vector.begin();
  }
  fb_general_vector::iterator end() {
    ASSERT_ALWAYS(only_general);
    return general_vector.end();
  }
  fb_general_vector::const_iterator cbegin() const {
    ASSERT_ALWAYS(only_general);
    return general_vector.cbegin();
  }
  fb_general_vector::const_iterator cend() const {
    ASSERT_ALWAYS(only_general);
    return general_vector.cend();
  }
  fb_general_vector *get_general_vector(){return &general_vector;}
  fb_vector_interface *get_n_roots_vector(const int n, const size_t slice) {
    return get_slices(n)->get_vector(slice);
  }
};

/* Splits the factor base for a polynomial into disjoint parts which are
   sieved over different sieve region sizes.
   For example,
   parts[0] will contain very small primes that are line-sieved,
   parts[1] contains bucket-sieved primes with 1 level of bucket sorting
   (i.e., hits get sorted into bucket regions of size 2^16)
*/
class fb_factorbase: public fb_interface, private NonCopyable {
  fb_part *parts[FB_MAX_PARTS];
  fbprime_t thresholds[FB_MAX_PARTS];
 public:
  fb_factorbase(const fbprime_t *thresholds, const size_t *nr_slices,
		const bool *only_general=NULL);
  ~fb_factorbase();
  void read(const char * const filename);
  void make_linear (const mpz_t *poly, fbprime_t powbound, bool do_projective);
  bool mmap_fbc(const char *) {return false;};
  void dump_fbc(const char *) {return;};
  unsigned long *extract_bycost(size_t &n, fbprime_t pmax, fbprime_t td_thresh);
  size_t size() const {return 0;}
  void fprint(FILE *) const;
  void append(const fb_general_entry &);
  void count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
  fb_part *get_part(const size_t n) {ASSERT_ALWAYS(n < FB_MAX_PARTS); return parts[n];}
};


unsigned char	fb_log (double, double, double);
fbprime_t       fb_pow (fbprime_t, unsigned long);
fbprime_t       fb_is_power (fbprime_t, unsigned long *);
unsigned char   fb_make_steps(fbprime_t *, fbprime_t, double);

#endif
