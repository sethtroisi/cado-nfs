#include "cado.h"
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cctype>
#include <gmp.h>
#include <pthread.h>
#ifdef HAVE_SYS_MMAN_H
#include <sys/mman.h>
#else
#define MAP_FAILED ((void *) -1)
#endif
#include "fb.hpp"
#include "mod_ul.h"
#include "verbose.h"
#include "getprime.h"
#include "gmp_aux.h"
#include "gzip.h"
#include "threadpool.hpp"
#include "misc.h"

/* {{{ fb_log fb_pow and friends */
/* Returns floor(log_2(n)) for n > 0, and 0 for n == 0 */
static unsigned int
fb_log_2 (fbprime_t n)
{
  unsigned int k;
  for (k = 0; n > 1; n /= 2, k++);
  return k;
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

/* the idea here is that we want to return the rounding of log(p^k) minus
 * the rounding of log(p^(k-1)). So it's close, in spirit, to log(p)...
 *
 * Note that the old fb_log function was taking log_scale as something
 * relative to natural logarithm. So while we had a precomputed scale so
 * that the-log-we-want(p) = log2(p) * scale, here we were feedint
 * scale/log(2), which is pretty awkward.
 *
 * (log2 is C99. I do recall that some bsd libm lacks it. It's a bug,
 * period).
 */
unsigned char
fb_log (double n, double log_scale, double offset)
{
  const long l = floor (log2 (n) * log_scale + offset + 0.5);
  return static_cast<unsigned char>(l);
}

unsigned char
fb_log_delta (const fbprime_t p, const unsigned long newexp,
              const unsigned long oldexp, const double log_scale)
{
  return fb_log (fb_pow(p, newexp), log_scale, 0.)
         - fb_log (fb_pow(p, oldexp), log_scale, 0.);
}
/* }}} */

static bool fb_linear_root (fbroot_t & root, cxx_mpz_poly const & poly, const fbprime_t q);

// Adapted from utils/ularith.h
// TODO: this function should go somewhere else...
static inline uint64_t
uint64_invmod(const uint64_t n)
{
    uint64_t r;
    ASSERT (n % UINT64_C(2) != UINT64_C(0));
    r = (UINT64_C(3) * n) ^ UINT64_C(2);
    r = UINT64_C(2) * r - (uint32_t) r * (uint32_t) r * (uint32_t) n;
    r = UINT64_C(2) * r - (uint32_t) r * (uint32_t) r * (uint32_t) n;
    r = UINT64_C(2) * r - (uint32_t) r * (uint32_t) r * (uint32_t) n;
    uint32_t k = (uint32_t)(r * n >> 32);
    k *= (uint32_t) r;
    r = r - ((uint64_t)k << 32);
    return r;
}

static inline redc_invp_t
compute_invq(fbprime_t q)
{
  if (q % 2 != 0) {
    if (sizeof(unsigned long) >= sizeof(redc_invp_t)) {
        return (redc_invp_t) (- ularith_invmod (q));
    } else {
        ASSERT(sizeof(redc_invp_t) == 8);
        return (redc_invp_t) (- uint64_invmod (q));
    }
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


/* Allow construction of a root from a linear polynomial and a prime (power) */
fb_general_root::fb_general_root (fbprime_t q, cxx_mpz_poly const & poly,
  const unsigned char nexp, const unsigned char oldexp)
  : exp(nexp), oldexp(oldexp)
{
  proj = fb_linear_root (r, poly, q);
}


/* Allow assignment-construction of general entries from simple entries */
template <int Nr_roots>
fb_entry_general::fb_entry_general (const fb_entry_x_roots<Nr_roots> &e) {
  p = q = e.p;
  k = 1;
  invq = e.invq;
  for (int i = 0; i != Nr_roots; i++) {
    /* Use simple constructor for root */
    roots[i] = e.roots[i];
  }
}


/* Return whether this is a simple factor base prime.
   It is simple if it is a prime (not a power) and all its roots are simple. */
bool
fb_entry_general::is_simple() const
{
  bool is_simple = (k == 1);
  for (unsigned char i = 0; i != nr_roots; i++) {
    is_simple &= roots[i].is_simple();
  }
  return is_simple;
}


/* Read roots from a factor base file line and store them in roots.
   line must point at the first character of the first root on the line.
   linenr is used only for printing error messages in case of parsing error.
   Returns the number of roots read. */
void
fb_entry_general::read_roots (const char *lineptr, const unsigned char nexp,
                              const unsigned char oldexp,
                              const unsigned long linenr)
{
    unsigned long long last_t = 0;

    nr_roots = 0;
    while (*lineptr != '\0')
    {
        if (nr_roots == MAX_DEGREE) {
            verbose_output_print (1, 0,
                    "# Error, too many roots for prime (power) %" FBPRIME_FORMAT
                    " in factor base line %lu\n", q, linenr);
            exit(EXIT_FAILURE);
        }
        /* Projective roots r, i.e., ar == b (mod q), are stored as q + r in
           the factor base file; since q can be a 32-bit value, we read the
           root as a 64-bit integer first and subtract q if necessary. */
        const unsigned long long t = strtoull_const (lineptr, &lineptr, 10);
        if (nr_roots > 0 && t <= last_t) {
            verbose_output_print (1, 0,
                "# Error, roots must be sorted in the fb file, line %lu\n",
                linenr);
            exit(EXIT_FAILURE);
        }
        last_t = t;

        roots[nr_roots++] = fb_general_root(t, q, nexp, oldexp);
        if (*lineptr != '\0' && *lineptr != ',') {
            verbose_output_print (1, 0,
                    "# Incorrect format in factor base file line %lu\n",
                    linenr);
            exit(EXIT_FAILURE);
        }
        if (*lineptr == ',')
            lineptr++;
    }

    if (nr_roots == 0) {
        verbose_output_print (1, 0, "# Error, no root for prime (power) %"
                FBPRIME_FORMAT " in factor base line %lu\n", q, linenr - 1);
        exit(EXIT_FAILURE);
    }
}

/* Parse a factor base line.
   Return 1 if the line could be parsed and was a "short version", i.e.,
   without explicit old and new exponent.
   Return 2 if the line could be parsed and was a "long version".
   Otherwise return 0. */
void
fb_entry_general::parse_line (const char * lineptr, const unsigned long linenr)
{
    q = strtoul_const (lineptr, &lineptr, 10);
    if (q == 0) {
        verbose_output_print (1, 0, "# fb_read: prime is not an integer on line %lu\n",
                              linenr);
        exit (EXIT_FAILURE);
    } else if (*lineptr != ':') {
        verbose_output_print (1, 0,
                "# fb_read: prime is not followed by colon on line %lu",
                linenr);
        exit (EXIT_FAILURE);
    }

    lineptr++; /* Skip colon after q */
    const bool longversion = (strchr(lineptr, ':') != NULL);

    /* NB: a short version is not permitted for a prime power, so we
     * do the test for prime powers only for long version */
    p = q;
    unsigned char nexp = 1, oldexp = 0;
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
        nexp = strtoul_const (lineptr, &lineptr, 10);

        if (nexp == 0) {
            verbose_output_print (1, 0, "# Error in fb_read: could not parse "
                "the integer after the colon of prime %" FBPRIME_FORMAT "\n",
                q);
            exit (EXIT_FAILURE);
        }
        if (*lineptr != ',') {
            verbose_output_print (1, 0,
		    "# fb_read: exp is not followed by comma on line %lu",
		    linenr);
            exit (EXIT_FAILURE);
        }
        lineptr++; /* skip comma */
        oldexp = strtoul_const (lineptr, &lineptr, 10);
        if (*lineptr != ':') {
            verbose_output_print (1, 0,
		    "# fb_read: oldlogp is not followed by colon on line %lu",
		    linenr);
            exit (EXIT_FAILURE);
        }
        ASSERT (nexp > oldexp);
        lineptr++; /* skip colon */
    }

    read_roots(lineptr, nexp, oldexp, linenr);

    /* exp and oldexp are a property of a root, not of a prime (power).
       The factor base file should specify them per root, but specifies
       them per prime instead - a bit of a design bug.
       For long version lines, we thus use the exp and oldexp values for all
       roots specified in that line. */
    for (unsigned char i = 1; i < nr_roots; i++) {
      roots[i].exp = roots[0].exp;
      roots[i].oldexp = roots[0].oldexp;
    }
}

void
fb_entry_general::fprint(FILE *out) const
{
  fprintf(out, "%" FBPRIME_FORMAT ": ", q);
  for (unsigned char i_root = 0; i_root != nr_roots; i_root++) {
    roots[i_root].fprint(out, q);
    if (i_root + 1 < nr_roots)
      fprintf(out, ",");
  }
  fprintf(out, "\n");
}

void
fb_entry_general::merge (const fb_entry_general &other)
{
  ASSERT_ALWAYS(p == other.p && q == other.q && k == other.k);
  for (unsigned char i_root = 0; i_root < other.nr_roots; i_root++) {
    ASSERT_ALWAYS(nr_roots < MAX_DEGREE);
    roots[nr_roots++] = other.roots[i_root];
  }
}

void
fb_entry_general::transform_roots(fb_entry_general::transformed_entry_t &result,
                                  const qlattice_basis &basis) const
{
  result.p = p;
  result.q = q;
  result.invq = invq;
  result.k = k;
  result.nr_roots = nr_roots;
  /* TODO: Use batch-inversion here */
  for (unsigned char i_root = 0; i_root != nr_roots; i_root++)
    roots[i_root].transform(result.roots[i_root], q, invq, basis);
}


template <int Nr_roots>
void
fb_entry_x_roots<Nr_roots>::transform_roots(fb_entry_x_roots<Nr_roots>::transformed_entry_t &result, const qlattice_basis &basis) const
{
  result.p = p;
  /* TODO: Use batch-inversion here */
  for (unsigned char i_root = 0; i_root != nr_roots; i_root++) {
    const unsigned long long t = fb_root_in_qlattice(p, roots[i_root], invq, basis);
    result.proj[i_root] = (t >= p);
    result.roots[i_root] = (t < p) ? t : (t - p);
  }
}

// FIXME: why do I have to make those instances explicit???
// If someone knows how to avoid that...

template void
fb_entry_x_roots<0>::transform_roots(fb_transformed_entry_x_roots<0> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<1>::transform_roots(fb_transformed_entry_x_roots<1> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<2>::transform_roots(fb_transformed_entry_x_roots<2> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<3>::transform_roots(fb_transformed_entry_x_roots<3> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<4>::transform_roots(fb_transformed_entry_x_roots<4> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<5>::transform_roots(fb_transformed_entry_x_roots<5> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<6>::transform_roots(fb_transformed_entry_x_roots<6> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<7>::transform_roots(fb_transformed_entry_x_roots<7> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<8>::transform_roots(fb_transformed_entry_x_roots<8> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<9>::transform_roots(fb_transformed_entry_x_roots<9> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<10>::transform_roots(fb_transformed_entry_x_roots<10> &, qlattice_basis const&) const; 



template <int Nr_roots>
void
fb_entry_x_roots<Nr_roots>::fprint(FILE *out) const
{
  fprintf(out, "%" FBPRIME_FORMAT ": ", p);
  for (int i = 0; i != Nr_roots; i++) {
    fprintf(out, "%" FBROOT_FORMAT "%s", roots[i],
	    (i + 1 < Nr_roots) ? "," : "");
  }
  fprintf(out, "\n");
}

/* Make one factor base entry for a linear polynomial poly[1] * x + poly[0]
   and the prime (power) q. We assume that poly[0] and poly[1] are coprime.
   Non-projective roots a/b such that poly[1] * a + poly[0] * b == 0 (mod q)
   with gcd(poly[1], q) = 1 are stored as a/b mod q.
   If do_projective != 0, also stores projective roots with gcd(q, f_1) > 1,
   but stores the reciprocal root.
   Returns true if the roots was projective, and false otherwise. */

static bool
fb_linear_root (fbroot_t & root, cxx_mpz_poly const & poly, const fbprime_t q)
{
  modulusul_t m;
  residueul_t r0, r1;
  bool is_projective;

  modul_initmod_ul (m, q);
  modul_init_noset0 (r0, m);
  modul_init_noset0 (r1, m);

  modul_set_ul_reduced (r0, mpz_fdiv_ui (poly->coeff[0], q), m);
  modul_set_ul_reduced (r1, mpz_fdiv_ui (poly->coeff[1], q), m);

  /* We want poly[1] * a + poly[0] * b == 0 <=>
     a/b == - poly[0] / poly[1] */
  is_projective = (modul_inv (r1, r1, m) == 0); /* r1 = 1 / poly[1] */

  if (is_projective)
    {
      ASSERT_ALWAYS(mpz_gcd_ui(NULL, poly->coeff[1], q) > 1);
      /* Set r1 = poly[0] % q, r0 = poly[1] (mod q) */
      modul_set (r1, r0, m);
      modul_set_ul_reduced (r0, mpz_fdiv_ui (poly->coeff[1], q), m);
      int rc = modul_inv (r1, r1, m);
      ASSERT_ALWAYS(rc != 0);
    }

  modul_mul (r1, r0, r1, m); /* r1 = poly[0] / poly[1] */
  modul_neg (r1, r1, m); /* r1 = - poly[0] / poly[1] */

  root = modul_get_ul (r1, m);

  modul_clear (r0, m);
  modul_clear (r1, m);
  modul_clearmod (m);

  return is_projective;
}

std::ostream& operator<<(std::ostream& o, fb_factorbase::key_type const & k)
{
    o << "scale=" << k.scale
      << ", "
      << "thresholds={";
    for(int i = 0 ; i < FB_MAX_PARTS ; i++) {
        if (i) o << ", ";
        o << k.thresholds[i];
    }
    o << "}";
    return o;
}

/* {{{ counting primes, prime ideals, and so on in the whole factor base.
 * In fact, these functions are not used presently.
 */
struct fb_factorbase::helper_functor_count_primes {
        template<typename T>
        size_t operator()(size_t t0, T const  & x) const {
            return t0 + x.size();
        }
};
struct fb_factorbase::helper_functor_count_prime_ideals {
        template<typename T>
        size_t operator()(size_t t0, T const  & x) const {
            if (T::is_general_type) {
                for(auto const & a : x)
                    t0 += a.get_nr_roots();
                return t0;
            } else {
                return t0 + T::fixed_nr_roots * x.size();
            }
        }
};
struct fb_factorbase::helper_functor_count_weight {
        template<typename T>
        double operator()(double t, T const  & x) const {
            for(auto const & e : x)
                t += e.weight();
            return t;
        }
};
struct fb_factorbase::helper_functor_count_combined {
    size_t & nprimes;
    size_t & nideals;
    double & weight;
    template<typename T>
    void operator()(T const  & x) const {
        for(auto const & e : x) {
            nprimes++;
            nideals += e.get_nr_roots();
            weight += e.weight();
        }
    }
};
struct fb_factorbase::helper_functor_count_primes_interval {
        fbprime_t pmin = 0;
        fbprime_t pmax = std::numeric_limits<fbprime_t>::max();
        template<typename T>
        size_t operator()(size_t t, T const  & x) const {
            for(auto const & e : x) {
                if (e.get_q() >= pmin && e.get_q() < pmax)
                    t++;
            }
            return t;
        }
};
struct fb_factorbase::helper_functor_count_prime_ideals_interval {
        fbprime_t pmin = 0;
        fbprime_t pmax = std::numeric_limits<fbprime_t>::max();
        template<typename T>
        size_t operator()(size_t t, T const  & x) const {
            for(auto const & e : x) {
                if (e.get_q() >= pmin && e.get_q() < pmax)
                    t += e.get_nr_roots();
            }
            return t;
        }
};
struct fb_factorbase::helper_functor_count_weight_interval {
        fbprime_t pmin = 0;
        fbprime_t pmax = std::numeric_limits<fbprime_t>::max();
        template<typename T>
        double operator()(double t, T const  & x) const {
            for(auto const & e : x)
                if (e.get_q() >= pmin && e.get_q() < pmax)
                    t += e.weight();
            return t;
        }
};
struct fb_factorbase::helper_functor_count_combined_interval {
    size_t & nprimes;
    size_t & nideals;
    double & weight;
    helper_functor_count_combined_interval(size_t & np, size_t & ni, double & w) :
        nprimes(np), nideals(ni), weight(w) {}
    fbprime_t pmin = 0;
    fbprime_t pmax = std::numeric_limits<fbprime_t>::max();
    template<typename T>
    void operator()(T const & x) const {
        for(auto const & e : x) {
            nprimes++;
            nideals += e.get_nr_roots();
            weight += e.weight();
        }
    }
};

/* outside visible interface */
size_t fb_factorbase::count_primes() const {
    return multityped_array_fold(helper_functor_count_primes(), 0, entries);
}
size_t fb_factorbase::count_prime_ideals() const {
    return multityped_array_fold(helper_functor_count_primes(), 0, entries);
}
size_t fb_factorbase::count_weight() const {
    return multityped_array_fold(helper_functor_count_weight(), 0, entries);
}
/* }}} */

/* {{{ append. */

/* FIXME: bizarrely, std::dequeue does not work, here. */
struct fb_factorbase::helper_functor_append {
    std::list<fb_entry_general> &pool;
    int deg;
    bool positive = true;
    helper_functor_append(std::list<fb_entry_general> &pool, int deg)
        : pool(pool), deg(deg) {}
    void switch_to_general_entries_only() { positive = false; }
    template<typename T>
        void operator()(T & x) {
            typedef typename T::value_type FB_ENTRY_TYPE;
            constexpr bool isG = FB_ENTRY_TYPE::is_general_type;
            constexpr unsigned char N = FB_ENTRY_TYPE::fixed_nr_roots;
            if (positive != !isG) return;
            /* normally we should not have primes with zero prime
             * ideals above them in the factor base */
            if (!isG && N == 0) return;
            /* likewise, no more roots than the degree of the number
             * field */
            if (!isG && N > deg) return;
            for(auto it = pool.begin(); it != pool.end(); ) {
                fb_entry_general E = std::move(*it);
                /* see above */
                ASSERT(E.nr_roots > 0 && E.nr_roots <= deg);
                bool must_go_to_general = !E.is_simple() || E.k > 1
                    /* || E.q < powlim TODO */;

                bool ok1 = isG && must_go_to_general;
                bool ok2 = !isG && !must_go_to_general && N == E.nr_roots;
                if (ok1 || ok2) {
                    auto it_next = it;
                    ++it_next;
                    x.push_back(std::move(E));
                    pool.erase(it);
                    it = it_next;
                } else {
                    ++it;
                }
            }
        }
};
void fb_factorbase::append(std::list<fb_entry_general> &pool) {
    /* The "positive" hack and the two passes are just here so that we
     * don't needlessly to a complete pass over the full list just to
     * trim a pocketful of special entries.
     */
    helper_functor_append A { pool, mpz_poly_degree(f) };
    multityped_array_foreach(A, entries);
    A.switch_to_general_entries_only();
    multityped_array_foreach(A, entries);
    ASSERT_ALWAYS(pool.empty());
}
/* }}} */

template<int n>
struct fb_entries_interval_factory {
    struct type {
        typedef typename fb_entries_factory<n>::type container_type;
        typename container_type::const_iterator begin, end;
    };
};

/* Goal: count the weight, the number of primes and so on, and find
 * iterators that point to the right portions.
 */
struct helper_functor_dispatch_weight_parts : public fb_factorbase::slicing::stats_type {
    fb_factorbase::key_type K;
    typedef multityped_array<fb_entries_interval_factory, -1, fb_factorbase::MAX_ROOTS+1> intervals_t;
    std::array<intervals_t, FB_MAX_PARTS> intervals;
    helper_functor_dispatch_weight_parts(fb_factorbase::key_type K) : K(K) {}

    /* TODO: factor base entries that made it only to the vector of
     * general entries will perhas deserve a special treatment.
     */
    template<typename T>
        int operator()(int toplevel, T const  & x) {
            /* T is fb_entries_factory<n>::type for some n */
            typedef typename T::const_iterator xit_t;
            typedef typename T::value_type FB_ENTRY_TYPE;
            xit_t it = x.begin();
            /* we're now processing entries with n roots. First make sure
             * that the begin pointer for the first part correctly points
             * to the beginning of our current vector of entries.
             */
            int i = 0;
            constexpr int n = FB_ENTRY_TYPE::is_general_type ? -1 : FB_ENTRY_TYPE::fixed_nr_roots;
            intervals[0].get<n>().begin = it;
            for( ; it != x.end() ; ++it) {
                fbprime_t q = it->get_q();
                /* This really really mandates that we have sorted
                 * vectors of entries in the factor base ! */
                for( ; i < FB_MAX_PARTS && q >= K.thresholds[i] ; ++i) {
                    intervals[i].get<n>().end = it;
                    ASSERT_ALWAYS(i + 1 < FB_MAX_PARTS);
                    intervals[i+1].get<n>().begin = it;
                }
                ASSERT_ALWAYS(i < FB_MAX_PARTS);
                primes[i]++;
                ideals[i]+= it->get_nr_roots(); /* == n, except when n==-1 */
                weight[i] += it->weight();
            }
            intervals[i].get<n>().end = it;
            /* for consistency, make sure all further parts are
             * understood as empty.
             */
            if (i > toplevel) toplevel = i;
            for( ; ++i < FB_MAX_PARTS ; ) {
                intervals[i].get<n>().begin = it;
                intervals[i].get<n>().end = it;
            }
            return toplevel;
        }
};

struct helper_functor_dispatch_small_sieved_primes {
    fb_factorbase::slicing & S;
    unsigned int td_thresh;
    /* TODO: it's a bit unsatisfactory that we do this comparison on
     * it->p for each prime.
     */
    template<typename T>
        void operator()(T const  & x) {
            /* T is fb_entries_interval_factory<n>::type for some n */
            /* entries between x.begin and x.end go to the vectors of
             * small-sieved primes */
            /* the product td_thresh * it->get_nr_roots() will be
             * folded outside the loop for the most common cases.
             */
            for(auto it = x.begin ; it != x.end ; ++it) {
                fb_entry_general G(*it);
                if (it->k > 1 || it->p <= td_thresh * it->get_nr_roots()) {
                    S.small_sieve_entries.rest.push_back(G);
                } else {
                    S.small_sieve_entries.resieved.push_back(G);
                }
            }
        }
};

struct helper_functor_subdivide_slices {
    fb_factorbase::slicing::part & dst;
    fb_factorbase::key_type K;
    slice_index_t index = 0;
    helper_functor_subdivide_slices(fb_factorbase::slicing::part & dst, fb_factorbase::key_type const & K) : dst(dst), K(K), index(0) {}
    template<typename T>
        void operator()(T const & x) {
            /* T is fb_entries_interval_factory<n>::type for some n */
            typedef T interval_t;
            typedef typename interval_t::container_type ventry_t;
            typedef typename ventry_t::value_type entry_t;
            typedef fb_slice<entry_t> slice_t;
            typedef std::vector<slice_t> vslice_t;

            constexpr int n = entry_t::is_general_type ? -1 : entry_t::fixed_nr_roots;
            // vslice_t & sdst(dst.get_slices_vector_for_nroots<n>());
            vslice_t & sdst(dst.slices.get<n>());
            size_t interval_width = x.end - x.begin;
            if (!interval_width) return;
            /* first scan to separate by values of logp */
            typename std::vector<entry_t>::const_iterator it = x.begin;
            typename std::vector<slice_t> pool;
            /* can't emplace_back here because private ctor */
            pool.push_back(slice_t(it, fb_log(it->get_q(), K.scale, 0)));
            for( ; it != x.end ; ++it) {
                unsigned char cur_logp = fb_log(it->get_q(), K.scale, 0);
                if (cur_logp == pool.back().logp) {
                    pool.back().weight += it->weight();
                    continue;
                }
                pool.back()._end = it;
                pool.push_back(slice_t(it, cur_logp));
                pool.back().weight += it->weight();
            }
            pool.back()._end = it;
            std::ostringstream n_eq;
            n_eq << "n=";
            if (n < 0) n_eq << "*"; else n_eq << n;
            verbose_output_print (0, 2, "# slices for %s roots: from %zu entries, we found %zu different logp values\n", n_eq.str().c_str(), interval_width, pool.size());
            for(auto const & s : pool) {
                verbose_output_print (0, 2, "#  %s logp=%d: %zu entries, weight=%f\n",
                        n_eq.str().c_str(),
                        (int) s.get_logp(),
                        s.end() - s.begin(),
                        s.get_weight());
            }

            /* bold move: all of them become slices... */
            /* of course that won't stay, it's obvious WIP.
             * TODO */
            sdst.swap(pool);

            /* And then we number all slices */
            for(auto & s : sdst) s.index = index++;
        }
};



fb_factorbase::slicing::slicing(fb_factorbase const & fb, fb_factorbase::key_type const & K) {

    /* parts[0] will be unused.
     * parts[1] will have primes between K.thresholds[1] and K.thresholds[2]
     * parts[2] will have primes above K.thresholds[2].
     */

    /* First thing we're going to do is count the weight of each part.
     * That used to be done in fb_part::_count_entries
     *
     * Now this is all done in the "foreach" below, which does a lot of
     * stuff using the helper_functor_dispatch_weight_parts structure above.
     */

    std::ostringstream os;
    os << K;
    verbose_output_print(0, 1, "# Creating new slicing on side %d for %s\n",
            fb.side, os.str().c_str());

    helper_functor_dispatch_weight_parts D { K };
    toplevel = multityped_array_fold(D, 0, fb.entries);
    for (int i = 0; i <= toplevel; i++) {
        size_t nr_primes = D.primes[i];
        size_t nr_roots = D.ideals[i];
        double weight = D.weight[i];
        int side = fb.side;
        if (!nr_primes) continue;
            verbose_output_print(0, 1, "# Number of primes in side-%d factor base part %d = %zu\n",
                    side, i, nr_primes);
            verbose_output_print(0, 1, "# Number of prime ideals in side-%d factor base part %d = %zu\n",
                    side, i, nr_roots);
            verbose_output_print(0, 1, "# Weight of primes in side-%d factor base part %d = %0.5g\n",
                    side, i, weight);
    }

    /* D.weight[i] is now what used to be called max_bucket_fill_ratio. We
     * will now make sure that slices are small enough so that a single
     * slice never ever exceeds some fraction of this weight  */
    stats = D;
     
    /* First, part 0 is treated in a special way. There's no slicing to
     * speak of. We simply populate the small_sieve_entries struct,
     * according to the rules that pertain to this data (which primes are
     * resieved, which are trial-divided, and so on).
     */

    multityped_array_foreach(helper_functor_dispatch_small_sieved_primes { *this, K.td_thresh }, D.intervals[0]);
    auto by_q = fb_entry_general::sort_byq();

    /* small sieve cares about this list being sorted by hit rate, I
     * believe. This is tricky because small sieve considers the hit rate
     * *per line*, and given some q here, the special-q lattice may
     * change the picture somewhat if the prime becomes projective: we
     * may have (p^k1,r1) and (p^k2,r2) two distinct roots above p with
     * k1 < k2, yet if only the latter becomes projective the hit rate
     * per line becomes p^k1 and p^(k2-1) (for example). So there's clear
     * potential for the ordering to be swapped.
     */
    std::sort(small_sieve_entries.resieved.begin(), small_sieve_entries.resieved.end(), by_q);
    std::sort(small_sieve_entries.rest.begin(), small_sieve_entries.rest.end(), by_q);

    /* Next, we have sets of begin and end pointers. We need to subdivide
     * them.
     */
    for (int i = 1; i < FB_MAX_PARTS; i++)
        multityped_array_foreach(helper_functor_subdivide_slices { parts[i], K }, D.intervals[i]);



    /* we're going to divide our vector in several
     * parts, and compute slices. That used to be done by many
     * functions:
     *
     * fb_factorbase::append
     * fb_part::append
     * fb_vector::append
     *
     * print_fb_statistics (yes, it's true)
     * fb_read
     * fb_factorbase::make_slices
     * sieve_info::update
     *
     * part0: will be small-sieved, so we don't care about
     * slices.
     * part1: will be bucket-sieved in N passes if toplevel==N
     * part2: if toplevel>=2: in N-1 passes if toplevel==N
     * and so on.
     */
    /* Append a factor base entry to the factor base.
   The new entry is inserted into the correct part, as determined by the
   size of the prime p, and within that part, into the correct fb_slice, as
   determined by the number of roots. */

#if 0
void
fb_factorbase::append(const fb_entry_general &fb_cur)
{
  int i;
  static bool printed_too_large_prime_warning = false;

  /* Find the smallest threshold t such that t >= q */
  for (i = 0; i < FB_MAX_PARTS && fb_cur.q > thresholds[i]; i++);
  /* No prime > largest threshold should ever be added */
  if (i == FB_MAX_PARTS) {
    if (!printed_too_large_prime_warning) {
      verbose_output_print(1, 0, "# Factor base entry %" FBPRIME_FORMAT
                           " is above factor base bound, skipping it "
                           "(and all other too large entries)\n", fb_cur.q);
      printed_too_large_prime_warning = true;
    }
    return; /* silently skip this entry */
  }
  if (i > toplevel) {
      toplevel = i;
  }
  parts[i]->append(fb_cur);
}



/* Append a factor base entry given in fb_cur to the
   correct vector, as determined by the number of roots.
   "Special" primes are added to the general-primes vector.
   This function is a de-multiplexer: Using "switch" to turn a run-time value
   (fb_cur->nr_roots) into a compile-time constant which can be used as a
   template specifier. */
void
fb_part::append(const fb_entry_general &fb_cur)
{
  /* Non-simple ones go in the general vector, or, if this is an only_general
     part, then all entries do */
  if (only_general || !fb_cur.is_simple()) {
    /* Append powers only up to this part's powlim */
    if (fb_cur.k == 1 || fb_cur.q < powlim)
      general_vector.append(fb_cur);
    return;
  }
  get_slices(fb_cur.nr_roots)->append(fb_cur);
}

    fb_vector::
          void append(const fb_entry_general &new_entry) {
    ASSERT_ALWAYS(!read_only);
    if (size >= alloc) {
      alloc += std::max(alloc, static_cast<size_t>(16));
      data   = static_cast<FB_ENTRY_TYPE *>(realloc(data, alloc * sizeof(FB_ENTRY_TYPE)));
      if (data == NULL)
        throw std::bad_alloc();
    }
    memset(data+size, 0, sizeof(FB_ENTRY_TYPE));
    data[size++] = new_entry;
  }

    void sieve_info::print_fb_statistics(int side)
{
    side_info & s(sides[side]);
    if (!s.fb) return;
    for (int i_part = 0; i_part < FB_MAX_PARTS; i_part++)
    {
        size_t nr_primes, nr_roots;
        double weight;
        s.fb->get_part(i_part)->count_entries(&nr_primes, &nr_roots, &weight);
        if (nr_primes != 0 || weight != 0.) {
            verbose_output_print(0, 1, "# Number of primes in side-%d factor base part %d = %zu\n",
                    side, i_part, nr_primes);
            verbose_output_print(0, 1, "# Number of prime ideals in side-%d factor base part %d = %zu\n",
                    side, i_part, nr_roots);
            verbose_output_print(0, 1, "# Weight of primes in side-%d factor base part %d = %0.5g\n",
                    side, i_part, weight);
        }
        s.max_bucket_fill_ratio[i_part] = weight;
    }
}

#endif
}

/* {{{ Generation of the factor base on the rational side */

struct fb_power_t {
    fbprime_t p, q;
    unsigned char k;
    inline bool operator<(fb_power_t const & o) const { return q < o.q; }
};

/* Create a list of prime powers (with exponent >1) up to powlim */
static std::vector<fb_power_t> fb_powers (fbprime_t powlim)
{
    std::vector<fb_power_t> powers;

    prime_info pi;
    prime_info_init(pi);
    for (fbprime_t p = 2; p <= powlim / p; p = getprime_mt(pi)) {
        unsigned char k = 2;
        for(fbprime_t q = p ; (q <= powlim / p) ; k++) {
            q *= p;
            powers.push_back(fb_power_t {p, q, k});
        }
    }
    prime_info_clear(pi);

    std::sort (powers.begin(), powers.end());
    return powers;
}

/* Generate a factor base with primes <= bound and prime powers <= powbound
 * for a linear polynomial. If projective != 0, adds projective roots
 * (for primes that divide leading coefficient).
 */

/*{{{ sequential code */
void fb_factorbase::make_linear (cxx_mpz_poly const & poly, unsigned long lim, unsigned long powlim)
{
    /* Prepare for computing powers up to that limit */
    std::vector<fb_power_t> powers(fb_powers(powlim));
    size_t next_pow = 0;

    std::string polystring = poly.print_poly("x");

    verbose_output_print(0, 1,
            "# Making factor base for polynomial g(x) = %s,\n"
            "# including primes up to %lu"
            " and prime powers up to %lu.\n",
            polystring.c_str(), lim, powlim);

    prime_info(pi);

    std::list<fb_entry_general> pool;
    size_t pool_size = 0;

    for (fbprime_t next_prime = 2; next_prime <= lim; ) {
        fb_entry_general fb_cur;
        /* Handle any prime powers that are smaller than next_prime */
        if (next_pow < powers.size() && powers[next_pow].q <= next_prime) {
            /* The list of powers must not include primes */
            ASSERT_ALWAYS(powers[next_pow].q < next_prime);
            fb_cur.q = powers[next_pow].q;
            fb_cur.p = powers[next_pow].p;
            fb_cur.k = powers[next_pow].k;
            next_pow++;
        } else {
            fb_cur.q = fb_cur.p = next_prime;
            fb_cur.k = 1;
            next_prime = getprime_mt(pi);
        }
        fb_cur.nr_roots = 1;
        fb_cur.roots[0].exp = fb_cur.k;
        fb_cur.roots[0].oldexp = fb_cur.k - 1U;
        fb_cur.roots[0].proj = fb_linear_root (fb_cur.roots[0].r, poly, fb_cur.q);
        fb_cur.invq = compute_invq(fb_cur.q);
        pool.push_back(fb_cur);
        if (++pool_size >= 1024) {
            /* enough to do a batch fill */
            append(pool);
            ASSERT_ALWAYS(pool.empty());
            pool_size = 0;
        }
    }
    append(pool);

    prime_info_clear(pi);

    finalize();
}
/*}}}*/
/* {{{ Parallel version, using thread pool. TODO: simplify ! */

#define GROUP 1024

// A task will handle 1024 entries before returning the result to the
// master thread.
// This task_info structure contains:
//   - general info (poly, number of valid entries)
//   - input for the computation
//   - output of the computation
typedef struct {
    mpz_poly_srcptr poly;
    unsigned int n;

    fbprime_t p[GROUP];
    fbprime_t q[GROUP];
    unsigned char k[GROUP];

    fbroot_t r[GROUP];
    bool proj[GROUP];
    redc_invp_t invq[GROUP];
} task_info_t;


class make_linear_thread_param: public task_parameters {
    public:
        task_info_t *T;
        make_linear_thread_param(task_info_t *_T) : T(_T) {}
        make_linear_thread_param() {}
};

class make_linear_thread_result: public task_result {
    public:
        task_info_t *T;
        const make_linear_thread_param *orig_param;
        make_linear_thread_result(task_info_t *_T, const make_linear_thread_param *_p)
            : T(_T), orig_param(_p) {
                ASSERT_ALWAYS(T == orig_param->T);
            }
};

    static task_result *
process_one_task(const worker_thread * worker MAYBE_UNUSED, const task_parameters *_param)
{
    const make_linear_thread_param *param =
        static_cast<const make_linear_thread_param *>(_param);
    task_info_t *T = param->T;
    for (unsigned int i = 0; i < T->n; ++i) {
        T->proj[i] = fb_linear_root (T->r[i], T->poly, T->q[i]);
        T->invq[i] = compute_invq(T->q[i]);
    }
    return new make_linear_thread_result(T, param);
}


// Prepare a new task. Return 0 if there are no new task to schedule.
// Otherwise, return the number of ideals put in the task.
static int get_new_task(task_info_t &T, fbprime_t &next_prime, prime_info& pi, const fbprime_t maxp, size_t &next_pow, std::vector<fb_power_t> const & powers)
{
    unsigned int i;
    for (i = 0; i < GROUP && next_prime <= maxp; ++i) {
        if (next_pow < powers.size() && powers[next_pow].q <= next_prime) {
            ASSERT_ALWAYS(powers[next_pow].q < next_prime);
            T.q[i] = powers[next_pow].q;
            T.p[i] = powers[next_pow].p;
            T.k[i] = powers[next_pow].k;
            next_pow++;
        } else {
            T.q[i] = T.p[i] = next_prime;
            T.k[i] = 1;
            next_prime = getprime_mt(pi);
        }
    }
    T.n = i;
    return i;
}

static void store_task_result(fb_factorbase &fb, task_info_t const & T)
{
    std::list<fb_entry_general> pool;
    for (unsigned int j = 0; j < T.n; ++j) {
        fb_entry_general fb_cur;
        fb_cur.q = T.q[j];
        fb_cur.p = T.p[j];
        fb_cur.k = T.k[j];
        fb_cur.nr_roots = 1;
        fb_cur.roots[0].exp = fb_cur.k;
        fb_cur.roots[0].oldexp = fb_cur.k - 1U;
        fb_cur.roots[0].proj = T.proj[j];
        fb_cur.roots[0].r = T.r[j];
        fb_cur.invq = T.invq[j];
        pool.push_back(fb_cur);
    }
    fb.append(pool);
}

void fb_factorbase::make_linear_threadpool (cxx_mpz_poly const & poly, unsigned long lim, unsigned long powlim, unsigned int nb_threads)
{
    /* Prepare for computing powers up to that limit */
    std::vector<fb_power_t> powers(fb_powers(powlim));
    size_t next_pow = 0;

    std::string polystring = poly.print_poly("x");

    verbose_output_print(0, 1,
            "# Making factor base for polynomial g(x) = %s,\n"
            "# including primes up to %lu"
            " and prime powers up to %lu"
            " using threadpool of %u threads.\n",
            polystring.c_str(), lim, powlim, nb_threads);

#define MARGIN 3
    // Prepare more tasks, so that threads keep being busy.
    unsigned int nb_tab = nb_threads + MARGIN;
    task_info_t * T = new task_info_t[nb_tab];
    make_linear_thread_param * params = new make_linear_thread_param[nb_tab];
    for (unsigned int i = 0; i < nb_tab; ++i) {
        T[i].poly = poly;
        params[i].T = &T[i];
    }

    fbprime_t maxp = lim;
    fbprime_t next_prime = 2;

    prime_info pi;
    prime_info_init(pi);

    thread_pool pool(nb_threads);

    // Stage 0: prepare tasks
    unsigned int active_task = 0;
    for (unsigned int i = 0; i < nb_tab; ++i) {
        int ret;
        ret = get_new_task(T[i], next_prime, pi, maxp, next_pow, powers);
        if (!ret)
            break;
        pool.add_task(process_one_task, &params[i], 0);
        active_task++;
    }

    // Stage 1: while there are still primes, wait for a result and
    // schedule a new task.
    for(int cont = 1 ; cont && active_task ; ) {
        task_result *result = pool.get_result();
        make_linear_thread_result *res =
            static_cast<make_linear_thread_result *>(result);
        active_task--;
        task_info_t * curr_T = res->T;
        store_task_result(*this, *curr_T);
        cont = get_new_task(*curr_T, next_prime, pi, maxp, next_pow, powers);
        if (cont) {
            active_task++;
            pool.add_task(process_one_task, res->orig_param, 0);
        }
        delete result;
    }

    // Stage 2: purge last tasks
    for (unsigned int i = 0; i < active_task; ++i) {
        task_result *result = pool.get_result();
        make_linear_thread_result *res =
            static_cast<make_linear_thread_result *>(result);
        task_info_t * curr_T = res->T;
        store_task_result(*this, *curr_T);
        delete result;
    }

    delete [] T;
    delete [] params;
    prime_info_clear(pi);
    finalize();
}
/* }}} */

/* }}} */

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

    int
fb_factorbase::read(const char * const filename, unsigned long lim, unsigned long powlim)
{
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
        verbose_output_print (1, 0, "# Could not open file %s for reading\n",
                filename);
        return 0;
    }

    std::list<fb_entry_general> pool;
    int pool_size = 0;
    size_t overflow = 0;
    while (!feof(fbfile)) {
        /* Sadly, the size parameter of fgets() is of type int */
        if (fgets (line, static_cast<int>(linesize), fbfile) == NULL)
            break;
        linenr++;
        if (read_strip_comment(line) == (size_t) 0) {
            /* Skip empty/comment lines */
            continue;
        }

        fb_entry_general C;
        C.parse_line (line, linenr);
        if (C.q > lim || (C.k > 1 && C.q > powlim)) {
            overflow++;
            continue;
        }
        C.invq = compute_invq(C.q);

        if (C.p > maxprime) maxprime = C.p;

        if (pool.empty() || C.q != pool.back().q) {
            pool.push_back(std::move(C));
            pool_size++;
        } else {
            pool.back().merge(C);
        }

        /* fb_fprint_entry (stdout, fb_cur); */
        nr_primes++;

        if (pool_size >= 1024) {
            /* enough to do a batch fill */
            append(pool);
            ASSERT_ALWAYS(pool.empty());
            pool_size = 0;
        }
    }

    append(pool);

    verbose_output_print (0, 2, "# Factor base successfully read, %lu primes, "
            "largest was %" FBPRIME_FORMAT "\n",
            nr_primes, maxprime);
    if (overflow) {
        verbose_output_print (0, 2, "# Note: %zu primes above limits (lim=%lu, powlim=%lu) were discarded\n", overflow, lim, powlim);
    }

    fclose_maybe_compressed (fbfile, filename);

    finalize();
    return 1;
}


/*  Factor base handling */
fb_factorbase::fb_factorbase(cxx_cado_poly const & cpoly, int side, unsigned long lim, unsigned long powlim, cxx_param_list & pl) : f(cpoly->pols[side]), side(side)
{
    char paramname[5];

    if (!lim) return;

    cxx_mpz_poly pol;
    mpz_poly_set(pol, cpoly->pols[side]);

    if (pol->deg > 1) {
        double tfb = seconds ();
        double tfb_wct = wct_seconds ();
        std::string polystring = pol.print_poly("x");
        verbose_output_print(0, 1,
                "# Reading side-%d factor base from disk"
                " for polynomial f%d(x) = %s\n",
                side, side, polystring.c_str());
        snprintf(paramname, sizeof(paramname), "fb%d", side);
        const char * fbfilename = param_list_lookup_string(pl, paramname);
        if (!fbfilename) {
            fprintf(stderr, "Error: factor base file for side %d is not given\n", side);
            exit(EXIT_FAILURE);
        }
        verbose_output_print(0, 1, "# Reading side-%d factor base from %s\n", side, fbfilename);
        if (!read(fbfilename, lim, powlim))
            exit(EXIT_FAILURE);
        tfb = seconds () - tfb;
        tfb_wct = wct_seconds () - tfb_wct;
        verbose_output_print(0, 1,
                "# Reading side-%d factor base took %1.1fs (%1.1fs real)\n",
                side, tfb, tfb_wct);
    } else {
        double tfb = seconds ();
        double tfb_wct = wct_seconds ();
        /* note: we parse again the -t option here -- it gets parsed
         * in the las_info ctor too */
        int nb_threads = 1;		/* default value */
        param_list_parse_int(pl, "t", &nb_threads);
        make_linear_threadpool (pol, lim, powlim, nb_threads);
        tfb = seconds () - tfb;
        tfb_wct = wct_seconds() - tfb_wct;
        verbose_output_print(0, 1,
                "# Creating side-%d rational factor base took %1.1fs (%1.1fs real)\n",
                side, tfb, tfb_wct);
    }
}

union fbc_header {
    char h[256];
    /* TODO: add a checksum of the polynomial, or maybe even the
     * polynomial itself ? */
    struct {
        unsigned long lim;
        unsigned long powlim;
    } sc;
};
fb_factorbase::fb_factorbase(cxx_cado_poly const & cpoly, int side, unsigned long lim, unsigned long powlim, FILE * fbc_filename) : f(cpoly->pols[side]), side(side)
{
    fbc_header header;
    int nr = fread(&header, 1, 256, fbc_filename);
    ASSERT_ALWAYS(nr == sizeof(header.h));
    if (lim != header.sc.lim) {
        fprintf(stderr, "Fatal error: cached factor base file is not consistent with the configured value lim%d=%lu\n", side, lim);
        exit(EXIT_FAILURE);
    }
    if (powlim != header.sc.powlim) {
        fprintf(stderr, "Fatal error: cached factor base file is not consistent with the configured value powlim%d=%lu\n", side, powlim);
        exit(EXIT_FAILURE);
    }

    /* Oh and this is all still TODO, of course ! */

    /* We assume that the cache file gets read in order, side 0 before
     * side 1.
     *
     * The difficulty is that we want the (allocator of the) internal
     * container for entries in the factor base to be mmap-able.  But
     * then, that means an immutable container, clearly !
     */

    /* TODO: re-enable */
#if 0
    const char * fbcfilename = param_list_lookup_string(pl, "fbc");

    if (fbcfilename != NULL) {
        /* Try to read the factor base cache file. If that fails, because
           the file does not exist or is not compatible with our parameters,
           it will be written after we generate the factor bases. */
        verbose_output_print(0, 1, "# Mapping memory image of factor base from file %s\n",
                fbcfilename);
        if (fb_mmap_fbc(fb, fbcfilename)) {
            verbose_output_print(0, 1, "# Finished mapping memory image of factor base\n");
            return;
        } else {
            verbose_output_print(0, 1, "# Could not map memory image of factor base\n");
        }
    }
#endif
#if 0

    if (fbcfilename != NULL) {
        verbose_output_print(0, 1, "# Writing memory image of factor base to file %s\n", fbcfilename);
        /* TODO: re-enable */
        fb_dump_fbc(fb, fbcfilename);
        verbose_output_print(0, 1, "# Finished writing memory image of factor base\n");
    }

    /* Note that max_bucket_fill_ratio and friends are set from within
     * print_fb_statistics, which is a bit ugly.
     *
     * (used to)
     */
#endif
}

template <class FB_ENTRY_TYPE>
plattices_vector_t
fb_slice<FB_ENTRY_TYPE>::make_lattice_bases(const qlattice_basis &basis,
    const int logI, const sublat_t &sublat) const
{
  typename FB_ENTRY_TYPE::transformed_entry_t transformed;
  /* Create a transformed vector and store the index of the fb_slice we currently
     transform */

  plattices_vector_t result(index);
  slice_offset_t i_entry = 0;
  for (auto it = _begin; it != _end; ++it, i_entry++) {
    if (!basis.is_coprime_to(it->p))
      continue;
    it->transform_roots(transformed, basis);
    for (unsigned char i_root = 0; i_root != transformed.nr_roots; i_root++) {
      const fbroot_t r = transformed.get_r(i_root);
      const bool proj = transformed.get_proj(i_root);
      /* If proj and r > 0, then r == 1/p (mod p^2), so all hits would be in
         locations with p | gcd(i,j). */
      if (LIKELY(!proj || r == 0)) {
        plattice_info_t pli = plattice_info_t(transformed.get_q(), r, proj, logI);
        plattice_enumerate_t ple = plattice_enumerate_t(pli, i_entry, logI, sublat);
        // Skip (0,0) unless we have sublattices.
        if (!sublat.m)
          ple.next();
        if (LIKELY(pli.a0 != 0)) {
          result.push_back(ple);
        }
      }
    }
  }
  /* This is moved, not copied */
  return result;
}

