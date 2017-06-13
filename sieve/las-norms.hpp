#ifndef LAS_NORMS_HPP_
#define LAS_NORMS_HPP_

#include <stdint.h>
#include "las-forwardtypes.hpp"
#include "las-siever-config.hpp"
#include "las-types.hpp"
#include "double_poly.h"


/* Initialize lognorms for the bucket_region number J. It's a wrapper.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime.
 * The sieve area S must be preallocated with at least (BUCKET_REGION +
 * MEMSET_MIN) space. Indeed, the algorithm that is used might write a
 * bit beyond the meaningful area.
 */
void init_norms_bucket_region (unsigned char *S, uint32_t J, sieve_info& si, unsigned int side, unsigned int smart);

double get_maxnorm_alg (double_poly_srcptr src_poly, const double X, const double Y);

void sieve_info_init_norm_data_sq (sieve_info& si, unsigned long q);

  
struct smart_norm_root {
  unsigned char derivative; /* 0 = root of F; 1 = root of F'; and so on */
  double value;             /* root of F(i,1) or derivatives of F. */
  smart_norm_root(unsigned char derivative = 0, double value = 0) : derivative(derivative), value(value) {}
};

/* These segments ((x, F(x)), (y, F(y))) are used in the smart normalization */
typedef struct sg_s {
  int begin, end;
  double f_begin, f_end;
} sg_t;

/* These functions are internals. Don't use them. Use the wrapper above.
   It's need to declare them here for units & coverage tests.
 */
void init_degree_one_norms_bucket_region_internal     (unsigned char *S, uint32_t J, uint32_t I, double scale, cxx_double_poly const &, double *cexp2);
void init_exact_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, cxx_double_poly const & fijd);
void init_smart_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, cxx_double_poly const & fijd, std::vector<smart_norm_root> const & roots);
void init_norms_roots_internal (cxx_double_poly const &, double max_abs_root, double precision, std::vector<smart_norm_root> & roots);
void init_norms_roots (sieve_info & si, unsigned int side);

struct sieve_range_adjust {/*{{{*/
    friend struct sieve_info;
private:
    las_todo_entry doing;
    siever_config conf;         /* This "conf" field is only used for a
                                 * few fields. In particular the
                                 * large prime bounds. We're specifically
                                 * *not* using the sieving fields, since
                                 * by design these can be decided *after*
                                 * the adjustment.  */
    cado_poly_srcptr cpoly;
    int nb_threads;
    cxx_double_poly fijd[2];
    int logA;
public:
    int logI;
    int J;
    qlattice_basis Q;

#if 0
    sieve_range_adjust(las_todo_entry const & doing, las_info const & las)
        : doing(doing), cpoly(las.cpoly), nb_threads(las.nb_threads)
    {
        /* See whether for this size of special-q, we have predefined
         * parameters (note: we're copying the default config, and then
         * we replace by an adjusted one if needed). */
        conf = las.get_config_for_q(doing);
        /* These two will be adjusted in the process */
        logA = conf.logA;
        logI = J = 0;
    }
#endif

    /* This is only for desperate cases. In las-duplicates, for the
     * moment it seems that we're lacking the las_info structure... */
    sieve_range_adjust(las_todo_entry const & doing, cado_poly_srcptr cpoly, siever_config const & conf, int nb_threads = 1)
        : doing(doing), conf(conf), cpoly(cpoly), nb_threads(nb_threads)
    {
        logA = conf.logA;
        logI = J = 0;
    }


    int SkewGauss() {
        int ret = ::SkewGauss(Q, doing.p, doing.r, cpoly->skew);
        Q.set_q(doing.p, doing.prime_sq);
        if (!Q.prime_sq)
            Q.prime_factors = Q.prime_factors;
        return ret;
    }

    /* There are three strategies to do a post-SkewGauss adjustment of
     * the q-lattice basis.  */

    /* implementation is in las-norms.cpp */
    // all these functions return 0 if they feel that the special-q
    // should be discarded.
    int sieve_info_adjust_IJ();    // "raw" J.
    int sieve_info_update_norm_data_Jmax(bool keep_logI = false);
    int adjust_with_estimated_yield();

    // a fall-back measure for desperate cases.
    // XXX when estimated_yield() wins, this will probably no longer be
    // necessary.
    void set_minimum_J_anyway();
    siever_config const& config() const { return conf; }
private:
    template<typename T> struct mat {
        T x[4];
        T const& operator()(int i, int j) const { return x[2*i+j]; }
        T & operator()(int i, int j) { return x[2*i+j]; }
        mat(T a, T b, T c, T d) { x[0]=a; x[1]=b; x[2]=c; x[3]=d; }
        mat(T y[4]) { x[0]=y[0]; x[1]=y[1]; x[2]=y[2]; x[3]=y[3]; }
    };
    template<typename T> struct vec {
        T x[2];
        vec(T a, T b) { x[0] = a; x[1] = b; }
        vec(T y[2]) { x[0] = y[0]; x[1] = y[1]; }
        T const& operator[](int i) const { return x[i]; }
        T & operator[](int i) { return x[i]; }
        T const& operator()(int i) const { return x[i]; }
        T & operator()(int i) { return x[i]; }
    };
    friend sieve_range_adjust::vec<double> operator*(sieve_range_adjust::vec<double> const& a, sieve_range_adjust::mat<int> const& m) ;
    friend qlattice_basis operator*(sieve_range_adjust::mat<int> const& m, qlattice_basis const& Q) ;
    void prepare_fijd();
    int adapt_threads(const char *);
    double estimate_yield_in_sieve_area(mat<int> const& shuffle, int squeeze, int N);
};/*}}}*/
#endif	/* LAS_NORMS_HPP_ */
