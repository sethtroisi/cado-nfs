#ifndef LAS_SIEVE_INFO_HPP_
#define LAS_SIEVE_INFO_HPP_


/* various includes that are needed to define the sieve_info type */
#include "cado_poly.h"
#include "las-siever-config.hpp"
#include "bucket.hpp"  /* for bkmult_specifier */
#include "las-todo-entry.hpp"
#include "las-qlattice.hpp"
#include "las-norms.hpp"
#include "trialdiv.hpp"
#include "fb.hpp"
#include "las-fill-in-buckets.hpp"
#include "las-smallsieve-types.hpp"
#include "las-unsieve.hpp"
#include "ecm/facul.hpp"

#include <memory>
#ifdef HAVE_BOOST_SHARED_PTR
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
namespace std { using boost::shared_ptr; using boost::make_shared; }
#endif

/*
 * A sieve_info struct is only a thin layer above the factor base. It
 * provides the factor base, and the possibility to access some data
 * structures that depend on the current configuration (e.g. on I).

 *
 * It must be regarded as a (mostly) const repository from where we fetch
 * data that is of interest to many special-q's.
 *
 * The data at the sieve_info level should not be attached to a
 * particular config, nor to particular q.
 *
 * In the multi-las context, sieve_info is *shared* between the different
 * sub-processes.
 *
 * TODO: at some point, the query mechanism for the cached entries will
 * depend on the memory binding.
 */
struct sieve_info {
    cxx_cado_poly cpoly;

    /* {{{ side_data */
    struct side_data {
        cxx_mpz_poly f;     /* handy. Alias to cpoly.pols[side] */

        /* fb depends on
         *      cpoly
         *      conf.sides[side].lim
         *      conf.sides[side].powlim
         *
         */
        fb_factorbase fb;

        inline bool no_fb() const { return fb.empty(); }

        fb_factorbase::slicing const * get_factorbase_slicing(fb_factorbase::key_type fbK) {
            /* implemented in fb_factorbase::operator[] (fb.hpp), with
             * locking.
             */
            return &fb[fbK];
        }


        /* trialdiv_data depends on
         *      fbK.td_thresh
         * and marginally on
         *      fbK.thresholds[0]
         */

        private:
        struct equivalent_fbK_for_td {
            bool operator()(fb_factorbase::key_type const & a, fb_factorbase::key_type const & b) const {
                if (a.thresholds[0] < b.thresholds[0]) return true;
                if (a.thresholds[0] > b.thresholds[0]) return false;
                return a.td_thresh < b.td_thresh;
            }
        };
        std::map<fb_factorbase::key_type, trialdiv_data, equivalent_fbK_for_td>
            trialdiv_data_cache;
        public:
        /* in las-trialdiv.cpp */
        trialdiv_data const * get_trialdiv_data(fb_factorbase::key_type fbK, fb_factorbase::slicing const * fbs);

        side_data(int side, cxx_cado_poly const & cpoly, cxx_param_list & pl, bool try_fbc);
    };

    /* }}} */

    side_data sides[2];

    /* Most of the member functions below which take a side argument
     * should be members of the side_data structure instead.
     *
     * Note that have have only init_* functions below, no clear_*. In
     * fact, it's mostly a misnomer, and I apologize for it. The
     * automatic dtor will work fine for all these fields, hence there is
     * no need for specific clear_* functions. The magic lies most of the
     * time in the use of shared_ptr's. Maybe the init_* ones should have
     * been named setup_* or something like this.
     */

    /* in las-unsieve.cpp */
    /* Data for unsieving locations where gcd(i,j) > 1 */
    /* This gets initialized only when I is finally decided */
    private:
    std::map<std::pair<int, int>, unsieve_data> us_cache;
    public:
    unsieve_data const * get_unsieve_data(siever_config const & conf);

    /* in las-unsieve.cpp */
    /* Data for divisibility tests p|i in lines where p|j */
    private:
    std::map<unsigned int, j_divisibility_helper> jdiv_cache;
    public:
    j_divisibility_helper const * get_j_divisibility_helper(int J);

    /* We use a shared_ptr only to provide the custom dtor.
     *
     * We have a cache because during the descent, several cofactoring
     * strategies are possible.
     */
    private:
    const char *cofactfilename;
    std::map<siever_config,
        std::shared_ptr<facul_strategies_t>,
        siever_config::has_same_cofactoring::comparison> facul_strategies_cache;
    public:
    facul_strategies_t const * get_strategies(siever_config const & conf);

    /* These functions must be called before actually sieving */
    void update_norm_data ();

    sieve_info(cxx_cado_poly const & cpoly, cxx_param_list & pl, bool try_fbc = false);
    ~sieve_info();
};

/*  */

#endif	/* LAS_SIEVE_INFO_HPP_ */