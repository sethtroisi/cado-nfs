#ifndef LAS_SIEVE_INFO_HPP_
#define LAS_SIEVE_INFO_HPP_


/* various includes that are needed to define the sieve_info type */
#include "cado_poly.h"
#include "las-siever-config.hpp"
#include "bucket.hpp"  /* for bkmult_specifier */
#include "las-todo-entry.hpp"
#include "las-qlattice.hpp"
#include "las-norms.hpp"
#include "trialdiv.h"
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

/* General information about the siever, based on some input-dependent
 * configuration data with an impact on the output (as opposed to e.g.
 * file names, or verbosity flags, which do not affect the output).
 */
struct sieve_info {
    cado_poly_srcptr cpoly; /* The polynomial pair */

    /* This conditions the validity of the sieve_info_side members
     * sides[0,1], as well as some other members */
    siever_config conf;

    /* This field gets initialized only via get_sieve_info_from_config */
    const bkmult_specifier * bk_multiplier = nullptr;

    las_todo_entry doing;

    // sieving area. Note that in the (conf) member struct, we find
    // logI_adjusted (will be renamed eventually), as well as logA.
    uint32_t J;
    uint32_t I;

    // description of the q-lattice. The values here should remain
    // compatible with those in ->conf (this concerns notably the bit
    // size as well as the special-q side).
    qlattice_basis qbasis;

    // parameters for bucket sieving
    /* Actual number of buckets at toplevel used by current special-q */
    uint32_t nb_buckets[FB_MAX_PARTS];

    /* Largest level for which the corresponding fb_part is not empty */
    int toplevel;

    /* {{{ sieve_side_info */
    struct side_info {
        std::shared_ptr<lognorm_base> lognorms;

        std::shared_ptr<trialdiv_divisor_t> trialdiv_data;
        std::shared_ptr<std::vector<fb_entry_general> > fb_smallsieved;
        /* Iterators that point into fb_smallsieved; the entries between these
           two iterators are to be small-sieved, the others are not. */
        size_t resieve_start_offset, resieve_end_offset;

        /* The reading, mapping or generating the factor base all create the
         * factor base in several pieces: small primes, and large primes split
         * into one piece for each thread.
         *
         * [§23.2.4.9 : references to an associative container remain valid]
         */
        std::shared_ptr<fb_factorbase> fb;
        fb_factorbase::slicing * fbs = NULL;

        /* Caching of the FK-basis in sublat mode */
        precomp_plattice_dense_t precomp_plattice_dense;

        /* When threads pick up this sieve_info structure, they should check
         * their bucket allocation */
        double max_bucket_fill_ratio[FB_MAX_PARTS];

        /* This is updated by applying the special-q lattice transform to
         * the factor base. */
        small_sieve_data_t ssd;

    };
    /* }}} */

    side_info sides[2];

    /* Most of the member functions below which take a side argument
     * should be members of the side_info structure instead.
     *
     * Note that have have only init_* functions below, no clear_*. In
     * fact, it's mostly a misnomer, and I apologize for it. The
     * automatic dtor will work fine for all these fields, hence there is
     * no need for specific clear_* functions. The magic lies most of the
     * time in the use of shared_ptr's. Maybe the init_* ones should have
     * been named setup_* or something like this.
     */

    /* in las-fb.cpp */
    void init_factor_bases( param_list_ptr pl);
    void share_factor_bases(sieve_info & other);

    void init_fb_smallsieved(int side);

    /* in las-trialdiv.cpp */
    void init_trialdiv(int side);

    /* in las-unsieve.cpp */
    /* Data for unsieving locations where gcd(i,j) > 1 */
    /* This gets initialized only when I is finally decided */
    unsieve_data us;
    void init_unsieve_data() { us = unsieve_data(conf.logI_adjusted, conf.logA); }

    /* in las-unsieve.cpp */
    /* Data for divisibility tests p|i in lines where p|j */
    /* This gets initialized only when J is finally decided. */
    j_divisibility_helper j_div;
    void init_j_div() { j_div = j_divisibility_helper(J); }

    /* in las-cofactor.cpp */
    std::shared_ptr<facul_strategies_t> strategies;
    void init_strategies(param_list_ptr pl);

    /* These functions must be called before actually sieving */
    void update_norm_data ();
    void update (unsigned int nr_workspaces);

    sieve_info(siever_config const & sc, cado_poly_srcptr cpoly, std::list<sieve_info> & sievers, cxx_param_list & pl, bool try_fbc = false);

    sieve_info() {
        cpoly = NULL;
        I = J = 0;
        memset(nb_buckets, 0, sizeof(nb_buckets));
        toplevel = 0;
    }

    /* This constructor is meant to cover the case where we want to store
     * all that to a common registry. It would be cleaner if the registry
     * were a hidden feature of the class, buried in las-types.cpp. The
     * danger is with the static initialization order fiasco which lurks
     * near, so I've been refraining so far.
     */
    static sieve_info & get_sieve_info_from_config(siever_config const & sc, cado_poly_srcptr cpoly, std::list<sieve_info> & registry, cxx_param_list & pl, bool try_fbc = false);

    void recover_per_sq_values(sieve_range_adjust const& Adj) {
        doing = Adj.doing;
        qbasis = Adj.Q;
        qbasis.set_q(doing.p, doing.prime_sq);
        if (!qbasis.prime_sq) {
            qbasis.prime_factors = doing.prime_factors;
        }
        ASSERT_ALWAYS(conf.logI_adjusted == Adj.logI);
        ASSERT_ALWAYS(I == (1UL << Adj.logI));
        J = Adj.J;
    }
};

/*  */

#endif	/* LAS_SIEVE_INFO_HPP_ */
