#ifndef LAS_SIEVER_CONFIG_HPP_
#define LAS_SIEVER_CONFIG_HPP_

#include "las-base.hpp"
#include "las-todo-entry.hpp"
#include "fb-types.h"
#include <string.h>
#include <map>
#include "params.h"

/* {{{ siever_config */
/* The following structure lists the fields with an impact on the siever.
 * Different values for these fields will correspond to different siever
 * structures.
 */
struct siever_config : public _padded_pod<siever_config> {
    /* The bit size of the special-q. Counting in bits is no necessity,
     * we could imagine being more accurate */
    unsigned int bitsize;  /* bitsize == 0 indicates end of table */
    int side;   /* special-q side */
    int logA;


    /* For a given logA, we may trigger configurations for various logI
     * values. logI_adjusted is a sieving-only parameter. */
    int logI_adjusted;
    sublat_t sublat;

    unsigned long bucket_thresh;    // bucket sieve primes >= bucket_thresh
    unsigned long bucket_thresh1;   // primes above are 2-level bucket-sieved
    unsigned int td_thresh;
    unsigned int skipped;           // don't sieve below this
    double bk_multiplier;           // how much margin when allocating buckets
    unsigned int unsieve_thresh;
    struct side_config {
        unsigned long lim; /* factor base bound */
        unsigned long powlim; /* bound on powers in the factor base */
        int lpb;           /* large prime bound is 2^lpb */
        int mfb;           /* bound for residuals is 2^mfb */
        int ncurves;       /* number of cofactorization curves */
        double lambda;     /* lambda sieve parameter */
        unsigned long qmin; /* smallest q sieved on this side, for dup sup */
        unsigned long qmax; /* largest q sieved on this side, for dup sup */
    };
    side_config sides[2];

    bool operator==(siever_config const & o) const { return memcmp(this, &o, sizeof(*this)) == 0; }

    /*{{{ has_same_config */
    struct has_same_config {
        siever_config const & sc;
        has_same_config(siever_config const & sc) : sc(sc) {}
        bool operator()(siever_config const& o) const { return o == sc; }
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
    };
    has_same_config same_config() const { return has_same_config(*this); }
    /*}}}*/
    /*{{{ has_same_config_q */
    struct has_same_config_q {
        siever_config const & sc;
        has_same_config_q(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const {
            return sc.side == o.side && sc.bitsize == o.bitsize;
        }
    };
    has_same_config_q same_config_q() const {
        return has_same_config_q(*this);
    }
    /*}}}*/
    /* {{{ has_same_fb_parameters */
    struct has_same_fb_parameters {
        siever_config const & sc;
        has_same_fb_parameters(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const {
            bool ok = true;
            // ok = ok && sc.logI_adjusted == o.logI_adjusted;
            ok = ok && sc.bucket_thresh == o.bucket_thresh;
            ok = ok && sc.bucket_thresh1 == o.bucket_thresh1;
            ok = ok && sc.td_thresh == o.td_thresh;
            ok = ok && sc.skipped == o.skipped;
            // ok = ok && sc.bk_multiplier == o.bk_multiplier;
            ok = ok && sc.unsieve_thresh == o.unsieve_thresh;
            for(int side = 0 ; side < 2 ; side++) {
                ok = ok && sc.sides[side].lim == o.sides[side].lim;
                ok = ok && sc.sides[side].powlim == o.sides[side].powlim;
            }
            return ok;
        }
    };
    has_same_fb_parameters same_fb_parameters() const { return has_same_fb_parameters(*this); }
    /*}}}*/
    /*{{{ has_same_sieving -- currently duplicates has_same_fb_parameters */
    struct has_same_sieving {
        siever_config const & sc;
        has_same_sieving(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const {
            return has_same_fb_parameters(sc)(o);
        }
    };
    has_same_sieving same_sieving() const { return has_same_sieving(*this); }
    /*}}}*/
    /*{{{ has_same_cofactoring */
    struct has_same_cofactoring {
        siever_config const & sc;
        has_same_cofactoring(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const {
            bool ok = true;
            for(int side = 0 ; side < 2 ; side++) {
                ok = ok && sc.sides[side].lambda == o.sides[side].lambda;
                ok = ok && sc.sides[side].lpb == o.sides[side].lpb;
                ok = ok && sc.sides[side].mfb == o.sides[side].mfb;
                ok = ok && sc.sides[side].ncurves == o.sides[side].ncurves;
            }
            return ok;
        }
    };
    has_same_cofactoring same_cofactoring() const { return has_same_cofactoring(*this); }
    /*}}}*/
};

/* }}} */

void siever_config_display(siever_config const & sc);

#endif	/* LAS_SIEVER_CONFIG_HPP_ */
