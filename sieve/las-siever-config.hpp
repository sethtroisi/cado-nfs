#ifndef LAS_SIEVER_CONFIG_HPP_
#define LAS_SIEVER_CONFIG_HPP_

#include "las-base.hpp"
#include "las-todo-entry.hpp"
#include "fb-types.h"
#include <string.h>
#include <string>
#include <map>
#include "params.h"

/* siever_config */
 
class bkmult_specifier {
    double base = 1.0;
    typedef std::map<std::pair<int, char>, double> dict_t;
    dict_t dict;
    public:
    typedef dict_t::key_type key_type;
    static std::string printkey(dict_t::key_type const& key) {
        char c[3] = { (char) ('0' + key.first), key.second, '\0' };
        return std::string(c);
    }
    template<typename T> static dict_t::key_type getkey() {
        return dict_t::key_type(T::level(), T::rtti[0]);
    }
    template<typename T> double get() const { return get(getkey<T>()); }
    double const & get(dict_t::key_type const& key) const {
        auto xx = dict.find(key);
        if (xx != dict.end()) return xx->second;
        return base;
    }
    double grow(dict_t::key_type const& key, double d) {
        double v = get(key) * d;
        return dict[key] = v;
    }
    template<typename T> double get(T const &) const { return get<T>(); }
    template<typename T> double operator()(T const &) const { return get<T>(); }
    template<typename T> double operator()() const { return get<T>(); }
    bkmult_specifier(double x) : base(x) {}
    bkmult_specifier(const char * specifier);
    std::string print_all() const;
};
/* The following structure lists the fields with an impact on the siever.
 * Different values for these fields will correspond to different siever
 * structures.
 */
struct siever_config {
    /* The bit size of the special-q. Counting in bits is no necessity,
     * we could imagine being more accurate */
    unsigned int bitsize=0;  /* bitsize == 0 indicates end of table */
    int side=0;   /* special-q side */
    int logA=0;


    /* For a given logA, we may trigger configurations for various logI
     * values. logI_adjusted is a sieving-only parameter. */
    int logI_adjusted=0;
    sublat_t sublat;
    unsigned long bucket_thresh=0;    // bucket sieve primes >= bucket_thresh
    unsigned long bucket_thresh1=0;   // primes above are 2-level bucket-sieved
    unsigned int td_thresh=0;
    unsigned int skipped=0;           // don't sieve below this
    bkmult_specifier bk_multiplier { 1.0 };     // how much margin when allocating buckets
    unsigned int unsieve_thresh=0;
    struct side_config {
        unsigned long lim=0;    /* factor base bound */
        unsigned long powlim=0; /* bound on powers in the factor base */
        int lpb=0;              /* large prime bound is 2^lpb */
        int mfb=0;              /* bound for residuals is 2^mfb */
        int ncurves=0;          /* number of cofactorization curves */
        double lambda=0;        /* lambda sieve parameter */
        unsigned long qmin=0;   /* smallest q sieved on this side, for dupsup */
        unsigned long qmax=0;   /* largest q sieved on this side, for dupsup */
    };
    side_config sides[2];

    void display() const;

    static void declare_usage(param_list_ptr pl);
    static bool parse_default(siever_config & sc, param_list_ptr pl);

    /*{{{ has_same_config */
    bool operator==(siever_config const & o) const { return memcmp(this, &o, sizeof(*this)) == 0; }

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


/* {{{ descent_hint
 *
 * This is used for the descent. For each factor size, we provide a
 * reasonable siever_config value
 *
 * We also provide, based on experience, info relative to how long it
 * takes to finish the smoothing process for a prime factor of this size.
 */
/* }}} */

struct siever_config_pool {
    typedef std::pair<int, unsigned int> key_type;

    struct descent_hint : public siever_config {
        double expected_time;
        double expected_success;
    };

    typedef std::map<key_type, descent_hint> hint_table_t;
    hint_table_t hints;

    descent_hint const * get_hint(int side, unsigned int bitsize) const {
        hint_table_t::const_iterator it = hints.find(key_type(side, bitsize));
        if (it == hints.end())
            return NULL;
        else
            return &it->second;
    }

    siever_config const * default_config_ptr;

    /* This needs not be complete. The default_config_ptr field points
     * here if it is complete. If not, the fields here are just used as a
     * base for initializing the other configurations */
    siever_config base;

    siever_config get_config_for_q(las_todo_entry const& doing) const;

    void grow_bk_multiplier(bkmult_specifier::key_type const& key, double d) {
        for(auto & c : hints)
            c.second.bk_multiplier.grow(key, d);
        base.bk_multiplier.grow(key, d);
    }


    siever_config_pool(cxx_param_list& pl);

    double hint_expected_time(key_type const &K) const {
        if (hints.find(K) == hints.end())
            return -1;
        return hints.at(K).expected_time;
    }
    double hint_expected_success(key_type const &K) const {
        if (hints.find(K) == hints.end())
            return -1;
        return hints.at(K).expected_success;
    }
};

#endif	/* LAS_SIEVER_CONFIG_HPP_ */
