#ifndef LAS_SMALLSIEVE_HPP_
#define LAS_SMALLSIEVE_HPP_

#include <stdarg.h>
#include <vector>
#include "fb-types.h"
#include "macros.h"

/* Structures for small sieve */

#define SSP_POW2        (1u<<0)
#define SSP_PROJ        (1u<<1)
#define SSP_DISCARD_SUBLAT     (1u<<2)
#define SSP_DISCARD_PROJ     (1u<<3)

/* spos_t is in las-forwardtypes.hpp */

class ssp_simple_t {
protected:
    fbprime_t p;
    fbprime_t r;
    fbprime_t offset;
public:
    unsigned char logp;

    ssp_simple_t() : p(0), r(0), offset(0), logp(0) {}
    ssp_simple_t(fbprime_t _p, fbprime_t _r, unsigned char _logp, unsigned int skip)
    : p(_p), r(_r), offset((r*skip)%p), logp(_logp)
    {}
    fbprime_t get_p() const {return p;}
    fbprime_t get_r() const {return r;}
    fbprime_t get_offset() const {return offset;}
    void set_p(const fbprime_t _p) {p = _p;}
    void set_r(const fbprime_t _r) {r = _r;}
    void set_offset(const fbprime_t _offset) {offset = _offset;}
    void print(FILE *) const;
};

class ssp_t : public ssp_simple_t {
    /* ordinary primes. Equation is (i-r*j) = 0 mod p */
    /* p may be a prime power, and even a power of two */
    /* Note that we need three fields for the projective primes anyway,
     * so we may have three here. Yet, admittedly the "offset" field is
     * rather useless. (see ssp_init_oa.)  We can recompute it on the fly
     * when needed.
     */
    unsigned char flags = 0;
public:

    /* Initialization procedures for the ssp data */
    /* Constructor for affine case */
    ssp_t() : ssp_simple_t(), flags(0) {}
    ssp_t(fbprime_t _p, fbprime_t _r, unsigned char _logp, unsigned int skip)
    : ssp_simple_t(_p, _r, _logp, skip)
    {}
    /* Constructor for affine or projective case */
    ssp_t(fbprime_t _p, fbprime_t _r, unsigned char _logp, unsigned int skip, bool proj)
    {
      if (proj || _p % 2 == 0) {
        init_proj(_p, _r, _logp, skip);
      } else {
        /* TODO: How to defer to the other constructor correctly? This here is bad */
        *this = ssp_t(_p, _r, _logp, skip);
      }
    }

    /* We could use the parent class' methods if we get rid of the ASSERT()s */
    fbprime_t get_p() const {ASSERT(!is_proj()); return p;}
    fbprime_t get_r() const {ASSERT(!is_proj()); return r;}
    fbprime_t get_offset() const {ASSERT(!is_proj()); return offset;}

    fbprime_t get_q() const {ASSERT(is_proj()); ASSERT(p > 0); return p;}
    fbprime_t get_g() const {ASSERT(is_proj()); ASSERT(r > 0); return r;}
    fbprime_t get_U() const {ASSERT(is_proj()); return offset;}

    void set_q(const fbprime_t q) {ASSERT(is_proj()); p = q;}
    void set_g(const fbprime_t g) {ASSERT(is_proj()); ASSERT(g > 0); r = g;}
    void set_U(const fbprime_t U) {ASSERT(is_proj()); offset = U;}

    bool is_pow2() const {return (flags & SSP_POW2) != 0;}
    bool is_proj() const {return (flags & SSP_PROJ) != 0;}
    bool is_discarded_sublat() const {return (flags & SSP_DISCARD_SUBLAT) != 0;}
    bool is_discarded_proj() const {return (flags & SSP_DISCARD_PROJ) != 0;}
    bool is_discarded() const {return is_discarded_proj() || is_discarded_sublat();}
    bool is_nice() const {return !is_pow2() && !is_proj() && !is_discarded();}

    void set_pow2() {flags |= SSP_POW2;}
    void set_proj() {flags |= SSP_PROJ;}
    void set_discarded_sublat() {flags |= SSP_DISCARD_SUBLAT;}
    void set_discarded() {flags |= SSP_DISCARD_PROJ;}
    
    void print(FILE *) const;
private:
    void init_proj(fbprime_t p, fbprime_t r, unsigned char _logp,
                   unsigned int skip MAYBE_UNUSED);
};

typedef struct {
    std::vector<ssp_simple_t> ssps;
    std::vector<ssp_t> ssp;
    /* These iterators tell which of the ssps entries should be used 
       for re-sieving */
    std::vector<ssp_simple_t>::const_iterator resieve_start, resieve_end;
} small_sieve_data_t;

/* Include this only now, as there's a cross dependency between the two
 * (our prototypes need sieve_info_t, which needs our datatypes...)
 */
#include "las-types.hpp"
#include "bucket.hpp"

extern void small_sieve_info(const char * what, int side, small_sieve_data_t const & r);
extern int small_sieve_dump(FILE *, const char *, va_list);
extern void small_sieve_clear(small_sieve_data_t & ssd);
extern void small_sieve_extract_interval(small_sieve_data_t & r, small_sieve_data_t const & s, int bounds[2]);
extern void small_sieve_init(small_sieve_data_t & ssd, unsigned int interleaving,
                      std::vector<fb_general_entry>::const_iterator fb_start,
                      std::vector<fb_general_entry>::const_iterator fb_end,
                      std::vector<fb_general_entry>::const_iterator resieve_start,
                      std::vector<fb_general_entry>::const_iterator resieve_end,
                      sieve_info const & si, int side);
extern void small_sieve_start(std::vector<spos_t> & ssdpos, small_sieve_data_t & ssd, unsigned int first_region_index, sieve_info const & si);
extern void small_sieve_copy_start(std::vector<spos_t>& res, std::vector<spos_t> const & base, int bounds[2]);
/*
extern void small_sieve_skip_stride(small_sieve_data_t *ssd, spos_t * ssdpos, unsigned int N, unsigned int interleaving, sieve_info const & si);
*/
extern void sieve_small_bucket_region(unsigned char *S, unsigned int N,
                               small_sieve_data_t & ssd,
                               std::vector<spos_t> & ssdpos,
                               sieve_info & si, int side MAYBE_UNUSED,
                               int nthreads,
                               where_am_I & w);

extern void
resieve_small_bucket_region (bucket_primes_t *BP, int N, unsigned char *S,
        small_sieve_data_t & ssd, 
        std::vector<spos_t> & ssdpos,
        sieve_info const & si,
        int interleaving,
        where_am_I& w MAYBE_UNUSED);


#endif	/* LAS_SMALLSIEVE_HPP_ */
