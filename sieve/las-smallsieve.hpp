#ifndef LAS_SMALLSIEVE_HPP_
#define LAS_SMALLSIEVE_HPP_

#include <stdarg.h>
#include "fb-types.h"

/* Structures for small sieve */

/* the two structures here are in the small sieve array. It's written as
 * "ordinary" versus "bad", but in reality one should rather see it as
 * "affine" versus "projective".
 */
class ssp_t {
    /* ordinary primes. Equation is (i-r*j) = 0 mod p */
    /* p may be a prime power, and even a power of two */
    /* Note that we need three fields for the projective primes anyway,
     * so we may have three here. Yet, admittedly the "offset" field is
     * rather useless. (see ssp_init_oa.)  We can recompute it on the fly
     * when needed.
     */
public:
    fbprime_t p;	// q if projective
    fbprime_t r;        // in [ 0, p [. g if projective
    fbprime_t offset;   // in [ 0, p [. U if projective

    fbprime_t get_q() const {return p;}
    fbprime_t get_g() const {return r;}
    fbprime_t get_U() const {return offset;}
    void set_q(fbprime_t q) {p=q;}
    void set_g(fbprime_t g) {r=g;}
    void set_U(fbprime_t U) {offset=U;}
};

#define SSP_POW2        (1u<<0)
#define SSP_PROJ        (1u<<1)
#define SSP_DISCARD_SUBLAT     (1u<<29)
#define SSP_DISCARD     (1u<<30)
#define SSP_END         (1u<<31)

typedef struct {
    int index;
    unsigned int event;
} ssp_marker_t;

typedef struct {
    ssp_marker_t * markers;
    ssp_t *ssp;
    int nb_ssp;
    unsigned char * logp;
} small_sieve_data_t;

/* Include this only now, as there's a cross dependency between the two
 * (our prototypes need sieve_info_t, which needs our datatypes...)
 */
#include "las-types.hpp"
#include "bucket.hpp"

extern void small_sieve_info(const char * what, int side, small_sieve_data_t * r);
extern int small_sieve_dump(FILE *, const char *, va_list);
extern void small_sieve_clear(small_sieve_data_t * ssd);
extern void small_sieve_extract_interval(small_sieve_data_t * r, small_sieve_data_t * s, int bounds[2]);
extern void small_sieve_init(small_sieve_data_t *ssd, unsigned int interleaving, const std::vector<fb_general_entry> *fb,
                      sieve_info const & si, int side);
extern int64_t * small_sieve_copy_start(int64_t * base, int bounds[2]);
extern int64_t * small_sieve_start(small_sieve_data_t *ssd, unsigned int N, sieve_info const & si);
/*
extern void small_sieve_skip_stride(small_sieve_data_t *ssd, int64_t * ssdpos, unsigned int N, unsigned int interleaving, sieve_info const & si);
*/
extern void sieve_small_bucket_region(unsigned char *S, int N,
			       small_sieve_data_t * ssd, int64_t * ssdpos, sieve_info& si,
                               int side,
                               int interleaving,
			       where_am_I& w MAYBE_UNUSED);

extern void
resieve_small_bucket_region (bucket_primes_t *BP, int N, unsigned char *S,
        small_sieve_data_t *ssd, int64_t * ssdpos,
        sieve_info const & si,
        int interleaving,
        where_am_I& w MAYBE_UNUSED);


#endif	/* LAS_SMALLSIEVE_HPP_ */
