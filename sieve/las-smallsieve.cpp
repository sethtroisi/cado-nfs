#include "cado.h"
#include <pthread.h>
#if defined(HAVE_SSE2)
#include <emmintrin.h>
#endif
#include "las-config.h"
#include "las-smallsieve.hpp"
#include "las-debug.hpp"
#include "las-qlattice.hpp"
#include "misc.h"
#include "portability.h"
#include "verbose.h"
#include "las-smallsieve-glue.hpp"


/* small sieve and resieving */

/* {{{ Some documentation first.
 *
 * Small primes or powers of small primes p^k with projective root.
 * These hit at 
 *   i*v == j*u (mod p^k) 
 * for some u,v in Z, but gcd(v, p^k) > 1.
 * We may assume gcd(u,p)==1, or we divide the entire equation by p.
 * XXX [ET]: we should also assume that v is a prime power, and that u
 * XXX [ET]: is within [0..p^k/v-1[ ; 
 * We store g = gcd(v, p^k), q = p^k / g, and U = u * (v/g)^(-1) (mod q).
 * XXX [ET]: which would then imply g==v, q=p^k/v, and U=u

 * Then we have
 *   i*v == j*u (mod p^k)  <==>  i == (j/g)*U (mod q)
 * with g|j. 
 * 
 * In other words, we can sieve this projective prime (power) much like a 
 * normal prime (power) q with root U, except that after sieving a line 
 * we don't advance by one line, but by g lines.
 * The case where g = q^k and thus q = 1 can be sieved more efficiently,
 * of course, since every entry in each g-th line will be hit, so that
 * the sieving should use long word transfers.

 * Just like for normal primes, the ssdpos value points at the first
 * position to sieve relative to the start of the current sieve region.

 * Within a line that starts at index line_start, for array element of 
 * index x, we have x - line_start = i+I/2. 
 * We skip j=0, as it contains only the single possible relation 
 * (i,j) = (1,0). 
 * For j=1*g, we want i=U (mod q), so x - line_start == I/2+U (mod q),
 * so we initialise 
 *   ssdpos = I*g + (I/2 + U) % q
 * to get the first array index in line j=g, 
 * then within a line sieve ssdpos + t*q < I, t in N,
 * and update 
 *   ssdpos = (ssdpos - line_start + U) % q + line_start + g*I 
 * to get the first position to sieve in the next suitable line.
 * }}} */

/* {{{ Some code for information purposes only */

void ssp_simple_t::print(FILE *f) const
{
    fprintf(f, "# p = %" FBPRIME_FORMAT ", r = %" FBROOT_FORMAT ", offset = %" FBPRIME_FORMAT ", logp = %hhu",
        this->p, this->r, this->offset, this->logp);
}

void ssp_t::print(FILE *f) const
{
    fprintf(f, "# p = %" FBPRIME_FORMAT ", r = %" FBROOT_FORMAT ", offset = %" FBPRIME_FORMAT ", logp = %hhu",
        this->p, this->r, this->offset, this->logp);
    if (this->is_pow2()) fprintf(f, " (power of 2)");
    if (this->is_ordinary3()) fprintf(f, " (ordinary p=3)");
    if (this->is_proj()) fprintf(f, " (projective root)");
    if (this->is_discarded_proj()) fprintf(f, "(discarded) because of projective root");
    if (this->is_discarded_sublat()) fprintf(f, "(discarded because not compatible with sub lattices)");
}

int
small_sieve_dump(FILE *f, const char *header, va_list va) /// XXX uh ? va_ ?
{
    const small_sieve_data_t * p_ssd = va_arg(va, const small_sieve_data_t *);

    fprintf(f, "%s", header);
    for(auto const & x : p_ssd->ssps) {
        x.print(f);
        fprintf(f, "\n");
    }
    for (auto const & x : p_ssd->ssp) {
        x.print(f);
        fprintf(f, "\n");
    }
    return 1;
}

static void small_sieve_print_contents(const char * prefix, small_sieve_data_t const & ssd)
{
    int nnice=ssd.ssps.size();
    int nproj=0;
    int npow2=0;
    int nordinary3=0;
    int ndiscard=0;
    for(auto const & ssp : ssd.ssp) {
        nproj += ssp.is_proj();
        npow2 += ssp.is_pow2();
        nordinary3 += ssp.is_ordinary3();
        ndiscard += ssp.is_discarded();
        nnice += ssp.is_nice();
    }

    verbose_output_start_batch();
    verbose_output_print(0, 2, "# %s: %d nice primes", prefix, nnice);
    /* Primes may be both even and projective... */
    if (npow2) verbose_output_print(0, 2, ", %d powers of 2", npow2);
    if (nordinary3) verbose_output_print(0, 2, ", %d pattern-sieved 3's", nordinary3);
    if (nproj) verbose_output_print(0, 2, ", and %d projective primes", nproj);
    verbose_output_print(0, 2, ".");
    if (ndiscard) verbose_output_print(0, 2, " %d discarded.", ndiscard);
    verbose_output_print(0, 2, "\n");
    /* With -v -v -v, dump all the small sieve data */
    verbose_output_vfprint (0, 4, small_sieve_dump, "# Dump of small sieve data:\n", &ssd);
    verbose_output_end_batch();
}


void small_sieve_info(const char * what, int side, small_sieve_data_t const & r)
{
    char * tmp;
    int rc = asprintf(&tmp, "%s(side %d)", what, side);
    ASSERT_ALWAYS(rc >= 0);
    small_sieve_print_contents(tmp, r);
    free(tmp);
}

/* }}} */

/* {{{ Sieve initialization / clearing : first the easy ones */
void small_sieve_clear(small_sieve_data_t & ssd)
{
    ssd.ssps.clear();
    ssd.ssp.clear();
}

void small_sieve_extract_interval(small_sieve_data_t & r, small_sieve_data_t const & s, int bounds[2])
{
    r = small_sieve_data_t();
    r.ssp.assign(s.ssp.begin() +  bounds[0], s.ssp.begin() + bounds[1]);
}

/* }}} */

/* {{{ Sieve initialization: now the real stuff */

void ssp_t::init_proj(fbprime_t p, fbprime_t r, unsigned char _logp, unsigned int skip MAYBE_UNUSED)/*{{{*/
{
    set_proj();
    if (p % 2 == 0)
      set_pow2();
    logp = _logp;
    unsigned int v = r; /* have consistent notations */
    unsigned int g = gcd_ul(p, v);
    fbprime_t q = p / g;
    set_q(q);
    set_g(g);
    if (q == 1) {
        ASSERT(r == 0);
        set_U(0);
    } else {
        int rc;
        uint32_t U = v / g; /* coprime to q */
        rc = invmod_32(&U, q);
        ASSERT_ALWAYS(rc != 0);
        set_U(U);
    }
}/*}}}*/

// Prepare sieving of small primes: initialize a small_sieve_data_t
// structure to be used thereafter during sieving each region.
// ssdpos points at the next position that will be hit by sieving,
// relative to the start of the next bucket region to sieve. It may exceed I 
// and even BUCKET_REGION
// It could actually be larger than 32 bits when I > 16.

void small_sieve_init(small_sieve_data_t & ssd,
                      unsigned int nthreads,
                      std::vector<fb_general_entry>::const_iterator fb_start,
                      std::vector<fb_general_entry>::const_iterator fb_end,
                      std::vector<fb_general_entry>::const_iterator resieve_start,
                      std::vector<fb_general_entry>::const_iterator resieve_end,
                      sieve_info const & si, const int side)
{
    const unsigned int thresh = si.conf.bucket_thresh;
    const int verbose = 0;
    where_am_I w MAYBE_UNUSED;

    ssd.ssp.clear();
    ssd.ssps.clear();

    // Do a pass on fb and projective primes, to fill in the data
    // while we have any regular primes or projective primes < thresh left
    
    // The processing of bucket region by nb_threads is interleaved.
    // It means that the positions for the small sieve must jump
    // over the (nb_threads - 1) regions after each region.
    // For typical primes, this jump is easily precomputed and goes into
    // the ssp struct.
    
    // If we are doing sublattices modulo m, then we jump virutally m
    // times faster.
    unsigned int sublatm = si.conf.sublat.m;
    const unsigned int skiprows = ((nthreads-1) << LOG_BUCKET_REGION) >> si.conf.logI_adjusted;
    bool saw_resieve_start = false, saw_resieve_end = false;
    
    for (std::vector<fb_general_entry>::const_iterator iter = fb_start;
         iter != fb_end; iter++) {
        /* p=pp^k, the prime or prime power in this entry, and pp is prime */
        const fbprime_t p = iter->q, pp = iter->p;
        WHERE_AM_I_UPDATE(w, p, p);

        if (iter == resieve_start) {
            ASSERT_ALWAYS(!saw_resieve_end);
            saw_resieve_start = true;
            ssd.resieve_start_offset = ssd.ssps.size();
        }
        if (iter == resieve_end) {
            ASSERT_ALWAYS(saw_resieve_start);
            saw_resieve_end = true;
            ssd.resieve_end_offset = ssd.ssps.size();
        }

        ASSERT_ALWAYS(p <= thresh);
        if (p > thresh) {
            continue;
        }

        for (int nr = 0; nr < iter->nr_roots; nr++) {
            const fb_general_root &root = iter->roots[nr];
            /* Convert into old format for projective roots by adding p if projective.
               FIXME, this sucks. */
            const fbroot_t r = root.r + (root.proj ? p : 0);

            if (!si.doing.is_coprime_to(pp)) {
                continue;
            }

            if (sublatm) {
                // In sublat mode, disable pattern sieving and primes
                // dividing m. (pp is the prime, here)
                //
                // FIXME. ok, they're certainly not "nice", but we should
                // sieve them nonetheless.
                if (pp == 2 || pp == 3 || (sublatm % pp) == 0) {
                    continue;
                }
            }

            const unsigned char logp = fb_log_delta (pp, root.exp, root.oldexp,
                si.sides[side].lognorms->scale);
            
            WHERE_AM_I_UPDATE(w, r, r);
            fbroot_t r_q = fb_root_in_qlattice(p, r, iter->invq, si.qbasis);
            /* If this root is somehow interesting (projective in (a,b) or
               in (i,j) plane), print a message */
            const bool is_proj_in_ij = r_q >= p;
            if (is_proj_in_ij) r_q -= p;
            if (verbose && (r > p || is_proj_in_ij))
                verbose_output_print(0, 1, "# small_sieve_init: side %d, prime %"
                        FBPRIME_FORMAT " root %s%" FBROOT_FORMAT " -> %s%" 
                        FBROOT_FORMAT "\n", side, p, 
                        r >= p ? "1/" : "", r % p,
                        is_proj_in_ij ? "1/" : "", r_q % p);

            ssp_t new_ssp(p, r_q, logp, skiprows, is_proj_in_ij);

            /* FIXME: what should we do with 3, given that it is pattern-sieved?
             */

            /* don't push projective primes that never hit. */
            if (new_ssp.is_proj()) {
                if (new_ssp.get_g() >= si.J) {
                    /* ... unless the number of lines to skip is >= J */
                    /* FIXME: we lose hits to (+-1,0) this way (the two locations
                       are equal up to sign, but we should sieve one of them!) */
                    if (verbose) {
                        verbose_output_print(0, 1,
                                "# small_sieve_init: not adding projective prime"
                                " (1:%" FBROOT_FORMAT ") mod %" FBPRIME_FORMAT ")"
                                " to small sieve  because g=%d >= si.J = %d\n",
                                r_q-p, p, new_ssp.get_g(), si.J);
                    }
                    continue;
                }
                /* projective primes of all sorts go to ssp anyway */
                ssd.ssp.push_back(new_ssp);
            } else if (new_ssp.is_ordinary3() || new_ssp.is_pow2()) {
                ssd.ssp.push_back(new_ssp);
            } else if (new_ssp.is_nice()) {
                ssd.ssps.push_back(new_ssp);
            } else {
                ASSERT_ALWAYS(0);
            }
        }
    }
    /* Now init_fb_smallsieved puts the resieved primes first anyway */
    ASSERT_ALWAYS(resieve_start == fb_start);
    if (resieve_start == fb_end) {
        ASSERT_ALWAYS(!saw_resieve_end);
        saw_resieve_start = true;
        ssd.resieve_start_offset = ssd.ssps.size();
    }
    if (resieve_end == fb_end) {
        ASSERT_ALWAYS(saw_resieve_start);
        saw_resieve_end = true;
        ssd.resieve_end_offset = ssd.ssps.size();
    }
    ASSERT_ALWAYS(saw_resieve_start && saw_resieve_end);

    /* arrange so that the small_sieve() ctor is happy */
    std::sort(ssd.ssps.begin(), ssd.ssps.begin() + ssd.resieve_end_offset);
    std::sort(ssd.ssps.begin() + ssd.resieve_end_offset, ssd.ssps.end());
}
/* }}} */

/* {{{ Creation of the ssdpos tables */

/*{{{ doc */
/* The places to be sieved are governed by the shape of the underlying
 * lattice of points.
 *
 * In the generic case, the lattice has determinant p (which may be a
 * prime power), and its basis in the (i,j) plane may be written as: (p
 * 0) (r 1) so that the lattice points are exactly the ones which satisfy
 * the equations i-r*j=0 mod p.  This means, in particular, that there
 * are hits on every line.
 *
 * With the extra assumption that p<=I and p<=BUCKET_REGION, this implies
 * that we'll have a starting point on the very first line of each bucket
 * region. However, versatility of the code may lead us to be slightly
 * tolerant about this assumption, and we'd like to see whether it is
 * possible to deal with larger p.
 *
 * The non-generic case is "projective". The lattice may be written as (q
 * 0) (U g) with q*g=p (which may be a prime power pp^k, and pp|g).
 *
 * What follows can be seen as a generalization of the generic case,
 * since it adapts by taking q=p and g=1.
 *
 * Note that when g>1, we no longer have a hit on every line. So the
 * first step is to try and see on which line we'll have a hit. Then the
 * (i,j) points will have to satisfy the equation i-(j/g)*U=0 mod q.
 *
 *
 * First task: the position of the first hit
 * -----------------------------------------
 *
 * We have a starting (current) bucket region with begin coordinates
 * (i0,j0). We want to find the first (i,j), >= (i0,j0) in processing
 * order, where p hits. Note that when this means j>j0, we may have i<i0.
 *
 * Specific for affine primes: in fact, because of the way we do the
 * sieve, we're rather going to compute the first i in the line j=j0,
 * relative to the current starting position i0, so that we have a hit.
 * The adjustments we make in the course of the computation will allow us
 * to catch up with the alignment constraints.
 *
 * Second task: actual sieve
 * -------------------------
 *
 * We have the first hit in ssdpos. We increment position by q (p in the
 * generic case), and update each point. The only difficult questions
 * are: - how do we deal with q even, given that we want to make sure we
 * don't sieve the (even,even) positions ?  - how do we deal, within the
 * sieve, with p reaching end of line ?  - when the position goes beyond
 * BUCKET_REGION, we need to move on to the next prime. We'll store this
 * "first off" position for further dealings with p. Do we need an
 * adjustment.  - the timing for end-of-bucket and end-of-line
 * adjustments is subtle and requires some care.
 *
 * Third task: reposition the sieve to jump over bucket regions
 * ------------------------------------------------------------
 *
 * In multithreaded mode, bucket regions are processed in an interleaved
 * fashion (it's not the only possible way to proceed). This means that
 * with N threads working, the position values per prime must be adjusted
 * after processing so that the N-1 bucket regions which have been sieved
 * by other threads are skipped.
 *
 * This adjustment may be made in a different manner depending on whether
 * we consider affine or projective primes.
 *
 *
 * end-of-bucket and end-of-line adjustments
 * -----------------------------------------
 *
 * This is a tricky part. When several bucket regions are in a line, once
 * we've marked hits in a bucket region and the sieve position has grown
 * past the end-of-bucket, the only needed adjustment is to subtract the
 * bucket region size from the position value: this will properly
 * represent the position, relative to the next bucket (if buckets are
 * processed in order). We do the same when there are no hits, so the
 * case "first hit is in the next bucket" does not need to be treated in
 * any specific way.
 *
 * On the other hand, at the end of a line, it's a bit different.
 * (i1-1)+j*I and i0+(j+1)*I are represented by consecutive bytes in
 * memory, so the coordinates (i1,j) and (i0,j+1) are equivalent
 * memory-wise. They're not equivalent with respect to the sieve
 * constraint that (X(i,j)=i-jr = 0 mod p). So for the "next position",
 * we have X(i,j+1)=-r if we don't pay attention. We fix that by adding r
 * to the sieved locations for line j+1, incorporating a reduction mod p
 * of the first sieved location.  This must be done for every line.
 *
 *
 */
/*}}}*/

/* This is copied from LOGNORM_FILL_COMMON_DEFS in las-norms.cpp ; from
 * logI, N, and LOG_BUCKET_REGION, define the integers i0, i1, j0, j1,
 * and I.
 */
#define SMALLSIEVE_COMMON_DEFS()                                         \
    const unsigned int log_lines_per_region = MAX(0, LOG_BUCKET_REGION - logI);\
    const unsigned int log_regions_per_line = MAX(0, logI - LOG_BUCKET_REGION);\
    const unsigned int regions_per_line = 1 << log_regions_per_line;           \
    const unsigned int region_rank_in_line = N & (regions_per_line - 1);       \
    const bool last_region_in_line MAYBE_UNUSED = region_rank_in_line == (regions_per_line - 1); \
    const unsigned int j0 = (N >> log_regions_per_line) << log_lines_per_region;    \
    const unsigned int j1 MAYBE_UNUSED = j0 + (1 << log_lines_per_region);    \
    const int I = 1 << logI;                                            \
    const int i0 = (region_rank_in_line << LOG_BUCKET_REGION) - I/2;          \
    const int i1 MAYBE_UNUSED = i0 + (1 << MIN(LOG_BUCKET_REGION, logI));     \
    /* those are (1,0,0) in the standard case */                        \
    const int sublatm MAYBE_UNUSED = si.conf.sublat.m ? si.conf.sublat.m : 1; \
    const unsigned int sublati0 MAYBE_UNUSED = si.conf.sublat.i0;       \
    const unsigned int sublatj0 MAYBE_UNUSED = si.conf.sublat.j0;       \
    const int row0_is_oddj MAYBE_UNUSED = (j0*sublatm + sublatj0) & 1;  \
    bool has_haxis = !j0;                                               \
    bool has_vaxis = region_rank_in_line == ((regions_per_line-1)/2);   \
    bool has_origin MAYBE_UNUSED = has_haxis && has_vaxis;              \
    do {} while (0)

/* Only compute the initial ssdpos fields. */
void small_sieve_start(std::vector<spos_t> & ssdpos,
        small_sieve_data_t & ssd,
        unsigned int first_region_index,
        sieve_info const & si)
{
#if 0
    int logI = si.conf.logI_adjusted;
    ssdpos.assign(ssd.ssp.size(), 0);

    bool is_fragment = logI > LOG_BUCKET_REGION;

    /* We want to compute the index of the "next" hit, counted from the
     * starting offset of the "current" bucket region at (i0,j0). The
     * next hit means that it has to be the first in the current bucket
     * region, or in the first bucket region that comes after it in
     * processing order.
     *
     * In processing order, the memory location for cell (i,j), relative
     * to (i0,j0), is (i-i0) + (j-j0)*I -- note that (i-i0) may be
     * negative if j>j0.
     *
     */

    if (is_fragment) {
        size_t index = 0;
        small_sieve_base<tribool_const<true> > C(si.conf.logI_adjusted, first_region_index, si.conf.sublat);
        for (std::vector<ssp_t>::const_iterator iter = ssd.ssp.begin();
             iter != ssd.ssp.end(); iter++) {
             ssp_t const & ssp = *iter;
            /* generic case only. For all other cases (which are rare enough
             * -- typically at most 20, counting powers of two and such), we
             *  compute the starting point from within
             *  sieve_small_bucket_region for each bucket region.
             */
            if (ssp.is_nice()) {
                ssdpos[index++] = C.first_position_ordinary_prime(ssp);
            } else if (ssp.is_proj()) {
                /* This also handles powers of 2 with projective root */
                ssdpos[index++] = C.first_position_projective_prime(ssp);
            } else if (ssp.is_pow2()) {
                /* Powers of 2 with affine root */
                ssdpos[index++] = C.first_position_power_of_two(ssp);
            } else {
                abort(); /* How did we get here? */
            }
        }
    } else {
        size_t index = 0;
        small_sieve_base<tribool_const<true> > C(si.conf.logI_adjusted, first_region_index, si.conf.sublat);
        for (ssp_t const & ssp : ssd.ssp) {
             if (ssp.is_nice())
                 ssdpos[index] = C.first_position_ordinary_prime(ssp);
             else if (ssp.is_proj())
                 ssdpos[index] = C.first_position_projective_prime(ssp);
             else if (ssp.is_pow2())
                 ssdpos[index] = C.first_position_power_of_two(ssp);
             else
                 abort();
             index++;
        }
    }
#endif

    /* We store start positions for the simple case only.
     * For all other cases (which are rare enough
     * -- typically at most 20, counting powers of two and such), we
     *  compute the starting point from within sieve_small_bucket_region
     *  for each bucket region.
     */
    small_sieve_base<> C(si.conf.logI_adjusted, first_region_index, si.conf.sublat);
    ssdpos.clear();
    ssdpos.reserve(ssd.ssps.size());
    for (ssp_simple_t const & ssp : ssd.ssps) {
        ssdpos.push_back(C.first_position_ordinary_prime(ssp));
    }
}
/* }}} */

/*{{{ ugliness that will go away soon. */
/* This adds extra logging for pattern sieving. Very slow.
 */
#define xxxUGLY_DEBUGGING

/* #define exactly one of these */
// #define xxxSMALLSIEVE_CRITICAL_UGLY_ASSEMBLY
// #define xxxSMALLSIEVE_CRITICAL_MANUAL_UNROLL
// #define SMALLSIEVE_CRITICAL_PLAIN

/* on an rsa155 test with I=14, the UGLY_ASSEMBLY version gains
 * practically nothing over the plain C version. The MANUAL_UNROLL
 * version is slower.
 *
 * On an rsa220 test with I=16, the UGLY_ASSEMBLY version gains about 3%
 * on the small sieve part strictly speaking, and then less than 1%
 * overall.
 */
#if 1
#if defined( HAVE_SSE2 ) && !defined( TRACK_CODE_PATH ) /* x86 optimized code */
#define SMALLSIEVE_CRITICAL_UGLY_ASSEMBLY
#else
#define SMALLSIEVE_CRITICAL_MANUAL_UNROLL
#endif
#else
#warning "using plain C code for critical part of small sieve"
#define SMALLSIEVE_CRITICAL_PLAIN
#endif


/* {{{ critical part of the small sieve */
#ifdef SMALLSIEVE_CRITICAL_UGLY_ASSEMBLY
      /* Some defines for the critical part of small sieve.
         0. The C code and the asm X86 code have the same algorithm.
         Read first the C code to understand easily the asm code.
         1. If there are less than 12 "T" in the line, the goal is to do
         only one jump (pipe-line breaks) and of course the instructions
         minimal number.
         2. 12 "T" seems the best in the critical loop. Before, gcc tries
         to optimize in a bad way the loop. For gcc generated code, the
         best is here a systematic code (12*2 instructions), like the C code.
         3. The asm X86 optimization uses addb logp,(pi+p_or_2p*[0,1,2])
         & three_p_or_2p = 3 * p_or_2p; the lea (add simulation) & real
         add interlace seems a bit interesting.
         So, the loop is smaller & faster (19 instructions versus 27 for
         gcc best X86 asm).
         4. Of course, the gain between the 2 versions is light, because
         the main problem is the access time of the L0 cache: read + write
         with sieve_increase(pi,logp,w), or *pi += logp in fact.
      */
#define U1                                                              \
        "addb %4,(%1)\n"                /* sieve_increase(pi,logp,w) */ \
        "addb %4,(%1,%3,1)\n"   /* sieve_increase(p_or_2p+pi,logp,w) */ \
        "addb %4,(%1,%3,2)\n" /* sieve_increase(p_or_2p*2+pi,logp,w) */ \
        "lea (%1,%2,1),%1\n"                    /* pi += p_or_2p * 3 */
#define U2                                                              \
        "cmp %5, %1\n"                    /* if (pi >= S1) break; */    \
        "jae 2f\n"                                                      \
        "addb %4,(%1)\n"                /* sieve_increase(pi,logp,w) */ \
        "lea (%1,%3,1),%1\n"                        /* pi += p_or_2p */
#define U do {                                                          \
        unsigned char *pi_end;                                          \
        size_t three_p_or_2p;                                           \
        __asm__ __volatile__ (                                          \
        "lea (%3,%3,2), %2\n"         /* three_p_or_2p = p_or_2p * 3 */ \
        "lea (%1,%2,4), %0\n"            /* pi_end = pi + p_or_2p*12 */ \
        "cmp %5, %0\n"                /* if (pi_end > S1) no loop */    \
        "jbe 0f\n"                                                      \
        "1:\n"                                                          \
        U2 U2 U2 U2 U2 U2 U2 U2 U2 U2 U2                                \
        "cmp %5, %1\n"                                                  \
        "jae 2f\n"                                                      \
        "addb %4,(%1)\n"                                                \
        "jmp 2f\n"                                                      \
        ".balign 8\n 0:\n"                              /* Main loop */ \
        U1 U1           /* sieve_increase(p_or_2p*[0..11]+pi,logp,w) */ \
        U1 U1                                    /* pi += p_or_2p*12 */ \
        "lea (%1,%2,4), %0\n"    /* if (pi+p_or_2p*12 > S1) break */    \
        "cmp %5, %0\n"                                                  \
        "jbe 0b\n"                                                      \
        "jmp 1b\n"                                                      \
        ".balign 8\n 2:\n"                                              \
        : "=&r"(pi_end), "+&r"(pi), "=&r"(three_p_or_2p)                \
        : "r"(p_or_2p), "q"(logp), "r"(S1) : "cc");                     \
        pi = pi_end;                                                    \
      } while (0)
#endif

#ifdef SMALLSIEVE_CRITICAL_MANUAL_UNROLL
#define T do {                                                          \
    WHERE_AM_I_UPDATE(w, x, x0 + pi - S0);                              \
    sieve_increase (pi, logp, w); pi += p_or_2p;                        \
} while(0)
#endif

inline size_t sieve_full_line(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#ifdef SMALLSIEVE_CRITICAL_UGLY_ASSEMBLY
    U;
#endif
#ifdef SMALLSIEVE_CRITICAL_MANUAL_UNROLL
    while (UNLIKELY(pi + p_or_2p * 12 <= S1))
    { T; T; T; T; T; T; T; T; T; T; T; T; }
    do {
        if (pi >= S1) break; T; if (pi >= S1) break; T;
        if (pi >= S1) break; T; if (pi >= S1) break; T;
        if (pi >= S1) break; T; if (pi >= S1) break; T;
        if (pi >= S1) break; T; if (pi >= S1) break; T;
        if (pi >= S1) break; T; if (pi >= S1) break; T;
        if (pi >= S1) break; T; if (pi >= S1) break; T;
    } while (0);
#endif
#ifdef SMALLSIEVE_CRITICAL_PLAIN
    for ( ; pi <= S1 ; pi += p_or_2p) {
        WHERE_AM_I_UPDATE(w, x, x0 + pi - S0);
        sieve_increase (pi, logp, w);
    }
#endif
    return pi - S1;
}

#ifdef SMALLSIEVE_CRITICAL_UGLY_ASSEMBLY
#undef U1
#undef U2
#undef U
#undef T
#endif
#ifdef SMALLSIEVE_CRITICAL_MANUAL_UNROLL
#undef T
#endif

/* }}} */
/*}}}*/


void sieve_one_nice_prime(unsigned char *S, const ssp_simple_t &ssps,
    int64_t & ssdpos, const sieve_info & si, const int i0, const int i1,
    const unsigned int j0, const unsigned int j1, const int I, const int logI,
    const int N, const int interleaving, where_am_I w)
{
    const fbprime_t p = ssps.get_p();
    const fbprime_t r = ssps.get_r();

    if (mpz_cmp_ui(si.qbasis.q, p) == 0)
        return;

    /* Don't sieve 3 again as it was pattern-sieved -- unless
     * it's projective, see TODO above. */
    if (p == 3)
      return;

    WHERE_AM_I_UPDATE(w, p, p);
    const unsigned char logp = ssps.logp;
    unsigned char *S0 = S;
    int pos = ssdpos;
#ifdef TRACE_K
    /* we're keeping track of it with ssdpos, but just in case,
     * make sure we get it right ! */
    small_sieve_base<> C(logI, N, si.conf.sublat);
    ASSERT_ALWAYS(pos == C.first_position_ordinary_prime (ssps));
#endif
    ASSERT(pos < (int) p);
    unsigned int i_compens_sublat = si.conf.sublat.i0 & 1;
    unsigned int j = j0;

    /* we sieve over the area [S0..S0+(i1-i0)], which may
     * actually be just a fragment of a line. After that, if
     * (i1-i0) is different from I, we'll break anyway. So
     * whether we add I or (i1-i0) to S0 does not matter much.
     */
    size_t overrun = 0; /* tame gcc */
    if (j & 1) {
        WHERE_AM_I_UPDATE(w, j, j - j0);
        overrun = sieve_full_line(S0, S0 + (i1 - i0), S0 - S,
                pos, p, logp, w);
        S0 += I;
        j++;
        pos += r; if (pos >= (int) p) pos -= p;
    }
    for( ; j < j1 ; ) {
        /* for j even, we sieve only odd pi, so step = 2p. */
        int xpos = ((i_compens_sublat + pos) & 1) ? pos : (pos+p);
        overrun = sieve_full_line(S0, S0 + (i1 - i0), S0 - S,
                xpos, p+p, logp, w);
        S0 += I;
        pos += r; if (pos >= (int) p) pos -= p; 
        if (++j >= j1) break;

        /* now j odd again */
        WHERE_AM_I_UPDATE(w, j, j - j0);
        overrun = sieve_full_line(S0, S0 + (i1 - i0), S0 - S,
                pos, p, logp, w);
        S0 += I;
        pos += r; if (pos >= (int) p) pos -= p;
        ++j;
    }

    if (logI > LOG_BUCKET_REGION) {
        /* quick notes for incremental adjustment in case I>B (B =
         * LOG_BUCKET_REGION).
         *
         * Let q = 2^(I-B).
         * Let N = a*q+b, and N'=N+interleaving=a'*q+b' ; N' is the
         * next bucket region we'll handle.
         *
         * Let interleaving = u*q+v
         *
         * The row increase is dj = (N' div q) - (N div q) = a'-a
         * The fragment increase is di = (N' mod q) - (N mod q) = b'-b
         *
         * Of course we have -q < b'-b < q
         *

         * dj can be written as (N'-N-(b'-b)) div q, which is an
         * exact division. we rewrite that as:
         *
         * dj = u + (v - (b'-b)) div q
         *
         * where the division is again exact. Now (v-(b'-b))
         * satisfies:
         * -q < v-(b'-b) < 2*q-1
         *
         * so that the quotient may only be 0 or 1.
         *
         * It is 1 if and only if v >= q + b'-b, which sounds like a
         * reasonable thing to check.
         */
        int N1 = N + interleaving;
        int Q = logI - LOG_BUCKET_REGION;
        int dj = (N1>>Q) - j0;
        int di = (N1&((1<<Q)-1)) - (N&((1<<Q)-1));
        /* Note that B_mod p is not reduced. It may be <0, and
         * may also be >= p if we sieved with 2p because of even j
         */
        int B_mod_p = overrun - pos;
        /* We may avoid some of the cost for the modular
         * reduction, here:
         *
         * dj is either always the same thing, or that same thing
         * + 1.
         *
         * di is within a small interval (albeit a centered one).
         * 
         * So it seems feasible to get by with a fixed number of
         * conditional subtractions.
         */
        pos = (pos + B_mod_p * di + dj * r) % p;
        if (pos < 0) pos += p;
#ifndef NDEBUG
        small_sieve_base<> C1(logI, N1, si.conf.sublat);
        ASSERT(pos == C1.first_position_ordinary_prime(ssps));
#endif
        ssdpos = pos;
    } else {
        /* skip stride */
        pos += ssps.get_offset();
        if (pos >= (int) p) pos -= p;
        ssdpos = pos;
    }
}


/* {{{ Pattern-sieve powers of 2 up to 2 * sizeof(long) */
template<bool is_fragment> void small_sieve<is_fragment>::pattern_sieve2(where_am_I & w MAYBE_UNUSED)
{
    /* TODO: use SSE2 (well, cost does not seem to be significant anyway).  */
    WHERE_AM_I_UPDATE(w, p, 2);
    /* First collect updates for powers of two in a pattern,
       then apply pattern to sieve line.
       Repeat for each line in bucket region. */
    for (unsigned int j = j0; j < j1; j++) {
        WHERE_AM_I_UPDATE(w, j, j - j0);
        unsigned long pattern[2] = {0,0};

        /* Prepare the pattern */

        /* Process only powers of 2. For now we do a full pass over the
         * ssp_t entries. It may seem annoying. However given that we're
         * now processing the lightweight ssp list and not the larger
         * list ssps that has all nice primes, it's perhaps not really
         * important to strive to minimize the cost of this processing
         * pass.
         */
        const bool verbose_pattern_2 = false;
        for(auto const & ssp : not_nice_primes) {
            if (!ssp.is_pow2())
                continue;
            if (ssp.is_proj()) {
                /* Do nothing. It's a projective power of 2, * but for the
                   moment these are not pattern-sieved. TODO: pattern sieve
                   them in those lines where they hit. */
                continue;
            }
            /* Powers of 2 greater than the pattern size get sieved below */
            if (ssp.get_p() > 2*sizeof(unsigned long))
                continue;
            /* *Affine* powers of two are relevant only for odd lines anyway:
             * indeed, if i-j*r=0 mod 2^k and j even, then i even too, so
             * useless.
             */
            if (j&1) {
                const fbprime_t p = ssp.get_p();
                /*
                if (mpz_cmp_ui(si.qbasis.q, p) == 0) {
                    continue;
                }
                */

                if (verbose_pattern_2) {
                    FILE *out;
                    verbose_output_start_batch();
                    for (size_t out_i=0; (out=verbose_output_get(TRACE_CHANNEL, 0, out_i)) != NULL; out_i++) {
                        fprintf(out, "# Adding ");
                        ssp.print(out);
                        fprintf(out, " to sieving pattern for powers of 2\n");
                    }
                    verbose_output_end_batch();
                }

                int pos = super::first_position_power_of_two(ssp, j - j0);

                if (pos < F()) {
                    // At the start of the bucket region, we are not sure
                    // that pos is reduced mod p when we enter here.
                    // (see comment at the end of small_sieve_skip_stride()).
                    // Hence, we have to do this reduction here:
                    pos &= (p-1);
                    ASSERT (pos < (spos_t) p);
                    ASSERT (j % 2);
                    for (int x = pos; x < (int) pattern2_size; x += p) {
                        if (verbose_pattern_2)
                            verbose_output_print(TRACE_CHANNEL, 0,
                                "# Hits at pattern[%d]\n", x);
                        ((unsigned char *)pattern)[x] += ssp.logp;
                    }
#ifdef UGLY_DEBUGGING
                    for (int x = pos; x < F() ; x+= p) {
                        WHERE_AM_I_UPDATE(w, x, (w.j << logI) + x);
                        sieve_increase(S + x, ssp.logp, w);
                        /* cancel the above action */
                        S[x] += ssp.logp;
                    }
#endif
#if 0
                    /* Skip two lines above, since we sieve only odd lines.
                     * Even lines would correspond to useless reports.
                     */
                    pos = ((pos + 2 * ssp.get_r()) & (p - 1)) + (2 << logI);
#endif
                }
#if 0
                /* In this loop, ssdpos gets updated to the first 
                   index to sieve relative to the start of the next line, 
                   but after all lines of this bucket region are processed, 
                   it will point to the first position to sieve relative
                   to the start of the next bucket region, as required */
                ssdpos[index] = pos - F();
#endif
            }
        }

        /* Apply the pattern */
        if (verbose_pattern_2)
            verbose_output_print(TRACE_CHANNEL, 0,
                "# Sieving pattern for powers of 2 is: %lx, %lx\n",
            pattern[0], pattern[1]);
        if (pattern[0] || pattern[1]) {
            unsigned long *S_ptr = (unsigned long *) (S + ((size_t) (j-j0) << logI));
            const unsigned long *S_end = (unsigned long *)(S + ((size_t) (j-j0) << logI) + F());

#ifdef TRACE_K /* {{{ */
            if (trace_on_range_Nx(w.N, w.j*I(), w.j*I()+I())) {
                unsigned int x = trace_Nx.x;
                unsigned int k = x % I();
                unsigned int v = (((unsigned char *)(pattern+((k/sizeof(unsigned long))&1)))[k%sizeof(unsigned long)]);
                if (v) {
                    WHERE_AM_I_UPDATE(w, x, x);
                    sieve_increase_logging(S + x, v, w);
                }
            }
#endif /* }}} */

            ASSERT((S_end - S_ptr) % 4 == 0);

            for ( ; S_ptr < S_end ; S_ptr += 4) {
                S_ptr[0] += pattern[0];
                S_ptr[1] += pattern[1];
                S_ptr[2] += pattern[0];
                S_ptr[3] += pattern[1];
            }
        }
    }
}
/* }}} */

/* {{{ Pattern-sieve 3 */
template<bool is_fragment> void small_sieve<is_fragment>::pattern_sieve3(where_am_I & w MAYBE_UNUSED)
{
    /* For the time being, it's really 3. But the only thing we care about is
     * the hit rate to be 1 every 3rd, on lines for which we hit. So the code
     * below could be improved to handle more cases, although presently it
     * does not. */
    const fbprime_t p = 3;
    WHERE_AM_I_UPDATE(w, p, p);
    /* First collect updates for powers of three in a pattern,
       then apply pattern to sieve line.
       Repeat for each line in bucket region. */
    for (unsigned int j = j0; j < j1; j++)
    {
        WHERE_AM_I_UPDATE(w, j, j - j0);
        unsigned long pattern[3];

        pattern[0] = pattern[1] = pattern[2] = 0UL;

        for(ssp_t const & ssp : not_nice_primes) {
            if (!ssp.is_ordinary3()) continue;

            ASSERT_ALWAYS(ssp.get_p() == 3);

            int pos = super::first_position_ordinary_prime(ssp, j - j0);

            ASSERT (pos < (spos_t) p);
            for (unsigned int x = pos; x < 3 * sizeof(unsigned long); x += p)
                ((unsigned char *)pattern)[x] += ssp.logp;
        }

        /* TODO sublat */
        if (!(j&1)) {
            /* We have an even j. There, we must not sieve even i's either ! */
            for (unsigned int x = 0; x < 3 * sizeof(unsigned long); x += 2)
                ((unsigned char *)pattern)[x] = 0;
        }

        /* We want to test if there is a non-zero entry in the pattern
         * within the first 6 entries (we sieve mod 3 and mod 2).
         * For compatibility with 32-bit machines, we test the first two
         * pattern unsigned longs */
        if (pattern[0] || pattern[1]) {
            unsigned long *S_ptr = (unsigned long *) (S + ((size_t) (j-j0) << logI));
            const unsigned long *S_end = (unsigned long *)(S + ((size_t) (j-j0) << logI) + F()) - 2;

#ifdef TRACE_K /* {{{ */
            if (trace_on_range_Nx(w.N, w.j*I(), w.j*I()+I())) {
                unsigned int x = trace_Nx.x;
                unsigned int k = x % I();
                unsigned int v = (((unsigned char *)(pattern+((k/sizeof(unsigned long))%3)))[k%sizeof(unsigned long)]);
                if (v) {
                    WHERE_AM_I_UPDATE(w, x, x);
                    sieve_increase_logging(S + x, v, w);
                }
            }
#endif /* }}} */

            for( ; S_ptr < S_end ; S_ptr += 3) {
                S_ptr[0] += pattern[0];
                S_ptr[1] += pattern[1];
                S_ptr[2] += pattern[2];
            }

            S_end += 2;
            if (S_ptr < S_end)
                *(S_ptr++) += pattern[0];
            if (S_ptr < S_end)
                *(S_ptr) += pattern[1];
        }
    }
} /* }}} */

/* {{{ Normal small sieve */
// Sieve small primes (up to p < bucket_thresh) of the factor base fb in the
// next sieve region S.
// Information about where we are is in ssd.
void sieve_small_bucket_region(unsigned char *S, unsigned int N,
                               small_sieve_data_t & ssd,
                               std::vector<spos_t> & ssdpos,
                               sieve_info & si, int side MAYBE_UNUSED,
                               int nthreads,
                               where_am_I & w)
{
    int logI = si.conf.logI_adjusted;
    bool is_fragment = logI > LOG_BUCKET_REGION;

#if 0
            if (fbprime_t(mpz_get_ui(si.qbasis.q)) == p) continue;

            /* Don't sieve 3 again as it was pattern-sieved -- unless
             * it's projective, see TODO above. */
            if (p == 3) continue;
#endif


    if (is_fragment) {
        small_sieve<true> SS(ssdpos, ssd.ssps, ssd.ssp, S, logI, N, si.conf.sublat, nthreads);
        SS.pattern_sieve2(w);
        SS.pattern_sieve3(w);
        SS.exceptional_sieve(w);
        SS.normal_sieve(w);
    } else {
        small_sieve<false> SS(ssdpos, ssd.ssps, ssd.ssp, S, logI, N, si.conf.sublat, nthreads);
        SS.pattern_sieve2(w);
        SS.pattern_sieve3(w);
        SS.exceptional_sieve(w);
        SS.normal_sieve(w);
    }
}
/*}}}*/


/* {{{ resieving. Different interface, since it plays with buckets as well.
 */

/* Sieve small primes (p < I) of the factor
   base fb in the next sieve region S, and add primes and the x position
   where they divide and where there's a sieve report to a bucket (rather
   than subtracting the log norm from S, as during sieving).
   Information about where we are is in ssd. */
void
resieve_small_bucket_region (bucket_primes_t *BP, int N, unsigned char *S,
        small_sieve_data_t & ssd, std::vector<spos_t> & ssdpos,
        sieve_info const & si, int nthreads MAYBE_UNUSED, where_am_I & w MAYBE_UNUSED)
{
    int logI = si.conf.logI_adjusted;
    SMALLSIEVE_COMMON_DEFS();
    small_sieve_base<> C(logI, N, si.conf.sublat);

    unsigned char *S_ptr;
    const int resieve_very_verbose = 0;

    unsigned int i_compens_sublat = sublati0 & 1;

    // Odd/even property of j is the same as for j+2, even with
    // sublat, unless sublat.m is even, which is not handled right
    // now. Same for i.
    ASSERT_ALWAYS(!sublatm || ((sublatm & 1) == 1));

    for(size_t index = ssd.resieve_start_offset ; index < ssd.resieve_end_offset ; index++) {
        auto const & ssps(ssd.ssps[index]);

        const fbprime_t p = ssps.get_p();
        fbprime_t r = ssps.get_r();
        WHERE_AM_I_UPDATE(w, p, p);
        int pos = ssdpos[index - ssd.resieve_start_offset];
        S_ptr = S;
        ASSERT(pos < (spos_t) p);
        /* for j even, we sieve only odd index. This translates into loops
         * which look as follows:
         *
         * j even: (sieve only odd index)
         *   for(index = pos + (p & -!(pos&1)) ; index < I ; index += p+p)
         *   (where (p & -!(pos&1)) is 0 if pos is odd, and p otherwise)
         * j odd: (sieve all values of index)
         *   for(index = pos                  ; index < I ; index += p)
         *
         * we may merge the two by setting q=p&-!((j&1)^row0_is_oddj)
         *
         * which, when (j+row0_is_oddj) is even, is p, and is 0
         * otherwise.
         *
         * In turn, since q changes for each j, 1 xor within the loop
         * is enough to make it alternate between 0 and p, once the
         * starting value is correct.
         *
         * TODO: ok, this is nice and good, but:
         *
         *  - I haven't seen this win for simple small sieve, so I
         *    doubt it's a good idea. Relying on the branch predictor
         *    or the compiler does not seem to be so stupid after all.
         *  - as present, the behaviour is obviously buggy for
         *    sublatm even.
         *  - we really want to have the same structure both for
         *    small sieve and resieving.
         */
        bool row0_even = (((j0&sublatm)+sublatj0) & 1) == 0;
        unsigned int q = row0_even ? p : 0;
        for (unsigned int j = j0; j < j1; j ++) {
            WHERE_AM_I_UPDATE(w, j, j);
            for (int i = pos + (q& -!((pos+i_compens_sublat)&1)) ; i < C.F() ; i += p+q) {
                if (LIKELY(S_ptr[i] == 255)) continue;
                bucket_update_t<1, primehint_t> prime;
                unsigned int x = ((size_t) (j-j0) << logI) + i;
                if (resieve_very_verbose) {
                    verbose_output_print(0, 1, "resieve_small_bucket_region: root %"
                            FBROOT_FORMAT ",%d divides at x = "
                            "%d = %u * %u + %d\n",
                            p, r, x, j, I, i);
                }
                prime.p = p;
                prime.x = x;
                ASSERT(prime.p >= si.conf.td_thresh);
                BP->push_update(prime);
            }
            pos += r;
            if (pos >= (spos_t) p) pos -= (spos_t) p;

            S_ptr += I;
            q ^= p;
        }
        pos += ssps.get_offset();
        if (pos >= (spos_t) p) pos -= p;

        ssdpos[index - ssd.resieve_start_offset] = pos;
    }

    for(auto const & ssp : ssd.ssp) {
        /* "not nice" cases are either projective or power of two. we
         * obviously won't resieve powers of two, so we're bound to deal
         * with only projective primes here.
         */
        ASSERT(ssp.is_pow2() || ssp.is_ordinary3() || ssp.is_proj());

        /* FIXME: I should not have to do this test */
        if (ssp.is_pow2() || ssp.is_ordinary3()) continue;

        /* TODO: it doesn't seem very smart to resieve projective primes
         */
        if (ssp.is_proj()) {
            const fbprime_t g = ssp.get_g();
            const fbprime_t q = ssp.get_q();
            const fbprime_t p = g * q;

            /* the code below definitely does not deal with projective
             * primes that do not sieve full lines.  */
            if (q > 1) continue;
            if (p < si.conf.td_thresh) continue;

            WHERE_AM_I_UPDATE(w, p, p);

            const uint64_t gI = (uint64_t)g << logI;

            /* Test every p-th line, starting at S[ssdpos] */
            spos_t pos = C.first_position_projective_prime(ssp);
            // This block is for the case where p divides at (1,0).
            if (UNLIKELY(has_origin && pos == (spos_t) gI)) {
                bucket_update_t<1, primehint_t> prime;
                prime.p = p;
                prime.x = 1 - i0;
                ASSERT(prime.p >= si.conf.td_thresh);
                BP->push_update(prime);
            }
            // Same as in sieving: we discard after checking for row 0.
            if (ssp.is_discarded_proj())
                continue;

            /* make sure ssdpos points at start of line or region when
             * we're sieving whole lines. */
            ASSERT(q > 1 || !(pos % (1 << MIN(logI, LOG_BUCKET_REGION))));

            if (resieve_very_verbose) {
                verbose_output_print(0, 1, "# resieving projective prime %" FBPRIME_FORMAT
                        ", i0 = %d\n", q, i0 + pos);
            }
            if (pos >> LOG_BUCKET_REGION)
                continue;
            unsigned int j = j0 + (pos >> logI);
            for (; j < j1; j += g) {
                unsigned char *S_ptr = S + pos;
                /* FIXME: sublat */
                if (!(j & 1)) {
                    /* Even j: test only odd ii-coordinates */
                    for (int ii = 1; ii < C.F(); ii += 2) {
                        if (S_ptr[ii] == 255) continue;
                        bucket_update_t<1, primehint_t> prime;
                        const unsigned int x = pos + ii;
                        if (resieve_very_verbose) {
                            verbose_output_print(0, 1, "# resieve_small_bucket_region even j: root %"
                                    FBROOT_FORMAT ",inf divides at x = %u\n",
                                    g, x);
                        }
                        prime.p = g;
                        prime.x = x;
                        BP->push_update(prime);
                    }
                } else {
                    /* Odd j: test all ii-coordinates */
                    for (int ii = 0; ii < C.F(); ii++) {
                        if (S_ptr[ii] == 255) continue;
                        bucket_update_t<1, primehint_t> prime;
                        const unsigned int x = pos + ii;
                        if (resieve_very_verbose) {
                            verbose_output_print(0, 1, "# resieve_small_bucket_region odd j: root %"
                                    FBROOT_FORMAT ",inf divides at x = %u\n",
                                    g, x);
                        }
                        prime.p = g;
                        prime.x = x;
                        BP->push_update(prime);
                    }
                }
                pos += gI;
            }
#if 0
            ssdpos[index] = pos - bucket_region;
            if (resieve_very_verbose) {
                verbose_output_print(0, 1, "# resieving: new pos = %" PRIu64
                        ", bucket_region = %d, "
                        "new ssdpos = %" PRIu64 "\n",
                        pos, bucket_region, ssdpos[index]);
            }
#endif
        }
    }
}
/* }}} */
