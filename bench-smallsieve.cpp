#include "cado.h"
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>
#include <time.h>
#include <vector>
#include <sstream>
#include <algorithm>

#define xxxLOG_BUCKET_REGION_IS_A_CONSTANT

#define HAVE_SSE2

#define WHERE_AM_I_UPDATE(a,b,c)  /**/
typedef void * where_am_I;
#define MAYBE_UNUSED __attribute__((unused))
#define ASSERT(x) assert(x)
#define croak__(x,y) do {						\
        fprintf(stderr,"%s in %s at %s:%d -- %s\n",			\
                (x),__func__,__FILE__,__LINE__,(y));			\
    } while (0)
#define ASSERT_ALWAYS(x)						\
    do {								\
        if (!(x)) {							\
            croak__("code BUG() : condition " #x " failed",		\
                    "Abort");						\
            abort();							\
        }								\
    } while (0)
#define EXPECT(x,val)	__builtin_expect(x,val)
#define LIKELY(x)	EXPECT(x,1)
#define UNLIKELY(x)	EXPECT(x,0)
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#define MAX(h,i) ((h) > (i) ? (h) : (i))
static inline void sieve_increase(unsigned char *S, const unsigned char logp, where_am_I& w MAYBE_UNUSED)
{
    *S += logp;
}


/* this does the small sieve loop, only for nice primes.
 * In theory, we could do this stuff for all primes, but we'll be happy
 * enough if we can bench code for the main branch that goes fast enough
 */

typedef uint32_t fbprime_t;

struct ssp_t {
    /* ordinary primes. Equation is (i-r*j) = 0 mod p */
    /* p may be a prime power, and even a power of two */
    /* Note that we need three fields for the projective primes anyway,
     * so we may have three here. Yet, admittedly the "offset" field is
     * rather useless. (see ssp_init_oa.)  We can recompute it on the fly
     * when needed.
     */
    fbprime_t p;	// q if projective
    fbprime_t r;        // in [ 0, p [. g if projective
    // unused:
    fbprime_t offset;   // in [ 0, p [. U if projective
    uint16_t flags=0;
    public:
    unsigned char logp;
    ssp_t(fbprime_t p, fbprime_t r, fbprime_t offset = 0)
        : p(p), r(r), offset(offset) {
            // we don't care if we put garbage.
            logp = p/3;
        }
    ssp_t() = default;
    /* this code does not deal with proj primes anyway, but the accessors
     * below are real */
    inline bool is_proj() const { return false; }

    fbprime_t get_p() const {ASSERT(!is_proj()); return p;}
    fbprime_t get_r() const {ASSERT(!is_proj()); return r;}
    fbprime_t get_offset() const {ASSERT(!is_proj()); return offset;}
};

#ifdef LOG_BUCKET_REGION_IS_A_CONSTANT
#define LOG_BUCKET_REGION 16
#else
int LOG_BUCKET_REGION = 16;
#endif

int interleaving = 1;   /* number of threads */

/* this is really a mock structure just for the fun of it. */
struct {
    struct {
        struct {
            uint32_t m = 0;
            uint32_t i0 = 0;
            uint32_t j0 = 0;
        } sublat;
        int logI_adjusted;
    } conf;
    struct {
        mpz_t q;
    } qbasis;
} si;

/* this is the exact same macro that is found in las-smallsieve.cpp and
 * (in part) in las-norms.cpp
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



/* {{{ preprocessor macro to generate assembly code for the small sieve */
/* {{{ Comments on the _OLDLOOP version (AF)
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
   }}} */
/* define sequences of operations rather systematically.
 *
 * ADD_NOINCR   *pi += logp.
 * ADD_INCR     *pi += logp; pi += incr
 * XADD_NOINCR  if (pi >= fence) goto exit; *pi += logp;
 * XADD_INCR    if (pi >= fence) goto exit; *pi += logp; pi += incr;
 */
/* {{{ *_ADD_NOINCR */
#define ONE_ADD_NOINCR(pi, logp)                                     	\
    "addb " logp ",(" pi ")\n"                /* pi[0] += logp */
#define TWO_ADDS_NOINCR(pi, incr, logp)                 		\
    ONE_ADD_NOINCR(pi, logp)               /* pi[0] += logp */       	\
    ONE_ADD_NOINCR(pi "," incr ",1", logp) /* pi[incr] += logp */
#define THREE_ADDS_NOINCR(pi, incr, logp)                 		\
    ONE_ADD_NOINCR(pi, logp)               /* pi[0] += logp */       	\
    ONE_ADD_NOINCR(pi "," incr ",1", logp) /* pi[incr] += logp */    	\
    ONE_ADD_NOINCR(pi "," incr ",2", logp) /* pi[2*incr] += logp */
/* }}} */
/* {{{ *_ADDS_INCR, including _PRECOMP variant */
#define ONE_ADD_INCR(pi, incr, logp)                             	\
    ONE_ADD_NOINCR(pi, logp)                                         	\
    "lea (" pi "," incr ",1)," pi "\n"        /* pi += incr */
#define TWO_ADDS_INCR(pi, incr, logp)                 			\
    ONE_ADD_NOINCR(pi, logp)               /* pi[0] += logp */       	\
    ONE_ADD_NOINCR(pi "," incr ",1", logp) /* pi[incr] += logp */    	\
    "lea (" pi "," incr ",2)," pi "\n"     /* pi += 2*incr */
/* The _PRECOMP variant uses a precomputed  bigincr = 3 * incr  */
#define THREE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)                \
     THREE_ADDS_NOINCR(pi, incr, logp)                                  \
    "lea (" pi "," bigincr ",1)," pi "\n"     /* pi += 3*incr */
/* {{{ expand to define new lengths */
#define FOUR_ADDS_INCR(pi, incr, logp)                                  \
    TWO_ADDS_INCR(pi, incr, logp)                                       \
    TWO_ADDS_INCR(pi, incr, logp)
#define EIGHT_ADDS_INCR(pi, incr, logp)                                 \
    FOUR_ADDS_INCR(pi, incr, logp)                                      \
    FOUR_ADDS_INCR(pi, incr, logp)
#define TWELVE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)               \
    THREE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)                    \
    THREE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)                    \
    THREE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)                    \
    THREE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)
#define SIXTEEN_ADDS_INCR(pi, incr, logp)                               \
    EIGHT_ADDS_INCR(pi, incr, logp)                                     \
    EIGHT_ADDS_INCR(pi, incr, logp)
/* }}} */
/* }}} */
/* {{{ caution: *_XADDS_INCR */
#define ONE_XADD_NOINCR(pi, fence, logp, exit)               	\
    "cmp " fence ", " pi "\n"                 /* if (pi >= S1) break; */\
    "jae " exit "\n"                                                    \
    ONE_ADD_NOINCR(pi, logp)
#define ONE_XADD_INCR(pi, incr, fence, logp, exit)       		\
    ONE_XADD_NOINCR(pi, fence, logp, exit)                   		\
    "lea (" pi "," incr ",1)," pi "\n"        /* pi += incr */
#define TWO_XADDS_INCR(pi, incr, fence, logp, exit)             	\
    ONE_XADD_INCR(pi, incr, fence, logp, exit)                  	\
    ONE_XADD_INCR(pi, incr, fence, logp, exit)
/* no need for THREE_XADDs */
/* {{{ expand to define new lengths */
#define FOUR_XADDS_INCR(pi, incr, fence, logp, exit)     		\
    TWO_XADDS_INCR(pi, incr, fence, logp, exit)          		\
    TWO_XADDS_INCR(pi, incr, fence, logp, exit)
#define EIGHT_XADDS_INCR(pi, incr, fence, logp, exit)    		\
    FOUR_XADDS_INCR(pi, incr, fence, logp, exit)         		\
    FOUR_XADDS_INCR(pi, incr, fence, logp, exit)
#define TWELVE_XADDS_INCR(pi, incr, fence, logp, exit) 			\
    EIGHT_XADDS_INCR(pi, incr, fence, logp, exit)        		\
    FOUR_XADDS_INCR(pi, incr, fence, logp, exit)
#define SIXTEEN_XADDS_INCR(pi, incr, fence, logp, exit)    		\
    EIGHT_XADDS_INCR(pi, incr, fence, logp, exit)         		\
    EIGHT_XADDS_INCR(pi, incr, fence, logp, exit)
/* }}} */
/* }}} */
/* {{{ We also have _NOLASTINCR variants. Those are used in the old
 * assembly loop. As a matter of fact, it is not a good idea to optimize
 * like this: we do need the end value of the pointer, and it's really
 * crucially important for the case there I>B. So using these functions
 * is part of the reason why the old code will not do in the longer term.
 */
#define TWO_XADDS_NOLASTINCR(pi, incr, fence, logp, exit)       	\
    ONE_XADD_INCR(pi, incr, fence, logp, exit)                  	\
    ONE_XADD_NOINCR(pi, fence, logp, exit)
#define FOUR_XADDS_NOLASTINCR(pi, incr, fence, logp, exit)   		\
    TWO_XADDS_INCR(pi, incr, fence, logp, exit)          		\
    TWO_XADDS_NOLASTINCR(pi, incr, fence, logp, exit)
#define TWELVE_XADDS_NOLASTINCR(pi, incr, fence, logp, exit) 		\
    EIGHT_XADDS_INCR(pi, incr, fence, logp, exit)        		\
    FOUR_XADDS_NOLASTINCR(pi, incr, fence, logp, exit)
#define SMALLSIEVE_ASSEMBLY_OLD(pi, incr, fence, logp) do {             \
    unsigned char *pi_end;                                              \
    size_t three_p_or_2p;                                               \
    __asm__ __volatile__ (                                              \
            "lea (%3,%3,2), %2\n"     /* three_p_or_2p = p_or_2p * 3 */ \
            "lea (%1,%2,4), %0\n"     /* pi_end = pi + p_or_2p*12 */    \
            "cmp %5, %0\n"            /* if (pi_end > S1) no loop */    \
            "jbe 0f\n"                                                  \
            "1:\n"                                                      \
            TWELVE_XADDS_NOLASTINCR("%1","%3","%5","%4","2f")           \
            "jmp 2f\n"                                                  \
            ".balign 8\n"                                               \
            "0:\n"                                                      \
            TWELVE_ADDS_INCR_PRECOMP("%1","%3","%2","%4")               \
            "lea (%1,%2,4), %0\n"    /* if (pi+p_or_2p*12 > S1) break */\
            "cmp %5, %0\n"                                              \
            "jbe 0b\n"                                                  \
            "jmp 1b\n"                                                  \
            ".balign 8\n"                                               \
            "2:\n"                                                      \
            : "=&r"(pi_end), "+&r"(pi), "=&r"(three_p_or_2p)            \
            : "r"(incr), "q"(logp), "r"(fence) : "cc");                 \
    pi = pi_end;                                                        \
} while (0)
/* }}} */
#if 0 /* {{{ unused GCC trick */
/* there is technically a way to drop the "volatile" qualifier, so that
 * gcc understands that we're touching memory, and which memory exactly.
 * Such a way can be found in the snippet below.
 * Alas, I don't know how it is possible to do so and at the same time
 * make sure that: 1 - the memory reference that is passed to the inline asm
 * is strictly of the form (%[register]), and 2 - we can actually get the
 * register back, or emit code that is similar to (%[register],%[displ])
 * from there.
 * And anyway, I don't really believe that it's important: this block is
 * really the expensive operation, what happens to the glue around it is
 * less important.
 */
#define SMALLSIEVE_ASSEMBLY_NEW(S0, pos, incr, fence, logp) do {        \
    /* see gcc 6.45.2.6 Clobbers */                                     \
    typedef unsigned char (*ss_area_t)[1<<16];                          \
    ss_area_t S0x = (ss_area_t) S0;                                     \
    asm (                                                               \
            "leaq (%[S0],%[pos]), %[pi]\n"                              \
            ONE_ADD_INCR("%[pi]","%[inc]", "%[logp]")            \
            ONE_XADD_INCR("%[pi]", "%[inc]", "%[ff]", "%[logp]", "0f") \
            ONE_XADD_INCR("%[pi]", "%[inc]", "%[ff]", "%[logp]", "0f") \
            ONE_XADD_INCR("%[pi]", "%[inc]", "%[ff]", "%[logp]", "0f") \
            "0:\n"                                                      \
            : [S0]"=m"(*S0x), [pi]"=&r"(pi)                             \
            : [pos]"r"(pos),                                            \
            [inc]"r"(incr), [logp]"q"(logp), [ff] "r" (fence));         \
} while (0)
#endif /* }}} */
/* {{{ asm volatile statements for all the variants */
/* {{{ unrolled statements, limited length */
#define SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp)          \
            "0:\n"                                                      \
            : [pi]"+r"(_pi)                                             \
            : [inc]"r"(_incr), [logp]"q"(_logp), [ff] "r" (_fence)
#define SMALLSIEVE_ASSEMBLY_NEW_1_TO_2(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            ONE_ADD_INCR("%[pi]","%[inc]","%[logp]")                    \
            ONE_XADD_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")      \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_2_TO_4(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            TWO_ADDS_INCR("%[pi]","%[inc]","%[logp]")                   \
            TWO_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_4_TO_8(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            FOUR_ADDS_INCR("%[pi]","%[inc]","%[logp]")                  \
            FOUR_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_8_TO_16(_pi, _incr, _fence, _logp)      \
    asm volatile (                                                      \
            EIGHT_ADDS_INCR("%[pi]","%[inc]","%[logp]")                 \
            EIGHT_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_16_TO_32(_pi, _incr, _fence, _logp)     \
    asm volatile (                                                      \
            SIXTEEN_ADDS_INCR("%[pi]","%[inc]","%[logp]")               \
            SIXTEEN_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
/* this is the point where the complete unrolling starts being excessive */
#define SMALLSIEVE_ASSEMBLY_NEW_32_TO_64(_pi, _incr, _fence, _logp)     \
    asm volatile (                                                      \
            SIXTEEN_ADDS_INCR("%[pi]","%[inc]","%[logp]")               \
            SIXTEEN_ADDS_INCR("%[pi]","%[inc]","%[logp]")               \
            SIXTEEN_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SIXTEEN_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
/* }}} */
/* {{{ loop8 */
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP8_COMMON(_pi, _incr, _fence, _logp, _PP) \
    asm volatile (                                                      \
            "leaq (%[pi], %[inc], 8), %%rdi\n"                          \
            "cmpq %[ff], %%rdi\n"                                       \
            "ja 1f\n"           /* if rdi > fence, no loop */           \
            ".balign 8\n"                                               \
            "2:\n"                                                      \
            _PP                                                         \
            EIGHT_ADDS_INCR("%[pi]","%[inc]","%[logp]")                 \
            "leaq (%%rdi, %[inc], 8), %%rdi\n"                          \
            "cmpq %[ff], %%rdi\n"                                       \
            "jbe 2b\n"                                                  \
            "1:\n"                                                      \
            EIGHT_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")   \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp)      \
            : "rdi")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP8(_pi, _incr, _fence, _logp)        \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP8_COMMON(_pi, _incr, _fence, _logp, \
                "")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP8P0(_pi, _incr, _fence, _logp)      \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP8_COMMON(_pi, _incr, _fence, _logp, \
                "prefetch  (%%rdi)\n")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP8P1(_pi, _incr, _fence, _logp)      \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP8_COMMON(_pi, _incr, _fence, _logp, \
                "prefetch  (%%rdi, %[inc], 8)\n")
/* }}} */
/* {{{ loop12 */
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP12_COMMON(_pi, _incr, _fence, _logp, _PP) \
    asm volatile (                                                      \
            "leaq (%[inc], %[inc], 2), %%r10\n"                         \
            "leaq (%[pi], %%r10, 4), %%rdi\n"                           \
            "cmpq %[ff], %%rdi\n"                                       \
            "ja 1f\n"           /* if rdi > fence, no loop */           \
            ".balign 8\n"                                               \
            "2:\n"                                                      \
            _PP                                                         \
            TWELVE_ADDS_INCR_PRECOMP("%[pi]","%[inc]","%%r10","%[logp]")\
            "leaq (%%rdi, %%r10, 4), %%rdi\n"                           \
            "cmpq %[ff], %%rdi\n"                                       \
            "jbe 2b\n"                                                  \
            "1:\n"                                                      \
            TWELVE_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")  \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp)      \
            : "rdi", "r10")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP12(_pi, _incr, _fence, _logp)       \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP12_COMMON(_pi, _incr, _fence, _logp,\
                "")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP12P0(_pi, _incr, _fence, _logp)     \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP12_COMMON(_pi, _incr, _fence, _logp,\
                "prefetch  (%%rdi)\n")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP12P1(_pi, _incr, _fence, _logp)    \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP12_COMMON(_pi, _incr, _fence, _logp,\
                "prefetch  (%%rdi, %%r10, 4)\n")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP12P2(_pi, _incr, _fence, _logp)    \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP12_COMMON(_pi, _incr, _fence, _logp,\
                "prefetch  (%%rdi, %%r10, 8)\n")
/* }}} */
/* {{{ loop16 */
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP16_COMMON(_pi, _incr, _fence, _logp, _PP) \
    asm volatile (                                                      \
            "leaq (%[pi], %[inc2], 8), %%rdi\n"                         \
            "cmpq %[ff], %%rdi\n"                                       \
            "ja 1f\n"           /* if rdi > fence, no loop */           \
            ".balign 8\n"                                               \
            "2:\n"                                                      \
            _PP                                                         \
            SIXTEEN_ADDS_INCR("%[pi]","%[inc]","%[logp]")               \
            "leaq (%%rdi, %[inc2], 8), %%rdi\n"                         \
            "cmpq %[ff], %%rdi\n"                                       \
            "jbe 2b\n"                                                  \
            "1:\n"                                                      \
            SIXTEEN_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp)      \
            , [inc2]"r"(2*(_incr)) : "rdi")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP16(_pi, _incr, _fence, _logp)        \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP16_COMMON(_pi, _incr, _fence, _logp, \
                "")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP16P0(_pi, _incr, _fence, _logp)      \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP16_COMMON(_pi, _incr, _fence, _logp, \
                "prefetch  (%%rdi)\n")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP16P1(_pi, _incr, _fence, _logp)      \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP16_COMMON(_pi, _incr, _fence, _logp, \
                "prefetch  (%%rdi, %[inc2], 8)\n")
/*}}}*/
/* {{{ the routines below are unused. The might be, in case we insist on
 * using the same code for both even and odd lines (but then it would
 * probably be a loop anyway) */
#define SMALLSIEVE_ASSEMBLY_NEW_0_TO_1(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            ONE_XADD_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")      \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_0_TO_2(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            TWO_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")     \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_1_TO_4(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            ONE_ADD_INCR("%[pi]","%[inc]","%[logp]")                    \
            TWO_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")     \
            ONE_XADD_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")      \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp)) 
/*}}}*/
/*}}}*/
/* }}} */

/*{{{ function objects for all the routines that we have */
#define BEGIN_FOBJ(name_)						\
    struct name_ {							\
        static constexpr const char * name = # name_;				\
        inline size_t operator()(					\
                unsigned char * S0,					\
                unsigned char * S1,					\
                size_t x0 MAYBE_UNUSED,					\
                size_t pos,						\
                size_t p_or_2p,						\
                unsigned char logp,					\
                where_am_I w MAYBE_UNUSED) const			\
        {								\
            unsigned char * pi = S0 + pos
#define END_FOBJ()							\
    return pi - S1;							\
}}
/* {{{ function objects for assembly code */
BEGIN_FOBJ(assembly_generic_loop8);
SMALLSIEVE_ASSEMBLY_NEW_LOOP8(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop12);
SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop16);
SMALLSIEVE_ASSEMBLY_NEW_LOOP16(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop8p0);
SMALLSIEVE_ASSEMBLY_NEW_LOOP8P0(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop12p0);
SMALLSIEVE_ASSEMBLY_NEW_LOOP12P0(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop16p0);
SMALLSIEVE_ASSEMBLY_NEW_LOOP16P0(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop8p1);
SMALLSIEVE_ASSEMBLY_NEW_LOOP8P1(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop12p1);
SMALLSIEVE_ASSEMBLY_NEW_LOOP12P1(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop16p1);
SMALLSIEVE_ASSEMBLY_NEW_LOOP16P1(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop12p2);
SMALLSIEVE_ASSEMBLY_NEW_LOOP12P2(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_oldloop);
SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly0);
SMALLSIEVE_ASSEMBLY_NEW_0_TO_1(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly1);
SMALLSIEVE_ASSEMBLY_NEW_1_TO_2(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly2);
SMALLSIEVE_ASSEMBLY_NEW_2_TO_4(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly3);
SMALLSIEVE_ASSEMBLY_NEW_4_TO_8(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly4);
SMALLSIEVE_ASSEMBLY_NEW_8_TO_16(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly5);
SMALLSIEVE_ASSEMBLY_NEW_16_TO_32(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly6);
SMALLSIEVE_ASSEMBLY_NEW_32_TO_64(pi, p_or_2p, S1, logp);
END_FOBJ();
/* }}} */
/* {{{ Function objects for C code, manually unrolled except for the
 * _loop versions. */
BEGIN_FOBJ(manual0);
if (pi < S1) { *pi += logp; pi += p_or_2p; }
END_FOBJ();
BEGIN_FOBJ(manual1);
do {
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; pi += p_or_2p;
} while (0);
END_FOBJ();
BEGIN_FOBJ(manual2);
do {
    *pi += logp; pi += p_or_2p;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; pi += p_or_2p;
} while (0);
END_FOBJ();
BEGIN_FOBJ(manual3);
do {
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; pi += p_or_2p;
} while (0);
END_FOBJ();
BEGIN_FOBJ(manual4);
do {
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; pi += p_or_2p;
} while (0);
END_FOBJ();
BEGIN_FOBJ(manual_oldloop);
#define T do {                                                          \
    WHERE_AM_I_UPDATE(w, x, x0 + pi - S0);                              \
    sieve_increase (pi, logp, w); pi += p_or_2p;                        \
} while(0)
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
#undef T
END_FOBJ();
BEGIN_FOBJ(manual_oldloop_nounroll);
for ( ; pi < S1 ; pi += p_or_2p) {
    WHERE_AM_I_UPDATE(w, x, x0 + pi - S0);
    sieve_increase (pi, logp, w);
}
END_FOBJ();
/*}}}*/
/*}}}*/

/* {{{ our selection of best routines indexed by ceil(log2(I/p)) */
template<int bit> struct best_evenline;
template<int bit> struct best_oddline ;
#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && !defined(TRACK_CODE_PATH)
template<int bit> struct best_evenline { typedef assembly_generic_loop16 type; };
template<int bit> struct best_oddline  { typedef assembly_generic_loop16 type; };
template<> struct best_evenline<1> { typedef manual0 type; };
template<> struct best_evenline<2> { typedef manual1 type; };
template<> struct best_evenline<3> { typedef manual2 type; };
template<> struct best_evenline<4> { typedef assembly3 type; };
template<> struct best_evenline<5> { typedef assembly4 type; };
template<> struct best_evenline<6> { typedef assembly5 type; };
template<> struct best_evenline<7> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<8> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<9> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<10> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<11> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<12> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<13> { typedef assembly_generic_loop16 type; };

template<> struct best_oddline<1> { typedef assembly1 type; };
template<> struct best_oddline<2> { typedef assembly2 type; };
template<> struct best_oddline<3> { typedef manual3 type; };
template<> struct best_oddline<4> { typedef assembly4 type; };
template<> struct best_oddline<5> { typedef assembly_generic_loop8 type; };
template<> struct best_oddline<6> { typedef assembly_generic_loop8 type; };
template<> struct best_oddline<7> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<8> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<9> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<10> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<11> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<12> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<13> { typedef assembly_generic_loop16 type; };
/* TODO: we could perhaps provide hints so that the generated code can
 * merge some rounds of the loop. Or maybe the compiler will do that ? */
#else
template<int bit> struct best_evenline { typedef manual_oldloop type; };
template<int bit> struct best_oddline  { typedef manual_oldloop type; };
#endif
/* }}} */


/******************************************************************/
/* now provides routines for testing */

/* we create lists of functions, because from all these tidbits we need
 * to generate the more complete sieving functions below.
 */
struct list_nil {};
template<typename T, typename U> struct list_car {};

struct all_generic_candidates {
    typedef list_car<assembly_generic_loop8,
            list_car<assembly_generic_loop12,
            list_car<assembly_generic_loop16,
            list_car<assembly_generic_loop8p0,
            list_car<assembly_generic_loop12p0,
            list_car<assembly_generic_loop16p0,
            list_car<assembly_generic_loop8p1,
            list_car<assembly_generic_loop12p1,
            list_car<assembly_generic_loop16p1,
            list_car<assembly_generic_loop12p2,
            list_car<assembly_generic_oldloop,
            list_car<manual_oldloop,
            list_car<manual_oldloop_nounroll,
            list_nil>>> >>> >>> >>> > type;
};

template<int bit> struct all_candidates_for_evenline {
    typedef all_generic_candidates::type type;
};
template<int bit> struct all_candidates_for_oddline {
    typedef all_generic_candidates::type type;
};
template<> struct all_candidates_for_evenline<1> {
    typedef list_car<assembly0,
            list_car<manual0,
            all_generic_candidates::type>> type;
};
template<> struct all_candidates_for_evenline<2> {
    typedef list_car<assembly1,
            list_car<manual1,
            all_generic_candidates::type>> type;
};
template<> struct all_candidates_for_evenline<3> {
    typedef list_car<assembly2,
            list_car<manual2,
            all_generic_candidates::type>> type;
};
template<> struct all_candidates_for_evenline<4> {
    typedef list_car<assembly3,
            list_car<manual3,
            all_generic_candidates::type>> type;
};
template<> struct all_candidates_for_evenline<5> {
    typedef list_car<assembly4,
            list_car<manual4,
            all_generic_candidates::type>> type;
};
template<> struct all_candidates_for_evenline<6> {
    typedef list_car<assembly5,
            all_generic_candidates::type> type;
};
template<> struct all_candidates_for_evenline<7> {
    typedef list_car<assembly6,
            all_generic_candidates::type> type;
};
template<> struct all_candidates_for_oddline<1> {
    typedef list_car<assembly1,
            list_car<manual1,
            all_generic_candidates::type>> type;
};
template<> struct all_candidates_for_oddline<2> {
    typedef list_car<assembly2,
            list_car<manual2,
            all_generic_candidates::type>> type;
};
template<> struct all_candidates_for_oddline<3> {
    typedef list_car<assembly3,
            list_car<manual3,
            all_generic_candidates::type>> type;
};
template<> struct all_candidates_for_oddline<4> {
    typedef list_car<assembly4,
            list_car<manual4,
            all_generic_candidates::type>> type;
};
template<> struct all_candidates_for_oddline<5> {
    typedef list_car<assembly5,
            all_generic_candidates::type> type;
};
template<> struct all_candidates_for_oddline<6> {
    typedef list_car<assembly6,
            all_generic_candidates::type> type;
};

/* For 2^(I-k) <= p < 2^(I-k+1), we have: 2^k >= 2^I/p > 2^(k-1), and
 *
 * so that the number of intervals of width p that fit within 2^I,
 * which is (2^(I-1))/p, is between:
 *
 * 2^k > (2^(I-1))/p >= 2^(k-1)
 *
 * The number of hits (endpoints of the intervals) is thus <= 2^k.
 *
 * For a lower bound, 2^(k-1) is attained because of the offset.
 *
 * Note that on even lines, p behaves as if it was doubled.
 */

#define NBITS_LESS 12

/*{{{ many functions, #if-protected with NBITS_LESS ...*/
#if 1
#if (NBITS_LESS == 1)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#if 1
    /* C or asm, that doesn't make a huge difference */
    if (pi < S1) *pi += logp;
    pi += p_or_2p;
#else
    SMALLSIEVE_ASSEMBLY_NEW_0_TO_1(pi, p_or_2p, S1, logp);
#endif
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#if 0
    /* more expensive than the inline asm, although it's a bit of a
     * mystery while that same code is a winner for even lines on 1-bit
     * larger primes... */
    do {
        *pi += logp;
        if ((pi += p_or_2p) >= S1) break;
        *pi += logp;
        pi += p_or_2p;
    } while (0);
#else
    SMALLSIEVE_ASSEMBLY_NEW_1_TO_2(pi, p_or_2p, S1, logp);
#endif
    return pi - S1;
}
#endif/*}}}*/
#if (NBITS_LESS == 2)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#if 1
    /* actually slightly better than the asm */
    do {
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; pi += p_or_2p;
    } while (0);
#else
    SMALLSIEVE_ASSEMBLY_NEW_1_TO_2(pi, p_or_2p, S1, logp);
#endif
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#if 0
    /* the two are more on less on par, with a minor advantage for the
     * assembly. 
     */
    do {
        *pi += logp; pi += p_or_2p;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; pi += p_or_2p;
    } while (0);
#else
    SMALLSIEVE_ASSEMBLY_NEW_2_TO_4(pi, p_or_2p, S1, logp);
#endif
    return pi - S1;
}
#endif/*}}}*/
#if (NBITS_LESS == 3)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#if 1
    /* actually slightly better than the asm */
    do {
        *pi += logp; pi += p_or_2p;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; pi += p_or_2p;
    } while (0);
#else
    SMALLSIEVE_ASSEMBLY_NEW_2_TO_4(pi, p_or_2p, S1, logp);
#endif
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#if 1
    /* very mildly better than the ASM, but I'm not too sure, to be
     * honest.
     */
    do {
        *pi += logp; pi += p_or_2p;
        *pi += logp; pi += p_or_2p;
        *pi += logp; pi += p_or_2p;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; pi += p_or_2p;
    } while (0);
#else
    SMALLSIEVE_ASSEMBLY_NEW_4_TO_8(pi, p_or_2p, S1, logp);
#endif
    return pi - S1;
}
#endif/*}}}*/
#if (NBITS_LESS == 4)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#if 0
    /* here, we seem to get considerably better performance with the asm
     * code.
     */
    do {
        *pi += logp; pi += p_or_2p;
        *pi += logp; pi += p_or_2p;
        *pi += logp; pi += p_or_2p;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; pi += p_or_2p;
    } while (0);
#else
    SMALLSIEVE_ASSEMBLY_NEW_4_TO_8(pi, p_or_2p, S1, logp);
#endif
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#if 0
    /* we get some improvement with the asm
     */
    do {
        *pi += logp; pi += p_or_2p;
        *pi += logp; pi += p_or_2p;
        *pi += logp; pi += p_or_2p;
        *pi += logp; pi += p_or_2p;
        *pi += logp; pi += p_or_2p;
        *pi += logp; pi += p_or_2p;
        *pi += logp; pi += p_or_2p;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; if ((pi += p_or_2p) >= S1) break;
        *pi += logp; pi += p_or_2p;
    } while (0);
#else
    SMALLSIEVE_ASSEMBLY_NEW_8_TO_16(pi, p_or_2p, S1, logp);
#endif
    return pi - S1;
}
#endif/*}}}*/
#if (NBITS_LESS == 5)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#if 0
    /* C code clearly bad here */
#else
    SMALLSIEVE_ASSEMBLY_NEW_8_TO_16(pi, p_or_2p, S1, logp);
#endif
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
#if 0
    /* C manually unrolled code is clearly catastrophic ! */
#else
    /* old code is actually a decent option, almost on par */
    SMALLSIEVE_ASSEMBLY_NEW_16_TO_32(pi, p_or_2p, S1, logp);
#endif
    return pi - S1;
}
#endif/*}}}*/
#if (NBITS_LESS == 6)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_16_TO_32(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    /* better to use a loop here */
    SMALLSIEVE_ASSEMBLY_NEW_LOOP8(pi, p_or_2p, S1, logp);
    return pi - S1;
}
#endif/*}}}*/
#if (NBITS_LESS == 7)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    /* better to use a loop here */
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
#endif/*}}}*/
/* at this point, choosing a LOOP12 or a LOOP16 is rather a matter of a
 * 0.1% difference.
 */
#if (NBITS_LESS == 8)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
#endif/*}}}*/
#if (NBITS_LESS == 9)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
#endif/*}}}*/
#if (NBITS_LESS == 10)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
#endif/*}}}*/
#if (NBITS_LESS == 11)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
#endif/*}}}*/
/* LOOP12, LOOP16: idem. Maybe a wee bit better for loop12 */
#if (NBITS_LESS == 12)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
#endif/*}}}*/
/* LOOP12, LOOP16: idem */
#if (NBITS_LESS == 13)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP16(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP16(pi, p_or_2p, S1, logp);
    return pi - S1;
}
#endif/*}}}*/
/* LOOP16 better than LOOP12 here */
#if (NBITS_LESS == 14)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP16(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP16(pi, p_or_2p, S1, logp);
    return pi - S1;
}
#endif/*}}}*/
#if (NBITS_LESS == 15)/*{{{*/
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP16(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP16(pi, p_or_2p, S1, logp);
    return pi - S1;
}
#endif/*}}}*/
#endif/*}}}*/

static inline size_t sieve_full_line(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S1, logp);
    return pi - S1;
}


void current_I18_branch(std::vector<int64_t> & positions, std::vector<ssp_t> primes, unsigned char * S, int logI, unsigned int N) /* {{{ */
{
    SMALLSIEVE_COMMON_DEFS();
    ASSERT_ALWAYS(positions.size() == primes.size());
    for(auto const & ssp : primes) {
        int64_t & p_pos(positions[&ssp - &primes.front()]);
        int pos = p_pos;

        const fbprime_t p = ssp.get_p();
        const fbprime_t r = ssp.get_r();
        const unsigned char logp = ssp.logp;
        unsigned char * S0 = S;
        where_am_I w MAYBE_UNUSED = 0;

        size_t overrun = 0; /* tame gcc */
        unsigned int i_compens_sublat = si.conf.sublat.i0 & 1;
        unsigned int j = j0;


        /* we sieve over the area [S0..S0+(i1-i0)], which may
         * actually be just a fragment of a line. After that, if
         * (i1-i0) is different from I, we'll break anyway. So
         * whether we add I or (i1-i0) to S0 does not matter much.
         */
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
        } else {
            /* skip stride */
            pos += ssp.get_offset();
            if (pos >= (int) p) pos -= p;
        }

        p_pos = pos;
    }
}/*}}}*/
void modified_I18_branch_C(std::vector<int64_t> & positions, std::vector<ssp_t> primes, unsigned char * S, int logI, unsigned int N) /* {{{ */
{
    SMALLSIEVE_COMMON_DEFS();
    ASSERT_ALWAYS(positions.size() == primes.size());
    for(auto const & ssp : primes) {
        int64_t & p_pos(positions[&ssp - &primes.front()]);
        int pos = p_pos;

        const fbprime_t p = ssp.get_p();
        const fbprime_t r = ssp.get_r();
        const unsigned char logp = ssp.logp;
        unsigned char * S0 = S;
        where_am_I w MAYBE_UNUSED = 0;

        size_t overrun = 0; /* tame gcc */

        /* we sieve over the area [S0..S0+(i1-i0)], which may
         * actually be just a fragment of a line. After that, if
         * (i1-i0) is different from I, we'll break anyway. So
         * whether we add I or (i1-i0) to S0 does not matter much.
         */
        fbprime_t spx = (p+p)^p;
        fbprime_t px = (j0&1)?p:(p+p);
        /* The increment is i1-i0, not I, to cover the case where 
         * B <= I; we then have i1-i0 = B
         * because there's a min. (and in that case, there's only a
         * single line).
         */
        for(unsigned int j = j0 ; j < j1 ; ++j) {
            /* for j even, we sieve only odd pi, so step = 2p. */
            WHERE_AM_I_UPDATE(w, j, j);
            overrun = sieve_full_line(S0, S0 + (i1 - i0), S0 - S,
                    pos + (((si.conf.sublat.i0^pos^1)&1&~j)?p:0), px, logp, w);
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p; 
            px ^= spx;
        }
        if (logI > LOG_BUCKET_REGION) {
            int N1 = N + interleaving;
            int Q = logI - LOG_BUCKET_REGION;
            int dj = (N1>>Q) - j0;
            int di = (N1&((1<<Q)-1)) - (N&((1<<Q)-1));
            int B_mod_p = overrun - pos;
            pos = (pos + B_mod_p * di + dj * r) % p;
            if (pos < 0) pos += p;
        } else {
            pos += ssp.get_offset();
            if (pos >= (int) p) pos -= p;
        }

        p_pos = pos;
    }
}/*}}}*/
void legacy_branch(std::vector<int64_t> & positions, std::vector<ssp_t> primes, unsigned char * S, int logI, unsigned int N) /*{{{*/
{
    SMALLSIEVE_COMMON_DEFS();
    const unsigned long bucket_region = (1UL << LOG_BUCKET_REGION);
    ASSERT_ALWAYS(positions.size() == primes.size());
    ASSERT_ALWAYS(logI <= LOG_BUCKET_REGION);
    where_am_I w MAYBE_UNUSED = 0;

    for(auto const & ssp : primes) {
        int64_t & p_pos(positions[&ssp - &primes.front()]);
        unsigned int pos = p_pos;

        const fbprime_t p = ssp.get_p();
        const fbprime_t r = ssp.get_r();
        const unsigned char logp = ssp.logp;

        unsigned long j;
        const int test_divisibility MAYBE_UNUSED = 0; /* very slow, but nice for debugging */
        const unsigned long nj = bucket_region >> si.conf.logI_adjusted; /* Nr. of lines 
                                                                            per bucket region */
        /* In order to check whether a j coordinate is even, we need to take
         * into account the bucket number, especially in case buckets are as
         * large as the sieve region. The row number corresponding to a given
         * i0 is i0/I, but we also need to add bucket_nr*bucket_size/I to
         * this, which is what this flag is for.
         * Sublat must also be taken into account.
         */
        int row0_is_oddj;
        if (si.conf.sublat.m == 0) {
            row0_is_oddj = (N << (LOG_BUCKET_REGION - si.conf.logI_adjusted)) & 1;
        } else {
            int row0 = (N << (LOG_BUCKET_REGION - si.conf.logI_adjusted));
            row0_is_oddj = (row0 * si.conf.sublat.m + si.conf.sublat.j0) & 1;
            // Odd/even property of j is the same as for j+2, even with
            // sublat, unless sublat.m is even, which is not handled right
            // now. Same for i.
            ASSERT_ALWAYS((si.conf.sublat.m & 1) == 1);
        }

        WHERE_AM_I_UPDATE(w, p, p);

        //XXX XXX XXX const unsigned char logp = ssd->logp[k];
        unsigned char *S_ptr = S;
        size_t p_or_2p = p;
        ASSERT(pos < p);
        unsigned int i_compens_sublat = 0;
        if (si.conf.sublat.m != 0) {
            i_compens_sublat = si.conf.sublat.i0 & 1;
        }
        j = 0;
        if (row0_is_oddj) goto j_odd;
        // it is possible to make this code work with sieve_full_line by
        // uncommenting all lines that match S0. Performance is
        // apparently identical, but it's not fair to say that this *is*
        // the reference code, as far as performance is concerned.
        // unsigned char * S0 MAYBE_UNUSED;
        do {
            unsigned char *pi;
            WHERE_AM_I_UPDATE(w, j, j);
            // S0 = S_ptr;
            pi = S_ptr + pos;
            S_ptr += I;
            /* for j even, we sieve only odd pi, so step = 2p. */
            p_or_2p += p_or_2p;
            if (!((i_compens_sublat + pos) & 1)) {
                pi += p;
            }
            SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S_ptr, logp);
            // sieve_full_line(S0, S_ptr, S0 - S, pi - S0, p_or_2p, logp, w);
            pos += r; if (pos >= p) pos -= p; 
            /* Next line */
            if (++j >= nj) break;
            p_or_2p >>= 1;
j_odd:
            WHERE_AM_I_UPDATE(w, j, j);
            // S0 = S_ptr;
            pi = S_ptr + pos;
            S_ptr += I;
            SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S_ptr, logp);
            // sieve_full_line(S0, S_ptr, S0 - S, pi - S0, p_or_2p, logp, w);
            pos += r; if (pos >= p) pos -= p;
        } while (++j < nj);
        //XXX XXX XXX XXX adding this here, as it is relevant to the comparison
        //with branch I18.
        pos += ssp.get_offset();
        if (pos >= (unsigned int) p)
            pos -= p;
        //XXX XXX XXX XXX

        p_pos = pos;
    }
}
/*}}}*/
void devel_branch(std::vector<int64_t> & positions, std::vector<ssp_t> primes, unsigned char * S, int logI, unsigned int N) /*{{{*/
{
    SMALLSIEVE_COMMON_DEFS();
    ASSERT_ALWAYS(positions.size() == primes.size());
    for(auto const & ssp : primes) {
        int64_t & p_pos(positions[&ssp - &primes.front()]);
        int pos = p_pos;

        const fbprime_t p = ssp.get_p();
        const fbprime_t r = ssp.get_r();
        const unsigned char logp = ssp.logp;
        unsigned char * S0 = S;
        where_am_I w MAYBE_UNUSED = 0;

        unsigned int j = j0;


        /* we sieve over the area [S0..S0+(i1-i0)], which may
         * actually be just a fragment of a line. After that, if
         * (i1-i0) is different from I, we'll break anyway. So
         * whether we add I or (i1-i0) to S0 does not matter much.
         */
        if (j0 & 1)
            goto j_odd_devel;

        for( ; j < j1 ; ) {
            /* for j even, we sieve only odd pi, so step = 2p. */
            {
            int xpos = ((si.conf.sublat.i0 + pos) & 1) ? pos : (pos+p);
            sieve_full_line_new_half(S0, S0 + (I), S0 - S,
                    xpos, p+p, logp, w);
            }
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p; 
            if (++j >= j1) break;

j_odd_devel:
            /* now j odd again */
            WHERE_AM_I_UPDATE(w, j, j - j0);
            sieve_full_line_new(S0, S0 + (I), S0 - S,
                    pos, p, logp, w);
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p;
            ++j;
        }
            /* skip stride */
            pos += ssp.get_offset();
            if (pos >= (int) p) pos -= p;

        p_pos = pos;
    }
}
/*}}}*/
template<typename even_code, typename odd_code> void devel_branch_meta(std::vector<int64_t> & positions, std::vector<ssp_t> primes, unsigned char * S, int logI, unsigned int N) /*{{{*/
{
    SMALLSIEVE_COMMON_DEFS();
    ASSERT_ALWAYS(positions.size() == primes.size());
    for(auto const & ssp : primes) {
        int64_t & p_pos(positions[&ssp - &primes.front()]);
        int pos = p_pos;

        const fbprime_t p = ssp.get_p();
        const fbprime_t r = ssp.get_r();
        const unsigned char logp = ssp.logp;
        unsigned char * S0 = S;
        where_am_I w MAYBE_UNUSED = 0;

        unsigned int j = j0;


        /* we sieve over the area [S0..S0+(i1-i0)], which may
         * actually be just a fragment of a line. After that, if
         * (i1-i0) is different from I, we'll break anyway. So
         * whether we add I or (i1-i0) to S0 does not matter much.
         */
        if (j0 & 1)
            goto j_odd_devel;

        for( ; j < j1 ; ) {
            /* for j even, we sieve only odd pi, so step = 2p. */
            {
            int xpos = ((si.conf.sublat.i0 + pos) & 1) ? pos : (pos+p);
            even_code()(S0, S0 + (I), S0 - S, xpos, p+p, logp, w);
            }
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p; 
            if (++j >= j1) break;

j_odd_devel:
            /* now j odd again */
            WHERE_AM_I_UPDATE(w, j, j - j0);
            odd_code()(S0, S0 + (I), S0 - S, pos, p, logp, w);
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p;
            ++j;
        }
            /* skip stride */
            pos += ssp.get_offset();
            if (pos >= (int) p) pos -= p;

        p_pos = pos;
    }
}/*}}}*/

typedef void (*ss_func)(std::vector<int64_t> & positions, std::vector<ssp_t> primes, unsigned char * S, int logI, unsigned int N);

struct cand_func {
    bool sel;
    ss_func f;
    const char * name;
    cand_func(bool sel, ss_func f, const char * name) : sel(sel), f(f), name(name) {}
};

typedef std::vector<cand_func> candidate_list;

/* first we create a list of functions that have the currently recorded
 * best code for even lines, and then we iterate on odd lines
 */

template<int bit> struct factory_for_bit_round1 {
    typedef typename best_evenline<bit>::type Be;
    typedef typename best_oddline<bit>::type Bo;
    template<typename T> struct iterator {
        void operator()(candidate_list &) {}
    };
    template<typename T, typename U> struct iterator<list_car<T, U>> {
        void operator()(candidate_list & res) const {
            bool sel = std::is_same<Bo, T>::value;
            res.push_back({sel, &devel_branch_meta<Be, T>, T::name});
            iterator<U>()(res);
        }
    };
    void operator()(candidate_list & res) const {
        iterator<typename all_candidates_for_oddline<bit>::type>()(res);
    }
};
template<int bit> struct factory_for_bit_round2 {
    typedef typename best_evenline<bit>::type Be;
    typedef typename best_oddline<bit>::type Bo;
    template<typename T> struct iterator {
        void operator()(candidate_list &) {}
    };
    template<typename T, typename U> struct iterator<list_car<T, U>> {
        void operator()(candidate_list & res) const {
            bool sel = std::is_same<Be, T>::value;
            res.push_back({sel, &devel_branch_meta<T, Bo>, T::name});
            iterator<U>()(res);
        }
    };
    void operator()(candidate_list & res) const {
        iterator<typename all_candidates_for_evenline<bit>::type>()(res);
    }
};

struct bench_base {
    std::vector<ssp_t> allprimes;
    std::vector<int64_t> positions;
    unsigned char * S;
    size_t B;
    int logI;
    int logA;

    bench_base(int logB, int logI, int logA) : logI(logI), logA(logA) {
        B = 1UL << logB;
        posix_memalign((void**)&S, B, B);
    }
    bench_base(bench_base const &) = delete;
    ~bench_base() { free(S); }
    void test(candidate_list const & cand, const char * pfx="", std::string const& goal="current selection") {

        if (cand.empty()) return;

        std::vector<unsigned char *> refS;
        size_t I = 1UL << logI;
        int Nmax = 1 << (logA - LOG_BUCKET_REGION);

        std::vector<double> timings;
        int sel_index = -1;

        for(auto const& bf : cand) {
            ss_func f = bf.f;
            if (bf.sel) sel_index = &bf-&cand.front();
            positions.clear();
            for(auto const & ssp : allprimes)
                positions.push_back((I/2)%ssp.get_p());
            memset(S, 0, B);
            clock_t tt = clock();
            for(int N = 0 ; N < Nmax ; N++) {
                (*f)(positions, allprimes, S, logI, N);
            }
            timings.push_back((double) (clock()-tt) / CLOCKS_PER_SEC);
            printf(".");
            fflush(stdout);
            /* FIXME: we're only checking the very last region, here ! */
            if (!refS.empty()) {
                if (memcmp(S, refS.front(), B) != 0) {
                    fprintf(stderr, "inconsistency between %s and %s\n",
                            bf.name, cand.front().name);
                    for(size_t i = 0, n = 0 ; i < B && n < 16 ; i++) {
                        if (S[i] != refS.front()[i]) {
                            fprintf(stderr, "%04x: %02x != %02x\n",
                                    (unsigned int) i,
                                    (unsigned int) S[i],
                                    (unsigned int) refS.front()[i]);
                            n++;
                        }
                    }
                }
            }
            unsigned char * Scopy;
            posix_memalign((void**)&Scopy, B, B);
            memcpy(Scopy, S, B);
            refS.push_back(Scopy);
        }
        printf("\n");

#define BBLACK(X) "\e[01;30m" X "\e[00;30m"
#define BRED(X)   "\e[01;31m" X "\e[00;30m"
#define BGREEN(X) "\e[01;32m" X "\e[00;30m"
#define BYELLOW(X)"\e[01;33m" X "\e[00;30m"
#define BBLUE(X)  "\e[01;34m" X "\e[00;30m"
#define BVIOLET(X)"\e[01;35m" X "\e[00;30m"
#define BNAVY(X)  "\e[01;36m" X "\e[00;30m"
#define BLACK(X) "\e[00;30m" X "\e[00;30m"
#define RED(X)   "\e[00;31m" X "\e[00;30m"
#define GREEN(X) "\e[00;32m" X "\e[00;30m"
#define YELLOW(X)"\e[00;33m" X "\e[00;30m"
#define BLUE(X)  "\e[00;34m" X "\e[00;30m"
#define VIOLET(X)"\e[00;35m" X "\e[00;30m"
#define NAVY(X)  "\e[00;36m" X "\e[00;30m"

        size_t best_index = std::min_element(timings.begin(), timings.end()) - timings.begin();
        double best_time = timings[best_index];
        if (sel_index >= 0) {
            bool sel_best = (size_t) sel_index == best_index;
            bool sel_tied = timings[sel_index] <= 1.05 * best_time;;
            for(auto const& bf : cand) {
                size_t index = &bf-&cand.front();
                double tt = timings[index];
                bool is_best = index == best_index;
                bool near_best = tt <= 1.05 * best_time;
                printf("%s%-24s:\t%.3f\t", pfx, bf.name, tt);
                if (bf.sel) {
                    if (sel_best) {
                        printf(BGREEN("%s") "\n", goal.c_str());
                    } else if (sel_tied) {
                        printf(GREEN("%s, tied") "\n", goal.c_str());
                    } else {
                        printf(BRED("%s, wrong ?") "\n", goal.c_str());
                    }
                } else if (is_best && !sel_best) {
                    if (sel_tied) {
                         printf(GREEN("current best, tied") "\n");
                    } else {
                        printf(BRED("current best") "\n");
                    }
                } else if (near_best) {
                    printf(GREEN("tied")"\n");
                } else {
                    printf("\n");
                }
            }
        } else {
            for(auto const& bf : cand) {
                size_t index = &bf-&cand.front();
                double tt = timings[index];
                bool is_best = index == best_index;
                bool near_best = tt <= 1.05 * best_time;
                printf("%s%-24s:\t%.3f\t", pfx, bf.name, tt);
                if (is_best) {
                    printf(GREEN("current best") "\n");
                } else if (near_best) {
                    printf(GREEN("tied")"\n");
                } else {
                    printf("\n");
                }
            }
        }

        for(unsigned char * Sp : refS)
            free(Sp);
    }
};

void store_primes(std::vector<ssp_t>& allprimes, int bmin, int bmax, gmp_randstate_t rstate)
{
    mpz_t pz;
    mpz_init(pz);
    mpz_set_ui(pz, 1u << bmin);
    for(;;) {
        mpz_nextprime(pz, pz);
        if (mpz_cmp_ui(pz, 1u << bmax) >= 0) break;
        unsigned long p = mpz_get_ui(pz);
        unsigned long r = gmp_urandomm_ui(rstate, p);
        allprimes.emplace_back(p, r, (r * (interleaving-1)) % p);
    }
    printf("created a list of %zu primes\n", allprimes.size());
    mpz_clear(pz);
}

candidate_list funcs {
    // modified_I18_branch_C,
        { false, current_I18_branch, "I18" },
        { false, legacy_branch, "legacy" },
        { false, devel_branch, "devel" },
};

int main(int argc, char * argv[])
{
    int logI = 16;
    int logA = 0;
    int bmin = 0;
    int bmax = 0;
    argv++,argc--;
    for( ; argc ; argv++, argc--) {
#ifndef LOG_BUCKET_REGION_IS_A_CONSTANT
        if (strcmp(*argv, "-B") == 0) {
            LOG_BUCKET_REGION = atoi(argv[1]);
            argv++,argc--;
            continue;
        }
#endif
        if (strcmp(*argv, "-I") == 0) {
            logI = atoi(argv[1]);
            argv++,argc--;
            continue;
        }
        if (strcmp(*argv, "-A") == 0) {
            logA = atoi(argv[1]);
            argv++,argc--;
            continue;
        }
        if (strcmp(*argv, "-bmax") == 0) {
            bmax = atoi(argv[1]);
            argv++,argc--;
            continue;
        }
        if (strcmp(*argv, "-bmin") == 0) {
            bmin = atoi(argv[1]);
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "unhandled arg: %s\n", *argv);
        exit(EXIT_FAILURE);
    }
    si.conf.logI_adjusted = logI;
    if (!bmax) bmax = logI;
    if (!logA) logA = 2*logI-1;

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    int nthreads=1;

    interleaving = nthreads;

    mpz_init_set_ui(si.qbasis.q, 4294967291);


    // std::vector<std::pair<size_t, size_t>> bounds_perbit(logI + 1, {0,0});
    for(int b = bmin ; b < bmax ; b++) {
        if (b == 0) continue;
#if 0
        /* primes within [1<<b, 1<<(b+1)] are (b+1)-bit primes. With respect to the
         * semantics that we use in the template code, those are counted as
         * (logI-bit) off the bound. So they will use the code in
         * all_candidates_for_evenline<logI-bit> and
         * all_candidates_for_oddline<logI-bit>.
         */
        size_t i0 = 0;
        for( ; i0 < allprimes.size() && !(allprimes[i0].get_p() >> b) ; i0++);
        size_t i1 = i0;
        for( ; i1 < allprimes.size() && (allprimes[i1].get_p() >> b) ; i1++);
        bounds_perbit[b] = { i0, i1 };
#endif

        printf("===== now doing specific tests for %d-bit primes using candidates<%d> =====\n", b+1, logI-b);
        bench_base bbase(LOG_BUCKET_REGION, logI, logA);
        // std::vector<ssp_t>& allprimes(bbase.allprimes);
        // std::vector<int64_t>& positions(bbase.positions);
        // unsigned char * & S (bbase.S);
        store_primes(bbase.allprimes, b, b + 1, rstate);


        // std::vector<ssp_t> lp(allprimes.begin() + i0, allprimes.begin() + i1);
        // std::vector<int64_t>& lx(positions.begin() + i0, positions.begin() + i1);
        candidate_list cand1, cand2;

        switch(logI-b) {
#define CASE(xxx)							\
            case xxx:							\
                factory_for_bit_round1<xxx>()(cand1);	\
                factory_for_bit_round2<xxx>()(cand2);	\
                break
            CASE(1); CASE(2); CASE(3); CASE(4);
            CASE(5); CASE(6); CASE(7); CASE(8);
            CASE(9); CASE(10); CASE(11); CASE(12);
            CASE(13); CASE(14); CASE(15);
            default:
            fprintf(stderr, "sorry, not handled\n");
        }

        std::vector<unsigned char *> refS;

        if (!cand1.empty()) {
            std::ostringstream goal_name;
            goal_name << "best_oddline<" << logI-b << ">";
            printf("  testing odd lines ");
            bbase.test(cand1, "    ", goal_name.str());
        }
        if (!cand2.empty()) {
            std::ostringstream goal_name;
            goal_name << "best_evenline<" << logI-b << ">";
            printf("  testing even lines ");
            bbase.test(cand2, "    ", goal_name.str());
        }
    }

    printf("===== now doing tests for complete functions =====\n");
    {
        bench_base bbase(LOG_BUCKET_REGION, logI, logA);
        // std::vector<ssp_t>& allprimes(bbase.allprimes);
        // std::vector<int64_t>& positions(bbase.positions);
        // unsigned char * & S (bbase.S);
        store_primes(bbase.allprimes, bmin, bmax, rstate);

        bbase.test(funcs);
    }

    mpz_clear(si.qbasis.q);
    gmp_randclear(rstate);
}

