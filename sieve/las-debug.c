#include "cado.h"

#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS

#include <limits.h>
#include <inttypes.h>
#include <string.h>

#include "las-config.h"
#include "las-types.h"
#include "las-debug.h"
#include "las-coordinates.h"
#include "mpz_poly.h"
#include "portability.h"

#if defined(__GLIBC__) && (defined(TRACE_K) || defined(CHECK_UNDERFLOW))
#include <execinfo.h>   /* For backtrace. Since glibc 2.1 */
#endif


#ifdef  TRACE_Nx
struct trace_Nx_t trace_Nx = TRACE_Nx;
#else
struct trace_Nx_t trace_Nx = { 0, UINT_MAX};
#endif

#ifdef TRACE_AB
struct trace_ab_t trace_ab = TRACE_AB;
#else
struct trace_ab_t trace_ab = { 0, 0 };
#endif

#ifdef  TRACE_IJ
struct trace_ij_t trace_ij = TRACE_IJ;
#else
struct trace_ij_t trace_ij = { 0, UINT_MAX, };
#endif

void trace_update_conditions(sieve_info_srcptr si MAYBE_UNUSED)
{
#ifdef TRACE_K
    if (trace_ab.a || trace_ab.b) {
        if (!ABToIJ(&trace_ij.i, &trace_ij.j, trace_ab.a, trace_ab.b, si)) {
            fprintf(stderr, "# Relation (%"PRId64",%"PRIu64") to be traced is outside of the current q-lattice\n", trace_ab.a, trace_ab.b);
            trace_ij.i=0;
            trace_ij.j=UINT_MAX;
            trace_Nx.N=0;
            trace_Nx.x=UINT_MAX;
        } else {
            IJToNx(&trace_Nx.N, &trace_Nx.x, trace_ij.i, trace_ij.j, si);
        }
    } else if (trace_ij.i || trace_ij.j < UINT_MAX) {
        IJToAB(&trace_ab.a, &trace_ab.b, trace_ij.i, trace_ij.j, si);
        IJToNx(&trace_Nx.N, &trace_Nx.x, trace_ij.i, trace_ij.j, si);
    } else if (trace_Nx.N || trace_Nx.x < UINT_MAX) {
        NxToIJ(&trace_ij.i, &trace_ij.j, trace_Nx.N, trace_Nx.x, si);
        IJToAB(&trace_ab.a, &trace_ab.b, trace_ij.i, trace_ij.j, si);
    }
    if (trace_ij.j < UINT_MAX && trace_ij.j >= si->J) {
        fprintf(stderr, "# Relation (%"PRId64",%"PRIu64") to be traced is outside of the current (i,j)-rectangle (j=%u)\n", trace_ab.a, trace_ab.b, trace_ij.j);
        trace_ij.i=0;
        trace_ij.j=UINT_MAX;
        trace_Nx.N=0;
        trace_Nx.x=UINT_MAX;
    }
    if (trace_ij.i || trace_ij.j < UINT_MAX) {
        fprintf(stderr, "# Tracing relation (a,b)=(%"PRId64",%"PRIu64") (i,j)=(%d,%u), (N,x)=(%u,%u)\n",
                trace_ab.a, trace_ab.b, trace_ij.i, trace_ij.j, trace_Nx.N, trace_Nx.x);
    }
#else
    return;
#endif
}

#if defined(TRACK_CODE_PATH) && defined(WANT_ASSERT_EXPENSIVE)
int test_divisible(where_am_I_ptr w)
{
    fbprime_t p = w->p;
    if (p==0) return 1;
    const unsigned int logI = w->si->conf->logI;
    const unsigned int I = 1U << logI;

    const unsigned long X = w->x + (w->N << LOG_BUCKET_REGION);
    long i = (long) (X & (I-1)) - (long) (I/2);
    unsigned long j = X >> logI;
    mpz_t v;

    mpz_init(v);
    mp_poly_homogeneous_eval_siui(v,
            w->si->sides[w->side]->fij,
            w->si->cpoly->pols[w->side]->degree, i, j);

    int rc = mpz_divisible_ui_p(v, (unsigned long) p);

    if (!rc) {
        gmp_fprintf(stderr, "# FAILED test_divisible(" FBPRIME_FORMAT
                ", %d, %lu, %.3s): i = %ld, j = %lu, norm = %Zd\n",
                w->p, w->N, w->x, sidenames[w->side], i, j, v);
    }
    mpz_clear(v);

    ASSERT(rc);

    return rc;
}
#endif

/* {{{ helper: sieve_decrease */
/* Decrease the sieve array entry *S by logp, with underflow checking 
   and tracing if desired. Variables x, bucket_nr, p, and si
   are used only for trace test and output */

#ifdef CHECK_UNDERFLOW
void sieve_decrease_underflow_trap(unsigned char *S, const unsigned char logp, where_am_I_ptr w)
{
    int i;
    unsigned int j;
    int64_t a;
    uint64_t b;
    static unsigned char maxdiff = 0;

    NxToIJ(&i, &j, w->N, w->x, w->si);
    IJToAB(&a, &b, i, j, w->si);
    if (logp - *S > maxdiff)
      {
        maxdiff = logp - *S;
        fprintf(stderr, "# Error, underflow at (N,x)=(%u, %u), "
                "(i,j)=(%d, %u), (a,b)=(%ld, %lu), S[x] = %hhu, log("
                FBPRIME_FORMAT ") = %hhu\n",
                w->N, w->x, i, j, a, b, *S, w->p, logp);
      }
    /* arrange so that the unconditional decrease which comes next
     * has the effect of taking the result to zero */
    *S = logp;
}
#endif

#ifdef TRACE_K
/* Do this so that the _real_ caller is always 2 floors up */
void sieve_decrease_logging_backend(unsigned char *S, const unsigned char logp, where_am_I_ptr w)
{
    ASSERT(test_divisible(w));
    if (!trace_on_spot_Nx(w->N, w->x))
        return;

#ifdef __GLIBC__
    void * callers_addresses[3];
    char ** callers = NULL;
    backtrace(callers_addresses, 3);
    callers = backtrace_symbols(callers_addresses, 3);
    char * freeme = strdup(callers[2]);
    char * caller = freeme;
    free(callers);
    char * opening = strchr(caller, '(');
    if (opening) {
        char * closing = strchr(opening + 1, ')');
        if (closing) {
            *closing='\0';
            caller = opening + 1;
        }
    }
#else
    const char * caller = "";
#endif
    fprintf(stderr, "# Subtract log(" FBPRIME_FORMAT ",%.3s) = %u from "
            "S[%u] = %hhu, from BA[%u] -> %hhu [%s]\n",
            w->p, sidenames[w->side], logp, w->x, *S, w->N, *S-logp, caller);
#ifdef __GLIBC__
    free(freeme);
#endif
}
void sieve_decrease_logging(unsigned char *S, const unsigned char logp, where_am_I_ptr w)
{
    sieve_decrease_logging_backend(S, logp, w);
}
#endif

#ifdef TRACE_K
void sieve_decrease(unsigned char *S, const unsigned char logp, where_am_I_ptr w)
{
    sieve_decrease_logging_backend(S, logp, w);
#ifdef CHECK_UNDERFLOW
    sieve_decrease_underflow_trap(S, logp, w);
#endif
    *S = (*S < logp) ? 0 : (*S - logp);
}
#endif
/* }}} */

