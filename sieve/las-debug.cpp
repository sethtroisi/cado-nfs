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
#include "portability.h"

#if defined(__GLIBC__) && (defined(TRACE_K) || defined(CHECK_UNDERFLOW))
#include <execinfo.h>   /* For backtrace. Since glibc 2.1 */
#endif

/* The trivial calls for when TRACE_K is *not* defined are inlines in
 * las-debug.h */


/* recall that TRACE_K requires TRACK_CODE_PATH ; so we may safely use
 * all where_am_I types here */

struct trace_Nx_t trace_Nx = { 0, UINT_MAX};
struct trace_ab_t trace_ab = { 0, 0 };
struct trace_ij_t trace_ij = { 0, UINT_MAX, };

/* two norms of the traced (a,b) pair */
mpz_t traced_norms[2];

/* This fills all the trace_* structures from the main one. The main
 * structure is the one for which a non-NULL pointer is passed.
 */
void trace_per_sq_init(sieve_info_srcptr si, const struct trace_Nx_t *Nx,
                       const struct trace_ab_t *ab,
                       const struct trace_ij_t *ij)
{
#ifndef TRACE_K
    if (Nx != NULL || ab != NULL || ij != NULL) {
        fprintf (stderr, "Error, relation tracing requested but this siever "
                 "was compiled without TRACE_K.\n");
        exit(EXIT_FAILURE);
    }
    return;
#endif
    /* At most one of the three coordinates must be specified */
    ASSERT_ALWAYS((Nx != NULL) + (ab != NULL) + (ij != NULL) <= 1);

    if (ab != NULL) {
      trace_ab = *ab;
      /* can possibly fall outside the q-lattice. We have to check for it */
      if (ABToIJ(&trace_ij.i, &trace_ij.j, trace_ab.a, trace_ab.b, si)) {
          IJToNx(&trace_Nx.N, &trace_Nx.x, trace_ij.i, trace_ij.j, si);
      } else {
          fprintf(stderr, "# Relation (%" PRId64 ",%" PRIu64 ") to be traced "
                  "is outside of the current q-lattice\n",
                  trace_ab.a, trace_ab.b);
          trace_ij.i=0;
          trace_ij.j=UINT_MAX;
          trace_Nx.N=0;
          trace_Nx.x=UINT_MAX;
      }
    } else if (ij != NULL) {
        trace_ij = *ij;
        IJToAB(&trace_ab.a, &trace_ab.b, trace_ij.i, trace_ij.j, si);
        IJToNx(&trace_Nx.N, &trace_Nx.x, trace_ij.i, trace_ij.j, si);
    } else if (Nx != NULL) {
        trace_Nx = *Nx;
        NxToIJ(&trace_ij.i, &trace_ij.j, trace_Nx.N, trace_Nx.x, si);
        IJToAB(&trace_ab.a, &trace_ab.b, trace_ij.i, trace_ij.j, si);
    }

    if (trace_ij.j < UINT_MAX && trace_ij.j >= si->J) {
        fprintf(stderr, "# Relation (%" PRId64 ",%" PRIu64 ") to be traced is "
                "outside of the current (i,j)-rectangle (j=%u)\n",
                trace_ab.a, trace_ab.b, trace_ij.j);
        trace_ij.i=0;
        trace_ij.j=UINT_MAX;
        trace_Nx.N=0;
        trace_Nx.x=UINT_MAX;
    }
    if (trace_ij.i || trace_ij.j < UINT_MAX) {
        fprintf(stderr, "# Tracing relation (a,b)=(%" PRId64 ",%" PRIu64 ") "
                "(i,j)=(%d,%u), (N,x)=(%u,%u)\n",
                trace_ab.a, trace_ab.b, trace_ij.i, trace_ij.j, trace_Nx.N,
                trace_Nx.x);
    }

    for(int side = 0 ; side < 2 ; side++) {
        mpz_init(traced_norms[side]);
        mpz_poly_homogeneous_eval_siui(traced_norms[side], 
                si->sides[side]->fij, trace_ij.i, trace_ij.j);
    }
}

void trace_per_sq_clear(sieve_info_srcptr si MAYBE_UNUSED)
{
    for(int side = 0 ; side < 2 ; side++)
        mpz_clear(traced_norms[side]);
}

#ifdef TRACE_K
int test_divisible(where_am_I_ptr w)
{
    /* we only check divisibility for the given (N,x) value */
    if (!trace_on_spot_Nx(w->N, w->x))
        return 1;

    /* Note that when we are reaching here through apply_one_bucket, we
     * do not know the prime number. */
    fbprime_t p = w->p;
    if (p==0) return 1;

    const unsigned int logI = w->si->conf->logI;
    const unsigned int I = 1U << logI;

    const unsigned long X = w->x + (w->N << LOG_BUCKET_REGION);
    long i = (long) (X & (I-1)) - (long) (I/2);
    unsigned long j = X >> logI;
    fbprime_t q;

    q = fb_is_power (p, NULL);
    if (q == 0)
        q = p;

    int rc = mpz_divisible_ui_p (traced_norms[w->side], (unsigned long) q);

    if (rc)
        mpz_divexact_ui (traced_norms[w->side], traced_norms[w->side], (unsigned long) q);
    else
        gmp_fprintf(stderr, "# FAILED test_divisible(p=%" FBPRIME_FORMAT
                ", N=%d, x=%lu, %.3s): i = %ld, j = %lu, norm = %Zd\n",
                w->p, w->N, w->x, sidenames[w->side], i, j, traced_norms[w->side]);

    return rc;
}

/* {{{ helper: sieve_increase */

/* Do this so that the _real_ caller is always 2 floors up */
void sieve_increase_logging_backend(unsigned char *S, const unsigned char logp, where_am_I_ptr w)
{
    if (!trace_on_spot_Nx(w->N, w->x))
        return;

    ASSERT_ALWAYS(test_divisible(w));

#ifdef __GLIBC__
    void * callers_addresses[3];
    char ** callers = NULL;
    backtrace(callers_addresses, 3);
    callers = backtrace_symbols(callers_addresses, 3);
    char * freeme = strdup(callers[2]);
    const char * caller = freeme;
    free(callers);
    char * opening = strchr(freeme, '(');
    if (opening) {
        char * closing = strchr(opening + 1, ')');
        if (closing) {
            *closing='\0';
            caller = opening + 1;
        }
    }
    if (!*caller) {
        caller="<no symbol (static?)>";
    }
#else
    const char * caller = "";
#endif
    if (w->p) 
        fprintf(stderr, "# Add log(%" FBPRIME_FORMAT ",%.3s) = %u to "
            "S[%u] = %hhu, from BA[%u] -> %hhu [%s]\n",
            w->p, sidenames[w->side], logp, w->x, *S, w->N, (unsigned char)(*S+logp), caller);
    else
        fprintf(stderr, "# Add log(hint=%lu,%.3s) = %u to "
            "S[%u] = %hhu, from BA[%u] -> %hhu [%s]\n",
            (unsigned long) w->h, sidenames[w->side], logp, w->x, *S, w->N, (unsigned char)(*S+logp), caller);
#ifdef __GLIBC__
    free(freeme);
#endif
}

void sieve_increase_logging(unsigned char *S, const unsigned char logp, where_am_I_ptr w)
{
    sieve_increase_logging_backend(S, logp, w);
}

/* Increase the sieve array entry *S by logp, with underflow checking
 * and tracing if desired. w is used only for trace test and output */

void sieve_increase(unsigned char *S, const unsigned char logp, where_am_I_ptr w)
{
    sieve_increase_logging_backend(S, logp, w);
#ifdef CHECK_UNDERFLOW
    sieve_increase_underflow_trap(S, logp, w);
#endif
    *S += logp;
}

#endif  /* TRACE_K */

/* This function is useful both with and without TRACE_K, as the flag
 * controlling it is CHECK_UNDERFLOW
 */
#ifdef CHECK_UNDERFLOW
void sieve_increase_underflow_trap(unsigned char *S, const unsigned char logp, where_am_I_ptr w)
{
    int i;
    unsigned int j;
    int64_t a;
    uint64_t b;
    static unsigned char maxdiff = ~0;

    NxToIJ(&i, &j, w->N, w->x, w->si);
    IJToAB(&a, &b, i, j, w->si);
    if ((unsigned int) logp + *S > maxdiff)
      {
        maxdiff = logp - *S;
        fprintf(stderr, "# Error, underflow at (N,x)=(%u, %u), "
                "(i,j)=(%d, %u), (a,b)=(%ld, %lu), S[x] = %hhu, log(%"
                FBPRIME_FORMAT ") = %hhu\n",
                w->N, w->x, i, j, a, b, *S, w->p, logp);
      }
    /* arrange so that the unconditional increase which comes next
     * has the effect of taking the result to maxdiff */
    *S = maxdiff - logp;
}
#endif

/* }}} */

