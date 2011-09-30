#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS

#include <limits.h>
#include <inttypes.h>

#include "las-debug.h"
#include "las-coordinates.h"

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

/* Test if entry x in bucket region n is divisible by p */
void test_divisible_x (const fbprime_t p, const unsigned long x, const int n,
		       sieve_info_srcptr si, int side);
