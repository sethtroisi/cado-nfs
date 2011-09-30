#ifndef LAS_DEBUG_H_
#define LAS_DEBUG_H_

#include <stdint.h>

#include "las-config.h"
#include "las-types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Define CHECK_UNDERFLOW to check for underflow when subtracting
   the rounded log(p) from sieve array locations */
//#define CHECK_UNDERFLOW


/* define TRACE_K, and exactly one of the TRACE_* values to something
 * non-zero, in order to get tracing information on a particular
 * relation.  In particular this traces the sieve array entry
 * corresponding to the relation. Upon startup, the three values below
 * are reconciled.
 *
 * (see Coordinate systems further down in this file)
 */
#define xxxTRACE_K
#define TRACE_AB { 2587866,6337 }
// #define TRACE_IJ
// #define TRACE_Nx


struct trace_Nx_t { unsigned int N; unsigned int x; };
struct trace_ab_t { int64_t a; uint64_t b; };
struct trace_ij_t { int i; unsigned int j; };

struct trace_Nx_t trace_Nx;
struct trace_ab_t trace_ab;
struct trace_ij_t trace_ij;

static inline int trace_on_spot_Nx(unsigned int N, unsigned int x) {
    if (trace_Nx.x == UINT_MAX) return 0;
    return N == trace_Nx.N && x == trace_Nx.x;
}

static inline int trace_on_range_Nx(unsigned int N, unsigned int x0, unsigned int x1) {
    if (trace_Nx.x == UINT_MAX) return 0;
    return N == trace_Nx.N && x0 <= trace_Nx.x && trace_Nx.x < x1;
}

static inline int trace_on_spot_x(unsigned int x) {
    return x == (trace_Nx.N << LOG_BUCKET_REGION) + trace_Nx.x;
}

static inline int trace_on_spot_ab(int a, unsigned int b) {
    return a == trace_ab.a && b == trace_ab.b;
}

static inline int trace_on_spot_ij(int i, unsigned int j) {
    return i == trace_ij.i && j == trace_ij.j;
}

extern void test_divisible_x (const fbprime_t p, const unsigned long x, const int n,
		       sieve_info_srcptr si, int side);
extern void trace_update_conditions(sieve_info_srcptr si);

#ifdef __cplusplus
}
#endif

#endif	/* LAS_DEBUG_H_ */
