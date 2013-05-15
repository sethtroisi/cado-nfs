#ifndef LAS_DEBUG_H_
#define LAS_DEBUG_H_

#include <stdint.h>

#include "las-config.h"
#include "las-types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct trace_Nx_t { unsigned int N; unsigned int x; };
struct trace_ab_t { int64_t a; uint64_t b; };
struct trace_ij_t { int i; unsigned int j; };

struct trace_Nx_t trace_Nx;
struct trace_ab_t trace_ab;
struct trace_ij_t trace_ij;

static inline int trace_on_spot_N(unsigned int N) {
    if (trace_Nx.x == UINT_MAX) return 0;
    return N == trace_Nx.N;
}

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

static inline int trace_on_spot_ab(int64_t a, uint64_t b) {
    return a == trace_ab.a && b == trace_ab.b;
}

static inline int trace_on_spot_ij(int i, unsigned int j) {
    return i == trace_ij.i && j == trace_ij.j;
}

extern void test_divisible_x (const fbprime_t p, const unsigned long x, const int n,
		       sieve_info_srcptr si, int side);
extern void trace_update_conditions(sieve_info_srcptr si);

#if defined(TRACK_CODE_PATH) && (defined(TRACE_AB) || defined(TRACE_IJ) || defined(TRACE_Nx))
extern int test_divisible(where_am_I_ptr w);
#else
static inline int test_divisible(where_am_I_ptr w MAYBE_UNUSED) { return 1; }
#endif

#ifdef TRACE_K
void sieve_decrease_logging(unsigned char *S, const unsigned char logp, where_am_I_ptr w);
void sieve_decrease(unsigned char *S, const unsigned char logp, where_am_I_ptr w);
#else   /* TRACE_K */
#ifdef CHECK_UNDERFLOW
void sieve_decrease_underflow_trap(unsigned char *S, const unsigned char logp, where_am_I_ptr w);
#endif  /* CHECK_UNDERFLOW */
static inline void sieve_decrease(unsigned char *S, const unsigned char logp, where_am_I_ptr w MAYBE_UNUSED)
{
#ifdef CHECK_UNDERFLOW
    if (*S < logp)
    sieve_decrease_underflow_trap(S, logp, w);
#endif  /* CHECK_UNDERFLOW */
    *S -= logp;
}
#endif  /* TRACE_K */

#ifdef __cplusplus
}
#endif

#endif	/* LAS_DEBUG_H_ */
