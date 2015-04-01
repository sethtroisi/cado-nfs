#include "cado.h"
#include <stdio.h>
#include "verbose.h"

#include "las-config.h"

void las_display_config_flags()
{
    verbose_output_print(0, 1, "# las.c flags:");
#ifdef SAFE_BUCKETS
    verbose_output_print(0, 1, " SAFE_BUCKETS");
#endif
#ifdef BUCKETS_IN_ONE_BIG_MALLOC
    verbose_output_print(0, 1, " BUCKETS_IN_ONE_BIG_MALLOC");
#endif
#ifdef PROFILE
    verbose_output_print(0, 1, " PROFILE");
#endif
#ifdef UNSIEVE_NOT_COPRIME
    verbose_output_print(0, 1, " UNSIEVE_NOT_COPRIME");
#endif
#ifdef WANT_ASSERT_EXPENSIVE
    verbose_output_print(0, 1, " WANT_ASSERT_EXPENSIVE");
#endif
#ifdef TRACE_K
    verbose_output_print(0, 1, " TRACE_K");
#endif
#ifdef TRACK_CODE_PATH
    verbose_output_print(0, 1, " TRACK_CODE_PATH");
#endif
#ifdef SUPPORT_I17
    verbose_output_print(0, 1, " SUPPORT_I17");
#endif
#ifdef ALG_LAZY
    verbose_output_print(0, 1, " ALG_LAZY");
    verbose_output_print(0, 1, " NORM_STRIDE=8 (locked)");
    verbose_output_print(0, 1, " VERT_NORM_STRIDE=%u (max)", VERT_NORM_STRIDE);
#endif
#ifdef ALG_RAT
    verbose_output_print(0, 1, " ALG_RAT");
#endif
#ifdef SUPPORT_LARGE_Q
    verbose_output_print(0, 1, " SUPPORT_LARGE_Q");
#endif
#ifdef SKIP_GCD3
    verbose_output_print(0, 1, " SKIP_GCD3");
#endif
#ifdef USE_CACHEBUFFER
    verbose_output_print(0, 1, " USE_CACHEBUFFER");
#endif
    verbose_output_print(0, 1, " NORM_BITS=%u", NORM_BITS);
    verbose_output_print(0, 1, " LOG_BUCKET_REGION=%u", LOG_BUCKET_REGION);
    verbose_output_print(0, 1, " GUARD=%1.2f", (double) GUARD);
    verbose_output_print(0, 1, " LOG_MAX=%.1f", LOG_MAX);
    verbose_output_print(0, 1, "\n");
}				/* }}} */
