#include "cado.h"
#include <stdio.h>

#include "las-config.h"

void las_display_config_flags(FILE * stream)
{
    fprintf(stream, "# las.c flags:");
#ifdef SSE_NORM_INIT
    fprintf(stream, " SSE_NORM_INIT");
#endif
#ifdef UGLY_HACK
    fprintf(stream, " UGLY_HACK");
#endif
#ifdef SAFE_BUCKETS
    fprintf(stream, " SAFE_BUCKETS");
#endif
#ifdef BUCKETS_IN_ONE_BIG_MALLOC
    fprintf(stream, " BUCKETS_IN_ONE_BIG_MALLOC");
#endif
#ifdef PROFILE
    fprintf(stream, " PROFILE");
#endif
#ifdef COFACTOR_TRICK
    fprintf(stream, " COFACTOR_TRICK");
#endif
#ifdef UNSIEVE_NOT_COPRIME
    fprintf(stream, " UNSIEVE_NOT_COPRIME");
#endif
#ifdef WANT_ASSERT_EXPENSIVE
    fprintf(stream, " WANT_ASSERT_EXPENSIVE");
#endif
#ifdef TRACE_K
    fprintf(stream, " TRACE_K");
#endif
#ifdef TRACK_CODE_PATH
    fprintf(stream, " TRACK_CODE_PATH");
#endif
#ifdef SUPPORT_I17
    fprintf(stream, " SUPPORT_I17");
#endif
#ifdef LAZY_NORMS
    fprintf(stream, " LAZY_NORMS");
    fprintf(stream, " NORM_STRIDE=%u", NORM_STRIDE);
    fprintf(stream, " VERT_NORM_STRIDE=%u", VERT_NORM_STRIDE);
#endif
    fprintf(stream, " NORM_BITS=%u", NORM_BITS);
    fprintf(stream, " LOG_BUCKET_REGION=%u", LOG_BUCKET_REGION);
    fprintf(stream, " BUCKET_LIMIT_FACTOR=%.1f",
            (double) BUCKET_LIMIT_FACTOR);
    fprintf(stream, " BUCKET_LIMIT_ADD=%u", BUCKET_LIMIT_ADD);
    fprintf(stream, " GUARD=%1.2f", GUARD);
    fprintf(stream, " LOG_MAX=%.1f", LOG_MAX);
    fprintf(stream, "\n");
}				/* }}} */

