#include "cado.h"
#include <stdio.h>

#include "las-config.h"

void las_display_config_flags(FILE * stream)
{
    fprintf(stream, "# las.c flags:");
#ifdef SAFE_BUCKETS
    fprintf(stream, " SAFE_BUCKETS");
#endif
#ifdef BUCKETS_IN_ONE_BIG_MALLOC
    fprintf(stream, " BUCKETS_IN_ONE_BIG_MALLOC");
#endif
#ifdef PROFILE
    fprintf(stream, " PROFILE");
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
#ifdef ALG_LAZY
    fprintf(stream, " ALG_LAZY");
    fprintf(stream, " NORM_STRIDE=8 (locked)");
    fprintf(stream, " VERT_NORM_STRIDE=%u (max)", VERT_NORM_STRIDE);
#endif
#ifdef ALG_RAT
    fprintf(stream, " ALG_RAT");
#endif
    fprintf(stream, " NORM_BITS=%u", NORM_BITS);
    fprintf(stream, " LOG_BUCKET_REGION=%u", LOG_BUCKET_REGION);
    fprintf(stream, " GUARD=%1.2f", (double) GUARD);
    fprintf(stream, " LOG_MAX=%.1f", LOG_MAX);
    fprintf(stream, " NB_CURVES=%u", NB_CURVES);
    fprintf(stream, "\n");
}				/* }}} */

