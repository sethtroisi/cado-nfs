#include "cado.h"
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "verbose.h"

#define G(X) CADO_VERBOSE_PRINT_ ## X
#define F(X) (UINT64_C(1) << G(X))


const char * verbose_flag_list[] = 
{
    [G(CMDLINE)]                = "print-cmdline",
    [G(MODIFIED_FILES)]         = "print-modified-files",
    [G(COMPILATION_INFO)]       = "print-compilation-info",
    [G(BWC_DISPATCH_SLAVES)]    = "bwc-dispatch-slaves",
    [G(BWC_DISPATCH_MASTER)]    = "bwc-dispatch-master",
    [G(BWC_TIMING_GRIDS)]       = "bwc-timing-grids",
    [G(BWC_ITERATION_TIMINGS)]  = "bwc-iteration-timings",
    [G(BWC_CACHE_BUILD)]        = "bwc-cache-build",
};

struct {
    const char * name;
    uint64_t mask;
} verbose_flag_groups[] = {
    { "all-cmdline",
            F(CMDLINE) |
            F(MODIFIED_FILES) |
            F(COMPILATION_INFO) },
    { "all-bwc-dispatch",
            F(BWC_DISPATCH_SLAVES) |
            F(BWC_DISPATCH_MASTER) |
            F(BWC_CACHE_BUILD) },
    { "all-bwc-sub-timings", 
            F(BWC_TIMING_GRIDS) |
            F(BWC_ITERATION_TIMINGS) },
};


uint64_t verbose_flag_word;

/* This must be called in single-threaded context, preferably at program
 * start */
void verbose_set_enabled_flags(param_list pl)
{
    verbose_flag_word = ~0UL;

    const char * v = param_list_lookup_string(pl, "verbose_flags");
    if (!v) return;

    char * w = strdup(v);
    char * p = w;
    char * q;
    for( ; *p != '\0' ; p = q) {
        q = strchr(p, ',');
        if (q) {
            *q++ = '\0';
        } else {
            q = p + strlen(p);
        }
        int enabled = 1;
        if (strncmp(p, "no-", 3) == 0) { enabled = 0; p += 3; }
        else if (strncmp(p, "no", 2) == 0) { enabled = 0; p += 2; }
        else if (*p == '^') { enabled = 0; p += 1; }
        else if (*p == '!') { enabled = 0; p += 1; }

        uint64_t mask = 0;
        for(size_t i = 0 ; i < sizeof(verbose_flag_list) / sizeof(verbose_flag_list[0]) ; i++) {
            if (strcmp(p, verbose_flag_list[i]) == 0) {
                mask = UINT64_C(1) << (int) i;
                break;
            }
        }
        for(size_t i = 0; i < sizeof(verbose_flag_groups) / sizeof(verbose_flag_groups[0]) ; i++) {
            if (strcmp(p, verbose_flag_groups[i].name) == 0) {
                mask = verbose_flag_groups[i].mask;
                break;
            }
        }
        if (!mask) {
            fprintf(stderr, "Verbose flag not recognized: %s\n", p);
            abort();
        }
        if (enabled) {
            verbose_flag_word |= mask;
        } else {
            verbose_flag_word &= ~mask;
        }
    }
    free(w);
}

