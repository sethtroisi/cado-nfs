#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <stdlib.h>

#include <sys/stat.h>

#include "macros.h"
#include "utils.h"

static size_t filesize(const char * path)
{
    struct stat sbuf[1];
    int rc = stat(path, sbuf);
    DIE_ERRNO_DIAG(rc < 0, "stat", path);
    return sbuf->st_size;
}

/* Here are several ways to apply a permutation to a vector.
 *
 * The ``direct'' way is the ``untwisting'' one. the permutation perm is
 * represented by the values perm(0) to perm(n-1) as 32-bit unsigned
 * integers within a file (correct endianness assumed). Then on input, a
 * vector (v_{\perm(0)},...,v_{perm(n-1)}) is transformed into
 * (v_0,...,v_{n-1}) ; in other words, entry i on input goes to position
 * perm(i) on output.
 *
 * The common assumption is that a sufficiently large memory area is
 * available to store two vectors, and two permutations. Note that these
 * functions are called within the core lineaer algebra routines at
 * checkpoint times. It follows that checkpointing is not an entirely
 * trivial task, in particular memory-wise. It is not so much of a
 * problem, since checkpoints are seldom performed. Furthermore, items in
 * core memory may be temporarily swapped out at checkpointing time for a
 * relatively low penalty.
 */

/* The current implementation is rather a memory hog, since it allocates
 * space for everyone prior to doing any work.
 */

void apply_perm_byline(void * out, void * in, unsigned int * sigma, size_t ds MAYBE_UNUSED, size_t n)
{
    unsigned int * ssizes = malloc(n * sizeof(unsigned int));
    unsigned int * doffsets = malloc(n * sizeof(unsigned int));
    unsigned int * doffsets2 = malloc(n * sizeof(unsigned int));
    
    ASSERT_ALWAYS(ssizes && doffsets2 && doffsets);

    for(size_t i = 0 ; i < n ; i++) { doffsets[i] = UINT_MAX; }

    size_t o;
    
    o = 0;
    for(size_t i = 0 ; i < n ; i++) {
        size_t s;
        for(s = 1 ; ((char*)in)[o++] != '\n' ; s++) {
            ASSERT(o < ds);
        }
        ssizes[i] = s;
        doffsets[sigma[i]] = s;
    }

    for(size_t i = 0 ; i < n ; i++) { ASSERT(doffsets[i] != UINT_MAX); }
    
    o = 0;
    for(size_t i = 0 ; i < n ; i++) {
        size_t no = o + doffsets[i];
        doffsets[i] = o;
        o = no;
    }
    for(size_t i = 0 ; i < n ; i++) {
        doffsets2[i] = doffsets[sigma[i]];
    }
    free(doffsets);

    size_t i;
    for(i = 0 ; i < n ; i++) {
        size_t s = ssizes[i];
        memcpy(out + doffsets2[i], in + i * s, s);
    }

    free(doffsets2);
    free(ssizes);
}

void apply_perm_binary(void * out, void * in, unsigned int * sigma,
        size_t ds, size_t n)
{
    size_t i;
    size_t bl = ds / n;
    for(i = 0 ; i < n ; i++) {
        memcpy(out + sigma[i] * bl, in + i * bl, bl);
    }
}

void usage()
{
    fprintf(stderr, "Usage: apply-perm [--byline] --in <input> --out <output> --perm <permutation>\n");
    exit(1);
}

int main(int argc, char * argv[])
{
    param_list pl;
    int byline = 0;
    int reverse = 0;

    param_list_init(pl);
    argv++,argc--;
    param_list_configure_knob(pl, "--byline", &byline);
    param_list_configure_knob(pl, "--reverse", &reverse);
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv))
            continue;
        fprintf(stderr, "Unparsed argument %s\n", argv[0]);
        usage();
    }

    const char * nameout = param_list_lookup_string(pl,"out");
    const char * namein = param_list_lookup_string(pl,"in");
    const char * nameperm = param_list_lookup_string(pl,"perm");

    if (!nameperm)
        usage();

    /* read the permutation */
    unsigned int * sigma;
    size_t n;
    {
        const char * path = nameperm;
        size_t ps = filesize(path);
        ASSERT_ALWAYS(ps && (ps % sizeof(unsigned int) == 0));
        n = ps / sizeof(unsigned int);
        FILE * fperm = fopen(path, "r");
        DIE_ERRNO_DIAG(fperm == NULL, "fopen", path);
        sigma = malloc(n * sizeof(unsigned int));
        ASSERT_ALWAYS(sigma);
        size_t nr = fread(sigma, sizeof(unsigned int), n, fperm);
        FATAL_ERROR_CHECK(nr != n, "sigma: short read");
        fclose(fperm);
    }
    if (reverse) {
        unsigned int * np = malloc(n * sizeof(unsigned int));
        ASSERT_ALWAYS(np);
        for(size_t i = 0 ; i < n ; i++) {
            np[sigma[i]] = i;
        }
        free(sigma);
        sigma = np;
    }

    void * in;
    size_t ds;
    {
        /* At the moment, we don't work as a filter, although given the
         * fact that we boldly allocate space for everything, this could
         * be done.
         */
        ds = filesize(namein);
        in = malloc(ds);
        ASSERT_ALWAYS(in);

        FILE * f = fopen(namein, "r");
        DIE_ERRNO_DIAG(f == NULL, "fopen", namein);

        if (!byline) {
            ASSERT_ALWAYS(ds % n == 0);
        }
        size_t nr = fread(in, 1, ds, f);
        FATAL_ERROR_CHECK(nr != ds, "input: short read");

        fclose(f);
    }

    void * out = malloc(ds);
    ASSERT_ALWAYS(out);
    if (byline) {
        apply_perm_byline(out, in, sigma, ds, n);
    } else {
        apply_perm_binary(out, in, sigma, ds, n);
    }

    free(in);
    free(sigma);

    if (nameout && strcmp(nameout, "-") != 0) {
        FILE * f = fopen(nameout, "w");
        DIE_ERRNO_DIAG(f == NULL, "fopen", nameout);
        size_t nr = fwrite(out, 1, ds, f);
        FATAL_ERROR_CHECK(nr != ds, "output: short write");
        fclose(f);
    } else {
        size_t nr = fwrite(out, 1, ds, stdout);
        FATAL_ERROR_CHECK(nr != ds, "output: short write");
    }
    free(out);

    param_list_clear(pl);

    return 0;
}
