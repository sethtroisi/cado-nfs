#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h> /* for PRId64 */
#include <ctype.h> /* for isxdigit */
#include <string.h>
#include <errno.h>
#include <algorithm>

#include "relation.h"
#include "gzip.h"
#include "timing.h"
#include "portability.h"
#include "relation-tools.h"
#include "relation.h"

using namespace std;

/*
 * Convention for I/O of rels:
 *   a and b are printed in decimal
 *   primes are printed in hexadecimal.
 *
 */

int
relation::parse(const char *line)
{
    int consumed;
    int side = 0;

    if (sscanf(line, "%" SCNd64 ",%" SCNu64 ":%n", &a, &b, &consumed) < 2)
        return 0;

    sides[0].clear();
    sides[1].clear();

    if (line[consumed] == ':') {
        side++;
        consumed++;
    }

    while(line[consumed] != '\0' && line[consumed] != '\n') {
        unsigned long p;
        int consumed_p;
        if (sscanf(line + consumed, "%lx%n", &p, &consumed_p) < 1)
            return 0;
        add(side, p);
        consumed += consumed_p;
        if (line[consumed] == ',')
            consumed++;
        else if (line[consumed] == ':') {
            side++;
            ASSERT_ALWAYS(side < 2);
            consumed++;
        }
    }
    compress();
    fixup_r();
    return 1;
}

void
relation::print (FILE *file, const char *prefix) const
{
    char buf[RELATION_MAX_BYTES], *p = buf;
    char * fence = buf + sizeof(buf);
    int c;

    c = strlcpy(p, prefix, fence - p);
    p += strnlen(prefix, sizeof(buf));

    c = snprintf(p, fence - p, "%" PRId64 ",%" PRIu64, a, b);
    p += c;

    for(int side = 0 ; side < 2 ; side++) {
        if (p + 1 < fence) *p++ = ':';
        char * op = p;
        for(unsigned int i = 0 ; i < sides[side].size() ; i++) {
            for(int e = sides[side][i].e ; e ; e--) {
                c = gmp_snprintf(p, fence - p, "%Zx,", sides[side][i].p);
                p += c;
            }
        }
        if (p > op) p--;
    }
    if (p + 1 < fence) *p++ = '\n';
    *p = '\0';
    /* print all in one go, this gives us a chance to minimize the risk
     * of mixing up output if we don't pay attention to I/O locking too
     * much in multithreaded context. */
    size_t written = fwrite (buf, 1, p - buf, file);
    if (written != (size_t) (p - buf)) {
        perror("Error writing relation");
        abort();
    }
}

void relation::add(int side, mpz_srcptr p)
{
    /* we have to compute a/b mod p. Since we're not critical here, don't
     * bother.
     */
    if (side == rational_side) {
        add(side, p, 0);
    } else {
        pr x;
        mpz_set(x.p, p);

        mpz_set_uint64(x.r, b);
        if (mpz_invert(x.r, x.r, x.p)) {
            mpz_mul_int64(x.r, x.r, a);
            mpz_mod(x.r, x.r, x.p);
        } else {
            mpz_set(x.r, x.p);
        }

        sides[side].push_back(x);
    }
}

void relation::add(int side, unsigned long p)
{
    if (side == rational_side) {
        add(side, p, 0);
    } else {
        /* use the function provided in relation-tools.c */
        add(side, p, relation_compute_r(a, b, p));
    }
}

void relation::fixup_r(bool also_rational)
{
    for(int side = 0 ; side < 2 ; side++) {
        if (!also_rational && side == rational_side)
            continue;
        for(unsigned int i = 0 ; i < sides[side].size() ; i++) {
            if (mpz_cmp_ui(sides[side][i].r,0) == 0) {
                pr & x(sides[side][i]);

                mpz_set_uint64(x.r, b);
                if (mpz_invert(x.r, x.r, x.p)) {
                    mpz_mul_int64(x.r, x.r, a);
                    mpz_mod(x.r, x.r, x.p);
                } else {
                    mpz_set(x.r, x.p);
                }
            }
        }
    }
}

struct pr_cmp {
    bool operator()(relation::pr const& a, relation::pr const& b) const {
        int c = mpz_cmp(a.p, b.p);
        if (c) { return c < 0; }
        c = mpz_cmp(a.r, b.r);
        return c < 0;
    }
};

inline bool operator==(relation::pr const& a, relation::pr const& b) {
    return mpz_cmp(a.p, b.p) == 0 && mpz_cmp(a.r, b.r) == 0;
}


void relation::compress()
{
    for(int side = 0 ; side < 2 ; side++) {
        vector<pr> & v(sides[side]);
        sort(v.begin(), v.end(), pr_cmp());
        unsigned int j = 0;
        for(unsigned int i = 0; i < v.size() ; j++) {
            if (j < i) {
                v[j] = v[i];
            }
            for(i++ ; i < v.size() && v[i] == v[j] ; i++) {
                v[j].e += v[i].e;
            }
        }
        v.erase(v.begin() + j, v.end());
    }
}

