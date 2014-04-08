#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "portability.h"
#include "filter_common.h"
#include "filter_badideals.h"
#include "mod_ul.h"

#ifdef FOR_FFS
#include "utils-ffs.h"
#endif

#define MAXLINE 1024

/*
  A line in badidealinfo has the form:
    p k r v0 v1 ... vs
  If means that if (a/b) mod p^k is r (with the usual convention for
  projective roots, see below), then the s corresponding columns should
  be filled with the values v0, v1, ... vs.

  More precisely, if vi is positive, then this is indeed the value, but if vi
  is negative, then the value (e - vi) should be put in the column, where e
  is the valuation of p in the norm.

  Remarks:
  - there should be a line of the form
     p,(r mod p):side: s
    in the .badideals file, in order to "declare" the appropriate number of
    columns for this (p,r).
  - the badidealinfo is supposed to cover all the cases, but not necessarily
    in a simple way (the set of congruences might involve some "mod p" and
    some "mod p^2" rules, for instance).

  Projective roots:
  If we are in the case where (b/a) == 0 mod p, we write "p^k + (1/r)"
  instead of r.
  E.g.
    2 3 10 -2 2
  means that we are dealing with the case (a:b) = (1:2) mod 2^3, and that in
  that case, we should write (e-2) and 2 in the two corresponding columns.

 */


void read_bad_ideals_info(const char *filename, allbad_info_t info)
{
    FILE *file = fopen(filename, "r");
    ASSERT_ALWAYS(file != NULL);

    info->n = 0;
    info->badid_info = NULL;

    char str[1024], *ptr, *nptr;
    while (fgets(str, MAXLINE, file)) {
        if (feof(file))
            break;
        if (str[0] == '#' || str[0] == '\n')
            continue;
        errno = 0;
        p_r_values_t p = strtoul(str, &ptr, 10);
        ASSERT_ALWAYS(errno == 0);
        unsigned long k = strtoul(ptr, &nptr, 10);
        ASSERT_ALWAYS(errno == 0);
        p_r_values_t rk = strtoul(nptr, &ptr, 10);
        ASSERT_ALWAYS(errno == 0);

        badid_info_struct_t item;
        item.p = p;
        item.k = k;
        item.rk = rk;
        item.pk = p;
        for (unsigned int i = 1; i < k; ++i)
            item.pk *= p;
        if (rk < item.pk)
            item.r = item.rk % p;
        else {
            p_r_values_t x = item.rk-item.pk;
            item.r = p + (x % p);
        }
        item.ncol = 0;
        do {
            long v = strtol(ptr, &nptr, 10);
            ASSERT_ALWAYS(errno == 0);
            if (ptr == nptr)
                break;
            item.val[item.ncol] = v;
            item.ncol++;
            ptr = nptr;
        } while (1);

        info->n++;
        info->badid_info = (badid_info_struct_t *) realloc(
                info->badid_info,
                (info->n)*sizeof(badid_info_struct_t));
        info->badid_info[info->n-1] = item;
    }
    fclose(file);
}

static inline p_r_values_t
compute_r_wrapper (int64_t x, int64_t y, p_r_values_t p)
{
#ifndef FOR_FFS
  if (y < 0)
  {
    p_r_values_t r = relation_compute_r (x, (uint64_t) (-y), p);
    return (r == 0) ? r : p - r;
  }
  else
    return relation_compute_r (x, (uint64_t) y, p);
#else
  return ffs_relation_compute_r (x, (uint64_t) y, p);
#endif
}

static inline int
is_divisible_wrapper (uint64_t b, p_r_values_t p)
{
#ifndef FOR_FFS
  return ((b % p) == 0);
#else
  return ffs_is_zero_mod (b, (uint64_t) p);
#endif
}

void
handle_bad_ideals (int *exp_above, int64_t a, uint64_t b, p_r_values_t p, int e,
                   allbad_info_t info)
{
    p_r_values_t r;
    if (is_divisible_wrapper(b, p)) /* same as ((b % p) == 0) */
        r = p + compute_r_wrapper (b, a, p);
    else
        r = compute_r_wrapper (a, b, p);
    for(int i = 0; i < info->n; ++i) {
        if (p != info->badid_info[i].p)
            continue;
        if (r != info->badid_info[i].r)
            continue;
        p_r_values_t rk;
        p_r_values_t pk = info->badid_info[i].pk;
        if (is_divisible_wrapper(b, p)) /* same as ((b % p) == 0) */
            rk = pk + compute_r_wrapper (b, a, pk);
        else
            rk = compute_r_wrapper (a, b, pk);
        if (rk != info->badid_info[i].rk)
            continue;
        for (unsigned int j = 0; j < info->badid_info[i].ncol; ++j) {
            int v = info->badid_info[i].val[j];
            if (v>=0)
                exp_above[j] = v;
            else {
                ASSERT_ALWAYS(e >= -v);
                exp_above[j] = e+v;
            }
        }
        return;
    }
    ASSERT_ALWAYS(0);
}

