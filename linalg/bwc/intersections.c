#include "bwc_config.h"
#include "macros.h"
#include "intersections.h"

unsigned int intersect_two_intervals(unsigned int * offset_a, unsigned int * offset_b, unsigned int a0, unsigned int a1, unsigned int b0, unsigned int b1)
{
    if (a0 >= b1 || b0 >= a1) {
        *offset_a = 0;
        *offset_b = 0;
        return 0;
    }

    int d0 = (int) (b0 - a0);
    int d1 = (int) (b1 - a1);

    unsigned int u0 = 0;
    unsigned int u1 = a1 - a0;
    unsigned int v0 = 0;
    unsigned int v1 = b1 - b0;

    if (d0 > 0) {
        u0 += d0;
        u1 -= d0;
    } else {
        v0 += -d0;
        v1 -= -d0;
    }

    if (d1 < 0) {
        u1 -= -d1;
    } else {
        v1 -= d1;
    }

    /* Our data will be of interest for indices in the
     * range [u0..u0+u1[. It goes to range [v0..v0+v1[
     */

    ASSERT_ALWAYS(u1 == v1);

    *offset_a = u0;
    *offset_b = v0;
    return u1;
}

#if 0
#include <vector>
/*  Computes the intersection of [i0..i1[ with the fences */
std::vector<isect_info>
intersect(unsigned int * fences, unsigned int x0, unsigned int x1, unsigned int nmax)
{
    std::vector<isect_info> res;

    x0 = MIN(x0, nmax);
    x1 = MIN(x1, nmax);

    if (x1 == x0) {
        return res;
    }

    unsigned int k;
    for(k = 0 ; fences[k + 1] <= x0 ; k++) ;
    // BUG_ON(k >= pi.nvjobs * pi.nvcores);
    /* We have fences[k] <= x0 < fences[k+1] */
    for( ; fences[k] < x1 ; k++) {
        ASSERT(x0 < fences[k+1]);
        ASSERT(fences[k] < x1);

        isect_info r;
        r.k = k;
        r.count = intersect_two_intervals(&r.offset_me, &r.offset_there,
                x0, x1, fences[k], fences[k+1]);

        res.push_back(r);
    }
    return res;
}

void intersect(unsigned int * plen, struct isect_info ** res, unsigned int * fences, unsigned int x0, unsigned int x1, unsigned int nmax)
{
    std::vector<isect_info> vv;
    vv = intersect(fences, x0, x1, nmax);
    *plen = vv.size();
    *res = (struct isect_info *) malloc(*plen * sizeof(struct isect_info));   
    std::copy(vv.begin(), vv.end(), *res);
}
#endif
