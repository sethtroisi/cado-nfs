#ifndef INTERSECTIONS_HPP_
#define INTERSECTIONS_HPP_


/*
 * Structure for the intersection between an interval [i0..i1[ and
 * two consecutives fences.
 */
struct isect_info {
    int k;                      // intersection with [fences[k]..fences[k+1][
    unsigned int offset_me;     // offset w.r.t. i0
    unsigned int offset_there;  // offset w.r.t. fences[k]
    unsigned int count;         // size of the intersection
};

/* This internal routine is used only in the case where the row and
 * column permutations are equal.
 */

/*  Computes the intersection of [i0..i1[ with the fences ; we clamp the
 *  interval [i0..i1[ to the maximum value nmax, which must correspond to
 *  the last defined fence.
 */
#ifdef __cplusplus
extern "C" {
#endif

unsigned int intersect_two_intervals(unsigned int * offset_a, unsigned int * offset_b, unsigned int a0, unsigned int a1, unsigned int b0, unsigned int b1);

void intersect(unsigned int * plen, struct isect_info ** res, unsigned int * fences, unsigned int i0, unsigned int i1, unsigned int nmax);

#ifdef __cplusplus
}
#endif

#endif	/* INTERSECTIONS_HPP_ */
