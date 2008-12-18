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

/*  Computes the intersection of [i0..i1[ with the fences */
#ifdef __cplusplus
extern "C" {
#endif

void intersect(unsigned int * plen, struct isect_info ** res, unsigned int * fences, unsigned int x0, unsigned int x1);

#ifdef __cplusplus
}
#endif

#endif	/* INTERSECTIONS_HPP_ */
