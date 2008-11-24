#ifndef INTERSECTIONS_HPP_
#define INTERSECTIONS_HPP_

struct isect_info {
    int k;
    unsigned int offset_me;
    unsigned int offset_there;
    unsigned int count;
};

/* This internal routine is used only in the case where the row and
 * column permutations are equal.
 */

/*  Computes the intersection of [i0..i1[ with the fences */
#ifdef __cplusplus

#include <vector>
std::vector<isect_info>
intersect(unsigned int * fences, unsigned int x0, unsigned int x1);
#endif


#endif	/* INTERSECTIONS_HPP_ */
