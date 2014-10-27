#ifndef CONVEX_HULL
#define CONVEX_HULL

#include "tab_point.h"

#define PI  3.141

//we suppose that the coordinates of pt2 are highter than those of pt1
double compute_angle(point_t * pt1, point_t * pt2);

//return index of the new point of convex hull
int select_next_point(tabular_point_t * t, point_t * pt);

int search_init_point(tabular_point_t * t);

tabular_point_t *convex_hull(tabular_point_t * t);

#endif				/* CONVEX_HULL */
