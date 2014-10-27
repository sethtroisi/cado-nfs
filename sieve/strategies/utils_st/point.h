#ifndef POINT
#define POINT

#include <stdio.h>

typedef struct point {
	int number;
	double x;
	double y;
} point_t;

point_t *point_create(int numero, double x, double y);

void point_free(point_t * t);

int point_get_number(point_t * t);

double point_get_x(point_t * t);

double point_get_y(point_t * t);

#endif				/* POINT */
