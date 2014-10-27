#include "point.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


point_t*
point_create (int numero, double x, double y)
{
  point_t* t = malloc(sizeof(*t));
  assert (t != NULL);
  t->number = numero;
  t->x = x;
  t->y = y;
  return t;
}

void
point_free (point_t* t)
{
  free (t);
}

int
point_get_number (point_t* t)
{
  return t->number;
}
double
point_get_x (point_t* t)
{
  return t->x;
}

double
point_get_y (point_t* t)
{
  return t->y ;
}
