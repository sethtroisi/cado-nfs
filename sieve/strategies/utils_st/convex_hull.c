#include <stdlib.h>

#include "convex_hull.h"
#include "math.h"


/*
  The following model has been simplified, because the first point is (0,0) : 
  the most left and down point.
*/

static double
convertEnDegre(const double angle){
  return angle*180/PI;
}

double CST_MUL = 1;

//we suppose that the coordinates of pt2 are highter than those of pt1
double
compute_angle (point_t *pt1, point_t *pt2)
{
  double op = pt2->y - pt1->y; 
  double adj = (pt2->x - pt1->x) * CST_MUL;
 
  if(adj==0){
    if(op>0)
      return 90;
    else       
      return 360;
  }
  else if(op==0){
    if(adj>0)
      return 0;
    else
      return 180;
  }
  double angle = convertEnDegre (atan (op/adj));
  return angle;
}

//return index of the new point of convex hull
int
select_next_point (tabular_point_t* t, point_t* pt)
{
  int res = -1;
  int len = tabular_point_get_index (t);
  point_t* elem;
  double angle_min = 90;
  double angle;

  for (int i = 0; i < len; i++)
    {
      elem = tabular_point_get_point (t, i);
      if ( pt->x < elem->x)
	{
	  angle = compute_angle (pt, elem);
	  if (angle < angle_min)
	    {
	      angle_min = angle;
	      res = i;
	    }
	}
    }
  return res;
}

int
search_init_point (tabular_point_t* t)
{
  int res = 0;
  int len = tabular_point_get_index (t);
  point_t* pt_min = tabular_point_get_point (t, res);
  point_t* elem;
  for (int i = 1; i < len; i++)
    {
      elem = tabular_point_get_point (t, i);
      if ( pt_min->x > elem->x && pt_min->y > elem->y)
	{
	  res = i;
	  pt_min = elem;
	}
    }
  return res;
}

tabular_point_t*
convex_hull (tabular_point_t* t)
{
  //max et min!
  int len = t->index;
  double minx = INFINITY;
  double maxx = 0;
  double miny = INFINITY;
  double maxy = 0;
  for (int i =0; i < len; i++)
    {
     point_t* p = tabular_point_get_point (t, i);
     if (p->x < minx)
       minx = p->x;
     if (p->x > maxx)
       maxx = p->x;
     if (p->y < miny)
       miny = p->y;
     if (p->y > maxy)
       maxy = p->y;
    }
  double scaley = log(maxy - miny)/log(10);
  double scalex = log(maxx - minx)/log(10);
  CST_MUL = pow(10, scaley - scalex);

  tabular_point_t* convex_hull = tabular_point_create ();
  int index = search_init_point (t);

  while (index != -1)
    {
      tabular_point_add_point (convex_hull, t->tab[index]);
      index = select_next_point (t, t->tab[index]);
    }

  return convex_hull;
}




/* /\* ______________ test convex hull __________________ *\/ */
/* int main() */
/* { */
/*   //srand(time(NULL)); */
/*   //init file */
/*   FILE *graph_file; */
/*   char init[6] = "lata1"; */
/*   graph_file = fopen (init, "w"); */

/*   tabular_point_t* list = tabular_point_create (); */

/*   int j = 0; */
/*   while (j < 100) */
/*     { */
/*       int num = rand() %100; */
/*       double y = rand() %1000; */
/*       double x = rand() %100; */

/*       tabular_point_add (list, num, x, y); */

/*       fprintf(graph_file,  "%f \t %f\n",x,y); */
/*       j++; */
/*     } */
/*   fclose(graph_file); */
/*   /\*________test _______*\/ */
/*   char init2[6] = "lata2"; */
/*   graph_file = fopen (init2, "w"); */

/*   printf ("LEN = %d\n", tabular_point_get_index (list)); */
/*   tabular_point_t* res = convex_hull (list); */

/*   int len = tabular_point_get_index (res); */
/*   point_t *el; */
/*   for (int i=0; i < len; i++) */
/*     { */
/*       el = tabular_point_get_point(res, i); */
/*       fprintf(graph_file,  "%f \t %f\n", el->x, el->y); */
/*     } */
/*   tabular_point_free (res); */
/*   tabular_point_free (list); */

/*   fclose(graph_file); */
/*   exit(1); */
/* } */
/*   /\*______________ test convex hull __________________ *\/ */
