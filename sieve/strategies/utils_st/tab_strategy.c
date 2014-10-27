#include "tab_strategy.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>


tabular_strategy_t*
tabular_strategy_create (void)
{
  tabular_strategy_t* t = malloc(sizeof(*t));
  assert (t != NULL);
  
  t->index = 0;
  t->size = 2;
    
  t->tab = malloc (t->size * sizeof(strategy_t*));
  assert (t->tab != NULL);

  return t;
}

void
tabular_strategy_free (tabular_strategy_t* t)
{
  for (int i =0; i < t->index; i++)
    strategy_free (t->tab[i]);
  free (t->tab);
  free (t);
}

void
tabular_strategy_realloc (tabular_strategy_t *t)
{
  t->tab = realloc(t->tab, t->size*2  *( sizeof(strategy_t*)));
  assert (t->tab!=NULL);
  t->size *=2;
}


tabular_strategy_t*
tabular_strategy_copy (tabular_strategy_t* t)
{
  tabular_strategy_t* res = tabular_strategy_create ();
  int len = t->index;
  for (int i = 0; i < len; i++)
    {
      tabular_strategy_add_strategy (res, t->tab[i]);
    }
  return res;
}


int
tabular_strategy_get_index (tabular_strategy_t* t)
{
  return t->index;
}


strategy_t*
tabular_strategy_get_strategy (tabular_strategy_t *t, int index)
{
  assert (index <= t->index);
  return t->tab[index];
}


void
tabular_strategy_add_strategy (tabular_strategy_t *t, strategy_t* strategy)
{
  if (t->index >= t->size)
    tabular_strategy_realloc (t);
  strategy_t* elem = strategy_copy (strategy);
  t->tab[t->index] = elem;
  t->index++;
}



void
tabular_strategy_concat (tabular_strategy_t* t1, tabular_strategy_t* t2)
{
  int len = t2->index;
  for (int i=0; i < len; i++)
    tabular_strategy_add_strategy (t1, t2->tab[i]);
}

tabular_strategy_t*
tabular_strategy_concat_st (tabular_strategy_t* t1, tabular_strategy_t* t2)
{
  tabular_strategy_t* t = tabular_strategy_create ();
  tabular_strategy_concat (t, t1);
  tabular_strategy_concat (t, t2);
  return t;
}




void
tabular_strategy_print_file (tabular_strategy_t* t, FILE* output_file)
{
  //fprintf (output_file, "\n\n"); 
  for (int i =0; i < t->index; i++)
    strategy_print_file (t->tab[i], output_file);
  //fprintf (output_file, "\n\n ");
}

void
tabular_strategy_print (tabular_strategy_t* t)
{
  tabular_strategy_print_file (t, stdout);
}

 /* //test strategy */
 /*  strategy_t* t = strategy_create (); */
 /*  fm_t* elem = fm_create (); */
 /*  unsigned long value[4] = {1,2,3,4}; */
 /*  double p = 0.5; */
 /*  double a = 2; */
 /*  fm_set_method (elem, value,4); */
 /*  fm_set_proba (elem, &p, 1); */
 /*  fm_set_time (elem, &a, 1); */

 /*  strategy_add_fm (t, elem); */
 /*  p = 2.5; */
 /*  a = 21; */
 /*  fm_set_method (elem, value,4); */
 /*  fm_set_proba (elem, &p, 1); */
 /*  fm_set_time (elem, &a, 1); */
 /*  strategy_add_fm (t, elem); */
 /*  unsigned long value2[6] = {6,5,4,3,2,1}; */
 /*  p = 4.5; */
 /*  a = 221; */
 /*  fm_set_method (elem, value2,6); */
 /*  fm_set_proba (elem, &p, 1); */
 /*  fm_set_time (elem, &a, 1); */

 /*  strategy_add_fm (t, elem); */
 /*  strategy_print (t); */

 /*  tabular_strategy_t* s = tabular_strategy_create (); */
 /*  tabular_strategy_add_strategy (s, t); */
 /*  /\* tabular_strategy_add_strategy (s, t); *\/ */
 /*  tabular_strategy_print (s); */
  
 /*  tabular_strategy_free (s); */
 /*  strategy_free (t); */
 /*  fm_free (elem); */
 /*  exit (1); */

