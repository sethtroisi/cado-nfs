#include "tab_strategy.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>

tabular_strategy_t *tabular_strategy_create(void)
{
    tabular_strategy_t *t = malloc(sizeof(*t));
    assert(t != NULL);

    t->index = 0;
    t->size = 2;

    t->tab = malloc(t->size * sizeof(strategy_t *));
    assert(t->tab != NULL);

    return t;
}

void tabular_strategy_free(tabular_strategy_t * t)
{
    if (t != NULL)
	{
	    for (int i = 0; i < t->index; i++)
		strategy_free(t->tab[i]);
	    free(t->tab);
	    free(t);
	}
}

void tabular_strategy_realloc(tabular_strategy_t * t)
{
    t->tab = realloc(t->tab, t->size * 2 * (sizeof(strategy_t *)));
    assert(t->tab != NULL);
    t->size *= 2;
}

tabular_strategy_t *tabular_strategy_copy(tabular_strategy_t * t)
{
    tabular_strategy_t *res = tabular_strategy_create();
    int len = t->index;
    for (int i = 0; i < len; i++) {
	tabular_strategy_add_strategy(res, t->tab[i]);
    }
    return res;
}

int tabular_strategy_get_index(tabular_strategy_t * t)
{
    return t->index;
}

strategy_t *tabular_strategy_get_strategy(tabular_strategy_t * t, int index)
{
    assert(index <= t->index);
    return t->tab[index];
}

void
tabular_strategy_add_strategy(tabular_strategy_t * t, strategy_t * strategy)
{
    if (t->index >= t->size)
	tabular_strategy_realloc(t);
    strategy_t *elem = strategy_copy(strategy);
    t->tab[t->index] = elem;
    t->index++;
}

void tabular_strategy_concat(tabular_strategy_t * t1, tabular_strategy_t * t2)
{
    int len = t2->index;
    for (int i = 0; i < len; i++)
	tabular_strategy_add_strategy(t1, t2->tab[i]);
}

tabular_strategy_t *tabular_strategy_concat_st(tabular_strategy_t * t1,
					       tabular_strategy_t * t2)
{
    tabular_strategy_t *t = tabular_strategy_create();
    tabular_strategy_concat(t, t1);
    tabular_strategy_concat(t, t2);
    return t;
}

/************************************************************************/
/*           PRINT AND SCAN OUR FILES OF FACTORING METHODS              */
/************************************************************************/

void tabular_strategy_fprint(FILE * output_file, tabular_strategy_t * t)
{
    for (int i = 0; i < t->index; i++)
	strategy_fprint (output_file, t->tab[i]);
}

void tabular_strategy_print(tabular_strategy_t * t)
{
    tabular_strategy_fprint (stdout, t);
}


static int
is_number (const char c)
{
  return (c >= 48 && c <= 57);
}

static void
next_number (FILE* file, char* current_char)
{
  //end the current number
  while (is_number (*current_char))
    *current_char = fgetc (file);
  
  //find the next number
  while (*current_char != EOF  && !is_number (*current_char))
    {
      *current_char = fgetc (file);
    }
  fseek(file, -1, SEEK_CUR);
}

static void
check_variables_matrix_strat (const unsigned long type, 
			      const unsigned long curve, 
			      const unsigned long b1, 
			      const unsigned long b2, 
			      const int side)
{
  int res = true;
  //side
  if (side > 2)
    {
      fprintf (stderr, "error with parameter side : %d\n", side);
      res = false;
    }
  //b1, b2
  if (b2 < b1)
    {
      fprintf (stderr, "error with parameter (b1,b2) : (%lu, %lu)\n", b1, b2);
      res = false;
    }
  //type
  if ((int)type > NB_METHOD)
    {
      fprintf (stderr, "error with parameter type : %lu\n", type);
      res = false;
    }
  //curve
  if ((int)curve > NB_CURVE)
    {
      fprintf (stderr, "error with parameter curve : %lu\n", curve);
      res = false;
    }
  if (!res)
      exit (EXIT_FAILURE);
} 


tabular_strategy_t*
tabular_strategy_fscan (FILE* file)
{
  if (file == NULL)
      return NULL;
  
  //allocate tabular_strategy
  tabular_strategy_t* tab = tabular_strategy_create ();

  //collect data
  char current_char = fgetc (file);
  int side;

  while (current_char != EOF)
    {
      //{{collect strat
      strategy_t* strat = strategy_create ();
      //{{{collect fm
      while (is_number (current_char))
	{
	  fm_t* elem = fm_create ();
	  fseek(file, -1, SEEK_CUR);
	  fscanf (file, "%lu", &elem->method[0]);
	  next_number (file, &current_char);

	  fscanf (file, "%lu", &elem->method[1]);
	  next_number (file, &current_char);

	  fscanf (file, "%lu", &elem->method[2]);
	  next_number (file, &current_char);

	  fscanf (file, "%lu", &elem->method[3]);
	  next_number (file, &current_char);

	  fscanf (file, "%d", &side);
	  check_variables_matrix_strat (elem->method[0], elem->method[1], elem->method[2],
					elem->method[3], side); 
	  //go to end of line: 10 = '\t'
	  while (current_char != 10)
	    current_char = fgetc (file);
	  current_char = fgetc (file);

	  strategy_add_fm_side (strat, elem, side);
	  fm_free (elem);
	}
      next_number (file, &current_char);
      fscanf (file, "%lf", &strat->proba);

      next_number (file, &current_char);
      fscanf (file, "%lf", &strat->time);

      next_number (file, &current_char);
      //}}
      //add to the tabular
      tabular_strategy_add_strategy (tab, strat);
      strategy_free (strat);
    }
  fclose(file);
  return tab;
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

//print and scan: test on ecrit et on relit et on doit avoir la meme chose
