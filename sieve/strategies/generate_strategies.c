#include "generate_strategies.h"
#include "facul.h"
#include "ecm.h"

#include <math.h>
#include <stdlib.h>
#include <float.h>

/*todo: I know :)! In few days, you will never see the words
  rationnal and algebraic :)
*/
#define RATIONNAL_SIDE 0
#define ALGEBRAIC_SIDE 1
static double EPSILON_DBL = LDBL_EPSILON;


static int
is_good_decomp (decomp_t * dec, int len_p_min, int len_p_max)
{
  int len = dec->len;
  for (int i = 0; i < len; i++)
    if (dec->tab[i] > len_p_max || dec->tab[i] < len_p_min)
      return false;
  return true;
}




/************************************************************************/
/*                      COLLECT DATA FOR ONLY ONE COFACTOR              */
/************************************************************************/



/*
  Compute the percent of good decomposition which is detected by the strategy.
 */
double
compute_proba_strategy (tabular_decomp_t * init_tab, strategy_t * strat,
			int len_p_min, int len_p_max)
{
  double all = 0.0;
  double nb_found_elem = 0;
  int nb_fm = tabular_fm_get_index (strat->tab_fm);
  int nb_decomp = init_tab->index;
  tabular_fm_t *tab_fm = strat->tab_fm;

  for (int index_decomp = 0; index_decomp < nb_decomp; index_decomp++)
    {
      decomp_t *dec = init_tab->tab[index_decomp];
      if (is_good_decomp (dec, len_p_min, len_p_max))
	{
	  double res = 1;
	  for (int index_fm = 0; index_fm < nb_fm; index_fm++)
	    {
	      fm_t *elem = tab_fm->tab[index_fm];
	      double *proba = fm_get_proba (elem);
	      int len_proba = fm_get_len_proba (elem);
	      double proba_method = 1;
	      for (int i = 0; i < dec->len; i++)
		{
		  int j = dec->tab[i] - elem->len_p_min;
		  if (j < 0)
		    {
		      printf ("%d, %d\n", dec->tab[i], elem->len_p_min);
		      fprintf (stderr, "ATTENTION A NOTRE VALEUR DE LB\n");
		      exit (1);
		    }

		  if (j < len_proba)
		    {
		      proba_method *= (1 - proba[j]);
		    }
		  //else //the probability seems closest to 0.
		}
	      res *= proba_method;
	      if (elem->method[0] == PM1_METHOD ||
		  elem->method[0] == PP1_27_METHOD ||
		  elem->method[0] == PP1_65_METHOD)
		//because if you chain PP1||PM1 to PM1||PP1-->they are
		//not independant.
		res = (proba_method + res) / 2;
	    }
	  nb_found_elem += (1 - res) * dec->nb_elem;
	}
      all += dec->nb_elem;
    }
  return nb_found_elem / all;
}


/*
Compute the mean cost of the strategy! 
 */
double
compute_time_strategy (tabular_decomp_t * init_tab, strategy_t * strat)
{

  double cost = 0;
  double all = 0;
  int nb_fm = tabular_fm_get_index (strat->tab_fm);
  int nb_decomp = init_tab->index;


  for (int index_decomp = 0; index_decomp < nb_decomp; index_decomp++)
    {
      decomp_t *dec = init_tab->tab[index_decomp];
      double t = 0, p = 0;
      double num = 0, denom = 0;
      double nb_elem = dec->nb_elem;	//the number of elements remaining.
      //compute the time of each decomposition
      for (int index_fm = 0; index_fm < nb_fm; index_fm++)
	{
	  num = 0;
	  denom = 0;
	  tabular_fm_t *tab_fm = strat->tab_fm;
	  fm_t *elem = tab_fm->tab[index_fm];
	  double *time = fm_get_time (elem);
	  int len_time = fm_get_len_time (elem);
	  double *proba = fm_get_proba (elem);
	  int len_proba = fm_get_len_proba (elem);
	  double moy_proba = 1;
	  //compute the time of each decomposition and for one FM
	  for (int i = 0; i < dec->len; i++)
	    {
	      //temps
	      int j = dec->tab[i] - elem->len_p_min;
	      if (j < 0)
		{
		  fprintf (stderr, "ATTENTION A NOTRE VALEUR DE LB\n");
		  exit (1);
		}
	      if (j >= len_time)
		t = time[len_time - 1];
	      else
		t = time[j];
	      if (j >= len_proba)
		p = 0;
	      else
		p = proba[j];
	      denom += p;
	      num += t * p;
	      //moy_proba
	      moy_proba *= (1 - p);
	    }

	  if (denom < EPSILON_DBL)	//denom == 0
	    {
	      cost += t * nb_elem;
	      //here, the method give zero succes percent for all 
	      //factors and the cost is the same for all of this. So cost = t
	    }
	  else
	    {
	      //num/denom = mean cost of the method[index_fm]
	      cost += num / denom * nb_elem;
	    }
	  nb_elem *= moy_proba;
	  if (elem->method[0] == PM1_METHOD ||
	      elem->method[0] == PP1_27_METHOD ||
	      elem->method[0] == PP1_65_METHOD)
	    {
	      nb_elem = (dec->nb_elem * moy_proba + nb_elem) / 2;
	    }
	}
      all += dec->nb_elem;
    }
  return cost / all;
}


/*
  As its name suggests, this function add one strategy to our array
  't' without the zero method.
*/
static void
tabular_strategy_add_strategy_without_zero (tabular_strategy_t * t,
					    strategy_t * strategy)
{
  if (t->index >= t->size)
    tabular_strategy_realloc (t);

  strategy_t *elem = strategy_create ();
  int len = strategy->tab_fm->index;
  int strat_is_zero = true;
  for (int i = 0; i < len; i++)
    {
      if (!tabular_fm_is_zero (strategy->tab_fm, i))
	{
	  strat_is_zero = false;
	  tabular_fm_add_fm (elem->tab_fm, strategy->tab_fm->tab[i]);
	}
    }
  if (strat_is_zero)
    tabular_fm_add_fm (elem->tab_fm, strategy->tab_fm->tab[0]);

  elem->proba = strategy->proba;
  elem->time = strategy->time;
  t->tab[t->index] = elem;
  t->index++;
}



/************************************************************************/
/*                   GENERATE MATRIX                                    */
/************************************************************************/

static void
generate_collect_iter_ecm_rc (fm_t * zero, tabular_fm_t * ecm_rc,
			      int ind_ecm_rc, strategy_t * strat, int ind_tab,
			      int index_iter, int len_iteration,
			      tabular_decomp_t * init_tab,
			      tabular_strategy_t * res, int lb, int ub)
{
  if (index_iter >= len_iteration)
    {
      tabular_strategy_add_strategy_without_zero (res, strat);
      int nb_strat = res->index - 1;
      double proba =
	compute_proba_strategy (init_tab, res->tab[nb_strat], lb, ub);
      double time = compute_time_strategy (init_tab, res->tab[nb_strat]);
      strategy_set_proba (res->tab[nb_strat], proba);
      strategy_set_time (res->tab[nb_strat], time);
    }
  else
    {
      for (int i = ind_ecm_rc; i <= ecm_rc->index; i++)
	{
	  if (i < ecm_rc->index)
	    tabular_fm_set_fm_index (strat->tab_fm, ecm_rc->tab[i], ind_tab);
	  else			// == ecm_rc->index
	    tabular_fm_set_fm_index (strat->tab_fm, zero, ind_tab);
	  generate_collect_iter_ecm_rc (zero, ecm_rc, i, strat, ind_tab + 1,
					index_iter + 1, len_iteration,
					init_tab, res, lb, ub);
	}
    }
}


/*
  This function allows to generate strategies, compute their data, and
  select the convex hull of these strategies for one size of cofactor
  'r'.  It allows to avoid to full the RAM when we generated
  strategies!
*/
tabular_strategy_t *
generate_strategies_oneside (fm_t * zero, tabular_fm_t * pm1,
			     tabular_fm_t * pp1, tabular_fm_t * ecm_m16,
			     tabular_fm_t * ecm_rc, int lb, int ub, int r)
{
  //tabular_decomp_t *init_tab = decomposition_of_cofactor (lb, r);
  char name_file[200];
  sprintf (name_file, "/localdisk/trichard/strategies/decomp_cofactor/decomp_%d_%d", r, lb);
  FILE *file = fopen (name_file, "r");

  tabular_decomp_t *init_tab = tabular_decomp_fscan (file);

  if (init_tab == NULL)
    {
      fprintf (stderr, "impossible to read '%s'\n", name_file);
      exit (EXIT_FAILURE);
    }
  fclose (file);

  tabular_strategy_t *all_strat = tabular_strategy_create ();
  //contains strategies which will be processed.

  tabular_strategy_t *res = tabular_strategy_create ();
  //contains the final result

  int len_pm1 = pm1->index;
  int len_pp1 = pp1->index;
  int len_ecm_m16 = ecm_m16->index;
  int len_ecm_rc = ecm_rc->index;
  int nb_iteration = 4;

  //init strat
  int len_strat = 3 + nb_iteration;
  strategy_t *strat = strategy_create ();
  for (int i = 0; i < len_strat; i++)
    strategy_add_fm (strat, zero);

  tabular_fm_t *tab_strat = strat->tab_fm;
  int nb_strat = 0;

  for (int ind_pm1 = -1; ind_pm1 < len_pm1; ind_pm1++)
    {
      double current_proba_pm1 =
	(ind_pm1 == -1) ? 0 : pm1->tab[ind_pm1]->proba[0];
      if (ind_pm1 != -1)
	tabular_fm_set_fm_index (tab_strat, pm1->tab[ind_pm1], 0);
      else
	tabular_fm_set_fm_index (tab_strat, zero, 0);
      int ind_pp1 = -1;
      while (ind_pp1 < len_pp1)
	{
	  if (ind_pp1 != -1)
	    {
	      while (ind_pp1 < len_pp1 &&
		     pp1->tab[ind_pp1]->proba[0] < current_proba_pm1)
		ind_pp1++;
	      if (ind_pp1 >= len_pp1)
		break;
	    }
	  if (ind_pp1 != -1)
	    tabular_fm_set_fm_index (tab_strat, pp1->tab[ind_pp1], 1);
	  else
	    tabular_fm_set_fm_index (tab_strat, zero, 1);

	  for (int ind_ecm_m16 = -1; ind_ecm_m16 < len_ecm_m16; ind_ecm_m16++)
	    {
	      if (ind_ecm_m16 != -1)
		tabular_fm_set_fm_index (tab_strat, ecm_m16->tab[ind_ecm_m16],
					 2);
	      else
		tabular_fm_set_fm_index (tab_strat, zero, 2);
	      double current_proba_ecm =
		(ind_ecm_m16 == -1) ? 0 : ecm_m16->tab[ind_ecm_m16]->proba[0];

	      int ind_ecm_rc = 0;
	      while (ind_ecm_rc < len_ecm_rc &&
		     ecm_rc->tab[ind_ecm_rc]->proba[0] < current_proba_ecm)
		ind_ecm_rc++;

	      generate_collect_iter_ecm_rc (zero, ecm_rc, ind_ecm_rc, strat,
					    3, 0, nb_iteration, init_tab,
					    all_strat, lb, ub);

	      //process data to avoid to full the RAM!!!
	      //{{
	      nb_strat = all_strat->index - 1;
	      if (nb_strat > 100000)
		{
		  //add strategies of the old convexhull and recompute
		  //the new convex hull
		  tabular_strategy_concat (all_strat, res);
		  tabular_strategy_free (res);
		  res = convex_hull_strategy (all_strat);
		  //clear previous collect and start a new collect.
		  tabular_strategy_free (all_strat);
		  all_strat = tabular_strategy_create ();
		  nb_strat = 0;
		}
	      //}}
	    }
	  ind_pp1++;
	}
    }
  if (nb_strat != 0)		//process data
    {
      //add strategies of the old convexhull and recompute the new
      //convex hull
      tabular_strategy_concat (all_strat, res);
      tabular_strategy_free (res);
      res = convex_hull_strategy (all_strat);
    }
  tabular_strategy_free (all_strat);
  strategy_free (strat);
  tabular_decomp_free (init_tab);
  return res;
}

/*
  As its name suggests, this function concatenates two strategies but
  take into account the side of each strategy: st1 will be the
  first_side and st2 the other side.
*/

static strategy_t *
first_concat_strategy (strategy_t * st1, strategy_t * st2, int first_side)
{
  int side = (first_side == 0) ? RATIONNAL_SIDE : ALGEBRAIC_SIDE;
  strategy_t *st = strategy_create ();
  int len1 = st1->tab_fm->index;
  int len2 = st2->tab_fm->index;
  st->len_side = len1 + len2;
  if (st->side == NULL)
    st->side = malloc (sizeof (int) * (st->len_side));
  else
    st->side = realloc (st->side, sizeof (int) * (st->len_side));

  for (int i = 0; i < len1; i++)
    {
      strategy_add_fm (st, st1->tab_fm->tab[i]);
      st->side[i] = side;
    }
  side = first_side ? RATIONNAL_SIDE : ALGEBRAIC_SIDE;
  for (int i = 0; i < len2; i++)
    {
      strategy_add_fm (st, st2->tab_fm->tab[i]);
      st->side[len1 + i] = side;
    }
  return st;
}


/*
  returns the best strategies to factor a couple (r1, r2), from a set
  of factoring methods for each side. Note that, the probability and
  the time o find a non-trivial for its side must be previously
  computed!
 */
tabular_strategy_t *
generate_strategy_r1_r2 (tabular_strategy_t * strat_r1,
			 tabular_strategy_t * strat_r2)
{
  int len_r1 = strat_r1->index;
  int len_r2 = strat_r2->index;

  tabular_strategy_t *strat_r1_r2 = tabular_strategy_create ();
  tabular_strategy_t *ch = tabular_strategy_create ();
  unsigned long nb_strat = 0;

  int r_init = 0;
  int a_init = 0;
  /*
     for each array of strategies, the first one is the zero strategy.
     There are two cases : 
     -first, we have a success at 0%.
     -Secondly, we have a success at 100%, because the cofactor is already
     prime and not too big.  The four lines below allow to avoid to
     have several unless methods with a zero probability.
   */

  if (strat_r1->tab[0]->proba < EPSILON_DBL)
    r_init = 1;
  if (strat_r2->tab[0]->proba < EPSILON_DBL)
    a_init = 1;

  // the first strategy is the zero strategy, so we add there directly
  // in the convex hull.
  strategy_t *st;
  if (strat_r1->tab[0]->proba < EPSILON_DBL ||
      strat_r2->tab[0]->proba < EPSILON_DBL)
    {
      st = first_concat_strategy (strat_r1->tab[0], strat_r2->tab[0],
				  RATIONNAL_SIDE);
      tabular_strategy_add_strategy (ch, st);
      strategy_free (st);
    }

  for (int r = r_init; r < len_r1; r++)	//rationnal side
    {
      double p1 = strat_r1->tab[r]->proba;
      double c1 = strat_r1->tab[r]->time;
      for (int a = a_init; a < len_r2; a++)	//algebraic side : len
	{
	  nb_strat++;
	  //compute success ans cost:
	  double p2 = strat_r2->tab[a]->proba;
	  double c2 = strat_r2->tab[a]->time;
	  double proba = p1 * p2;
	  double tps1 = c1 + p1 * c2;
	  //mean time when we begin by the RATIONNAL_SIDE
	  double tps2 = c2 + p2 * c1;
	  //mean time when we begin by the ALGEBRAIC_SIDE
	  if (tps1 < tps2)
	    {
	      st = first_concat_strategy (strat_r1->tab[r], strat_r2->tab[a],
					  RATIONNAL_SIDE);
	      strategy_set_proba (st, proba);
	      strategy_set_time (st, tps1);
	    }
	  else
	    {
	      st =
		first_concat_strategy (strat_r2->tab[a], strat_r1->tab[r],
				       ALGEBRAIC_SIDE);
	      strategy_set_proba (st, proba);
	      strategy_set_time (st, tps2);
	    }
	  tabular_strategy_add_strategy (strat_r1_r2, st);
	  strategy_free (st);
	}

      //process data to avoid to full the RAM!!!
      if (nb_strat > 100000 || r == (len_r1 - 1))
	{
	  //add strategies of the old convexhull and recompute the new
	  //convex hull

	  tabular_strategy_concat (strat_r1_r2, ch);
	  tabular_strategy_free (ch);
	  /* printf( "avant CH\n"); */
	  /* tabular_strategy_print (strat_r1_r2); */

	  ch = convex_hull_strategy (strat_r1_r2);
	  //clear previous collect and start a new collect.
	  tabular_strategy_free (strat_r1_r2);
	  strat_r1_r2 = tabular_strategy_create ();
	  nb_strat = 0;
	}
    }

//free
  tabular_strategy_free (strat_r1_r2);

  return ch;
}



tabular_strategy_t ***
generate_matrix (tabular_fm_t * pm1, tabular_fm_t * pp1,
		 tabular_fm_t * ecm_m16, tabular_fm_t * ecm_rc,
		 int rlb, int rub, int rmfb, int alb, int aub, int amfb)
{
  /*
     allocate the matrix res which contains all good strategies for each couple (r_1, r_2).
     r_1 is the lenght of the cofactor in the rationnal side.
     r_2 is the lenght of the cofactor in the algebraic side.
   */
  tabular_strategy_t ***matrix = malloc (sizeof (*matrix) * (rmfb + 1));
  assert (matrix != NULL);
  for (int r1 = 0; r1 <= rmfb; r1++)
    {
      matrix[r1] = malloc (sizeof (*matrix[r1]) * (amfb + 1));
      assert (matrix[r1] != NULL);
    }

  printf ("\n COLLECT DATA\n\n");
  /*
     for each r_1 in [lbr..rmfb], we precompute the data for all
     strategy.  for each r_2, the values will be computed
     sequentially.  define two zero strategies to manage the case
     where r is prime or the case r is impossible (r<lb): - st_zero_0
     : process a zero method and the success is null because r<lb or r
     is a big prime (ub < r < lb^2) - st_zero_1 : process a zero
     method and the success is equal to 100% because r=1, and so this
     cofactor is already factorised.
   */
  fm_t *zero = fm_create ();
  unsigned long method_zero[4] = { 0, 0, 0, 0 };
  fm_set_method (zero, method_zero, 4);

  strategy_t *st_zero_0 = strategy_create ();
  strategy_add_fm (st_zero_0, zero);
  strategy_set_proba (st_zero_0, 0.0);
  strategy_set_time (st_zero_0, 0.0);

  strategy_t *st_zero_1 = strategy_create ();
  strategy_add_fm (st_zero_1, zero);
  strategy_set_proba (st_zero_1, 1.0);
  strategy_set_time (st_zero_1, 0.0);


  tabular_strategy_t **data_rat = malloc (sizeof (*data_rat) * (rmfb + 1));
  int lim_rat = rlb * 2 - 1;
  for (int r1 = 0; r1 <= rmfb; r1++)
    {
      /*
         process the collect for the different cases.
       */
      tabular_strategy_t *strat_r1;
      if (r1 < lim_rat)		//in this case, r1 is prime or equal
	//to 1. So no cofactorization is
	//necessary.
	{
	  strat_r1 = tabular_strategy_create ();
	  if (r1 != 1 && (r1 < rlb || r1 > rub))
	    tabular_strategy_add_strategy (strat_r1, st_zero_0);
	  else			// r1==1 or rlb<= r1 <= rub 
	    tabular_strategy_add_strategy (strat_r1, st_zero_1);
	}
      else
	{			//in this case, the cofactor is
	  //composite. So let's go to work!
	  strat_r1 =
	    generate_strategies_oneside (zero, pm1, pp1, ecm_m16,
					 ecm_rc, rlb, rub, r1);
	}
      data_rat[r1] = strat_r1;
    }
  //todo: warning duplicate code!!!
  printf ("\n COLLECT DATA : step 2\n\n");
  /*
     read good elements for r_2 in the array data_r2 and compute the
     data for each r_1. So :
   */
  int lim_alg = alb * 2 - 1;
  for (int r2 = 0; r2 <= amfb; r2++)
    {
      tabular_strategy_t *strat_r2;
      if (r2 < lim_alg)		//in this case, r2 is prime or equal
	//to 1. So no cofactorization is
	//necessary.
	{
	  strat_r2 = tabular_strategy_create ();
	  if (r2 != 1 && (r2 < alb || r2 > aub))
	    tabular_strategy_add_strategy (strat_r2, st_zero_0);
	  else			// r2==1 or alb<= r2 < aub
	    tabular_strategy_add_strategy (strat_r2, st_zero_1);
	}
      else
	{			//in this case the cofactor is composite.
	  strat_r2 =
	    generate_strategies_oneside (zero, pm1, pp1, ecm_m16,
					 ecm_rc, alb, aub, r2);
	}
      for (int r1 = 0; r1 <= rmfb; r1++)
	{
	  tabular_strategy_t *res =
	    generate_strategy_r1_r2 (data_rat[r1], strat_r2);
	  matrix[r1][r2] = res;
	}
      tabular_strategy_free (strat_r2);
    }

  //free
  for (int r1 = 0; r1 <= rmfb; r1++)
    tabular_strategy_free (data_rat[r1]);
  free (data_rat);

  strategy_free (st_zero_0);
  strategy_free (st_zero_1);
  fm_free (zero);
  return matrix;
}







/************************************************************************/
/*                      CONVEX_HULL_FM                                  */
/************************************************************************/


/*
  These following functions allow to use the module convex_hull, to
  compute the convex hull of a set of strategies.
 */


tabular_point_t *
convert_tab_point_to_tab_strategy (tabular_strategy_t * t)
{
  tabular_point_t *res = tabular_point_create ();
  int len = t->index;
  strategy_t *elem;
  for (int i = 0; i < len; i++)
    {
      elem = t->tab[i];
      tabular_point_add (res, i, elem->proba, elem->time);
    }
  return res;
}

tabular_strategy_t *
convert_tab_strategy_to_tab_point (tabular_point_t * t,
				   tabular_strategy_t * init)
{
  tabular_strategy_t *res = tabular_strategy_create ();
  int len = tabular_point_get_index (t);
  for (int i = 0; i < len; i++)
    {
      int index = point_get_number (tabular_point_get_point (t, i));
      tabular_strategy_add_strategy (res, init->tab[index]);
    }
  return res;
}

tabular_strategy_t *
convex_hull_strategy (tabular_strategy_t * t)
{
  tabular_point_t *tmp = convert_tab_point_to_tab_strategy (t);
  tabular_point_t *res = convex_hull (tmp);
  tabular_strategy_t *res_strat = convert_tab_strategy_to_tab_point (res, t);
  tabular_point_free (tmp);
  tabular_point_free (res);
  return res_strat;
}
