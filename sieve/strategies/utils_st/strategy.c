#include "cado.h"
#include "portability.h"
#include "utils.h"

#include "strategy.h"

#include <stdio.h>
#include <stdlib.h>

strategy_t *strategy_create()
{
    strategy_t *t = malloc(sizeof(*t));
    ASSERT_ALWAYS(t != NULL);
    t->tab_fm = tabular_fm_create();
    t->proba = 0;
    t->time = 0;
    t->side = NULL;
    t->len_side = 0;
    return t;
}

void strategy_free(strategy_t * t)
{
    if (t != NULL) {
	tabular_fm_free(t->tab_fm);
	if (t->side != NULL)
	    free(t->side);
	free(t);
    }
}

tabular_fm_t *strategy_get_tab_fm(strategy_t * t)
{
    return t->tab_fm;
}

double strategy_get_proba(strategy_t * t)
{
    return t->proba;
}

double strategy_get_time(strategy_t * t)
{
    return t->time;
}

void strategy_set_proba(strategy_t * t, double proba)
{
    t->proba = proba;
}

void strategy_set_time(strategy_t * t, double time)
{
    t->time = time;
}

void strategy_add_fm(strategy_t * t, fm_t * elem)
{
    tabular_fm_add_fm(t->tab_fm, elem);
}

void strategy_add_fm_side(strategy_t * t, fm_t * elem, int side)
{
    tabular_fm_add_fm(t->tab_fm, elem);
    if (t->side == NULL) {
	t->len_side = 1;
	t->side = malloc(sizeof(int) * (t->len_side));
	ASSERT_ALWAYS(t->side != NULL);
    } else {
	t->len_side++;
	t->side = realloc(t->side, sizeof(int) * (t->len_side));
    }
    t->side[t->len_side - 1] = side;
}

strategy_t *strategy_copy(strategy_t * t)
{
    strategy_t *elem = strategy_create();
    tabular_fm_concat(elem->tab_fm, t->tab_fm);
    elem->proba = t->proba;
    elem->time = t->time;
    //side
    if (t->side != NULL) {
	elem->len_side = t->len_side;
	elem->side = malloc(sizeof(int) * (elem->len_side));
	for (int i = 0; i < elem->len_side; i++)
	    elem->side[i] = t->side[i];
    }
    return elem;
}

void strategy_print_file(strategy_t * t, FILE * output_file)
{
    tabular_fm_t *tmp = t->tab_fm;
    //test if the varaible side is used!
    int is_alloced_side = false;
    if (t->side != NULL && t->len_side == tmp->index)
	is_alloced_side = true;

    for (int i = 0; i < tmp->index; i++) {
	fm_t *fm = tmp->tab[i];
	for (int j = 0; j < fm->len_method; j++)
	    fprintf(output_file, "%lu\t", fm->method[j]);
	if (is_alloced_side)
	    fprintf(output_file, "side = %d", t->side[i]);
	fprintf(output_file, "\n");
    }
    fprintf(output_file, "\t Probability: %1.10lf\n", t->proba);
    fprintf(output_file, "\t Time: %lf\n", t->time);
}

void strategy_print(strategy_t * t)
{
    strategy_print_file(t, stdout);
}

/* //test strategy */
/* strategy_t* t = strategy_create (); */
/* fm_t* elem = fm_create (); */
/* unsigned long value[4] = {1,2,3,4}; */
/* double p = 0.5; */
/* double a = 2; */
/* fm_set_method (elem, value,4); */
/* fm_set_proba (elem, &p, 1); */
/* fm_set_time (elem, &a, 1); */

/* strategy_add_fm (t, elem); */
/* p = 2.5; */
/* a = 21; */
/* fm_set_method (elem, value,4); */
/* fm_set_proba (elem, &p, 1); */
/* fm_set_time (elem, &a, 1); */
/* strategy_add_fm (t, elem); */
/* unsigned long value2[6] = {6,5,4,3,2,1}; */
/* p = 4.5; */
/* a = 221; */
/* fm_set_method (elem, value2,6); */
/* fm_set_proba (elem, &p, 1); */
/* fm_set_time (elem, &a, 1); */

/* strategy_add_fm (t, elem); */
/* strategy_print (t); */

/* strategy_free (t); */
/* fm_free (elem); */
/* exit (1); */
