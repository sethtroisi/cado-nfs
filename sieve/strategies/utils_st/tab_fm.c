#include "cado.h"
#include "portability.h"
#include "utils.h"
#include "tab_fm.h"

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "facul.h"

static const double EPSILON_DBL = 0.000001;

tabular_fm_t *tabular_fm_create(void)
{
    tabular_fm_t *t = malloc(sizeof(*t));
    ASSERT_ALWAYS(t != NULL);

    t->index = 0;
    t->size = 2;

    t->tab = malloc(t->size * sizeof(fm_t *));
    ASSERT_ALWAYS(t->tab != NULL);

    return t;
}

void tabular_fm_free(tabular_fm_t * t)
{
    for (int i = 0; i < t->index; i++)	//size
	fm_free(t->tab[i]);
    free(t->tab);
    free(t);
}

void tabular_fm_realloc(tabular_fm_t * t)
{
    t->tab = realloc(t->tab, t->size * 2 * (sizeof(fm_t *)));
    ASSERT_ALWAYS(t->tab != NULL);
    t->size *= 2;
}

int tabular_fm_get_index(tabular_fm_t * t)
{
    return t->index;
}

void tabular_fm_add_fm(tabular_fm_t * t, fm_t * fm)
{
    tabular_fm_set_fm_index(t, fm, t->index);
}

void
tabular_fm_add(tabular_fm_t * t, unsigned long *method, int len_method,
	       double *proba, int len_proba, double *time, int len_time,
	       int len_p_min)
{
    tabular_fm_set_index(t, method, len_method, proba, len_proba, time,
			 len_time, len_p_min, t->index);
}

void tabular_fm_set_fm_index(tabular_fm_t * t, fm_t * fm, int ind)
{
    tabular_fm_set_index(t, fm->method, fm->len_method, fm->proba,
			 fm->len_proba, fm->time, fm->len_time, fm->len_p_min,
			 ind);
}

void
tabular_fm_set_index(tabular_fm_t * t, unsigned long *method, int len_method,
		     double *proba, int len_proba, double *time, int len_time,
		     int len_p_min, int ind)
{
    if (ind >= t->size)
	tabular_fm_realloc(t);

    if (ind >= t->index) {
	t->tab[ind] = fm_create();
	ASSERT_ALWAYS(t->tab[ind] != NULL);
    }

    fm_set_method(t->tab[ind], method, len_method);
    fm_set_proba(t->tab[ind], proba, len_proba, len_p_min);
    fm_set_time(t->tab[ind], time, len_time);
    if (ind >= t->index)
	t->index++;
}

fm_t *tabular_fm_get_fm(tabular_fm_t * t, int index)
{
    ASSERT_ALWAYS(index <= t->index);
    return t->tab[index];
}

void tabular_fm_concat(tabular_fm_t * t1, tabular_fm_t * t2)
{
    int len = t2->index;
    for (int i = 0; i < len; i++)
	tabular_fm_add_fm(t1, t2->tab[i]);
}

void tabular_fm_put_zero(tabular_fm_t * t, int index)
{
    ASSERT_ALWAYS(index <= t->index);
    fm_put_zero(t->tab[index]);
}

bool tabular_fm_is_zero(tabular_fm_t * t, int index)
{
    ASSERT_ALWAYS(index <= t->index);
    return fm_is_zero(t->tab[index]);
}

tabular_fm_t *extract_fm_method(tabular_fm_t * t, int method, int curve)
{
    tabular_fm_t *res = tabular_fm_create();
    int len = t->index;
    for (int i = 0; i < len; i++) {
	fm_t *el = t->tab[i];
	if ((int)el->method[0] == method) {
	    if (method == EC_METHOD) {
		if ((int)el->method[1] == curve)
		    tabular_fm_add_fm(res, el);
	    } else
		tabular_fm_add_fm(res, el);
	}
    }
    return res;
}

int tab_fm_is_exist(tabular_fm_t * c, fm_t * f)
{
    int len = c->index;
    for (int i = 0; i < len; i++) {
	if (fm_is_equal(c->tab[i], f))
	    return true;
    }
    return false;
}

/************************************************************************/
/*           PRINT AND SCAN OUR FILES OF FACTORING METHODS              */
/************************************************************************/

int tabular_fm_print(tabular_fm_t * t)
{
    return tabular_fm_fprint(stdout, t);
}

int tabular_fm_fprint(FILE * file, tabular_fm_t * t)
{
    int len = t->index;
    for (int i = 0; i < len; i++) {
	fm_t *elem = tabular_fm_get_fm(t, i);
	if (fm_fprint(file, elem) < 0)
	    return -1;
    }
    return 0;
}

static int is_elem(const char c)
{
    return (c >= 48 && c <= 57) || c == '|';
}

static void next_elem(FILE * file, char *current_char)
{
    //end the current element
    while (*current_char != EOF && is_elem(*current_char))
	*current_char = fgetc(file);

    //find the next element
    while (*current_char != EOF && !is_elem(*current_char)) {
	*current_char = fgetc(file);
    }
    fseek(file, -1, SEEK_CUR);
}

static fm_t *sub_routine_fm_fscanf(FILE * file, char *current_char)
{
    fm_t *fm = fm_create();
    //One knows the number of times and methods but not for the probabilities.
    int len_method = 4;
    int len_proba = 200;
    int len_time = 4;

    unsigned long method[len_method];
    double proba[len_proba];
    double time[len_time];

    int ind = 0;
    while (ind < len_method && *current_char != '|') {
	fscanf(file, "%lu", &method[ind++]);
	next_elem(file, current_char);
    }
    next_elem(file, current_char);

    fm_set_method(fm, method, ind);

    fscanf(file, "%d", &fm->len_p_min);
    next_elem(file, current_char);
    ind = 0;
    while (ind < len_proba && *current_char != '|') {
	fscanf(file, "%lf", &proba[ind++]);
	next_elem(file, current_char);
    }
    next_elem(file, current_char);

    fm_set_proba(fm, proba, ind, fm->len_p_min);

    ind = 0;
    while (ind < len_time && *current_char != '|') {
	fscanf(file, "%lf", &time[ind++]);
	next_elem(file, current_char);
    }
    next_elem(file, current_char);

    fm_set_time(fm, time, ind);

    return fm;
}

int fm_fscan(FILE * file, tabular_fm_t * res)
{
    if (file == NULL)
	return -1;

    char current_char = fgetc(file);
    fseek(file, -1, SEEK_CUR);
    while (current_char != EOF) {
	fm_t *fm = sub_routine_fm_fscanf(file, &current_char);
	tabular_fm_add_fm(res, fm);
	fm_free(fm);
    }

    return 0;
}

/************************************************************************/
/*                      SORT_TAB_FM                                     */
/************************************************************************/

//return a positive value if el1 is greater than el2 and a negative value otherwise.
static int fm_cmp(fm_t * el1, fm_t * el2)
{
    //compare the mean value of time/succes
    //we will suppose that len_proba equal to len_time

    int len1 = el1->len_proba;
    int len2 = el2->len_proba;
    int len_min = (len1 < len2) ? len1 : len2;
    //  double compromise1, compromise2;
    int tot1 = 0, tot2 = 0;
    int non_zero;
    for (int i = 0; i < len_min; i++) {
	non_zero = true;
	if (el1->proba[i] < EPSILON_DBL) {
	    tot2++;
	    non_zero = false;
	}
	if (el2->proba[i] < EPSILON_DBL) {
	    tot1++;
	    non_zero = false;
	}
	if (non_zero) {
	    /* compromise1 = el1->time[i]/el1->proba[i]; */
	    /* compromise2 = el2->time[i]/el2->proba[i]; */
	    //double diff_compromise = (compromise1 - compromise2);
	    double diff_compromise = el1->proba[i] - el2->proba[i];
	    if (diff_compromise < EPSILON_DBL)
		tot2++;
	    else if (diff_compromise > EPSILON_DBL)
		tot1++;
	}
    }
    return tot1 - tot2;
}

void fm_swap(tabular_fm_t * t, int index1, int index2)
{
    fm_t *c = t->tab[index1];
    t->tab[index1] = t->tab[index2];
    t->tab[index2] = c;
}

static void tabular_fm_sort_rec(tabular_fm_t * t, int begin, int end)
{
    int index_max = begin;
    for (int i = begin; i < end; i++) {
	if (fm_cmp(t->tab[i], t->tab[index_max]) > 0) {
	    index_max = i;
	}
    }
    fm_swap(t, end - 1, index_max);
}

void tabular_fm_sort(tabular_fm_t * t)
{
    int max = t->index;
    while (max > 0) {
	tabular_fm_sort_rec(t, 0, max);
	max--;
    }
}

/* /\* gcc -c -std=c99 -Wall fm.c -o fm.o; gcc -std=c99 -Wall -g -c tab_fm.c -lm  -o tab_fm.o; gcc -g tab_fm.o fm.o -lm -o tab_fm*\/ */
/* int */
/* main () */
/* { */
/*   fm_t* t = fm_create(); */
/*   unsigned long elem[4]; */
/*   for (int i=0; i < 4; i++) */
/*     { */
/*       elem[i] = i; */
/*     } */
/*   fm_set_method (t, elem, 4); */
/*   double* elem2 = malloc (sizeof(double)*10); */

/*   for (int i=0; i < 10; i++) */
/*     { */
/*       elem2[i] = i*10; */
/*     } */
/*   fm_set_proba (t, elem2, 10, 0); */

/*   double elem3[10]; */
/*   for (int i=0; i < 10; i++) */
/*     { */
/*       elem3[i] = i*100; */
/*     } */
/*   fm_set_time (t, elem3, 10, 0); */

/*   fm_print (t); */

/*   //test tab_fm_t */
/*   tabular_fm_t* test = tabular_fm_create(); */
/*   tabular_fm_add_fm (test, t); */
/*   tabular_fm_print (test); */
/*   tabular_fm_add (test, elem, 4, elem2, 10, elem3, 10); */
/*   tabular_fm_print (test); */
/*   tabular_fm_add (test, elem, 4, elem2, 10, elem3, 10); */
/*   tabular_fm_print (test); */
/*   unsigned long a = 2; */
/*   double p = 0.5; */
/*   double t1 = 4.2; */
/*   tabular_fm_add (test, &a, 1, &p, 1, &t1, 1); */
/*   tabular_fm_print (test); */
/*   fm_set_proba(test->tab[3], elem2, 10, 0); */
/*   tabular_fm_print (test); */

/*   printf("test2\n"); */
/*   tabular_fm_t* test2 = tabular_fm_create(); */
/*   tabular_fm_add (test2, &a, 1, &p, 1, &t1, 1); */
/*   tabular_fm_print (test2); */
/*   printf("concat\n"); */
/*   tabular_fm_concat (test, test2); */
/*   tabular_fm_print (test); */

/*   //free */
/*   free (elem2); */
/*   fm_free (t); */
/*   tabular_fm_free (test); */
/*   tabular_fm_free (test2); */
/* } */
