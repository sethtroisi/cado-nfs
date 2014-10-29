#include "cado.h"
#include "portability.h"
#include "utils.h"

#include "decomp.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

decomp_t *decomp_create(int len, int *tab, double nb_elem)
{
    decomp_t *t = malloc(sizeof(*t));
    ASSERT_ALWAYS(t != NULL);
    t->len = len;
    t->tab = malloc(t->len * sizeof(int));
    ASSERT_ALWAYS(t->tab != NULL);
    for (int i = 0; i < t->len; i++)
	t->tab[i] = tab[i];
    t->nb_elem = nb_elem;
    return t;
}

void decomp_free(decomp_t * t)
{
    if (t != NULL)
	{
	    free(t->tab);
	    free(t);
	}
}

double decomp_get_nb_elem(decomp_t * t)
{
    return t->nb_elem;
}

int *decomp_get_decomp(decomp_t * t)
{
    return t->tab;
}

int decomp_get_len(decomp_t * t)
{
    return t->len;
}

void decomp_set_decomp(decomp_t * t, int *tab, int len)
{
    t->len = len;
    for (int i = 0; i < t->len; i++)
	t->tab[i] = tab[i];
}

void decomp_set_nb_elem(decomp_t * t, double nb_elem)
{
    t->nb_elem = nb_elem;
}

void decomp_print_file(decomp_t * t, FILE * output_file)
{
    if (t == NULL || t->len == 0)
	return;
    fprintf(output_file, "[ ");
    for (int l = 0; l < t->len - 1; l++)
	fprintf(output_file, "%d, ", t->tab[l]);
    fprintf(output_file, "%d ] : %f\n", t->tab[t->len - 1], t->nb_elem);
}

void decomp_print(decomp_t * t)
{
    decomp_print_file(t, stdout);
}
