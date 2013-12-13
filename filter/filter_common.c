#include "cado.h"
#include <stdio.h>
#include <stdlib.h>

#include "portability.h"
#include "utils.h"
#include "filter_common.h"

static inline unsigned int earlyparsed_relation_nb_above_min_index(earlyparsed_relation_srcptr rel, index_t min_index)
{
    unsigned int nb_above_min_index = 0;
    for (unsigned int i = 0; i < rel->nb; i++) {
        nb_above_min_index += rel->primes[i].h >= min_index;
    }
    return nb_above_min_index;
}

/* Write relation j from buffer to rel_compact
 * We put in rel_compact only primes such that their index h is greater or
 * equal to min_index
 * Return the number of new primes
 */

inline unsigned int
insert_rel_in_table_no_e(earlyparsed_relation_ptr my_br, index_t min_index,
			 index_t ** rel_compact,
			 weight_t * ideals_weight)
{
    unsigned int nprimes = 0;
    unsigned int itmp;
    index_t *my_tmp;
    index_t h;

    itmp = 0;

    unsigned int nb_above_min_index = earlyparsed_relation_nb_above_min_index(my_br, min_index);

    my_tmp = index_my_malloc(1 + nb_above_min_index);

    for (unsigned int i = 0; i < my_br->nb; i++) {
	h = my_br->primes[i].h;
	if (ideals_weight[h] == 0) {
	    ideals_weight[h] = 1;
	    nprimes++;
	} else if (ideals_weight[h] != UMAX(weight_t))
	    ideals_weight[h]++;

	if (h >= min_index)
	    my_tmp[itmp++] = h;
    }

    my_tmp[itmp] = UMAX(*my_tmp);	/* sentinel */
    rel_compact[my_br->num] = my_tmp;

    return nprimes;
}

inline unsigned int
insert_rel_in_table_with_e(earlyparsed_relation_ptr my_br, index_t min_index,
                           ideal_merge_t **rel_compact, int32_t *ideals_weight)
{
    unsigned int nprimes = 0;
    unsigned int i, itmp, nb_above_min;
    ideal_merge_t *my_tmp;
    index_t h;

    itmp = 0;
    nb_above_min = earlyparsed_relation_nb_above_min_index(my_br, min_index);
    my_tmp = idealmerge_my_malloc(1 + nb_above_min);

    for (i = 0; i < my_br->nb; i++) {
	h = my_br->primes[i].h;
	if (ideals_weight[h] == 0) {
	    ideals_weight[h] = 1;
	    nprimes++;
	} else if (ideals_weight[h] != SMAX(int32_t))
	    ideals_weight[h]++;

	if (my_tmp && h >= min_index) {
	    my_tmp[itmp].id = h;
	    my_tmp[itmp].e = my_br->primes[i].e;
	    itmp++;
	}
    }

    if (my_tmp) {
	my_tmp[itmp].id = UMAX(my_tmp[itmp].id);	/* sentinel */
	rel_compact[my_br->num] = my_tmp;
    }

    return nprimes;
}

int
cmp_index (const void *p, const void *q)
{
  index_t x = *((index_t *)p);
  index_t y = *((index_t *)q);
  return (x <= y ? -1 : 1);
}

int
cmp_ideal_merge (const void *p, const void *q)
{
  ideal_merge_t x = *((ideal_merge_t *)p);
  ideal_merge_t y = *((ideal_merge_t *)q);
  return (x.id <= y.id ? -1 : 1);
}

int
cmp_int (const void *p, const void *q)
{
  int x = *((int *)p);
  int y = *((int *)q);
  return (x <= y ? -1 : 1);
}

/* We also compare x[1] and y[1] to make the code deterministic
   since in case x[0] = y[0], qsort() may give different results on
   different machines */
int
cmp_int2 (const void *p, const void *q)
{
  int *x = (int*) p;
  int *y = (int*) q;

  if (x[0] < y[0])
    return -1;
  else if (x[0] > y[0])
    return 1;
  else
    return (x[1] < y[1]) ? 1 : -1;
}
