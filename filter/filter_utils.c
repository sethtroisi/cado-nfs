#include "cado.h"
#include <stdio.h>
#include <stdlib.h>

#include "portability.h"
#include "utils.h"
#include "filter_utils.h"

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

inline index_t
insert_rel_in_table_no_e(earlyparsed_relation_ptr my_br, index_t min_index,
			 index_t ** rel_compact,
			 weight_t * ideals_weight)
{
    index_t nprimes = 0;
    unsigned int itmp;
    index_t *my_tmp;
    index_t h;

    itmp = 0;

    unsigned int nb_above_min_index = earlyparsed_relation_nb_above_min_index(my_br, min_index);

    my_tmp = my_malloc(1 + nb_above_min_index);

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

inline index_t
insert_rel_in_table_with_e(earlyparsed_relation_ptr my_br, index_t min_index,
			   uint8_t no_storage, ideal_merge_t ** rel_compact,
			   weight_t * ideals_weight)
{
    index_t nprimes = 0;
    unsigned int i, itmp;
    ideal_merge_t *my_tmp;
    index_t h;

    itmp = 0;
    if (no_storage)
	my_tmp = NULL;
    else {
	//FIXME for now can't use my malloc, because it expected an index_t table
        unsigned int nb_above_min_index = earlyparsed_relation_nb_above_min_index(my_br, min_index);

	my_tmp =
	    (ideal_merge_t *) malloc((1 + nb_above_min_index) *
				     sizeof(ideal_merge_t));
	ASSERT_ALWAYS(my_tmp != NULL);
    }

    for (i = 0; i < my_br->nb; i++) {
	h = my_br->primes[i].h;
	if (ideals_weight[h] == 0) {
	    ideals_weight[h] = 1;
	    nprimes++;
	} else if (ideals_weight[h] != UMAX(weight_t))
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
