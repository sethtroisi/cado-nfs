#include "cado.h"
#include <stdio.h>
#include <stdlib.h>

#include "portability.h"
#include "utils.h"
#include "filter_utils.h"

/* Print the relation 'rel' in matrix format, i.e., a line of the form:
 * a,b:t_1,t_2,...,t_k

 * a (signed decimal) is a
 * b (nonnegative decimal) is b
 * t_1 ... t_k (hexadecimal) are the indices of the ideals (starting at 0)
 */
//TODO try with buf_rel_t instead of buf_rel_t *
inline void print_relation(FILE * file, buf_rel_t * rel)
{
    char buf[1 << 12], *p;
    char *op;
    size_t t;
    unsigned int i, j;

    p = d64toa16(buf, rel->a);
    *p++ = ',';
    p = u64toa16(p, rel->b);
    *p++ = ':';

    for (i = 0; i < rel->nb; i++) {
	op = p;
	p = u64toa16(p, (uint64_t) rel->primes[i].h);
	*p++ = ',';
	t = p - op;
	for (j = (unsigned int) ((rel->primes[i].e) - 1); j--; p += t)
	    memcpy(p, op, t);
    }

    *(--p) = '\n';
    p[1] = 0;
    fputs(buf, file);
}

 /* Write relation j from buffer to rel_compact
    We put in rel_compact only primes such that their index h is greater or
    equal to min_index
    Return the number of new primes
  */
//TODO try with buf_rel_t instead of buf_rel_t *
inline index_t
insert_rel_in_table_no_e(buf_rel_t * my_br, index_t min_index,
			 uint8_t no_storage, index_t ** rel_compact,
			 weight_t * ideals_weight)
{
    index_t nprimes = 0;
    unsigned int i, itmp;
    index_t *my_tmp;
    index_t h;

    itmp = 0;
    my_tmp = no_storage ? NULL : my_malloc(my_br->nb_above_min_index);

    for (i = 0; i < my_br->nb; i++) {
	h = my_br->primes[i].h;
	if (ideals_weight[h] == 0) {
	    ideals_weight[h] = 1;
	    nprimes++;
	} else if (ideals_weight[h] != UMAX(weight_t))
	    ideals_weight[h]++;

	if (!no_storage && h >= min_index)
	    my_tmp[itmp++] = h;
    }

    if (!no_storage) {
	my_tmp[itmp] = UMAX(*my_tmp);	/* sentinel */
	rel_compact[my_br->num] = my_tmp;
    }

    return nprimes;
}

inline index_t
insert_rel_in_table_with_e(buf_rel_t * my_br, index_t min_index,
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
	my_tmp =
	    (ideal_merge_t *) malloc(my_br->nb_above_min_index *
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

	if (!no_storage && h >= min_index) {
	    my_tmp[itmp].id = h;
	    my_tmp[itmp].e = my_br->primes[i].e;
	    itmp++;
	}
    }

    if (!no_storage) {
	my_tmp[itmp].id = UMAX(my_tmp[itmp].id);	/* sentinel */
	rel_compact[my_br->num] = my_tmp;
    }

    return nprimes;
}

/* if duplicate is_dup = 1, else is_dup = 0 
 * return i for sanity check
 */
inline uint32_t
insert_relation_in_dup_hashtable(uint32_t * H, unsigned long K,
				 buf_rel_t * rel, double *cost,
				 unsigned int *is_dup)
{
    uint64_t h;
    uint32_t i, j;

    h = CA_DUP2 * (uint64_t) rel->a + CB_DUP2 * rel->b;
    i = h % K;
    j = (uint32_t) (h >> 32);
    /* Note: in the case where K > 2^32, i and j share some bits.
     * The high bits of i are in j. These bits correspond therefore to far away
     * positions in the tables, and keeping them in j can only help.
     * FIXME: TODO: that's wrong!!! it would be better do take i from high
     bits instead!
     */
    while (H[i] != 0 && H[i] != j) {
	i++;
	if (UNLIKELY(i == K))
	    i = 0;
	(*cost)++;
    }

    if (H[i] == j)
	*is_dup = 1;
    else {
	H[i] = j;
	*is_dup = 0;
    }

    return i;
}


inline void compute_index_rel(renumber_t tab, buf_rel_t * rel)
{
    unsigned int i;
    p_r_values_t r;
    prime_t *pr = rel->primes;
    int side;

    for (i = 0; i < rel->nb; i++) {
	/* HACK a relation is on the form a,b:side0:side1
	 * we put the side in primes[i].h
	 */
	side = (int) rel->primes[i].h;
	if (side != tab->rat) {
#ifndef FOR_FFS
	    r = (p_r_values_t) findroot(rel->a, rel->b, pr[i].p);
#else
	    r = (p_r_values_t) findroot_ffs(rel->a, rel->b, pr[i].p);
#endif
	    pr[i].h = renumber_get_index_from_p_r(tab, pr[i].p, r, side);
	} else
	    pr[i].h = renumber_get_index_from_p_r(tab, pr[i].p, 0, side);
    }
}
