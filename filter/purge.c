/* purge --- remove singletons
 * 
 * Copyright 2008, 2009, 2010, 2011, 2012 Alain Filbois, Francois Morain,
 *                                        Paul Zimmermann
 * 
 * This file is part of CADO-NFS.
 * 
 * CADO-NFS is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * CADO-NFS is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with CADO-NFS; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
*/

/* References:
 * On the Number Field Sieve Integer Factorisation Algorithm,
 * Stefania Cavallar, PhD Thesis, University of Leiden, 2002.
 */

/*
 * This program works in two passes over the relation files:
 * - the first pass loads in memory only rational primes >= minpr and algebraic
 *   ideals >= minpa, but stores all ideals in the hash table, and keeps a count
 *   of the number of non-stored ideals for each relation.
 *   By default minpr and minpa are taken as rlim and alim respectively.
 *   Then simultaneously singleton removal is performed, and heavy relations
 *   are discarded, until the final excess is 'keep'.
 * - the second pass goes through the relations again, and dumps the remaining
 *   ones in the format needed by 'merge'.

 * This program uses the following data structures:
 * rel_used[i]    - non-zero iff relation i is kept (so far)
 * rel_compact[i] - list of 'h' indices in H table of considered (p,r) for row i
 *                  (terminated by a -1 sentinel)
 * ideals_weight [h] - number of occurrences of h in current relations

 * Exit value:
 * - 0 if enough relations
 * - 1 if an error occurred (then we should abort the factorization)
 * - 2 if not enough relations
 */

/*
 * index_t : (32 + 32 * need64) bits. Signed.
 * p_r_values_t : 32 if H.hm < 2^32-1, otherwise 64 bits.
 * Note : sizeof(p_r_values_t) >= sizeof(index_t)
*/

#include "cado.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>		/* for _O_BINARY */
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>
#ifdef HAVE_LIBGEN_H
#include <libgen.h>
#endif

#include "mod_ul.c"
#include "portability.h"
#include "utils.h"
#include "typedefs.h"

#include "filter_utils.h"

//#define STAT
//#define STAT_VALUES_COEFF //STAT must be defined. Interesting only DL

//#define USE_CAVALLAR_WEIGHT_FUNCTION

#define STR(s) XSTR(s)
#define XSTR(s) #s
#define DEFAULT_NPASS 50
#define DEFAULT_KEEP 160
#define DEFAULT_REQUIRED_EXCESS 0.1
#define DEFAULT_NPT 4

/* Main variables */

/* This one is passed to all functions, so it's morally a global */
struct purge_data_s {
    /* for minimal changes, don't put these here _yet_ */
    // uint64_t min_index;
    // index_t **rel_compact;	/* see main documentation */
    // weight_t *ideals_weight;
    info_mat_t info;
    /* fd[0]: printed relations */
    /* fd[1]: deleted relations */
    FILE * fd[2];
};
typedef struct purge_data_s purge_data[1];
typedef struct purge_data_s * purge_data_ptr;
typedef const struct purge_data_s * purge_data_srcptr;

uint64_t min_index;
index_t **rel_compact;	/* see main documentation */
weight_t *ideals_weight;


static bit_vector rel_used, Tbv;
static index_t *sum2_index = NULL;	/*sum of rows index for primes of weight 2 */

static uint64_t nrelmax = 0, nprimemax = 0;
static int64_t keep = DEFAULT_KEEP;	/* maximun final excess */
static unsigned int npass = DEFAULT_NPASS;
static double required_excess = DEFAULT_REQUIRED_EXCESS;
static unsigned int npt = DEFAULT_NPT;

static float w_ccc;

#ifdef STAT
uint64_t __stat_weight;
#ifdef STAT_VALUES_COEFF
#define STAT_VALUES_COEFF_LEN 10
uint64_t __stat_nbcoeffofvalue[STAT_VALUES_COEFF_LEN + 1];
#endif
#endif

/*****************************************************************************/


/* Delete a relation: set rel_used[i] to 0, update the count of primes
 * in that relation.
 * Warning: we only update the count of primes that we consider, i.e.,
 * primes with index >= min_index.
 */
static unsigned int delete_relation (uint64_t i)
{
    index_t *tab;
    unsigned int nremoveprimes = 0;
    weight_t *o;

    for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++) {
	o = &(ideals_weight[*tab]);
	ASSERT(*o);
	if (*o < UMAX(*o) && !(--(*o)))
	    nremoveprimes++;
    }

    /* We do not free rel_compact[i] as it is freed with my_malloc_free_all */
    bit_vector_clearbit(rel_used, (size_t) i);

    return nremoveprimes;
}



/*****************************************************************************/
/* Code for clique removal.
 * A clique is a connected components of the relation R, where R(i1,i2) iff
 * i1 and i2 share a prime of weight 2.
 * We remove the heaviest cliques.
 * Each ideal h contributes to the weight of the cliques depending on its
 * weight (see weight_function_clique).
 */

typedef struct {
    float w;
    index_t i;
} comp_t;

static int compare(const void *v1, const void *v2)
{
    comp_t *w1 = (comp_t *) v1;
    comp_t *w2 = (comp_t *) v2;

    return (w1->w >= w2->w) ? -1 : 1;
}

float weight_function_clique(weight_t w)
{
#ifdef USE_CAVALLAR_WEIGHT_FUNCTION
    if (w >= 3)
	return ldexpf(1, -(w - 1));
    else if (w == 2)
	return 0.25;
    else
	return 0.0;
#else
    if (w >= 3)
	return powf(0.8, (float) (w - 2.0));
    else if (w == 2)
	return 0.125;
    else
	return 0.0;
#endif
}

/* Compute connected component of row i for the relation R(i1,i2) if rows
 * i1 and i2 share a prime of weight 2.
 * Return number of rows of components, and put in w the total weight. */
static index_t compute_connected_component(index_t i)
{
    index_t *myrel_compact = rel_compact[i], h, k, n = 1;

    bit_vector_setbit(Tbv, (size_t) i);	/* mark row as visited */
    while ((h = *myrel_compact++) != UMAX(h)) {
	if (ideals_weight[h] == 2) {
	    k = sum2_index[h] - i;	/* other row where prime of index h appears */
	    if (!bit_vector_getbit(Tbv, (size_t) k))	/* row k not visited yet */
		n += compute_connected_component(k);
	}
	/* we use the multiplier 5 here, so that the average weight (assumed to
	   be 1) is in the middle of the Count[10] array */
	w_ccc += 5.0 * weight_function_clique(ideals_weight[h]);
    }
    return n;
}

/* Delete connected component of row i, assuming the bit-vector is set.
 * Warning: we might have some H->hashcount[h] = 3, which is decreased
 * to 2, but we don't want to treat that case. Thus we check in addition
 * that sum2_index[h] <> 0, which only occurs when H->hashcount[h] = 2
 * initially.
 */
static index_t delete_connected_component(index_t i, uint64_t * nprimes)
{
    index_t *myrel_compact = rel_compact[i], h, k, w = 1;

    bit_vector_clearbit(Tbv, (size_t) i);	/* mark row as visited */
    /* bit i of rel_used is cleared in delete_relation below */
    while ((h = *myrel_compact++) != UMAX(h)) {
	if (ideals_weight[h] == 2 && sum2_index[h]) {	/* first row that contains ideal of index h */
	    k = sum2_index[h] - i;	/* other row where prime of index h appears */
	    if (bit_vector_getbit(Tbv, (size_t) k) == 1)	/* row k was not visited yet */
		w += delete_connected_component(k, nprimes);
	}
    }
    *nprimes -= delete_relation(i);
    return w;
}

static void
cliques_removal(int64_t target_excess, uint64_t * nrels, uint64_t * nprimes)
{
    int64_t excess = (((int64_t) * nrels) - *nprimes);
    int64_t chunk;
    uint64_t i;
    comp_t *tmp = NULL;		/* (weight, index) */
    index_t *myrelcompact, h;
    double N = 0.0;		/* number of rows */
    unsigned int wceil, j, ltmp = 0, alloctmp = 0xff;
#define MAX_WEIGHT 10
    unsigned int Count[MAX_WEIGHT], *pc;	/* Count[w] is the # of components of weight >= w */

    if (excess <= keep || excess <= target_excess)
	return;

    chunk = excess - target_excess;

    /* first collect sums for primes with weight 2, and compute total weight */
    memset(sum2_index, 0, nprimemax * sizeof(index_t));
    for (i = 0; i < nrelmax; i++)
	if (bit_vector_getbit(rel_used, (size_t) i)) {
	    for (myrelcompact = rel_compact[i];
		 (h = *myrelcompact++) != UMAX(h);)
		if (ideals_weight[h] == 2)
		    sum2_index[h] += i;
	    N += 1.0;
	}

    ASSERT_ALWAYS(N == (double) *nrels);

    /* now initialize bit table for relations used */
    bit_vector_neg(Tbv, rel_used);
    memset(Count, 0, sizeof(unsigned int) * MAX_WEIGHT);
    tmp = (comp_t *) malloc(alloctmp * sizeof(comp_t));
    ASSERT_ALWAYS(tmp != NULL);

    for (i = 0; i < nrelmax; i++) {
	if (!bit_vector_getbit(Tbv, (size_t) i)) {
	    w_ccc = 0.0;
	    h = compute_connected_component(i);
	    wceil = (unsigned int) w_ccc + 1;
	    /* don't consider weight-0 components */
	    if ((w_ccc > 0.0)
		&& ((wceil >= MAX_WEIGHT) || (Count[wceil] < chunk))) {
		if (ltmp >= alloctmp) {
		    alloctmp = 2 * alloctmp + 1;
		    tmp = (comp_t *) realloc(tmp, alloctmp * sizeof(comp_t));
		}
		tmp[ltmp++] = (comp_t) {
		w_ccc, i};
		if (wceil > MAX_WEIGHT)
		    wceil = MAX_WEIGHT;
		for (pc = &(Count[wceil]); pc-- != Count; (*pc)++);
	    }
	}
    }
    qsort(tmp, ltmp, sizeof(comp_t), compare);

    for (j = 0; j < ltmp && *nrels > target_excess + *nprimes; j++)
	*nrels -= delete_connected_component(tmp[j].i, nprimes);

    fprintf(stderr, "    deleted %u heavier connected components at %2.2lf\n",
	    j, seconds());

#if DEBUG >= 1
    fprintf(stderr, "    DEBUG: ltmp=%u chunk=%u target=%u\n",
	    ltmp, chunk, target_excess);
#endif

    free(tmp);
}

/*****************************************************************************/
/* Code for singletons removal.
 * Exist in multithread if __sync_sub_and_fetch exists.
 */

#ifndef HAVE_SYNC_FETCH
static void onepass_singleton_removal(uint64_t * nrels, uint64_t * nprimes)
{
    index_t *tab;
    uint64_t i, nremoverels = 0, nremoveprimes = 0;

    for (i = 0; i < nrelmax; i++) {
	if (bit_vector_getbit(rel_used, (size_t) i)) {
	    for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++) {
		if (ideals_weight[*tab] == 1) {
		    nremoveprimes += delete_relation(i);
		    nremoverels++;
		    break;
		}
	    }
	}
    }

    *nrels -= nremoverels;
    *nprimes -= nremoveprimes;
}

#else				/* ifndef HAVE_SYNC_FETCH */

typedef struct {
    unsigned int nb;
    pthread_t mt;
    uint64_t begin, end, sup_nrel, sup_npri;
} ti_t;
static ti_t *ti;

/* Hightest criticality for performance. I inline all myself. */
static void onepass_thread_singleton_removal(ti_t * mti)
{
    index_t *tab;
    uint64_t i;
    weight_t *o;
    bv_t j;

    mti->sup_nrel = mti->sup_npri = 0;
    for (i = mti->begin; i < mti->end; i++) {
	j = (((bv_t) 1) << (i & (BV_BITS - 1)));

	if (rel_used->p[i >> LN2_BV_BITS] & j)
	    for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++)
		if (ideals_weight[*tab] == 1) {
		    for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++) {
			o = &(ideals_weight[*tab]);
			ASSERT(*o);
			if (*o < UMAX(*o) && !__sync_sub_and_fetch(o, 1))
			    (mti->sup_npri)++;
		    }
		    /* rel_compact[i] is not freed , it is freed with my_malloc_free_all */
		    rel_used->p[i >> LN2_BV_BITS] &= ~j;
		    (mti->sup_nrel)++;
		    break;
		}
    }
    pthread_exit(NULL);
}

static void
onepass_singleton_parallel_removal(unsigned int nb_thread, uint64_t * nrels,
				   uint64_t * nprimes)
{
    pthread_attr_t attr;
    uint64_t pas, incr;
    unsigned int i;
    int err;

    ti = (ti_t *) malloc(nb_thread * sizeof(ti_t));
    ASSERT_ALWAYS(ti != NULL);
    ti[0].begin = 0;
    pas = (nrelmax / nb_thread) & ((uint64_t) ~ (BV_BITS - 1));
    incr = 0;
    for (i = 0, incr = 0; i < nb_thread - 1;) {
	incr += pas;
	ti[i].nb = i;
	ti[i++].end = incr;
	ti[i].begin = incr;
    }
    ti[i].nb = i;
    ti[i].end = nrelmax;
    pthread_attr_init(&attr);
    pthread_attr_setstacksize(&attr, 1 << 16);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    for (i = 0; i < nb_thread; i++)
	if ((err =
	     pthread_create(&ti[i].mt, &attr,
			    (void *) onepass_thread_singleton_removal,
			    &ti[i]))) {
	    fprintf(stderr,
		    "onepass_singleton_parallel_removal : pthread_create error 1: %d. %s\n",
		    err, strerror(errno));
	    exit(1);
	}
    for (i = 0; i < nb_thread; i++) {
	pthread_join(ti[i].mt, NULL);
	*nrels -= ti[i].sup_nrel;
	*nprimes -= ti[i].sup_npri;
    }
    pthread_attr_destroy(&attr);
    if (ti != NULL)
	free(ti);
    ti = NULL;
}
#endif				/* ifdef HAVE_SYNC_FETCH */

static void
remove_all_singletons(uint64_t * nrels, uint64_t * nprimes, int64_t * excess)
{
    uint64_t oldnrels;
    *excess = (((int64_t) * nrels) - *nprimes);
    fprintf(stderr,
	    "  nrels=%" PRIu64 " nprimes=%" PRIu64 " excess=%" PRId64 "\n",
	    *nrels, *nprimes, *excess);
    do {
	oldnrels = *nrels;
#ifdef HAVE_SYNC_FETCH
	ASSERT(npt);
	onepass_singleton_parallel_removal(npt, nrels, nprimes);
#else
	onepass_singleton_removal(nrels, nprimes);
#endif
	*excess = (((int64_t) * nrels) - *nprimes);
	fprintf(stderr,
		"  new_nrels=%" PRIu64 " new_nprimes=%" PRIu64 " excess=%" PRId64
		"" " at %2.2lf\n", *nrels, *nprimes, *excess, seconds());
    } while (oldnrels != *nrels);
}

static void singletons_and_cliques_removal(uint64_t * nrels, uint64_t * nprimes)
{
    uint64_t oldnrels = 0;
    int64_t oldexcess, excess, target_excess;
    unsigned int count;

    //First step of singletons removal
    remove_all_singletons(nrels, nprimes, &excess);

    if (excess <= 0) {		/* covers case nrel = nprimes = 0 */
	fprintf(stderr, "number of relations <= number of ideals\n");
	exit(2);
    }

    if ((double) excess < required_excess * ((double) *nprimes)) {
	fprintf(stderr,
		"(excess / nprimes) = %.2f < %.2f. See -required_excess "
		"argument.\n", ((double) excess / (double) *nprimes),
		required_excess);
	exit(2);
    }

    int64_t chunk = excess / npass;

    //npass pass of clique removal + singletons removal
    for (count = 0; count < npass && excess > 0; count++) {
	oldnrels = *nrels;
	oldexcess = excess;
	target_excess = excess - chunk;
	if (target_excess < keep)
	    target_excess = keep;
	fprintf(stderr, "Step %u on %u: target excess is %" PRId64 "\n",
		count, npass, target_excess);
	cliques_removal(target_excess, nrels, nprimes);

	remove_all_singletons(nrels, nprimes, &excess);
	fprintf(stderr, "  [each excess row deleted %2.2lf rows]\n",
		(double) (oldnrels - *nrels) / (double) (oldexcess - excess));
    }


    /* May need an extra pass of clique removal + singletons removal if excess is
       still larger than keep. It may happen due to the fact that each clique does
       not make the excess go down by one but can (rarely) left the excess
       unchanged. */
    if (excess > keep) {
	oldnrels = *nrels;
	oldexcess = excess;
	target_excess = excess - chunk;
	target_excess = keep;

	fprintf(stderr, "Step extra: target excess is %" PRId64 "\n",
		target_excess);
	cliques_removal(target_excess, nrels, nprimes);

	remove_all_singletons(nrels, nprimes, &excess);
	fprintf(stderr, "  [each excess row deleted %2.2lf rows]\n",
		(double) (oldnrels - *nrels) / (double) (oldexcess - excess));
    }
}

/*****************************************************************************/
/* I/O functions */

/* Callback functions called by filter_rels */
void *insert_rel_into_table(purge_data_ptr arg, earlyparsed_relation_ptr rel)
{
    ASSERT_ALWAYS(rel->num < nrelmax);

    arg->info.nprimes +=
        insert_rel_in_table_no_e(rel, min_index, 
                rel_compact, ideals_weight);
#ifdef STAT
    /* here we also used to accumulate the number of primes above
     * min_index in arg->info.W */
    arg->info.W += earlyparsed_relation_nb_above_min_index(rel, min_index);
#endif

    return NULL;
}

void *thread_print(purge_data_ptr arg, earlyparsed_relation_ptr rel)
{
    if (bit_vector_getbit(rel_used, rel->num)) {
        arg->info.W += rel->nb;
        fputs(rel->line, arg->fd[0]);
    } else if (arg->fd[1] != NULL) {
        fputs(rel->line, arg->fd[1]);
    }
    return NULL;
}


/*********** utility functions for purge binary ****************/



  /* Build the file list (ugly). It is the concatenation of all
   *  b s p
   * where:
   *    b is the basepath (empty if not given)
   *    s ranges over all subdirs listed in the subdirlist (empty if no
   *    such list)
   *    p ranges over all paths listed in the filelist.
   */
char **filelist_from_file_with_subdirlist(const char *basepath,
					  const char *filelist,
					  const char *subdirlist)
{
    /* count the number of files in the filelist */
    int nfiles = 0;
    int nsubdirs = 0;
    char **fl = filelist_from_file(NULL, filelist, 0);
    for (char **p = fl; *p; p++, nfiles++);

    char **sl = filelist_from_file(basepath, subdirlist, 1);
    for (char **p = sl; *p; p++, nsubdirs++);

    char ** fic = (char **) malloc((nsubdirs * nfiles + 1) * sizeof(char *));
    ASSERT_ALWAYS(fic != NULL);

    char **full = fic;
    for (char **f = fl; *f; f++) {
	for (char **s = sl; *s; s++, full++) {
	    int ret = asprintf(full, "%s/%s", *s, *f);
	    ASSERT_ALWAYS(ret >= 0);
	}
    }
    *full = NULL;
    filelist_clear(fl);
    filelist_clear(sl);
    return fic;
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "filelist", "file containing a list of input files");
  param_list_decl_usage(pl, "subdirlist",
                               "file containing a list of subdirectories");
  param_list_decl_usage(pl, "basepath", "path added to all file in filelist");
  param_list_decl_usage(pl, "out", "outfile for remaining relations");
  param_list_decl_usage(pl, "nrels", "number of initial relations");
  param_list_decl_usage(pl, "nprimes",
                                  "number of prime ideals in renumber table");
  param_list_decl_usage(pl, "minindex", "index of the first considered prime");
  param_list_decl_usage(pl, "keep", "wanted excess at the end of purge "
                                    "(default " STR(DEFAULT_KEEP) ")");
  param_list_decl_usage(pl, "npass", "number of step of clique removal "
                                     "(default " STR(DEFAULT_NPASS) ")");
  param_list_decl_usage(pl, "required_excess", "\% of excess required at the "
                            "end of the 1st singleton removal step (default "
                            STR(DEFAULT_REQUIRED_EXCESS) ")");
  param_list_decl_usage(pl, "outdel", "outfile for deleted relations (for DL)");
  param_list_decl_usage(pl, "npthr", "number of threads used for singletons "
                                     "removal (default " STR(DEFAULT_NPT) ")");
  param_list_decl_usage(pl, "force-posix-threads", "(switch)");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}


/*************************** main ********************************************/

int main(int argc, char **argv)
{
    char * argv0 = argv[0];
    int k;
    param_list pl;
    min_index = UMAX(uint64_t);
    uint64_t nrels, nprimes;
    size_t tot_alloc_bytes = 0;
    char ** input_files;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;		/* Binary open for all files */
#endif

    double wct0 = wct_seconds();

    param_list_init(pl);
    declare_usage(pl);
    argv++,argc--;

    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

    if (argc == 0)
      usage (pl, argv0);

    /* read all command-line parameters */
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
    }
    /* print command-line arguments */
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    /* read command-line parameters */
    param_list_parse_uint64(pl, "nrels", &nrelmax);
    param_list_parse_uint64(pl, "nprimes", &nprimemax);
    param_list_parse_int64(pl, "keep", &keep);

    /* Only look at relations of index above minindex */
    param_list_parse_uint64(pl, "minindex", &min_index);


    param_list_parse_uint(pl, "npthr", &npt);
    param_list_parse_uint(pl, "npass", &npass);
    param_list_parse_double(pl, "required_excess", &required_excess);

    /* These three parameters specify the set of input files, of the form
     * <base path>/<one of the possible subdirs>/<one of the possible
     * file names>
     *
     * possible subdirs are lister in the file passed as subdirlist.
     * Ditto for possible file names.
     *
     * file names need not be basenames, i.e. they may contain directory
     * components. subdirlist and basepath are optional.
     */
    const char *basepath = param_list_lookup_string(pl, "basepath");
    const char *subdirlist = param_list_lookup_string(pl, "subdirlist");
    const char *filelist = param_list_lookup_string(pl, "filelist");
    const char *purgedname = param_list_lookup_string(pl, "out");
    const char *deletedname = param_list_lookup_string(pl, "outdel");

    if (param_list_warn_unused(pl))
    {
      fprintf(stderr, "Error, unused parameters are given\n");
      usage(pl, argv0);
    }

    /* }}} */


    /*{{{ argument checking, and some statistics for things related to
     * the hash table. It needs several static parameters on the command
     * line. This is cumbersome, but while it can probably be avoided, it
     * also hard to do so efficiently */
    if ((basepath || subdirlist) && !filelist)
    {
      fprintf(stderr, "Error, -basepath / -subdirlist only valid with -filelist\n");
      usage(pl, argv0);
    }
    if ((filelist != NULL) + (argc != 0) != 1) {
      fprintf(stderr, "Error, provide either -filelist or freeform file names\n");
      usage(pl, argv0);
    }
    if (nrelmax == 0)
    {
      fprintf(stderr, "Error, missing -nrels command line argument "
                      "(or nrels = 0)\n");
      usage(pl, argv0);
    }
    if (nprimemax == 0)
    {
      fprintf(stderr, "Error, missing -nprimes command line argument "
                      "(or nprimes = 0)\n");
      usage(pl, argv0);
    }
    if (min_index > nprimemax)
    {
      fprintf(stderr, "Error, missing -minindex command line argument "
                      "or (minindex > nprimes)\n");
      usage(pl, argv0);
    }
    /* If nrels or nprimes > 2^32, then we need index_t to be 64-bit */
    if (((nprimemax >> 32) != 0 || (nrelmax >> 32) != 0) && sizeof(index_t) < 8)
    {
      fprintf(stderr, "Error, -nrels or -nprimes is too large for a 32-bit "
                      "program\nSee #define index_size in typedefs.h\n");
      exit(EXIT_FAILURE);
    }

    /* Printing relevant information */
    fprintf(stderr, "Weight function used during clique removal:\n"
                    "  0     1     2     3     4     5     6     7\n");
    for (k = 0; k < 8; k++)
      fprintf(stderr, "%0.3f ", weight_function_clique((weight_t) k));
    fprintf(stderr, "\n");

    fprintf(stderr, "Number of relations is %" PRIu64 "\n", nrelmax);
    fprintf(stderr, "Number of prime ideals below the two large prime bounds: "
                    "%" PRIu64 "\n", nprimemax);
    /*}}}*/

    /* {{{ Allocating memory. We are keeping track of the total
     * malloc()'ed amount only for informational purposes */

    /* {{{ Some macros for tracking memory-consuming variables */
#define ALLOC_VERBOSE_MALLOC(type_, variable_, amount_) do {			\
    variable_ = (type_ *) malloc(amount_ * sizeof(type_));		\
    ASSERT_ALWAYS(variable_ != NULL);					\
    size_t cur_alloc = amount_ * sizeof(type_);                         \
    tot_alloc_bytes += cur_alloc;                                       \
    fprintf(stderr,							\
            "Allocated " #variable_ " of %zuMB (total %zuMB so far)\n",	\
	    cur_alloc >> 20, tot_alloc_bytes >> 20);                  	\
} while (0)

#define ALLOC_VERBOSE_CALLOC(type_, variable_, amount_) do {		\
    ALLOC_VERBOSE_MALLOC(type_, variable_, amount_);                    \
    /* Do this now so that we crash early if kernel overcommitted memory\
     */									\
    memset(variable_, 0, amount_);					\
} while (0)

#define ALLOC_VERBOSE_BIT_VECTOR(variable_, amount_) do {			\
    bit_vector_init(variable_, amount_);				\
    ASSERT_ALWAYS(variable_->p != NULL);				\
    size_t cur_alloc = bit_vector_memory_footprint(variable_);		\
    tot_alloc_bytes += cur_alloc;					\
    fprintf(stderr,                                                     \
            "Allocated rel_used of %zuMB (total %zuMB so far)\n",       \
	    cur_alloc >> 20, tot_alloc_bytes >> 20);			\
} while (0)
    /* }}} */

    ALLOC_VERBOSE_BIT_VECTOR(Tbv, nrelmax);
    ALLOC_VERBOSE_MALLOC(index_t, sum2_index, nprimemax);


    set_antebuffer_path(argv0, param_list_lookup_string(pl, "path_antebuffer"));
    /* }}} */

    /*{{{ build the list of input files from the given args
     * If no filelist is given, files are on the command-line.
     * If no subdirlist is given, files are easily construct from basepath and
     * filelist.
     * If subdirlist is given, it is a little bit trickier, see
     * filelist_from_file_with_subdirlist for more details. */
    if (!filelist)
      input_files = argv;
    else if (!subdirlist)
      input_files = filelist_from_file(basepath, filelist, 0);
    else
      input_files = filelist_from_file_with_subdirlist(basepath, filelist, 
                                                       subdirlist);
    /*}}}*/

    /****************** Begin interesting stuff *************************/
    purge_data pd;

    memset(pd, 0, sizeof(purge_data));

    ALLOC_VERBOSE_MALLOC(index_t*, rel_compact, nrelmax);
    ALLOC_VERBOSE_CALLOC(weight_t, ideals_weight, nprimemax);

    fprintf(stderr, "Pass 1, reading and storing ideals with index h >= "
            "%" PRIu64 "\n", min_index);

    nrels = nrelmax;

    /* first pass over relations in files */
    /* Note: Now that we no longer take a bitmap on input, all
     * relations are considered active at this point, so that we
     * do not need to pass a bitmap to filter_rels */
    pd->info.nrels = filter_rels(
            input_files,
            (filter_rels_callback_t) &insert_rel_into_table, pd,
            EARLYPARSE_NEED_PRIMES | EARLYPARSE_NEED_NB,
            NULL, NULL);

    if (pd->info.nrels != nrels) {
	fprintf(stderr,
		"Error, -nrels value should match the number of scanned "
		"relations\nexpected %" PRIu64 " relations, found %" PRIu64
		"\n", nrelmax, pd->info.nrels);
        abort();
    }

    nprimes = pd->info.nprimes;

    tot_alloc_bytes += get_my_malloc_bytes();
    fprintf(stderr, "Allocated rel_compact[i] %zuMB (total %zuMB so far)\n",
	    get_my_malloc_bytes() >> 20, tot_alloc_bytes >> 20);
#ifdef STAT
    {
	size_t tmp = ((uint64_t) buf_arg.W + nrelmax) * sizeof(index_t);
	double ratio = 100.0 * (double) (((double) tmp) /
					 ((double) get_my_malloc_bytes()));
	fprintf(stderr,
		"STAT: W_active=%1.0f\nSTAT: Should take %zuMB in memory, "
		"take %zuMB (%.2f %%)\n", buf_arg.W, tmp >> 20,
		get_my_malloc_bytes() >> 20, ratio);
    }
#endif

    ALLOC_VERBOSE_BIT_VECTOR(rel_used, nrels);
    bit_vector_set(rel_used, 1);

    singletons_and_cliques_removal(&nrels, &nprimes);
    if (nrels < nprimes) {
      fprintf(stderr, "number of relations < number of ideals\n");
      exit(2);
    }
    if (nrels == 0 || nprimes == 0) {
      fprintf(stderr, "number of relations or number of ideals is 0\n");
      exit(2);
    }

    /* free rel_compact[i] and rel_compact. We no longer need them */
    my_malloc_free_all();
    tot_alloc_bytes -= get_my_malloc_bytes();
    fprintf(stderr, "Freed rel_compact[i] %zuMB (total %zuMB so far)\n",
            get_my_malloc_bytes() >> 20, tot_alloc_bytes >> 20);

    free(rel_compact);
    size_t tmp = (nrelmax * sizeof(index_t *));
    tot_alloc_bytes -= tmp;
    fprintf(stderr, "Freed rel_compact %zuMB (total %zuMB so far)\n",
            tmp >> 20, tot_alloc_bytes >> 20);

    /* reread the relation files and convert them to the new coding */
    fprintf(stderr, "Storing remaining relations...\n");

    if (!(pd->fd[0] = fopen_maybe_compressed(purgedname, "w"))) {
	fprintf(stderr, "Error, cannot open file %s for writing.\n",
		purgedname);
	exit(1);
    }

    if (deletedname != NULL) {
	if (!(pd->fd[1] = fopen_maybe_compressed(deletedname, "w"))) {
	    fprintf(stderr, "Error, cannot open file %s for writing.\n",
		    deletedname);
	    exit(1);
	}
    }

    /* Write the header line for the file of remaining relations:
     * compute last index i such that ideals_weight[i] != 0
     */
    {
	uint64_t last_used = nprimemax - 1;
	while (ideals_weight[last_used] == 0)
	    last_used--;

	fprintf(pd->fd[0], "# %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", nrels,
		last_used + 1, nprimes);
    }

    /* second pass over relations in files */
    filter_rels(
            input_files,
            (filter_rels_callback_t) &thread_print, pd,
            EARLYPARSE_NEED_LINE | EARLYPARSE_NEED_NB,
            rel_used, NULL);

    /* write final values to stdout */
    fprintf(stdout, "Final values:\nnrels=%" PRIu64 " nprimes=%" PRIu64 " "
            "excess=%" PRId64 "\nweight=%1.0f weight*nrels=%1.2e\n",
            nrels, nprimes, ((int64_t) nrels) - nprimes, pd->info.W,
            pd->info.W * (double) nrels);
    fflush(stdout);

    /* Free allocated stuff */
    if (ideals_weight != NULL)
	free(ideals_weight);
    ideals_weight = NULL;
    if (sum2_index != NULL)
	free(sum2_index);
    sum2_index = NULL;

    bit_vector_clear(rel_used);
    bit_vector_clear(Tbv);

    if (filelist)
        filelist_clear(input_files);

    fclose_maybe_compressed(pd->fd[0], purgedname);
    if (pd->fd[1]) fclose_maybe_compressed(pd->fd[1], deletedname);

    /* print usage of time and memory */
    print_timing_and_memory(wct0);

    param_list_clear(pl);

    return 0;
}
