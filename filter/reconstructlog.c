#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>   /* for _O_BINARY */

#include "filter_common.h"

#ifdef FOR_FFS
#include "utils-ffs.h"
#endif
#define DEBUG 0

/* to use in the GF(p^n) computations, uncomment the following */
//#define FOR_GFPN

stats_data_t stats; /* struct for printing progress */

/*********************** mutex for multi threaded version ********************/
/* used as mutual exclusion lock for reading the status of logarithms */
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static void
mutex_lock(pthread_mutex_t *lock)
{
  pthread_mutex_lock (lock);
}

static void
mutex_unlock(pthread_mutex_t *lock)
{
  pthread_mutex_unlock (lock);
}

/**** Relations structure used for computing the logarithms from the rels ****/
typedef struct
{
  weight_t nb_unknown;
  ideal_merge_t *unknown;
  mpz_t log_known_part;
} log_rel_t;

/* Init a table of log_rel_t of size nrels (malloc + all mpz_init) */
static log_rel_t *
log_rel_init (uint64_t nrels)
{
  log_rel_t *rels;
  uint64_t i;
  rels = (log_rel_t *) malloc (nrels * sizeof (log_rel_t));
  FATAL_ERROR_CHECK(rels == NULL, "Cannot allocate memory");
  memset(rels, 0, nrels * sizeof(log_rel_t));
  for (i = 0; i < nrels; i++)
    mpz_init(rels[i].log_known_part);
  return rels;
}

/* Free what is allocated by log_rel_init */
static void
log_rel_free (log_rel_t *rels, uint64_t nrels)
{
  uint64_t i;
  my_malloc_free_all();
  for (i = 0; i < nrels; i++)
    mpz_clear (rels[i].log_known_part);
  free(rels);
}

/***** Light relations structure (for constructing the dependancy graph) *****/
typedef struct
{
  weight_t len;
  index_t *needed;
} light_rel;

typedef light_rel *light_rels_t;

/* Init a light_rels_t of size nrels */
static light_rels_t
light_rels_init (uint64_t nrels)
{
  light_rels_t rels;
  rels = (light_rels_t) malloc (nrels * sizeof (light_rel));
  FATAL_ERROR_CHECK(rels == NULL, "Cannot allocate memory");
  memset(rels, 0, nrels * sizeof(light_rel));
  return rels;
}

/* Free what is allocated by light_rel_init (and during reading for the
 * "unknown" array in the light_rel structure */
static void
light_rels_free (light_rels_t rels)
{
  my_malloc_free_all();
  free(rels);
}

/****************** Struct for the tab of log (mpz_t *) *********************/
struct logtab_struct
{
  uint64_t nprimes;
  uint64_t nknown;
  unsigned int nbsm;
  mpz_ptr q;
  mpz_t *tab;
};

typedef struct logtab_struct logtab_t[1];
typedef struct logtab_struct * logtab_ptr;

static void
logtab_init (logtab_t log, uint64_t nprimes, unsigned int nbsm, mpz_t q)
{
  log->tab = NULL;
  uint64_t size = nprimes + nbsm;
  size_t q_nbits = mpz_size(q) * GMP_LIMB_BITS;

  log->nprimes = nprimes;
  log->nknown = 0;
  log->nbsm = nbsm;
  log->q = q;
  log->tab = (mpz_t *) malloc (size * sizeof(mpz_t));
  FATAL_ERROR_CHECK(log->tab == NULL, "Cannot allocate memory");
  for (uint64_t i = 0; i < size; i++)
  {
    mpz_init2 (log->tab[i], q_nbits);
    mpz_set_si(log->tab[i], -1);
  }
}

static void
logtab_insert (logtab_t log, index_t h, mpz_t logvalue)
{
  if (mpz_cmp_ui (logvalue, 0) < 0)
  {
    fprintf (stderr, "Warning, log is negative for h = %" PRid "\n", h);
    mpz_mod (logvalue, logvalue, log->q);
  }
  else if (mpz_cmp (logvalue, log->q) >= 0)
  {
    fprintf (stderr, "Warning, log >= q for h = %" PRid "\n", h);
    mpz_mod (logvalue, logvalue, log->q);
  }
  if (mpz_cmp_ui (logvalue, 0) == 0)
    fprintf (stderr, "Warning, log is zero for h = %" PRid "\n", h);

  if (mpz_sgn (log->tab[h]) >= 0) //already known
  {
    ASSERT_ALWAYS (mpz_cmp (log->tab[h], logvalue) == 0);
  }
  else
  {
    mpz_set (log->tab[h], logvalue);
    if (h < log->nprimes) // log of SM columns are not taken into account
      log->nknown++;
  }
}

static void
logtab_free (logtab_t log)
{
  uint64_t size = log->nprimes + log->nbsm;
  for (uint64_t i = 0; i < size; i++)
    mpz_clear (log->tab[i]);
  free (log->tab);
}

/************ Struct used for reading rels files with process_rels ***********/
typedef struct
{
  log_rel_t *rels;
  logtab_ptr log;
  const char * abunitsfilename;
} read_data_t;

/* Init a read_data_t structure for nrels rels and nprimes primes. Assume the
 * table of log is already allocated, but not the table of log_rel_t     */
static void
read_data_init (read_data_t *data, logtab_t log, uint64_t nrels, 
		const char * abunitsfilename)
{
  data->rels = log_rel_init(nrels);
  data->log = log;
  data->abunitsfilename = abunitsfilename;
}

static void
read_data_free (read_data_t *data, uint64_t nrels)
{
  log_rel_free (data->rels, nrels);
}

/***************** utils functions for adding logarithms *********************/
/* Compute a <- (a + l*e) mod q  , with e being a int32_t */
static inline void
mpz_add_log_mod_si (mpz_t a, mpz_t l, int32_t e, mpz_t q)
{
  if (e > 0)
    mpz_addmul_ui (a, l, e);
  else
    mpz_submul_ui (a, l, -e);
  mpz_mod (a, a, q);
}

/* Compute a <- (a + l*e) mod q  , with e being a mpz_t */
static inline void
mpz_add_log_mod_mpz (mpz_ptr a, mpz_t l, mpz_t e, mpz_t q)
{
  mpz_t t;
  mpz_init(t);
  mpz_mul(t, l, e);  // t <- l*e
  mpz_mod(t, t, q);  // t <- l*e mod q
  mpz_add(t, a, t);  // t <- l*e mod q + a
  mpz_mod (a, t, q); // a <- (a+l*e) mod q
  mpz_clear(t);
}

/************************ Handling of the SMs *******************************/
/* number of SM that must be used for side 1. Must be 0 for FFS */
unsigned int nbsm1 = 0;
/* for side 0 */
#ifdef FOR_GFPN
unsigned int nbsm0 = 0;
#endif

#ifndef FOR_FFS /* Not needed for FFS */
mpz_t smexp1; /* exponent for SM */
# ifdef FOR_GFPN
mpz_t smexp0; /* exponent for SM on the 0 side */
mpz_poly_ptr F0;
# endif
mpz_poly_ptr F1;
mpz_t q2;    /* q^2 */
mpz_t invq2;
mpz_t *smlog;

/* Init the data needed for SM computation, given the polynomial,
 * the group order and the exponent */
static void
sm_data_init (mpz_t q, logtab_t log)
{
  mpz_init(q2);
  mpz_mul(q2, q, q);

  mpz_init(invq2);
  barrett_init(invq2, q2);

  smlog = &(log->tab[log->nprimes]);
}

static void
sm_data_free ()
{
  mpz_clear(q2);
  mpz_clear(invq2);
}

#ifdef FOR_GFPN
/* Very naive code for finding the valuations of the units. */
static int
get_units_for_ab(int* u, int64_t a, uint64_t b, const char * abunitsfilename)
{
    FILE *in = fopen(abunitsfilename, "r");
    int64_t aa;
    uint64_t bb;
    int32_t ok = 0;
    while(fscanf(in, "%" PRId64 " %" PRIu64, &aa, &bb) != EOF){
	for(int i = 0; i < (int)nbsm0; i++)
	    fscanf(in, "%d", u+i);
	if(aa == a && bb == b){
	    printf("# GOTCHA a and b\n");
	    ok = 1;
	    break;
	}
    }
    fclose(in);
    if(ok == 0){
	fprintf(stderr, "GASP: no eps found ");
	fprintf(stderr, "for a=%" PRId64 " b=%" PRIu64 "\n", a, b);
	exit(-1);
    }
    return ok;
}

/* we need retrieve the exponents of the units in a-b*x. */
static inline void
add_unit_contribution (mpz_ptr l, int64_t a, uint64_t b, mpz_t q,
		       const char * abunitsfilename)
{
    int u[10];
    if(get_units_for_ab(u, a, b, abunitsfilename)){
	/* add contrib */
	for(int i = 0; i < (int)nbsm0; i++)
	    mpz_add_log_mod_si (l, smlog[nbsm1+i], u[i], q);
    }
}
#endif // FOR_GFPN

/* Given a and b, compute the SM and add the contribution to l */
static inline void
add_sm_contribution (mpz_ptr l, int64_t a, uint64_t b, mpz_t q,
		     mpz_poly_ptr F, mpz_t smexp, unsigned int nbsm,
		     mpz_t *smlg)
{
  mpz_poly_t SMres;
  int degF = F->deg;
  mpz_poly_init(SMres, degF);
  SMres->deg = 0;
  mpz_poly_setcoeff_si(SMres, 0, 1);
  sm_single_rel(SMres, a, b, F, smexp, q, q2, invq2);
  unsigned int i;
  for (i = degF - SMres->deg - 1; i < nbsm; i++)
    mpz_add_log_mod_mpz (l, smlg[i], SMres->coeff[degF-1-i], q);
  mpz_poly_clear(SMres);
}

/* TODO: everything is easily adaptable for sm/sm or sm/unit or unit/unit.
   We need just throw away global variables smexp, F, etc. At least replacing
   them by smexp1 and F1, so as to function-ify the above function.
*/
static inline void
add_sm_contributions (mpz_ptr l, int64_t a, uint64_t b, mpz_t q,
		      const char * abunitsfilename MAYBE_UNUSED)
{
  add_sm_contribution(l, a, b, q, F1, smexp1, nbsm1, smlog);
#ifdef FOR_GFPN
  if(mpz_sgn(smexp0) == 0)
    add_unit_contribution(l, a, b, q, abunitsfilename);
  else
    add_sm_contribution(l, a, b, q, F0, smexp0, nbsm0, smlog+nbsm1);
#endif // FOR_GFPN
}

/* Callback function called by filter_rels in compute_log_from_rels */
void *
thread_sm (void * context_data, earlyparsed_relation_ptr rel)
{
  read_data_t *data = (read_data_t *) context_data;
  log_rel_t *lrel = &(data->rels[rel->num]);

  if (nbsm1 > 0)
      add_sm_contributions (lrel->log_known_part, rel->a, rel->b, data->log->q,
			    data->abunitsfilename);

  return NULL;
}
#endif /* ifndef FOR_FFS */

/****************** Computation of missing logarithms ************************/
/* Callback function called by filter_rels in compute_log_from_rels */
void *
thread_insert (void * context_data, earlyparsed_relation_ptr rel)
{
  read_data_t *data = (read_data_t *) context_data;
  log_rel_t *lrel = &(data->rels[rel->num]);
  unsigned int next = 0;
  ideal_merge_t buf[REL_MAX_SIZE];

  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    exponent_t e = rel->primes[i].e;

	  if (mpz_sgn(data->log->tab[h]) >= 0)
      mpz_add_log_mod_si (lrel->log_known_part, data->log->tab[h], e,
                          data->log->q);
    else
      buf[next++] = (ideal_merge_t) {.id = h, .e = e};
  }

  lrel->unknown = idealmerge_my_malloc (next);
  lrel->nb_unknown = next;
  memcpy(lrel->unknown, buf, next * sizeof(ideal_merge_t));

  return NULL;
}

/* Return the number of unknown logarithms in the relation.
 * rels[i].nb_unknown may not be up-to-date (can only be greater than the actual
 * value) */
static inline weight_t
nb_unknown_log (read_data_t *data, uint64_t i)
{
  weight_t j, k, len = data->rels[i].nb_unknown;
  ideal_merge_t *p = data->rels[i].unknown;

  for (j = 0, k = 0; k < len; k++)
  {
    mutex_lock (&lock);
    int sgn = mpz_sgn(data->log->tab[p[k].id]);
    mutex_unlock (&lock);

	  if (sgn < 0) // we do not know the log if this ideal
    {
      if (j != k)
        p[j] = p[k];
      j++;
    }
    else // We know this log, add it to log_know_part
      mpz_add_log_mod_si (data->rels[i].log_known_part, data->log->tab[p[k].id],
                          p[k].e, data->log->q);
  }

  data->rels[i].nb_unknown = j;
  return j;
}

/* In a relation with 1 missing logarithm of exponent e, compute its values,
 * i.e. compute   dest <- (-vlog / e) mod q
 * Return 0 if dest was already known (i.e. computed between the call to
 * nb_unknown_log and the call to this functions), return 1 otherwise.
 */
static inline unsigned int
compute_missing_log (mpz_t dest, mpz_t vlog, int32_t e, mpz_t q)
{
  unsigned int ret;
  mpz_t tmp;
  mpz_init_set_si (tmp, e);
  mpz_invert (tmp, tmp, q);
  mpz_neg (vlog, vlog);
  mpz_mul (vlog, vlog, tmp);
  mpz_mod (tmp, vlog, q);

  mutex_lock (&lock);

  if (mpz_sgn(dest) < 0) // we do not already know the log
  {
    mpz_set (dest, tmp);
    ret = 1;
  }
  else // log was already computed by another thread
  {
    ret = 0;
  }
  mutex_unlock (&lock);

  mpz_clear (tmp);
  return ret;
}

/* Compute all missing logarithms for relations in [start,end[.
 * Return the number of computed logarithms */
static uint64_t
log_do_one_part_of_iter (read_data_t *data, bit_vector not_used, uint64_t start,
                         uint64_t end)
{
  uint64_t i, computed = 0;

  for (i = start; i < end; i++)
  {
    if (bit_vector_getbit(not_used, (size_t) i))
    {
      weight_t nb = nb_unknown_log(data, i);
      if (nb <= 1)
      {
        bit_vector_clearbit(not_used, (size_t) i);
        mpz_ptr vlog = data->rels[i].log_known_part;
        if (nb == 0 && mpz_cmp_ui (vlog, 0) != 0)
        {
          gmp_fprintf (stderr, "Error, no unknown log in rel %" PRIu64 " and sum"
                       " of log is not zero, sum is:\n%Zd\n", i, vlog);
          exit (EXIT_FAILURE);
        }
        else if (nb == 1)
        {
          ideal_merge_t ideal = data->rels[i].unknown[0];
          computed +=
             compute_missing_log (data->log->tab[ideal.id], vlog, ideal.e,
                                  data->log->q);
        }
      }
    }
  }

  return computed;
}

#define log_do_one_iter_mono(d, bv, n) log_do_one_part_of_iter (d, bv, 0, n)

/************************** Dependancy graph *********************************/
typedef struct
{
  uint8_t state;
  uint64_t i;
} node_dep_t;

typedef struct
{
  uint64_t size;
  node_dep_t *tab;
} graph_dep_t;

/* Macro for state of node_dep_t. UNKNOWN must be 0. */
#define NODE_DEP_LOG_UNKNOWN 0
#define NODE_DEP_LOG_KNOWN_FROM_LOGFILE 1
#define NODE_DEP_LOG_RECONSTRUCTED 2

#define GRAPH_DEP_IS_LOG_UNKNOWN(G, h) (G.tab[h].state == NODE_DEP_LOG_UNKNOWN)

graph_dep_t
graph_dep_init (uint64_t size)
{
  node_dep_t *tab = NULL;
  tab = (node_dep_t *) malloc (sizeof(node_dep_t) * size);
  FATAL_ERROR_CHECK(tab == NULL, "Cannot allocate memory");
  memset (tab, 0, sizeof(node_dep_t) * size);
  return (graph_dep_t) {.size = size, .tab = tab};
}

void
graph_dep_clear (graph_dep_t G)
{
  free(G.tab);
}

/* Set G[h].state accordingly to log[h] values */
void
graph_dep_set_log_already_known (graph_dep_t G, logtab_t log)
{
  for (uint64_t h = 0; h < log->nprimes; h++)
  {
    if (mpz_sgn(log->tab[h]) >= 0)
      G.tab[h].state = NODE_DEP_LOG_KNOWN_FROM_LOGFILE;
  }
}

uint64_t
graph_dep_needed_rels_from_index (graph_dep_t G, index_t h, light_rels_t rels,
                                  bit_vector needed_rels)
{
  if (G.tab[h].state == NODE_DEP_LOG_UNKNOWN)
  {
    fprintf (stderr, "Error: logarithms of %" PRid" cannot be reconstructed "
                     "from this set of relations. Abort!\n", h);
    abort();
  }
  else if (G.tab[h].state == NODE_DEP_LOG_KNOWN_FROM_LOGFILE)
  {
    /* We know the wanted logarithm from linear algebra, no new relation is
     * necessary */
#if DEBUG >= 1
    fprintf (stderr, "DEBUG: h = %" PRid " is known from logfile\n", h);
#endif
    return 0;
  }
  else
  {
    uint64_t relnum = G.tab[h].i;
    bit_vector_setbit (needed_rels, relnum);
    uint64_t nadded = 1;
    weight_t nb_needed = rels[relnum].len;

#if DEBUG >= 1
    fprintf (stderr, "DEBUG: h = %" PRid " can be reconstructed\n", h);
    fprintf (stderr, "DEBUG:     relation %" PRIu64 " added\n", relnum);
    fprintf (stderr, "DEBUG:     depends of %u others logs\n", nb_needed-1);
#endif

    for (weight_t j = 0; j < nb_needed; j++)
    {
      index_t hh = rels[relnum].needed[j];
      if (!bit_vector_getbit (needed_rels, G.tab[hh].i))
      {
        nadded += graph_dep_needed_rels_from_index (G, hh, rels, needed_rels);
      }
    }
    return nadded;
  }
}

/* Structure containing the data necessary for reading rels for dep. graph */
typedef struct
{
  light_rels_t rels;
  graph_dep_t G;
} dep_read_data_t;

/* Callback function called by filter_rels in compute_needed_rels */
void *
dep_thread_insert (void * context_data, earlyparsed_relation_ptr rel)
{
  dep_read_data_t *data = (dep_read_data_t *) context_data;
  light_rels_t lrel = &(data->rels[rel->num]);
  unsigned int next = 0;
  index_t buf[REL_MAX_SIZE];

  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
	  if (GRAPH_DEP_IS_LOG_UNKNOWN(data->G, h))
      buf[next++] = h;
  }

  lrel->needed = index_my_malloc (next);
  lrel->len = next;
  memcpy(lrel->needed, buf, next * sizeof(ideal_merge_t));

  return NULL;
}

/* Return the number of unknown logarithms in the relation and put in h the last
 * unknown logarithm of the relations (usefull when the number of unknown
 * logarithms is 1) */
static inline weight_t
dep_nb_unknown_log (dep_read_data_t *data, uint64_t i, index_t *h)
{
  weight_t k, nb;
  index_t *p = data->rels[i].needed;

  for (nb = 0, k = 0; k < data->rels[i].len; k++)
  {
    mutex_lock (&lock);
    int unknow = GRAPH_DEP_IS_LOG_UNKNOWN (data->G, p[k]);
    mutex_unlock (&lock);

	  if (unknow) // we do not know the log if this ideal
    {
      nb++;
      *h = p[k];
    }
  }
  return nb;
}

/* Compute all dependancies for relations in [start,end[.
 * Return the number of dependancies found */
static uint64_t
dep_do_one_part_of_iter (dep_read_data_t *data, bit_vector not_used,
                         uint64_t start, uint64_t end)
{
  uint64_t i, computed = 0;

  for (i = start; i < end; i++)
  {
    if (bit_vector_getbit(not_used, (size_t) i))
    {
      index_t h;
      weight_t nb = dep_nb_unknown_log(data, i, &h);
      if (nb <= 1)
      {
        bit_vector_clearbit(not_used, (size_t) i);
        if (nb == 1)
        {
          data->G.tab[h].state = NODE_DEP_LOG_RECONSTRUCTED;
          data->G.tab[h].i = i;
          computed++;
        }
      }
    }
  }

  return computed;
}

#define dep_do_one_iter_mono(d, bv, n) dep_do_one_part_of_iter (d, bv, 0, n)

/******************** Code for multi thread version **************************/
#define SIZE_BLOCK 1024

typedef struct {
  void *data;
  bit_vector_ptr not_used;
  uint64_t offset;
  uint64_t nb;
  uint64_t computed;
  int version; /* 0 means call dep_* , 1 means call log_* */
} thread_info;

void * thread_start(void *arg)
{
  thread_info *ti = (thread_info *) arg;
  bit_vector_ptr not_used = ti->not_used;
  uint64_t start = ti->offset;
  uint64_t end = start + ti->nb;

  if (ti->version == 0)
  {
    dep_read_data_t *data = (dep_read_data_t *) ti->data;
    ti->computed = dep_do_one_part_of_iter (data, not_used, start, end);
  }
  else
  {
    read_data_t *data = (read_data_t *) ti->data;
    ti->computed = log_do_one_part_of_iter (data, not_used, start, end);
  }
  return NULL;
}

static uint64_t
do_one_iter_mt (void *data, bit_vector not_used, int nt, uint64_t nrels,
                int version)
{
  // We'll use a rotating buffer of thread id.
  pthread_t *threads;
  threads = (pthread_t *) malloc( nt * sizeof(pthread_t));
  int active_threads = 0;  // number of running threads
  int threads_head = 0;    // next thread to wait / restart.

  // Prepare the main loop
  uint64_t i = 0; // counter of relation.
  uint64_t computed = 0;

  // Arguments for threads
  thread_info *tis;
  tis = (thread_info *) malloc( nt * sizeof(thread_info));
  for (int i = 0; i < nt; ++i)
  {
    tis[i].data = data;
    tis[i].not_used = not_used;
    tis[i].version = version;
    // offset and nb must be adjusted.
  }

  // Main loop
  while ((i < nrels) || (active_threads > 0))
  {
    // Start / restart as many threads as allowed
    if ((active_threads < nt) && (i < nrels))
    {
      tis[threads_head].offset = i;
      tis[threads_head].nb = MIN(SIZE_BLOCK, nrels-i);
      pthread_create(&threads[threads_head], NULL,
          &thread_start, (void *)(&tis[threads_head]));
      i += SIZE_BLOCK;
      active_threads++;
      threads_head++;
      if (threads_head == nt)
        threads_head = 0;
      continue;
    }
    // Wait for the next thread to finish in order to print result.
    pthread_join(threads[threads_head], NULL);
    active_threads--;
    computed += tis[threads_head].computed;

    // If we are at the end, no job will be restarted, but head still
    // must be incremented.
    if (i >= nrels)
    {
      threads_head++;
      if (threads_head == nt)
        threads_head = 0;
    }
  }

  free(tis);
  free(threads);
  return computed;
}

/* Compute all missing logarithms possible. Run through all the relations once.
 * Multi thread version
 * Return the number of computed logarithms */
static inline uint64_t
log_do_one_iter_mt (read_data_t *d, bit_vector bv, int nt, uint64_t nrels)
{
  return do_one_iter_mt ((void *) d, bv, nt, nrels, 1);
}

/* Compute all missing logarithms possible. Run through all the relations once.
 * Multi thread version
 * Return the number of computed logarithms */
static inline uint64_t
dep_do_one_iter_mt (dep_read_data_t *d, bit_vector bv, int nt, uint64_t nrels)
{
  return do_one_iter_mt ((void *) d, bv, nt, nrels, 0);
}

/***************** Important functions called by main ************************/
/* Read the logarithms computed by the linear algebra */
static void
read_log_format_LA (logtab_t log, const char *logfile, const char *idealsfile)
{
  uint64_t i, ncols, col;
  index_t h;
  mpz_t tmp_log;
  FILE *flog = NULL, *fid = NULL;

  printf ("# Reading logarithms in LA format from %s\n", logfile);
  printf ("# Reading links between matrix columns and ideals from %s\n",
                                                                idealsfile);
  fflush(stdout);
  flog = fopen_maybe_compressed (logfile, "r");
  FATAL_ERROR_CHECK(flog == NULL, "Cannot open file for reading logarithms");
  fid = fopen_maybe_compressed (idealsfile, "r");
  FATAL_ERROR_CHECK(flog == NULL, "Cannot open ideals file");

  if (fscanf (fid, "# %" SCNu64 "\n", &ncols) != 1)
  {
    fprintf(stderr, "Error while reading first line of %s\n", idealsfile);
    abort();
  }

  mpz_init (tmp_log);
  i = 0;
  stats_init (stats, stdout, &i, nbits(ncols)-5, "Read", "logarithms", "", "logs");
  while (fscanf (fid, "%" SCNu64 " %" PRid "\n", &col, &h) == 2)
  {
    FATAL_ERROR_CHECK (col >= ncols, "Too big value of column number");
    FATAL_ERROR_CHECK (h >= log->nprimes, "Too big value of index");

    int ret = gmp_fscanf (flog, "%Zd\n", tmp_log);
    FATAL_ERROR_CHECK (ret != 1, "Error in file containing logarithms values");

    ASSERT_ALWAYS (col == i);
    logtab_insert (log, h, tmp_log);
    i++;
    if (stats_test_progress (stats))
      stats_print_progress (stats, i, 0, 0, 0);
  }
  stats_print_progress (stats, i, 0, 0, 1);
  ASSERT_ALWAYS (feof(fid));
  ASSERT_ALWAYS (i == ncols);

  for (unsigned int nsm = 0; nsm < nbsm1; nsm++)
  {
    int ret = gmp_fscanf (flog, "%Zd\n", tmp_log);
    FATAL_ERROR_CHECK (ret != 1, "Error in file containing logarithms values");
    logtab_insert (log, log->nprimes+nsm, tmp_log);
  }
#ifdef FOR_GFPN
  for (unsigned int nu = 0; nu < nbsm0; nu++)
  {
    int ret = gmp_fscanf (flog, "%Zd\n", tmp_log);
    FATAL_ERROR_CHECK (ret != 1, "Error in file containing logarithms values");
    logtab_insert (log, log->nprimes+nbsm1+nu, tmp_log);
  }
#endif
  while (gmp_fscanf (flog, "%Zd\n", tmp_log) == 1)
    printf ("Warning, line %" PRIu64 " is ignored\n", i++);
  ASSERT_ALWAYS (feof(flog));

  if (nbsm1)
    printf ("# Logarithms for %u SM1 columns were also read\n", nbsm1);
#ifdef FOR_GFPN
  if (nbsm0)
    printf ("# Logarithms for %u SM0 columns were also read\n", nbsm0);
#endif
  mpz_clear (tmp_log);
  fclose_maybe_compressed (flog, logfile);
  fclose_maybe_compressed (fid, idealsfile);
}

/* Read the logarithms in output format of reconstructlog */
static void
read_log_format_reconstruct (logtab_t log, MAYBE_UNUSED renumber_t renumb,
                             const char *filename)
{
  uint64_t nread = 0;
  index_t h;
  mpz_t tmp_log;
  FILE *f = NULL;
  int ret;

  printf ("# Reading logarithms in reconstruct format from %s\n", filename);
  fflush(stdout);
  f = fopen_maybe_compressed (filename, "r");
  FATAL_ERROR_CHECK(f == NULL, "Cannot open file for reading");

  mpz_init (tmp_log);
  stats_init (stats, stdout, &nread, nbits(renumb->size)-5, "Read", "logarithms", "",
              "logs");
  if (renumb->add_full_col)
  {
    ret = gmp_fscanf (f, "%" SCNid " added column %Zd\n", &h, tmp_log);
    ASSERT_ALWAYS (ret == 2 && h == 0);
    nread++;
    logtab_insert (log, h, tmp_log);
  }
  for (int i = 0; i < renumb->bad_ideals.n; i++)
  {
    for (int k = 0; k < renumb->bad_ideals.nb[i]; k++)
    {
      ret = gmp_fscanf (f, "%" SCNid " bad ideals %Zd\n", &h, tmp_log);
      ASSERT_ALWAYS (ret == 2);
      nread++;
      logtab_insert (log, h, tmp_log);
    }
  }
  while (gmp_fscanf (f, "%" SCNid " %*" SCNpr " %*d %*s %Zd\n", &h, tmp_log)
          == 2)
  {
    nread++;
    logtab_insert (log, h, tmp_log);
    if (stats_test_progress (stats))
      stats_print_progress (stats, nread, 0, 0, 0);
  }
  stats_print_progress (stats, nread, 0, 0, 1);

  for (unsigned int nsm = 0; nsm < log->nbsm; nsm++)
  {
    unsigned int n;
    if (nsm == 0) /* h was already red by previous gmp_fscanf */
    {
      ret = gmp_fscanf (f, "SM col %u %Zd\n", &n, tmp_log);
      ASSERT_ALWAYS (ret == 2);
    }
    else
    {
      ret = gmp_fscanf (f, "%" SCNid " SM col %u %Zd\n", &h, &n, tmp_log);
      ASSERT_ALWAYS (ret == 3);
    }
    ASSERT_ALWAYS (n == nsm);
    ASSERT_ALWAYS (h == (index_t) nsm + log->nprimes);
    logtab_insert (log, h, tmp_log);
  }
  ASSERT_ALWAYS (feof(f));

  if (nbsm1)
    printf ("# Logarithms for %u SM columns were also read\n", nbsm1);
  mpz_clear (tmp_log);
  fclose_maybe_compressed (f, filename);
}

/* Write values of the known logarithms. */
static void
write_log (const char *filename, logtab_t log, renumber_t tab, cado_poly poly)
{
  uint64_t i;
  FILE *f = NULL;

  printf ("# Opening %s for writing logarithms\n", filename);
  fflush(stdout);
  f = fopen_maybe_compressed (filename, "w");
  FATAL_ERROR_CHECK(f == NULL, "Cannot open file for writing");

  uint64_t nknown = 0;
  stats_init (stats, stdout, &nknown, nbits(tab->size)-5, "Wrote",
              "known logarithms", "ideals", "logs");
  for (i = 0; i < tab->size; i++)
  {
	  if (mpz_sgn(log->tab[i]) >= 0) // we know the log if this ideal
    {
      nknown++;
      if (tab->table[i] == RENUMBER_SPECIAL_VALUE)
      {
        ASSERT_ALWAYS (mpz_cmp (log->tab[i], log->q) < 0);
        if (i == 0 && tab->add_full_col)
          gmp_fprintf (f, "%" PRid " added column %Zd\n", i, log->tab[i]);
        else
          gmp_fprintf (f, "%" PRid " bad ideals %Zd\n", i, log->tab[i]);
      }
      else
      {
        ASSERT_ALWAYS (mpz_cmp (log->tab[i], log->q) < 0);
        p_r_values_t p, r;
        int side;
        renumber_get_p_r_from_index (tab, &p, &r, &side, i, poly);
        if (side != tab->rat)
          gmp_fprintf (f, "%" PRid " %" PRpr " %d %" PRpr " %Zd\n", i, p, side,
                                                               r, log->tab[i]);
        else
          gmp_fprintf (f, "%" PRid " %" PRpr " %d rat %Zd\n", i, p, side,
                                                                  log->tab[i]);
      }
      if (stats_test_progress (stats))
        stats_print_progress (stats, nknown, i+1, 0, 0);
    }
  }
  stats_print_progress (stats, nknown, tab->size, 0, 1);
  for (unsigned int nsm = 0, i = tab->size; nsm < log->nbsm; nsm++)
  {
    ASSERT_ALWAYS (mpz_sgn (log->tab[i+nsm]) >= 0);
    ASSERT_ALWAYS (mpz_cmp (log->tab[i+nsm], log->q) < 0);
    gmp_fprintf (f, "%" PRid " SM col %u %Zd\n", i+nsm, nsm, log->tab[i+nsm]);
  }

  uint64_t missing = tab->size - nknown;
  printf ("# factor base contains %" PRIu64 " elements\n"
          "# logarithms of %" PRIu64 " elements are known (%.1f%%)\n"
          "# logarithms of %" PRIu64 " elements are missing (%.1f%%)\n",
          tab->size, nknown, 100.0 * nknown / (double) tab->size,
          missing, 100.0 * missing / (double) tab->size);
  fclose_maybe_compressed (f, filename);
  ASSERT_ALWAYS (log->nknown == nknown);
}

/* Given a filename, compute all the possible logarithms of ideals appearing in
 * the file. Return the number of computed logarithms.
 * Modify its first argument bit_vector needed_rels */
static uint64_t
compute_log_from_rels (bit_vector needed_rels,
                       const char *relspfilename, uint64_t nrels_purged,
                       const char *relsdfilename, uint64_t nrels_del,
                       uint64_t nrels_needed, logtab_t log, int nt,
		       const char *abunitsfilename)
{
  double wct_tt0, wct_tt;
  uint64_t total_computed = 0, iter = 0, computed;
  uint64_t nrels = nrels_purged + nrels_del;
  ASSERT_ALWAYS (nrels_needed > 0);
  read_data_t data;
  read_data_init(&data, log, nrels, abunitsfilename);

  /* Reading all relations */
  printf ("# Reading relations from %s and %s\n", relspfilename, relsdfilename);
  if (nrels_needed != nrels)
    printf ("# Parsing only %" PRIu64 " needed relations out of %" PRIu64 "\n",
            nrels_needed, nrels);
#if ! defined (FOR_FFS) && DEBUG >= 1
  printf ("# DEBUG: Using %d thread(s) for thread_sm\n", nt);
#endif
  fflush(stdout);
  char *fic[3] = {(char *) relspfilename, (char *) relsdfilename, NULL};
  struct filter_rels_description desc[3] = {
                   { .f = thread_insert, .arg=&data, .n=1},
#ifndef FOR_FFS
                   { .f = thread_sm,     .arg=&data, .n=nt},
#endif
                   { .f = NULL,          .arg=0,     .n=0}
      };
  filter_rels2 (fic, desc, EARLYPARSE_NEED_AB_HEXA | EARLYPARSE_NEED_INDEX,
                needed_rels, NULL);

  /* computing missing log */
  printf ("# Starting to computing missing logarithms from rels\n");

  /* adjust the number of threads based on the number of needed relations */
  double ntm = ceil((nrels_needed + 0.0)/SIZE_BLOCK);
  if (nt > ntm)
    nt = (int) ntm;

  if (nt > 1)
    printf("# Using multithread version with %d threads\n", nt);
  else
    printf("# Using monothread version\n");

  wct_tt0 = wct_seconds();
  do
  {
    printf ("# Iteration %" PRIu64 ": starting...\n", iter);
    fflush(stdout);
    wct_tt = wct_seconds();

    if (nt > 1)
      computed = log_do_one_iter_mt (&data, needed_rels, nt, nrels);
    else
      computed = log_do_one_iter_mono (&data, needed_rels, nrels);
    total_computed += computed;

    printf ("# Iteration %" PRIu64 ": %" PRIu64 " new logarithms computed\n",
            iter, computed);
    printf ("# Iteration %" PRIu64 " took %.1fs (wall-clock time).\n",
            iter, wct_seconds() - wct_tt);

    iter++;
  } while (computed);

  printf ("# Computing %" PRIu64 " new logarithms took %.1fs (wall-clock "
          "time)\n", total_computed, wct_seconds() - wct_tt0);

  size_t c = bit_vector_popcount(needed_rels);
  if (c != 0)
    fprintf(stderr, "### Warning, %zu relations were not used\n", c);

  read_data_free(&data, nrels);
  return total_computed;
}

/* Given a filename, compute all the relations needed to compute the logarithms
 * appearing in the file.
 * needed_rels should be initialized before calling this function. Its size
 * must be nrels_purged + nrels_del.
 * Output:
 *    bit_vector needed_rels, where bits of needed rels are set.
 *    Return the number of needed_rels.*/
static uint64_t
compute_needed_rels (bit_vector needed_rels,
                     const char *relspfilename, uint64_t nrels_purged,
                     const char *relsdfilename, uint64_t nrels_del,
                     logtab_t log, const char *wanted_filename, int nt)
{
  double wct_tt0, wct_tt;
  uint64_t total_computed = 0, iter = 0, computed;
  uint64_t nrels = nrels_purged + nrels_del;
  graph_dep_t dep_graph = graph_dep_init (log->nprimes);
  light_rels_t rels = light_rels_init (nrels);

  graph_dep_set_log_already_known (dep_graph, log);

  dep_read_data_t data = (dep_read_data_t) {.rels = rels, .G = dep_graph};

  /* Init bit_vector to remember which relations were already used */
  bit_vector_set (needed_rels, 1);

  /* Reading all relations */
  printf ("# Reading relations from %s and %s\n", relspfilename, relsdfilename);
  fflush(stdout);
  char *fic[3] = {(char *) relspfilename, (char *) relsdfilename, NULL};
  filter_rels (fic, (filter_rels_callback_t) &dep_thread_insert, (void *) &data,
               EARLYPARSE_NEED_INDEX, NULL, NULL);

  /* computing dependancies */
  printf ("# Starting to compute dependancies from rels\n");

  /* adjust the number of threads based on the number of relations */
  double ntm = ceil((nrels + 0.0)/SIZE_BLOCK);
  if (nt > ntm)
    nt = (int) ntm;

  if (nt > 1)
    printf("# Using multithread version with %d threads\n", nt);
  else
    printf("# Using monothread version\n");

  wct_tt0 = wct_seconds();
  do
  {
    printf ("# Iteration %" PRIu64 ": starting...\n", iter);
    fflush(stdout);
    wct_tt = wct_seconds();

    if (nt > 1)
      computed = dep_do_one_iter_mt (&data, needed_rels, nt, nrels);
    else
      computed = dep_do_one_iter_mono (&data, needed_rels, nrels);
    total_computed += computed;

    printf ("# Iteration %" PRIu64 ": %" PRIu64 " new dependancies computed\n",
            iter, computed);
    printf ("# Iteration %" PRIu64 " took %.1fs (wall-clock time).\n",
            iter, wct_seconds() - wct_tt);

    iter++;
  } while (computed);

  printf ("# Computing dependancies took %.1fs (wall-clock time)\n",
          wct_seconds() - wct_tt0);

  FILE *f = NULL;
  printf ("# Reading wanted logarithms from %s\n", wanted_filename);
  fflush(stdout);
  f = fopen_maybe_compressed (wanted_filename, "r");
  FATAL_ERROR_CHECK(f == NULL, "Cannot open file for reading");

  bit_vector_set (needed_rels, 0);
  index_t h;
  uint64_t nadded, nrels_necessary = 0, nwanted_log = 0;
  wct_tt = wct_seconds();
  while (fscanf (f, "%" SCNid "\n", &h) == 1)
  {
    FATAL_ERROR_CHECK (h >= log->nprimes, "Too big value of index");
    printf ("# Computing rels necessary for wanted log %" PRid "\n", h);
    fflush(stdout);
    nadded = graph_dep_needed_rels_from_index (dep_graph, h, rels, needed_rels);
    nrels_necessary += nadded;
    printf ("-> %" PRIu64 " needed relations were added (%" PRIu64 " so far)\n",
            nadded, nrels_necessary);
    nwanted_log++;
  }

  fclose_maybe_compressed (f, wanted_filename);
  printf ("# Reading %" PRIu64 " wanted logarithms took %.1fs\n", nwanted_log,
          wct_seconds() - wct_tt);
  printf ("# %" PRIu64 " relations are needed to compute these logarithms\n",
          nrels_necessary);
  ASSERT_ALWAYS (nrels_necessary == bit_vector_popcount (needed_rels));

  light_rels_free (rels);
  graph_dep_clear (dep_graph);
  return nrels_necessary;
}

/********************* usage functions and main ******************************/
static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "log", "input file containing known logarithms");
  param_list_decl_usage(pl, "logformat", "format of input log file: 'LA' or "
                                         "'reconstruct' (default is 'LA')");
  param_list_decl_usage(pl, "gorder", "group order (see sm -gorder parameter)");
  param_list_decl_usage(pl, "out", "output file for logarithms");
  param_list_decl_usage(pl, "renumber", "input file for renumbering table");
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "ideals", "link between matrix cols and ideals "
                                      "(see replay -ideals parameter)");
  param_list_decl_usage(pl, "purged", "file with purged relations "
                                      "(see purge -out parameter)");
  param_list_decl_usage(pl, "relsdel", "file with relations deleted by purge "
                                      "(see purge -outdel parameter)");
  param_list_decl_usage(pl, "nrels", "number of relations (same as purge "
                                     "-nrels parameter)");
  param_list_decl_usage(pl, "partial", "(switch) do not reconstruct everything "
                                       "that can be reconstructed");
#ifndef FOR_FFS
  param_list_decl_usage(pl, "sm", "number of SM to add to relations");
  // TODO: replace smexp by smexp1
  param_list_decl_usage(pl, "smexp", "sm exponent (see sm -smexp1 parameter)");
# ifdef FOR_GFPN
  param_list_decl_usage(pl, "sm0", "number of SM to add to relations");
  param_list_decl_usage(pl, "smexp0", "sm0 exponent (see sm -smexp parameter)");
  param_list_decl_usage(pl, "abunits", "units for all (a, b) pairs from purged and relsdels");
# endif
#endif
  param_list_decl_usage(pl, "mt", "number of threads (default 1)");
  param_list_decl_usage(pl, "wanted", "file containing list of wanted logs");
  param_list_decl_usage(pl, "force-posix-threads", "(switch)");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}


int
main(int argc, char *argv[])
{
  char *argv0 = argv[0];

  renumber_t renumber_table;
  uint64_t nrels_tot = 0, nrels_purged, nideals_purged, nrels_del, nrels_needed;
  uint64_t nprimes;
  int mt = 1;
  int partial = 0;

  mpz_t q;
  logtab_t log;
  cado_poly poly;

  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  argv++,argc--;

  param_list_configure_switch(pl, "partial", &partial);
  param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

#ifdef HAVE_MINGW
  _fmode = _O_BINARY;     /* Binary open for all files */
#endif

  if (argc == 0)
    usage(pl, argv0);

  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
    fprintf (stderr, "Unknown option: %s\n", argv[0]);
    usage(pl, argv0);
  }
  /* print command-line arguments */
  param_list_print_command_line (stdout, pl);
  fflush(stdout);

  mpz_init (q);
  const char * logfilename = param_list_lookup_string(pl, "log");
  const char * logformat = param_list_lookup_string(pl, "logformat");
  const char * idealsfilename = param_list_lookup_string(pl, "ideals");
  const char * relsdfilename = param_list_lookup_string(pl, "relsdel");
  const char * relspfilename = param_list_lookup_string(pl, "purged");
  const char * outfilename = param_list_lookup_string(pl, "out");
  const char * renumberfilename = param_list_lookup_string(pl, "renumber");
  const char * polyfilename = param_list_lookup_string(pl, "poly");
  const char * wantedfilename = param_list_lookup_string(pl, "wanted");
  param_list_parse_uint64(pl, "nrels", &nrels_tot);
  param_list_parse_mpz(pl, "gorder", q);
#ifndef FOR_FFS
  param_list_parse_uint(pl, "sm", &nbsm1);
  mpz_init (smexp1);
  param_list_parse_mpz(pl, "smexp", smexp1);
# ifdef FOR_GFPN
  param_list_parse_uint(pl, "sm0", &nbsm0);
  mpz_init (smexp0);
  param_list_parse_mpz(pl, "smexp0", smexp0);
  const char * abunitsfilename = param_list_lookup_string(pl, "abunits");
# else
  const char * abunitsfilename = NULL;
# endif
#endif
  param_list_parse_int(pl, "mt", &mt);
  const char *path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");

  /* Some checks on command line arguments */
  if (param_list_warn_unused(pl))
  {
    fprintf(stderr, "Error, unused parameters are given\n");
    usage(pl, argv0);
  }

  if (mpz_cmp_ui (q, 0) <= 0)
  {
    fprintf(stderr, "Error, missing -gorder command line argument "
                    "(or gorder <= 0)\n");
    usage (pl, argv0);
  }
#ifndef FOR_FFS
  if (nbsm1 != 0 && mpz_cmp_ui (smexp1, 0) <= 0)
  {
    fprintf(stderr, "Error, missing -smexp1 command line argument "
                    "(or smexp1 <= 0)\n");
    usage (pl, argv0);
  }
#endif
  if (nrels_tot == 0)
  {
    fprintf(stderr, "Error, missing -nrels command line argument "
                    "(or nrels = 0)\n");
    usage (pl, argv0);
  }
  if (logfilename == NULL)
  {
    fprintf(stderr, "Error, missing -log command line argument\n");
    usage (pl, argv0);
  }
  if (relspfilename == NULL)
  {
    fprintf(stderr, "Error, missing -purged command line argument\n");
    usage (pl, argv0);
  }
  if (relsdfilename == NULL)
  {
    fprintf(stderr, "Error, missing -relsdel command line argument\n");
    usage (pl, argv0);
  }
  if (outfilename == NULL)
  {
    fprintf(stderr, "Error, missing -out command line argument\n");
    usage (pl, argv0);
  }
  if (renumberfilename == NULL)
  {
    fprintf(stderr, "Error, missing -renumber command line argument\n");
    usage (pl, argv0);
  }
  if (polyfilename == NULL)
  {
    fprintf(stderr, "Error, missing -poly command line argument\n");
    usage (pl, argv0);
  }
  if (mt < 1)
  {
    fprintf(stderr, "Error: parameter mt must be at least 1\n");
    usage (pl, argv0);
  }

  if (logformat != NULL)
  {
    if (strcmp(logformat, "LA") != 0 && strcmp(logformat, "reconstruct") != 0)
    {
      fprintf(stderr, "Error, unknown -formatlog argument. Must be 'LA' or "
                      "'reconstruct'\n");
      usage (pl, argv0);
    }
  }
  if ((logformat == NULL || strcmp(logformat, "LA") == 0) &&
      idealsfilename == NULL)
  {
    fprintf(stderr, "Error, missing -ideals command line argument\n");
    usage (pl, argv0);
  }

  if (wantedfilename != NULL && !partial)
  {
    fprintf(stderr, "Warning, -wanted command line argument is ignored if "
                    "-partial is not set\n");
  }

  cado_poly_init (poly);
#ifndef FOR_FFS
  if (!cado_poly_read (poly, polyfilename))
#else
  if (!ffs_poly_read (poly, polyfilename))
#endif
  {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

#ifndef FOR_FFS
  /* Get mpz_poly_t F from cado_poly pol (algebraic side) */
  F1 = poly->pols[ALGEBRAIC_SIDE];
# ifndef FOR_GFPN // humf, this is quite temporary...!
  FATAL_ERROR_CHECK(nbsm1 > (unsigned int) poly->alg->deg, "Too many SM");
# else
  F0 = poly->pols[RATIONAL_SIDE];
  if (abunitsfilename == NULL)
  {
    fprintf(stderr, "Error, missing -abunits command line argument\n");
    usage (pl, argv0);
  }
# endif // FOR_GFPN
#else
  FATAL_ERROR_CHECK(nbsm1 != 0, "sm should be 0 for FFS");
#endif

  set_antebuffer_path (argv0, path_antebuffer);

  /* Reading renumber file */
  printf ("\n###### Reading renumber file ######\n");
  renumber_init_for_reading (renumber_table);
  renumber_read_table (renumber_table, renumberfilename);
  nprimes = renumber_table->size;

  /* Read number of rows and cols on first line of purged file */
  purgedfile_read_firstline (relspfilename, &nrels_purged, &nideals_purged);
  nrels_del = nrels_tot - nrels_purged;

  /* Malloc'ing log tab and reading values of log */
  printf ("\n###### Reading known logarithms ######\n");
  fflush(stdout);
#ifndef FOR_GFPN
  logtab_init (log, nprimes, nbsm1, q);
#else
  logtab_init (log, nprimes, nbsm1+nbsm0, q);
#endif
  if (logformat == NULL || strcmp(logformat, "LA") == 0)
    read_log_format_LA (log, logfilename, idealsfilename);
  else
    read_log_format_reconstruct (log, renumber_table, logfilename);

  /* Init bit_vector of rels that must be process by compute_log_from_rels */
  bit_vector rels_to_process;
  bit_vector_init (rels_to_process, nrels_tot);

  if (partial)
  {
    if (wantedfilename == NULL)
    {
      bit_vector_set (rels_to_process, 0);
      nrels_needed = 0;
    }
    else /* We compute needed rels for logarithms in wantedfilename */
    {
      printf ("\n###### Computing needed rels ######\n");
      nrels_needed =
        compute_needed_rels (rels_to_process, relspfilename, nrels_purged,
                             relsdfilename, nrels_del, log, wantedfilename, mt);
    }
  }
  else
  {
    bit_vector_set (rels_to_process, 1);
    nrels_needed = nrels_tot;
  }

#ifndef FOR_FFS
  sm_data_init(q, log); /* Init data needed to compute SM */
#endif

  /* Computing logs using rels in purged file */
  printf ("\n###### Computing logarithms using rels ######\n");
  if (nrels_needed > 0)
  {
    log->nknown += compute_log_from_rels (rels_to_process, relspfilename,
                                          nrels_purged, relsdfilename,
                                          nrels_del, nrels_needed, log, mt,
					  abunitsfilename);
    printf ("# %" PRIu64 " logarithms are known.\n", log->nknown);
  }
  else
    printf ("# All wanted logarithms are already known, skipping this step\n");
  fflush(stdout);

  /* Writing all the logs in outfile */
  printf ("\n###### Writing logarithms in a file ######\n");
  write_log (outfilename, log, renumber_table, poly);

  /* freeing and closing */
#ifndef FOR_FFS
  sm_data_free();
  mpz_clear(smexp1);
# ifdef FOR_GFPM
  mpz_clear(smexp0);
# endif
#endif
  logtab_free (log);
  mpz_clear(q);

  renumber_clear (renumber_table);
  bit_vector_clear(rels_to_process);
  cado_poly_clear (poly);
  param_list_clear (pl);
  return EXIT_SUCCESS;
}
