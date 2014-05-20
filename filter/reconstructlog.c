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

unsigned int nbsm = 0; /* number of SM that should be used. Must be 0 for FFS */


/* used as mutual exclusion lock for the array of log */
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

/* 2 functions to compute a <- (a + l*e) mod q (e is either int32_t or mpz_t) */
static inline void
mpz_add_log_mod_si (mpz_t a, mpz_t l, int32_t e, mpz_t q)
{
  if (e > 0)
    mpz_addmul_ui (a, l, e);
  else
    mpz_submul_ui (a, l, -e);
  mpz_mod (a, a, q);
}

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

/* Structure for a relation */
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

/* Structure containing the data for the reading of a relation file */
typedef struct
{
  log_rel_t *rels;
  int32_t *w;
  mpz_t *log;
  mpz_ptr q;
  uint64_t nideals;

} read_data_t;

/* Init a read_data_t structure for nrels rels and nprimes primes. Assume the
 * table of log is already allocated, but not the table of log_rel_t     */
static void
read_data_init (read_data_t *data, mpz_t *log, mpz_t q, uint64_t nprimes,
                uint64_t nrels)
{
  data->w = (int32_t *) malloc (nprimes * sizeof(int32_t));
  FATAL_ERROR_CHECK (data->w == NULL, "Cannot allocate memory");
  memset(data->w, 0, nprimes * sizeof(int32_t));
  data->rels = log_rel_init(nrels);
  data->log = log;
  data->q = q;
  data->nideals = 0;
}

static void
read_data_free (read_data_t *data, uint64_t nrels)
{
  free(data->w);
  log_rel_free (data->rels, nrels);
}

#ifndef FOR_FFS /* Not needed for FFS */
mpz_t smexp; /* exponent for SM */
mpz_poly_ptr F;
mpz_t q2;    /* q^2 */
mpz_t invq2;
mpz_t *smlog;

/* Init the data needed for SM computation, given the polynomial,
 * the group order and the exponent */
static void
sm_data_init (mpz_t q, mpz_t *log, uint64_t nprimes)
{
  mpz_init(q2);
  mpz_mul(q2, q, q);

  mpz_init(invq2);
  barrett_init(invq2, q2);

  smlog = &(log[nprimes]);
}

static void
sm_data_free ()
{
  mpz_clear(q2);
  mpz_clear(invq2);
}

/* Given a and b, compute the SM and add the contribution to l */
static inline void
add_sm_contribution (mpz_ptr l, int64_t a, uint64_t b, mpz_t q)
{
  mpz_poly_t SMres;
  int degF = F->deg;
  mpz_poly_init(SMres, degF);
  SMres->deg = 0;
  mpz_poly_setcoeff_si(SMres, 0, 1);
  sm_single_rel(SMres, a, b, F, smexp, q, q2, invq2);
  unsigned int i;
  for (i = degF - SMres->deg - 1; i < nbsm; i++)
    mpz_add_log_mod_mpz (l, smlog[i], SMres->coeff[degF-1-i], q);
  mpz_poly_clear(SMres);
}

/* Callback function called by filter_rels in compute_log_from_relfile */
void *
thread_sm (void * context_data, earlyparsed_relation_ptr rel)
{
  read_data_t *data = (read_data_t *) context_data;
  log_rel_t *lrel = &(data->rels[rel->num]);

  if (nbsm > 0)
    add_sm_contribution (lrel->log_known_part, rel->a, rel->b, data->q);

  return NULL;
}
#endif /* ifndef FOR_FFS */

/* Read the index file produced by replay (needed to read the logarithms)*/
static index_t *
read_index (const char *filename, uint64_t *ncols, uint64_t nprimes)
{
  FILE *f = NULL;
  uint64_t nbread = 0, i, j;
  index_t *tab = NULL;

  double tt = seconds();
  printf ("# Reading ideals file from %s\n", filename);
  fflush(stdout);
  f = fopen_maybe_compressed (filename, "r");
  FATAL_ERROR_CHECK(f == NULL, "Cannot open file for reading");

  if (fscanf (f, "# %" SCNu64 "\n", ncols) != 1)
  {
    fprintf(stderr, "Error while reading first line of %s\n", filename);
    exit(EXIT_FAILURE);
  }

  tab = (index_t *) malloc (*ncols * sizeof (index_t));
  FATAL_ERROR_CHECK(tab == NULL, "Cannot allocate memory");

  while (fscanf (f, "%" SCNu64 " %" SCNx64 "\n", &i, &j) == 2)
  {
    FATAL_ERROR_CHECK (i >= *ncols, "Too big value of column number");
    FATAL_ERROR_CHECK (j >= nprimes, "Too big value of index");
    tab[i] = j;
    nbread++;
  }

  FATAL_ERROR_CHECK (nbread != *ncols, "Not enough or too many index read");
  fclose_maybe_compressed (f, filename);
  printf ("# Reading %" PRIu64 " index took %.1fs\n", nbread, seconds() - tt);
  return tab;
}

/* Read the logarithms computed by the linear algebra */
static mpz_t *
read_log (index_t *mat_renum, const char *filename, mpz_t q, unsigned int nbsm,
          uint64_t ncols, uint64_t nprimes)
{
  uint64_t i;
  mpz_t tmp_log, *log = NULL;
  FILE *f = NULL;

  double tt = seconds();
  printf ("# Reading logarithms computed by LA from %s\n", filename);
  fflush(stdout);
  f = fopen_maybe_compressed (filename, "r");
  FATAL_ERROR_CHECK(f == NULL, "Cannot open file for reading");

  log = (mpz_t *) malloc ((nprimes + nbsm) * sizeof(mpz_t));
  FATAL_ERROR_CHECK(log == NULL, "Cannot allocate memory");
  size_t q_nbits = mpz_size(q) * GMP_LIMB_BITS;
  for (i = 0; i < nprimes + nbsm; i++)
  {
    mpz_init2 (log[i], q_nbits);
    mpz_set_si(log[i], -1);
  }

  mpz_init (tmp_log);
  i = 0;
  while (i < ncols + nbsm)
  {
    int ret = gmp_fscanf (f, "%Zd\n", tmp_log);
    FATAL_ERROR_CHECK (ret != 1, "Error in file containing logarithms values");

    if (mpz_cmp_ui (tmp_log, 0) < 0)
    {
      fprintf (stderr, "Warning, log is negative for cols %" PRIu64 "\n", i);
      mpz_mod (tmp_log, tmp_log, q);
    }
    else if (mpz_cmp (tmp_log, q) >= 0)
    {
      fprintf (stderr, "Warning, log >= q for cols %" PRIu64 "\n", i);
      mpz_mod (tmp_log, tmp_log, q);
    }
    if (mpz_cmp_ui (tmp_log, 0) == 0)
      fprintf (stderr, "Warning, log is zero for cols %" PRIu64 "\n", i);

    if (i < ncols)
      mpz_set (log[mat_renum[i]], tmp_log);
    else // log that corresponds to a SM columns in the matrix
      mpz_set (log[nprimes+(i-ncols)], tmp_log);
    i++;
  }

  while (gmp_fscanf (f, "%Zd\n", tmp_log) == 1)
    printf ("Warning, line %" PRIu64 " is ignored\n", i++);

  printf ("# Reading %" PRIu64 " logs took %.1fs\n", ncols, seconds() - tt);
  if (nbsm)
    printf ("# logs for %u SM columns were also read\n", nbsm);
  mpz_clear (tmp_log);
  fclose_maybe_compressed (f, filename);
  return log;
}

/* Write values of the known logarithms. Return the number of missing values */
static uint64_t
write_log (const char *filename, mpz_t *log, mpz_t q, renumber_t tab,
           cado_poly poly, uint64_t known_log)
{
  uint64_t i, missing = 0;
  double tt = seconds();
  FILE *f = NULL;

  printf ("# Opening %s for writing logarithms\n", filename);
  fflush(stdout);
  f = fopen_maybe_compressed (filename, "w");
  FATAL_ERROR_CHECK(f == NULL, "Cannot open file for writing");

  for (i = 0; i < tab->size; i++)
  {
	  if (mpz_sgn(log[i]) < 0) // we do not know the log if this ideal
      missing++;
    else if (tab->table[i] == RENUMBER_SPECIAL_VALUE)
    {
      ASSERT_ALWAYS (mpz_cmp (log[i], q) < 0);
      if (i == 0 && tab->add_full_col)
        gmp_fprintf (f, "%" PRid " added column %Zd\n", i, log[i]);
      else
        gmp_fprintf (f, "%" PRid " bad ideals %Zd\n", i, log[i]);
    }
    else
    {
      ASSERT_ALWAYS (mpz_cmp (log[i], q) < 0);
      p_r_values_t p, r;
      int side;
      renumber_get_p_r_from_index (tab, &p, &r, &side, i, poly);
      if (side != tab->rat)
        gmp_fprintf (f, "%" PRid " %" PRpr " %d %" PRpr " %Zd\n", i, p, side,
                                                                  r, log[i]);
      else
        gmp_fprintf (f, "%" PRid " %" PRpr " %d rat %Zd\n", i, p, side, log[i]);
    }
  }

  printf ("# Writing logarithms took %.1fs\n", seconds()-tt);
  printf ("# %" PRIu64 " logarithms of elements of the factor base are known, "
          "%" PRIu64 " are missing.\n", known_log, missing);
  fclose_maybe_compressed (f, filename);
  ASSERT_ALWAYS (known_log + missing == tab->size);
  return missing;
}

/* Callback function called by filter_rels in compute_log_from_relfile */
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

    if (data->w[h] == 0)
    {
      data->w[h] = 1;
      data->nideals++;
    }
    else if (data->w[h] != SMAX(int32_t))
      data->w[h]++;

	  if (mpz_sgn(data->log[h]) >= 0)
      mpz_add_log_mod_si (lrel->log_known_part, data->log[h], e, data->q);
    else
      buf[next++] = (ideal_merge_t) {.id = h, .e = e};
  }

  lrel->unknown = idealmerge_my_malloc (next);
  lrel->nb_unknown = next;
  memcpy(lrel->unknown, buf, next * sizeof(ideal_merge_t));

  return NULL;
}

/* Debug functions */
#if DEBUG == 1
/* Debug functions */ /* List all unused relations.*/
static void
check_unused_rel (bit_vector not_used, log_rel_t *rels, uint64_t nrels)
{
  for (uint64_t i = 0; i < nrels; i++)
  {
    if (bit_vector_getbit(not_used, (size_t) i))
    {
      fprintf (stderr, "DEBUG: rel i=%" PRIu64 ": %u unknown logarithms: ",
                       i, rels[i].nb_unknown);
      for (weight_t k = 0; k < rels[i].nb_unknown; k++)
        fprintf (stderr, "%" PRid " ", rels[i].unknown[k].id);
      fprintf (stderr, "\n");
    }
  }
}

/* If an ideal had a non zero weight, it logarithm should be known at the end of
 * compute_log_from_relfile */
static void
check_unknown_log (read_data_t *data, uint64_t nprimes)
{
  uint64_t c = 0;
  int32_t *weight = data->w;
  mpz_t *log = data->log;
  for (index_t k = 0; k < nprimes; k++)
    if (weight[k])
	    if (mpz_sgn(log[k]) < 0)
      {
        c++;
        fprintf (stderr, "DEBUG: ideal %" PRid " had weight %u but its "
                         "logarithm is unknown\n", k, weight[k]);
      }
  fprintf (stderr, "DEBUG: %" PRIu64 " more logarithms should be known\n", c);
}
#endif

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
    int sgn = mpz_sgn(data->log[p[k].id]);
    mutex_unlock (&lock);

	  if (sgn < 0) // we do not know the log if this ideal
    {
      if (j != k)
        p[j] = p[k];
      j++;
    }
    else // We know this log, add it to log_know_part
      mpz_add_log_mod_si (data->rels[i].log_known_part, data->log[p[k].id],
                          p[k].e, data->q);
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
do_one_part_of_iter (read_data_t *data, bit_vector not_used, uint64_t start,
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
          gmp_fprintf (stderr, "Error, no unknow log in rel %" PRIu64 " and sum"
                       " of log is not zero, sum is:\n%Zd\n", i, vlog);
          exit (EXIT_FAILURE);
        }
        else if (nb == 1)
        {
          ideal_merge_t ideal = data->rels[i].unknown[0];
          computed +=
             compute_missing_log (data->log[ideal.id], vlog, ideal.e, data->q);
        }
      }
    }
  }

  return computed;
}

/* Compute all missing logarithms possible. Run through all the relations once.
 * Mono thread version
 * Return the number of computed logarithms */
static uint64_t
do_one_iter_mono (read_data_t *data, bit_vector not_used, uint64_t nrels)
{
  return do_one_part_of_iter (data, not_used, 0, nrels);
}

/* Code for multi thread version */
#define SIZE_BLOCK 1024

typedef struct {
  read_data_t *data;
  bit_vector_ptr not_used;
  uint64_t offset;
  uint64_t nb;
  uint64_t computed;
} thread_info;

void * thread_start(void *arg)
{
  thread_info *ti = (thread_info *) arg;
  read_data_t *data = ti->data;
  bit_vector_ptr not_used = ti->not_used;
  uint64_t start = ti->offset;
  uint64_t end = start + ti->nb;

  ti->computed = do_one_part_of_iter (data, not_used, start, end);

  return NULL;
}

/* Compute all missing logarithms possible. Run through all the relations once.
 * Multi thread version
 * Return the number of computed logarithms */
static uint64_t
do_one_iter_mt (read_data_t *data, bit_vector not_used, int nt, uint64_t nrels)
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



/* Given a filename, compute all the possible logarithms of ideals appearing in
 * the file. Return the number of computed logarithms. */
static uint64_t
compute_log_from_relfile (const char *filename, uint64_t nrels, mpz_t q,
                          mpz_t *log, uint64_t nprimes, int nt)
{
  double tt0, tt, wct_tt0, wct_tt;
  uint64_t total_computed = 0, iter = 0, computed;
  read_data_t data;
  read_data_init(&data, log, q, nprimes, nrels);

  /* Init bit_vector to remember which relations were already used */
  bit_vector not_used;
  bit_vector_init(not_used, nrels);
  bit_vector_set(not_used, 1);

  /* adjust the number of threads based on the number of relations */
  double ntm = ceil((nrels + 0.0)/SIZE_BLOCK);
  if (nt > ntm)
    nt = (int) ntm;

  if (nt > 1)
    printf("# Using multithread version with %d threads\n", nt);
  else
    printf("# Using monothread version\n");

  /* Reading all relations */
  printf ("# Reading relations from %s\n", filename);
  fflush(stdout);
  char *fic[2] = {(char *) filename, NULL};
  struct filter_rels_description desc[3] = {
                   { .f = thread_insert, .arg=&data, .n=1},
#ifndef FOR_FFS
                   { .f = thread_sm,     .arg=&data, .n=nt},
#endif
                   { .f = NULL,          .arg=0,     .n=0}
      };
  filter_rels2 (fic, desc, EARLYPARSE_NEED_AB_HEXA | EARLYPARSE_NEED_INDEX,
                NULL, NULL);

  /* computing missing log */
  printf ("# Starting to computing missing logarithms from rels\n");
  tt0 = seconds();
  wct_tt0 = wct_seconds();
  do
  {
    printf ("# Iteration %" PRIu64 ": starting...\n", iter);
    fflush(stdout);
    tt = seconds();
    wct_tt = wct_seconds();

    if (nt > 1)
      computed = do_one_iter_mt (&data, not_used, nt, nrels);
    else
      computed = do_one_iter_mono (&data, not_used, nrels);
    total_computed += computed;

    printf ("# Iteration %" PRIu64 ": %" PRIu64 " new logarithms computed\n",
            iter, computed);
    if (nt > 1)
      printf ("# Iteration %" PRIu64 " took %.1fs (CPU time), %.1fs (wall-clock"
              " time).\n", iter, seconds() - tt, wct_seconds() - wct_tt);
    else
      printf ("# Iteration %" PRIu64 " took %.1fs\n", iter, seconds() - tt);

    iter++;
  } while (computed);

  if (nt > 1)
    printf ("# Computing %" PRIu64 " new logarithms took %.1fs (CPU time), "
            "%.1fs (wall-clock time)\n", total_computed, seconds() - tt0,
            wct_seconds() - wct_tt0);
  else
    printf ("# Computing %" PRIu64 " new logarithms took %.1fs\n",
            total_computed, seconds() - tt0);

  size_t c = bit_vector_popcount(not_used);
  if (c != 0)
    fprintf(stderr, "### Warning, %zu relations were not used\n", c);
#if DEBUG == 1
  if (c != 0)
    check_unused_rel(not_used, data->rels, nrels);
  check_unknown_log (&data, nprimes);
#endif

  read_data_free(&data, nrels);
  bit_vector_clear(not_used);
  return total_computed;
}


static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "log", "input file containing known logarithms");
  param_list_decl_usage(pl, "gorder", "group order (see sm -gorder parameter");
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
  param_list_decl_usage(pl, "smexp", "sm exponent (see sm -smexp parameter)");
#endif
  param_list_decl_usage(pl, "mt", "number of threads (default 1)");
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
  uint64_t ncols_matrix, nrels_tot = 0, nrels_purged, nideals_purged, nrels_del;
  uint64_t nprimes, i, known_log;
  index_t *matrix_indexing = NULL;
  int mt = 1;
  int partial = 0;

  mpz_t q, *log = NULL;
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
  const char * idealsfilename = param_list_lookup_string(pl, "ideals");
  const char * relsdfilename = param_list_lookup_string(pl, "relsdel");
  const char * relspfilename = param_list_lookup_string(pl, "purged");
  const char * outfilename = param_list_lookup_string(pl, "out");
  const char * renumberfilename = param_list_lookup_string(pl, "renumber");
  const char * polyfilename = param_list_lookup_string(pl, "poly");
  param_list_parse_uint64(pl, "nrels", &nrels_tot);
  param_list_parse_mpz(pl, "gorder", q);
#ifndef FOR_FFS
  param_list_parse_uint(pl, "sm", &nbsm);
  mpz_init (smexp);
  param_list_parse_mpz(pl, "smexp", smexp);
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
  if (nbsm != 0 && mpz_cmp_ui (smexp, 0) <= 0)
  {
    fprintf(stderr, "Error, missing -smexp command line argument "
                    "(or smexp <= 0)\n");
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
  if (idealsfilename == NULL)
  {
    fprintf(stderr, "Error, missing -ideals command line argument\n");
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
  F = poly->pols[ALGEBRAIC_SIDE];
  FATAL_ERROR_CHECK(nbsm > (unsigned int) poly->alg->deg, "Too many SM");
#else
  FATAL_ERROR_CHECK(nbsm != 0, "sm should be 0 for FFS");
#endif


  set_antebuffer_path (argv0, path_antebuffer);

  /* Reading renumber file */
  renumber_init (renumber_table, poly, NULL);
  renumber_read_table (renumber_table, renumberfilename);
  nprimes = renumber_table->size;

  /* Read number of rows and cols on first line of purged file */
  purgedfile_read_firstline (relspfilename, &nrels_purged, &nideals_purged);
  nrels_del = nrels_tot - nrels_purged;

  /* Opening ideals file, malloc'ing matrix_indexing and reading ideals file */
  matrix_indexing = read_index (idealsfilename, &ncols_matrix, nprimes);

  /* Opening log file, malloc'ing log tab and reading values of log */
  log = read_log (matrix_indexing, logfilename, q, nbsm, ncols_matrix, nprimes);
  known_log = ncols_matrix;
  free (matrix_indexing);

  if (!partial) {

    /* Init data needed to compute SM */
#ifndef FOR_FFS
    sm_data_init(q, log, nprimes);
#endif

    /* Computing logs using rels in purged file */
    printf ("\n###### Computing logs using rels in purged file ######\n");
    known_log +=
     compute_log_from_relfile (relspfilename, nrels_purged, q, log, nprimes, mt);
    printf ("# %" PRIu64 " known logarithms so far.\n", known_log);
    fflush(stdout);

    /* Computing logs using rels in deleted file */
    printf ("\n###### Computing log using rels in deleted file ######\n");
    known_log +=
      compute_log_from_relfile (relsdfilename, nrels_del, q, log, nprimes, mt);
    printf ("# %" PRIu64 " known logarithms at the end.\n\n", known_log);
    fflush(stdout);

#ifndef FOR_FFS
    sm_data_free();
#endif
  }

  /* Writing all the logs in outfile */
  write_log (outfilename, log, q, renumber_table, poly, known_log);

  /* freeing and closing */
  for (i = 0; i < nprimes + nbsm; i++)
    mpz_clear(log[i]);
  free(log);
  mpz_clear(q);
#ifndef FOR_FFS
  mpz_clear(smexp);
#endif

  renumber_free (renumber_table);
  cado_poly_clear (poly);
  param_list_clear (pl);
  return EXIT_SUCCESS;
}
