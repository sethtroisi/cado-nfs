#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "filter_utils.h"

#ifdef FOR_FFS
#include "fppol.h"
#include "fq.h"
#include "utils_ffs.h"
#endif
#define DEBUG 0

static ideal_merge_t **rel_purged;
static weight_t *ideals_weight;
mpz_t vlog, invert_coeff;

void
read_matrix_indexing (index_t *tab, FILE *file, index_t ncols, index_t nprimes)
{
  index_t nbread = 0;
  index_t index;
  index_t j;

  while (fscanf (file, "%u %x\n", &index, &j) == 2)
  {
    ASSERT_ALWAYS (index < ncols);
    ASSERT_ALWAYS (j < nprimes);
    tab[index] = j;
    nbread++;
  }

  ASSERT_ALWAYS (nbread == ncols);
}

void
read_log (mpz_t *log, index_t *mat_renum, FILE *logfile, mpz_t q, index_t ncols)
{
  index_t nbread = 0;
  index_t i;
  mpz_t vlog;

  mpz_init (vlog);

  while (gmp_fscanf (logfile, "%"PRid" %Zd\n", &i, vlog) == 2)
  {
    if (mpz_cmp_ui (vlog, 0) < 0)
    {
      fprintf (stderr, "  Warning, log is negative for cols %" PRid "\n", i);
      mpz_mod (vlog, vlog, q);
    }
    else if (mpz_cmp_ui (vlog, 0) == 0)
      fprintf (stderr, "  Warning, log is zero for cols %" PRid "\n", i);
    else if (mpz_cmp (vlog, q) >= 0)
    {
      fprintf (stderr, "  Warning, log >= q for cols %" PRid "\n", i);
      mpz_mod (vlog, vlog, q);
    }

    mpz_set (log[mat_renum[i]], vlog);
    nbread++;
  }

  ASSERT_ALWAYS (nbread == ncols);

  mpz_clear (vlog);
}

static int 
usable (index_t i, mpz_t *log)
{
  unsigned int count = 0;
  ideal_merge_t *p = rel_purged[i];

#if DEBUG >= 2
  fprintf (stderr, "    rel i=%u unknown logarithms: ", i);
#endif

  for (; p->id != UMAX(p->id); p++)
  {
    if (mpz_cmp_si(log[p->id], -1) <= 0) // we do not know the log if this ideal
    {
#if DEBUG >= 2
  fprintf (stderr, "%" PRxid " ", p->id);
#endif
        count++;
    }
  }

#if DEBUG >= 2
  fprintf (stderr, "(count=%u)\n", count);
#endif
  return (count <= 1);
}

// return the number of computed log (0 or 1)
static unsigned int
compute_log (index_t i, mpz_t *log, mpz_t q)
{
  int32_t coeff = 0;
  index_t unknown_ideal, h;
  ideal_merge_t *p = rel_purged[i];

  mpz_set_ui (vlog, 0); 

  for (; p->id != UMAX(p->id); p++)
  {
    h = p->id;
    if (mpz_cmp_si(log[h], -1) <= 0) // we do not know the log if this ideal
    {
      if (coeff == 0) //first time we see an unknown ideal
      {
        coeff = p->e;
        unknown_ideal = h;
      }
      else
      {
        fprintf (stderr, "Error, too much unknown ideals in relation %u\n", i);
        exit (1);
      }
    }
    else
      mpz_addmul_ui (vlog, log[h], (unsigned int) p->e);
  }

  if (coeff == 0) //no unknown ideal, sum should be zero
  {
    mpz_mod (vlog, vlog, q);
    if (mpz_cmp_ui (vlog, 0) != 0)
    {
      gmp_fprintf (stderr, "    No unknow log in rel %u and sum of log is not "
                           "zero (sum is %Zd), error!\n", i, vlog);
      exit (1);
    }
#if DEBUG >= 2
    else
    {
      fprintf (stderr, "    No unknow log in rel %u but sum of log is zero,"
                       " continue.\n", i);
    }
#endif
    return 0;
  }
  else
  { 
    mpz_set_si (invert_coeff, coeff);
    mpz_invert (invert_coeff, invert_coeff, q);
    mpz_neg (vlog, vlog);
    mpz_mul (vlog, vlog, invert_coeff);
    mpz_mod (log[unknown_ideal], vlog, q);
#if DEBUG >= 2
    fprintf (stderr, "    New logarithm computed for index %" PRxid ".\n",
                                                                unknown_ideal);
#endif
    return 1;
  }
}

#if DEBUG == 1
static void
check_unused_rel (bit_vector not_used, mpz_t *log, index_t nrels)
{
  unsigned int count;
  ideal_merge_t *p;

  for (index_t i = 0; i < nrels; i++)
  {
    if (bit_vector_getbit(not_used, (size_t) i))
    {
      count = 0;
      p = rel_purged[i];
      fprintf (stderr, "  rel i=%u unknown logarithm: ", i);
      for (; p->id != UMAX(p->id); p++)
        if (mpz_cmp_si(log[p->id], -1) <= 0)
        {
          fprintf (stderr, "%" PRxid " ", p->id);
          count++;
        }
      
      fprintf (stderr, "(count=%u)\n", count);
    }
  }
}
#endif

static inline unsigned int
count_bits (uint64_t n)
{
  unsigned int c = 0;
  for (; n ; c++)
    n &= (n-1);

  return c;
}

static index_t
compute_missing_log (mpz_t *log, mpz_t q, index_t nrels)
{
  index_t i, change, computed = 0, total_computed = 0;
  unsigned int iter = 1;
  double tt;
  bit_vector not_used;
  mpz_init (vlog);
  mpz_init (invert_coeff);

  bit_vector_init(not_used, nrels);
  ASSERT_ALWAYS (not_used->p != NULL);
  bit_vector_set(not_used, 1);
  if (nrels & (BV_BITS - 1))
    not_used->p[nrels>>LN2_BV_BITS] &= (((bv_t) 1)<<(nrels & (BV_BITS - 1))) - 1;

  do
  {
    computed = 0;
    change = 0;
    tt = seconds();
    fprintf (stderr, "  Iteration %u: begin\n", iter);
    for (i = 0; i < nrels; i++)
    {
      if (bit_vector_getbit(not_used, (size_t) i) && usable(i, log))
      {
        computed += compute_log(i, log, q);
        bit_vector_clearbit(not_used, (size_t) i);
      }
      if (i >> 18 != change >> 18)
      {
        fprintf(stderr, "  Iteration %u: %" PRid " lines read, %" PRid " new "
                        "logarithms computed\n", iter, i, computed);
        change = i;
      }
    }
    total_computed += computed;
    fprintf (stderr, "  Iteration %u: end with %" PRid " new logarithms "
                     "computed in %.1fs.\n", iter, computed, seconds()-tt);
    iter++;
  } while (computed);

  index_t c = 0; 
  size_t j;
  for (j = 0; j < iceildiv(not_used->n, BV_BITS); j++)
    c += count_bits(not_used->p[j]);
  if (c != 0)
    fprintf(stderr, "Warning, %" PRid " relations were not used\n", c);
#if DEBUG == 1
  if (c != 0)
    check_unused_rel(not_used, log, nrels);
#endif

  bit_vector_clear(not_used);
  mpz_clear (vlog);
  mpz_clear (invert_coeff);
  return total_computed;
}

static index_t
write_log (FILE *outfile, mpz_t *log, mpz_t q, renumber_t tab, cado_poly poly)
{
  index_t i, missing = 0;

  for (i = 0; i < tab->size; i++)
  {
    if (mpz_cmp_si (log[i], -1) <= 0)
    {
#if DEBUG >= 2
      gmp_fprintf (stderr, "log[%" PRxid "] is negative: %Zd\n", i, log[i]);
#endif
      missing++;
    }
    else if (tab->table[i] == RENUMBER_SPECIAL_VALUE)
    {
      ASSERT_ALWAYS (mpz_cmp (log[i], q) < 0);
      if (i == 0 && tab->add_full_col)
        gmp_fprintf (outfile, "%" PRid " added column %Zd\n", i, log[i]);
      else
        gmp_fprintf (outfile, "%" PRid " bad ideals %Zd\n", i, log[i]);
    }
    else
    {
      ASSERT_ALWAYS (mpz_cmp (log[i], q) < 0);
      p_r_values_t p, r;
      int side;
      renumber_get_p_r_from_index (tab, &p, &r, &side, i, poly);
      if (side != tab->rat)
        gmp_fprintf (outfile, "%" PRid " %" PRpr " %d %" PRpr " %Zd\n", i, p,
                              side, r, log[i]);
      else
        gmp_fprintf (outfile, "%" PRid " %" PRpr " %d rat %Zd\n", i, p, side,
                              log[i]);
    }
  }

  return missing;
}

#if DEBUG >= 1
// if weight[k] != 0 => k appear in a rel => log should be known
static void
check_unknown_log (mpz_t *log, index_t nprimes)
{
  index_t c = 0;
  for (index_t k = 0; k < nprimes; k++)
    if (ideals_weight[k])
      if (mpz_cmp_si(log[k], -1) <= 0)
      {
        c++;
        fprintf (stderr, "ideal %" PRxid " appear in rels (with weight %u) "
                         "but its logarithm is unknown\n", k, ideals_weight[k]);
      }

  fprintf (stderr, "%" PRid " more logarithms of ideals should be known\n", c);
}
#endif

/* Callback function called by preempt_scan_relations */
void *
thread_insert (buf_arg_t *arg)
{
  unsigned int j;
  unsigned long cpy_cpt_rel_b;
  buf_rel_t *my_rel;

  cpy_cpt_rel_b = cpt_rel_b;
  for ( ; ; )
  {
    while (cpt_rel_a == cpy_cpt_rel_b)
    {
      if (!is_finish())
        NANOSLEEP();
      else if (cpt_rel_a == cpy_cpt_rel_b)
        pthread_exit(NULL);
    }

    if (cpt_rel_a == cpy_cpt_rel_b + 1)
      NANOSLEEP();

    j = (unsigned int) (cpy_cpt_rel_b & (SIZE_BUF_REL - 1));
    my_rel = &(arg->rels[j]);

    arg->info.nprimes += insert_rel_in_table_with_e (my_rel, 0, 0, rel_purged,
                                                           ideals_weight);

    test_and_print_progress_now ();
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}

static void
usage (const char *argv0)
{
  fprintf (stderr, "Usage: %s [options]\n", argv0);
  fprintf (stderr, "# Mandatory command-line options:\n");
  fprintf (stderr, "    -log <f>         file containing known logarithms\n");
  fprintf (stderr, "    -ideals <f>      link between matrix cols and ideals\n");
  fprintf (stderr, "    -relsdel <f>     file with rels deleted by purge\n");
  fprintf (stderr, "    -relspurged <f>  file with rels that survived purge\n");
  fprintf (stderr, "    -out <f>         ouput file for logarithms\n");
  fprintf (stderr, "    -renumber <f>    renumbering table\n");
  fprintf (stderr, "    -poly <f>        poly file\n");
  fprintf (stderr, "    -nrels <n>       number of rels (same as purge)\n");
  fprintf (stderr, "    -q <n>           computations are done modulo q\n");
  exit (1);
}

int
main(int argc, char *argv[])
{
  char *argv0 = argv[0];

  uint64_t ncols_matrix = 0; /* number of column in the matrix */
  renumber_t renumber_table;
  index_t nrels_purged = 0, nrels_del;
  uint64_t nrels_tot = 0;
  index_t *matrix_indexing = NULL;

  int ret;
  FILE *logfile = NULL, *outfile = NULL, *idealsreplayfile = NULL;
  mpz_t *log = NULL;
  mpz_t q;
  double tt = 0.0;
  index_t i, ncomputed;
  index_t count_missing_log, known_log;
  cado_poly poly;

  param_list pl;
  param_list_init(pl);
  argv++,argc--;

  if (argc == 0)
    usage(argv0);

  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
    fprintf (stderr, "Unknown option: %s\n", argv[0]);
    usage(argv0);
  }

  /* Update parameter list at least once to register argc/argv pointers. */
  param_list_update_cmdline (pl, &argc, &argv);
  /* print command-line arguments */
  param_list_print_command_line (stdout, pl);
  fflush(stdout);

  const char * logfilename = param_list_lookup_string(pl, "log");
  const char * idealsreplayfilename = param_list_lookup_string(pl, "ideals");
  const char * relsdfilename = param_list_lookup_string(pl, "relsdel");
  const char * relspfilename = param_list_lookup_string(pl, "relspurged");
  const char * outfilename = param_list_lookup_string(pl, "out");
  const char * renumberfilename = param_list_lookup_string(pl, "renumber");
  const char * polyfilename = param_list_lookup_string(pl, "poly");
  param_list_parse_uint64(pl, "nrels", &nrels_tot);

  mpz_init (q);
  param_list_parse_mpz(pl, "q", q);
  if (mpz_cmp_ui (q, 0) <= 0)
  {
    fprintf (stderr, "Error, -q should be positive.");
    usage (argv0);
  }

  if (nrels_tot == 0)
  {
    fprintf (stderr, "Error, missing -nrels ... option (or nrels=0)\n");
    usage (argv0);
  }
  /* If nrels_tot > 2^32, then we need index_t to be 64-bit */
  if (((nrels_tot >> 32) != 0) && sizeof(index_t) < 8)
  {
    fprintf (stderr, "Error, -nrels is too large for a 32-bit program\n"
                     "See #define index_size in typedefs.h\n");
    exit(1);
  }
  if (logfilename == NULL)
  {
    fprintf (stderr, "Error, -log is missing.\n");
    usage (argv0);
  }
  if (idealsreplayfilename == NULL)
  {
    fprintf (stderr, "Error, -ideals is missing.\n");
    usage (argv0);
  }
  if (relspfilename == NULL)
  {
    fprintf (stderr, "Error, -relspurged is missing.\n");
    usage (argv0);
  }
  if (relsdfilename == NULL)
  {
    fprintf (stderr, "Error, -relsdel is missing.\n");
    usage (argv0);
  }
  if (outfilename == NULL)
  {
    fprintf (stderr, "Error, -out is missing.\n");
    usage (argv0);
  }
  if (renumberfilename == NULL)
  {
    fprintf (stderr, "Error, -renumber is missing.\n");
    usage (argv0);
  }
  if (polyfilename == NULL)
  {
    fprintf (stderr, "Error, -poly is missing.\n");
    usage (argv0);
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

  if (param_list_warn_unused (pl))
    usage (argv0);

  /* Opening out file */
  outfile = fopen_maybe_compressed (outfilename, "w");
  ASSERT_ALWAYS (outfile != NULL);

  /* Opening rels.purged file. Reading number of rels on first line */
  purgedfile_stream ps;
  purgedfile_stream_init(ps);
  purgedfile_stream_openfile (ps, relspfilename);
  nrels_purged = ps->nrows;
  nrels_del = nrels_tot - nrels_purged;
  purgedfile_stream_closefile (ps);

  /* Reading renumber file */
  renumber_init (renumber_table, poly);
  renumber_read_table (renumber_table, renumberfilename);

  /* Opening logfile (produced by linalg). Reading ncols_matrix on first line */
  logfile = fopen_maybe_compressed (logfilename, "r");
  ASSERT_ALWAYS (logfile != NULL);
  ret = fscanf (logfile, "%" SCNu64 "\n", &ncols_matrix);
  ASSERT_ALWAYS (ret == 1);

  /* Malloc'ing matrix_indexing and log mpz table */
  matrix_indexing = (index_t *) malloc (ncols_matrix * sizeof (index_t));
  ASSERT_ALWAYS (matrix_indexing != NULL);

  log = (mpz_t *) malloc (renumber_table->size * sizeof(mpz_t));
  for (i = 0; i < renumber_table->size; i++)
    mpz_init_set_si(log[i], -1);

  /* Opening ideals file, malloc'ing matrix_indexing and reading ideals file */
  fprintf (stderr, "Reading ideals file from %s\n", idealsreplayfilename);
  tt = seconds();
  idealsreplayfile = fopen_maybe_compressed (idealsreplayfilename, "r");
  ASSERT_ALWAYS (idealsreplayfile != NULL);
  read_matrix_indexing (matrix_indexing, idealsreplayfile, ncols_matrix,
                                                           renumber_table->size);
  fprintf (stderr, "Reading index file took %.1fs\n", seconds() - tt);

  /* Reading log file */
  fprintf (stderr, "Reading logarithms computed by LA from %s\n", logfilename);
  tt = seconds();
  read_log (log, matrix_indexing, logfile, q, ncols_matrix);
  fprintf (stderr, "Reading %" PRIu64 " logs took %.1fs\n", ncols_matrix,
                                                                seconds() - tt);
  known_log = ncols_matrix;


  /*********** for rels.purged file ***********************************/

  /* Malloc'ing rel_purged and ideals_weight */
  rel_purged = (ideal_merge_t **) malloc (nrels_purged * sizeof(ideal_merge_t*));
  ASSERT_ALWAYS (rel_purged != NULL);

  ideals_weight = (weight_t *) malloc (renumber_table->size * sizeof(weight_t));
  ASSERT_ALWAYS (ideals_weight != NULL);
  memset(ideals_weight, 0, renumber_table->size * sizeof(weight_t));

  /* Reading all relations */
  fprintf (stderr, "Reading relations from %s\n", relspfilename);
  set_antebuffer_path (argv0, NULL);
  char *fic[2] = {(char *) relspfilename, NULL};
  process_rels (fic, &thread_insert, NULL, 0, NULL, NULL, STEP_RECONSTRUCT);

  /* computing missing log */
  fprintf (stderr, "Computing missing logarithms from rels left "
                   "after purge...\n");
  tt = seconds();
  ncomputed = compute_missing_log (log, q, nrels_purged);
  known_log += ncomputed;
  fprintf (stderr, "Computing %"PRid" new logarithms took %.1fs "
                   "(%" PRid " known logarithms so far)\n", ncomputed, 
                   seconds() - tt, known_log);

#if DEBUG >= 1
  check_unknown_log (log, renumber_table->size);
#endif

  /*********** for rels.deleted file ***********************************/
  /* Malloc'ing rel_purged and memset ideals_weight to zero */
  free(rel_purged);
  rel_purged = (ideal_merge_t **) malloc (nrels_del * sizeof(ideal_merge_t*));
  ASSERT_ALWAYS (rel_purged != NULL);
  memset(ideals_weight, 0, renumber_table->size * sizeof(weight_t));

  /* Reading all relations */
  fprintf (stderr, "Reading relations from %s\n", relsdfilename);
  char *fic2[2] = {(char *) relsdfilename, NULL};
  process_rels (fic2, &thread_insert, NULL, 0, NULL, NULL, STEP_RECONSTRUCT);

  /* computing missing log */
  fprintf (stderr, "Computing missing logarithms from all rels...\n");
  tt = seconds();
  ncomputed = compute_missing_log (log, q, nrels_del);
  known_log += ncomputed;
  fprintf (stderr, "Computing %"PRid" new logarithms took %.1fs "
                   "(%" PRid " known logarithms so far)\n", ncomputed, 
                   seconds() - tt, known_log);

#if DEBUG >= 1
  check_unknown_log (log, renumber_table->size);
#endif

  /*********************** end *****************************************/
  /* writing all the logs in outfile */
  fprintf (stderr, "Writing logarithms in %s\n", outfilename);
  tt = seconds();
  count_missing_log = write_log (outfile, log, q, renumber_table, poly);
  fprintf (stderr, "  %" PRid " logarithms are known, %" PRid " are missing\n",
                   known_log, count_missing_log);
  fprintf (stderr, "Writing logarithms took %.1fs\n", seconds()-tt);
  ASSERT_ALWAYS (known_log + count_missing_log == renumber_table->size);

  /* freeing and closing */
  fclose_maybe_compressed (outfile, outfilename);
  fclose_maybe_compressed (logfile, logfilename);
  fclose_maybe_compressed (idealsreplayfile, idealsreplayfilename);

  free (matrix_indexing);
  free (rel_purged);
  free (ideals_weight);

  for (i = 0; i < renumber_table->size; i++)
    mpz_clear(log[i]);
  free(log);
  mpz_clear(q);

  renumber_free (renumber_table);
  cado_poly_clear (poly);
  return 0;
}
