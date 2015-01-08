/* Shirokauer maps 
   
Input:

* A list of the (npurged) a,b pairs. This is obtained from the
  purgedfile.
* A matrix of (small_nrows) rows and (npurged) cols, which indicates
  the contents of each relation-set. This is obtained from the
  indexfile.
* The sub-group order (ell) such that ell | p-1
  Note: All computations are done mod ell^2.
* (eps): the exponent used in the computation of the Shirokauer maps.
  Note: eps = ppcm(eps_i), where eps_i = ell^(deg(f_i)) - 1 and f = f_1 ... f_k mod ell
  
Output

* A matrix of (small_nrows) rows and (nmaps)=deg(f) cols (mpz_t).  For each
  relation (rel) the (nmaps) Shirokauer maps are computed as the second
  least-significant digit of the ell-adic representation of the polynomial 
  equal to (rel^eps - 1) / ell.

  In case of two algebraic sides, SM's are computed for sides 0..1 in that
  order.

*/

#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pthread.h>
#include <errno.h>

#include "macros.h"
#include "filter_common.h"

stats_data_t stats; /* struct for printing progress */

void *
thread_sm (void * context_data, earlyparsed_relation_ptr rel)
{
  mpz_poly_t * abpolys = (mpz_poly_t *) context_data;
  mpz_poly_init_set_ab(abpolys[rel->num], rel->a, rel->b);

  return NULL;
}

sm_relset_ptr build_rel_sets(const char * purgedname, const char * indexname,
                             uint64_t * small_nrows, mpz_poly_ptr *F,
                             const mpz_t ell2)
{
  uint64_t nrows, ncols, small_ncols, len_relset;
  uint64_t r[MAX_LEN_RELSET];
  int64_t e[MAX_LEN_RELSET];
  int ret;

  /* array of (a,b) pairs from (purgedname) file */
  mpz_poly_t *pairs;

  purgedfile_read_firstline (purgedname, &nrows, &ncols);
  pairs = (mpz_poly_t *) malloc (nrows * sizeof(mpz_poly_t));
  ASSERT_ALWAYS (pairs != NULL);
  /* For each rel, read the a,b-pair and init the corresponding poly pairs[] */
  fprintf(stdout, "\n# Reading %" PRIu64 " (a,b) pairs\n", nrows);
  fflush(stdout);
  char *fic[2] = {(char *) purgedname, NULL};
  filter_rels (fic, (filter_rels_callback_t) thread_sm, pairs,
          EARLYPARSE_NEED_AB_HEXA, NULL, NULL);


  /* Array of (small_nrows) relation-sets built from array (pairs) and
     (indexname) file  */
  sm_relset_ptr rels;
  FILE * ix = fopen_maybe_compressed(indexname, "r");

  /* small_ncols isn't used here: we don't care */
  ret = fscanf(ix, "%" SCNu64 " %" SCNu64 "", small_nrows, &small_ncols);
  ASSERT(ret == 2);

  rels = (sm_relset_ptr) malloc (*small_nrows * sizeof(sm_relset_t));
  ASSERT_ALWAYS (rels != NULL);

  fprintf(stdout, "\n# Building %" PRIu64 " relation-sets\n", *small_nrows);
  fflush(stdout);
  uint64_t i;
  stats_init (stats, stdout, &i, nbits(*small_nrows)-5, "Computed",
              "relation-sets", "", "relsets");
  for(i = 0 ; i < *small_nrows ; i++)
  {
    ret = fscanf(ix, "%" SCNu64 "", &len_relset);
    ASSERT_ALWAYS(ret == 1 && len_relset < MAX_LEN_RELSET);

    for (uint64_t k = 0 ; k < len_relset ; k++)
    {
      ret = fscanf(ix, " %" SCNx64 ":%" SCNd64 "", &r[k], &e[k]); 
      ASSERT_ALWAYS(ret == 2);
    }
    
    int dF[2] = {0, 0};
    for (int s = 0; s < 2; ++s) {
        if (F[s] == NULL) continue;
        dF[s] = F[s]->deg;
    }
    sm_relset_init (&rels[i], dF);
    sm_build_one_relset (&rels[i], r, e, len_relset, pairs, F, ell2);

    if (stats_test_progress(stats))
      stats_print_progress (stats, i, 0, 0, 0);
  }
  stats_print_progress (stats, *small_nrows, 0, 0, 1);
  fclose_maybe_compressed(ix, indexname);

  for (uint64_t i = 0; i < nrows; i++)
    mpz_poly_clear (pairs[i]);
  free (pairs);
  
  return rels;
}

struct thread_info {
  int offset;
  int nb;
  sm_relset_ptr rels;
  mpz_poly_ptr *F;
  mpz_ptr *eps;
  mpz_srcptr ell;
  mpz_srcptr ell2;
  mpz_srcptr invl2;
  mpz_poly_t *sm0;
  mpz_poly_t *sm1;
  int *nsm;
};

void * thread_start(void *arg) {
  struct thread_info *ti = (struct thread_info *) arg;
  sm_relset_ptr rels = ti->rels;
  mpz_poly_ptr *F = ti->F;
  mpz_ptr *eps = ti->eps;
  mpz_srcptr ell = ti->ell;
  mpz_srcptr ell2 = ti->ell2;
  mpz_srcptr invl2 = ti->invl2;
  mpz_poly_t *sm0 = ti->sm0;
  mpz_poly_t *sm1 = ti->sm1;
  int *nsm = ti->nsm;
  int offset = ti->offset;

  for (int i = 0; i < ti->nb; i++) {
    if (nsm[0] > 0) {
      mpz_poly_reduce_frac_mod_f_mod_mpz(rels[offset+i].num[0],
              rels[offset+i].denom[0], F[0], ell2);
      compute_sm (sm0[i], rels[offset+i].num[0], F[0], ell, eps[0], ell2, invl2);
    }
    if (nsm[1] > 0) {
      mpz_poly_reduce_frac_mod_f_mod_mpz(rels[offset+i].num[1],
              rels[offset+i].denom[1], F[1], ell2);
      compute_sm (sm1[i], rels[offset+i].num[1], F[1], ell, eps[1], ell2, invl2);
    }
  }

  return NULL;
}

#define SM_BLOCK 512

void mt_sm (int nt, const char * outname, sm_relset_ptr rels, uint64_t sr,
            mpz_poly_ptr *F, mpz_ptr *eps,
            const mpz_t ell, const mpz_t ell2,
            int *nsm)
{
  // allocate space for results of threads
  mpz_poly_t **SM[2];
  for (int side = 0; side < 2; side++) {
      SM[side] = (mpz_poly_t **) malloc(nt*sizeof(mpz_poly_t *));
      for (int i = 0; i < nt; ++i) {
          SM[side][i] = (mpz_poly_t *) malloc(SM_BLOCK*sizeof(mpz_poly_t));
          for (int j = 0; j < SM_BLOCK; ++j)
              if (F[side] != 0)
                  mpz_poly_init(SM[side][i][j], F[side]->deg);
      }
  }

  // We'll use a rotating buffer of thread id.
  pthread_t *threads;
  threads = (pthread_t *) malloc(nt*sizeof(pthread_t));
  int active_threads = 0;  // number of running threads
  int threads_head = 0;    // next thread to wait / restart.
  
  // Prepare the main loop
  uint64_t i = 0; // counter of relation-sets.
  uint64_t out_cpt = 0; // counter of already printed relation-sets;
  FILE * out = fopen(outname, "w");
  int nsm_total=0;
  for (int side = 0; side < 2; side++) {
      nsm_total += nsm[side];
  }
  gmp_fprintf(out, "%" PRIu64 " %d %Zd\n", sr, nsm_total, ell);

  mpz_t invl2;
  mpz_init(invl2);
  barrett_init(invl2, ell2);

  // Arguments for threads
  struct thread_info *tis;
  tis = (struct thread_info*) malloc(nt*sizeof(struct thread_info));
  for (int i = 0; i < nt; ++i) {
    tis[i].rels = rels;
    tis[i].F = F;
    tis[i].eps = eps;
    tis[i].ell = ell;
    tis[i].ell2 = ell2;
    tis[i].invl2 = invl2;
    tis[i].sm0 = SM[0][i];
    tis[i].sm1 = SM[1][i];
    tis[i].nsm = nsm;
    // offset and nb must be adjusted.
  }

  // Main loop
  stats_init (stats, stdout, &out_cpt, nbits(sr)-5, "Computed", "SMs", "", "SMs");
  while ((i < sr) || (active_threads > 0)) {
    // Start / restart as many threads as allowed
    if ((active_threads < nt) && (i < sr)) { 
      tis[threads_head].offset = i;
      tis[threads_head].nb = MIN(SM_BLOCK, sr-i);
      pthread_create(&threads[threads_head], NULL, 
          &thread_start, (void *)(&tis[threads_head]));
      i += SM_BLOCK;
      active_threads++;
      threads_head++; 
      if (threads_head == nt) 
        threads_head = 0;
      continue;
    }
    // Wait for the next thread to finish in order to print result.
    pthread_join(threads[threads_head], NULL);
    active_threads--;
    for (int k = 0; k < SM_BLOCK && out_cpt < sr; ++k, ++out_cpt) {
      if (F[0] != NULL)
        print_sm (out, SM[0][threads_head][k], nsm[0], F[0]->deg);
      if (F[1] != NULL) {
        if (F[0] != NULL)
          fprintf(out, " ");
        print_sm (out, SM[1][threads_head][k], nsm[1], F[1]->deg);
      }
      fprintf(out, "\n");
    }

    // report
    if (stats_test_progress(stats))
      stats_print_progress (stats, out_cpt, 0, 0, 0);

    // If we are at the end, no job will be restarted, but head still
    // must be incremented.
    if (i >= sr) { 
      threads_head++;
      if (threads_head == nt) 
        threads_head = 0;
    }
  }
  stats_print_progress (stats, sr, 0, 0, 1);

  mpz_clear(invl2);
  fclose(out);
  free(tis);
  free(threads);
  for (int side = 0; side < 2; side++) {
      for (int i = 0; i < nt; ++i) {
          for (int j = 0; j < SM_BLOCK; ++j) {
              if (F[side] != NULL)
                  mpz_poly_clear(SM[side][i][j]);
          }
          free(SM[side][i]);
      }
      free(SM[side]);
  }
}


void sm (const char * outname, sm_relset_ptr rels, uint64_t sr,
        mpz_poly_ptr *F, mpz_ptr *eps,
        const mpz_t ell, const mpz_t ell2, int *nsm)
{
  FILE * out = fopen(outname, "w");
  DIE_ERRNO_DIAG(out==NULL, "fopen", outname);
  mpz_poly_t SM;
  mpz_t invl2;

  mpz_init(invl2);
  barrett_init(invl2, ell2);

  int dd = 0;
  if (F[0] != NULL)
      dd = MAX(dd, F[0]->deg);
  if (F[1] != NULL)
      dd = MAX(dd, F[1]->deg);
  mpz_poly_init(SM, dd);

  int nsm_total=0;
  for (int side = 0; side < 2; side++) {
      nsm_total += nsm[side];
  }
  gmp_fprintf(out, "%" PRIu64 " %d %Zd\n", sr, nsm_total, ell);

  uint64_t i;
  stats_init (stats, stdout, &i, nbits(sr)-5, "Computed", "SMs", "", "SMs");
  for (i = 0; i < sr; i++) {
    for (int side = 0; side < 2; side++) {
      if (nsm[side] == 0) continue;
      mpz_poly_reduce_frac_mod_f_mod_mpz (rels[i].num[side], rels[i].denom[side],
              F[side], ell2);
      compute_sm (SM, rels[i].num[side], F[side], ell, eps[side], ell2, invl2);
      print_sm (out, SM, nsm[side], F[side]->deg);
      if ((side == 0) && (nsm[1] != 0))
          fprintf(out, " ");
    }
    fprintf(out, "\n");
    // report
    if (stats_test_progress(stats))
      stats_print_progress (stats, i, 0, 0, 0);
  }
  stats_print_progress (stats, sr, 0, 0, 1);

  mpz_poly_clear (SM);
  mpz_clear(invl2);
  fclose(out);
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "(required) poly file");
  param_list_decl_usage(pl, "purged", "(required) purged file");
  param_list_decl_usage(pl, "index", "(required) index file");
  param_list_decl_usage(pl, "out", "output file");
  param_list_decl_usage(pl, "gorder", "(required) group order");
  param_list_decl_usage(pl, "smexp0", "(required) sm-exponent");
  param_list_decl_usage(pl, "nsm0", "number of SM on the 0-side, default deg(polynomial))");
  param_list_decl_usage(pl, "smexp1", "(required) sm-exponent");
  param_list_decl_usage(pl, "nsm1", "number of SM on the 1-side, default deg(polynomial))");
  param_list_decl_usage(pl, "t", "number of threads (default 1)");
  verbose_decl_usage(pl);
}

static void usage (const char *argv, const char * missing, param_list pl)
{
  if (missing) {
    fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
        missing);
  }
  param_list_print_usage(pl, argv, stderr);
  exit (EXIT_FAILURE);
}


/* -------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
  char *argv0 = argv[0];

  const char *polyfile = NULL;
  const char *purgedfile = NULL;
  const char *indexfile = NULL;
  const char *outfile = NULL;

  param_list pl;
  cado_poly pol;
  mpz_poly_ptr F[2];
  sm_relset_ptr rels = NULL;
  uint64_t sr;
  mpz_t ell, ell2, epsilon[2];
  mpz_ptr eps[2];
  eps[0] = &epsilon[0][0];
  eps[1] = &epsilon[1][0];
  int mt = 1;
  double t0;

  /* read params */
  param_list_init(pl);
  declare_usage(pl);

  if (argc == 1)
    usage (argv[0], NULL, pl);

  argc--,argv++;
  for ( ; argc ; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0, NULL, pl);
  }

  /* Read poly filename from command line */
  if ((polyfile = param_list_lookup_string(pl, "poly")) == NULL) {
      fprintf(stderr, "Error: parameter -poly is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  /* Read purged filename from command line */
  if ((purgedfile = param_list_lookup_string(pl, "purged")) == NULL) {
      fprintf(stderr, "Error: parameter -purged is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  /* Read index filename from command line */
  if ((indexfile = param_list_lookup_string(pl, "index")) == NULL) {
      fprintf(stderr, "Error: parameter -index is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  /* Read outfile filename from command line */
  if ((outfile = param_list_lookup_string(pl, "out")) == NULL) {
      fprintf(stderr, "Error: parameter -out is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  /* Read ell from command line (assuming radix 10) */
  mpz_init (ell);
  if (!param_list_parse_mpz(pl, "gorder", ell)) {
      fprintf(stderr, "Error: parameter -gorder is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  /* Read sm exponent from command line (assuming radix 10) */
  for (int side = 0; side < 2; side++) {
      mpz_init (eps[side]);
      char str[10];
      sprintf(str, "smexp%c", '0'+side);
      if (!param_list_parse_mpz(pl, str, eps[side])) {
          fprintf(stderr, "Error: parameter -%s is mandatory\n", str);
          param_list_print_usage(pl, argv0, stderr);
          exit(EXIT_FAILURE);
      }
  }
  param_list_parse_int(pl, "t", &mt);
  if (mt < 1) {
    fprintf(stderr, "Error: parameter mt must be at least 1\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  /* Init polynomial */
  cado_poly_init (pol);
  cado_poly_read(pol, polyfile);
  F[0] = pol->pols[RATIONAL_SIDE];
  F[1] = pol->pols[ALGEBRAIC_SIDE];

  /* Read number of sm to be printed from command line */
  int nsm[2];
  nsm[0] = F[0]->deg;
  nsm[1] = F[1]->deg;
  param_list_parse_int(pl, "nsm0", &nsm[0]);
  param_list_parse_int(pl, "nsm1", &nsm[1]);
  if (nsm[0] + nsm[1] == 0) {
      fprintf(stderr, "Error: no SM to compute!\n");
      exit(EXIT_FAILURE);
  }
  if (nsm[0] > F[0]->deg || nsm[1] > F[1]->deg)
  {
    fprintf(stderr, "Error: nsm can not exceed the degree\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  // If nsm is 0 on one side, then set F[side] to NULL to desactivate the
  // corresponding computations.
  for (int side = 0; side < 2; ++side) {
      if (nsm[side] == 0)
          F[side] = NULL;
  }
  
  if (param_list_warn_unused(pl))
    usage (argv0, NULL, pl);
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);

  /* Print F, ell, smexp and ell2 */
  mpz_init(ell2);
  mpz_mul(ell2, ell, ell);

  for (int side = 0; side < 2; side++) {
      if (F[side] == NULL) continue;
      fprintf(stdout, "\n# Polynomial on side %d:\nF[%d] = ", side, side);
      mpz_poly_fprintf(stdout, F[side]);
      gmp_fprintf(stdout,
              "# Sub-group order:\nell = %Zi\n# Computation is done "
              "modulo ell2 = ell^2:\nell2 = %Zi\n# Shirokauer maps' "
              "exponent:\neps[%d] = %Zi\n", ell, ell2, side, eps[side]);
      fflush(stdout);
  }

  t0 = seconds();
  rels = build_rel_sets(purgedfile, indexfile, &sr, F, ell2);

  /* adjust the number of threads based on the number of relations */
  double ntm = ceil((sr + 0.0)/SM_BLOCK);
  if (mt > ntm)
    mt = (int) ntm;

  fprintf(stdout, "\n# Computing Shirokauer maps for %" PRIu64 " relations "
                  "using %d threads\n", sr, mt);
  fflush(stdout);

  if (mt == 1)
    sm(outfile, rels, sr, F, eps, ell, ell2, nsm);
  else
    mt_sm(mt, outfile, rels, sr, F, eps, ell, ell2, nsm);

  fprintf(stdout, "\n# sm completed in %2.2lf seconds\n", seconds() - t0);
  fflush(stdout);

  for (uint64_t i = 0; i < sr; i++)
    sm_relset_clear (&rels[i]);
  free(rels);
  mpz_clear(eps[0]);
  mpz_clear(eps[1]);
  mpz_clear(ell);
  mpz_clear(ell2);
  cado_poly_clear(pol);
  param_list_clear(pl);

  return 0;
}
