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
  least-significant digit of the ell-adic representation of the polynomial equal
  to (rel^eps - 1) / ell.  
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

#include "filter_common.h"

void *
thread_sm (void * context_data, earlyparsed_relation_ptr rel)
{
  mpz_poly_t * abpolys = (mpz_poly_t *) context_data;
  mpz_poly_init_set_ab(abpolys[rel->num], rel->a, rel->b);

  return NULL;
}

sm_relset_ptr build_rel_sets(const char * purgedname, const char * indexname,
                             uint64_t * small_nrows, mpz_poly_t F,
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

  fprintf(stdout, "\n# Building relation-sets\n");
  for(uint64_t i = 0 ; i < *small_nrows ; i++)
  {
    ret = fscanf(ix, "%" SCNu64 "", &len_relset);
    ASSERT_ALWAYS(ret == 1 && len_relset < MAX_LEN_RELSET);

    for (uint64_t k = 0 ; k < len_relset ; k++)
    {
      ret = fscanf(ix, " %" SCNx64 ":%" SCNd64 "", &r[k], &e[k]); 
      ASSERT_ALWAYS(ret == 2);
    }
     
    sm_relset_init (&rels[i], F->deg);
    sm_build_one_relset (&rels[i], r, e, len_relset, pairs, F, ell2);
  }

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
  mpz_poly_srcptr F;
  mpz_srcptr eps;
  mpz_srcptr ell;
  mpz_srcptr ell2;
  mpz_srcptr invl2;
  mpz_poly_t *sm;
};

void * thread_start(void *arg) {
  struct thread_info *ti = (struct thread_info *) arg;
  sm_relset_ptr rels = ti->rels;
  mpz_poly_srcptr F = ti->F;
  mpz_srcptr eps = ti->eps;
  mpz_srcptr ell = ti->ell;
  mpz_srcptr ell2 = ti->ell2;
  mpz_srcptr invl2 = ti->invl2;
  mpz_poly_t *sm = ti->sm;
  int offset = ti->offset;

  for (int i = 0; i < ti->nb; i++) {
    mpz_poly_reduce_frac_mod_f_mod_mpz(rels[offset+i].num, rels[offset+i].denom,
                                       F, ell2);
    compute_sm (sm[i], rels[offset+i].num, F, ell, eps, ell2, invl2);
  }

  return NULL;
}

#define SM_BLOCK 500

void mt_sm (int nt, const char * outname, sm_relset_ptr rels, uint64_t sr,
            mpz_poly_t F, const mpz_t eps, const mpz_t ell, const mpz_t ell2,
            int nsm)
{
  // allocate space for results of threads
  mpz_poly_t **SM;
  SM = (mpz_poly_t **) malloc(nt*sizeof(mpz_poly_t *));
  for (int i = 0; i < nt; ++i) {
    SM[i] = (mpz_poly_t *) malloc(SM_BLOCK*sizeof(mpz_poly_t));
    for (int j = 0; j < SM_BLOCK; ++j)
      mpz_poly_init(SM[i][j], F->deg);
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
  fprintf(out, "%" PRIu64 "\n", sr);
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
    tis[i].sm = SM[i];
    // offset and nb must be adjusted.
  }

  // Main loop
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
    for (int k = 0; k < SM_BLOCK && out_cpt < sr; ++k, ++out_cpt)
      print_sm (out, SM[threads_head][k], nsm);

    // If we are at the end, no job will be restarted, but head still
    // must be incremented.
    if (i >= sr) { 
      threads_head++;
      if (threads_head == nt) 
        threads_head = 0;
    }
  }

  mpz_clear(invl2);
  fclose(out);
  free(tis);
  free(threads);
  for (int i = 0; i < nt; ++i) {
    for (int j = 0; j < SM_BLOCK; ++j) {
      mpz_poly_clear(SM[i][j]);
    }
    free(SM[i]);
  }
  free(SM);
}


void sm (const char * outname, sm_relset_ptr rels, uint64_t sr, mpz_poly_t F,
         const mpz_t eps, const mpz_t ell, const mpz_t ell2, int nsm)
{
  FILE * out = fopen(outname, "w");
  mpz_poly_t SM;
  mpz_t invl2;

  mpz_init(invl2);
  barrett_init(invl2, ell2);

  mpz_poly_init(SM, F->deg);

  fprintf(out, "%" PRIu64 "\n", sr);

  for (uint64_t i=0; i<sr; i++) {
    mpz_poly_reduce_frac_mod_f_mod_mpz (rels[i].num, rels[i].denom, F, ell2);
    compute_sm (SM, rels[i].num, F, ell, eps, ell2, invl2);
    print_sm (out, SM, nsm);
  }

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
  param_list_decl_usage(pl, "smexp", "(required) sm-exponent");
  param_list_decl_usage(pl, "mt", "number of threads (default 1)");
  param_list_decl_usage(pl, "nsm", "number of SM (default deg(alg polynomial))");
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
  mpz_poly_ptr F;
  sm_relset_ptr rels = NULL;
  uint64_t sr;
  mpz_t ell, ell2, eps;
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
  mpz_init (eps);
  if (!param_list_parse_mpz(pl, "smexp", eps)) {
      fprintf(stderr, "Error: parameter -smexp is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_int(pl, "mt", &mt);
  if (mt < 1) {
    fprintf(stderr, "Error: parameter mt must be at least 1\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  /* Init polynomial */
  cado_poly_init (pol);
  cado_poly_read(pol, polyfile);
  F = pol->pols[ALGEBRAIC_SIDE];

  /* Read number of sm to be printed from command line */
  int nsm = pol->alg->deg;
  param_list_parse_int(pl, "nsm", &nsm);
  if (nsm < 1 || nsm > pol->alg->deg)
  {
    fprintf(stderr, "Error: parameter nsm must be at least 1 and at most "
                    "the degree of the rational polynomial\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }
  
  if (param_list_warn_unused(pl))
    usage (argv0, NULL, pl);
  param_list_print_command_line (stdout, pl);

  /* Print F, ell, smexp and ell2 */
  mpz_init(ell2);
  mpz_mul(ell2, ell, ell);

  fprintf(stdout, "\n# Algebraic polynomial:\nF = ");
  mpz_poly_fprintf(stdout, F);
  gmp_fprintf(stdout, "# Sub-group order:\nell = %Zi\n# Computation is done "
                      "modulo ell2 = ell^2:\nell2 = %Zi\n# Shirokauer maps' "
                      "exponent:\neps = %Zi\n", ell, ell2, eps);

  t0 = seconds();
  rels = build_rel_sets(purgedfile, indexfile, &sr, F, ell2);

  /* adjust the number of threads based on the number of relations */
  double ntm = ceil((sr + 0.0)/SM_BLOCK);
  if (mt > ntm)
    mt = (int) ntm;

  fprintf(stdout, "\n# Computing Shirokauer maps for %" PRIu64 " relations "
                  "using %d threads\n", sr, mt);

  if (mt == 1)
    sm(outfile, rels, sr, F, eps, ell, ell2, nsm);
  else
    mt_sm(mt, outfile, rels, sr, F, eps, ell, ell2, nsm);

  fprintf(stdout, "\n# sm completed in %2.2lf seconds\n", seconds() - t0);

  for (uint64_t i = 0; i < sr; i++)
    sm_relset_clear (&rels[i]);
  free(rels);
  mpz_clear(eps);
  mpz_clear(ell);
  mpz_clear(ell2);
  cado_poly_clear(pol);
  param_list_clear(pl);

  return 0;
}
