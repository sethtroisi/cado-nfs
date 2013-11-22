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

#include "utils.h"
#include "timing.h"
#include "filter_utils.h"


sm_relset_ptr build_rel_sets(const char * purgedname, const char * indexname,
			  int * small_nrows, poly_t F, const mpz_t ell2)
{
  purgedfile_stream ps;
  FILE * ix = fopen_maybe_compressed(indexname, "r");

  /* array of (a,b) pairs from (purgedname) file */
  poly_t *pairs;
  
  /* Array of (small_nrows) relation sets built from array (pairs) and (indexname) file  */
  sm_relset_ptr rels;

  purgedfile_stream_init(ps);
  purgedfile_stream_openfile(ps, purgedname);

  pairs = (poly_t *) malloc(ps->nrows * sizeof(poly_t));

  /* Parse purgedfile (a,b)-pairs only */
  ps->parse_only_ab = 1;
  uint64_t npairs;
  for(npairs = 0 ; purgedfile_stream_get(ps, NULL) >= 0 ; npairs++) {
    ASSERT_ALWAYS(npairs < ps->nrows);
    poly_alloc_and_set_from_ab(pairs[npairs], ps->a, ps->b);
  }

  /* small_ncols isn't used here: we don't care */
  int small_ncols;
  int ret = fscanf(ix, "%d %d", small_nrows, &small_ncols);
  ASSERT(ret == 2);

  rels = malloc(*small_nrows * sizeof(sm_relset_t));

  for (int k = 0 ; k < *small_nrows ; k++) {
    poly_alloc(rels[k].num, F->deg);
    poly_alloc(rels[k].denom, F->deg);
  }
    
  unsigned int ridx;
  long e, nc;
  poly_t tmp;
  
  mpz_t ee;
  mpz_init(ee);  

  poly_alloc(tmp, F->deg);

  for(int i = 0 ; i < *small_nrows ; i++) {
    ret = fscanf(ix, "%ld", &nc); 
    ASSERT_ALWAYS(ret == 1);

    (rels[i].num)->deg = 0;
    (rels[i].denom)->deg = 0;
    poly_setcoeff_si(rels[i].num, 0, 1);      /* rels[i].num = 1   */
    poly_setcoeff_si(rels[i].denom, 0, 1);    /* rels[i].denom = 1 */

    for(int k = 0 ; k < nc ; k++) {
      ret = fscanf(ix, "%x:%ld", &ridx, &e); 
      ASSERT_ALWAYS(ret == 2);

      /* Should never happen! */
      ASSERT_ALWAYS(e != 0);

      if (e > 0) {
	  mpz_set_si(ee, e);
	  /* TODO: poly_long_power_mod_f_mod_mpz */
	  poly_power_mod_f_mod_mpz(tmp, pairs[ridx], F, ee, ell2);
	  poly_mul_mod_f_mod_mpz(rels[i].num, rels[i].num, tmp, F, ell2, NULL);
      }
      else {
	  mpz_set_si(ee, -e);
	  /* TODO: poly_long_power_mod_f_mod_mpz */
	  poly_power_mod_f_mod_mpz(tmp, pairs[ridx], F, ee, ell2);
	  poly_mul_mod_f_mod_mpz(rels[i].denom, rels[i].denom, tmp, F, ell2, NULL);
      }
    }
    cleandeg(rels[i].num, F->deg);
    cleandeg(rels[i].denom, F->deg);

  }
  poly_free(tmp);

  fclose_maybe_compressed(ix, indexname);
  purgedfile_stream_closefile(ps);
  purgedfile_stream_clear(ps);

  while (npairs > 0)
    poly_free (pairs[--npairs]);
  free (pairs);
  mpz_clear (ee);
  
  return rels;
}

struct thread_info {
  int offset;
  int nb;
  sm_relset_ptr rels;
  poly_srcptr F;
  mpz_srcptr eps;
  mpz_srcptr ell;
  mpz_srcptr ell2;
  mpz_srcptr invl2;
  poly_t *sm;
};

void * thread_start(void *arg) {
  struct thread_info *ti = (struct thread_info *) arg;
  sm_relset_ptr rels = ti->rels;
  poly_srcptr F = ti->F;
  mpz_srcptr eps = ti->eps;
  mpz_srcptr ell = ti->ell;
  mpz_srcptr ell2 = ti->ell2;
  mpz_srcptr invl2 = ti->invl2;
  poly_t *sm = ti->sm;
  int offset = ti->offset;
  mpz_t tmp;

  mpz_init(tmp);

  poly_t g, U, V;
  poly_alloc(g, 0);
  poly_alloc (U,0);
  poly_alloc (V,0);


  for (int i = 0; i < ti->nb; i++) {
    poly_reduce_frac_mod_f_mod_mpz (&rels[offset+i], F, ell2, tmp, g, U, V);
    compute_sm (sm[i], rels[offset+i].num, F, ell, eps, ell2, invl2);
  }

  mpz_clear(tmp);
  poly_free(g);
  poly_free(U);
  poly_free(V);
  return NULL;
}

#define SM_BLOCK 500

void mt_sm(int nt, const char * outname, sm_relset_ptr rels, int sr, poly_t F,
    const mpz_t eps, const mpz_t ell, const mpz_t ell2)
{
  // allocate space for results of threads
  poly_t **SM;
  SM = (poly_t **) malloc(nt*sizeof(poly_t *));
  for (int i = 0; i < nt; ++i) {
    SM[i] = (poly_t *) malloc(SM_BLOCK*sizeof(poly_t));
    for (int j = 0; j < SM_BLOCK; ++j)
      poly_alloc(SM[i][j], F->deg);
  }

  // We'll use a rotating buffer of thread id.
  pthread_t *threads;
  threads = (pthread_t *) malloc(nt*sizeof(pthread_t));
  int active_threads = 0;  // number of running threads
  int threads_head = 0;    // next thread to wait / restart.
  
  // Prepare the main loop
  int i = 0; // counter of relation-sets.
  int out_cpt = 0; // counter of already printed relation-sets;
  FILE * out = fopen(outname, "w");
  fprintf(out, "%d\n", sr);
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
      print_sm (out, SM[threads_head][k], F->deg);

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
      poly_free(SM[i][j]);
    }
    free(SM[i]);
  }
  free(SM);
}


void sm(const char * outname, sm_relset_ptr rels, int sr, poly_t F,
	const mpz_t eps, const mpz_t ell, const mpz_t ell2)
{
  FILE * out = fopen(outname, "w");
  poly_t SM;
  mpz_t invl2;
  mpz_t tmp;
  
  mpz_init(tmp);

  mpz_init(invl2);
  barrett_init(invl2, ell2);


  fprintf(stderr, "\tBuilding %d relation sets mod ell^2:\n", sr);
  fprintf(stderr, "\tell^2 = ");
  mpz_out_str(stderr, 10, ell2);
  fprintf(stderr, "\n");

  poly_alloc(SM, F->deg);
  SM->deg = 0;
  poly_setcoeff_si(SM, 0, 1);

  poly_t g, U, V;
  poly_alloc(g, 0);
  poly_alloc (U,0);
  poly_alloc (V,0);
  
  fprintf(out, "%d\n", sr);

  for (int i=0; i<sr; i++) {
    poly_reduce_frac_mod_f_mod_mpz (&rels[i], F, ell2, tmp, g, U, V);
    compute_sm (SM, rels[i].num, F, ell, eps, ell2, invl2);
    print_sm (out, SM, F->deg);
  }

  poly_free(SM);
  poly_free(U);
  poly_free(V);
  mpz_clear(invl2);
  mpz_clear(tmp);
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
  const char *group_order = NULL;
  const char *sm_exponent = NULL;

  param_list pl;
  cado_poly pol;
  poly_t F;
  sm_relset_ptr rels = NULL;
  int sr;
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


  if ((polyfile = param_list_lookup_string(pl, "poly")) == NULL) {
      fprintf(stderr, "Error: parameter -poly is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  if ((purgedfile = param_list_lookup_string(pl, "purged")) == NULL) {
      fprintf(stderr, "Error: parameter -purged is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  if ((indexfile = param_list_lookup_string(pl, "index")) == NULL) {
      fprintf(stderr, "Error: parameter -index is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  if ((group_order = param_list_lookup_string(pl, "gorder")) == NULL) {
      fprintf(stderr, "Error: parameter -gorder is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  if ((sm_exponent = param_list_lookup_string(pl, "smexp")) == NULL) {
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

  outfile = param_list_lookup_string(pl, "out");

  cado_poly_init (pol);
  cado_poly_read(pol, polyfile);
  
  if (param_list_warn_unused(pl))
    usage (argv0, NULL, pl);
  param_list_print_command_line (stdout, pl);

  /* Construct poly_t F from cado_poly pol (algebraic side) */
  poly_t_from_cado_poly_alg(F, pol);
  fprintf(stderr, "F = ");
  poly_print(F);

  /* read ell from command line (assuming radix 10) */
  mpz_init_set_str(ell, group_order, 10);

  mpz_init(ell2);
  mpz_mul(ell2, ell, ell);

  fprintf(stderr, "\nSub-group order:\n\tell = ");
  mpz_out_str(stderr, 10, ell);
  fprintf(stderr, "\n");

  mpz_init_set_str(eps, sm_exponent, 10);

  fprintf(stderr, "\nShirokauer maps' exponent:\n\teps = ");
  mpz_out_str(stderr, 10, eps);
  fprintf(stderr, "\n");

  t0 = seconds();
  rels = build_rel_sets(purgedfile, indexfile, &sr, F, ell2);

  fprintf(stderr, "\nComputing Shirokauer maps for %d relations\n", sr);

  /* adjust the number of threads based on the number of relations */
  double ntm = ceil((sr + 0.0)/SM_BLOCK);
  if (mt > ntm)
    mt = (int) ntm;

  fprintf(stderr, "using %d threads\n", mt);

  if (mt == 1)
    sm(outfile, rels, sr, F, eps, ell, ell2);
  else
    mt_sm(mt, outfile, rels, sr, F, eps, ell, ell2);

  fprintf(stderr, "\nsm completed in %2.2lf seconds\n", seconds() - t0);

  mpz_clear(eps);
  mpz_clear(ell);
  mpz_clear(ell2);
  poly_free(F);
  cado_poly_clear(pol);
  param_list_clear(pl);

  return 0;
}
