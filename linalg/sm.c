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
  to (rel^eps - 1) / ell.  Note: In the very unlikely case where the second lsd
  is zero, the program stops!
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


typedef struct {
  poly_t num;
  poly_t denom;
} relset_struct_t;


typedef relset_struct_t relset_t[1];
typedef relset_struct_t * relset_ptr;
typedef const relset_struct_t * relset_srcptr;


/* Q = P^a mod f, mod p. Note, p is mpz_t */
void
poly_power_mod_f_mod_mpz_Barrett (poly_t Q, const poly_t P, const poly_t f,
				  const mpz_t a, const mpz_t p,
                                  MAYBE_UNUSED const mpz_t invp)
{
  int k = mpz_sizeinbase(a, 2);
  poly_t R;

  if (mpz_cmp_ui(a, 0) == 0) {
    Q->deg = 0;
    mpz_set_ui(Q->coeff[0], 1);
    return;
  }

  poly_alloc(R, 2*f->deg);

  // Initialize R to P
  poly_copy(R, P);

  // Horner
  for (k -= 2; k >= 0; k--)
  {
    poly_sqr_mod_f_mod_mpz(R, R, f, p, invp);  // R <- R^2
    if (mpz_tstbit(a, k))
      poly_mul_mod_f_mod_mpz(R, R, P, f, p, invp);  // R <- R*P
  }

  poly_copy(Q, R);
  poly_free(R);
}


relset_ptr build_rel_sets(const char * purgedname, const char * indexname,
			  int * small_nrows, poly_t F, const mpz_t ell2)
{
  purgedfile_stream ps;
  FILE * ix = fopen(indexname, "r");

  /* array of (a,b) pairs from (purgedname) file */
  /* change to ab_pair *pairs, where ab_pair is uint64_t[2] */
  poly_t *pairs;
  
  /* Array of (small_nrows) relation sets built from array (pairs) and (indexname) file
     rels->num is used for (a,b)-pairs with positive multiplicities, whereas
     rels->denom for those with negative ones
  */
  relset_ptr rels;

  purgedfile_stream_init(ps);
  purgedfile_stream_openfile(ps, purgedname);

  pairs = (poly_t *) malloc(ps->nrows * sizeof(poly_t));

  /* Parse purgedfile (a,b)-pairs only*/
  ps->parse_only_ab = 1;
  int npairs;
  for(npairs = 0 ; purgedfile_stream_get(ps, NULL) >= 0 ; npairs++) {
    ASSERT_ALWAYS(npairs < ps->nrows);
    if (ps->b == 0) {
	/* freerels */
	poly_alloc(pairs[npairs], 0);
	poly_setcoeff_si(pairs[npairs], 0, ps->a);
      }
    else {
      /* (a,b)-pair is a degree-1 poly */
      poly_alloc(pairs[npairs], 1);
      poly_setcoeff_si(pairs[npairs], 0, ps->a);
      poly_setcoeff_si(pairs[npairs], 1, -ps->b);
    }
  }

  /* small_ncols isn't used here: we don't care. */
  int small_ncols;
  int ret = fscanf(ix, "%d %d", small_nrows, &small_ncols);
  ASSERT(ret == 2);

  rels = malloc(*small_nrows * sizeof(relset_t));

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
  }
  poly_free(tmp);

  fclose(ix);
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
  relset_srcptr rels;
  poly_srcptr F;
  mpz_srcptr eps;
  mpz_srcptr ell2;
  mpz_srcptr invl2;
  poly_t *sm;
};

void * thread_start(void *arg) {
  struct thread_info *ti = (struct thread_info *) arg;
  relset_srcptr rels = ti->rels;
  poly_srcptr F = ti->F;
  mpz_srcptr eps = ti->eps;
  mpz_srcptr ell2 = ti->ell2;
  mpz_srcptr invl2 = ti->invl2;
  poly_t *sm = ti->sm;
  int offset = ti->offset;

  poly_t SMn, SMd;
  poly_alloc(SMn, F->deg);
  poly_alloc(SMd, F->deg);

  for (int i = 0; i < ti->nb; i++) {
    poly_power_mod_f_mod_mpz_Barrett(SMn, rels[offset+i].num,
        F, eps, ell2, invl2);
    poly_sub_ui(SMn, 1);

    poly_power_mod_f_mod_mpz_Barrett(SMd, rels[offset+i].denom,
        F, eps, ell2, invl2);
    poly_sub_ui(SMd, 1);

    poly_sub_mod_mpz(sm[i], SMn, SMd, ell2);
  }
  poly_free(SMn);
  poly_free(SMd);
  return NULL;
}

#define SM_BLOCK 500

void mt_sm(int nt, const char * outname, relset_srcptr rels, int sr, poly_t F,
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
    tis[i].ell2 = ell2;
    tis[i].invl2 = invl2;
    tis[i].sm = SM[i];
    // offset and nb must be adjusted.
  }

  // Main loop
  while ((i < sr) || (active_threads > 0)) {
    // Start / restart threads as many threads as allowed
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
    for (int k = 0; k < SM_BLOCK; ++k) {
      if (out_cpt >= sr)
        break;
      poly_ptr sm = SM[threads_head][k];
      for(int j=0; j<F->deg; j++) {
        if (j > sm->deg) {
          fprintf(out, "0 ");
          continue;
        }
        ASSERT_ALWAYS(mpz_divisible_p(sm->coeff[j], ell));
        mpz_divexact(sm->coeff[j], sm->coeff[j], ell);
        ASSERT_ALWAYS(mpz_cmp(ell, sm->coeff[j])>0);

        mpz_out_str(out, 10, sm->coeff[j]);
        fprintf(out, " ");
      }
      fprintf(out, "\n");
      out_cpt++;
    }
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


void shirokauer_maps(const char * outname, relset_srcptr rels, int sr, poly_t F, const mpz_t eps, const mpz_t ell, const mpz_t ell2)
{
  FILE * out = fopen(outname, "w");
  poly_t SMn, SMd, SM;
  mpz_t invl2;
  
  /* mpz_init(l2); */
  mpz_init(invl2);
  //  mpz_init(tmp);

  barrett_init(invl2, ell2);

  mpz_out_str(stderr, 10, invl2);
  fprintf(stderr, "\n");

  fprintf(stderr, "\tBuilding %d relation sets mod ell^2:\n\tell^2 = ", sr);
  mpz_out_str(stderr, 10, ell2);
  fprintf(stderr, "\n");

  poly_alloc(SMn, F->deg);
  poly_alloc(SMd, F->deg);
  poly_alloc(SM, F->deg);

  fprintf(out, "%d\n", sr);

  for (int i=0; i<sr; i++) {

    poly_power_mod_f_mod_mpz_Barrett(SMn, rels[i].num, F, eps, ell2, invl2);
    poly_sub_ui(SMn, 1);

    /* fprintf(stderr, "SMn: "); */
    /* poly_print(SMn); */

    poly_power_mod_f_mod_mpz_Barrett(SMd, rels[i].denom, F, eps, ell2, invl2);
    poly_sub_ui(SMd, 1);

    /* fprintf(stderr, "SMd: "); */
    /* poly_print(SMd); */

    poly_sub_mod_mpz(SM, SMn, SMd, ell2);

    for(int j=0; j<F->deg; j++) {
      if (j > SM->deg) {
          fprintf(out, "0 ");
          continue;
      }
      ASSERT_ALWAYS(mpz_divisible_p(SM->coeff[j], ell));
      mpz_divexact(SM->coeff[j], SM->coeff[j], ell);
      ASSERT_ALWAYS(mpz_cmp(ell, SM->coeff[j])>0);

      mpz_out_str(out, 10, SM->coeff[j]);
      fprintf(out, " ");
    }
    fprintf(out, "\n");

  }

  poly_free(SMn);
  poly_free(SMd);
  poly_free(SM);
  mpz_clear(invl2);
  fclose(out);
}


void usage(const char * me)
{
  fprintf(stderr, "Usage: %s --poly polyname --purged purgedname "
      "--index indexname --out outname --gorder group-order "
      "--smexp sm-exponent [-mt nb_thread]\n", me);
}

/* -------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
  
  const char *purgedname = NULL;
  const char *indexname = NULL;
  const char *outname = NULL;
  const char *group_order = NULL;
  const char *sm_exponent = NULL;

  param_list pl;
  cado_poly pol;
  poly_t F;
  int deg;
  mpz_t *f;
  relset_ptr rels;
  int sr;
  mpz_t ell, ell2, eps;
  int mt = 1;

  double t0;

  char * me = *argv;

  /* print the command line */
  fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
  for (int i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf(stderr, "\n");

  param_list_init(pl);
  argc--,argv++;

  if (!argc) {
      usage(me);
      exit(1);
  }
  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) continue;
    if (strcmp(*argv, "--help") == 0) {
      usage(me);
      exit(0);
    } else {
      fprintf(stderr, "unexpected argument: %s\n", *argv);
      usage(me);
      exit(1);
    }
  }
  
  purgedname = param_list_lookup_string(pl, "purged");
  indexname = param_list_lookup_string(pl, "index");
  outname = param_list_lookup_string(pl, "out");
  group_order = param_list_lookup_string(pl, "gorder");
  sm_exponent = param_list_lookup_string(pl, "smexp");
  param_list_parse_int(pl, "mt", &mt);

  cado_poly_init (pol);

  const char * tmp;

  ASSERT_ALWAYS((tmp = param_list_lookup_string(pl, "poly")) != NULL);
  cado_poly_read(pol, tmp);
  
  if (param_list_warn_unused(pl))
    exit(1);

  /* Construct poly_t F from cado_poly pol (algebraic side) */
  deg = pol->pols[ALGEBRAIC_SIDE]->degree;
  f = pol->pols[ALGEBRAIC_SIDE]->f;
  ASSERT_ALWAYS(deg > 1);
  poly_alloc (F, deg);
  for (int i = deg; i >= 0; --i)
    poly_setcoeff (F, i, f[i]);

  ASSERT_ALWAYS(purgedname != NULL);
  ASSERT_ALWAYS(indexname != NULL);
  ASSERT_ALWAYS(group_order != NULL);
  ASSERT_ALWAYS(sm_exponent != NULL);
  ASSERT_ALWAYS(mt >= 1);

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
  rels = build_rel_sets(purgedname, indexname, &sr, F, ell2);

  fprintf(stderr, "\nComputing Shirokauer maps for %d relations\n", sr);

  if (mt == 1)
    shirokauer_maps(outname, rels, sr, F, eps, ell, ell2);
  else
    mt_sm(mt, outname, rels, sr, F, eps, ell, ell2);

  fprintf(stderr, "\nsm completed in %2.2lf seconds\n", seconds() - t0);

  mpz_clear(eps);
  mpz_clear(ell);
  mpz_clear(ell2);
  poly_free(F);
  cado_poly_clear(pol);
  param_list_clear(pl);

  return 0;
}
