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

* A matrix of (small_nrows) rows and (nmaps)=deg(f) cols (mpz_t).
  For each relation (rel) the (nmaps) Shirokauer maps are computed as the second least-significant digit
  of the ell-adic representation of the polynomial equal to (rel^eps - 1) / ell.
  Note: In the very unlikely case where the second lsd is zero, the program stops!
*/

#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "utils.h"
#include "timing.h"


typedef struct {
  poly_t num;
  poly_t denom;
} relset_struct_t;


typedef relset_struct_t relset_t[1];
typedef relset_struct_t * relset_ptr;
typedef const relset_struct_t * relset_srcptr;


relset_ptr
build_rel_sets (const char * purgedname, const char * indexname ,
                int * small_nrows, poly_t F, mpz_t ell)
{
  purgedfile_stream ps;
  FILE * ix = fopen(indexname, "r");

  /* array of (a,b) pairs from (purgedname) file */
  poly_t *pairs;

  /* Array of (small_nrows) relation sets built from array (pairs) and (indexname) file
     rels->num is used for (a,b)-pairs with positive multiplicities, whereas
     rels->denom for those with negative ones
  */
  relset_ptr rels;

  /* We build relation sets mod ell^2 because Shirokaueur maps will be computed
     mod ell^2
  */
  mpz_t ell2;
  mpz_init(ell2);
  mpz_mul(ell2, ell, ell);

  purgedfile_stream_init(ps);
  purgedfile_stream_openfile(ps, purgedname);

  pairs = (poly_t *) malloc(ps->nrows * sizeof(poly_t));

  /* Parse purgedfile (a,b)-pairs only*/
  ps->parse_only_ab = 1;
  int npairs;
  for(npairs = 0 ; purgedfile_stream_get(ps, NULL) >= 0 ; npairs++) {
    ASSERT_ALWAYS(npairs < ps->nrows);
    /* (a,b)-pair is a degree-1 poly */
    poly_alloc(pairs[npairs], 1);
    poly_setcoeff_si(pairs[npairs], 0, ps->a);
    poly_setcoeff_si(pairs[npairs], 1, -ps->b);
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

  poly_alloc (tmp, F->deg);
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

      mpz_set_si(ee, (e > 0) ? e : -e);
      /* TODO: poly_long_power_mod_f_mod_mpz */
      poly_power_mod_f_mod_mpz(tmp, pairs[ridx], F, ee, ell2);
      poly_mul_mod_f_mod_mpz(rels[i].num, rels[i].num, tmp, F, ell2, NULL);
    }
  }
  poly_free (tmp);

  fclose(ix);
  mpz_clear (ell2);
  purgedfile_stream_closefile(ps);
  purgedfile_stream_clear(ps);
  while (npairs > 0)
    poly_free (pairs[--npairs]);
  free (pairs);
  mpz_clear (ee);
  
  return rels;
}


void shirokauer_maps(const char * outname, relset_srcptr rels, int sr, poly_t F, mpz_t eps, mpz_t ell)
{
  FILE * out = fopen(outname, "w");
  poly_t SMn, SMd, SM;
  mpz_t l2;
  
  mpz_init(l2);

  /* All computations are done mod l^2 */
  mpz_mul(l2, ell, ell);

  poly_alloc(SMn, F->deg);
  poly_alloc(SMd, F->deg);
  poly_alloc(SM, F->deg);

  fprintf(out, "%d\n", sr);

  for (int i=0; i<sr; i++) {

    poly_power_mod_f_mod_mpz(SMn, rels[i].num, F, eps, l2);
    poly_sub_ui(SMn, 1);
  
    poly_power_mod_f_mod_mpz(SMd, rels[i].denom, F, eps, l2);
    poly_sub_ui(SMd, 1);

    poly_sub_mod_mpz(SM, SMn, SMd, l2);

    for(int j=0; j<=SM->deg; j++) {
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
  mpz_clear (l2);
  fclose(out);
}


void usage(const char * me)
{
  fprintf(stderr, "Usage: %s --poly polyname --purged purgedname --index indexname --out outname --gorder group-order --smexp sm-exponent\n", me);
}

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
  mpz_t ell, eps;

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

  cado_poly_init (pol);

  const char * tmp;

  ASSERT_ALWAYS((tmp = param_list_lookup_string(pl, "poly")) != NULL);

  if (param_list_warn_unused(pl))
    exit(1);

  cado_poly_read(pol, tmp);
  
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

  t0 = seconds();

  /* read ell from command line (assuming radix 10) */
  mpz_init_set_str(ell, group_order, 10);

  rels = build_rel_sets(purgedname, indexname, &sr, F, ell);

  fprintf(stderr, "Building %d relation sets took %lf seconds\n", sr, seconds() - t0);

  /* SM exponent from command line (assuming radix 10) */
  mpz_init_set_str(eps, sm_exponent, 10);

  shirokauer_maps(outname, rels, sr, F, eps, ell);

  fprintf(stderr, "Computing Shirokauer maps for %d relations took %2.2lf seconds\n", sr, seconds() - t0);

  mpz_clear (eps);
  poly_free (F);
  cado_poly_clear(pol);
  param_list_clear(pl);

  return 0;
}
