/**
 * @file ropt_quadratic.c
 * Main function for quadratic rotation.
 * Called by ropt.c and then call stage1.c and stage2.c.
 */


#include "ropt_quadratic.h"
#include "ropt.h"


/**
 * Tune stage 1 for quadratic rotation.
 */
static void
ropt_quadratic_stage1 ( ropt_poly_t poly,
                        ropt_bound_t bound,
                        ropt_s1param_t s1param,
                        ropt_param_t param,
                        alpha_pq *alpha_pqueue )
{
  int i, j, k, r, w, old_i, old_nbest_sl, used, *w_good, w_top[3];
  const unsigned int size_tmp_alpha_pqueue = 64;
  double score;
  mpz_t m, u, v, mod;
  alpha_pq *tmp_alpha_pqueue;

  mpz_init (u);
  mpz_init (v);
  mpz_init (mod);
  mpz_init_set (m, poly->g[0]);
  mpz_neg (m, m);
  w_good = (int *) malloc (size_tmp_alpha_pqueue * sizeof(int));

  /* parameters for tunning quadratic rotation */
  j = param->verbose;
  param->verbose = 0;
  old_nbest_sl = s1param->nbest_sl;
  s1param->nbest_sl = 16; // Note: used in ropt_stage1()
  s1param->nbest_sl_tunemode = 1; // Note: used in ropt_stage1()

  /* tmp queue for tunning stage 1 */
  new_alpha_pq (&tmp_alpha_pqueue, size_tmp_alpha_pqueue);

  /* Step 1: for each quadratic rotation, find sublattices quickly */
  k = 0;
  old_i = 0;
  for (i = bound->global_w_boundl; i <= bound->global_w_boundr; i++) {

    if (j >= 2 && k % 10 == 0)
      fprintf (stderr, "# Info: quadratic rotation range %d*x^2\n", i);

    old_i = rotate_aux (poly->f, poly->g[1], m, old_i, i, 2);

    ropt_poly_setup (poly);

    r = ropt_stage1 (poly, bound, s1param, param, tmp_alpha_pqueue, i);

    k ++;
  }

  /* get/rotate back */
  param->verbose = j;
  rotate_aux (poly->f, poly->g[1], m, old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

  /* fetch all w */
  used = tmp_alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {

    /* sublattice in w, u, v */
    extract_alpha_pq (tmp_alpha_pqueue, &w, u, v, mod, &score);
    w_good[i] = w;

#if RANK_SUBLATTICE_BY_E
    if (param->verbose >= 3) {
      gmp_fprintf (stderr, "# Info: got %4d sublattice (%d, %Zd, %Zd) "
                   "(mod %Zd), E %.2f\n", i + 1, w, u, v, mod, score);
    }
#else
    if (param->verbose >= 3) {
      gmp_fprintf (stderr, "# Info: got %4d sublattice (%d, %Zd, %Zd) "
                   "(mod %Zd), alpha: %.2f\n", i + 1, w, u, v, mod, score);
    }
#endif
  }

  /* Step 2: fetch top three w */
  for (i = 0; i < 3; i ++) w_top[i] = bound->global_w_boundr + 1;
  k = 0;
  for (i = used - 1; i > 0; i --) {
    if (k > 2)
      break;    
    r = 0;
    for (j = 0; j < 3; j ++) {
      if (w_top[j] == w_good[i])
        r = 1;
    }
    if (r != 1)
      w_top[k++] = w_good[i];
  }
  
  /* Step 3: for each top w, re-detect (u, v) harder */
  s1param->nbest_sl = old_nbest_sl; // recover the larger nbest_sl
  s1param->nbest_sl_tunemode = 0; // Note: used in ropt_stage1()
  old_i = 0;
  for (i = 0; i < 3; i ++) {
    w = w_top[i];
    
    /* w_top[3] is not full? */
    if (w > bound->global_w_boundr + 1)
      continue;
    
    if (param->verbose >= 2)
      fprintf (stderr, "# Info: find sublattice on quadratic rotation "
               "by %d*x^2\n", w);

    old_i = rotate_aux (poly->f, poly->g[1], m, old_i, w, 2);

    ropt_poly_setup (poly);

    ropt_bound_reset (poly, bound, param);

    r = ropt_stage1 (poly, bound, s1param, param, alpha_pqueue, w);
  }

  /* rotate back */
  rotate_aux (poly->f, poly->g[1], m, old_i, 0, 2);
  ropt_poly_setup (poly);

  old_i = 0;

  /* clear */
  free_alpha_pq (&tmp_alpha_pqueue);
  free (w_good);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);
  mpz_clear (m);
}


/**
 * Rank found sublattices by test sieving.
 */
static void
ropt_quadratic_tune ( ropt_poly_t poly,
                      ropt_bound_t bound,
                      ropt_s1param_t s1param,
                      ropt_param_t param,
                      ropt_info_t info,
                      alpha_pq *alpha_pqueue,
                      MurphyE_pq *global_E_pqueue )
{
  /* tune mode */
  info->mode = 1;

  int i, j, w, used, old_i;
  double score, old_MurphyE;
  mpz_t m, u, tmpu, v, mod, old_mod;
  ropt_s2param_t s2param;
  alpha_pq *tmp_alpha_pqueue;

  mpz_init_set (m, poly->g[0]);
  mpz_neg (m, m);

#ifdef ROPT_LINEAR_TUNE_HARDER
  int k;
#endif

  mpz_init (u);
  mpz_init (tmpu);
  mpz_init (v);
  mpz_init (mod);
  mpz_init (old_mod);

  ropt_s2param_init (poly, s2param);
  new_alpha_pq (&tmp_alpha_pqueue, s1param->nbest_sl);

  /* Step 1: test sieve on alpha_pqueue */
  used = alpha_pqueue->used - 1;
  old_i = 0;

  for (i = 0; i < used; i ++) {

    /* sublattice in w, u, v */
    extract_alpha_pq (alpha_pqueue, &w, u, v, mod, &score);

    old_i = rotate_aux (poly->f, poly->g[1], m, old_i, w, 2);

    ropt_poly_setup (poly);

    // ropt_bound_reset (poly, bound, param);  //may have this

    /* rank by alpha or by E */
#if RANK_SUBLATTICE_BY_E
    if (param->verbose >= 1) {
      gmp_fprintf (stderr, "# Info: tune %4d sublattice (%d, %Zd, %Zd) "
                   "(mod %Zd), E: %.2f\n", i + 1, w, u, v, mod, score);
    }
#else
    if (param->verbose >= 1) {
      gmp_fprintf (stderr, "# Info: tune %4d sublattice (%d, %Zd, %Zd) "
                   "(mod %Zd), alpha: %.2f\n", i + 1, w, u, v, mod, score);
    }
#endif

    /* detect positive good u (may also detect good mod with more time) */
    mpz_set (tmpu, u);
    j = 0;
    while (mpz_cmp_si(tmpu, bound->global_u_boundr) < 0) {

      if (j > 16) break;

      ropt_s2param_setup_tune (s1param, s2param, tmpu, v, mod,
                               0, size_tune_sievearray, NP-1);

      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

      insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                       -info->ave_MurphyE);

      if (param->verbose >= 3) {
        gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                     "sublattice (%d, %Zd, %Zd) (mod %Zd) (%d, %d)\n",
                     info->ave_MurphyE, info->best_MurphyE,
                     w, tmpu, v, mod, i + 1, j);
      }

      mpz_add (tmpu, tmpu, mod); // consider original u itself.
      j ++;
    }
    
    /* detect negative good u */
    mpz_set (tmpu, u);
    j = 0;
    while (mpz_cmp_si (tmpu, bound->global_u_boundl) > 0) {

      if (j >= 16)
        break;
      mpz_sub (tmpu, tmpu, mod);

      ropt_s2param_setup_tune (s1param, s2param, tmpu, v, mod,
                               0, size_tune_sievearray, NP-1);

      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

      insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                       -info->ave_MurphyE);

      if (param->verbose >= 3) {
        gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                     "sublattice (%d, %Zd, %Zd) (mod %Zd) (%d, %d)\n",
                     info->ave_MurphyE, info->best_MurphyE,
                     w, tmpu, v, mod, i + 1, j);
      }
      j ++;
    }
  }

  /* Step 2: slight larger range sieve on best sublattices */
  reset_alpha_pq (alpha_pqueue);
  used =  tmp_alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {

    extract_alpha_pq (tmp_alpha_pqueue, &w, u, v, mod, &score);

    if (param->verbose >= 1) {
      gmp_fprintf ( stderr, "# Info: ave. E: %.2e, tune (#%4d ) sublattice "
                    "(%d, %Zd, %Zd) (mod %Zd)\n",
                    -score, i + 1, w, u, v, mod );
    }

    insert_alpha_pq (alpha_pqueue, w, u, v, mod, score);

    j = 0;
    old_MurphyE = -score;
    while (j < 8) {

      mpz_mul_ui (mod, mod, primes[s1param->tlen_e_sl + j]);

      ropt_s2param_setup_tune (s1param, s2param, u, v, mod,
                               0, size_tune_sievearray, 20);
      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

      if (old_MurphyE > info->ave_MurphyE)
        break;

      old_MurphyE = info->ave_MurphyE;

      insert_alpha_pq (alpha_pqueue, w, u, v, mod, -info->ave_MurphyE);

      j ++;
    }
  }

  /* free s2param */
  free_alpha_pq (&tmp_alpha_pqueue);
  ropt_s2param_free (poly, s2param);
  mpz_clear (m);
  mpz_clear (u);
  mpz_clear (tmpu);
  mpz_clear (v);
  mpz_clear (mod);
  mpz_clear (old_mod);

  /* normal mode (only control stderr) */
  info->mode = 0;
}


/**
 * Quadratic rotation.
 */
void
ropt_quadratic ( ropt_poly_t poly,
                 ropt_bestpoly_t bestpoly,
                 ropt_param_t param,
                 ropt_info_t info )
{
  ropt_bound_t bound;
  ropt_s1param_t s1param;
  alpha_pq *alpha_pqueue;
  MurphyE_pq *global_E_pqueue;

  /* setup bound, s1param, alpha_pqueue, tsieve_E_pqueue */
  ropt_bound_init (bound);
  ropt_bound_setup (poly, bound, param);
  ropt_s1param_init (s1param);
  ropt_s1param_setup (poly, s1param, bound, param);
  new_alpha_pq (&alpha_pqueue, s1param->nbest_sl);
  new_MurphyE_pq (&global_E_pqueue, s1param->nbest_sl);

  /* Step 1: find good sublattice */
  ropt_quadratic_stage1 (poly, bound, s1param, param, alpha_pqueue);


  /* Step 2: rank/tune above found sublattices by short sieving */
  ropt_quadratic_tune (poly, bound, s1param, param, info, alpha_pqueue,
                       global_E_pqueue);

  /* Step 3, root sieve */
  ropt_linear_sieve (poly, bound, s1param, param, info, alpha_pqueue,
                        global_E_pqueue);

  /* Step 4, return best poly */
  ropt_get_bestpoly (poly, global_E_pqueue, bestpoly);

  /* free */
  free_MurphyE_pq (&global_E_pqueue);
  free_alpha_pq (&alpha_pqueue);
  ropt_s1param_free (s1param);
  ropt_bound_free (bound);
}
