/**
 * @file ropt_linear.c
 * Main function for linear rotation.
 * Called by ropt.c and call stage1.c and stage2.c.
 */


#include "ropt_linear.h"
//#define ROPT_LINEAR_TUNE_HARDER


/**
 * Rank found sublattices by test sieving.
 */
static void
ropt_linear_tune ( ropt_poly_t poly,
                   ropt_bound_t bound,
                   ropt_s1param_t s1param,
                   ropt_param_t param,
                   ropt_info_t info,
                   alpha_pq *alpha_pqueue,
                   MurphyE_pq *global_E_pqueue )
{
  /* tune mode */
  info->mode = 1;

  int i, j, w, used;
  double score, old_MurphyE;
  mpz_t u, tmpu, v, mod, old_mod;
  ropt_s2param_t s2param;
  alpha_pq *tmp_alpha_pqueue;

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
  for (i = 0; i < used; i ++) {

    /* sublattice in w, u, v */
    extract_alpha_pq (alpha_pqueue, &w, u, v, mod, &score);

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

      if (j > 32) break;

#ifdef ROPT_LINEAR_TUNE_HARDER /* slow tunning process*/
      k = 0;
      old_MurphyE = 0.0;
      mpz_set (old_mod, mod);
      while (k < 3) {

        ropt_s2param_setup_tune (s1param, s2param, tmpu, v, mod,
                                 0, size_tune_sievearray, NP-1);
        ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

        if (old_MurphyE > info->best_MurphyE)
          break;
        else
          old_MurphyE = info->best_MurphyE;

        insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                         -info->ave_MurphyE);

        if (param->verbose >= 3) {
          gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                       "sublattice (%d, %Zd, %Zd) (mod %Zd) "
                       "(%d, %d, %d)\n",
                       info->ave_MurphyE, info->best_MurphyE,
                       w, tmpu, v, mod, i + 1, j, k);
        }

        mpz_mul_ui (mod, mod, primes[s1param->tlen_e_sl + k]);

        k ++;
      }
      mpz_set (mod, old_mod);

#else

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

#endif

      mpz_add (tmpu, tmpu, mod); // consider original u itself.
      j ++;
    }

    /* detect negative good u */
    mpz_set (tmpu, u);
    j = 0;
    while (mpz_cmp_si (tmpu, bound->global_u_boundl) > 0) {

      if (j >= 32)
        break;
      mpz_sub (tmpu, tmpu, mod);

#ifdef ROPT_LINEAR_TUNE_HARDER /* slow tunning process*/

      k = 0;
      old_MurphyE = 0.0;
      mpz_set (old_mod, mod);
      while (k < 3) {

        ropt_s2param_setup_tune (s1param, s2param, tmpu, v, mod,
                                 0, size_tune_sievearray, NP-1);
        ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

        if (old_MurphyE > info->best_MurphyE)
          break;
        else
          old_MurphyE = info->best_MurphyE;

        insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                         -info->ave_MurphyE);

        if (param->verbose >= 3) {
          gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                       "sublattice (%d, %Zd, %Zd) (mod %Zd) "
                       "(%d, %d, %d)\n",
                       info->ave_MurphyE, info->best_MurphyE,
                       w, tmpu, v, mod, i + 1, j, k);
        }

        mpz_mul_ui (mod, mod, primes[s1param->tlen_e_sl + k]);

        k ++;
      }
      mpz_set (mod, old_mod);

#else

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

#endif

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
  mpz_clear (u);
  mpz_clear (tmpu);
  mpz_clear (v);
  mpz_clear (mod);
  mpz_clear (old_mod);

  /* normal mode (only control stderr) */
  info->mode = 0;
}


/**
 * Call root sieve
 */
void
ropt_linear_sieve ( ropt_poly_t poly,
                    ropt_bound_t bound,
                    ropt_s1param_t s1param,
                    ropt_param_t param,
                    ropt_info_t info,
                    alpha_pq *alpha_pqueue,
                    MurphyE_pq *global_E_pqueue )
{
  int i, w, w_old = 0, used;
  double score;
  mpz_t u, v, mod, u_old, v_old, mod_old;
  ropt_s2param_t s2param;
  MurphyE_pq *tmp_E_pqueue;

  mpz_init (u);
  mpz_init (v);
  mpz_init (mod);
  mpz_init (u_old);
  mpz_init (v_old);
  mpz_init (mod_old);

  new_MurphyE_pq (&tmp_E_pqueue, s1param->nbest_sl);
  ropt_s2param_init (poly, s2param);

  /* alpha_pqueue contains best ave. E sublattices.
     global_E_pqueue contains best E polynomials. We
     form a final queue using both, with a focus on
     the second. */

  /* Step 1, extract global_E_pqueue */
  used = global_E_pqueue->used - 1;
  for (i = 0; i < used; i ++) {

    extract_MurphyE_pq (global_E_pqueue, &w, u, v, mod, &score);

    if (i == 0) {

      insert_MurphyE_pq (tmp_E_pqueue, w, u, v, mod, score);

      w_old = w;
      mpz_set (u_old, u);
      mpz_set (v_old, v);
      mpz_set (mod_old, mod);

      if (param->verbose >= 2) {
        gmp_fprintf (stderr, "# Info: found best E: %.2e on (#%4d) "
                     "sublattice (%d, %Zd, %Zd) (mod %Zd)\n",
                     score, i + 1, w, u, v, mod);
      }
    }
    else {

      if ( mpz_cmp (u, u_old) == 0 && mpz_cmp (v, v_old) == 0 &&
           mpz_cmp (mod, mod_old) == 0 && w_old == w )
        continue;
      else {

        insert_MurphyE_pq (tmp_E_pqueue, w, u, v, mod, score);

        w_old = w;
        mpz_set (u_old, u);
        mpz_set (v_old, v);
        mpz_set (mod_old, mod);

        if (param->verbose >= 2) {
          gmp_fprintf (stderr, "# Info: found best E: %.2e on (#%4d) "
                       "sublattice (%d, %Zd, %Zd) (mod %Zd)\n",
                       score, i + 1, w, u, v, mod);
        }
      }
    }
  }

  /* re-insert global_E_pqueue as it may contains some top polynomials found
     in sieving test, which may not be found again due to size property. */
  used = tmp_E_pqueue->used - 1;
  for (i = 1; i < used; i ++) {
    insert_MurphyE_pq (global_E_pqueue, tmp_E_pqueue->w[i],
                       tmp_E_pqueue->u[i], tmp_E_pqueue->v[i], 
                       tmp_E_pqueue->modulus[i], tmp_E_pqueue->E[i]);
  }

  /* Step 2, extract alpha_pqueue */
  used = alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {

    /* Note: score here is negative*/
    extract_alpha_pq (alpha_pqueue, &w, u, v, mod, &score);

    if (i == 0) {

      insert_MurphyE_pq (tmp_E_pqueue, w, u, v, mod, -score);

      w_old = w;
      mpz_set (u_old, u);
      mpz_set (v_old, v);
      mpz_set (mod_old, mod);

      if (param->verbose >= 2) {
        gmp_fprintf (stderr, "# Info: found best ave. E: %.2e on (#%4d) "
                     "sublattice (%d, %Zd, %Zd) (mod %Zd)\n",
                     -score, i + 1, w, u, v, mod);
      }
    }
    else {
      if ( mpz_cmp (u, u_old) == 0 && mpz_cmp (v, v_old) == 0 
           && mpz_cmp (mod, mod_old) == 0 && w_old == w )
        continue;

      insert_MurphyE_pq (tmp_E_pqueue, w, u, v, mod, -score);

      w_old = w;
      mpz_set (u_old, u);
      mpz_set (v_old, v);
      mpz_set (mod_old, mod);

      if (param->verbose >= 2) {
        gmp_fprintf (stderr, "# Info: found best ave. E: %.2e on (#%4d) "
                     "sublattice (%d, %Zd, %Zd) (mod %Zd)\n",
                     -score, i + 1, w, u, v, mod);
      }
    }
  }


  /* Step3, final root sieve */
  mpz_t m;
  int old_i = 0;
  mpz_init_set (m, poly->g[0]);
  mpz_neg (m, m);
  used = tmp_E_pqueue->used - 1;
  for (i = 0; i < used; i ++) {

    extract_MurphyE_pq (tmp_E_pqueue, &w, u, v, mod,
                        &score);

    /* rotate */
    old_i = rotate_aux (poly->f, poly->g[1], m, old_i, w, 2);
    ropt_poly_setup (poly);

    if (param->verbose >= 2) {
      gmp_fprintf (stderr, "# Info: Sieve sublattice (# %2d), "
                   "(w, u, v): (%d, %Zd, %Zd) (mod %Zd), "
                   "tsieve. E: %.2e\n",
                   i + 1, w, u, v, mod, score);
    }

    ropt_s2param_setup (bound, s1param, s2param, param, u, v, mod);

    ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

    /* extra sieve for small skewness */
    if (score > info->best_MurphyE) {
      info->mode = 0; // normal mode

      if (param->verbose >= 2) {
        fprintf (stderr, "# Info: Re-sieve sublattice (# %2d) in "
                 "range %u.\n", i + 1, size_tune_sievearray * 2);
      }
     
      ropt_s2param_setup_tune (s1param, s2param, u, v, mod,
                               0, size_tune_sievearray * 2, NP - 1);
      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);
    }

#if 0
    /* a bit longer for the top 3 sublattices */
    if ( used <= (i + 3) ) {
      int j = 0;
      double old_MurphyE = 0.0;
      while (j < 8) {

        mpz_set_ui (mod, default_sublattice_prod[j]);

        if (param->verbose >= 2) {
          gmp_fprintf (stderr, "# Info: Re-sieve sublattice (# %2d) "
                       "for mod=%Zd\n", i + 1, mod);
        }
        ropt_s2param_setup (bound, s1param, s2param, param,
                            u, v, mod);

        ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

        if (old_MurphyE > info->best_MurphyE)
          break;
        else
          old_MurphyE = info->best_MurphyE;
        j ++;
      }
    }
#endif
  }
  rotate_aux (poly->f, poly->g[1], m, old_i, 0, 2);
  ropt_poly_setup (poly);

  /* free */
  free_MurphyE_pq (&tmp_E_pqueue);
  ropt_s2param_free (poly, s2param);
  mpz_clear (m);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);
  mpz_clear (u_old);
  mpz_clear (v_old);
  mpz_clear (mod_old);
}


/**
 * Linear rotation.
 */
void
ropt_linear ( ropt_poly_t poly,
              ropt_bestpoly_t bestpoly,
              ropt_param_t param,
              ropt_info_t info )
{
  int r;
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
  
  /* Step 1:, find good sublattices */
  r = ropt_stage1 (poly, bound, s1param, param, alpha_pqueue, 0);
  if (r == -1) return;
  
  /* Step 2: rank/tune above found sublattices by short sieving */
  ropt_linear_tune (poly, bound, s1param, param, info, alpha_pqueue,
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
