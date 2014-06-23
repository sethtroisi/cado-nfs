/**
 * @file ropt_quadratic.c
 * Main function for quadratic rotation.
 * Called by ropt.c and then call stage1.c and stage2.c.
 */


#include "cado.h"
#include "ropt_quadratic.h"
#include "ropt.h"
#include "portability.h"

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
  const int numw = 3;
  int i, j, k, r, w, old_i, old_nbest_sl, used, *w_good, w_top[numw];
  const unsigned int size_alpha_pqueue_all_w =
      s1param->nbest_sl * TUNE_NUM_SUBLATTICE;
  double score;
  mpz_t m, u, v, mod;
  alpha_pq *alpha_pqueue_all_w;
  
  mpz_init (u);
  mpz_init (v);
  mpz_init (mod);
  mpz_init_set (m, poly->g[0]);
  mpz_neg (m, m);
  w_good = (int *) malloc (size_alpha_pqueue_all_w * sizeof(int));

  /* tmp queue for tunning stage 1 */
  new_alpha_pq (&alpha_pqueue_all_w, size_alpha_pqueue_all_w);

  /* parameters for tunning quadratic rotation */
  j = param->verbose;
  param->verbose = 0;
  old_nbest_sl = s1param->nbest_sl;

  /* used in ropt_stage1(), the larger the slower for each
     quadratic rotation w. */
  s1param->nbest_sl = 16;

  /* used in ropt_stage1(), uses this to reduce the length of 
     individual sublattice queue */
  s1param->nbest_sl_tunemode = 1;

  /* Step 1: for each quadratic rotation, find sublattices in a
     quick, but approximated manner */
  k = 0;
  old_i = 0;
  for (i = bound->global_w_boundl; i <= bound->global_w_boundr; i++) {
    if (j >= 2 && k % 10 == 0)
      fprintf (stderr, "# Info: quadratic rotation range %d*x^2\n", i);
    old_i = rotate_aux (poly->f, poly->g[1], m, old_i, i, 2);
    ropt_poly_setup (poly);
    r = ropt_stage1 (poly, bound, s1param, param, alpha_pqueue_all_w, i);
    k ++;
  }

  /* get/rotate back */
  param->verbose = j;
  rotate_aux (poly->f, poly->g[1], m, old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

  /* fetch all w */
  used = alpha_pqueue_all_w->used - 1;
  for (i = 0; i < used; i ++) {

    /* sublattice in w, u, v */
    extract_alpha_pq (alpha_pqueue_all_w, &w, u, v, mod, &score);
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
  for (i = 0; i < numw; i ++)
    w_top[i] = bound->global_w_boundr + 1;
  k = 0;
  for (i = used - 1; i > 0; i --) {
    if (k > 2)
      break;    
    r = 0;
    for (j = 0; j < numw; j ++) {
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
  for (i = 0; i < numw; i ++) {
    w = w_top[i];
    
    /* w_top[3] is not full? */
    if (w >= bound->global_w_boundr + 1)
      continue;
    
    if (param->verbose >= 2)
      fprintf (stderr, "\n# Info: find sublattice on quadratic rotation "
               "by %d*x^2\n", w);

    old_i = rotate_aux (poly->f, poly->g[1], m, old_i, w, 2);

    ropt_poly_setup (poly);

    ropt_bound_reset (poly, bound, param);

    r = ropt_stage1 (poly, bound, s1param, param, alpha_pqueue, w);
  }
  if (param->verbose >= 2) fprintf (stderr, "\n");

  /* rotate back */
  rotate_aux (poly->f, poly->g[1], m, old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

  /* clear */
  free_alpha_pq (&alpha_pqueue_all_w);
  free (w_good);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);
  mpz_clear (m);
}


/**
 * Rank found sublattices in alpha_pqueue by test sieving.
 * Note: alpha_pqueue contains more (double/quad) sublattices than 
 * nbest_sl (final sublattices to be root-sieved). In this function,
 * - the first pass of tunning reduces #alpha_pqueue to nbest_sl.
 * - the second pass of tunning identifies good u and v.
 * - the third pass of tunning identifies good mod.
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
  /* tune mode, prevent from fprint */
  info->mode = 1;

  int i, j, w, used, old_i;
  double score, old_MurphyE;
  mpz_t m, u, tmpu, v, mod, old_mod;
  ropt_s2param_t s2param;
  alpha_pq *tmp_alpha_pqueue;

  mpz_init_set (m, poly->g[0]);
  mpz_neg (m, m);
  mpz_init (u);
  mpz_init (tmpu);
  mpz_init (v);
  mpz_init (mod);
  mpz_init (old_mod);

  ropt_s2param_init (poly, s2param);
  new_alpha_pq (&tmp_alpha_pqueue, s1param->nbest_sl);

  /* Step 1: test sieve on alpha_pqueue, reduce the sublattice number
     to the s1param->nbest_sl, and record them to tmp_alpha_pqueue. */
  used = alpha_pqueue->used - 1;
  old_i = 0;
  for (i = 0; i < used; i ++) {

    /* sublattice in w, u, v */
    extract_alpha_pq (alpha_pqueue, &w, u, v, mod, &score);

    old_i = rotate_aux (poly->f, poly->g[1], m, old_i, w, 2);

    ropt_poly_setup (poly);

    //ropt_bound_reset (poly, bound, param); // not really necessary

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

    ropt_s2param_setup_tune (s1param, s2param, u, v, mod,
                             0, size_tune_sievearray, NP-1);

    ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

    insert_alpha_pq (tmp_alpha_pqueue, w, u, v, mod,
                     -info->ave_MurphyE);

    if (param->verbose >= 3) {
      gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                   "sublattice (%d, %Zd, %Zd) (mod %Zd)\n",
                   info->ave_MurphyE, info->best_MurphyE,
                   w, u, v, mod);
    }
  }
  /* rotate back */
  rotate_aux (poly->f, poly->g[1], m, old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;


  /* Step 2: detect good (u, v). This is the slow step. */
  reset_alpha_pq (alpha_pqueue);
  used = tmp_alpha_pqueue->used - 1;
  old_i = 0;
  for (i = 0; i < used; i ++) {

    /* sublattice in w, u, v */
    extract_alpha_pq (tmp_alpha_pqueue, &w, u, v, mod, &score);

    old_i = rotate_aux (poly->f, poly->g[1], m, old_i, w, 2);

    ropt_poly_setup (poly);

    // ropt_bound_reset (poly, bound, param);  //may have this
    if (param->verbose >= 1) {
      gmp_fprintf (stderr, "# Info: tune %4d sublattice (%d, %Zd, %Zd) "
                   "(mod %Zd), E: %.2e\n", i + 1, w, u, v, mod, -score);
    }

    /* detect positive good u (may also detect good mod with more time) */
    mpz_set (tmpu, u);
    j = 0;
    while (mpz_cmp_si(tmpu, bound->global_u_boundr) < 0) {

      if (j > 16) break;

      ropt_s2param_setup_tune (s1param, s2param, tmpu, v, mod,
                               0, size_tune_sievearray, NP-1);

      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

      insert_alpha_pq (alpha_pqueue, w, tmpu, v, mod,
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

      insert_alpha_pq (alpha_pqueue, w, tmpu, v, mod,
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
  /* rotate back */
  rotate_aux (poly->f, poly->g[1], m, old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

  /* put alpha_pqueue to tmp_alpha_pqueue */
  used = alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {
    extract_alpha_pq (alpha_pqueue, &w, u, v, mod, &score);
    insert_alpha_pq (tmp_alpha_pqueue, w, u, v, mod, score);
  }

  /* Step 3: detect good mod */
  reset_alpha_pq (alpha_pqueue);
  used = tmp_alpha_pqueue->used - 1;
  old_i = 0;
  for (i = 0; i < used; i ++) {

    extract_alpha_pq (tmp_alpha_pqueue, &w, u, v, mod, &score);

    old_i = rotate_aux (poly->f, poly->g[1], m, old_i, w, 2);

    ropt_poly_setup (poly);

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
  /* rotate back */
  rotate_aux (poly->f, poly->g[1], m, old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

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
 * Call root sieve. It will root-sieve two passes of the following
 * two queues:
 * - alpha_pqueue which records 'average best E' sublattices.
 * - global_E_pqueue which records 'top E' sublattices.
 */
void
ropt_quadratic_sieve ( ropt_poly_t poly,
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
  }
  /* rotate back */
  rotate_aux (poly->f, poly->g[1], m, old_i, 0, 2);
  ropt_poly_setup (poly);
  old_i = 0;

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
 * Quadratic rotation.
 */
void
ropt_quadratic ( ropt_poly_t poly,
                 ropt_bestpoly_t bestpoly,
                 ropt_param_t param,
                 ropt_info_t info )
{
  unsigned long t1, t2, t3;
  ropt_bound_t bound;
  ropt_s1param_t s1param;
  alpha_pq *alpha_pqueue;
  MurphyE_pq *global_E_pqueue;

  /* setup bound, s1param, alpha_pqueue, tsieve_E_pqueue */
  ropt_bound_init (bound);
  ropt_bound_setup (poly, bound, param);
  ropt_s1param_init (s1param);
  ropt_s1param_setup (poly, s1param, bound, param);

  /* nbest_sl is the number of sublattice to be finally root-sieved.
     Here we use some larger value since alpha_pqueue records more 
     sublattices to be tunned */
  new_alpha_pq (&alpha_pqueue, s1param->nbest_sl * TUNE_NUM_SUBLATTICE);
  new_MurphyE_pq (&global_E_pqueue, s1param->nbest_sl);

  /* Step 1: find good sublattice */
  t1 = milliseconds ();
  ropt_quadratic_stage1 (poly, bound, s1param, param, alpha_pqueue);
  t1 = milliseconds () - t1;
  
  /* Step 2: rank/tune above found sublattices by short sieving */
  t2 = milliseconds ();
  ropt_quadratic_tune (poly, bound, s1param, param, info, alpha_pqueue,
                       global_E_pqueue);
  t2 = milliseconds () - t2;
  
  /* Step 3, root sieve */
  t3 = milliseconds ();
  ropt_quadratic_sieve (poly, bound, s1param, param, info, alpha_pqueue,
                        global_E_pqueue);
  t3 = milliseconds () - t3;
  
  if (param->verbose >= 2) {
    fprintf ( stderr, "# Stat: tot (stage 1) took %lums\n", t1 );
    fprintf ( stderr, "# Stat: tot (tuning ) took %lums\n", t2 );
    fprintf ( stderr, "# Stat: tot (stage 2) took %lums\n", t3 );
  }

  /* Step 4, return best poly */
  ropt_get_bestpoly (poly, global_E_pqueue, bestpoly);
  
  /* free */
  free_MurphyE_pq (&global_E_pqueue);
  free_alpha_pq (&alpha_pqueue);
  ropt_s1param_free (s1param);
  ropt_bound_free (bound);
}
