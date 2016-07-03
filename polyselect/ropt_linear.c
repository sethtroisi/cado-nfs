/**
 * @file ropt_linear.c
 * Main function for linear rotation.
 * Called by ropt.c and call stage1.c and stage2.c.
 *
 * This is description of the flow:
 * - first we call rotp_stage1() to get N lats based only on alpha;
 * - then fast_tune() gets N' << N good lats based on crude test-sieve;
 * - then slow_tune() gets N'' ~ N lats based on finer and recursieve test-sieve;
 * - finally we root-sieve the top few candidates from above.
 *
 * For the above figures:
 *  N  = size_total_sublattices[][1] * ropteffort * TUNE_NUM_SUBLATTICE_STAGE1
 *  N' = size_total_sublattices[][1] * ropteffort
 * 
 *
 */ 


#include "cado.h"
#include "ropt_linear.h"
#include "portability.h"
//#define ROPT_LINEAR_TUNE_HARDER


/**
 * This is a fast tuning process which deal with 
 */
static double
ropt_linear_tune_fast ( ropt_poly_t poly,
                        ropt_s1param_t s1param,
                        ropt_param_t param,
                        ropt_info_t info,
                        alpha_pq *alpha_pqueue,
                        MurphyE_pq *global_E_pqueue,
                        unsigned int curr_size_tune )
{
  /* tune mode */
  info->mode = 1;

  int i, w, used;
  double score, max_score;
  mpz_t u, v, mod;
  ropt_s2param_t s2param;
  alpha_pq *tmp_alpha_pqueue;
  
  mpz_init (u);
  mpz_init (v);
  mpz_init (mod);
  ropt_s2param_init (poly, s2param);
  new_alpha_pq (&tmp_alpha_pqueue, s1param->nbest_sl);
  
  /* fast test-sieve on lats from stage 1 */
  used = alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {

    /* sublattice in w, u, v */
    extract_alpha_pq (alpha_pqueue, &w, u, v, mod, &score);

#if RANK_SUBLATTICE_BY_E
    char stmp[] = "E";  /* RANK_SUBLATTICE_BY_E */
#else
    char stmp[] = "alpha";
#endif
    if (param->verbose >= 2) {
      gmp_fprintf (stderr, "# Info: tune [%4d], %s: %.2f, lat (%d, %Zd, %Zd) "
                   "(mod %Zd)\n", i+1,  stmp, score, w, u, v, mod);
    }
    
    ropt_s2param_setup_tune (s1param, s2param, u, v, mod,
                             0, curr_size_tune, NP-1);
    ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);
    insert_alpha_pq (tmp_alpha_pqueue, w, u, v, mod, -info->ave_MurphyE);
    if (param->verbose >= 4)
      gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e, "
                   "lat (%d, %Zd, %Zd) (mod %Zd)\n",
                   info->ave_MurphyE, info->best_MurphyE,
                   w, u, v, mod);
  }
  
  /* End: save to alpha_pqueue */
  reset_alpha_pq (alpha_pqueue);
  used =  tmp_alpha_pqueue->used - 1;
  max_score = 0;
  for (i = 0; i < used; i ++) {
    extract_alpha_pq (tmp_alpha_pqueue, &w, u, v, mod, &score);
    if (param->verbose >= 3)
      gmp_fprintf (stderr, "# Info:  got [%4d], E: %.2e, "
                   "lat (%d, %Zd, %Zd) (mod %Zd)\n",
                   i+1, -score, w, u, v, mod);
    insert_alpha_pq (alpha_pqueue, w, u, v, mod, score);
    if (max_score < -score)
      max_score = score;
  }
  max_score = -max_score;
  
  /* free s2param */
  free_alpha_pq (&tmp_alpha_pqueue);
  ropt_s2param_free (poly, s2param);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);

  /* normal mode (only control stderr) */
  info->mode = 0;
  return max_score;
}


/**
 * Rank found sublattices by test sieving.
 */
static void
ropt_linear_tune_slow ( ropt_poly_t poly,
                        ropt_bound_t bound,
                        ropt_s1param_t s1param,
                        ropt_param_t param,
                        ropt_info_t info,
                        alpha_pq *alpha_pqueue,
                        MurphyE_pq *global_E_pqueue,
                        unsigned int curr_size_tune,
                        int tune_round,
                        unsigned int curr_nbest )
{
  /* tune mode */
  info->mode = 1;

  int i, j, w, used;
  double score, old_MurphyE;
  mpz_t u, tmpu, v, mod, old_mod;
  ropt_s2param_t s2param;
  alpha_pq *tmp_alpha_pqueue, *tmp2_alpha_pqueue;
#if TUNE_EARLY_ABORT
  double ave_E, ave_bestE, new_ave_E, new_ave_bestE;
  int acc;
#endif  
  
#ifdef ROPT_LINEAR_TUNE_HARDER
  int k;
#endif

  mpz_init (u);
  mpz_init (tmpu);
  mpz_init (v);
  mpz_init (mod);
  mpz_init (old_mod);

  ropt_s2param_init (poly, s2param);
  if (curr_nbest > (unsigned) alpha_pqueue->used-1)
    curr_nbest = alpha_pqueue->used-1;
  if (curr_nbest < 2) curr_nbest = 2;
  new_alpha_pq (&tmp_alpha_pqueue, curr_nbest);
  new_alpha_pq (&tmp2_alpha_pqueue, curr_nbest);
  
  /* Step 1: test sieve on alpha_pqueue */
  used = alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {

    /* sublattice in w, u, v */
    extract_alpha_pq (alpha_pqueue, &w, u, v, mod, &score);

#if RANK_SUBLATTICE_BY_E
    char stmp[] = "E";  /* RANK_SUBLATTICE_BY_E */
#else
    char stmp[] = "alpha";
#endif
    if (param->verbose >= 1) {
      if (tune_round==1) {
        /* first tune round, score are in E, not murphyE */
        gmp_fprintf (stderr, "# Info: tune [%4d], %s: %.2f, lat (%d, %Zd, %Zd) "
                     "(mod %Zd)\n", i+1,  stmp, score, w, u, v, mod);
      }
      else {
        gmp_fprintf (stderr, "# Info: tune [%4d], E: %.2e, "
                     "lat (%d, %Zd, %Zd) (mod %Zd)\n",
                     -score, i + 1, w, u, v, mod);
      }
    }
    
    /* detect positive good u (may also detect good mod with more time) */
    mpz_set (tmpu, u);
    j = 0;
#if TUNE_EARLY_ABORT
    ave_E = 0.0;
    ave_bestE = 0.0;
    new_ave_E = 0.0;
    new_ave_bestE = 0.0;
    acc = 0;
#endif
    
    while (mpz_cmp_si(tmpu, bound->global_u_boundr) <= 0) {
      
      if (j > TUNE_BOUND_ON_UV_TRIALS)
        break;

#ifdef ROPT_LINEAR_TUNE_HARDER /* slow tuning process */
      k = 0;
      old_MurphyE = 0.0;
      mpz_set (old_mod, mod);
      while (k < 3) {

        ropt_s2param_setup_tune (s1param, s2param, tmpu, v, mod,
                                 0, curr_size_tune, NP-1);

        ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

        if (old_MurphyE > info->best_MurphyE)
          break;
        else
          old_MurphyE = info->best_MurphyE;

        insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                         -info->ave_MurphyE);

        if (param->verbose >= 3) {
          gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                       "lat (%d, %Zd, %Zd) (mod %Zd) "
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
                               0, curr_size_tune, NP-1);

      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

      insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                       -info->ave_MurphyE);

      if (param->verbose >= 3) {
        gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                     "lat (%d, %Zd, %Zd) (mod %Zd) (%d, %d)\n",
                     info->ave_MurphyE, info->best_MurphyE,
                     w, tmpu, v, mod, i + 1, j);
      }
      
#endif

      mpz_add (tmpu, tmpu, mod); // consider original u itself.
#if TUNE_EARLY_ABORT
      new_ave_E = ave_E + info->ave_MurphyE;
      new_ave_bestE = ave_bestE + info->best_MurphyE;
      if (j>0) {
        if ((ave_E/j>new_ave_E/(j+1)) && (ave_bestE/j > new_ave_bestE/(j+1))) {
          acc ++;
          
          if (acc >= TUNE_EARLY_ABORT_THR) {
            //printf (" early abort after acc %d\n", acc);
            break;
          }
          
        }
      }
      ave_E = new_ave_E;
      ave_bestE = new_ave_bestE;
#endif      
      j ++;
    }

    /* detect negative good u */
    mpz_set (tmpu, u);
    j = 0;
#if TUNE_EARLY_ABORT
    ave_E = 0.0;
    ave_bestE = 0.0;
    new_ave_E = 0.0;
    new_ave_bestE = 0.0;
    acc = 0;
#endif

    while (mpz_cmp_si (tmpu, bound->global_u_boundl) > 0) {

      if (j > TUNE_BOUND_ON_UV_TRIALS)
        break;
      mpz_sub (tmpu, tmpu, mod);

#ifdef ROPT_LINEAR_TUNE_HARDER /* slow tuning process */

      k = 0;
      old_MurphyE = 0.0;
      mpz_set (old_mod, mod);
      while (k < 3) {

        ropt_s2param_setup_tune (s1param, s2param, tmpu, v, mod,
                                 0, curr_size_tune, NP-1);
        ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

        if (old_MurphyE > info->best_MurphyE)
          break;
        else
          old_MurphyE = info->best_MurphyE;

        insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                         -info->ave_MurphyE);

        if (param->verbose >= 3) {
          gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                       "lat (%d, %Zd, %Zd) (mod %Zd) "
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
                               0, curr_size_tune, NP-1);

      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

      insert_alpha_pq (tmp_alpha_pqueue, w, tmpu, v, mod,
                       -info->ave_MurphyE);

      if (param->verbose >= 3) {
        gmp_fprintf (stderr, "# Info: ave. E: %.2e, best E: %.2e on "
                     "lat (%d, %Zd, %Zd) (mod %Zd) (%d, %d)\n",
                     info->ave_MurphyE, info->best_MurphyE,
                     w, tmpu, v, mod, i + 1, j);
      }

#endif

#if TUNE_EARLY_ABORT
      new_ave_E = ave_E + info->ave_MurphyE;
      new_ave_bestE = ave_bestE + info->best_MurphyE;
      if (j>0) {
        if ((ave_E/j>new_ave_E/(j+1)) && (ave_bestE/j > new_ave_bestE/(j+1))) {
          acc ++;
          if (acc >= TUNE_EARLY_ABORT_THR)
            break;
        }
      }
      ave_E = new_ave_E;
      ave_bestE = new_ave_bestE;
#endif      

      j ++;

    }
  }
  
  /* Step 2: slight larger range sieve on best sublattices */
  used =  tmp_alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {
    extract_alpha_pq (tmp_alpha_pqueue, &w, u, v, mod, &score);
    if (param->verbose >= 1) {
      gmp_fprintf (stderr, "# Info: tune [%4d], E: %.2e, "
                   "lat (%d, %Zd, %Zd) (mod %Zd)\n",
                   -score, i + 1, w, u, v, mod);
    }

    insert_alpha_pq (tmp2_alpha_pqueue, w, u, v, mod, score);

    j = 0;
    old_MurphyE = -score;
    while (j < 8) {

      mpz_mul_ui (mod, mod, primes[s1param->tlen_e_sl + j]);

      ropt_s2param_setup_tune (s1param, s2param, u, v, mod,
                               0, curr_size_tune, 20);
      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);

      /* abort if E is becoming worse */
      if (old_MurphyE > info->ave_MurphyE)
        break;

      old_MurphyE = info->ave_MurphyE;
      insert_alpha_pq (tmp2_alpha_pqueue, w, u, v, mod, -info->ave_MurphyE);
      j ++;
    }
  }

  /* End: save to alpha_pqueue */
  reset_alpha_pq (alpha_pqueue);
  used =  tmp2_alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {
    extract_alpha_pq (tmp2_alpha_pqueue, &w, u, v, mod, &score);
    insert_alpha_pq (alpha_pqueue, w, u, v, mod, score);
  }
  
  /* free s2param */
  free_alpha_pq (&tmp_alpha_pqueue);
  free_alpha_pq (&tmp2_alpha_pqueue);
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
 * Recursive tuning
 */
static void
ropt_linear_tune ( ropt_poly_t poly,
                   ropt_bound_t bound,
                   ropt_s1param_t s1param,
                   ropt_param_t param,
                   ropt_info_t info,
                   alpha_pq *alpha_pqueue,
#if TUNE_LOGNORM_INCR
                   alpha_pq *pre_alpha_pqueue,
#endif
                   MurphyE_pq *global_E_pqueue )
{
#if TUNE_LOGNORM_INCR
  int i, w, used;
  double score;
  mpz_t u, v, mod;
  mpz_init (u);
  mpz_init (v);
  mpz_init (mod);
#endif
  unsigned int old_nbest = s1param->nbest_sl;
  unsigned int curr_size_tune = size_tune_sievearray;
  int r = 1;

  /* step 2: fast tune for lats from stage 1 */
  if (param->verbose >= 1) {
    printf ("# Info: tune (round %d), tunesize: %u, queuesize: %u, new_queuesize: %u\n",
            r, curr_size_tune, alpha_pqueue->used-1, s1param->nbest_sl);
  }
  ropt_linear_tune_fast (poly, s1param, param, info, alpha_pqueue,
                         global_E_pqueue, curr_size_tune);
  r ++;


#if TUNE_LOGNORM_INCR
  /* aux: in case pretune not empty, get it */
  used = pre_alpha_pqueue->used - 1;
  for (i = 0; i < used; i ++) {
    extract_alpha_pq (pre_alpha_pqueue, &w, u, v, mod, &score);
    if (param->verbose >= 3) {
      gmp_fprintf (stderr, "# Info: pre-tune [%4d], E: %.2e, lat (%d, %Zd, %Zd) "
                   "(mod %Zd)\n", i+1, -score, w, u, v, mod);
    }
    insert_alpha_pq (alpha_pqueue, w, u, v, mod, score);
   }
#endif

  /* step 3: recursieve (slower but finer) tune */
  while (1) {
    if (param->verbose >= 1) {
      printf ("# Info: tune (round %d), tunesize: %u, nbest: %u\n",
              r, curr_size_tune, s1param->nbest_sl);
    }
    curr_size_tune = curr_size_tune*2;
    s1param->nbest_sl = s1param->nbest_sl/2;
    ropt_linear_tune_slow (poly, bound, s1param, param, info, alpha_pqueue,
      global_E_pqueue, curr_size_tune, r, s1param->nbest_sl);
    r ++;
    remove_rep_alpha (alpha_pqueue);
    /* doing 2 rounds of tuning seems enough */
    if (s1param->nbest_sl < old_nbest/2) break;
  }
  s1param->nbest_sl = old_nbest;

#if TUNE_LOGNORM_INCR
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);
#endif
  return;
}


/**
 *a Call root sieve
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

  new_MurphyE_pq (&tmp_E_pqueue, s1param->nbest_sieve);
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
                     "lat (%d, %Zd, %Zd) (mod %Zd)\n",
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
                       "lat (%d, %Zd, %Zd) (mod %Zd)\n",
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
                     "lat (%d, %Zd, %Zd) (mod %Zd)\n",
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
                     "lat (%d, %Zd, %Zd) (mod %Zd)\n",
                     -score, i + 1, w, u, v, mod);
      }
    }
  }


  /* Step3, final root sieve */
  used = tmp_E_pqueue->used - 1;
  for (i = 0; i < used; i ++) {

    extract_MurphyE_pq (tmp_E_pqueue, &w, u, v, mod,
                        &score);

    if (param->verbose >= 2) {
      gmp_fprintf (stderr, "# Info: Sieve lat (# %2d), "
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
        fprintf (stderr, "# Info: Re-sieve lat (# %2d) in "
                 "range %u.\n", i + 1, size_tune_sievearray * 2);
      }
     
      ropt_s2param_setup_tune (s1param, s2param, u, v, mod,
                               0, size_tune_sievearray * 2, NP - 1);
      ropt_stage2 (poly, s2param, param, info, global_E_pqueue, w);
    }
  }

  /* free */
  free_MurphyE_pq (&tmp_E_pqueue);
  ropt_s2param_free (poly, s2param);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);
  mpz_clear (u_old);
  mpz_clear (v_old);
  mpz_clear (mod_old);
}


#if TUNE_LOGNORM_INCR
/**
 * Tune the BOUND_LOGNORM_INCR_MAX
 */
static double
ropt_linear_pre_tune ( ropt_poly_t poly,
                       ropt_param_t param,
                       alpha_pq *alpha_pqueue,
                       ropt_info_t info,
                       MurphyE_pq *global_E_pqueue )
{
  info->mode = 0;
  unsigned int pqueue_size = 32;
  unsigned int s1_size = 16;

  /* setup bound, pqueue */
  int i, r, w, used, k = -1;
  ropt_bound_t bound;
  alpha_pq *pqueue;
  ropt_s1param_t s1param;
  mpz_t u, v, mod;
  mpz_init (u);
  mpz_init (v);
  mpz_init (mod);
  ropt_bound_init (bound);
  new_alpha_pq (&pqueue, pqueue_size);
  ropt_s1param_init (s1param);
  double score, maxscore,
    incr = BOUND_LOGNORM_INCR_MAX,
    best_incr = BOUND_LOGNORM_INCR_MAX;
  maxscore = 0;
  while ((double)k<=param->effort) {

    /* change bound and test-sieve */
    ropt_bound_setup_incr (poly, bound, param, incr);
    ropt_s1param_resetup (poly, s1param, bound, param, s1_size);
    r = ropt_stage1 (poly, bound, s1param, param, pqueue, 0);
    if (r == -1) break;
    score = ropt_linear_tune_fast (poly, s1param, param, info, pqueue,
                                   global_E_pqueue, size_tune_sievearray);
    /* update maxscore */
    if (maxscore < score) {
      maxscore = score;
      best_incr = incr;
    }

    /* better not wasting sieved slots */
    used = pqueue->used - 1;
    for (i = 0; i < used; i ++) {
      extract_alpha_pq (pqueue, &w, u, v, mod, &score);
      insert_alpha_pq (alpha_pqueue, w, u, v, mod, score);
    }

    /* reset queue */
    reset_alpha_pq (pqueue);
    incr += BOUND_LOGNORM_INCR_MAX_TUNESTEP;
    k ++;
  }
  
  if (param->verbose >= 1) {
    gmp_fprintf (stderr, "# Info: tune (incr), best_lognorm_incr: %f\n",
                 best_incr);
  }

  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);
  free_alpha_pq (&pqueue);
  ropt_s1param_free (s1param);
  ropt_bound_free (bound);
  info->mode = 0;
  return best_incr;
}
#endif

/**
 * Linear rotation: deg 5.
 */
void
ropt_linear_deg5 ( ropt_poly_t poly,
                   ropt_bestpoly_t bestpoly,
                   ropt_param_t param,
                   ropt_info_t info )
{

  /* setup bound, s1param, alpha_pqueue, global_E_pqueue */
  int r, old_nbest_sl;
  double t1, t2, t3;
#if TUNE_LOGNORM_INCR
  double incr;
#endif
  ropt_bound_t bound;
  ropt_s1param_t s1param;
  alpha_pq *alpha_pqueue;
#if TUNE_LOGNORM_INCR
  alpha_pq *pre_alpha_pqueue;
#endif
  MurphyE_pq *global_E_pqueue;
  ropt_bound_init (bound);
  ropt_bound_setup (poly, bound, param, BOUND_LOGNORM_INCR_MAX);
  ropt_s1param_init (s1param);
  ropt_s1param_setup (poly, s1param, bound, param);
  new_MurphyE_pq (&global_E_pqueue, s1param->nbest_sl);
#if TUNE_LOGNORM_INCR
  new_alpha_pq (&pre_alpha_pqueue, s1param->nbest_sl);
#endif  
  new_alpha_pq (&alpha_pqueue, s1param->nbest_sl*TUNE_NUM_SUBLATTICE);

  /* Step 1: tune ropt_bound first */
#if TUNE_LOGNORM_INCR
  incr = ropt_linear_pre_tune (poly, param, pre_alpha_pqueue,
                               info, global_E_pqueue);
  ropt_bound_setup_incr (poly, bound, param, incr);
  ropt_s1param_resetup (poly, s1param, bound, param, s1param->nbest_sl);
#endif
  
  /* Step 2:, tuning to find good sublattices */
  t1 = seconds_thread ();
  old_nbest_sl = s1param->nbest_sl;
  s1param->nbest_sl = s1param->nbest_sl*TUNE_NUM_SUBLATTICE_STAGE1;
  r = ropt_stage1 (poly, bound, s1param, param, alpha_pqueue, 0);
  s1param->nbest_sl = old_nbest_sl;
  t1 = seconds_thread () - t1;
  if (r == -1) return;
  
  /* Step 3: rank/tune above found sublattices by short sieving */
  t2 = seconds_thread ();
#if TUNE_LOGNORM_INCR
  ropt_linear_tune (poly, bound, s1param, param, info, alpha_pqueue,
                    pre_alpha_pqueue, global_E_pqueue);
#else
  ropt_linear_tune (poly, bound, s1param, param, info, alpha_pqueue,
                    global_E_pqueue);
#endif
  t2 = seconds_thread () - t2;

  /* Step 3, root sieve */
  t3 = seconds_thread ();
  ropt_linear_sieve (poly, bound, s1param, param, info, alpha_pqueue,
                     global_E_pqueue);
  t3 = seconds_thread () - t3;

  info->ropt_time_stage1 = t1;
  info->ropt_time_tuning = t2;
  info->ropt_time_stage2 = t3;

  if (param->verbose >= 1) {
    fprintf ( stderr, "# Stat: this polynomial (stage 1) took %.2fs\n", t1 );
    fprintf ( stderr, "# Stat: this polynomial (tuning ) took %.2fs\n", t2 );
    fprintf ( stderr, "# Stat: this polynomial (stage 2) took %.2fs\n", t3 );
  }

  /* Step 4, return best poly */
  ropt_get_bestpoly (poly, global_E_pqueue, bestpoly);

  /* free */
  free_MurphyE_pq (&global_E_pqueue);
#if TUNE_LOGNORM_INCR
  free_alpha_pq (&pre_alpha_pqueue);
#endif
  free_alpha_pq (&alpha_pqueue);
  ropt_s1param_free (s1param);
  ropt_bound_free (bound);
}


/**
 * Linear rotation: deg 4.
 * No tuning at all.
 */
void
ropt_linear_deg34 ( ropt_poly_t poly,
                    ropt_bestpoly_t bestpoly,
                    ropt_param_t param,
                    ropt_info_t info )
{
  unsigned long ub, vb;
  ropt_bound_t bound;
  ropt_s1param_t s1param;
  ropt_s2param_t s2param;
  MurphyE_pq *global_E_pqueue;
  mpz_t u, v, mod;

  mpz_init_set_ui (u, 0);
  mpz_init_set_ui (v, 0);
  mpz_init_set_ui (mod, 1);
  
  /* setup bound, s1param */
  ropt_bound_init (bound);
  ropt_bound_setup (poly, bound, param, BOUND_LOGNORM_INCR_MAX);
  ropt_s1param_init (s1param);
  ropt_s1param_setup (poly, s1param, bound, param);
  new_MurphyE_pq (&global_E_pqueue, s1param->nbest_sl);
  ropt_s2param_init (poly, s2param);
  
  /* Step 1, set up lattice (0, 0) (mod 1) and run root sieve */
  ub = (bound->global_u_boundr > 0) ?
    bound->global_u_boundr : (-bound->global_u_boundr);
  vb = mpz_get_ui (bound->global_v_boundr);
  if (vb < 128) vb = 128;
  
  ropt_s2param_setup_tune (s1param, s2param, u, v, mod,
                           ub, vb, NP - 1);

  ropt_stage2 (poly, s2param, param, info, global_E_pqueue, 0);

  /* Step 2, return best poly */
  ropt_get_bestpoly (poly, global_E_pqueue, bestpoly);

  /* free */
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (mod);  
  free_MurphyE_pq (&global_E_pqueue);
  ropt_s1param_free (s1param);
  ropt_s2param_free (poly, s2param);
  ropt_bound_free (bound);
}


/**
 * Linear rotation: for deg 4 and deg 5
 */
void
ropt_linear ( ropt_poly_t poly,
              ropt_bestpoly_t bestpoly,
              ropt_param_t param,
              ropt_info_t info )
{
  if (poly->d == 3)
    ropt_linear_deg34 (poly, bestpoly, param, info);
  else if (poly->d == 5 || poly->d == 4)
    ropt_linear_deg5 (poly, bestpoly, param, info);
  else 
    fprintf (stderr, "Error: ropt_linear() only support deg 4 and 5.");
}
