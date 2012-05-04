/**
 * @file ropt_quadratic.c
 * Main function for quadratic rotation.
 * Called by ropt.c and then call stage1.c and stage2.c.
 */

#include "ropt_quadratic.h"
#include "ropt.h"

/*
  Find good sublattices and then linear rotation.
*/
double
ropt_quadratic ( ropt_poly_t rs,
                      bestpoly_t bestpoly,
                      param_t param,
                      int verbose )
{
  int i, k, old_i, re, used, *w;
  double ave_MurphyE, ave2_MurphyE = 0.0, *score;
  mpz_t m, *u, *v, *mod;
  ropt_param_t rsparam;
  sub_alpha_pq *alpha_pqueue;

  mpz_init_set (m, rs->g[0]);
  mpz_neg (m, m);
  rsparam_init (rsparam, rs, param);
  rsparam_setup (rsparam, rs, param, verbose);
  new_sub_alpha_pq (&alpha_pqueue, rsparam->nbest_sl);

  // if needed, read e_sl from input.
  if (param->flag == 1) {
    for (i = 0; i < LEN_SUBLATTICE_PRIMES; i ++)
      rsparam->e_sl[i] = param->s1_e_sl[i];
  }
  
  /* STAGE 1: for each quadratic rotation i, find sublattice */
  re = 0;
  old_i = 0;
  for (i = param->w_left_bound; i < param->w_length + param->w_left_bound; i++) {

    if (verbose == 2)
      fprintf (stderr, "# Info: quadratic rotation by %d*x^2\n", i);

    old_i = rotate_aux (rs->f, rs->g[1], m, old_i, i, 2);
    rsstr_setup (rs);

    // either use input e_sl[] or tune only once.
    rsparam_reset_bounds (rsparam, rs, param, verbose);
    
#if TUNE_FIND_SUBLATTICE
    if (re == 0) {
      // tune sublattice parameters?
      ropt_param_tune_findlat (rs, rsparam, param, 6, i, verbose);
    }
#endif

    re ++;
    
    // stage 1. rs is already rotated, i here is only for book-keeping purpose.
    k = ropt_stage1 ( rs,
                      rsparam,
                      alpha_pqueue,
                      verbose,
                      i );
    if (k == -1) {
      // rotate back
      rotate_aux (rs->f, rs->g[1], m, old_i, 0, 2);
      old_i = 0;
      rsparam_free (rsparam);
      mpz_clear (m);
      free_sub_alpha_pq (&alpha_pqueue);
      return -1;
    }
  }
  // rotate back.
  rotate_aux (rs->f, rs->g[1], m, old_i, 0, 2);
  old_i = 0;

  // use another array to save the priority queue.
  used = alpha_pqueue->used - 1;

  w = (int *) malloc ( used * sizeof (int) );
  score = (double *) malloc ( used * sizeof (double) );
  u = (mpz_t *) malloc ( used * sizeof (mpz_t) );
  v = (mpz_t *) malloc ( used * sizeof (mpz_t) );
  mod = (mpz_t *) malloc ( used * sizeof (mpz_t) );
  for (i = 0; i < used; i ++) {
    mpz_init (u[i]);
    mpz_init (v[i]);
    mpz_init (mod[i]);
  }

/* --------------------------------------- */
/* rank found sublattices by test sieving. */
/* --------------------------------------- */
#if TUNE_RANK_SUBLATTICE

  sub_alpha_pq *tsieve_MurphyE_pqueue;
  new_sub_alpha_pq ( &tsieve_MurphyE_pqueue,
                     rsparam->nbest_sl );
  ropt_param_tune_ranklat ( rs,
                         rsparam,
                         param,
                         alpha_pqueue,
                         used,
                         tsieve_MurphyE_pqueue,
                         verbose );
  for (i = 0; i < used; i ++) {
    // put all sublattices into another array.
    extract_sub_alpha_pq ( tsieve_MurphyE_pqueue,
                           &(w[i]),
                           u[i],
                           v[i],
                           mod[i],
                           &(score[i]) );
    if (verbose == 2) {
      gmp_fprintf ( stderr, "# Info: %4d sublattice (w, u, v): (%d, %Zd, %Zd)"
                    " (mod %Zd), tsieve E: %.2e\n",
                    i + 1,
                    w[i],
                    u[i],
                    v[i],
                    mod[i],
                    -score[i] );
    }
  }
  free_sub_alpha_pq (&tsieve_MurphyE_pqueue);

  // re-detect sublattices with the best quadratic roation.
  if (param->w_length > 1) {
    
    i = used - 1;
    old_i = 0;
    if (verbose == 2)
      fprintf (stderr, "# Info: Found best quadratic rotation by %d*x^2\n", w[i]);
    old_i = rotate_aux (rs->f, rs->g[1], m, old_i, w[i], 2);
    rsstr_setup (rs);
    rsparam_reset_bounds (rsparam, rs, param, verbose);
    k = ropt_stage1 ( rs,
                      rsparam,
                      alpha_pqueue,
                      verbose,
                      w[i] );
    // if failed in ropt_stage1, return.
    if (k == -1) {
      rotate_aux (rs->f, rs->g[1], m, old_i, 0, 2);
      old_i = 0;
      rsparam_free (rsparam);
      mpz_clear (m);
      free_sub_alpha_pq (&alpha_pqueue);
      return -1;
    }
    // rotate back.
    rotate_aux (rs->f, rs->g[1], m, old_i, 0, 2);

    old_i = 0;
    // use another array to save the priority queue.
    used = alpha_pqueue->used - 1;
  }
#else

  /* ----------------------------------------------- */
  /* rank found sublattices by partial alpha values. */
  /* ----------------------------------------------- */
  for (i = 0; i < used; i ++) {
    // put all sublattices into another array.
    extract_sub_alpha_pq ( alpha_pqueue,
                           &(w[i]),
                           u[i],
                           v[i],
                           mod[i],
                           &(score[i]));

    if (verbose == 2) {
      gmp_fprintf ( stderr, "# Info: %4d sublattice (w, u, v): (%d, %Zd, %Zd) (mod %Zd), alpha: %.2f\n",
                    i + 1,
                    w[i],
                    u[i],
                    v[i],
                    mod[i],
                    score[i] );
    }
  }
  free_sub_alpha_pq (&alpha_pqueue); // free alpha queue.
#endif


  // E priority queue for all sublattice, we rank
  // the top three polynomials from each sublattice.
  MurphyE_pq *global_E_pqueue;
  new_MurphyE_pq (&global_E_pqueue, 4);

  /* STAGE2: for each sublattice, do the root sieve */
  re = 0;
  for (i = used - 1; i >= 0; i --) {

    /* for polyselect2.c, we don't want to spent too much time. */
    if (verbose == 0) {
      if (re > SHORT_NUM_SIEVE_SUBLATTICE)
        break;
    }

    /* don't output in tune mode */
    if (verbose == 2) {
      gmp_fprintf ( stderr,
                    "\n# Info: Sieve on sublattice (# %2d), (w, u, v): (%d, %Zd, %Zd) (mod %Zd) \n# Info: alpha: %.2f, proj_alpha: %.2f, exp_min_alpha: %.2f\n",
                    i + 1,
                    w[i],
                    u[i],
                    v[i],
                    mod[i],
                    score[i],
                    rs->alpha_proj,
                    rsparam->exp_min_alpha_rs );
    }

    // rotate polynomial by f + rot*x^2 for various rot.
    old_i = rotate_aux (rs->f, rs->g[1], m, old_i, w[i], 2);
    rsstr_setup (rs);

    rsparam_reset_bounds (rsparam, rs, param, verbose);

    // print_poly_fg (rs->f, rs->g, rs->d, rs->n, rs->m, 1);
    ave_MurphyE = ropt_stage2 ( rs,
                                rsparam,
                                param,
                                global_E_pqueue,
                                w[i],
                                u[i],
                                v[i],
                                mod[i],
                                verbose );
    ave2_MurphyE += ave_MurphyE;
    re ++;
  }

  // rotate back.
  rotate_aux (rs->f, rs->g[1], m, old_i, 0, 2);
  old_i = 0;
  rsparam_free (rsparam);

  // record the best polynomial.
  if (verbose == 0 || verbose == 2) {
    ropt_return_bestpoly ( rs,
                           global_E_pqueue,
                           bestpoly );

    if (verbose == 2) {
      fprintf (stderr, "\n# Info: Best E is:\n");
      print_poly_fg (bestpoly->f, bestpoly->g, rs->d, rs->n, 2);
    }
  }

  free_MurphyE_pq (&global_E_pqueue);
  mpz_clear (m);

  for (i = 0; i < used; i ++) {
    mpz_clear (u[i]);
    mpz_clear (v[i]);
    mpz_clear (mod[i]);
  }
  free (u);
  free (v);
  free (mod);
  free (w);
  free (score);

  ave2_MurphyE /= (double) re;
  return ave2_MurphyE;
}
