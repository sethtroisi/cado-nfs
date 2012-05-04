/**
 * @file ropt.c
 * Called by polyselect2l.c or ropt_main.c and then call ropt_small.c or ropt_large.c
 */

#include "ropt.h"


/**
 * Find best poly. This is somehow redundant.
 */
void
ropt_return_bestpoly ( ropt_poly_t rs,
                       MurphyE_pq *global_E_pqueue,
                       bestpoly_t bestpoly )
{
  double ave_MurphyE = 0.0, best_E = 0.0;
  int i, old_i, k;
  mpz_t m, *fuv, *guv;

  mpz_init_set (m, rs->g[0]);
  mpz_neg (m, m);

  /* var for computing E */
  fuv = (mpz_t*) malloc ((rs->d + 1) * sizeof (mpz_t));
  guv = (mpz_t*) malloc (2 * sizeof (mpz_t));
  if (fuv == NULL || guv == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in ropt_return_bestpoly().\n");
    exit (1);
  }
  for (i = 0; i <= rs->d; i++)
    mpz_init_set (fuv[i], rs->f[i]);
  for (i = 0; i < 2; i++)
    mpz_init_set (guv[i], rs->g[i]);

  /* output all polys in the global queue */
  for (i = 1; i < global_E_pqueue->used; i ++) {

    old_i = 0;
    old_i = rotate_aux (rs->f, rs->g[1], m, old_i, global_E_pqueue->w[i], 2);
    for (k = 0; k <= rs->d; k++)
      mpz_set (fuv[k], rs->f[k]);

    for (k = 0; k < 2; k++)
      mpz_set (guv[k], rs->g[k]);

    compute_fuv_mp (fuv, rs->f, rs->g, rs->d, global_E_pqueue->u[i], global_E_pqueue->v[i]);
    optimize_aux (fuv, rs->d, guv, 0, 0, CIRCULAR);
    ave_MurphyE = print_poly_fg (fuv, guv, rs->d, rs->n, 0); // only output when verbose == 2.

    if (ave_MurphyE > best_E) {
      best_E = ave_MurphyE;
      for (k = 0; k <= rs->d; k++)
        mpz_set (bestpoly->f[k], fuv[k]);
      for (k = 0; k < 2; k++)
        mpz_set (bestpoly->g[k], guv[k]);
    }
    rotate_aux (rs->f, rs->g[1], m, old_i, 0, 2);
  }

  for (i = 0; i <= rs->d; i++)
    mpz_clear (fuv[i]);
  for (i = 0; i < 2; i++)
    mpz_clear (guv[i]);
  mpz_clear (m);
  free (fuv);
  free (guv);
}


/*
  sieve only
*/
static inline double
ropt_run_stage2only ( ropt_poly_t rs,
                      param_t param )
{
  int i, old_i;
  double ave_MurphyE, alpha_lat;
  mpz_t m,  *fuv, *guv;
  ropt_param_t rsparam;

  mpz_init_set (m, rs->g[0]);
  mpz_neg (m, m);
  fuv = (mpz_t*) malloc ((rs->d + 1) * sizeof (mpz_t));
  guv = (mpz_t*) malloc (2 * sizeof (mpz_t));
  if (fuv == NULL || guv == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in ropt_run_stage2only().\n");
    exit (1);
  }
  for (i = 0; i <= rs->d; i++)
    mpz_init_set (fuv[i], rs->f[i]);
  for (i = 0; i < 2; i++)
    mpz_init_set (guv[i], rs->g[i]);

  rsparam_init (rsparam, rs, param);
  rsparam_setup (rsparam, rs, param, 2);

  /* rotate polynomial by f + rot*x^2 */
  old_i = 0;
  old_i = rotate_aux (rs->f, rs->g[1], m, old_i, param->s2_w, 2);
  rsstr_setup (rs);
  rsparam_reset_bounds (rsparam, rs, param, 2);

  /* alpha values on the subllatice primes */
  compute_fuv_mp (fuv, rs->f, rs->g, rs->d, param->s2_u, param->s2_v);
  alpha_lat = get_alpha (fuv, rs->d, 2000);
  gmp_fprintf ( stderr,
                "\n# Info: Sieve on sublattice, (w, u, v): (%d, %Zd, %Zd) (mod %Zd) \n# Info: alpha: %.2f, proj_alpha: %.2f, exp_min_alpha: %.2f\n",
                param->s2_w,
                param->s2_u,
                param->s2_v,
                param->s2_mod,
                alpha_lat,
                rs->alpha_proj,
                rsparam->exp_min_alpha_rs );

  //print_poly_fg (rs->f, rs->g, rs->d, rs->n, 2);
  MurphyE_pq *global_E_pqueue;
  new_MurphyE_pq (&global_E_pqueue, 4);

  ave_MurphyE = ropt_stage2 ( rs,
                              rsparam,
                              param,
                              global_E_pqueue,
                              param->s2_w,
                              param->s2_u,
                              param->s2_v,
                              param->s2_mod,
                              2 );

  /* rotate back */
  rotate_aux (rs->f, rs->g[1], m, old_i, 0, 2);
  old_i = 0;

  /* Record the best polynomial */
  bestpoly_t bestpoly;
  bestpoly_init (bestpoly, rs->d);
  for (i = 0; i <= rs->d; i++)
  {
    mpz_set (bestpoly->f[i], rs->f[i]);
  }
  mpz_neg (bestpoly->g[0], rs->g[0]);
  mpz_set (bestpoly->g[1], rs->g[1]);

  ropt_return_bestpoly ( rs,
                           global_E_pqueue,
                           bestpoly );

  fprintf (stderr, "\n# Info: Best E is:\n");
  print_poly_fg (bestpoly->f, bestpoly->g, rs->d, rs->n, 2);

  bestpoly_free (bestpoly, rs->d);
  rsparam_free (rsparam);
  free_MurphyE_pq (&global_E_pqueue);

  mpz_clear (m);
  for (i = 0; i <= rs->d; i++)
    mpz_clear (fuv[i]);
  for (i = 0; i < 2; i++)
    mpz_clear (guv[i]);
  free (fuv);
  free (guv);

  return ave_MurphyE;
}


/*
  Called by polyselect2 or polyselect2l.
*/
void
ropt_polyselect ( mpz_t *f,
                  int d,
                  mpz_t m,
                  mpz_t l,
                  mpz_t N,
                  int max_k,
                  int verbose )
{
  /* rootsieve_struct */
  int i;
  ropt_poly_t rs;
  param_t param;

  rsstr_init (rs);
  mpz_set (rs->g[1], l);
  mpz_neg (rs->g[0], m);
  for (i = 0; i <=d; i ++)
    mpz_set (rs->f[i], f[i]);
  mpz_set (rs->n, N);
  rsstr_setup (rs);
  param_init (param);

  /* bestpoly for return */
  bestpoly_t bestpoly;
  bestpoly_init (bestpoly, d);
  for (i = 0; i <= d; i++)
  {
    mpz_set (bestpoly->f[i], f[i]);
  }
  mpz_neg (bestpoly->g[0], m);
  mpz_set (bestpoly->g[1], l);

  /* auto-mode */
  param->flag = 0;
  /* upper bound for sieve */
  param->s2_Vmax = max_k;

  /* start main rootsieve function */
  ropt (rs, bestpoly, param, verbose);

  /* set and free */
  for (i = 0; i <= d; i++) {
    mpz_set (f[i], bestpoly->f[i]);
  }

  mpz_neg (m, bestpoly->g[0]);
  mpz_set (l, bestpoly->g[1]);

  bestpoly_free (bestpoly, d);
  param_clear (param);
  rsstr_free (rs);
}


/**
 * For the polynomial in rs, start the two-stage root optimization.
 */
void
ropt ( ropt_poly_t rs,
       bestpoly_t bestpoly,
       param_t param,
       int verbose )
{
  /* L1 cache size */
  int ret = cachesize_cpuid (0);
  if ( (2048 <= ret) && (ret <= (1 << 20)) ) {
    L1_CACHESIZE = ret / 2;
    TUNE_SIEVEARRAY_SIZE = L1_CACHESIZE / 2;
  }
  else {
    ret = cachesize_guess (0);
    if ( (2048 <= ret)  && (ret <= (1 << 20)) ) {
      L1_CACHESIZE = ret / 2;
      TUNE_SIEVEARRAY_SIZE = L1_CACHESIZE / 2;
    }
  }

  if (verbose == 2)
    fprintf ( stderr, "# Info: L1_CACHESIZE/2: %d, TUNE_SIEVEARRAY_SIZE: %d.\n",
              L1_CACHESIZE, TUNE_SIEVEARRAY_SIZE );
  // finishring cache detection

  
  /* stage 2 (sieve) only */
  if (param->flag == 2) {
    ropt_run_stage2only (rs, param);
  }
  else {

    if (rs->d == 5 || rs->d == 4) {
      /* not quadratic rot for deg 5 polynomial */
      param->w_left_bound = 0;
      param->w_length = 1;
      ropt_linear (rs, bestpoly, param, verbose);
    }
    else if (rs->d == 6) {
      ropt_quadratic (rs, bestpoly, param, verbose);
    }
    else {
      fprintf (stderr, "Error: only support deg 4, 5 or 6.\n");
      exit(1);
    }
  }
}

