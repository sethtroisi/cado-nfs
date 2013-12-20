/**
 * @file ropt.c
 * This is the main interface, which is called from
 * either polyselect2*.c or ropt.c. There are three
 * main functions:
 * - ropt_do_both_stages()
 * - ropt_do_stage2()
 * - ropt_polyselect()
 */


#include "cado.h"
#include "ropt.h"
#include "portability.h"


/**
 * Find best poly. This is somehow redundant.
 */
void
ropt_get_bestpoly ( ropt_poly_t poly,
                    MurphyE_pq *global_E_pqueue,
                    ropt_bestpoly_t bestpoly )
{
  double ave_MurphyE = 0.0, best_E = 0.0;
  int i, old_i, k;
  mpz_t m, *fuv, *guv;
  mpz_poly_t Fuv;

  mpz_init_set (m, poly->g[0]);
  mpz_neg (m, m);

  /* var for computing E */
  mpz_poly_init (Fuv, poly->d);
  Fuv->deg = poly->d;
  fuv = Fuv->coeff;
  guv = (mpz_t*) malloc (2 * sizeof (mpz_t));
  if (guv == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in "
             "ropt_get_bestpoly().\n");
    exit (1);
  }
  for (i = 0; i <= poly->d; i++)
    mpz_set (fuv[i], poly->f[i]);
  for (i = 0; i < 2; i++)
    mpz_init_set (guv[i], poly->g[i]);

  /* output all polys in the global queue */
  for (i = 1; i < global_E_pqueue->used; i ++) {

    old_i = 0;
    old_i = rotate_aux (poly->f, poly->g[1], m, old_i,
                        global_E_pqueue->w[i], 2);
    for (k = 0; k <= poly->d; k++)
      mpz_set (fuv[k], poly->f[k]);

    for (k = 0; k < 2; k++)
      mpz_set (guv[k], poly->g[k]);

    compute_fuv_mp (fuv, poly->f, poly->g, poly->d,
                    global_E_pqueue->u[i], global_E_pqueue->v[i]);
    optimize_aux (Fuv, guv, 0, 0, CIRCULAR);
    ave_MurphyE = print_poly_fg (Fuv, guv, poly->n, 0);

    if (ave_MurphyE > best_E) {
      best_E = ave_MurphyE;
      for (k = 0; k <= poly->d; k++)
        mpz_set (bestpoly->f[k], fuv[k]);
      for (k = 0; k < 2; k++)
        mpz_set (bestpoly->g[k], guv[k]);
    }
    rotate_aux (poly->f, poly->g[1], m, old_i, 0, 2);
  }


  mpz_poly_clear (Fuv);
  for (i = 0; i < 2; i++)
    mpz_clear (guv[i]);
  mpz_clear (m);
  free (guv);
}


/**
 * Root sieve only.
 */
static void
ropt_do_stage2 (ropt_poly_t poly,
                ropt_bestpoly_t bestpoly,
                ropt_param_t param,
                ropt_info_t info )
{

  int i, old_i;
  double alpha_lat;
  mpz_t m,  *fuv, *guv;
  ropt_bound_t bound;
  ropt_s2param_t s2param;
  MurphyE_pq *global_E_pqueue;
  mpz_poly_t Fuv;

  ropt_bound_init (bound);
  new_MurphyE_pq (&global_E_pqueue, 4);

  mpz_init_set (m, poly->g[0]);
  mpz_neg (m, m);
  mpz_poly_init (Fuv, poly->d);
  Fuv->deg = poly->d;
  fuv = Fuv->coeff;
  guv = (mpz_t*) malloc (2 * sizeof (mpz_t));
  if (guv == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in "
             "ropt_do_stage2().\n");
    exit (1);
  }
  for (i = 0; i <= poly->d; i++)
    mpz_set (fuv[i], poly->f[i]);
  for (i = 0; i < 2; i++)
    mpz_init_set (guv[i], poly->g[i]);

  /* rotate polynomial by f + rot*x^2 */
  old_i = 0;
  old_i = rotate_aux (poly->f, poly->g[1], m, old_i, param->s2_w, 2);
  //print_poly_fg (poly->f, poly->g, poly->d, poly->n, 1);

  /* reset after rotation */
  ropt_poly_setup (poly);
  ropt_bound_setup (poly, bound, param);

  /* print some basic information */
  compute_fuv_mp (fuv, poly->f, poly->g, poly->d, param->s2_u, param->s2_v);
  alpha_lat = get_alpha (Fuv, 2000);
  gmp_fprintf ( stderr,
                "\n# Info: Sieve on sublattice, (w, u, v): (%d, %Zd, %Zd) "
                "(mod %Zd)\n"
                "# Info: alpha: %.2f, proj_alpha: %.2f, "
                "exp_alpha: %.2f\n",
                param->s2_w,
                param->s2_u,
                param->s2_v,
                param->s2_mod,
                alpha_lat,
                poly->alpha_proj,
                bound->exp_min_alpha );

  /* root sieve */
  ropt_s2param_init (poly, s2param);
  ropt_s2param_setup_stage2_only (bound, s2param, param,
                                  param->s2_u, param->s2_v, param->s2_mod);
  info->w = param->s2_w;
  ropt_stage2 (poly, s2param, param, info, global_E_pqueue, param->s2_w);

  /* rotate back */
  rotate_aux (poly->f, poly->g[1], m, old_i, 0, 2);
  old_i = 0;

  /* return best poly */
  ropt_get_bestpoly (poly, global_E_pqueue, bestpoly);

  /* free */
  free_MurphyE_pq (&global_E_pqueue);
  ropt_bound_free (bound);
  ropt_s2param_free (poly, s2param);
  mpz_clear (m);
  mpz_poly_clear (Fuv);
  for (i = 0; i < 2; i++)
    mpz_clear (guv[i]);
  free (guv);

}


/**
 * Ropt linear or quadratic.
 */
static void
ropt_do_both_stages ( ropt_poly_t poly,
                      ropt_bestpoly_t bestpoly,
                      ropt_param_t param,
                      ropt_info_t info )
{
  if (poly->d == 5 || poly->d == 4 || poly->d == 3)
    ropt_linear (poly, bestpoly, param, info);
  else if (poly->d == 6)
    ropt_quadratic (poly, bestpoly, param, info);
  else {
    fprintf (stderr, "Error: only support deg 4, 5 or 6.\n");
    exit(1);
  }
}


/**
 * Start the two-stage root optimization on poly.
 * @param poly contains polynomial.
 * @param param input parameters.
 * @return top polynomial in bestpoly.
 */
void
ropt ( ropt_poly_t poly,
       ropt_bestpoly_t bestpoly,
       ropt_param_t param,
       ropt_info_t info )
{

  /* print cache size */
  if (param->verbose == 2)
    fprintf ( stderr, "# Info: L1_cachesize/2: %d, "
              "size_tune_sievearray: %d\n",
              L1_cachesize, size_tune_sievearray );

  if (param->stage_flag == 2)
    ropt_do_stage2 (poly, bestpoly, param, info);
  else
    ropt_do_both_stages (poly, bestpoly, param, info);

}


/**
 * Called by polyselect2 or polyselect2l.
 */
void
ropt_polyselect ( mpz_t *f,
                  int d,
                  mpz_t m,
                  mpz_t l,
                  mpz_t N ,
                  int verbose )
{
  int i;
  ropt_poly_t poly;
  ropt_poly_init (poly);

  /* setup poly */
  mpz_set (poly->g[1], l);
  mpz_neg (poly->g[0], m);
  for (i = 0; i <=d; i ++)
    mpz_set (poly->f[i], f[i]);
  mpz_set (poly->n, N);
  ropt_poly_setup (poly);

  ropt_info_t info;
  ropt_info_init (info);

  ropt_param_t param;
  ropt_param_init (param);
  param->verbose = verbose;

  ropt_bestpoly_t bestpoly;
  ropt_bestpoly_init (bestpoly, poly->d);
  ropt_bestpoly_setup (bestpoly, poly->f, poly->g, poly->d);

  /* cal main function */
  ropt_do_both_stages (poly, bestpoly, param, info);
  
  /* bring bestpoly back to polyselect2* */
  for (i = 0; i <= d; i++)
    mpz_set (f[i], bestpoly->f[i]);
  mpz_neg (m, bestpoly->g[0]);
  mpz_set (l, bestpoly->g[1]);

  /* free */
  ropt_param_free (param);
  ropt_bestpoly_free (bestpoly, poly->d);
  ropt_info_free (info);
  ropt_poly_free (poly);
}
