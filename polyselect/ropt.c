/**
 * @file ropt.c
 * Called by polyselect2l.c or ropt_main.c and then call stage1.c and stage2.c.
 */

#include "ropt.h"

/* 
   Record best poly
*/
static inline void
ropt_main_run_bestpoly ( rsstr_t rs,
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
    fprintf (stderr, "Error, cannot allocate memory in ropt_main_run_bestpoly().\n");
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
  Run: find good sublattices and then sieve.
*/
static inline double
ropt_main_run ( rsstr_t rs,
                bestpoly_t bestpoly,
                param_t param,
                int verbose )
{
  int i, k, old_i, re, used;
  double ave_MurphyE, ave2_MurphyE = 0.0;
  mpz_t m;
  rsparam_t rsparam;
  sub_alpha_pq *alpha_pqueue;

  mpz_init_set (m, rs->g[0]);
  mpz_neg (m, m);
  rsparam_init (rsparam, rs, param);
  rsparam_setup (rsparam, rs, param, verbose);
  new_sub_alpha_pq (&alpha_pqueue, rsparam->nbest_sl);

  /* read e_sl from input */
  if (param->flag == 1) {
    for (i = 0; i < LEN_SUBLATTICE_PRIMES; i ++)
      rsparam->e_sl[i] = param->s1_e_sl[i];
  }

  /* For each qudratic rotation i, find sublattice */
  re = 0;
  old_i = 0;
  for (i = param->w_left_bound; i < param->w_length + param->w_left_bound; i++) {

    if (verbose == 2)
      fprintf (stderr, "# Info: quadratic rotation by %d*x^2\n", i);

    old_i = rotate_aux (rs->f, rs->g[1], m, old_i, i, 2);
    rsstr_setup (rs);

    verbose = 0;
    /* either use input e_sl[] or tune only once */
    rsparam_reset_bounds (rsparam, rs, param, verbose);
    verbose = 2;
      
#if WANT_TUNE
    if (re == 0) {
      /* do we really want to tune? slow */
      rsparam_tune (rs, rsparam, param, 6, i, verbose);
    }
#endif

    re ++;
    k = ropt_stage1 ( rs,
                      rsparam,
                      alpha_pqueue,
                      verbose,
                      i );
    if (k == -1) {
      /* rotate back */
      rotate_aux (rs->f, rs->g[1], m, old_i, 0, 2);
      old_i = 0;
      rsparam_free (rsparam);
      mpz_clear (m);
      free_sub_alpha_pq (&alpha_pqueue);
      return -1;
    }
  }

  /* rotate back */
  rotate_aux (rs->f, rs->g[1], m, old_i, 0, 2);
  old_i = 0;

  /* use another array to save the priority queue */
  used = alpha_pqueue->used - 1;
  mpz_t *u, *v, *mod;
  int *w;
  double *sub_alpha;
  w = (int *) malloc ( used * sizeof (int));
  sub_alpha = (double *) malloc ( used * sizeof (double));
  u = (mpz_t *) malloc ( used * sizeof (mpz_t));
  v = (mpz_t *) malloc ( used * sizeof (mpz_t));
  mod = (mpz_t *) malloc ( used * sizeof (mpz_t));
  for (i = 0; i < used; i ++) {
    mpz_init (u[i]);
    mpz_init (v[i]);
    mpz_init (mod[i]);
  }

  /* put all sublattices into another array */
  for (i = 0; i < used; i ++) {

    extract_sub_alpha_pq ( alpha_pqueue,
                           &(w[i]),
                           u[i],
                           v[i],
                           mod[i],
                           &(sub_alpha[i]) );

    if (verbose != 0) {
      gmp_fprintf ( stderr, "# Info: %4d sublattice (w, u, v): (%d, %Zd, %Zd) (mod %Zd), alpha: %.2f\n",
                    i + 1,
                    w[i],
                    u[i],
                    v[i],
                    mod[i],
                    sub_alpha[i] );

      // This is only for pbs submission purpose.
      // Note: 1. Change N and path (appeared twice)
      //       2. Change sieve length -umax and -vmax
      /*
        fprintf (stderr, "#!/bin/bash\n#PBS -N rsa_%d\n#PBS -l nodes=1,walltime=160:00:00\n#PBS -q route\n#PBS -m ae\n", i);
        gmp_fprintf ( stderr, "/home/bai/cado-nfs/build/orac/polyselect/ropt -fm /home/bai/cado-nfs/build/orac/polyselect/rsa768.poly --s2 -n %Zd -d %d -w %d -u %Zd -v %Zd -umax 64 -vmax 317325312 -mod %Zd > /home/bai/cado-nfs/build/orac/polyselect/rsa704_%d.out 2>&1\n\n",
        rs->n,
        rs->d,
        w[i],
        u[i],
        v[i],
        mod[i],
        i );
      */
    }
  }
  /* free alpha queue */
  free_sub_alpha_pq (&alpha_pqueue);

  /* E priority queue for all sublattice, only consider the top three polynomials */
  MurphyE_pq *global_E_pqueue;
  new_MurphyE_pq (&global_E_pqueue, 4);

  /* For each sublattice, do the root sieve */
  re = 0;
  for (i = used - 1; i >= 0; i --) {

    /* for polyselect2.c */
    if (verbose == 0) {
      if (re > 3)
        break;
    }

    /* don't output in tune mode */
    if (verbose != 0) {
      gmp_fprintf ( stderr,
                    "\n# Info: Sieve on sublattice (# %2d), (w, u, v): (%d, %Zd, %Zd) (mod %Zd) \n# Info: alpha: %.2f, proj_alpha: %.2f, exp_min_alpha: %.2f\n",
                    i + 1,
                    w[i],
                    u[i],
                    v[i],
                    mod[i],
                    sub_alpha[i],
                    rs->alpha_proj,
                    rsparam->exp_min_alpha_rs );
    }

    /* rotate polynomial by f + rot*x^2 for various rot */
    old_i = rotate_aux (rs->f, rs->g[1], m, old_i, w[i], 2);
    rsstr_setup (rs);
    rsparam_reset_bounds (rsparam, rs, param, verbose);

    //print_poly_fg (rs->f, rs->g, rs->d, rs->n, rs->m, 1);
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

  /* rotate back */
  rotate_aux (rs->f, rs->g[1], m, old_i, 0, 2);
  old_i = 0;
  rsparam_free (rsparam);

  /* Record the best polynomial */
  if (verbose == 0 || verbose == 2) {
    ropt_main_run_bestpoly ( rs,
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
  free (sub_alpha);

  ave2_MurphyE /= (double) re;
  return ave2_MurphyE;
}


/*
  sieve only
*/
static inline double
ropt_main_run_stage2only ( rsstr_t rs,
                           param_t param )
{
  int i, old_i;
  double ave_MurphyE, alpha_lat;
  mpz_t m,  *fuv, *guv;
  rsparam_t rsparam;

  mpz_init_set (m, rs->g[0]);
  mpz_neg (m, m);
  fuv = (mpz_t*) malloc ((rs->d + 1) * sizeof (mpz_t));
  guv = (mpz_t*) malloc (2 * sizeof (mpz_t));
  if (fuv == NULL || guv == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in ropt_main_run_stage2only().\n");
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

  ropt_main_run_bestpoly ( rs,
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
  For the current polynomial (rs), start two-stage root sieve.
*/
void
ropt_main ( rsstr_t rs,
            bestpoly_t bestpoly,
            param_t param,
            int verbose )
{

  /* stage 2 (sieve) only */
  if (param->flag == 2) {
    if (rs->d == 5) {
      fprintf (stderr, "Error: sieve-only mode support deg 6 only.\n");
      exit(1);
    }
    ropt_main_run_stage2only (rs, param);
  }
  else {
    if (rs->d == 5 || rs->d == 4) {
      /* not qudratci rot for deg 5 polynomial */
      param->w_left_bound = 0;
      param->w_length = 1;
      ropt_main_run (rs, bestpoly, param, verbose);
    }
    else if (rs->d == 6) {
      ropt_main_run (rs, bestpoly, param, verbose);
    }
    else {
      fprintf (stderr, "Error: only support deg 4, 5 or 6.\n");
      exit(1);
    }
  }
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
  rsstr_t rs;
  param_t param;

  rsstr_init (rs);
  mpz_set (rs->g[1], l);
  mpz_neg (rs->g[0], m);
  for (i = 0; i <=d; i ++)
    mpz_set (rs->f[i], f[i]);
  mpz_set (rs->n, N);
  rsstr_setup (rs);
  param_init (param);

  /* this should be fast to tune */
  if (d == 6) {
    param->w_left_bound = 1;
    param->w_length = 1;
  }

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
  ropt_main (rs, bestpoly, param, verbose);

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
