/**
 * @file ropt_str.c
 * Basic structs used in ropt.
 */


#include "cado.h"
#include "ropt_str.h"
#include "portability.h"


/* -----------------*/
/* Static Functions */
/* -----------------*/


/**
 * Replace f + k0 * x^t * (b*x - m) by f + k * x^t * (b*x - m), 
 * and return k to k0
 */
static inline void
rotate_aux_mpz ( mpz_t *f,
                 mpz_t b,
                 mpz_t m,
                 mpz_t k0,
                 mpz_t k,
                 unsigned int t )
{
  mpz_t tmp;
  mpz_init (tmp);
  mpz_sub (tmp, k, k0);
  mpz_addmul (f[t + 1], b, tmp);
  mpz_submul (f[t], m, tmp);
  mpz_set (k0, k);
  mpz_clear (tmp);
}


/**
 * Find bound V for constant rotation.
 */
static inline void
rotate_bounds_V_mpz ( ropt_poly_t poly,
                      ropt_bound_t bound )
{

  int i, j;
  double skewness, lognorm;
  mpz_t *f, *g, b, m, tmpv, V;

  f = (mpz_t *) malloc ( (poly->d + 1)* sizeof (mpz_t));
  g = (mpz_t *) malloc ( 2 * sizeof (mpz_t));
  for (j = 0; j <= poly->d; j ++)
    mpz_init_set (f[j], poly->f[j]);
  for (j = 0; j < 2; j ++)
    mpz_init_set (g[j], poly->g[j]);
  mpz_init_set (b, poly->g[1]);
  mpz_init_set (m, poly->g[0]);
  mpz_neg (m, m);
  mpz_init_set_si (V, 1);
  mpz_init_set_ui (tmpv, 0);

  /* look for positive V: 2, 4, 8, ... */
  for (i = 0; i < 150; i++, mpz_mul_si (V, V, 2) )
  {
    /* rotate by w*x */
    rotate_aux_mpz (poly->f, b, m, tmpv, V, 0);

    /* translation-optimize the rotated polynomial */
    for (j = 0; j <= poly->d; j ++)
      mpz_set (f[j], poly->f[j]);
    for (j = 0; j < 2; j ++)
      mpz_set (g[j], poly->g[j]);

    optimize_aux (f, poly->d, g, 0, 0, DEFAULT_L2_METHOD);

    skewness = L2_skewness ( f, poly->d, SKEWNESS_DEFAULT_PREC,
                             DEFAULT_L2_METHOD );

    lognorm = L2_lognorm (f, poly->d, skewness, DEFAULT_L2_METHOD);

    if (lognorm > bound->bound_lognorm)
      break;
  }

  mpz_set (bound->global_v_boundr, V);
  mpz_set_si (V, 0);

  /* rotate back */
  rotate_aux_mpz (poly->f, b, m, tmpv, V, 0);

  /* look for negative k: -2, -4, -8, ... */
  mpz_set_si (V, -1);
  for (i = 0; i < 150; i++, mpz_mul_si (V, V, 2))
  {
    /* rotate by w*x */
    rotate_aux_mpz (poly->f, b, m, tmpv, V, 0);

    /* translation-optimize the rotated polynomial */
    for (j = 0; j <= poly->d; j ++)
      mpz_set (f[j], poly->f[j]);
    for (j = 0; j < 2; j ++)
      mpz_set (g[j], poly->g[j]);

    optimize_aux (f, poly->d, g, 0, 0, DEFAULT_L2_METHOD);

    /* get lognorm */
    skewness = L2_skewness ( f, poly->d, SKEWNESS_DEFAULT_PREC,
                             DEFAULT_L2_METHOD );

    lognorm = L2_lognorm (f, poly->d, skewness, DEFAULT_L2_METHOD);

    if (lognorm > bound->bound_lognorm)
      break;
  }

  mpz_set (bound->global_v_boundl, V);

  /* rotate back */
  mpz_set_ui (V, 0);
  rotate_aux_mpz (poly->f, b, m, tmpv, V, 0);

  for (j = 0; j <= poly->d; j ++)
    mpz_clear (f[j]);
  for (j = 0; j < 2; j ++)
    mpz_clear (g[j]);
  free (f);
  free (g);
  mpz_clear (b);
  mpz_clear (m);
  mpz_clear (V);
  mpz_clear (tmpv);
}


/**
 * Find bound U for linear rotation.
 */
static inline void
rotate_bounds_U_lu ( ropt_poly_t poly,
                     ropt_bound_t bound )
{
  unsigned int i;
  int j;
  double skewness, lognorm;
  mpz_t *f, *g, b, m;
  f = (mpz_t *) malloc ( (poly->d + 1)* sizeof (mpz_t));
  g = (mpz_t *) malloc ( 2 * sizeof (mpz_t));
  for (j = 0; j <= poly->d; j ++)
    mpz_init_set (f[j], poly->f[j]);
  for (j = 0; j < 2; j ++)
    mpz_init_set (g[j], poly->g[j]);
  mpz_init_set (b, poly->g[1]);
  mpz_init_set (m, poly->g[0]);
  mpz_neg (m, m);

  /* Look for positive w: 1, 2, 4, 8, ... Note (sizeof (long) * 8 - 3)
     to prevent overflow in the rotate_aux(). */
  long w0 = 0, w = 1;
  for (i = 0; i < (sizeof (long) * 8 - 3); i++, w *= 2)
  {
    /* rotate by w*x */
    w0 = rotate_aux (poly->f, b, m, w0, w, 1);

    /* translation-optimize the rotated polynomial */
    for (j = 0; j <= poly->d; j ++)
      mpz_set (f[j], poly->f[j]);
    for (j = 0; j < 2; j ++)
      mpz_set (g[j], poly->g[j]);

    optimize_aux (f, poly->d, g, 0, 0, CIRCULAR);

    skewness = L2_skewness ( f, poly->d, SKEWNESS_DEFAULT_PREC,
                             DEFAULT_L2_METHOD );

    lognorm = L2_lognorm (f, poly->d, skewness, DEFAULT_L2_METHOD);
    /*
      fprintf (stderr, "# DEBUG --- [%d-th] U: %ld, lognorm: %f, "
      "bound_lognorm: %f\n", i, w, lognorm,  bound->bound_lognorm);
    */
    if (lognorm > bound->bound_lognorm)
      break;
  }

  bound->global_u_boundr = w;

  /* go back to w=0 */
  rotate_aux (poly->f, b, m, w0, 0, 1);

  /* Look for negative w: -1, -2, -4, -8, ... Note (sizeof (long) * 8 - 3)
     to prevent overflow in the rotate_aux(). */
  w0 = 0;
  w = -1;
  for (i = 0; i < (sizeof (long int) * 8 - 3); i++, w *= 2)
  {
    /* rotate by w*x */
    w0 = rotate_aux (poly->f, b, m, w0, w, 1);

    /* translation-optimize the rotated polynomial */
    for (j = 0; j <= poly->d; j ++)
      mpz_set (f[j], poly->f[j]);
    for (j = 0; j < 2; j ++)
      mpz_set (g[j], poly->g[j]);

    optimize_aux (f, poly->d, g, 0, 0, CIRCULAR);

    skewness = L2_skewness ( f, poly->d, SKEWNESS_DEFAULT_PREC,
                             DEFAULT_L2_METHOD );

    lognorm = L2_lognorm (f, poly->d, skewness, DEFAULT_L2_METHOD);
    /*
      fprintf (stderr, "# DEBUG --- [%d-th] U: %ld, lognorm: %f, "
      "bound_lognorm: %f\n", i, w, lognorm,  bound->bound_lognorm);
    */
    if (lognorm > bound->bound_lognorm)
      break;
  }

  bound->global_u_boundl = w;

  /* go back to w=0 */
  rotate_aux (poly->f, b, m, w0, 0, 1);

  for (j = 0; j <= poly->d; j ++)
    mpz_clear (f[j]);
  for (j = 0; j < 2; j ++)
    mpz_clear (g[j]);
  free (f);
  free (g);
  mpz_clear (b);
  mpz_clear (m);
}


/**
 * Find bound W for quadratic rotation.
 */
static inline void
rotate_bounds_W_lu ( ropt_poly_t poly,
                     ropt_bound_t bound )
{
  int i, j;
  double skewness, lognorm;
  mpz_t *f, *g, b, m;
  f = (mpz_t *) malloc ( (poly->d + 1)* sizeof (mpz_t));
  g = (mpz_t *) malloc ( 2 * sizeof (mpz_t));
  for (i = 0; i <= poly->d; i ++)
    mpz_init_set (f[i], poly->f[i]);
  for (i = 0; i < 2; i ++)
    mpz_init_set (g[i], poly->g[i]);
  mpz_init_set (b, poly->g[1]);
  mpz_init_set (m, poly->g[0]);
  mpz_neg (m, m);

  /* look for positive w: , ... 0, 1, 2 */
  long w0 = 0, w = 0;
  for (i = 0; i < 4096; i++, w++)
  {
    /* rotate by w*x */
    w0 = rotate_aux (poly->f, b, m, w0, w, 2);

    /* translation-optimize the rotated polynomial */
    for (j = 0; j <= poly->d; j ++)
      mpz_set (f[j], poly->f[j]);
    for (j = 0; j < 2; j ++)
      mpz_set (g[j], poly->g[j]);

    optimize_aux (f, poly->d, g, 0, 0, CIRCULAR);

    skewness = L2_skewness (f, poly->d, SKEWNESS_DEFAULT_PREC,
                            DEFAULT_L2_METHOD);

    lognorm = L2_lognorm (f, poly->d, skewness, DEFAULT_L2_METHOD);

    /*
      fprintf (stderr, "# DEBUG --- [%d-th] W: %ld, lognorm: %f, "
      "bound_lognorm: %f\n", i, w, lognorm,  bound->bound_lognorm);
    */       
    if (lognorm > bound->bound_lognorm)
      break;
  }

  bound->global_w_boundr = w;

  /* go back to w=0 */
  rotate_aux (poly->f, b, m, w0, 0, 2);

  /* look for negative w: , ... 0, -1, -2 */
  w0 = 0;
  w = 0;
  for (i = 0; i < 4096; i++, w--)
  {
    /* rotate by w*x */
    w0 = rotate_aux (poly->f, b, m, w0, w, 2);

    /* translation-optimize the rotated polynomial */
    for (j = 0; j <= poly->d; j ++)
      mpz_set (f[j], poly->f[j]);
    for (j = 0; j < 2; j ++)
      mpz_set (g[j], poly->g[j]);

    optimize_aux (f, poly->d, g, 0, 0, CIRCULAR);

    skewness = L2_skewness (f, poly->d, SKEWNESS_DEFAULT_PREC,
                            DEFAULT_L2_METHOD);

    lognorm = L2_lognorm (f, poly->d, skewness, DEFAULT_L2_METHOD);

    /*
      fprintf (stderr, "# DEBUG --- [%d-th] W: %ld, lognorm: %f,
      bound_lognorm: %f\n", i, w, lognorm,  bound->bound_lognorm);
    */
    if (lognorm > bound->bound_lognorm)
      break;
  }

  bound->global_w_boundl = w;
  
  /* go back to w=0 */
  rotate_aux (poly->f, b, m, w0, 0, 2);

  for (i = 0; i <= poly->d; i ++)
    mpz_clear (f[i]);
  for (i = 0; i < 2; i ++)
    mpz_clear (g[i]);
  free (f);
  free (g);
  mpz_clear (b);
  mpz_clear (m);
}


/* -----------------*/
/* Public Functions */
/* -----------------*/


/**
 * Init ropt_poly_t. 
 */
void
ropt_poly_init ( ropt_poly_t poly )
{
  unsigned int i;
  mpz_init (poly->n);
  mpz_init (poly->m);

  /* fx, gx holds pre-computed values f(r), g(r) where 0 <= r < p. */
  poly->f = (mpz_t*) malloc ((MAX_DEGREE + 1) * sizeof (mpz_t));
  poly->g = (mpz_t*) malloc ((MAX_DEGREE + 1) * sizeof (mpz_t));
  (poly->fx) = (mpz_t *) malloc ((primes[NP-1]+1) * sizeof (mpz_t));
  (poly->gx) = (mpz_t *) malloc ((primes[NP-1]+1) * sizeof (mpz_t));
  (poly->numerator) = (mpz_t *) malloc ((primes[NP-1]+1) * sizeof (mpz_t));

  if ( (poly->f == NULL) || (poly->g == NULL) ||
       ((poly->fx) == NULL) || ((poly->gx) == NULL) ||
       ((poly->numerator) == NULL) ) {
    fprintf (stderr, "Error, cannot allocate memory for polynomial"
             " coefficients in ropt_poly_init().\n");
    exit(1);
  }

  for (i = 0; i <= MAX_DEGREE; i++) {
    mpz_init (poly->f[i]);
    mpz_init (poly->g[i]);
  }

  for (i = 0; i <= primes[NP-1]; i++) {
    mpz_init (poly->fx[i]);
    mpz_init (poly->gx[i]);
    mpz_init (poly->numerator[i]);
  }
}


/**
 * Evaluation polynomials at many points. 
 */
static inline void
ropt_poly_setup_eval ( mpz_t *f,
                       mpz_t *g,
                       mpz_t *fr,
                       mpz_t *gr,
                       mpz_t *numerator,
                       const unsigned int *p,
                       int d )
{
  unsigned long i;
  mpz_t tmp;
  mpz_init (tmp);

  for ( i = 0; i <= p[NP-1]; i ++ ) {
    mpz_set_ui (fr[i], 0);
    mpz_set_ui (gr[i], 0);
    eval_poly_ui (fr[i], f, d, i);
    eval_poly_ui (gr[i], g, 1, i);
    eval_poly_diff_ui (numerator[i], f, d, i);
    mpz_mul (numerator[i], numerator[i], gr[i]);
    mpz_neg (numerator[i], numerator[i]);
    mpz_mul (tmp, fr[i], g[1]);
    mpz_add (numerator[i], tmp, numerator[i]);
  }

  mpz_clear (tmp);
}


/**
 * Precompute fx, gx and numerator in ropt_poly_t. Note: poly->f, 
 * poly->g, poly->d, poly->n must be set in prior.
 * This function can be called to reset poly after rotation.
 */
void
ropt_poly_setup ( ropt_poly_t poly )
{
  int i;
  mpz_t t;
  mpz_init (t);

  /* degree */
  for ( (poly->d) = MAX_DEGREE; 
        (poly->d) > 0 && mpz_cmp_ui ((poly->f[poly->d]), 0) == 0;
        poly->d -- );

  /* m = -Y0/Y1 mod n */
  mpz_invert (poly->m, poly->g[1], poly->n);
  mpz_neg (poly->m, poly->m);
  mpz_mul (poly->m, poly->m, poly->g[0]);
  mpz_mod (poly->m, poly->m, poly->n);

  /* check if m is a root of f mod n */
  mpz_set (t, poly->f[poly->d]);

  for (i = poly->d - 1; i >= 0; i --) {
    mpz_mul (t, t, poly->m);
    mpz_add (t, t, poly->f[i]);
    mpz_mod (t, t, poly->n);
  }

  if (mpz_cmp_ui (t, 0) != 0) {
    fprintf (stderr, "ERROR: The following polynomial have no common"
             " root. \n");
    print_cadopoly_fg (stderr, poly->f, poly->g, poly->d, poly->n);
    exit (1);
  }

  /* pre-compute f(r) for all r < B */
  ropt_poly_setup_eval ( poly->f, poly->g, poly->fx, poly->gx, 
                         poly->numerator, primes, poly->d );

  /* projective alpha */
  poly->alpha_proj = get_biased_alpha_projective (poly->f, poly->d, 2000);

  mpz_clear (t);
}


/**
 * Free ropt_poly_t.
 */
void
ropt_poly_free ( ropt_poly_t poly )
{
  unsigned int i;

  for (i = 0; i <= primes[NP-1]; i ++) {
    mpz_clear(poly->fx[i]);
    mpz_clear(poly->gx[i]);
    mpz_clear(poly->numerator[i]);
  }

  for (i = 0; i <= MAX_DEGREE; i++) {
    mpz_clear (poly->f[i]);
    mpz_clear (poly->g[i]);
  }

  mpz_clear (poly->n);
  mpz_clear (poly->m);

  free (poly->fx);
  free (poly->gx);
  free (poly->numerator);
  free (poly->f);
  free (poly->g);
}


/**
 * Init ropt_bound_t. 
 */
void
ropt_bound_init ( ropt_bound_t bound )
{
  bound->global_w_boundl = 0;
  bound->global_w_boundr = 0;
  bound->global_u_boundl = 0;
  bound->global_u_boundr = 0;
  mpz_init_set_ui (bound->global_v_boundl, 0UL);
  mpz_init_set_ui (bound->global_v_boundr, 0UL);
  bound->init_lognorm = 0.0;
  bound->bound_lognorm = 0.0;
  bound->bound_lognorm_ratio = 0.0;
  bound->exp_min_alpha = 0.0;
}


/**
 * Subroutine for ropt_bound_setup().
 */
static inline void
ropt_bound_setup_normbound ( ropt_poly_t poly,
                             ropt_bound_t bound,
                             ropt_param_t param )
{
  double skewness = L2_skewness ( poly->f, poly->d, SKEWNESS_DEFAULT_PREC,
                                  DEFAULT_L2_METHOD );
  bound->init_lognorm = L2_lognorm ( poly->f, poly->d, skewness,
                                     DEFAULT_L2_METHOD );

  /* setup lognorm bound, either from input or by default. */
  if (param->bound_lognorm > 0) {
    bound->bound_lognorm = param->bound_lognorm;
  }
  else {
    /* The higher, the more margin in computing the sieving bound
       w, u and v, hence the larger the sieving bound, and hence
       larger individual sublattices. */
    bound->bound_lognorm_ratio = BOUND_LOGNORM_RATIO;
    bound->bound_lognorm = bound->init_lognorm * bound->bound_lognorm_ratio;
  }
}


/**
 * Subroutine for ropt_bound_setup().
 */
static inline void
ropt_bound_setup_globalbound ( ropt_poly_t poly,
                               ropt_bound_t bound,
                               ropt_param_t param )
{
  /* w bound */
  if (poly->d == 6) {
    if (param->w_length > 0) {
      bound->global_w_boundl = param->w_left_bound;
      bound->global_w_boundr = param->w_left_bound + param->w_length - 1;
    }
    else
      rotate_bounds_W_lu (poly, bound);
  }
  else {
    bound->global_w_boundr = 0;
    bound->global_w_boundl = 0;
    param->w_left_bound = 0;
    param->w_length = 1;
  }

  /* u bound */
  rotate_bounds_U_lu (poly, bound);

  /* v bound */
  rotate_bounds_V_mpz (poly, bound);

}


/**
 * Subroutine for ropt_bound_setup().
 */
static inline void
ropt_bound_setup_others ( ropt_bound_t bound )
{
  long size, sizel, sizer;
  sizel = mpz_sizeinbase (bound->global_v_boundl, 2);
  sizer = mpz_sizeinbase (bound->global_v_boundr, 2);
  size = (sizel > sizer) ? sizel : sizer;

  if (bound->global_u_boundr == 0)
    bound->global_u_boundr = 1;
  if (bound->global_u_boundl == 0)
    bound->global_u_boundr = 1;

  /*
    printf ("size: %ld, uboundr: %ld, uboundl: %ld\n", 
            size, bound->global_u_boundr, bound->global_u_boundl);
  */

  size += (int) log2 ((double) (bound->global_u_boundr -
                                bound->global_u_boundl));
  if (size > 150)
    size = 150;
  
  bound->exp_min_alpha = exp_alpha[size-1];
}


/**
 * Setup ropt_bound_t (independent of manually-input param).
 * Note, this function should be called in the very beginning
 * before doing any rotation, since the init_lognorm parameter
 * will be set to decide the rotation range in the rest.
 */
void
ropt_bound_setup ( ropt_poly_t poly,
                   ropt_bound_t bound,
                   ropt_param_t param )
{
  /* set bound->bound_lognorm */
  ropt_bound_setup_normbound (poly, bound, param);

  /* set w, u, v bounds */
  ropt_bound_setup_globalbound (poly, bound, param);

  /* set exp_min_alpha */
  ropt_bound_setup_others (bound);

  if (param->verbose >= 1) {
    gmp_fprintf ( stderr, "# Info: global bounds (%d:%d, %ld:%ld, %Zd:%Zd)"
                  " gives:\n",
                  bound->global_w_boundl,
                  bound->global_w_boundr,
                  bound->global_u_boundl,
                  bound->global_u_boundr,
                  bound->global_v_boundl,
                  bound->global_v_boundr );
    gmp_fprintf ( stderr, "# Info: exp_alpha: %.3f, norm bound: %.3f, "
                  "initial norm: %.3f\n",
                  bound->exp_min_alpha,
                  bound->bound_lognorm,
                  bound->init_lognorm );
  }
}


/**
 * For existing rsparam and when rs is changed,
 */
void
ropt_bound_reset ( ropt_poly_t poly,
                   ropt_bound_t bound,
                   ropt_param_t param )
{

  /* u bound */
  rotate_bounds_U_lu (poly, bound);

  /* v bound */
  rotate_bounds_V_mpz (poly, bound);

  /* set exp_min_alpha */
  ropt_bound_setup_others (bound);

  if (param->verbose >= 2) {
    gmp_fprintf ( stderr, "# Info: reset bounds (%d:%d, %ld:%ld, %Zd:%Zd)"
                  " gives:\n",
                  bound->global_w_boundl,
                  bound->global_w_boundr,
                  bound->global_u_boundl,
                  bound->global_u_boundr,
                  bound->global_v_boundl,
                  bound->global_v_boundr );
    gmp_fprintf ( stderr, "# Info: exp_alpha: %.3f, L2 bound: %.3f, "
                  "initial L2: %.3f\n",
                  bound->exp_min_alpha,
                  bound->bound_lognorm,
                  bound->init_lognorm );
  }
}


/**
 * Free ropt_bound_t. 
 */
void
ropt_bound_free ( ropt_bound_t bound )
{
  mpz_clear (bound->global_v_boundl);
  mpz_clear (bound->global_v_boundr);
}


/**
 * Init stage 1 parameters. The default customisation/parameters
 * happens in ropt_s1param_setup().
 */
void
ropt_s1param_init ( ropt_s1param_t s1param )
{
  int i;

  /* will be set either from param (stdin) or by default */
  s1param->len_e_sl = 0;
  s1param->tlen_e_sl = 0;
  s1param->nbest_sl = 0;

  /* set to 1 for using quicker, smaller nbest_sl for tunning
     sublattices */
  s1param->nbest_sl_tunemode = 0;
  
  mpz_init_set_ui (s1param->modulus, 1UL);

  s1param->e_sl = (unsigned int*) 
    malloc ( NUM_SUBLATTICE_PRIMES * sizeof (unsigned int) );

  s1param->individual_nbest_sl = (unsigned int*) 
    malloc ( NUM_SUBLATTICE_PRIMES * sizeof (unsigned int) );

  if (s1param->e_sl == NULL || s1param->individual_nbest_sl == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in "
             "ropt_s1param_init().\n");
    exit (1);
  }

  for (i = 0; i < NUM_SUBLATTICE_PRIMES; i ++) {
    s1param->e_sl[i] = 0;
    s1param->individual_nbest_sl[i] = 0;
  }
}


/**
 * Helper function to setup "e_sl[]".
 */
static inline void
ropt_s1param_setup_e_sl ( ropt_poly_t poly,
                          ropt_s1param_t s1param,
                          ropt_bound_t bound,
                          ropt_param_t param )
{
  unsigned int i, j;
  unsigned long sublattice;
  mpz_t bound_by_v;
  mpz_init (bound_by_v);

  /* e_sl */
  unsigned long bound_by_u = (unsigned long)
    (bound->global_u_boundr < (-bound->global_u_boundl)) ?
    bound->global_u_boundr : (-bound->global_u_boundl);

  mpz_fdiv_q_ui (bound_by_v, bound->global_v_boundr,
                 SIZE_SIEVEARRAY_V_MIN);

  /* fix for small skewnss */
  if ( mpz_cmp_ui (bound_by_v, bound_by_u) < 0 ) {
    sublattice = mpz_get_ui (bound_by_v);
    if (sublattice > bound_by_u / 16)
      sublattice = bound_by_u / 16;
  }
  else 
    sublattice = bound_by_u / 16;

  double skew = L2_skewness (poly->f, poly->d, SKEWNESS_DEFAULT_PREC, 
                             DEFAULT_L2_METHOD);

  /* fix for small skewness but large bound */
  if ( (double) sublattice > (skew / 16.0) )
    sublattice = (unsigned long) (skew / 16.0);

  for (i = 0; i < NUM_DEFAULT_SUBLATTICE; i ++)
    if (default_sublattice_prod[i] > sublattice)
      break;
  /* fix if i too large */
  i = (i >= NUM_DEFAULT_SUBLATTICE) ? 
    (NUM_DEFAULT_SUBLATTICE - 1) : i;

  /* set e_sl[] from default array */
  for (j = 0; j < NUM_SUBLATTICE_PRIMES; j++) {
    s1param->e_sl[j] = default_sublattice_pe[i][j];
  }

  /* overwrite e_sl[] from from stdin, if needed */
  if (param->s1_num_e_sl != 0)
    for (i = 0; i < NUM_SUBLATTICE_PRIMES; i++)
      s1param->e_sl[i] = param->s1_e_sl[i];

  mpz_clear (bound_by_v);
}


/**
 *  Function to setup "individual_nbest_sl[]".
 */
void
ropt_s1param_setup_individual_nbest_sl (ropt_s1param_t s1param)
{
  unsigned int i;
  for (i = 0; i < s1param->len_e_sl; i ++)
    {
      s1param->individual_nbest_sl[i] =
        size_each_sublattice[s1param->tlen_e_sl - 1][i];
#ifdef HAVE_MINGW
      fprintf (stderr, "# ropt_str:792 individual_nbest_sl[%u] = %d\n", i,
               s1param->individual_nbest_sl[i]);
#endif
    }
}


/**
 *  Function to setup shorter "individual_nbest_sl[]".
 */
void
ropt_s1param_setup_individual_nbest_sl_tune (ropt_s1param_t s1param)
{
  unsigned int i;
  for (i = 0; i < s1param->len_e_sl; i ++)
    {
      s1param->individual_nbest_sl[i] = size_each_sublattice_tune[i];
#ifdef HAVE_MINGW
      fprintf (stderr, "# ropt_str:810: individual_nbest_sl[%u] = %d\n", i,
               s1param->individual_nbest_sl[i]);
#endif
    }
}


/**
 * Setup s1param by default parameters and/or param (stdin).
 */
void
ropt_s1param_setup ( ropt_poly_t poly,
                     ropt_s1param_t s1param,
                     ropt_bound_t bound,
                     ropt_param_t param )
{
  unsigned int i, j;
  
  /* Set 1: "len_e_sl" and "tlen_e_sl" */
  s1param->len_e_sl = NUM_SUBLATTICE_PRIMES;
  s1param->tlen_e_sl = s1param->len_e_sl;

  /* Set 2: "nbest_sl", depending on the size of n */
  j = mpz_sizeinbase (poly->n, 10);
  for (i = 0; i < 7; i ++)
    if (size_total_sublattices[i][0] > j)
      break;
  
  s1param->nbest_sl = size_total_sublattices[i][1];

  //printf ("s1param->nbest_sl: %u\n", s1param->nbest_sl);
  
  /* Set 3: set "e_sl[]" */
  ropt_s1param_setup_e_sl (poly, s1param, bound, param);

  /* Set 4: set "individual_nbest_sl[]" */
  ropt_s1param_setup_individual_nbest_sl (s1param);
}


/**
 * Free s1param.
 */
void
ropt_s1param_free ( ropt_s1param_t s1param )
{
  free(s1param->e_sl);
  free(s1param->individual_nbest_sl);
  mpz_clear (s1param->modulus);
}


/**
 * Init ropt_s2param.
 */
void
ropt_s2param_init ( ropt_poly_t poly,
                    ropt_s2param_t s2param )
{
  int i;
  s2param->len_p_rs = NP - 1;
  s2param->Amax = s2param->Amin = 0;
  s2param->Bmax = s2param->Bmin = 0;

  mpz_init_set_ui (s2param->Umax, 0UL);
  mpz_init_set_ui (s2param->Umin, 0UL);
  mpz_init_set_ui (s2param->Vmax, 0UL);
  mpz_init_set_ui (s2param->Vmin, 0UL);
  mpz_init_set_ui (s2param->A, 0UL);
  mpz_init_set_ui (s2param->B, 0UL);
  mpz_init_set_ui (s2param->MOD, 0UL);

  s2param->f = (mpz_t*) malloc ((poly->d + 1) * sizeof (mpz_t));
  s2param->g = (mpz_t*) malloc (2 * sizeof (mpz_t));
  if (s2param->f == NULL || s2param->g == NULL) {
    fprintf ( stderr, "Error, cannot allocate memory in "
              "ropt_s2param_init().\n" );
    exit (1);
  }

  for (i = 0; i <= poly->d; i++)
    mpz_init (s2param->f[i]);
  mpz_init (s2param->g[0]);
  mpz_init (s2param->g[1]);
}


/**
 * Free ropt_s2param.
 */
void
ropt_s2param_free ( ropt_poly_t poly,
                    ropt_s2param_t s2param )
{
  int i;
  
  mpz_clear (s2param->Umax);
  mpz_clear (s2param->Umin);
  mpz_clear (s2param->Vmax);
  mpz_clear (s2param->Vmin);
  mpz_clear (s2param->A);
  mpz_clear (s2param->B);
  mpz_clear (s2param->MOD);
  for (i = 0; i <= poly->d; i++)
    mpz_clear (s2param->f[i]);
  mpz_clear (s2param->g[0]);
  mpz_clear (s2param->g[1]);
  free (s2param->f);
  free (s2param->g);
}


/**
 * Setup sublattice (u, v), mod and Umax, Umin, Vmax, Vmin.
 */
static inline void
ropt_s2param_setup_sublattice ( ropt_s2param_t s2param,
                                mpz_t A,
                                mpz_t B,
                                mpz_t MOD )
{
  mpz_set (s2param->A, A);
  mpz_set (s2param->B, B);
  /* the s1param->modulus might be not true since rsparam might be
     changed for different quadratic rotations. Instead, the true mod
     is recorded in the priority queue, now in 'MOD'. */
  mpz_set (s2param->MOD, MOD);
  ab2uv (s2param->A, s2param->MOD, s2param->Amax, s2param->Umax);
  ab2uv (s2param->A, s2param->MOD, s2param->Amin, s2param->Umin);
  ab2uv (s2param->B, s2param->MOD, s2param->Bmax, s2param->Vmax);
  ab2uv (s2param->B, s2param->MOD, s2param->Bmin, s2param->Vmin);
}


/**
 * Set sieving region for s2param -> Amax, Amin, Amax, Amin.
 */
static inline void
ropt_s2param_setup_range ( ropt_bound_t bound,
                           ropt_s2param_t s2param,
                           ropt_param_t param,
                           mpz_t mod )
{
  /* read sieve length from param (stdin) */
  if (param->s2_Amax >= 0 && param->s2_Bmax > 0) {
    s2param->Amax = param->s2_Amax;
    s2param->Bmax = param->s2_Bmax;
  }
  /* or compute sieve V length */
  else {
    s2param->Amax = 0;
    unsigned long len;
    mpz_t q;
    mpz_init (q);
    mpz_fdiv_q (q, bound->global_v_boundr, mod);
    len =  mpz_get_ui (q);

    /* upper bound */
    s2param->Bmax = ( (len > SIZE_SIEVEARRAY_V_MAX) ?
                      SIZE_SIEVEARRAY_V_MAX : (long) len );

    /* fix if len is too small */
    if (s2param->Bmax == 0)
      s2param->Bmax = 128;

    mpz_clear (q);
  }

  s2param->Bmin = -s2param->Bmax;
  s2param->Amin = -s2param->Amax;
}


/**
 * Set up sieving length and sublattice for s2param.
 */
void
ropt_s2param_setup ( ropt_bound_t bound,
                     ropt_s1param_t s1param,
                     ropt_s2param_t s2param,
                     ropt_param_t param,
                     mpz_t A,
                     mpz_t B,
                     mpz_t MOD )
{
  /* normally we should have more primes than those in the sublattice */
  if (s2param->len_p_rs < s1param->len_e_sl) {
    fprintf ( stderr, "# Warning: number of primes considered "
              "in stage 2 is smaller than that in "
              "stage 1. See ropt_s2param_setup().\n" );
  }

  /* set sieving length */
  ropt_s2param_setup_range (bound, s2param, param, MOD);

  /* set sublattice */
  ropt_s2param_setup_sublattice (s2param, A, B, MOD);
}


/**
 * Set up s2param without s1param.
 */
void
ropt_s2param_setup_stage2_only ( ropt_bound_t bound,
                                 ropt_s2param_t s2param,
                                 ropt_param_t param,
                                 mpz_t A,
                                 mpz_t B,
                                 mpz_t MOD )
{
  s2param->len_p_rs = NP - 1;

  /* set sieving length */
  ropt_s2param_setup_range (bound, s2param, param, MOD);

  /* set sublattice */
  ropt_s2param_setup_sublattice (s2param, A, B, MOD);
}


/**
 * Set tune (short) sieving region for s2param -> Amax, Amin, Amax, Amin.
 */
static inline void
ropt_s2param_setup_tune_range  ( ropt_s2param_t s2param,
                                 unsigned long Amax,
                                 unsigned long Bmax )
{
  s2param->Amax = (long) Amax;
  s2param->Bmax = (long) Bmax;
  s2param->Bmin = -s2param->Bmax;
  s2param->Amin = -s2param->Amax;
}


/**
 * Set up tune sieving length and sublattice for s2param.
 */
void
ropt_s2param_setup_tune ( ropt_s1param_t s1param,
                          ropt_s2param_t s2param,
                          mpz_t A,
                          mpz_t B,
                          mpz_t MOD,
                          unsigned long Amax,
                          unsigned long Bmax,
                          unsigned int len_p_rs )
{
  /* setup s2param->len_p_rs */
  if (len_p_rs < s1param->tlen_e_sl) {
    fprintf ( stderr, "# Warning: number of primes considered "
              "in stage 2 (%d) is smaller than that (%d) in "
              "stage 1. See ropt_s2param_setup().\n",
              len_p_rs, s1param->tlen_e_sl );
  }
  s2param->len_p_rs = len_p_rs;
  
  /* set tune sieving length */
  ropt_s2param_setup_tune_range (s2param, Amax, Bmax);

  /* set sublattice */
  ropt_s2param_setup_sublattice (s2param, A, B, MOD);
}


/**
 * Print s2param.
 */
void
ropt_s2param_print ( ropt_s2param_t s2param )
{
  fprintf ( stderr, "# Info: sieving matrix: "
            "[%ld, %ld] x [%ld, %ld], len. prime = %d.\n",
            s2param->Amin, s2param->Amax,
            s2param->Bmin, s2param->Bmax,
            s2param->len_p_rs );

  gmp_fprintf ( stderr,
                "# Info: (u, v) = (%Zd + i * %Zd, %Zd + j * %Zd)\n",
                s2param->A, s2param->MOD, s2param->B,
                s2param->MOD );

  gmp_fprintf ( stderr,
                "# Info: (Amin: %4ld, Amax: %4ld) -> "
                "(Umin: %6Zd, Umax: %6Zd)\n",
                s2param->Amin, s2param->Amax,
                s2param->Umin, s2param->Umax );

  gmp_fprintf ( stderr,
                "# Info: (Bmin: %4ld, Bmax: %4ld) -> "
                "(Vmin: %6Zd, Vmax: %6Zd)\n",
                s2param->Bmin, s2param->Bmax,
                s2param->Vmin, s2param->Vmax );
}


/**
 * Init bestpoly.
 */
void
ropt_bestpoly_init ( ropt_bestpoly_t bestpoly,
                     int d )
{
  int i;
  bestpoly->f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  bestpoly->g = (mpz_t*) malloc (2 * sizeof (mpz_t));

  if ((bestpoly->f == NULL) || (bestpoly->g == NULL)) {
    fprintf ( stderr, "Error, cannot allocate memory for"
              " ropt_bestpoly_init()\n" );
    exit (1);
  }

  for (i = 0; i <= d; i++)
    mpz_init (bestpoly->f[i]);
  mpz_init (bestpoly->g[0]);
  mpz_init (bestpoly->g[1]);
}


/**
 * Setup bestpoly.
 */
void
ropt_bestpoly_setup ( ropt_bestpoly_t bestpoly,
                      mpz_t *f,
                      mpz_t *g,
                      int d )
{
  int i;
  for (i = 0; i <= d; i++)
    mpz_set (bestpoly->f[i], f[i]);
  mpz_set (bestpoly->g[0], g[0]);
  mpz_set (bestpoly->g[1], g[1]);
}


/**
 * Free bestpoly.
 */
void
ropt_bestpoly_free ( ropt_bestpoly_t bestpoly,
                     int d )
{
  int i;
  for (i = 0; i <= d; i++)
    mpz_clear (bestpoly->f[i]);

  mpz_clear (bestpoly->g[0]);
  mpz_clear (bestpoly->g[1]);
  free (bestpoly->f);
  free (bestpoly->g);
}


/**
 * Init param.
 */
void
ropt_param_init ( ropt_param_t param )
{
  mpz_init (param->s2_u);
  mpz_init (param->s2_v);
  mpz_init (param->s2_mod);
  mpz_init (param->n);
  mpz_set_ui (param->s2_u, 0);
  mpz_set_ui (param->s2_v, 0);
  mpz_set_ui (param->s2_mod, 0);
  mpz_set_ui (param->n, 0);
  param->w_left_bound = 0;
  param->w_length = 0;
  param->s1_num_e_sl = 0;
  param->s2_Amax = 0;
  param->s2_Bmax = 0;
  param->s2_w = 0;
  param->bound_lognorm = 0;
  param->s1_e_sl = (unsigned short*)
    malloc ( NUM_SUBLATTICE_PRIMES * sizeof (unsigned short) );
  int i;
  for (i = 0; i < NUM_SUBLATTICE_PRIMES; i ++)
    param->s1_e_sl[i] = 0;
  param->d = 0;
  param->verbose = 0;
}


/**
 * Free param.
 */
void
ropt_param_free ( ropt_param_t param )
{
  mpz_clear (param->s2_u);
  mpz_clear (param->s2_v);
  mpz_clear (param->s2_mod);
  mpz_clear (param->n);
  free (param->s1_e_sl);
}


/**
 * Init info.
 */
void
ropt_info_init ( ropt_info_t info )
{
  info->ave_MurphyE = 0.0;
  info->best_MurphyE = 0.0;
  info->mode = 0;
  info->w = 0;
}


/**
 * Free info.
 */
void
ropt_info_free ( ropt_info_t info )
{
  info->ave_MurphyE = 0.0;
  info->best_MurphyE = 0.0;
  info->mode = 0;
  info->w = 0;
}
