/**
 * @file ropt_stage2.c
 * Called by ropt.c to root sieve.
 */

#include "ropt_stage2.h"

/*
  Change coordinate from (i, j) to (u, v).
*/
#if DEBUG
void
print_sievearray ( double **A,
                   long A0,
                   long A1,
                   long B0,
                   long B1,
                   unsigned long K_ST,
                   unsigned long J_ST,
                   unsigned long MOD )
{
  long i, j;

  for (i = 0; i < B1 - B0 + 1; i ++) {
    for (j = 0; j < A1 - A0 + 1; j ++) {
      if (j % 16 == 0) {
        printf ("\n");
        printf ("[%3ld*%lu+%lu, %4ld*%lu+%lu]=[%4ld, %4ld]: ",
                i + B0, MOD, K_ST, j + A0, MOD, J_ST, (i + B0)*MOD+K_ST, (j + A0)*MOD+J_ST);
      }
      if (j % 4 == 0)
        printf (" ");
      printf ("%1.1f ", A[i][j]);
    }
    printf ("\n");
  }
  printf ("\n");
}
#endif

/*
  init the best poly (for the return)
*/
static inline void
sievearray_init ( sievearray_t sievearray,
                  unsigned int i,
                  unsigned int j )
{
  /* Init array A
     sage: sum ([p/(p^2-1)*log(p) for p in prime_range(200)])
     4.2774597204802731 */
  float tmpf = SUP_ALPHA;
  int16_t tmpu = (int16_t) ceil (tmpf * 1000.0);
  unsigned long len = i * j, k;

  /* allocate matrix A. */
  sievearray->len_i = i;
  sievearray->len_j = j;
  sievearray->array = (int16_t *) malloc ( len * sizeof (int16_t));

  if ( sievearray->array != NULL )
    for (k = 0; k < len; k++)
      sievearray->array[k] = tmpu;
  else {
    fprintf (stderr, "Error, cannot allocate memory in sievearray_init().\n");
    exit (1);
  }
}


/*
  re-set root sieve array with biased alpha.
*/
static inline void
sievearray_reset ( sievearray_t sievearray )
{
  float tmpf = SUP_ALPHA;
  int16_t tmpu = (int16_t) ceil (tmpf * 1000.0);
  unsigned long k, len = sievearray->len_i * sievearray->len_j;

  if (sievearray->array != NULL)
    for (k = 0; k < len; k++)
      sievearray->array[k] = tmpu;
  else {
    fprintf (stderr, "Error, null memory in sievearray_reset().\n");
    exit (1);
  }
}


/*
  re-set root sieve array with biased alpha.
*/
static inline void
sievearray_free ( sievearray_t sievearray )
{
  free (sievearray->array);
}


/*
  root sieve over ARRAY with initial point
  (i, j) (mod pe), where i is fixed, so only
  need to sieve along the line j.
*/
static inline long
rootsieve_run_line ( int16_t *ARRAY,
                     long V,
                     long j,
                     unsigned int pe,
                     int16_t sub )
{
  long pel = (long) pe;

  /* Bounding j with B0 < tmp + p*j < B1 */
  while ( j <= V ) {
    ARRAY[j] = ARRAY[j] - sub;
    j += pel;
  }
  return j;
}


/*
  Compute the indices for some special cases
  Note that pl < pe
*/
static inline long
rootsieve_run_multroot_lift_idx ( unsigned int v,
                                  long Bmin,
                                  mpz_t B,
                                  mpz_t MOD,
                                  unsigned int pl,
                                  unsigned int pe )
{
  if (pl >= pe) {
    fprintf (stderr, "Non-invertible element in rootsieve_run_multroot_lift_idx().\n");
    exit(1);
  }
  /* findthe val_p (MOD), which must be < e */
  unsigned int a, b, c;
  long j_idx;

  /* (B-v)/p^l + (MOD/p^l)*k = 0 (mod p^{e-l})
     Note (B-v)/p^l is exact due to our constructions. */
  a = (unsigned int) mpz_fdiv_ui (B, pe);
  a = (a + pe - v) / pl; // +pe ensure it is positive
  b = (unsigned int) mpz_fdiv_ui (MOD, pe);
  b = b /pl;
  pl = pe / pl;  // it is now p^{e-l}
  c = (unsigned int) solve_lineq (a, b, 0, pl);
  j_idx = (long) ceil (((double) Bmin - (double) c) / (double) pl)
    * (long) pl + (long) c;
  j_idx = ab2ij (Bmin, j_idx);
  return j_idx;
}

#define DEBUG_MULTROOT_LIFT 0
/*
  rootsieve_run_multroot_lift for the multiple root. Since we are
  sure that level 1 node only contains one multiple root
  r, we don't need to (and can't, otherwise, may count
  repeatedly) consider single root at all.
*/
static inline void
rootsieve_run_multroot_lift ( node *currnode,
                              sievearray_t sa,
                              unsigned int *f_ui,
                              unsigned int *g_ui,
                              unsigned int *fuv_ui,
                              int d,
                              rsbound_t rsbound,
                              unsigned int p,
                              unsigned int max_e,
                              unsigned int curr_e,
                              int16_t sub )
{
  /* recursion end */
  if (currnode == NULL || curr_e > max_e || sub == 0)
    return;

  /* variables */
  int16_t subtmp;
  unsigned int i, j, k, l, nroots, pe, pl, pem1, fr, gr, step;
  node *tmpnode = NULL, *tmpnode2 = NULL;
  long i_idx, j_idx;

  /* compute p^e */
  pe = pem1 = 1;
  for (k = 0; k < curr_e - 1; k ++)
    pem1 = pem1 * p;
  pe = pem1 * p;

  /* val of p in MOD */
  l = 0;
  pl = p;
  while (mpz_divisible_ui_p(rsbound->MOD, pl) != 0) {
    l ++;
    pl *= p;
  } // pl is p^{l+1}
  pl /= p;

  /* loop until all siblings are checked. */
  while ( currnode != NULL ) {

    /* loop all (lifted multiple) roots */
    for (nroots = 0; nroots < currnode->nr; nroots++) {

      /* compute g(r) */
      gr = eval_poly_ui_mod (g_ui, 1, currnode->r[nroots], pe);
      if (gr % p == 0)
        continue;

      /* compute f_uv(x) and then evaluate it at r. */
      compute_fuv_ui (fuv_ui, f_ui, g_ui, d, currnode->u, currnode->v, pe);
      fr = eval_poly_ui_mod (fuv_ui, d, currnode->r[nroots], pe);
      fr = fr / pem1;

      /* solve on fr + gr*x = 0 (mod p), where x = uu*r + vv. */
      fr = (unsigned int) solve_lineq (fr, gr, 0, p); //fr = fr % p;

      /* solve (i, j) in  fr = i*r + j (mod p).
         - if r is not invertible, loop i with fixed j;
         - otherwise, loop j to solve i. */
      if (currnode->r[nroots] % p == 0) {

        for (i = 0; i < p; i ++) {

          /* if p|MOD and current pe|MOD:
             B + MOD*k = v (mod pe), so B = v (mod pe)
             A + MOD*k = u (mod pe), so A = u (mod pe)
             otherwise, no need to record this (u, v) at all. */
          if (mpz_divisible_ui_p (rsbound->MOD, p) != 0) {
            if (mpz_divisible_ui_p (rsbound->MOD, pe) != 0) {
              if ( mpz_fdiv_ui (rsbound->B, pe) != (currnode->v + fr * pem1) ||
                   mpz_fdiv_ui (rsbound->A, pe) != (currnode->u + i * pem1) )
              {
#if DEBUG_MULTROOT_LIFT
                fprintf (stderr, "skip found (%u, %u) in level %u.\n",
                         currnode->u + i * pem1, currnode->v + fr * pem1, curr_e);
#endif
                continue;
              }
            }
          }

          /* r is a multiple root, add r + k * p^{e-1}. */
          for (k = 0; k < p; k ++) {

#if DEBUG_MULTROOT_LIFT
            fprintf ( stderr, "level %u, (%u, %u), r: %u\n",
                      curr_e, currnode->u + i * pem1, currnode->v + fr * pem1,
                      currnode->r[nroots] + k * pem1 );
#endif

            insert_node ( currnode, &tmpnode,
                          currnode->u + i * pem1,
                          currnode->v + fr * pem1,
                          currnode->r[nroots] + k * pem1,
                          curr_e, p, pe, 0 );
          }

          /* r is a multiple root, add r + k * p^{e-1}. */
          for (k = 0; k < p; k ++) {
            insert_node (currnode, &tmpnode, currnode->u + pem1 * i,
                         currnode->v + pem1 * fr, currnode->r[nroots] + k * pem1, curr_e, p, pe, 2);
          }
        }
      }
      /* p -|- (currnode->r[nroots]), loop j to solve i. */
      else {
        for (j = 0; j < p; j ++) {

          /* given j, solve i in  fr = i*r + j (mod p). */
          i = solve_lineq (j, currnode->r[nroots], fr, p);

          /* if p|MOD and current pe|MOD:
             B + MOD*k = v (mod pe), so B = v (mod pe)
             A + MOD*k = u (mod pe), so A = u (mod pe)
             otherwise, no need to record this (u, v) at all. */
          if (mpz_divisible_ui_p (rsbound->MOD, p) != 0) {
            if (mpz_divisible_ui_p (rsbound->MOD, pe) != 0) {
              if ( mpz_fdiv_ui (rsbound->B, pe) != (currnode->v + j * pem1) ||
                   mpz_fdiv_ui (rsbound->A, pe) != (currnode->u + i * pem1) )
              {
#if DEBUG_MULTROOT_LIFT
                fprintf (stderr, "skip found (%u, %u) in level %u.\n",
                         currnode->u + i * pem1, currnode->v + fr * pem1, curr_e);
#endif
                continue;
              }
            }
          }

          /* r is a multiple root, add r + k * p^{e-1}. */
          for (k = 0; k < p; k ++) {
#if DEBUG_MULTROOT_LIFT
            fprintf ( stderr, "level %u, (%u, %u), r: %u\n",
                      curr_e, currnode->u + i * pem1, currnode->v + j * pem1,
                      currnode->r[nroots] + k * pem1 );
#endif
            insert_node (currnode, &tmpnode, currnode->u + pem1 * i,
                         currnode->v + pem1 * j, currnode->r[nroots] + k * pem1, curr_e, p, pe, 2);
          }
        }
      }
    }  // next root of current (u, v)

    /* recursieve to next level, curr_e + 1 */
    rootsieve_run_multroot_lift ( currnode->firstchild,
                                  sa,
                                  f_ui,
                                  g_ui,
                                  fuv_ui,
                                  d,
                                  rsbound,
                                  p,
                                  max_e,
                                  curr_e + 1,
                                  sub );

    /* we are in the second level from bottom, consider all
       children of this node (bottom nodes) and sieve */
    if (curr_e == max_e) {

      tmpnode = currnode->firstchild;

      while (tmpnode != NULL) {

        /* p|MOD, the inverse doesn't exist */
        if (mpz_divisible_ui_p (rsbound->MOD, p) != 0) {

          /* This shouldn't be called since we didn't record them
             during the lift. Only for safty purpose. */
          if (mpz_divisible_ui_p (rsbound->MOD, pe) != 0) {
            tmpnode2 = tmpnode;
            tmpnode = tmpnode->nextsibling;
#if DEBUG_MULTROOT_LIFT
            fprintf (stderr, "deleting bottomnode (%u, %u) with %u roots in level %u.\n",
                     tmpnode2->u, tmpnode2->v, tmpnode2->nr, curr_e);
#endif
            free_node (&tmpnode2);
            continue;
          }
          /* If current pe-|-MOD, then in B + MOD*x = v (mod pe),
             it is still possible to solve such x. */
          else {
            j_idx = rootsieve_run_multroot_lift_idx
              ( tmpnode->v, rsbound->Bmin,
                rsbound->B, rsbound->MOD, pl, pe );

            i_idx = rootsieve_run_multroot_lift_idx
              ( tmpnode->u, rsbound->Amin,
                rsbound->A, rsbound->MOD, pl, pe );
            step = pe / pl;
          }
        }
        /* p-|-MOD, the inverse exists */
        else {
          j_idx = uv2ij_mod (rsbound->B, rsbound->Bmin, rsbound->MOD, tmpnode->v, pe);
          i_idx = uv2ij_mod (rsbound->A, rsbound->Amin, rsbound->MOD, tmpnode->u, pe);
          step = pe;
        }

        /* be careful about this */
        subtmp = (int16_t) ceil ( (double) sub / (double) pem1  * tmpnode->nr);

        /* For each i block */
        long old_j = j_idx;
        while (i_idx < (long) sa->len_i) {
          j_idx = old_j;
          rootsieve_run_line ( sa->array,
                               i_idx * sa->len_j + sa->len_j - 1,
                               i_idx * sa->len_j + j_idx,
                               step,
                               subtmp );
          i_idx += step;
        }

        tmpnode2 = tmpnode;
        tmpnode = tmpnode->nextsibling;

#if DEBUG_MULTROOT_LIFT
        fprintf (stderr, "deleting bottomnode (%u, %u) with %u roots in level %u AND ",
                 tmpnode2->u, tmpnode2->v, tmpnode2->nr, curr_e);
        fprintf (stderr, "sieving bottomnode %d in steps %u.\n", subtmp, step); // pe
#endif

        free_node (&tmpnode2);
      }
    }

    /* delete current node and move to next sibling. */
    tmpnode = currnode;
    currnode = currnode->nextsibling;
    if (currnode != NULL)
      (currnode->parent)->firstchild = currnode;

    /* p|MOD, the inverse doesn't exist */
    if (mpz_divisible_ui_p (rsbound->MOD, p) != 0) {

      /* This shouldn't be called since we didn't record them
         during the lift. Only for safty purpose. */
      if (mpz_divisible_ui_p (rsbound->MOD, pem1) != 0) {
#if DEBUG_MULTROOT_LIFT
        fprintf (stderr, "deleting (%u, %u) with %u roots in level %u.\n",
                 tmpnode->u, tmpnode->v, tmpnode->nr, curr_e - 1);
#endif
        free_node (&tmpnode);
        continue;
      }
      /* If current pem1-|-MOD, then in B + MOD*x = v (mod pem1),
         it is still possible to solve such x. */
      else {
        j_idx = rootsieve_run_multroot_lift_idx
          ( tmpnode->v, rsbound->Bmin,
            rsbound->B, rsbound->MOD, pl, pem1 );
        i_idx = rootsieve_run_multroot_lift_idx
          ( tmpnode->u, rsbound->Amin,
            rsbound->A, rsbound->MOD, pl, pem1 );
        step = pem1 / pl;
      }
    }
    /* p-|-MOD, the inverse exists */
    else {
      j_idx = uv2ij_mod (rsbound->B, rsbound->Bmin, rsbound->MOD, tmpnode->v, pem1);
      i_idx = uv2ij_mod (rsbound->A, rsbound->Amin, rsbound->MOD, tmpnode->u, pem1);
      step = pem1;
    }

    /* be careful about this */
    subtmp = (int16_t) ceil ( (double) sub / (double) pem1 * (double) p  * tmpnode->nr );

    /* For each i block */
    long old_j = j_idx;
    while (i_idx < (long) sa->len_i) {

      j_idx = old_j;
      rootsieve_run_line ( sa->array,
                           i_idx * sa->len_j + sa->len_j - 1,
                           i_idx * sa->len_j + j_idx,
                           step, subtmp );
      i_idx += step;
    }

#if DEBUG_MULTROOT_LIFT
    fprintf (stderr, "deleting (%u, %u) with %u roots in level %u AND ",
             tmpnode->u, tmpnode->v, tmpnode->nr, curr_e - 1);
    fprintf (stderr, "sieving %d in steps %u.\n", subtmp, step); //pem1
#endif

    free_node (&tmpnode);
  }

  return;
}


/*
  For multiple root r of the (i, j) (mod p)
*/
static inline void
rootsieve_run_multroot ( sievearray_t sa,
                         ropt_poly_t rs,
                         rsbound_t rsbound,
                         unsigned int u,
                         unsigned int i,
                         long j,
                         unsigned int r,
                         unsigned int p,
                         unsigned int e,
                         int16_t sub )
{
  /* sieve f(r) + ((u+i*p)*r + (B + MOD * (j+k*p)))*g(r) for i, k,
     we sieve (mod p) case in the beginning since r would
     be a multiple root for all (j + k*p).
     In the rootsieve_v(), it must be that
     roottype_flag = 1, e = 1; We can sieve starting from j.
     Note if roottype_flag = 2, e must > 1. We are not here. */
  if (e <= 1) {

    /* i is kept the same, next_i denotes the u-dim sieve */
    unsigned int next_i = i;
    long old_j = j;

    /* For each i block */
    while (next_i < sa->len_i) {

      j = old_j;
      rootsieve_run_line ( sa->array,
                           next_i * sa->len_j + sa->len_j - 1,
                           next_i * sa->len_j + j, p, sub );
      next_i += p;
    }
    return;
  }

  /* some variables */
  unsigned int v, pe, *f_ui, *g_ui, *fuv_ui;
  mpz_t tmpz;

  mpz_init_set_ui (tmpz, 0UL);

  /* compute p^e */
  pe = 1;
  for (v = 0; v < e ; v ++)
    pe = pe * p;

  /* use s.p instead of m.p */
  f_ui = (unsigned int*) malloc ((rs->d + 1) * sizeof (unsigned int));
  fuv_ui = (unsigned int*) malloc ((rs->d + 1) * sizeof (unsigned int));
  g_ui = (unsigned int*) malloc ((2) * sizeof (unsigned int));
  if ((f_ui == NULL) || (g_ui == NULL) || (fuv_ui == NULL)) {
    fprintf (stderr, "Error, cannot allocate memory in rootsieve_run_multroot(). \n");
    exit (1);
  }
  reduce_poly_ul (f_ui, rs->f, rs->d, pe);
  reduce_poly_ul (g_ui, rs->g, 1, pe);

  /* j -> v (mod p) */
  ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, j, tmpz);
  v = (unsigned int) mpz_fdiv_ui (tmpz, p);
  // v = (unsigned int) ( (tmp >= 0) ? tmp : tmp + (int) p );

#if DEBUG_MULTROOT_LIFT
  fprintf (stderr, "\n f(%u) + ( %u * %u + %u ) * g(%u) = 0 (mod %u)\n",
           r, u, r, v, r, p);
#endif

  /* we've already known that r is a multiple root for f_{u, v}. */
  node *tmpnode, *root;
  new_tree (&root);
  root = new_node ();
  insert_node (root, &tmpnode, u, v, r, 1, p, p, 2);

  /* lift to higher p^e */
  rootsieve_run_multroot_lift ( root->firstchild,
                                sa,
                                f_ui,
                                g_ui,
                                fuv_ui,
                                rs->d,
                                rsbound,
                                p,
                                e,
                                2,
                                sub );

  /* free, either root itself or root with a level 1 node.
     Note, if e>1, freeing will be done in the lift function;
     However, if e=1, we need to manaully free root->firstchild. */
  free_node (&root);
  tmpnode = NULL;

  free (f_ui);
  free (fuv_ui);
  free (g_ui);
  mpz_clear (tmpz);
}


#define DEBUG_ROOTSIEVE_UV_ONEBLOCK 0
/*
  Root sieve for f + u*x*g(x) + v*g(x)
*/
static inline void
rootsieve_uv_oneblock ( sievearray_t sa,
                        ropt_poly_t rs,
                        rsbound_t rsbound,
                        ropt_param_t rsparam )
{
  unsigned int nbb, p, r, u, v, k, max_e, totnbb, fr_ui, gr_ui, pe = 1;
  unsigned short np;
  long i, bblock_size, i_idx;
  mpz_t tmpz1;
  mpz_init (tmpz1);

  /* some arrays for holding index */
  long j_idx[primes[rsparam->len_p_rs]], j_idx_i0[primes[rsparam->len_p_rs]];
  int16_t subsgl[rsparam->len_p_rs], submul[rsparam->len_p_rs];
  char roottype_flag[rsparam->len_p_rs];
  bblock_size = L1_CACHESIZE;
  totnbb = (unsigned int) ((rsbound->Bmax - rsbound->Bmin + 1) / bblock_size);

#if DEBUG_ROOTSIEVE_UV_ONEBLOCK
  fprintf ( stderr,
            "# Stat: totnbb: %u, bblock_size: %ld, total_size: %ld\n",
            totnbb, bblock_size, (long) totnbb * bblock_size );
#endif

  /* Init index array */
  float subf;
  for ( np = 0; np < rsparam->len_p_rs; np ++ ) {
    p = primes[np];
    subf = (float) p * log ( (float) p) / ((float) p * (float) p - 1.0);
    subsgl[np] = (int16_t) ceil (subf * 1000.0);
    subf = log ( (float) p) / ( (float) p + 1.0);
    submul[np] = (int16_t) ceil (subf * 1000.0);
    roottype_flag[np] = 0;
  }

  /* Idx holders for each r < B */
  for ( r = 0; r < primes[rsparam->len_p_rs]; r ++ ) {
    j_idx[r] = 0;
    j_idx_i0[r] = 0;
  }

  /* For each p < p_bound*/
  for (np = 0; np < rsparam->len_p_rs; np ++) {

    p = primes[np];

    /* e depends on p */
    max_e = (unsigned int) (log (200.0) / log ((double) p));
    pe = 1;
    for (k = 0; k < max_e; k ++)
      pe = pe * p;

    /* pe | MOD in A + MOD*i = u (mod pe), don't need to sieve */
    if (mpz_divisible_ui_p (rsbound->MOD, pe) != 0)
      continue;

    /* For each u < p */
    for (u = 0; u < p; u++) {

      i = 0;

      /* compute u->i in A + MOD*i = u (mod p) */
      if (mpz_divisible_ui_p (rsbound->MOD, p) != 0) {
        if (mpz_fdiv_ui (rsbound->A, p) != u)
          continue;
        else {
          i = -1;
        }
      }
      else {
        /* i >= 0 */
        i = uv2ij_mod (rsbound->A, rsbound->Amin, rsbound->MOD, u, p);
      }

      /* i is kept the same, i_idx denotes the u-dim sieve */
      i_idx = i;

      /* For each i block */
      while (i_idx < (long) sa->len_i) {

        /* For each j block on a i-line */
        for (nbb = 0; nbb < totnbb + 1; nbb ++) {

          /* For each r < p_bound */
          for (r = 0; r < p; r++) {

            /* skip */
            if (mpz_divisible_ui_p(rs->gx[r], p) != 0)
              continue;

            /* u*g(r)^2 - f(r)g'(r) + f'(r)g(r) */
            mpz_mul (tmpz1, rs->gx[r], rs->gx[r]);
            mpz_mul_ui (tmpz1, tmpz1, u);
            mpz_sub (tmpz1, tmpz1, rs->numerator[r]);

            /* if nbb == 0:
               -- if i_idx == i, compute indices;
               -- elseif i_idx != 0, retrieve;
               if nbb != 0:
               -- next; */
            if (nbb == 0) {
              if (i_idx == i) {

                /* reset for all r < p for the first block */
                roottype_flag[r] = 0;
                j_idx[r] = 0;
                j_idx_i0[r] = 0;

                /* compute v in f(r) + u*r*g(r) + v*g(r) = 0 (mod p) */
                fr_ui = mpz_fdiv_ui (rs->fx[r], p);
                gr_ui = mpz_fdiv_ui (rs->gx[r], p);
                v = compute_v_ui (fr_ui, gr_ui, r, u, p);

                /* compute v->j in B + MOD*j = v (mod p) */
                if (mpz_divisible_ui_p (rsbound->MOD, p) != 0) {

                  /* now A + MOD*i = u (mod p)
                     and B + MOD*j = v (mod p) */
                  if (mpz_fdiv_ui (rsbound->B, p) == v) {

                    /* case where MOD has p-valuation smaller than e.
                       - If r is multiple, do the lift;
                       - If r is simple, all slots on the array are equivalent - ignore. */
                    if (mpz_divisible_ui_p(tmpz1, p) == 0)
                      continue;
                    else
                      roottype_flag[r] = 2;
#if DEBUG_ROOTSIEVE_UV_ONEBLOCK
                    fprintf (stderr, "[mult] f(%u) + ( %u * %u + %u ) * g(%u) = 0 (mod %u)\n",
                             r, u, r, v, r, p);
                    fprintf (stderr, " j_idx: %ld\n", j_idx[r]);
#endif
                  }
                  else
                    continue;
                }
                else {
                  roottype_flag[r] = 1;

                  /* u -> i in A + MOD*i = u (mod p) has been computed. */
                  /* v -> j in B + MOD*j = v (mod p) */
                  j_idx_i0[r] = uv2ij_mod (rsbound->B, rsbound->Bmin, rsbound->MOD, v, p);
                  j_idx[r] = j_idx_i0[r];
#if DEBUG_ROOTSIEVE_UV_ONEBLOCK
                  fprintf (stderr, "[simple] f(%u) + ( %u * %u + %u ) * g(%u) = 0 (mod %u)\n",
                           r, u, r, v, r, p);
                  fprintf (stderr, " i_idx: %ld, j_idx: %ld\n", i_idx, j_idx[r]);
#endif
                }
              }
              /* if nbb == 0, but i_idx != i, retrieve info from nbb = 0 */
              else {
                /* don't change roottype_flag[np] */
                j_idx[r] = j_idx_i0[r];
              }
            }
            /* nbb != 0, need to find j_idx[r]'s j-position with respect to a line */
            else {
              j_idx[r] %= sa->len_j;
            }

            /* cases where we don't want to sieve at all */
            if (roottype_flag[r] == 0)
              continue;

            /* r is a simple root for current (u, v, p). Then r
               is simple root for (u, v+i*p, p) for all i. */
            if (mpz_divisible_ui_p(tmpz1, p) == 0) {

              if (nbb != totnbb) {
                j_idx[r] = rootsieve_run_line (
                  sa->array,
                  i_idx * sa->len_j + bblock_size * (nbb + 1),
                  i_idx * sa->len_j + j_idx[r],
                  p,
                  subsgl[np] );
              }
              else {
                j_idx[r] = rootsieve_run_line (
                  sa->array,
                  i_idx * sa->len_j + (rsbound->Bmax - rsbound->Bmin),
                  i_idx * sa->len_j + j_idx[r],
                  p,
                  subsgl[np] );
              }
            }
            /* r is multiple root for current u, p */
            else {

              /* don't sieve in block, assume cm << cs */
              if (nbb == 0 && i_idx == i) {

                /* Two cases:
                   If roottype_flag = 1, then MOD != 0 (mod p), everything is normal;
                   If MOD = 0 (mod p), then at nbb = 0, roottype_flag = 2;
                   Also note, must use u (mod pe) instead of u (mod p) */
                rootsieve_run_multroot ( sa,
                                         rs,
                                         rsbound,
                                         u,
                                         i,
                                         j_idx[r],
                                         r,
                                         p,
                                         max_e,
                                         submul[np] );
              }
            }
          } // next r
        } // next j-block

        if (i_idx == -1)
          break;
        else {
          i_idx += p;
        }

      } // next i-block
    } // next u
  } // next p

  mpz_clear (tmpz1);
}


/*
  Rootsieve for f + (u*x +v)*g for each sublattice
*/
static inline double
rootsieve_uv ( ropt_poly_t rs,
               rsbound_t rsbound,
               ropt_param_t rsparam,
               MurphyE_pq *global_E_pqueue,
               mpz_t *fuv,
               mpz_t *guv,
               int w,
               int verbose )
{
  /* for each sieving array, we first look at the E
     of #TOPALPHA polynomials which have best alpha */
  const int TOPALPHA = TOPALPHA_EACH_SIEVEARRAY;
  const int TOPE = TOPE_EACH_SUBLATTICE;
  int l;
  long tmpBmax, tmpBmin;
  unsigned long j, array_mem_size, B_block_size, len_B, len_A;
  double MurphyE;
  mpz_t tmpv, tmpu;
  mpz_init (tmpv);
  mpz_init (tmpu);

  len_A = (unsigned long) (rsbound->Amax - rsbound->Amin + 1);
  len_B = (unsigned long) (rsbound->Bmax - rsbound->Bmin + 1);

  /* tune sieve mode */
  if (verbose == 1) {
    array_mem_size = TUNE_SIEVEARRAY_SIZE; // len_A should == 1, see ropt_param.c
  }
  /* polyselect2* mode */
  else if (verbose == 0) {
    array_mem_size = MAX_SHORT_SIEVEARRAY_SIZE;
    /* if this is all ready too long */
    if ( len_B < (double) array_mem_size / len_A ) {
      array_mem_size = len_B * len_A;
    }
  }
  /* ropt_main mode */
  else {
    array_mem_size = MAX_LONG_SIEVEARRAY_SIZE;
    /* if this is all ready too long */
    if ( len_B < (double) array_mem_size / len_A ) {
      array_mem_size = len_B * len_A;
    }
  }

  if (array_mem_size <= (unsigned long) len_A) {
    fprintf (stderr, "Error: Amax - Amin + 1 = %lu is too long. This happend when A is too large while B is too small. To be fixed in future.\n", rsbound->Amax - rsbound->Amin);
    exit(1);
  }

  B_block_size = (unsigned long) ( (double) array_mem_size / (double) len_A );

  /* E priority queue for each sublattice */
  MurphyE_pq *E_pqueue;
  new_MurphyE_pq (&E_pqueue, TOPE);
  rootscore_pq *alpha_pqueue;
  new_rootscore_pq (&alpha_pqueue, TOPALPHA);
  sievearray_t sa;
  sievearray_init (sa, len_A, B_block_size);
  /* get the correct value */
  array_mem_size = len_A * B_block_size;

  
  /* for each i -> each u = A + MOD * i */
  tmpBmax = rsbound->Bmax;
  tmpBmin = rsbound->Bmin;

  long k = 0;
  int st = cputime ();
  do {
    k ++;

    /* fo reach block of size B_block_size */
    rsbound->Bmax = rsbound->Bmin + B_block_size - 1;
    if (rsbound->Bmax > tmpBmax)
      rsbound->Bmax = tmpBmax;

#if DEBUG
    fprintf (stderr, "[%ld, %ld] ", rsbound->Bmin, rsbound->Bmax);
#endif

    /* root sieve for v. Note, rsbound with changed Bmin and Bmax will be used in rootsieve_v(). */
    rootsieve_uv_oneblock (sa, rs, rsbound, rsparam);

    /* record top n slots */
    for (j = 0; j < array_mem_size; j++) {
      insert_rootscore_pq ( alpha_pqueue,
                            (j - j % (B_block_size)) / B_block_size,
                            j % (B_block_size),
                            sa->array[j] );

      // LEAVE THIS FOR DEBUG //
      /* if (MAT[j] < -2900) { */
      /*  ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, i, tmpu); */
      /*  ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, j, tmpv); */
      /*  gmp_fprintf (stderr, "MAT[]: %d, u: %Zd, v: %Zd, i: %lu, j: %lu\n", MAT[j], tmpu, tmpv, i, j); */
      /* } */
    }

    /* output polynomials (put them into the MurphyE priority queue) */
    for (l = 1; l < alpha_pqueue->used; l ++) {

      ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, alpha_pqueue->i[l], tmpu);
      ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, alpha_pqueue->j[l], tmpv);
      compute_fuv_mp (fuv, rs->f, rs->g, rs->d, tmpu, tmpv);
      mpz_set (guv[0], rs->g[0]);
      mpz_set (guv[1], rs->g[1]);

      //optimize (fuv, rs->d, guv, 0, 0); not correct for deg 6 polynomial
      optimize_aux (fuv, rs->d, guv, 0, 0, CIRCULAR);
#if DEBUG
      gmp_fprintf ( stderr,
                    "\n# Found #%2d (u=%Zd, v=%Zd)",
                    l, tmpu, tmpv );
#endif

      MurphyE = print_poly_fg (fuv, guv, rs->d, rs->n, 0);
      insert_MurphyE_pq ( E_pqueue, w, tmpu, tmpv, MurphyE );
    }

    /* next j */
    rsbound->Bmin = rsbound->Bmax + 1;

    /* reset sieve array and priority queue */
    sievearray_reset (sa);
    reset_rootscore_pq (alpha_pqueue);

  } while (rsbound->Bmin < tmpBmax);

  /* recover */
  rsbound->Bmax = tmpBmax;
  rsbound->Bmin = tmpBmin;

  /* free priority queue and sieving array */
  sievearray_free (sa);
  free_rootscore_pq (&alpha_pqueue);

  /* output polynomials */
  MurphyE = 0.0;

  for (l = 1; l < E_pqueue->used; l ++) {

    if (verbose == 2) {
      gmp_fprintf ( stderr, "\n# Found (%dth) (w=%d, u=%Zd, v=%Zd) gives E = %1.2e",
                    l, E_pqueue->w[l], E_pqueue->u[l], E_pqueue->v[l], E_pqueue->E[l]);
    }

    /* insert E scores to a global queue (before optimization) */
    insert_MurphyE_pq ( global_E_pqueue, E_pqueue->w[l], E_pqueue->u[l], E_pqueue->v[l], E_pqueue->E[l] );
    compute_fuv_mp (fuv, rs->f, rs->g, rs->d, E_pqueue->u[l], E_pqueue->v[l]);
    mpz_set (guv[0], rs->g[0]);
    mpz_set (guv[1], rs->g[1]);
    optimize_aux (fuv, rs->d, guv, 0, 0, CIRCULAR);
    MurphyE += print_poly_fg (fuv, guv, rs->d, rs->n, verbose);
  }

  MurphyE = MurphyE / (double) (E_pqueue->used - 1);

  if (verbose == 2) {
    gmp_fprintf ( stderr,
                  "\n# Stat: ave_MurphyE of top %ld polynomials: %1.2e (on sublattice %Zd, %Zd)\n",
                  E_pqueue->used - 1,
                  MurphyE,
                  rsbound->A, rsbound->B );

    fprintf ( stderr, "# Stat: root sieve took %dms\n",
              cputime () - st );
  }

  free_MurphyE_pq (&E_pqueue);
  mpz_clear (tmpv);
  mpz_clear (tmpu);

  return MurphyE;
}


/*
  Stage 2: root sieve on the sublattice points;
  For each sublattice, call rootsieve_uv().
*/
static inline double
ropt_stage2_run ( ropt_poly_t rs,
                  rsbound_t rsbound,
                  ropt_param_t rsparam,
                  MurphyE_pq *global_E_pqueue,
                  int w,
                  int verbose )
{
  int i;
  mpz_t *fuv, *guv;
  double ave_MurphyE = 0.0;

  /* fuv is f+(u*x+v)*g */
  fuv = (mpz_t*) malloc ((rs->d + 1) * sizeof (mpz_t));
  guv = (mpz_t*) malloc (2 * sizeof (mpz_t));
  if (fuv == NULL || guv == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in ropt_stage2_run().\n");
    exit (1);
  }
  for (i = 0; i <= rs->d; i++)
    mpz_init_set (fuv[i], rs->f[i]);
  for (i = 0; i < 2; i++)
    mpz_init_set (guv[i], rs->g[i]);

  /* alpha values on the subllatice primes */
  compute_fuv_mp (fuv, rs->f, rs->g, rs->d, rsbound->A, rsbound->B);

  /* On this sublattice, sieve (i, j) */
  ave_MurphyE = rootsieve_uv ( rs,
                               rsbound,
                               rsparam,
                               global_E_pqueue,
                               fuv,
                               guv,
                               w,
                               verbose );

  for (i = 0; i <= rs->d; i++)
    mpz_clear (fuv[i]);
  for (i = 0; i < 2; i++)
    mpz_clear (guv[i]);
  free (fuv);
  free (guv);

  return ave_MurphyE;
}


/*
  Stage 2: root sieve on the sublattice points.
*/
double
ropt_stage2 ( ropt_poly_t rs,
              ropt_param_t rsparam,
              param_t param,
              MurphyE_pq *global_E_pqueue,
              int w,
              mpz_t u,
              mpz_t v,
              mpz_t mod,
              int verbose )
{
  double ave_MurphyE = 0.0;
  rsbound_t rsbound;

  /* init rsbound */
  rsbound_init (rsbound);

  /* set root sieve bounds, depending on U, V bounds */
  rsbound_setup_AB_bound (rsbound, rsparam, param, mod, verbose);

  if (verbose == 2)
    fprintf ( stderr, "# Info: sieving matrix size: [%ld, %ld] x [%ld, %ld]\n",
              rsbound->Amin, rsbound->Amax,
              rsbound->Bmin, rsbound->Bmax );

  rsbound_setup_sublattice ( rsbound,
                             u,
                             v,
                             mod );

  /* compute exact sieving bounds UV given size AB depending
     on current A, B, MOD */
  if (verbose == 2)
    rsbound_print (rsbound);

  /* root sieve */
  ave_MurphyE = ropt_stage2_run ( rs,
                                  rsbound,
                                  rsparam,
                                  global_E_pqueue,
                                  w,
                                  verbose );
  /* free rsbound */
  rsbound_free (rsbound);
  return ave_MurphyE;
}
