/**
 * @file ropt_stage1.c
 * Called by ropt.c to find congruence classes.
 */

#include "ropt_stage1.h"

/*
  Find good sublattice, the lifted cases.
*/
static inline void
find_sublattice_lift ( node *firstchild,
                       listnode **top,
                       unsigned int * f_ui,
                       unsigned int * g_ui,
                       unsigned int * fuv_ui,
                       int d,
                       unsigned int p,
                       char e,
                       char curr_e )
{
  /* recursion end */
  if (firstchild == NULL || curr_e > e)
    return;

  char l;
  unsigned int i, j, k, nroots, pe, pem1, fr, gr;
  node *currnode, *tmpnode = NULL, *tmpnode2 = NULL;

  /* compute p^e */
  pem1 = 1;
  for (l = 0; l < curr_e - 1; l ++)
    pem1 = pem1 * p;
  pe = pem1 * p;

  /* loop until all siblings are checked. */
  currnode = firstchild;
  while (currnode != NULL) {

    /*
      printf("-----\n");
      printf("p: %u, e: %d -> %d\n", p, curr_e - 1, curr_e);
      printf("(u, v) pair: (%u, %u) has roots: \n",
      currnode->u, currnode->v);
      for (i = 0; i < (currnode->nr); i++)
      printf("\tr[%u]: %u\n", i, currnode->r[i]);
      printf("-----\n");
    */

    /* compute f_uv(x) and then evaluate it at r. */
    compute_fuv_ui (fuv_ui, f_ui, g_ui, d, currnode->u, currnode->v, pe);

    /* loop all roots */
    for (nroots = 0; nroots < currnode->nr; nroots++) {

      gr = eval_poly_ui_mod (g_ui, 1, currnode->r[nroots], pe);
      if (gr % p == 0)
        continue;

      /* If the root is multiple */
      if (currnode->roottype[nroots] == 2) {

        fr = eval_poly_ui_mod (fuv_ui, d, currnode->r[nroots], pe);
        fr = fr / pem1;

        /* solve on fr + gr*A = 0 (mod p), where A = i*r + j. */
        fr = (unsigned int) solve_lineq (fr, gr, 0, p);

        /* we want to solve (i, j) in  fr = i*r + j (mod p).
           - if r is not invertible, loop i with fixed j;
           - otherwise, loop j to solve i. */
        if (currnode->r[nroots] % p == 0) {

          for (i = 0; i < p; i ++) {

#if DEBUG_FIND_SUBLATTICE
            fprintf (stderr, "fr: %u, r: %u,  (i, j): (%u, %u) -> (%u, %u) (non-invertible, multiple)\n",
                     fr, currnode->r[nroots], i, fr, currnode->u + pem1 * i,
                     currnode->v + pem1 * fr);
#endif
            /* r is a multiple root, add r + k * p^{e-1}. */
            for (k = 0; k < p; k ++) {
              insert_node (currnode, &tmpnode, currnode->u + pem1 * i,
                           currnode->v + pem1 * fr, currnode->r[nroots] + k * pem1, curr_e, p, pe, 2);
            }

            /* count the lifted single roots for any lifted pairs (u, v). Note
               the lifted single roots will not be computed actually. */
            for (k = 0; k < (currnode->nr); k++) {
              if (currnode->roottype[k] == 1) {
                insert_node (currnode, &tmpnode, currnode->u + pem1 * i,
                             currnode->v + pem1 * fr, currnode->r[k], curr_e, p, pe, 1);
              }
            }
          }
        }
        else {
          for (j = 0; j < p; j ++) {

            /* given j, solve i in  fr = i*r + j (mod p). */
            i = solve_lineq (j, currnode->r[nroots], fr, p);
#if DEBUG_FIND_SUBLATTICE
            fprintf (stderr, "fr: %u, r: %u,  (i, j): (%u, %u) -> (%u, %u) (invertible, multiple)\n",
                     fr, currnode->r[nroots], i, j, currnode->u + pem1 * i,
                     currnode->v + pem1 * j);
#endif
            /* r is a multiple root, add r + k * p^{e-1}. */
            for (k = 0; k < p; k ++) {
              insert_node (currnode, &tmpnode, currnode->u + pem1 * i,
                           currnode->v + pem1 * j, currnode->r[nroots] + k * pem1, curr_e, p, pe, 2);
            }

            /* count the lifted single roots for any lifted pairs (u, v). Note
               the lifted single roots will not be computed actually. */
            for (k = 0; k < (currnode->nr); k++) {
              if (currnode->roottype[k] == 1) {
                insert_node (currnode, &tmpnode, currnode->u + pem1 * i,
                             currnode->v + pem1 * j, currnode->r[k], curr_e, p, pe, 1);
              }
            }
          }
        }
      }
    }  // next root of current (u, v)

    /* recursieve to next level, curr_e + 1 */
    find_sublattice_lift ( currnode->firstchild,
                           top,
                           f_ui,
                           g_ui,
                           fuv_ui,
                           d,
                           p,
                           e,
                           curr_e + 1 );

    /* If current node is the 2nd bottom leave, add the bottom level leaves
       with highest valuations to the list and delete them. */
    if (curr_e == e) {
      tmpnode = currnode->firstchild;
      while (tmpnode != NULL) {
        insert_listnode (top, tmpnode->u, tmpnode->v, tmpnode->val, e);
        tmpnode2 = tmpnode;
        tmpnode = tmpnode->nextsibling;

#if DEBUG_FIND_SUBLATTICE
        fprintf (stderr, "DEBUG_FIND_SUBLATTICE (2nd bottom): p: %u, (%u, %u), val: %f, e: %d, max_e: %d\n",
                 p, tmpnode2->u, tmpnode2->v, tmpnode2->val, curr_e, e);
#endif

        free_node (&tmpnode2);
      }
    }

    /* delete current node and move to next sibling. */
    tmpnode = currnode;
    currnode = currnode->nextsibling;
    if (currnode != NULL)
      (currnode->parent)->firstchild = currnode;

#if DEBUG_FIND_SUBLATTICE
    fprintf (stderr, "DEBUG_FIND_SUBLATTICE (bottom): p: %u, (%u, %u), val: %f, e: %d, max_e: %d\n",
             p, tmpnode->u, tmpnode->v, tmpnode->val, curr_e, e);
#endif

    free_node (&tmpnode);
  } // next sibling of current node

  return;
}


/*
  Find sublattices, the base case. Note, for the (u, v)
  pairs, we consider all the simple + multi roots.
*/
static inline void
find_sublattice ( listnode **top,
                  rsstr_t rs,
                  unsigned int p,
                  char e )
{
  unsigned int pe, r, u, v, fx_ui, gx_ui;
  char i;
  node *currnode, *root;
  mpz_t tmp;
  mpz_init (tmp);

  /* compute p^e */
  pe = 1;
  for (i = 0; i < e; i ++)
    pe = pe * p;

  /* new (u, v, val) tree */
  new_tree (&root);
  root = new_node ();

  /* for each root 0 <= r < p  */
  for (r = 0; r < p; r ++) {

    /* skip these */
    if (mpz_divisible_ui_p(rs->gx[r], p) != 0)
      continue;

    /* use single precision */
    fx_ui = (unsigned int) mpz_fdiv_ui (rs->fx[r], p);
    gx_ui = (unsigned int) mpz_fdiv_ui (rs->gx[r], p);

    for (u = 0; u < p; u ++) {

      /* u*g(r)^2 - f(r)g'(r) + f'(r)g(r) */
      mpz_mul (tmp, rs->gx[r], rs->gx[r]);
      mpz_mul_ui (tmp, tmp, u);
      mpz_sub (tmp, tmp, rs->numerator[r]);

      /* compute v in f(r) + u*r*g(r) + v*g(r) = 0 (mod p) */
      v =  compute_v_ui (fx_ui, gx_ui, r, u, p);

      /* simple root */
      if (mpz_divisible_ui_p(tmp, p) == 0) {

        insert_node (root, &currnode, u, v, r, 1, p, p, 1);
      }
      else {
        insert_node (root, &currnode, u, v, r, 1, p, p, 2);
      }
    }
  }

  /* If e == 1, stop */
  if (e == 1) {
    node *tmpnode;
    currnode = root->firstchild;
    while (currnode != NULL) {
      insert_listnode (top, currnode->u, currnode->v, currnode->val, 1);
      tmpnode = currnode;
      currnode = currnode->nextsibling;
      free_node (&tmpnode);
    }
    tmpnode = NULL;
  }
  /* if e > 1, lift to higher p^e */
  else {
    /* find_sublattice_lift() only consider those pairs with at least one multiple
       , hence we need to consider those (u,v) which solely have single roots */
    node *tmpnode = NULL, *lastnode = NULL;
    unsigned int j, c;
    currnode = root->firstchild;
    while (currnode != NULL) {

      c = 1;
      for (j = 0; j < currnode->nr; j ++)
        if (currnode->roottype[j] == 2)
          c = 0;

      /* case when (u, v) only has single roots */
      if (c == 1) {

        insert_listnode (top, currnode->u, currnode->v, currnode->nr / ( (double) p - 1), 1);

        /* delete this node */
        tmpnode = currnode;
        currnode = currnode->nextsibling;
        if (lastnode != NULL)
          lastnode->nextsibling = currnode;
        else
          (currnode->parent)->firstchild = currnode;
        free_node (&tmpnode);
      }
      else {
        lastnode = currnode;
        currnode = currnode->nextsibling;
      }
    }

    /* data struct for the lift */
    unsigned int *f_ui, *g_ui, *fuv_ui;
    f_ui = (unsigned int*) malloc ( (rs->d + 1) * sizeof (unsigned int) );
    fuv_ui = (unsigned int*) malloc ( (rs->d + 1) * sizeof (unsigned int) );
    g_ui = (unsigned int*) malloc ( 2 * sizeof (unsigned int) );
    if ( (f_ui == NULL) ||
         (g_ui == NULL) ||
         (fuv_ui == NULL) ) {
      fprintf (stderr, "Error, cannot allocate memory in find_sublattice(). \n");
      exit (1);
    }
    /* compute f (mod pe) */
    reduce_poly_ul (f_ui, rs->f, rs->d, pe);
    reduce_poly_ul (g_ui, rs->g, 1, pe);

    find_sublattice_lift ( root->firstchild,
                           top,
                           f_ui,
                           g_ui,
                           fuv_ui,
                           rs->d,
                           p,
                           e,
                           2 );

    /* clear */
    free (f_ui);
    free (fuv_ui);
    free (g_ui);
  }

  mpz_clear (tmp);
  free_node(&root);
}


#if 0
/*
  partition.
*/
static long
quick_sort_2d_ld_partition  ( long **array,
                              const unsigned short dim,
                              const long l,
                              const long h )
{
  long pivot_index = l + (h - l) / 2, tmp = 0, i, first_large = l;

  tmp = array [h][dim];
  array [h][dim] = array [pivot_index][dim];
  array [pivot_index][dim] = tmp;
  tmp = array [h][1-dim];
  array [h][1-dim] = array [pivot_index][1-dim];
  array [pivot_index][1-dim] = tmp;

  for (i = l; i <= h; i ++) {
    if ( abs(array[i][dim]) < abs(array[h][dim]) ) {

      tmp = array [i][dim];
      array [i][dim] = array [first_large][dim];
      array [first_large][dim] = tmp;
      tmp = array [i][1-dim];
      array [i][1-dim] = array [first_large][1-dim];
      array [first_large][1-dim] = tmp;

      first_large ++;
    }
  }

  tmp = array [h][dim];
  array [h][dim] = array [first_large][dim];
  array [first_large][dim] = tmp;
  tmp = array [h][1-dim];
  array [h][1-dim] = array [first_large][1-dim];
  array [first_large][1-dim] = tmp;

  return first_large;
}


/*
  Do a quick_selection since we only want the nbest.
*/
static void
quick_sort_2d_ld ( long **array,
                   const unsigned short dim,
                   const long l,
                   const long h,
                   const long nbest )
{
  /* quick selection */
  if (l < h) {
    long pivot = quick_sort_2d_ld_partition  (array, dim, l, h);
    if (pivot > nbest)
      quick_sort_2d_ld (array, dim, l, pivot - 1, nbest);
    if (pivot < nbest)
      quick_sort_2d_ld (array, dim, pivot + 1, h, nbest);
  }
}
#endif


/*
  Compute crt and add (u, v) to queue.
*/
static inline void
return_all_sublattices_crt ( rsparam_t rsparam,
                             unsigned int *ind,
                             unsigned int ***individual_sublattices,
                             sublattice_pq *pqueue )
{
  int i;
  unsigned int pe[rsparam->tlen_e_sl], j, e;
  mpz_t tmpp1, tmpp2, tmpu1, tmpu2, tmpv, re;
  mpz_init (tmpp1);
  mpz_init (tmpp2);
  mpz_init (tmpu1);
  mpz_init (tmpu2);
  mpz_init (tmpv);
  mpz_init (re);

  /* compute pe[] */
  for (j = 0; j < rsparam->tlen_e_sl; j ++) {
    pe[j] = 1;
    /* note, the actual e can be different from the rsparam->e_sl[] */
    for (e = 0; e < individual_sublattices[j][ind[j]][2]; e ++)
      pe[j] *= primes[j];
  }

  /* compute modulus based on pe[], put it temperorily to rsparam->modulus */
  mpz_set_ui (rsparam->modulus, 1);
  for (j = 0; j < rsparam->tlen_e_sl; j ++) {
    mpz_mul_ui (rsparam->modulus, rsparam->modulus, pe[j]);
  }

  /* compute u */
  mpz_set_ui (tmpu1, individual_sublattices[0][ind[0]][0]);
  mpz_set_ui (tmpp1, pe[0]);

  for (i = 1; i < rsparam->tlen_e_sl; i ++) {

    mpz_set_ui (tmpu2, individual_sublattices[i][ind[i]][0]);
    mpz_set_ui (tmpp2, pe[i]);

    crt_pair_mp ( tmpu1,
                  tmpp1,
                  tmpu2,
                  tmpp2,
                  re );

    mpz_mul_ui (tmpp1, tmpp1, pe[i]);
    mpz_set (tmpu1, re);
  }

  /* re < 0 by construction */
  mpz_sub (re, tmpu1, rsparam->modulus);

  /* if u is good, compute v */
  if ( mpz_cmp_ui (tmpu1, rsparam->global_u_bound_rs) <= 0 ||
       mpz_cmpabs_ui (re, rsparam->global_u_bound_rs) <= 0 ) {

    /* compute v */
    mpz_set_ui (tmpu2, individual_sublattices[0][ind[0]][1]);
    mpz_set_ui (tmpp1, pe[0]);
    for (i = 1; i < rsparam->tlen_e_sl; i ++) {

      mpz_set_ui (tmpv, individual_sublattices[i][ind[i]][1]);
      mpz_set_ui (tmpp2, pe[i]);

      crt_pair_mp ( tmpu2,
                    tmpp1,
                    tmpv,
                    tmpp2,
                    re );

      mpz_mul_ui (tmpp1, tmpp1, pe[i]);
      mpz_set (tmpu2, re);
    }

    /* (u, v) pair in (tmpu1, tmpu2) */
    if (mpz_cmp_ui (tmpu1, rsparam->global_u_bound_rs) > 0)
      mpz_sub (tmpu1, tmpu1, rsparam->modulus);

    /* insert this node */
    insert_sublattice_pq ( pqueue, tmpu1, tmpu2, rsparam->modulus );
  }

  mpz_clear (tmpp1);
  mpz_clear (tmpp2);
  mpz_clear (tmpu1);
  mpz_clear (tmpu2);
  mpz_clear (tmpv);
  mpz_clear (re);
}

/*
  static void
  SORT ( listnode **top,
  rsstr_t rs,
  unsigned int p,
  unsigned int e )
*/


/*
  Return all sublattices by calling CRT, where for each sublattice,
  the seperate (mod p) valuations are the best.
*/
static int
return_all_sublattices ( rsstr_t rs,
                         rsparam_t rsparam,
                         sublattice_pq *pqueue,
                         int verbose )
{
  /* At least consider three primes */
  if (rsparam->tlen_e_sl < 2) {
    fprintf ( stderr,
              "Error: At least consider two primes (2, 3) in return_all_sublattice. \n" );
    exit(1);
  }

  unsigned short i, k, tmp_e_sl;
  unsigned int j, size[rsparam->tlen_e_sl], tsize[rsparam->tlen_e_sl],
    ind[rsparam->tlen_e_sl], ***individual_sublattices;
  unsigned long count;
  listnode *top, *tmp;

  /* Sublattice[i][length][] save (u, v) for prime[i] */
  individual_sublattices = (unsigned int ***)
    malloc ( rsparam->tlen_e_sl * sizeof (unsigned int **) );
  if ( individual_sublattices == NULL ) {
    fprintf (stderr,
             "Error, cannot allocate memory in return_all_sublattices(). \n");
    exit (1);
  }

  /* For each prime[i], lift the polynomial
     and return (u, v) with best score. */
  for (i = 0; i < rsparam->tlen_e_sl; i ++) {

    new_list (&top);
    /* find individual sublattices */
    find_sublattice ( &top, rs,
                      primes[i], rsparam->e_sl[i] );

    tsize[i] = count_list (top);

    /* If individual list has 0 len, this set of parameters fails. */
    if (tsize[i] == 0) {
      for (k = 0; k < rsparam->tlen_e_sl; k ++) {
        for (j = 0; j < size[i]; j ++) {
          if (individual_sublattices[k][j] != NULL)
            free (individual_sublattices[k][j]);
        }
        if (individual_sublattices[k] != NULL)
          free (individual_sublattices[k]);
      }
      if (individual_sublattices != NULL)
        free (individual_sublattices);
      return 0;
    }

    /* If the list is too long (> ncrts_sl), consider lower the e
       in p^e. This seems beneficial since we don't want to throw
       away any good sublattices. Might be slower but safer */
    tmp_e_sl = rsparam->e_sl[i];
    while (tsize[i] > rsparam->ncrts_sl) {

      free_list (&top);
      new_list (&top);
      tmp_e_sl --;

      /* find individual sublattices */
      find_sublattice ( &top, rs,
                        primes[i], tmp_e_sl );

      tsize[i] = count_list (top);
    }

    /* now, tsize[i] must be <= rsparam->ncrts_sl and >0 */
    size[i] = tsize[i];

    /* allocate array for the i-th prime. */
    (individual_sublattices)[i] = (unsigned int **)
      malloc ( size[i] * sizeof(unsigned int *) );

    if ( (individual_sublattices)[i] == NULL ) {
      fprintf (stderr, "Error, cannot allocate memory in return_all_sublattices(). \n");
      exit (1);
    }

    tmp = top;
    for (j = 0; j < size[i]; j ++) {
      (individual_sublattices)[i][j] =
        (unsigned int *) malloc ( 3 * sizeof(unsigned int));
      if ( individual_sublattices[i][j] == NULL ) {
        fprintf (stderr, "Error, cannot allocate memory in return_all_sublattices(). \n");
        exit (1);
      }
      individual_sublattices[i][j][0] = tmp->u;
      individual_sublattices[i][j][1] = tmp->v;
      individual_sublattices[i][j][2] = (unsigned int) tmp->e;
      tmp = tmp->next;
      //fprintf (stderr, "SUBLATTICE: #%lu, (%lu, %lu)\n", j, individual_sublattices[i][j][0], individual_sublattices[i][j][1]);
    }

    //print_list (top);

    if (verbose == 2)
      fprintf ( stderr,
                "# Info: p: %2u, max_e: %2u (old_e: %2u), list_size: %6u, size_cutoff: %6u\n",
                primes[i], tmp_e_sl, rsparam->e_sl[i], tsize[i], size[i] );

    free_list (&top);
  }

  /* Loop over combinations of all arrays. This is awkward.
     We could map 0 ... \prod pe[i] to the indices of the
     arrays in the price of using more arithmetic. */

  /* 2, 3, 5, 7 */
  if (rsparam->tlen_e_sl == 4) {
    for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < size[3]; ind[3] ++)
            return_all_sublattices_crt ( rsparam,
                                         ind,
                                         individual_sublattices,
                                         pqueue );
  }
  /* 2, 3, 5, 7, 11 */
  else if (rsparam->tlen_e_sl == 5) {
    for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < size[4]; ind[4] ++)
              return_all_sublattices_crt ( rsparam,
                                           ind,
                                           individual_sublattices,
                                           pqueue );
  }
  /* 2, 3, 5, 7, 11, 13 */
  else if (rsparam->tlen_e_sl == 6) {
    for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < size[4]; ind[4] ++)
              for (ind[5] = 0; ind[5] < size[5]; ind[5] ++)
                return_all_sublattices_crt ( rsparam,
                                             ind,
                                             individual_sublattices,
                                             pqueue );
  }
  /* 2, 3, 5, 7, 11, 13, 17 */
  else if (rsparam->tlen_e_sl == 7) {
    for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < size[4]; ind[4] ++)
              for (ind[5] = 0; ind[5] < size[5]; ind[5] ++)
                for (ind[6] = 0; ind[6] < size[6]; ind[6] ++)
                  return_all_sublattices_crt ( rsparam,
                                               ind,
                                               individual_sublattices,
                                               pqueue );
  }
  /* 2, 3, 5, 7, 11, 13, 17, 19 */
  else if (rsparam->tlen_e_sl == 8) {
    for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < size[4]; ind[4] ++)
              for (ind[5] = 0; ind[5] < size[5]; ind[5] ++)
                for (ind[6] = 0; ind[6] < size[6]; ind[6] ++)
                  for (ind[7] = 0; ind[7] < size[7]; ind[7] ++)
                    return_all_sublattices_crt ( rsparam,
                                                 ind,
                                                 individual_sublattices,
                                                 pqueue );
  }
  /* 2, 3, 5, 7, 11, 13, 17, 19, 23 */
  else if (rsparam->tlen_e_sl == 9) {
    for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < size[4]; ind[4] ++)
              for (ind[5] = 0; ind[5] < size[5]; ind[5] ++)
                for (ind[6] = 0; ind[6] < size[6]; ind[6] ++)
                  for (ind[7] = 0; ind[7] < size[7]; ind[7] ++)
                    for (ind[8] = 0; ind[8] < size[8]; ind[8] ++)
                      return_all_sublattices_crt ( rsparam,
                                                   ind,
                                                   individual_sublattices,
                                                   pqueue );
  }
  /* 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 */
  else if (rsparam->tlen_e_sl == 10) {
    for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
          for (ind[3] = 0; ind[3] < size[3]; ind[3] ++)
            for (ind[4] = 0; ind[4] < size[4]; ind[4] ++)
              for (ind[5] = 0; ind[5] < size[5]; ind[5] ++)
                for (ind[6] = 0; ind[6] < size[6]; ind[6] ++)
                  for (ind[7] = 0; ind[7] < size[7]; ind[7] ++)
                    for (ind[8] = 0; ind[8] < size[8]; ind[8] ++)
                      for (ind[9] = 0; ind[9] < size[9]; ind[9] ++)
                        return_all_sublattices_crt ( rsparam,
                                                     ind,
                                                     individual_sublattices,
                                                     pqueue );
  }
  /* 2, 3, 5 */
  else if (rsparam->tlen_e_sl == 3) {
    for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
        for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
          return_all_sublattices_crt ( rsparam,
                                       ind,
                                       individual_sublattices,
                                       pqueue );
  }
  /* 2, 3 */
  else if (rsparam->tlen_e_sl == 2) {
    for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
        return_all_sublattices_crt ( rsparam,
                                     ind,
                                     individual_sublattices,
                                     pqueue );
  }
  /* too aggressive */
  else {
    fprintf (stderr, "Error, number of primes in \"-e\" (len_e_sl) should be between 2 to 10\n");
    exit (1);
  }


  /* info */
  if (verbose == 2) {
    /* Compute rsparam->modulus */
    count = 1;
    for (i = 0; i < rsparam->tlen_e_sl; i ++) {
      count *= size[i];
    }
    fprintf (stderr, "# Info: computed %lu CRTs\n", count);
  }

  /* clearence */
  for (i = 0; i < rsparam->tlen_e_sl; i ++) {
    for (j = 0; j < size[i]; j ++) {
      free (individual_sublattices[i][j]);
    }
    free (individual_sublattices[i]);
  }
  free (individual_sublattices);

  return 1;
}


/*
  Call return_all_sublattice() to return good sublattices.

  Choose the "rsbound->nbest_sl" best ones among them if there
  are more quantities than it. The "best" property is ranked
  by the size of u.
*/
static inline int
return_best_sublattice ( rsstr_t rs,
                         rsparam_t rsparam,
                         sublattice_pq *pqueue,
                         int verbose )
{
  unsigned long i = 0UL, global_u_bound_rs_tmp;

  /* find the actual rsparam->len_e_sl, excluding 0's */
  rsparam->tlen_e_sl = rsparam->len_e_sl;
  global_u_bound_rs_tmp = rsparam->global_u_bound_rs;
  for (i = 0; i < rsparam->len_e_sl; i ++) {
    if (rsparam->e_sl[i] == 0)
      rsparam->tlen_e_sl --;
  }

  int ret = return_all_sublattices ( rs,
                                     rsparam,
                                     pqueue,
                                     verbose );

  /* If failed, return. Some individual sublattices has length 0 */
  if (ret == 0) {
    return -1;
  }

  /* If no sublattice is found with u < rsparam->global_u_bound_rs,
     then we try to enlarge u bound. However, it might be better
     to enlarge e_sl[] to allow to check more sublattices. */
  int count = 1;
  while (pqueue->used == 1) {

    if (rsparam->global_u_bound_rs < LONG_MAX) {
      rsparam->global_u_bound_rs *= 2;
      if (verbose == 2) {
        fprintf ( stderr,
                  "# Warn: Not enough sublattice classes. Reset \"rsparam->global_u_bound_rs = %lu\" (#%d)\n",
                  rsparam->global_u_bound_rs, count );
      }
    }
    else
      return -1;

    return_all_sublattices ( rs,
                             rsparam,
                             pqueue,
                             verbose );
    count ++;
  }

  /* info */
  if (verbose == 2) {
    fprintf ( stderr,
              "# Info: found %d sublattices with small u, where |u| < %lu (\"rsparam->global_u_bound_rs\")\n",
              pqueue->used - 1, rsparam->global_u_bound_rs);
  }

  /* recover, for deg 6 poly */
  rsparam->global_u_bound_rs = global_u_bound_rs_tmp;

  return 1;
}


/*
  Stage 1: record good sublattices in "alpha_pqueue".
*/
int
ropt_stage1 ( rsstr_t rs,
              rsparam_t rsparam,
              sub_alpha_pq *alpha_pqueue,
              int verbose,
              int w )
{
  int st, i, re;
  mpz_t *fuv;
  double alpha_lat;
#if DEBUG
  double skew, logmu;
#endif
  sublattice_pq *pqueue;

  /* fuv is f+(u*x+v)*g */
  fuv = (mpz_t*) malloc ((rs->d + 1) * sizeof (mpz_t));
  if (fuv == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in ropt_stage1.\n");
    exit (1);
  }
  for (i = 0; i <= rs->d; i++)
    mpz_init_set (fuv[i], rs->f[i]);

  /* priority queue */
  new_sublattice_pq (&pqueue, rsparam->nbest_sl);

  /* return the nbest good sublattices to queue ranked by the size of u */
  st = cputime ();
  re = return_best_sublattice ( rs,
                                rsparam,
                                pqueue,
                                verbose );

  /* failed, free queue and ret */
  if (re == -1) {
    free_sublattice_pq (&pqueue);
    for (i = 0; i <= rs->d; i++) {
      mpz_clear (fuv[i]);
    }
    free (fuv);
    return -1;
  }

  if (verbose == 2)
    gmp_fprintf ( stderr, "# Info: find best sublattices took %dms\n",
                  cputime () - st );

  /* filter all sublattices into another, global queue ranked by the
     sublattice alphas. */
  for (i = 1; i < pqueue->used; i ++) {

    compute_fuv_mp (fuv, rs->f, rs->g, rs->d, pqueue->u[i], pqueue->v[i]);
    //alpha_lat = get_biased_alpha_affine (fuv, rs->d, primes[rsparam->tlen_e_sl - 1]);
    alpha_lat = get_alpha (fuv, rs->d, 2000);
#if DEBUG
    skew = L2_skewness (fuv, rs->d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    logmu = L2_lognorm (fuv, rs->d, skew, DEFAULT_L2_METHOD);
    gmp_fprintf ( stderr, "# Info: insert sublattice #%4d, (w, u, v): (%d, %Zd, %Zd), alpha: %.2f, logmu: %.2f\n",
                  i,
                  w,
                  pqueue->u[i],
                  pqueue->v[i],
                  alpha_lat,
                  logmu );
#endif
    /* insert to a global priority queue */
    insert_sub_alpha_pq ( alpha_pqueue,
                          w,
                          pqueue->u[i],
                          pqueue->v[i],
                          pqueue->modulus[i],
                          alpha_lat );
  }

  /* free priority queue */
  free_sublattice_pq (&pqueue);
  for (i = 0; i <= rs->d; i++) {
    mpz_clear (fuv[i]);
  }
  free (fuv);

  return 0;
}


#if WANT_TUNE
/*
  auxiliary for tune
*/
static inline double
rsparam_tune_aux ( rsstr_t rs,
                   rsparam_t rsparam,
                   param_t param,
                   sub_alpha_pq *alpha_pqueue,
                   int nbest_sl_tunecut,
                   int rot_w )
{
  int i, re;
  double ave_MurphyE, ave2_MurphyE = 0.0;

  /* re-init the pqeueu */
  alpha_pqueue->used = 1;

  re = ropt_stage1 ( rs,
                               rsparam,
                               alpha_pqueue,
                               0,
                               rot_w );
  if (re == -1) {
    return -1;
  }

  /* use another array to save the priority queue */
  int used = alpha_pqueue->used - 1;
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

    /*
      gmp_fprintf ( stderr, "# Tune: #%4d sublattice (w, u, v): (%d, %Zd, %Zd) (mod %Zd), alpha: %.2f\n",
      i + 1,
      w[i],
      u[i],
      v[i],
      mod[i],
      sub_alpha[i] );
    */
  }

  /* For each sublattice, do the root sieve */
  MurphyE_pq *global_E_pqueue;
  new_MurphyE_pq (&global_E_pqueue, 4);

  re = 0;
  for (i = used - 1; i >= 0; i --) {

    if (re > nbest_sl_tunecut)
      break;

    /* TBC
       gmp_fprintf ( stderr,
       "\n# Tune: Sieve on sublattice (# %2d), (w, u, v): (%d, %Zd, %Zd)  (mod %Zd)\n# Info: affine_alpha: %.2f, proj_alpha: %.2f, exp_min_alpha: %.2f\n",
       i + 1,
       w[i],
       u[i],
       v[i],
       mod[i],
       sub_alpha[i],
       rs->alpha_proj,
       rsparam->exp_min_alpha_
    */
    ave_MurphyE = ropt_stage2 ( rs,
                                rsparam,
                                param,
                                global_E_pqueue,
                                w[i],
                                u[i],
                                v[i],
                                mod[i],
                                1 ); // tune mode, this is 1 in order to to set the correct length (short) of sieve array.
    ave2_MurphyE += ave_MurphyE;
    re ++;
  }

  free_MurphyE_pq (&global_E_pqueue);
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
  Tune parameters for find_sublattice().
*/
double
rsparam_tune ( rsstr_t rs,
               rsparam_t rsparam,
               param_t param,
               int num_trials,
               int w,
               int verbose )
{
  unsigned short i, j, best_j = 0, tmp_e_sl[rsparam->len_e_sl], k;
  unsigned int p, pearr[rsparam->len_e_sl];
  double ave_MurphyE = 0, best_MurphyE = 0;

  for (i = 0; i < rsparam->len_e_sl; i ++) {
    p = primes[i];
    pearr[i] = 1;
    tmp_e_sl[i] = rsparam->e_sl[i];
    for (j = 0; j < rsparam->e_sl[i]; j ++) {
      pearr[i] *= p;
    }
  }

  if (verbose != 0) {
    /* first set of parameters */
    fprintf (stderr, "# Tune:");
    for (i = 0; i < rsparam->len_e_sl; i ++) {
      fprintf (stderr, " %u^%u=%u ", primes[i],
               rsparam->e_sl[i], pearr[i]);
    }
  }
  /* queue */
  sub_alpha_pq *alpha_pqueue;
  new_sub_alpha_pq (&alpha_pqueue, rsparam->nbest_sl);

  best_MurphyE = rsparam_tune_aux ( rs,
                                    rsparam,
                                    param,
                                    alpha_pqueue,
                                    3, // nbest_sl_tunecut
                                    w );
  alpha_pqueue->used = 1;

  if (verbose != 0) {
    if (best_MurphyE == -1)
      fprintf (stderr, " ave_MurphyE: failed\n");
    else
      fprintf (stderr, " ave_MurphyE: %.3e\n", best_MurphyE);
  }
  /* Other sets of parameters, try next "num_trials" sets
     of parameters and keep the best one. */
  char flag;
  for (j = 0; j < num_trials; j ++) {

    flag = 0;

    /* loops */
    for (i = 1; i < rsparam->len_e_sl; i ++) {
      if (pearr[i] * primes[i] < pearr[i-1]) {
        pearr[i] *= primes[i];
        rsparam->e_sl[i] += 1;
        flag = 1;
        break;
      }
    }

    if (flag == 0) {
      pearr[0] *= primes[0];
      rsparam->e_sl[0] += 1;
    }

    if (verbose != 0) {
      /* some info and tune */
      fprintf (stderr, "# Tune:");
      for (i = 0; i < rsparam->len_e_sl; i ++) {
        fprintf (stderr, " %u^%u=%u ", primes[i],
                 rsparam->e_sl[i], pearr[i]);
      }
    }
    /* test sieve */
    ave_MurphyE = rsparam_tune_aux ( rs,
                                     rsparam,
                                     param,
                                     alpha_pqueue,
                                     8, // nbest_sl_tunecut
                                     w );
    alpha_pqueue->used = 1;

    if (verbose != 0) {
      if (ave_MurphyE == -1)
        fprintf (stderr, " ave_MurphyE: failed\n");
      else
        fprintf (stderr, " ave_MurphyE: %.3e\n", ave_MurphyE);
    }
    if (ave_MurphyE >= best_MurphyE) {
      best_MurphyE = ave_MurphyE;
      //printf ("current best: %.6e\n", best_MurphyE);
      best_j = j;
      for (k = 0; k < rsparam->len_e_sl; k ++) {
        //printf ("exp: %u\n", rsparam->e_sl[k]);
        tmp_e_sl[k] = rsparam->e_sl[k];
      }
    }
  }

  free_sub_alpha_pq (&alpha_pqueue);

  /* finally, save best parameters back */
  for (k = 0; k < rsparam->len_e_sl; k ++) {
    rsparam->e_sl[k] = tmp_e_sl[k];
    //printf ("e: %u\n", rsparam->e_sl[k]);
  }

  if (verbose != 0) {
    /* output best parameters */
    fprintf (stderr, "# Tune: best parameters ");
    for (i = 0; i < rsparam->len_e_sl; i ++) {
      fprintf (stderr, "%u:%u ", primes[i], rsparam->e_sl[i]);
    }
    fprintf (stderr, "\n");
  }
  return best_MurphyE;
}

#endif
