/**
 * @file ropt_tree.c
 * Liftings and priority queue structs.
 */

#include "ropt_tree.h"

#if 0
/*
  Print the info for the node
*/
static inline void
print_node ( node *pnode )
{
  fprintf(stderr, "(%u,%u):%u:%d:(%.2f)\n",
          pnode->u, pnode->v, pnode->nr, pnode->e, pnode->val);
  /* int i; */
  /* for (i = 0; i < pnode->nr; i++) */
  /* printf ("pnode->r[%d]: %lu\n", i, pnode->r[i]); */
}


/*
  Print a tree, non-recursive. Two styles.
  - If level == 0, print the whole tree;
  - If level == i, print nodes at height i.
*/
static inline void
print_tree ( node *root,
             int level )
{
  if (level < 0) {
    fprintf(stderr,"Error: level >= 0. \n");
    exit(1);
  }

  int i, currlevel = 0;
  node *ptr = root;
  printf ("\n");

  if (!ptr)  /* if empty */
    return;

  /* find leftmost */
  while (ptr->firstchild) {
    ptr = ptr->firstchild;
    ++ currlevel;
  }
  while (ptr) {
    if (level != 0) {
      if (currlevel == level) {
        for (i = 0; i < currlevel; ++i)
          printf("-   ");
        print_node(ptr);
      }
    }
    else {
      for (i = 0; i < currlevel; ++i)
        printf("-   ");
      print_node(ptr);
    }

    if (ptr->nextsibling) {
      ptr = ptr->nextsibling;
      //++ level; // aligned layout, uncommont for sloped layout
      while (ptr->firstchild) {
        ptr = ptr->firstchild;
        ++ currlevel;
      }
    }
    else {
      ptr = ptr->parent;  /* move up */
      -- currlevel;
    }
  }
  printf ("\n");
  return;
}
#endif


/*
  Create new empty node.
*/
node *
new_node ( void )
{
  node *pnode;
  pnode = (node *) malloc(sizeof(node));
  if (pnode == NULL) {
    fprintf(stderr,"Error: malloc failed in new_node()\n");
    exit(1);
  }
  pnode->firstchild = pnode->nextsibling = pnode->parent = NULL;
  pnode->r = NULL;
  pnode->roottype = NULL;
  pnode->u = pnode->v = pnode->nr = pnode->e = 0;
  pnode->val = 0.0;
  return pnode;
}


/*
  Initialise new empty tree.
*/
void
new_tree ( node **root )
{
  *root = NULL;
}


/*
  Free current node ptr.
*/
void
free_node ( node **ptr )
{
  if (*ptr) {
    if ((*ptr)->r)
      free ((*ptr)->r);
    if ((*ptr)->roottype)
      free ((*ptr)->roottype);
    (*ptr)->parent = (*ptr)->firstchild =
      (*ptr)->nextsibling = NULL;
    free (*ptr);
  }
}


#if 0
/*
  Free tree given by root.
*/
void
free_tree ( node *root )
{
  node *ptr = root;

  /* if empty */
  if (!ptr)
    return;

  if (ptr->firstchild)
    free_tree (ptr->firstchild);

  if (ptr->nextsibling)
    free_tree (ptr->nextsibling);

  //print_node (ptr);
  free_node (&ptr);
  return;
}
#endif


/*
  Malloc or realloc (double) the memory space of pnode->r[] and pnode->roottype[].
*/
void
alloc_r_node ( node *pnode )
{
  //assert (pnode->alloc <= pnode->nr);

  if (pnode->nr == 0) {
    pnode->r = (unsigned int *)
      malloc ( sizeof (unsigned int) );
    pnode->roottype = (char *)
      malloc ( sizeof (char) );
    pnode->alloc = 1;
  }
  else {
    pnode->r = (unsigned int *)
      realloc ( pnode->r, 2 * pnode->nr * sizeof (unsigned int) );
    pnode->roottype = (char *)
      realloc ( pnode->roottype, 2 * pnode->nr * sizeof (char) );
    pnode->alloc = (pnode->alloc * 2);
  }

  if (pnode->r == NULL) {
    fprintf (stderr, "Error, cannot reallocate memory in alloc_r_node().\n");
    exit (1);
  }
}


/*
  Insert child node to the current node, where *parent is
  the ptr to the parent of the current node.

  - If (u, v) exits in any childrent of the parent;
  -- Check if r exists;
  --- If not, add r and/or "is_multiple=k".
  - If (u, v) doesnot exit;
  -- Add a node with (u, v, r, curr_e, val) and/or "is_multiple=k".
*/
#define DEBUG_INSERT_NODE 0
void
insert_node ( node *parent,
              node **currnode,
              unsigned int u,
              unsigned int v,
              unsigned int r,
              char curr_e,
              unsigned int p,
              unsigned int pe,
              char k )
{
  unsigned int i;
  node *lastnode = NULL;
  node *nextnode = NULL;
  nextnode = parent->firstchild;
  lastnode = nextnode;

  /* add r to some existing node (u, v) */
  while ((nextnode != NULL)) {

    if ((nextnode->u == u) && (nextnode->v == v)) {

      /* if (u,v) pair exists, see whether
         the root r already exists in the node. */
      for (i = 0; i < nextnode->nr; i ++) {
        if ( nextnode->r[i] == r )
          return;
      }

      /* otherwise, insert this root */
      if (nextnode->alloc <= nextnode->nr)
        alloc_r_node (nextnode);
      nextnode->roottype[nextnode->nr] = k;
      nextnode->r[nextnode->nr] = r;
      nextnode->nr += 1;

      if (k == 2)
        nextnode->val += 1.0 / (double) pe;
      else {
        if (curr_e == 1)
          nextnode->val += 1.0 / ((double)p - 1.0);
      }

#if DEBUG_INSERT_NODE
      if (pe == 7)
        fprintf ( stderr, "add to (u: %u, v: %u), nr: %u, r: %u, e: %d, val: %f\n",
                  u, v, nextnode->nr, nextnode->r[nextnode->nr-1], curr_e, nextnode->val );
#endif

      break;
    }
    lastnode = nextnode;
    nextnode = nextnode->nextsibling;
  }

  /* add new (u, v) node */
  if (nextnode == NULL) {

    *currnode = new_node();

    /* has left sibling? */
    if (lastnode != NULL)
      lastnode->nextsibling = (*currnode);
    else
      parent->firstchild = (*currnode);

    (*currnode)->parent = parent;
    (*currnode)->u = u;
    (*currnode)->v = v;

    alloc_r_node (*currnode);
    (*currnode)->r[0] = r;
    (*currnode)->roottype[0] = k;
    (*currnode)->nr += 1;
    (*currnode)->e = curr_e;
    /* such (u, v) is new, hence we inheritate the val from its parent. */

    if (k == 2)
      (*currnode)->val += 1.0 / (double) pe + (*currnode)->parent->val;
    else {
      if (curr_e == 1)
        (*currnode)->val += 1.0 / ((double)p - 1.0);
      else
        (*currnode)->val = (*currnode)->parent->val;
    }

#if DEBUG_INSERT_NODE
    if (pe == 7)
      fprintf ( stderr, "creating (u: %u, v: %u), nr: %u, r: %u, e: %d, val: %f\n", u, v,
                (*currnode)->nr, (*currnode)->r[0], curr_e, (*currnode)->val );
#endif
  }
  //print_tree(parent, 0);
}

/*
  Initialise new list for (u, v) p-valuations.
*/
void
new_list ( listnode **top )
{
  *top = NULL;
}


/*
  free non-empty listnode.
*/
void
free_listnode ( listnode **pplistnode )
{
  if (*pplistnode)
    free (*pplistnode);
}


/*
  Del the list pointed by top.
*/
void
free_list ( listnode **pptop )
{
  if (*pptop != NULL) {
    listnode *ptr;
    while ( (*pptop)->next != NULL) {
      ptr = (*pptop);
      (*pptop) = (*pptop)->next;
      free_listnode (&ptr);
    }
    free_listnode (pptop);
  }
}


#if 0
/*
  Print the info for the listnode
*/
void
print_listnode ( listnode *plistnode )
{
  printf("(u, v): (%u,%u), val: %.2f, e: %d\n", plistnode->u, plistnode->v, plistnode->val, plistnode->e);
}


/*
  Print the list pointed by top.
*/
void
print_list ( listnode *ptop )
{
  if (ptop) {
    while ( ptop->next ) {
      print_listnode (ptop);
      ptop = ptop->next;
    }
    print_listnode (ptop);
  }
}
#endif

/*
  Return the length of the list.
*/
unsigned long
count_list ( listnode *ptop)
{
  unsigned long s = 0UL;
  if (ptop) {
    while ( ptop->next ) {
      ptop = ptop->next;
      s ++;
    }
    //print_listnode (ptop);
    s ++;
  }
  return s;
}


/*
  Create new empty listnode.
*/
listnode *
new_listnode ( unsigned int u,
               unsigned int v,
               double val,
               char e )
{
  listnode *plistnode = NULL;
  plistnode = (listnode *) malloc (sizeof(listnode));
  if (plistnode == NULL) {
    fprintf(stderr,"Error: malloc failed in new_listnode()\n");
    exit(1);
  }
  plistnode->u = u;
  plistnode->v = v;
  plistnode->val = val;
  plistnode->e = e;
  plistnode->next = NULL;
  return plistnode;
}


/*
  Insert a listnode to current list <- top, only record best one.
*/
void
insert_listnode ( listnode **top,
                  unsigned int u,
                  unsigned int v,
                  double val,
                  char e )
{
  listnode *newlistnode;

  /* if empty, create node */
  if ( (*top) == NULL ) {
    newlistnode = new_listnode (u, v, val, e);
    (*top) = newlistnode;
    (*top)->next = NULL;
  }
  else {
    /* if income has better val, delete and add */
    if ( (*top)->val < val) {
      newlistnode = new_listnode (u, v, val, e);
      free_list (top);
      (*top) = newlistnode;
      (*top)->next = NULL;
    }
    /* if income has equal val, add */
    else if ( (*top)->val == val) {
      newlistnode = new_listnode (u, v, val, e);
      newlistnode->next = (*top);
      (*top) = newlistnode;
    }
  }
  newlistnode = NULL;
}


/*
  Some indices, note the queue is shifted by 1.
*/
inline int
pq_parent ( int i )
{
  return (i>>1);
}

inline int
pq_leftchild ( int i )
{
  return (i<<1);
}

inline int
pq_rightchild ( int i )
{
  return ( (i << 1) + 1 );
}


/*
  Create priority queue for best sublattices of length len.
*/
void
new_sublattice_pq ( sublattice_pq **ppqueue,
                    unsigned long len )
{
  if ( len < 3 ) {
    fprintf(stderr,"Error: len < 3 in new_sublattice_pq()\n");
    exit(1);
  }

  (*ppqueue) = (sublattice_pq *) malloc (sizeof (sublattice_pq));
  if ( (*ppqueue) == NULL) {
    fprintf(stderr,"Error: malloc failed in new_sublattice_pq()\n");
    exit(1);
  }

  (*ppqueue)->len = len;

  (*ppqueue)->u = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->v = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->modulus = (mpz_t *) malloc (len* sizeof (mpz_t));
  if ( (*ppqueue)->u == NULL
       || (*ppqueue)->v == NULL
       ||  (*ppqueue)->v == NULL ) {
    fprintf(stderr,"Error: malloc failed in new_sublattice_pq()\n");
    exit(1);
  }

  int i;
  for (i = 0; i < (*ppqueue)->len; i++)
  {
    mpz_init ( (*ppqueue)->u[i] );
    mpz_init ( (*ppqueue)->v[i] );
    mpz_init ( (*ppqueue)->modulus[i] );
  }

  mpz_set_str ( (*ppqueue)->u[0], "340282366920938463463374607431768211456", 10 ); // 2^128
  mpz_set_ui ( (*ppqueue)->v[0], 0 );
  mpz_set_ui ( (*ppqueue)->modulus[0], 0 );

  (*ppqueue)->used = 1; // u[0] and v[0] are null elements

}


/*
  Create priority queue for best sublattices of length len.
*/
void
free_sublattice_pq ( sublattice_pq **ppqueue )
{
  int i;
  for (i = 0; i < (*ppqueue)->len; i++)
  {
    mpz_clear ( (*ppqueue)->u[i] );
    mpz_clear ( (*ppqueue)->v[i] );
    mpz_clear ( (*ppqueue)->modulus[i] );
  }

  free ( (*ppqueue)->u );
  free ( (*ppqueue)->v );
  free ( (*ppqueue)->modulus );
  free ( *ppqueue );
}


/*
  Sift-up to add, if the queue is not full.
*/
static inline void
insert_sublattice_pq_up ( sublattice_pq *pqueue,
                          mpz_t u,
                          mpz_t v,
                          mpz_t mod )
{
  int i;

  for ( i = pqueue->used;
        mpz_cmpabs (u, pqueue->u[pq_parent(i)]) > 0; i /= 2 ) {

    mpz_set ( pqueue->u[i], pqueue->u[pq_parent(i)] );
    mpz_set ( pqueue->v[i], pqueue->v[pq_parent(i)] );
    mpz_set ( pqueue->modulus[i], pqueue->modulus[pq_parent(i)] );
  }

  mpz_set (pqueue->u[i], u);
  mpz_set (pqueue->v[i], v);
  mpz_set (pqueue->modulus[i], mod);
  pqueue->used ++;
}


/*
  Sift-down, if the heap is full.
*/
static inline void
insert_sublattice_pq_down ( sublattice_pq *pqueue,
                            mpz_t u,
                            mpz_t v,
                            mpz_t mod )
{
  int i, l;

  for (i = 1; i*2 < pqueue->used; i = l) {

    l = (i << 1);

    /* right > left ? */
    if ( (l+1) < pqueue->used &&
         mpz_cmpabs (pqueue->u[l+1], pqueue->u[l]) > 0 )
      l ++;

    /* switch larger child with parent */
    if ( mpz_cmpabs (pqueue->u[l], u) > 0 ) {
      mpz_set (pqueue->u[i], pqueue->u[l]);
      mpz_set (pqueue->v[i], pqueue->v[l]);
      mpz_set (pqueue->modulus[i], pqueue->modulus[l]);
    }
    else
      break;
  }

  mpz_set (pqueue->u[i], u);
  mpz_set (pqueue->v[i], v);
  mpz_set (pqueue->modulus[i], mod);
}


/*
  Insert to the priority queue.
*/
void
insert_sublattice_pq ( sublattice_pq *pqueue,
                       mpz_t u,
                       mpz_t v,
                       mpz_t mod )
{
  //gmp_fprintf (stderr, "# Debug: inserting (%Zd, %Zd), used: %d, len: %d\n", u, v, pqueue->used, pqueue->len);

  /* queue is full,  */
  if (pqueue->len == pqueue->used) {
    if ( mpz_cmpabs (pqueue->u[1], u) > 0 ) {
      insert_sublattice_pq_down (pqueue, u, v, mod);
    }
  }

  /* queue is not full, sift-up */
  else if (pqueue->len > pqueue->used) {
    insert_sublattice_pq_up (pqueue, u, v, mod);
  }
  else {
    fprintf(stderr,"Error: error (pqueue->len < pqueue->used) in insert_sublattice_pq()\n");
    exit(1);
  }
}


/*
  Create priority queue for best root scores (among MAT[] for a sublattice).
*/
void
new_rootscore_pq ( rootscore_pq **ppqueue,
                   unsigned long len )
{
  if ( len < 3 ) {
    fprintf(stderr,"Error: len < 3 in new_rootscore_pq()\n");
    exit(1);
  }

  (*ppqueue) = (rootscore_pq *) malloc (sizeof (rootscore_pq));
  if ( (*ppqueue) == NULL) {
    fprintf(stderr,"Error: malloc failed in new_rootscore_pq()\n");
    exit(1);
  }

  (*ppqueue)->len = len;

  (*ppqueue)->i = (long *) malloc (len* sizeof (long));
  (*ppqueue)->j = (long *) malloc (len* sizeof (long));
  (*ppqueue)->alpha = (int16_t *) malloc (len* sizeof (int16_t));

  if ( (*ppqueue)->i == NULL ||
       (*ppqueue)->j == NULL ||
       (*ppqueue)->alpha == NULL ) {
    fprintf(stderr,"Error: malloc failed in new_rootscore_pq()\n");
    exit(1);
  }

  /* i[0] and j[0] are null elements */
  (*ppqueue)->i[0] = 0;
  (*ppqueue)->j[0] = 0;
  (*ppqueue)->alpha[0] = INT16_MAX;
  (*ppqueue)->used = 1;

}


/*
  Reset priority queue for best root scores (among MAT[] for a sublattice).
*/
void
reset_rootscore_pq ( rootscore_pq *pqueue )
{

  if ( pqueue == NULL) {
    fprintf(stderr,"Error: malloc failed in new_rootscore_pq()\n");
    exit(1);
  }

  if ( pqueue->i == NULL ||
       pqueue->j == NULL ||
       pqueue->alpha == NULL ) {
    fprintf(stderr,"Error: malloc failed in new_rootscore_pq()\n");
    exit(1);
  }

  /* i[0] and j[0] are null elements */
  pqueue->i[0] = 0;
  pqueue->j[0] = 0;
  pqueue->alpha[0] = INT16_MAX;
  pqueue->used = 1;
}


/*
  Free
*/
void
free_rootscore_pq ( rootscore_pq **ppqueue )
{
  free ( (*ppqueue)->i );
  free ( (*ppqueue)->j );
  free ( (*ppqueue)->alpha );
  free ( *ppqueue );
}


/*
  Sift-up to add, if the queue is not full.
*/
inline void
insert_rootscore_pq_up ( rootscore_pq *pqueue,
                         long i,
                         long j,
                         int16_t alpha )
{
  int k;

  for ( k = pqueue->used;
        alpha > pqueue->alpha[pq_parent(k)];
        k /= 2 )
  {
    pqueue->i[k] = pqueue->i[pq_parent(k)];
    pqueue->j[k] = pqueue->j[pq_parent(k)];
    pqueue->alpha[k] = pqueue->alpha[pq_parent(k)];
  }

  pqueue->i[k] = i;
  pqueue->j[k] = j;
  pqueue->alpha[k] = alpha;

  pqueue->used ++;
}


/*
  Sift-down, if the heap is full.
*/
inline void
insert_rootscore_pq_down ( rootscore_pq *pqueue,
                           long i,
                           long j,
                           int16_t alpha )
{
  int k, l;

  for (k = 1; k*2 < pqueue->used; k = l) {

    l = (k << 1);

    /* right > left ? */
    if ( (l+1) < pqueue->used &&
         (pqueue->alpha[l+1] > pqueue->alpha[l]) )
      l ++;

    /* switch larger child with parent */
    if ( pqueue->alpha[l] > alpha ) {
      pqueue->i[k] = pqueue->i[l];
      pqueue->j[k] = pqueue->j[l];
      pqueue->alpha[k] = pqueue->alpha[l];
    }
    else
      break;
  }

  pqueue->i[k] = i;
  pqueue->j[k] = j;
  pqueue->alpha[k] = alpha;
}


/*
  Insert to the priority queue.
*/
void
insert_rootscore_pq ( rootscore_pq *pqueue,
                      long i,
                      long j,
                      int16_t alpha )
{
  /* fprintf (stderr, "# Debug: inserting (%ld, %ld), alpha: %d, used: %d, len: %d\n", */
  /*   i, j, alpha, pqueue->used, pqueue->len); */

  /* queue is full,  */
  if (pqueue->len == pqueue->used) {
    if ( alpha < pqueue->alpha[1] ) {

      insert_rootscore_pq_down (pqueue, i, j, alpha);
    }
  }
  /* queue is not full, sift-up */
  else if (pqueue->len > pqueue->used) {
    insert_rootscore_pq_up (pqueue, i, j, alpha);
  }
  else {
    fprintf(stderr,"Error: error (pqueue->len < pqueue->used) in insert_sublattice_pq()\n");
    exit(1);
  }
}


/*
  Create priority queue for MurphyE
*/
void
new_MurphyE_pq ( MurphyE_pq **ppqueue,
                 unsigned long len )
{
  if ( len < 3 ) {
    fprintf(stderr,"Error: len < 3 in new_MurphyE_pq()\n");
    exit(1);
  }

  (*ppqueue) = (MurphyE_pq *) malloc (sizeof (MurphyE_pq));
  if ( (*ppqueue) == NULL) {
    fprintf(stderr,"Error: malloc failed in new_MurphyE_pq()\n");
    exit(1);
  }

  (*ppqueue)->len = len;
  (*ppqueue)->w = (int *) malloc (len* sizeof (int));
  (*ppqueue)->u = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->v = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->E = (double *) malloc (len* sizeof (double));

  if ( (*ppqueue)->u == NULL ||
       (*ppqueue)->v == NULL ||
       (*ppqueue)->w == NULL ||
       (*ppqueue)->E == NULL ) {
    fprintf(stderr,"Error: malloc failed in new_MurphyE_pq()\n");
    exit(1);
  }

  int i;
  for (i = 0; i < (*ppqueue)->len; i++)
  {
    mpz_init ( (*ppqueue)->u[i] );
    mpz_init ( (*ppqueue)->v[i] );
  }

  mpz_set_ui ( (*ppqueue)->u[0], 0 );
  mpz_set_ui ( (*ppqueue)->v[0], 0 );
  (*ppqueue)->w[0] = 0;
  (*ppqueue)->E[0] = -DBL_MAX; // E should be positive, so larger than 0 anyway.

  (*ppqueue)->used = 1;
}


/*
  Free
*/
void
free_MurphyE_pq ( MurphyE_pq **ppqueue )
{

  int i;
  for (i = 0; i < (*ppqueue)->len; i++)
  {
    mpz_clear ( (*ppqueue)->u[i] );
    mpz_clear ( (*ppqueue)->v[i] );
  }

  free ( (*ppqueue)->w );
  free ( (*ppqueue)->u );
  free ( (*ppqueue)->v );
  free ( (*ppqueue)->E );
  free ( *ppqueue );

}


/*
  Sift-up to add, if the queue is not full.
*/
inline void
insert_MurphyE_pq_up ( MurphyE_pq *pqueue,
                       int w,
                       mpz_t u,
                       mpz_t v,
                       double E )
{
  int k;

  for ( k = pqueue->used;
        E < pqueue->E[pq_parent(k)];
        k /= 2 )
  {

    mpz_set ( pqueue->u[k], pqueue->u[pq_parent(k)] );
    mpz_set ( pqueue->v[k], pqueue->v[pq_parent(k)] );
    pqueue->w[k] = pqueue->w[pq_parent(k)];
    pqueue->E[k] = pqueue->E[pq_parent(k)];
  }

  mpz_set (pqueue->u[k], u);
  mpz_set (pqueue->v[k], v);
  pqueue->w[k] = w;
  pqueue->E[k] = E;

  pqueue->used ++;
}


/*
  Sift-down, if the heap is full.
*/
inline void
insert_MurphyE_pq_down ( MurphyE_pq *pqueue,
                         int w,
                         mpz_t u,
                         mpz_t v,
                         double E )
{
  int k, l;

  for (k = 1; k*2 < pqueue->used; k = l) {

    l = (k << 1);

    /* right < left ? */
    if ( (l+1) < pqueue->used &&
         (pqueue->E[l+1] < pqueue->E[l]) )
      l ++;

    /* switch smaller child with parent */
    if ( pqueue->E[l] < E ) {

      mpz_set (pqueue->u[k], pqueue->u[l]);
      mpz_set (pqueue->v[k], pqueue->v[l]);
      pqueue->w[k] = pqueue->w[l];
      pqueue->E[k] = pqueue->E[l];
    }
    else
      break;
  }

  mpz_set (pqueue->u[k], u);
  mpz_set (pqueue->v[k], v);
  pqueue->w[k] = w;
  pqueue->E[k] = E;
}


/*
  Insert to the priority queue.
*/
void
insert_MurphyE_pq ( MurphyE_pq *pqueue,
                    int w,
                    mpz_t u,
                    mpz_t v,
                    double E )
{

  /* queue is full,  */
  if (pqueue->len == pqueue->used) {
    if ( E > pqueue->E[1] ) {
      insert_MurphyE_pq_down (pqueue, w, u, v, E);
    }
  }
  /* queue is not full, sift-up */
  else if (pqueue->len > pqueue->used) {
    insert_MurphyE_pq_up (pqueue, w, u, v, E);
  }
  else {
    fprintf(stderr,"Error: error (pqueue->len < pqueue->used) in insert_MurphyE_pq()\n");
    exit(1);
  }
}



/*
  Create priority queue for sublattices with best alpha.
*/
void
new_sub_alpha_pq ( sub_alpha_pq **ppqueue,
                   unsigned long len )
{
  if ( len < 3 ) {
    fprintf(stderr,"Error: len < 3 in new_sub_alpha_pq()\n");
    exit(1);
  }

  (*ppqueue) = (sub_alpha_pq *) malloc (sizeof (sub_alpha_pq));
  if ( (*ppqueue) == NULL) {
    fprintf(stderr,"Error: malloc failed in new_sub_alpha_pq()\n");
    exit(1);
  }

  (*ppqueue)->len = len;

  (*ppqueue)->w = (int *) malloc (len* sizeof (int));
  (*ppqueue)->sub_alpha = (double *) malloc (len* sizeof (double));
  (*ppqueue)->u = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->v = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->modulus = (mpz_t *) malloc (len* sizeof (mpz_t));

  if ( (*ppqueue)->w == NULL ||
       (*ppqueue)->sub_alpha == NULL ||
       (*ppqueue)->modulus == NULL ||
       (*ppqueue)->u == NULL ||
       (*ppqueue)->v == NULL ) {
    fprintf(stderr,"Error: malloc failed in new_sub_slpha_pq()\n");
    exit(1);
  }

  int i;
  for (i = 0; i < (*ppqueue)->len; i++)
  {
    mpz_init ( (*ppqueue)->u[i] );
    mpz_init ( (*ppqueue)->v[i] );
    mpz_init ( (*ppqueue)->modulus[i] );
  }

  mpz_set_ui ( (*ppqueue)->u[0], 0 );
  mpz_set_ui ( (*ppqueue)->v[0], 0 );
  mpz_set_ui ( (*ppqueue)->modulus[0], 0 );
  (*ppqueue)->sub_alpha[0] = DBL_MAX;
  (*ppqueue)->w[0] = 0;

  (*ppqueue)->used = 1;
}


/*
  free
*/
void
free_sub_alpha_pq ( sub_alpha_pq **ppqueue )
{
  int i;
  for (i = 0; i < (*ppqueue)->len; i++)
  {
    mpz_clear ( (*ppqueue)->u[i] );
    mpz_clear ( (*ppqueue)->v[i] );
    mpz_clear ( (*ppqueue)->modulus[i] );
  }

  free ( (*ppqueue)->u );
  free ( (*ppqueue)->v );
  free ( (*ppqueue)->modulus );
  free ( (*ppqueue)->w );
  free ( (*ppqueue)->sub_alpha );
  free ( *ppqueue );
}


/*
  Sift-up to add, if the queue is not full.
*/
inline void
insert_sub_alpha_pq_up ( sub_alpha_pq *pqueue,
                         int w,
                         mpz_t u,
                         mpz_t v,
                         mpz_t modulus,
                         double alpha )
{
  int k;

  for ( k = pqueue->used;
        alpha > pqueue->sub_alpha[pq_parent(k)];
        k /= 2 )
  {
    mpz_set ( pqueue->u[k], pqueue->u[pq_parent(k)] );
    mpz_set ( pqueue->v[k], pqueue->v[pq_parent(k)] );
    mpz_set ( pqueue->modulus[k], pqueue->modulus[pq_parent(k)] );
    pqueue->w[k] = pqueue->w[pq_parent(k)];
    pqueue->sub_alpha[k] = pqueue->sub_alpha[pq_parent(k)];
  }

  mpz_set (pqueue->u[k], u);
  mpz_set (pqueue->v[k], v);
  mpz_set (pqueue->modulus[k], modulus);
  pqueue->w[k] = w;
  pqueue->sub_alpha[k] = alpha;

  pqueue->used ++;
}


/*
  Sift-down, if the heap is full.
*/
inline void
insert_sub_alpha_pq_down ( sub_alpha_pq *pqueue,
                           int w,
                           mpz_t u,
                           mpz_t v,
                           mpz_t modulus,
                           double alpha )
{
  int k, l;

  for (k = 1; k*2 < pqueue->used; k = l) {

    l = (k << 1);

    /* right > left ? */
    if ( (l+1) < pqueue->used &&
         (pqueue->sub_alpha[l+1] > pqueue->sub_alpha[l]) )
      l ++;

    /* switch larger child with parent */
    if ( pqueue->sub_alpha[l] > alpha ) {
      mpz_set (pqueue->u[k], pqueue->u[l]);
      mpz_set (pqueue->v[k], pqueue->v[l]);
      mpz_set (pqueue->modulus[k], pqueue->modulus[l]);
      pqueue->w[k] = pqueue->w[l];
      pqueue->sub_alpha[k] = pqueue->sub_alpha[l];
    }
    else
      break;
  }

  mpz_set (pqueue->u[k], u);
  mpz_set (pqueue->v[k], v);
  mpz_set (pqueue->modulus[k], modulus);
  pqueue->w[k] = w;
  pqueue->sub_alpha[k] = alpha;
}


/*
  Extract the max of the priority queue.
*/
void
extract_sub_alpha_pq ( sub_alpha_pq *pqueue,
                       int *w,
                       mpz_t u,
                       mpz_t v,
                       mpz_t modulus,
                       double *alpha )
{
  // don't extract u[0] since it is just a placeholder.
  pqueue->used --;
  mpz_set (u, pqueue->u[1]);
  mpz_set (v, pqueue->v[1]);
  mpz_set (modulus, pqueue->modulus[1]);
  *alpha = pqueue->sub_alpha[1];
  *w = pqueue->w[1];

  insert_sub_alpha_pq_down ( pqueue,
                             pqueue->w[pqueue->used],
                             pqueue->u[pqueue->used],
                             pqueue->v[pqueue->used],
                             pqueue->modulus[pqueue->used],
                             pqueue->sub_alpha[pqueue->used] );
}


/*
  Insert to the priority queue.
*/
void
insert_sub_alpha_pq ( sub_alpha_pq *pqueue,
                      int w,
                      mpz_t u,
                      mpz_t v,
                      mpz_t modulus,
                      double alpha )
{

  /* queue is full,  */
  if (pqueue->len == pqueue->used) {
    if ( pqueue->sub_alpha[1] > alpha) {
      insert_sub_alpha_pq_down (pqueue, w, u, v, modulus, alpha);
    }
  }

  /* queue is not full, sift-up */
  else if (pqueue->len > pqueue->used) {
    insert_sub_alpha_pq_up (pqueue, w, u, v, modulus, alpha);
  }
  else {
    fprintf(stderr,"Error: error (pqueue->len < pqueue->used) in insert_sub_alpha_pq()\n");
    exit(1);
  }
}
