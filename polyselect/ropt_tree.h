#ifndef ROPT_TREE_H
#define ROPT_TREE_H

#include "auxiliary.h"

/* Struct for the lift. Note we could rely on stack
   in the recursive calls. But this is more convenient
   as long as the memory is not a problem. */
typedef struct node_t {
  unsigned int u;
  unsigned int v;
  unsigned int *r;
  char *roottype;
  unsigned int alloc;
  unsigned int nr;
  char e;
  float val;
  struct node_t *firstchild;
  struct node_t *nextsibling;
  // might remove parent in future.
  struct node_t *parent;
} node;

/* Priority queue to record sublattices (u, v) */
typedef struct sublattice_pq_t {
  mpz_t *u;
  mpz_t *v;
  mpz_t *modulus;
  int len;
  int used;
} sublattice_pq;

/* Priority queue to record alpha values */
typedef struct rootscore_pq_t {
  long *i;
  long *j;
  int16_t *alpha;
  int len;
  int used;
} rootscore_pq;

/* Priority queue to record E */
typedef struct MurphyE_pq_t {
  int *w;
  mpz_t *u;
  mpz_t *v;
  double *E;
  int len;
  int used;
} MurphyE_pq;

/* Priority queue to record sublattices (w, u, v)'s alpha's */
typedef struct sub_alpha_pq_t {
  int *w;
  mpz_t *u;
  mpz_t *v;
  mpz_t *modulus;
  double *sub_alpha;
  int len;
  int used;
} sub_alpha_pq;

/* Priority queue for the (u, v) with e */
typedef struct single_sub_alpha_pq_t {
  unsigned int *u;
  unsigned int *v;
  char *e;
  double *val;
  int used;
  int len;
} single_sub_alpha_pq;


/* --- declarations --- */

void new_tree ( node **root );

node* new_node ( void );

void insert_node ( node *parent,
                   node **currnode,
                   unsigned int u,
                   unsigned int v,
                   unsigned int r,
                   char curr_e,
                   unsigned int p,
                   unsigned int pe,
                   char k );

void free_node ( node **ptr );


/* sublattice_pq */
void new_sublattice_pq ( sublattice_pq **ppqueue,
                         unsigned long len );

void insert_sublattice_pq ( sublattice_pq *pqueue,
                            mpz_t u,
                            mpz_t v,
                            mpz_t mod );

void free_sublattice_pq ( sublattice_pq **ppqueue );


/* sub_alpha_pq */
void new_sub_alpha_pq ( sub_alpha_pq **ppqueue,
                        unsigned long len );

void insert_sub_alpha_pq ( sub_alpha_pq *pqueue, 
                           int w,
                           mpz_t u,
                           mpz_t v,
                           mpz_t modulus,
                           double alpha );

void insert_sub_alpha_pq_up ( sub_alpha_pq *pqueue,
                              int w,
                              mpz_t u,
                              mpz_t v,
                              mpz_t modulus,
                              double alpha );

void insert_sub_alpha_pq_down ( sub_alpha_pq *pqueue,
                                int w,
                                mpz_t u,
                                mpz_t v,
                                mpz_t modulus,
                                double alpha );

void extract_sub_alpha_pq ( sub_alpha_pq *pqueue,
                            int *w,
                            mpz_t u,
                            mpz_t v,
                            mpz_t modulus,
                            double *alpha );

void free_sub_alpha_pq ( sub_alpha_pq **ppqueue );


/* MurphyE_pq */
void new_MurphyE_pq ( MurphyE_pq **ppqueue,
                      unsigned long len );

void insert_MurphyE_pq ( MurphyE_pq *pqueue,
                         int w,
                         mpz_t u,
                         mpz_t v,
                         double E );

void insert_MurphyE_pq_up ( MurphyE_pq *pqueue,
                            int w,
                            mpz_t u,
                            mpz_t v,
                            double E );

void insert_MurphyE_pq_down ( MurphyE_pq *pqueue,
                              int w,
                              mpz_t u,
                              mpz_t v,
                              double E );

void free_MurphyE_pq ( MurphyE_pq **ppqueue );


/* single_sub_alpha_pq */
void new_single_sub_alpha_pq ( single_sub_alpha_pq **top,
                               unsigned long len );

void insert_single_sub_alpha_pq ( single_sub_alpha_pq *top,
                                  unsigned int u,
                                  unsigned int v,
                                  double val,
                                  char e );

void insert_single_sub_alpha_pq_up ( single_sub_alpha_pq *top,
                                     unsigned int u,
                                     unsigned int v,
                                     double val,
                                     char e );

void insert_single_sub_alpha_down ( single_sub_alpha_pq *top,
                                    unsigned int u,
                                    unsigned int v,
                                    double val,
                                    char e );

void free_single_sub_alpha_pq ( single_sub_alpha_pq **top );

/* rootscore_pq */
void new_rootscore_pq ( rootscore_pq **ppqueue,
                        unsigned long len );

void insert_rootscore_pq ( rootscore_pq *pqueue,
                           long i,
                           long j,
                           int16_t alpha );

void insert_rootscore_pq_down ( rootscore_pq *pqueue,
                                long i,
                                long j,
                                int16_t alpha );

void insert_rootscore_pq_up ( rootscore_pq *pqueue,
                              long i,
                              long j,
                              int16_t alpha );

void reset_rootscore_pq ( rootscore_pq *pqueue );

void free_rootscore_pq ( rootscore_pq **ppqueue );

int pq_parent ( int i );

#endif /* ROPT_TREE_H */
