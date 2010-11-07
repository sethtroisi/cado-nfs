#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "utils.h"
#include "auxiliary.h"
#include "murphyE.h"
#include "rho.h"
#include "cado.h"

/* At most consider so much */
#define NP 46
static unsigned int primes[NP] = {
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
	31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
	73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
	127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
	179, 181, 191, 193, 197, 199
};

/* Sieve bound for (A + MOD*i)*x + (B + MOD*j) */
typedef struct {
	 long Umax;
	 long Umin;
	 long Vmax;
	 long Vmin;
	 long Amax;
	 long Amin;
	 long Bmax;
	 long Bmin;
	 unsigned long A;
	 unsigned long B;
	 unsigned long MOD;
} _rsbound_t;
typedef _rsbound_t rsbound_t[1];

/* Sieving data struct */
typedef struct {
     mpz_t n;
     mpz_t m;
     mpz_t *f;
     mpz_t *g;
	 mpz_t *fx;
	 mpz_t *gx;
	 mpz_t *numerator;
     double skew;
     int d;
} _rsstr_t;
typedef _rsstr_t rsstr_t[1];

/* Parameters for root sieve (customize) */
typedef struct {
	 /* regarding find_sublattice() functions. */
	 unsigned short len_e_sl;
	 unsigned short *e_sl;
	 unsigned long nbest_sl;
	 unsigned long modulus;

	 /* regarding rootsieve_run() functions. */
	 float sizebound_ratio_rs;
	 unsigned long sizebound_u_rs;
	 mpz_t sizebound_v_rs;
	 unsigned short len_p_rs;
	 float alpha_bound_rs;
} _rsparam_t;
typedef _rsparam_t rsparam_t[1];


/* Struct for the lift. Note we could rely on stack
   in the recursive calls. But this is more convenient
   as long as the memory is not a problem. */
typedef struct node_t {
	 unsigned long u;
	 unsigned long v;
	 unsigned long *r;
	 unsigned short *roottype;
	 unsigned int alloc;
	 unsigned int nr;
	 unsigned int e;
	 double val;
	 struct node_t *firstchild;
	 struct node_t *nextsibling;
	 // might remove parent in future.
	 struct node_t *parent;
} node;


/* List struct for the (u, v) with valuations;
   it is similar to a stack but different in
   permittable operations. */
typedef struct listnode_t {
	 unsigned long u;
	 unsigned long v;
	 double val;
	 struct listnode_t *next;
} listnode;


/* -- ADD FUNCTION DECLARES -- */
static inline unsigned long solve_lineq ( unsigned long a,
										  unsigned long b,
										  unsigned long c,
										  unsigned long p );
