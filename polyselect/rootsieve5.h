#include <gmp.h>

/* Things you don't want to change */
#define MAX_LINE_LENGTH 4096
#define PI 3.14159265358979324
#define MAX_DEGREE 6
#define DEBUG 0
#define SUP_ALPHA 4.843
#define NP 46

/* Sieve bound for (A + MOD*i)*x + (B + MOD*j), mainly
   used for stage 2 */
typedef struct {
	 int Wbound;
	 mpz_t Umax;
	 mpz_t Umin;
	 mpz_t Vmax;
	 mpz_t Vmin;
	 long Amax;
	 long Amin;
	 long Bmax;
	 long Bmin;
	 long Blocksize;
	 mpz_t A;
	 mpz_t B;
	 mpz_t MOD;
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
	 double alpha_proj;
     int d;
} _rsstr_t;
typedef _rsstr_t rsstr_t[1];


/* parameters for stage 1 */
typedef struct {

	 /* regarding find_sublattice() functions. */
	 unsigned short len_e_sl;
	 unsigned short tlen_e_sl;
	 unsigned short *e_sl;
	 unsigned short ncrts_sl;
	 unsigned short len_p_rs;
	 unsigned long nbest_sl;

	 /* only those u < the following in sublattices().
		also choose e_sl and len_e_sl such that \prod p_i^{e_i}
		is larger than global_u_bound_rs, then we only need to
		do a line sieving 1x[-U, U]  */
	 unsigned long global_w_bound_rs;
	 unsigned long global_u_bound_rs;
	 /* regarding rootsieve_run() functions. */
	 float sizebound_ratio_rs;
	 float exp_min_alpha_rs;

	 mpz_t global_v_bound_rs;
	 mpz_t modulus;

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


/* Priority queue to record sublattices (u, v) */
typedef struct sublattice_pq_t {
	 mpz_t *u;
	 mpz_t *v;
	 int len;
	 int used;
} sublattice_pq;


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


/* Sieving data struct */
typedef struct {
     mpz_t *f;
     mpz_t *g;
} _bestpoly_t;
typedef _bestpoly_t bestpoly_t[1];


/* Sieving data struct */
typedef struct{
	 int w_left_bound;
	 int w_length;
	 unsigned short s1_num_e_sl;
	 unsigned short *s1_e_sl;

	 int s2_w;
	 mpz_t s2_u;
	 mpz_t s2_v;
	 mpz_t s2_mod;
	 long s2_Amax;
	 long s2_Bmax;

	 /* flag == 0, only use -w and -l
	    flag == 1, use stage 1 params
		flag == 2, use stage 2 params */
	 char flag;
} _param_t;
typedef _param_t param_t[1];

/* -- ADD FUNCTION DECLARES -- */
