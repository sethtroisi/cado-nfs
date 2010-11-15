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

/* 2^i for i <= 150. generated by the asymptotic expansion
   of the first order statistics as pointed by Emmanuel Thom\'e */
static double exp_alpha[] = { 0 ,
							  -1.15365373215 ,-1.52232116258 ,-1.84599014685, -2.13527943823, -2.39830159659, /* 2^2, ...*/
							  -2.64067276005 ,-2.86635204676 ,-3.0782200169 ,-3.27843684076 ,-3.46866583492 ,
							  -3.65021705343 ,-3.8241424823 ,-3.99130105332 ,-4.15240425248 ,-4.30804888816 ,
							  -4.45874113786 ,-4.60491453048 ,-4.74694362117 ,-4.88515454833 ,-5.01983329465 ,
							  -5.15123223117 ,-5.27957535905 ,-5.40506255088 ,-5.52787301451 ,-5.64816814584 ,
							  -5.76609389702 ,-5.88178275645 ,-5.99535541547 ,-6.10692217992 ,-6.21658417274 ,
							  -6.32443436393 ,-6.43055845718 ,-6.53503565678 ,-6.63793933391 ,-6.73933760789 ,
							  -6.83929385532 ,-6.93786715768 ,-7.03511269625 ,-7.13108210169 ,-7.22582376443 ,
							  -7.3193831112 ,-7.41180285201 ,-7.50312320134 ,-7.59338207681 ,-7.68261527801 ,
							  -7.77085664785 ,-7.85813821855 ,-7.94449034388 ,-8.02994181933 ,-8.11451999143 ,
							  -8.19825085747 ,-8.28115915656 ,-8.36326845294 ,-8.44460121241 ,-8.52517887245 ,
							  -8.60502190671 ,-8.68414988441 ,-8.76258152512 ,-8.84033474938 ,-8.9174267255 ,
							  -8.99387391288 ,-9.0696921022 ,-9.14489645276 ,-9.21950152716 ,-9.29352132356 ,
							  -9.36696930575 ,-9.43985843123 ,-9.51220117738 ,-9.58400956592 ,-9.65529518581 ,
							  -9.72606921473 ,-9.79634243911 ,-9.86612527306 ,-9.93542777601 ,-10.0042596694 ,
							  -10.0726303522 ,-10.140548916 ,-10.2080241583 ,-10.2750645962 ,-10.3416784783 ,
							  -10.4078737968 ,-10.4736582979 ,-10.5390394928 ,-10.6040246674 ,-10.6686208913 ,
							  -10.7328350271 ,-10.7966737386 ,-10.8601434986 ,-10.9232505968 ,-10.9860011466 ,
							  -11.0484010921 ,-11.1104562145 ,-11.1721721386 ,-11.233554338 ,-11.2946081414 ,
							  -11.3553387375 ,-11.4157511801 ,-11.4758503931 ,-11.535641175 ,-11.5951282035 ,
							  -11.6543160394 ,-11.713209131 ,-11.7718118176 ,-11.8301283334 ,-11.888162811 ,
							  -11.9459192847 ,-12.0034016939 ,-12.0606138859 ,-12.1175596192 ,-12.1742425661 ,
							  -12.2306663156 ,-12.2868343759 ,-12.342750177 ,-12.3984170731 ,-12.4538383449 ,
							  -12.5090172017 ,-12.5639567839 ,-12.6186601648 ,-12.6731303524 ,-12.7273702918 ,
							  -12.7813828666 ,-12.8351709011 ,-12.8887371614 ,-12.9420843577 ,-12.9952151453 ,
							  -13.0481321268 ,-13.1008378527 ,-13.1533348236 ,-13.2056254911 ,-13.2577122596 ,
							  -13.3095974869 ,-13.3612834859 ,-13.4127725259 ,-13.4640668332 ,-13.5151685928 ,
							  -13.5660799493 ,-13.6168030075 ,-13.6673398341 ,-13.7176924584 ,-13.7678628729 ,
							  -13.8178530347 ,-13.8676648661 ,-13.9173002556 ,-13.9667610588 ,-14.0160490987 ,
							  -14.0651661673 ,-14.1141140255 ,-14.1628944044 ,-14.211509006 /* 2^2, 2^3 ... */ };


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
