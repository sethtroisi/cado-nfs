/*
  Root sieve for degree 5 and degree 6 polynomials.

  [1. Run]
  The format of the input polynomial and output polynomials
  are in CADO format.

  rootsieve5 < POLY -w 0 -l 8 > OUTPUT 2> ERRINFO
  rootsieve5 -f POLY -w 0 -l 8 > OUTPUT 2> ERRINFO

  Note that, "-w" and "-l" will be ignored for deg 5 polynomials.
  For deg 6 polynomial, "-w" is the leftmost of a qudratic rotation
  and "-l" is the step from the "-w".

  [2. Algorithm]

  There are two steps in the root sieve.

  Given a polynomial, the code first tries to find good
  (u, v) (mod (p1^e1 * p2^e2* \cdots *pn^en)) for small prime powers.
  This step is done by hensel-lift-like procedure in a p^2-ary tree data
  structure (for each such p) and discard any bad/impossible (u, v)
  during the tree-building. The multiple roots are updated in a way
  following Emmanuel Thome's idea. Finally, Finally we use CRT to
  find some good (u, v) pairs, those with small u. Note if there are
  too many pi^ei, the CRTs will domimated the running time.

  In the second step, we do the actually root sieve for all primes
  up to some bound. We only consider the (u, v) points lying on the
  sublattice within the sieving range (U, V).

  [3. Degree 6]

  For degree 6 polynomial, the following processes are called in order,

  (a) - For each qudratic rotation w;
  (b) -- Tune parameters;
  (c) -- Find good sublattices (w, u, v)
  (d) - Compare (priority queue) all good sublattices and pick up top ones.
  (e) - For each such sublattice (w, u, v)
  (f) -- Do the root sieve.

  Details:
  (a) In the following command,

  rootsieve5 < POLY -w 0 -l 8 > OUTPUT 2> ERRINFO

  "-w" defines the leftmost point of qudratic rotation.
  f'(x) = f(x) + w*x^2*g(x)

  "-l" defines the steps for quadratic rotation, (w+l-1)

  (b, c)
  The code starts by looking each qudratic rotation, say
  [w, \cdots, w+l-1]. For each rotated polynomial, we will
  have to find a suitable set of parameters (p1^e1 * p2^e2* \cdots *pn^en)
  to produce the sublattice. Note that, if it is larger, then
  we could find polynomials with better alpha, but probably
  worse size; reversly, a small parameter gives worse alpha, but
  probably better size. It is hard to prebuilt a universel parameter.
  Hence  we tune the parameters (p1^e1 * p2^e2* \cdots *pn^en)
  using a trial sieving. The starting point for the parameter
  tunning is in rsparam_setup(). In general, there is no need
  to change this.

  (d)
  After good sublattices (w, u, v) are found for all w in the
  permitted range (as you set by "-w" and "-l"), we will compare
  the alpha values between all these sublattices. At the moment,
  we only pick up the top "rsparam->nbest_sl = 128" sublattices
  for sieving. You may change this in rsparam().

  (e, f)
  For each survived sublattice (w, u, v), we do the root sieve.
  The permitted bound for "v" could be huge which runs forever.
  Therefore, function "rsbound_setup_AB_bound()" limits the sieving
  range for "v".

  Even the range "v" is limited, it may take much memory.
  Therefore, we divide it into "SIEVEARRAY_SIZE = 2097152" in #define.
  In my laptop, it is similar to L2 cache. Setting it larger take more
  memory and can reduce (setup) time.

  For each such "SIEVEARRAY_SIZE", we actually sieving in blocks
  of "L1_SIZE 12288".

  One important parameter:  "TOPALPHA_EACH_SIEVEARRAY 8"

  For each sieve array of "SIEVEARRAY_SIZE", compute the MurphyE of
  8 polynomials whose alpha values are the best among this array.
  These poynomials will be then filtered into another priority queue
  which records "TOPE_EACH_SUBLATTICE=8" polynomials with top MurphyE.

  Note that, for polynomials of small skewness. Size can be more
  important, hence you may want to set "TOPALPHA_EACH_SIEVEARRAY 8"
  larger. However, this may reduce the performance since MurphyE
  computation is slow.

  [4. Tuning]

  The "#define" in the beginning and functions:
  "rsbound_setup_AB_bound()" and  "rsparam_setup ()"
  might be tuned to fit your computation.

  [5. TODO and history]

  -- (Dec) block sieving.
  -- (Dec 30) reduced memory usage in return_all_sublattices().
  -- (Jan 2) addded priority queue, changed precion in return_all_sublattices().
  -- (Feb) some tunnings, correct bugs in rootsieve_v().
  -- (Mar) readme.

  [6. Bugs]

  Please report bugs to Shi Bai (shi.bai AT anu.edu.au).
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "auxiliary.h"
#include "murphyE.h"
#include "rho.h"
#include "rootsieve5.h"

/* Things for tuning */
#define L1_SIZE 12288 // ~ l1 cache
#define SIEVEARRAY_SIZE 20971520 // ~ l2 cache = 2097152
#define TUNE_SIEVEARRAY_SIZE L1_SIZE / 2
#define TOPALPHA_EACH_SIEVEARRAY 8 // For each "SIEVEARRAY_SIZE", record top 8 poly's alpha avalues.
#define TOPE_EACH_SUBLATTICE 8 // For sublattice, record top 8 poly's E avalues.
#define LEN_SUBLATTICE_PRIMES 9

/* Define primes, exp_alpha. */
unsigned int primes[NP] = {
	 2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
	 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
	 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
	 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
	 179, 181, 191, 193, 197, 199
};

/* next_prime_idx[3] = ind(5) = 2 in above table */
unsigned char next_prime_idx[] = {
	 0, /* 0 */
	 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, /*1 - 10*/
	 5, 5, 6, 6, 6, 6, 7, 7, 8, 8,
	 8, 8, 9, 9, 9, 9, 9, 9, 10, 10,
	 11, 11, 11, 11, 11, 11, 12, 12, 12, 12,
	 13, 13, 14, 14, 14, 14, 15, 15, 15, 15,
	 15, 15, 16, 16, 16, 16, 16, 16, 17, 17,
	 18, 18, 18, 18, 18, 18, 19, 19, 19, 19,
	 20, 20, 21, 21, 21, 21, 21, 21, 22, 22,
	 22, 22, 23, 23, 23, 23, 23, 23, 24, 24,
	 24, 24, 24, 24, 24, 24, 25, 25, 25, 25,
	 26, 26, 27, 27, 27, 27, 28, 28, 29, 29,
	 29, 29, 30, 30, 30, 30, 30, 30, 30, 30,
	 30, 30, 30, 30, 30, 30, 31, 31, 31, 31,
	 32, 32, 32, 32, 32, 32, 33, 33, 34, 34,
	 34, 34, 34, 34, 34, 34, 34, 34, 35, 35,
	 36, 36, 36, 36, 36, 36, 37, 37, 37, 37,
	 37, 37, 38, 38, 38, 38, 39, 39, 39, 39,
	 39, 39, 40, 40, 40, 40, 40, 40, 41, 41,
	 42, 42, 42, 42, 42, 42, 42, 42, 42, 42,
	 43, 43, 44, 44, 44, 44, 45, 45 /* 191 - 198 */
};


/* 2^i for i <= 150. generated by the asymptotic expansion
   of the first order statistics as pointed by Emmanuel Thom\'e */
double exp_alpha[] = {
	 0 ,
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
	 -14.0651661673 ,-14.1141140255 ,-14.1628944044 ,-14.211509006 };

/* local declare */
static inline unsigned long solve_lineq ( unsigned long a,
										  unsigned long b,
										  unsigned long c,
										  unsigned long p );
static double get_biased_alpha_projective ( mpz_t *f,
											const int d,
											unsigned long B );

/*-----------------------------*/
/*   @Input-output functions.  */
/*-----------------------------*/

/*
  Some information-related routines
*/
static int
readline ( char *s )
{
	 int c;

	 while ((c = getchar ()) != '\n' && c != EOF)
		  *s++ = c;
	 *s = '\0';
	 fprintf (stderr, "%s", s);
	 return c;
}

static void
read_ggnfs ( mpz_t N,
			 mpz_t *f,
			 mpz_t *g,
			 mpz_t M )
{
	 fflush(stdin);
	 int i, ret, d;
	 char s[MAX_LINE_LENGTH]; /* input buffer */

	 while (feof (stdin) == 0) {
		  ret = readline (s);
		  if (ret == EOF)
			   break;
		  if (strlen (s) + 1 >= MAX_LINE_LENGTH) {
			   fprintf (stderr, "Error, buffer overflow\n");
			   exit (1);
		  }
		  if (strncmp (s, "n:", 2) == 0) { /* n: input number */
			   if (mpz_set_str (N, s + 2, 0) != 0) {
					fprintf (stderr, "Error while reading N\n`");
					exit (1);
			   }
		  }
		  else if (sscanf (s, "c%d:", &i) == 1) {/* ci: coeff of degree i */
			   if (i > MAX_DEGREE) {
					fprintf (stderr, "Error, too large degree %d\n", i);
					exit (1);
			   }
			   if (mpz_set_str (f[i], s + 3, 0) != 0) {
					fprintf (stderr, "Error while reading f[%d]", i);
					exit (1);
			   }
		  }
		  else if (strncmp (s, "Y1:", 3) == 0) {
			   if (mpz_set_str (g[1], s + 3, 0) != 0) {
					fprintf (stderr, "Error while reading Y1");
					exit (1);
			   }
		  }
		  else if (strncmp (s, "Y0:", 3) == 0) {
			   if (mpz_set_str (g[0], s + 3, 0) != 0)
			   {
					fprintf (stderr, "Error while reading Y0");
					exit (1);
			   }
		  }
		  else if (strncmp (s, "# M", 3) == 0 || strncmp (s, "m:", 2) == 0) {
			   if (mpz_set_str (M, s + 2 + (s[0] == '#'), 0) != 0) {
					fprintf (stderr, "Error while reading M or m");
					exit (1);
			   }
		  }
	 }

	 for (d = MAX_DEGREE; d > 0 && mpz_cmp_ui (f[d], 0) == 0; d --);
	 if (mpz_cmp_ui (M, 0) == 0) {
		  mpz_t t;
		  /* M = -Y0/Y1 mod N */
		  mpz_invert (M, g[1], N);
		  mpz_neg (M, M);
		  mpz_mul (M, M, g[0]);
		  mpz_mod (M, M, N);
		  /* check M is also a root of the algebraic polynomial mod N */
		  mpz_init_set (t, f[d]);
		  for (i = d - 1; i >= 0; i --) {
			   mpz_mul (t, t, M);
			   mpz_add (t, t, f[i]);
			   mpz_mod (t, t, N);
		  }
		  if (mpz_cmp_ui (t, 0) != 0) {
			   fprintf (stderr, "Polynomials have no common root mod N\n");
			   exit (1);
		  }
		  mpz_clear (t);
	 }
}


/* Care: poly_print, print_poly are defined elsewhere. */
static void
print_poly2 ( mpz_t *f,
			  mpz_t *g,
			  int deg,
			  mpz_t N )
{
	 int i;
	 gmp_printf ("\nn: %Zd\n", N);
	 for (i = deg; i >= 0; i --)
	 {
		  gmp_printf ("c%d: %Zd\n", i, f[i]);
	 }
	 for (i = 1; i >= 0; i --)
	 {
		  gmp_printf ("Y%d: %Zd\n", i, g[i]);
	 }
}


/*
  Print polynomial and info: lognorm, skew, alpha.
*/
static double
print_poly_info ( mpz_t *f,
				  mpz_t *g,
				  int d,
				  mpz_t N,
				  mpz_t M,
				  int verbose )
{
	 /* print info about the polynomial */
	 unsigned int nroots = 0;
	 double skew, logmu, alpha, e, alpha_proj;
	 int i;

	 /* initlize cado_poly for Murphy E */
	 cado_poly cpoly;
	 cado_poly_init(cpoly);

	 for (i = 0; i < (d + 1); i++) {
		  mpz_set(cpoly->f[i], f[i]);
	 }
	 for (i = 0; i < 2; i++) {
		  mpz_set(cpoly->g[i], g[i]);
	 }

	 if (verbose == 2) {

		  /* output original poly */
		  gmp_printf ("\nn: %Zd\n", N);
		  for (i = d; i >= 0; i --) {
			   gmp_printf ("c%d: %Zd\n", i, f[i]);
		  }
		  for (i = 1; i >= 0; i --) {
			   gmp_printf ("Y%d: %Zd\n", i, g[i]);
		  }
		  if (verbose == 3) // don't want m in general, and this m might be wrong
		  	   gmp_printf ("m: %Zd\n", M);
	 }

	 /* compute skew, logmu, nroots */
	 nroots = numberOfRealRoots (f, d, 0, 0);
	 skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
	 logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);
	 alpha = get_alpha (f, d, ALPHA_BOUND);
	 alpha_proj = get_biased_alpha_projective (f, d, ALPHA_BOUND);

	 mpz_set (cpoly->n, N);
	 cpoly->degree = d;
	 cpoly->degreeg = 2;
	 cpoly->skew = skew;
	 e = MurphyE (cpoly, BOUND_F, BOUND_G, AREA, MURPHY_K);

	 if (verbose == 2) {
		  printf ("# skew: %.2f, ", skew);
		  printf ("lognorm: %.2f, alpha: %.2f, (alpha_proj: %.2f) E: %.2f, nr: %u \n# MurphyE: %1.2e (Bf=%.0f, Bg=%.0f, area=%1.2e)\n",
				  logmu,
				  alpha,
				  alpha_proj,
				  logmu + alpha,
				  nroots,
				  e,
				  BOUND_F,
				  BOUND_G,
				  AREA );
	 }

     fflush( stdout );
	 cado_poly_clear (cpoly);

	 return e;
}


/*-----------------------------*/
/*   @Tree-related functions.  */
/*-----------------------------*/


#if 0
/*
  Print the info for the node
*/
static void
print_node ( node *pnode )
{
	 printf("(%lu,%lu):%d:%d:(%.2f)\n",pnode->u, pnode->v, pnode->nr, pnode->e, pnode->val);
	 /* int i; */
	 /* for (i = 0; i < pnode->nr; i++) */
	 /*  	  printf ("pnode->r[%d]: %lu\n", i, pnode->r[i]); */
}


/*
  Print a tree, non-recursive. Two styles.
  - If level == 0, print the whole tree;
  - If level == i, print nodes at height i.
*/
static void
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
static node *
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
static void
new_tree ( node **root )
{
	 *root = NULL;
}


/*
  Free current node ptr.
*/
static void
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
static void
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
  Double the memory space of pnode->r[]. If need_roottype
  is 1, we also want to allocate memory for pnode->roottype[].
  This is only used for base case in the find_sublattice().

  TBA: This function is buggy. You should always call
  alloc_r_node(ptr, 1) and then alloc_r_node(ptr, 0) without
  any overlapping. In our case, we call alloc_r_node(ptr, 1)
  for all the base case in find_sublattice(), and call
  alloc_r_node(ptr, 0) for all the lifted cases.
*/
static void
alloc_r_node ( node *pnode,
			   int need_roottype )
{
	 //assert (pnode->alloc <= pnode->nr);

	 if (pnode->nr == 0) {
		  pnode->r = (unsigned long *) realloc (pnode->r, (sizeof (unsigned long)));
		  if (need_roottype == 1)
			   pnode->roottype = (unsigned short *) realloc (pnode->roottype, (sizeof (unsigned short)));
		  pnode->alloc = 1;
	 }
	 else {
		  pnode->r = (unsigned long *) realloc (pnode->r, (unsigned long) (2 * (pnode->nr)) * sizeof (unsigned long));
		  if (need_roottype == 1)
			   pnode->roottype = (unsigned short *) realloc (pnode->roottype, (unsigned long) (2 * (pnode->nr)) * sizeof (unsigned short));
		  pnode->alloc = (unsigned long) (pnode->alloc * 2);
	 }

	 if (pnode->r == NULL) {
		  fprintf (stderr, "Error, cannot reallocate memory in alloc_r_node\n");
		  exit (1);
	 }
	 if ( (need_roottype == 1) && (pnode->roottype == NULL) ) {
		  fprintf (stderr, "Error, cannot reallocate memory in alloc_r_node\n");
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

  k = 0 means don't allocate for currnode->roottype[]
  k = 1/2 means currnode->roottype[] = 1 or 2
*/
static void
insert_node ( node *parent,
			  node **currnode,
			  unsigned long u,
			  unsigned long v,
			  unsigned long r,
			  unsigned int curr_e,
			  unsigned long pe,
			  int k )
{
	 unsigned int i, r_exit = 0;
	 node *lastnode = NULL;
	 node *nextnode = NULL;
	 nextnode = parent->firstchild;
	 lastnode = nextnode;

	 /* if the (u, v) pair already exists, don't
		create new node, but add possibly new root. */
	 while ((nextnode != NULL)) {

		  if ((nextnode->u == u) && (nextnode->v == v)) {
			   /* if (u,v) pair exists, see whether
				  the root r already exists in the node. */
			   for (i = 0; i < nextnode->nr; i ++) {
					if ((nextnode->r[i]) == r)
						 r_exit = 1;
			   }
			   /* r already exits, do nothing and return */
			   if (r_exit == 1)
					return;

			   /* otherwise, insert this root */
			   if (nextnode->alloc <= nextnode->nr) {
					/* rellocate memory for the new root without allocating nextnode->roottype[]. */
					if (k == 0)
						 alloc_r_node (nextnode, 0);
					/* rellocate memory for the new root and allocate nextnode->roottype[]. */
					else
						 alloc_r_node (nextnode, 1);
					//printf ("alloc: %lu, nr: %lu\n", nextnode->alloc, nextnode->nr);
			   }

			   if (k != 0)
					nextnode->roottype[nextnode->nr] = k;
			   nextnode->r[nextnode->nr] = r;
			   nextnode->nr += 1;
			   /* such (u, v) exists, hence we increase the valuation. */
			   nextnode->val += 1.0 / (double) pe;
			   /* printf ("(u: %lu, v: %lu), nr: %lu, r: %lu, e: %lu, val: %f\n", u, v,
				  nextnode->nr, nextnode->r[nextnode->nr], curr_e, nextnode->val); */
			   break;
		  }
		  lastnode = nextnode;
		  nextnode = nextnode->nextsibling;
	 }

	 /* add new node */
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
		  /* rellocate memory for the new root without allocating nextnode->roottype[]. */

		  if (k == 0)
			   alloc_r_node (*currnode, 0);
		  /* rellocate memory for the new root and allocate nextnode->roottype[]. */
		  else {
			   alloc_r_node (*currnode, 1);
			   (*currnode)->roottype[0] = k;
		  }
		  (*currnode)->r[0] = r;
		  (*currnode)->nr += 1;
		  (*currnode)->e = curr_e;
		  /* such (u, v) is new, hence we inheritate the val from its parent. */
		  (*currnode)->val = 1.0 / (double) pe + (*currnode)->parent->val;
		  //printf ("(u: %lu, v: %lu), e: %lu, val: %f\n", u, v, curr_e, (*currnode)->val);
	 }

	 //print_tree(parent, 0);
}


/*-----------------------------*/
/* @Listnode (u, v) valuations */
/*-----------------------------*/


/*
  Initialise new list for (u, v) p-valuations.
*/
static void
new_list ( listnode **top )
{
	 *top = NULL;
}


/*
  free non-empty listnode.
*/
static void
free_listnode ( listnode **pplistnode )
{
	 if (*pplistnode)
		  free (*pplistnode);
}


/*
  Del the list pointed by top.
*/
static void
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
static void
print_listnode ( listnode *plistnode )
{
	 printf("(u, v): (%lu,%lu), val: %.2f\n", plistnode->u, plistnode->v, plistnode->val);
}


/*
  Print the list pointed by top.
*/
static void
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
static unsigned long
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
static listnode *
new_listnode ( unsigned long u,
			   unsigned long v,
			   double val )
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
	 plistnode->next = NULL;
	 return plistnode;
}


/*
  Insert a listnode to current list <- top, only record best one.
*/
static void
insert_listnode ( listnode **top,
				  unsigned long u,
				  unsigned long v,
				  double val )
{
	 listnode *newlistnode;

	 /* if empty, create node */
	 if ( (*top) == NULL ) {
		  newlistnode = new_listnode (u, v, val);
		  (*top) = newlistnode;
		  (*top)->next = NULL;
	 }
	 else {
		  /* if income has better val, delete and add */
		  if ( (*top)->val < val) {
			   newlistnode = new_listnode (u, v, val);
			   free_list (top);
			   (*top) = newlistnode;
			   (*top)->next = NULL;
		  }
		  /* if income has equal val, add */
		  else if ( (*top)->val == val) {
			   newlistnode = new_listnode (u, v, val);
			   newlistnode->next = (*top);
			   (*top) = newlistnode;
		  }
	 }
	 newlistnode = NULL;
}


/*
  Insert a listnode to current list, record all. Line 2323
*/
static void
insert_listnode_plain ( listnode **top,
						unsigned long u,
						unsigned long v,
						double val )
{
	 listnode *newlistnode;

	 /* if empty, create node */
	 if ( (*top) == NULL ) {
		  newlistnode = new_listnode (u, v, val);
		  (*top) = newlistnode;
		  (*top)->next = NULL;
	 }
	 else {
		  /* if income has equal val, add */
		  newlistnode = new_listnode (u, v, val);
		  newlistnode->next = (*top);
		  (*top) = newlistnode;
	 }
	 newlistnode = NULL;
}

#if 0
/*
  Scan a tree, return the best (u, v).
  //scan_tree (root, e[i], &top);
  */
static void
scan_tree ( node *root,
			int level,
			listnode **top )
{

	 if (level <= 0) {
		  fprintf(stderr,"Error: level > 0. \n");
		  exit(1);
	 }

	 int currlevel = 0;
	 node *ptr = root;

	 if (!ptr)  /* if empty */
		  return;

     /* find leftmost */
	 while (ptr->firstchild) {
		  ptr = ptr->firstchild;
		  ++ currlevel;
	 }
	 while (ptr) {
		  if (currlevel == level) {
			   insert_listnode (top, ptr->u, ptr->v, ptr->val);
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
	 return;
}
#endif


/*
  Some indices, note the queue is shifted by 1.
*/
static inline int
pq_parent ( int i )
{
	 return (i>>1);
}

static inline int
pq_leftchild ( int i )
{
	 return (i<<1);
}

static inline int
pq_rightchild ( int i )
{
	 return ( (i << 1) + 1 );
}


/*
  Create priority queue for best sublattices of length len.
*/
static void
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
	 if ( (*ppqueue)->u == NULL || (*ppqueue)->v == NULL ) {
		  fprintf(stderr,"Error: malloc failed in new_sublattice_pq()\n");
		  exit(1);
	 }

	 int i;
	 for (i = 0; i < (*ppqueue)->len; i++)
	 {
		  mpz_init ( (*ppqueue)->u[i] );
		  mpz_init ( (*ppqueue)->v[i] );
	 }

	 mpz_set_str ( (*ppqueue)->u[0], "340282366920938463463374607431768211456", 10 ); // 2^128
	 mpz_set_ui ( (*ppqueue)->v[0], 0 );

	 (*ppqueue)->used = 1UL; // u[0] and v[0] are null elements

}


/*
  Create priority queue for best sublattices of length len.
*/
static void
free_sublattice_pq ( sublattice_pq **ppqueue )
{
	 int i;
	 for (i = 0; i < (*ppqueue)->len; i++)
	 {
		  mpz_clear ( (*ppqueue)->u[i] );
		  mpz_clear ( (*ppqueue)->v[i] );
	 }

	 free ( (*ppqueue)->u );
	 free ( (*ppqueue)->v );
	 free ( *ppqueue );
}


/*
  Sift-up to add, if the queue is not full.
*/
static inline void
insert_sublattice_pq_up ( sublattice_pq *pqueue,
						  mpz_t u,
						  mpz_t v )
{
	 int i;

	 for ( i = pqueue->used;
		   mpz_cmpabs (u, pqueue->u[pq_parent(i)]) > 0; i /= 2 ) {

		  mpz_set ( pqueue->u[i], pqueue->u[pq_parent(i)] );
		  mpz_set ( pqueue->v[i], pqueue->v[pq_parent(i)] );
	 }

	 mpz_set (pqueue->u[i], u);
	 mpz_set (pqueue->v[i], v);
	 pqueue->used ++;
}


/*
  Sift-down, if the heap is full.
*/
static inline void
insert_sublattice_pq_down ( sublattice_pq *pqueue,
							mpz_t u,
							mpz_t v )
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
		  }
		  else
			   break;
	 }

	 mpz_set (pqueue->u[i], u);
	 mpz_set (pqueue->v[i], v);
}

#if 0
/*
  Extract the max of the priority queue.
*/
static void
extract_sublattice_pq ( sublattice_pq *pqueue,
						mpz_t u,
						mpz_t v )
{
	 // don't extract u[0] since it is just a placeholder.
	 pqueue->used --;
	 mpz_set (u, pqueue->u[1]);
	 mpz_set (v, pqueue->v[1]);

	 insert_sublattice_pq_down ( pqueue,
								 pqueue->u[pqueue->used],
								 pqueue->v[pqueue->used] );
}
#endif

/*
  Insert to the priority queue.
*/
static void
insert_sublattice_pq ( sublattice_pq *pqueue,
					   mpz_t u,
					   mpz_t v )
{
	 //gmp_fprintf (stderr, "# Debug: inserting (%Zd, %Zd), used: %d, len: %d\n", u, v, pqueue->used, pqueue->len);

	 /* queue is full,  */
	 if (pqueue->len == pqueue->used) {
		  if ( mpz_cmpabs (pqueue->u[1], u) > 0 ) {
			   insert_sublattice_pq_down (pqueue, u, v);
		  }
	 }

	 /* queue is not full, sift-up */
	 else if (pqueue->len > pqueue->used) {
		  insert_sublattice_pq_up (pqueue, u, v);
	 }
	 else {
		  fprintf(stderr,"Error: error (pqueue->len < pqueue->used) in insert_sublattice_pq()\n");
		  exit(1);
	 }
}


/*
  Create priority queue for best root scores (among MAT[] for a sublattice).
*/
static void
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
	 (*ppqueue)->used = 1UL;

}


/*
  Reset priority queue for best root scores (among MAT[] for a sublattice).
*/
static void
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
	 pqueue->used = 1UL;
}


/*
  Free
*/
static void
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
static inline void
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
static inline void
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
static void
insert_rootscore_pq ( rootscore_pq *pqueue,
					  long i,
					  long j,
					  int16_t alpha )
{
	 /* fprintf (stderr, "# Debug: inserting (%ld, %ld), alpha: %d, used: %d, len: %d\n", */
	 /* 		  i, j, alpha, pqueue->used, pqueue->len); */

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
static void
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

	 (*ppqueue)->used = 1UL;
}


/*
  Free
*/
static void
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
static inline void
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
static inline void
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
static void
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
static void
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

	 (*ppqueue)->used = 1UL;
}


/*
  free
*/
static void
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
static inline void
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
static inline void
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
static void
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
static void
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


static void
crt_pair_mp ( mpz_t a,
			  mpz_t p1,
			  mpz_t b,
			  mpz_t p2,
			  mpz_t re )
{
	 mpz_t tmp;
	 mpz_init (tmp);
	 /* 1/p1 (mod p2) */
	 mpz_invert (tmp, p1, p2);
	 mpz_sub (re, b, a);
	 mpz_mul (re, re, tmp);
	 mpz_mod (re, re, p2);
	 mpz_mul (re, re, p1);
	 mpz_add (re, re, a);
	 mpz_clear (tmp);
}


#if 0


/*
  Recover x (mod p1*p2) in
  x = a (mod p1)
  x = b (mod p2)

  - p1 is probably larger than p2.
*/
static unsigned long
crt_pair ( unsigned long a,
		   unsigned long p1,
		   unsigned long b,
		   unsigned long p2 )
{
	 unsigned long tmp;
	 /* solve x in a + p1*x = b (mod p2) */
	 tmp = solve_lineq (a, p1, b, p2);
	 tmp = tmp * p1 + a;
	 return tmp;
}


/*
  Do crt for each element between lists <- top1 and top2
  Return the new list pointed to top3.
*/
static void
crt_list ( listnode *top1,
		   unsigned long pe1,
		   listnode *top2,
		   unsigned long pe2,
		   listnode **top3 )
{
	 if ( (!top1) || (!top2) ) {
		  fprintf (stderr, "Error, empty list in crt_list().\n");
		  exit (1);
	 }

	 unsigned long tmp1, tmp2;
	 listnode *holdtop2 = top2, *holdtop3 = NULL;

	 do {
		  top2 = holdtop2;
		  do {
			   tmp1 = crt_pair(top1->u, pe1, top2->u, pe2);
			   tmp2 = crt_pair(top1->v, pe1, top2->v, pe2);

			   /* printf ("(%lu, %lu) in %lu\t", top1->u, top1->v, pe1);
				  printf ("(%lu, %lu) in %lu\n", top2->u, top2->v, pe2); */

			   /* add to list <- (*top3) */
			   holdtop3 = (*top3);
			   (*top3) = new_listnode (tmp1, tmp2, 0.0);
			   (*top3)->next = holdtop3;

			   //printf ("Sol (u, v): %lu, %lu\n", tmp1, tmp2);
			   top2 = top2->next;
		  }
		  while (top2);
		  top1 = top1->next;
	 }
	 while (top1);
}
#endif

/*-----------------------------*/
/*   @Some arithmetics.        */
/*-----------------------------*/


/*
  Solve x in a + b*x = c (mod p)
*/
static inline unsigned long
solve_lineq ( unsigned long a,
			  unsigned long b,
			  unsigned long c,
			  unsigned long p )
{
	 /* in general, we should know that gcd(b, p) = 1 */
	 if (b % p == 0) {
		  fprintf (stderr, "Error, impossible inverse in solve_lineq().\n");
		  exit (1);
	 }

	 unsigned long tmp;
	 modulusul_t mod;
	 residueul_t tmpr, ar, cr;
	 modul_initmod_ul (mod, p);
	 modul_init (tmpr, mod);
	 modul_init (cr, mod);
	 modul_init (ar, mod);
	 modul_set_ul (cr, c, mod);
	 modul_set_ul (ar, a, mod);
	 modul_sub (cr, cr, ar, mod);
	 modul_set_ul (tmpr, b, mod);
	 modul_inv (tmpr, tmpr, mod);
	 modul_mul (tmpr, tmpr, cr, mod);
	 tmp = modul_get_ul(tmpr, mod);
	 modul_clear (tmpr, mod);
	 modul_clear (cr, mod);
	 modul_clear (ar, mod);
	 modul_clearmod (mod);
	 return tmp;
}

#if 0
/*
  Could call from auxiliary.c.
*/
static void
content_poly ( mpz_t g,
			   mpz_t *f,
			   int d )
{
	 int i;

	 ASSERT(d >= 1);

	 mpz_gcd (g, f[0], f[1]);
	 for (i = 2; i <= d; i++)
		  mpz_gcd (g, g, f[i]);
}
#endif

/*
  Affine part of the special valution for polynomial f over p.
*/
static double
special_valuation_affine ( mpz_t * f,
						   int d,
						   unsigned long p,
						   mpz_t disc )
{
	 double v;
	 int pvaluation_disc = 0;
	 double pd = (double) p;

	 if (mpz_divisible_ui_p(disc, p)) {
		  mpz_t t;
		  pvaluation_disc++;
		  mpz_init(t);
		  mpz_divexact_ui(t, disc, p);
		  if (mpz_divisible_ui_p(t, p))
			   pvaluation_disc++;
		  mpz_clear(t);
	 }

	 if (pvaluation_disc == 0) {
		  /* case 1: root must be simple*/
		  int e = 0;
		  e = poly_roots_ulong(NULL, f, d, p);

		  return (pd * e) / (pd * pd - 1);
	 }
	 /* else if (pvaluation_disc == 1) { */
	 /* 	  /\* case 2: special case where p^2 does not divide disc *\/ */
	 /* 	  int e = 0; */
	 /* 	  e = poly_roots_ulong(NULL, f, d, p); */

	 /* 	  /\* something special here. *\/ */
	 /* 	  return (pd * e - 1) / (pd * pd - 1); */

	 /* } */
	 else {
		  v = special_val0(f, d, p) * pd;
		  v /= pd + 1.0;
		  return v;
	 }
}


/*
  Find biased alpha_projective for a poly f. It uses
  some hacks here which need to be changed in future.
  Until now, since this will only be done several
  times, hence the speed is not critical.

  Note that, the returned alpha is the  -val * log(p)
  biased part in the alpha. Hence, we can just add
  this to our affine part.
*/
double
get_biased_alpha_projective ( mpz_t *f,
							  const int d,
							  unsigned long B )
{
	 double alpha, e;
	 unsigned long p;
	 mpz_t disc;

	 mpz_init (disc);
	 discriminant (disc, f, d);

	 /* prime p=2 */
	 e = special_valuation (f, d, 2, disc) - special_valuation_affine (f, d, 2, disc);

	 /* 1/(p-1) is counted in the affine part */
	 alpha =  (- e) * log (2.0);

	 /* FIXME: generate all primes up to B and pass them to get_alpha */
	 for (p = 3; p <= B; p += 2)
		  if (isprime (p)) {
			   e = special_valuation(f, d, p, disc) - special_valuation_affine (f, d, p, disc);
			   alpha += (- e) * log ((double) p);
		  }

	 mpz_clear (disc);

	 return alpha;
}

#if 0
/*
  Similar to above, but for affine part.
*/
static double
get_biased_alpha_affine ( mpz_t *f,
						  const int d,
						  unsigned long B )
{
	 double alpha, e;
	 unsigned long p;
	 mpz_t disc;

	 mpz_init (disc);
	 discriminant (disc, f, d);

	 /* prime p=2 */
	 e = special_valuation_affine (f, d, 2, disc);
	 alpha =  (1.0 - e) * log (2.0);

	 //printf ("\np: %u, val: %f, alpha: %f\n", 2, e, alpha);

	 /* FIXME: generate all primes up to B and pass them to get_alpha */
	 for (p = 3; p <= B; p += 2)
		  if (isprime (p)) {
			   e = special_valuation_affine (f, d, p, disc);
			   alpha += (1.0 / (double) (p - 1) - e) * log ((double) p);
			   //printf ("\np: %u, val: %f, alpha: %f\n", p, e, alpha);

		  }
	 mpz_clear (disc);
	 return alpha;
}


/*
  Contribution from a particular multiple root r of the polynomial f
  over p. Note, r must also be a double root of f mod p.
*/
static double
average_valuation_affine_root ( mpz_t *f,
								int d,
								unsigned long p,
								unsigned long r )
{
	 unsigned long v = 0UL;
	 int i, j;
	 mpz_t c, *fv;
	 double val;

	 mpz_init (c);

	 /* init fv */
	 fv = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
	 if (fv == NULL) {
		  fprintf (stderr, "Error, cannot allocate memory in average_valuation_affine_root.\n");
		  exit (1);
	 }

	 for (i = 0; i <= d; i++)
		  mpz_init_set (fv[i], f[i]);

	 /* remove the p-valuations from fv */
	 content_poly (c, f, d);
	 while (mpz_divisible_ui_p(c, p)) {
		  v += 1;
		  for (i = 0; i <= d; i ++) {
			   mpz_fdiv_q_ui (fv[i], f[i], p);
		  }
	 }

	 /* first translate, then scale */
	 for (i = d - 1; i >= 0; i--)
		  for (j = i; j < d; j++)
			   mpz_addmul_ui (fv[j], fv[j+1], r);
	 /* t is p^i */
	 mpz_set_ui(c, 1);
	 for (i = 0; i <= d; i++) {
		  mpz_mul(fv[i], fv[i], c);
		  mpz_mul_ui(c, c, p);
	 }

	 /* now c is disc. */
	 discriminant (c, fv, d);
	 val = special_valuation_affine (fv, d, p, c);
	 val = val / (double) p;

	 /* clear */
	 for (i = 0; i <= d; i++) {
		  mpz_clear (fv[i]);
	 }

	 /* !!! REMEMBER THIS !!! */
	 free (fv);
	 mpz_clear(c);
	 return val;
}
#endif

/*
  Change coordinate from (u, v) to (a, b),
  where A + MOD*a = u.
*/
static inline long
uv2ab ( long A,
		long MOD,
		long u )
{
	 return ( (u - A) / MOD);
}


/*
  Find coordinate a such that
  A + MOD*a = u (mod p).
*/
static inline unsigned long
uv2ab_mod ( mpz_t A,
			mpz_t MOD,
			unsigned long U,
			unsigned long p )
{
	 unsigned long a = mpz_fdiv_ui (A, p);
	 unsigned long mod = mpz_fdiv_ui (MOD, p);
	 unsigned long u = U % p;

	 /* compute the A + MOD * a = u (mod p) */
	 return solve_lineq(a, mod, u, p);
}


/*
  Change coordinate from (a, b) to (u, v),
  where A + MOD*a = u.
*/
static inline void
ab2uv ( mpz_t A,
		mpz_t MOD,
		long a,
		mpz_t u )
{
	 mpz_mul_si (u, MOD, a);
	 mpz_add (u, u, A);
}


/*
  Change coordinate from (a, b) to the index of
  the sieving array, where index = a - Amin,
  where Amin is negative.
*/
static inline long
ab2ij ( long Amin,
		long a )
{
	 return ( a - Amin );
}


/*
  Change coordinate from (i, j) to (a, b).
*/
static inline long
ij2ab ( long Amin,
		long i )
{
	 return ( i + Amin );
}


/*
  Change coordinate from (i, j) to (u, v).
*/
static inline void
ij2uv ( mpz_t A,
		mpz_t MOD,
		long Amin,
		long i,
		mpz_t u )
{
	 ab2uv(A, MOD, ij2ab(Amin, i), u);
}


#if DEBUG
/*
  Change coordinate from (i, j) to (u, v).
*/
static void
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
  Given numerator = f(x)*g'(x) - f'(x)*g(x) and g(x),
  compute special_u = numerator/g^2(x) (mod p).
*/
static inline unsigned long
compute_special_u ( unsigned long numerator,
					unsigned long gx,
					unsigned long p )
{
	 modulusul_t mod;
	 residueul_t tmp, tmp1;
	 unsigned long special_u;
	 modul_initmod_ul (mod, p);
	 modul_init (tmp, mod);

	 modul_set_ul (tmp1, numerator, mod);
	 gx = gx * gx; /* This should be ok since gx < p < 2000 */
	 modul_set_ul (tmp, gx, mod);

	 modul_inv (tmp, tmp, mod);
	 modul_mul (tmp, tmp, tmp1, mod);
	 special_u = modul_get_ul (tmp, mod);

	 modul_clear (tmp, mod);
	 modul_clear (tmp1, mod);
	 modul_clearmod (mod);

	 return special_u;
}


/*
  Similar to above, but test whether r is a multiple root.
  Given numerator = f(x)*g'(x) - f'(x)*g(x) and g(x),
  test  u*g^2(x) == numerator (mod p).
*/
static inline unsigned long
test_special_u ( unsigned long numerator,
				 unsigned long gx,
				 unsigned long u,
				 unsigned long p )
{
	 modulusul_t mod;
	 residueul_t tmp, tmp1;
	 unsigned long re;
	 modul_initmod_ul (mod, p);
	 modul_init (tmp, mod);
	 modul_init (tmp1, mod);

	 gx = gx * gx; /* This should be ok since gx < p < 2000 */
	 modul_set_ul (tmp, gx, mod);
	 modul_set_ul (tmp1, u, mod);
	 modul_mul (tmp, tmp, tmp1, mod);
	 re = (modul_get_ul (tmp, mod) == numerator);

	 modul_clear (tmp, mod);
	 modul_clear (tmp1, mod);
	 modul_clearmod (mod);
	 return re;
}


/*
  Compute fuv = f+(u*x+v)*g,
  f(r) + u*r*g(r) + v*g(r) = 0
  The inputs for f and g are mpz.
*/
static inline void
compute_fuv_mp ( mpz_t *fuv,
				 mpz_t *f,
				 mpz_t *g,
				 int d,
				 mpz_t u,
				 mpz_t v )
{
	 mpz_t tmp, tmp1;
	 mpz_init (tmp);
	 mpz_init (tmp1);
	 int i = 0;

	 for (i = 3; i <= d; i ++)
		  mpz_set (fuv[i], f[i]);

	 /* f + u*g1*x^2
		+ (g0*u* + v*g1)*x
		+ v*g0 */

	 /* Note, u, v are signed long! */
	 /* u*g1*x^2 */
	 mpz_mul (tmp, g[1], u);
	 mpz_add (fuv[2], f[2], tmp);

	 /* (g0*u* + v*g1)*x */
	 mpz_mul (tmp, g[0], u);
	 mpz_mul (tmp1, g[1], v);
	 mpz_add (tmp, tmp, tmp1);
	 mpz_add (fuv[1], f[1], tmp);

	 /* v*g0 */
	 mpz_mul (tmp, g[0], v);
	 mpz_add (fuv[0], f[0], tmp);

	 mpz_clear (tmp);
	 mpz_clear (tmp1);
}


/*
  Compute fuv = f+(u*x+v)*g,
  f(r) + u*r*g(r) + v*g(r) = 0
  The inputs for f and g are mpz.
*/
static inline void
compute_fuv ( mpz_t *fuv,
			  mpz_t *f,
			  mpz_t *g,
			  int d,
			  long u,
			  long v)
{
	 mpz_t tmp, tmp1;
	 mpz_init (tmp);
	 mpz_init (tmp1);
	 int i = 0;

	 for (i = 3; i <= d; i ++)
		  mpz_set (fuv[i], f[i]);

	 /* f + u*g1*x^2
		+ (g0*u* + v*g1)*x
		+ v*g0 */

	 /* Note, u, v are signed long! */
	 /* u*g1*x^2 */
	 mpz_mul_si (tmp, g[1], u);
	 mpz_add (fuv[2], f[2], tmp);

	 /* (g0*u* + v*g1)*x */
	 mpz_mul_si (tmp, g[0], u);
	 mpz_mul_si (tmp1, g[1], v);
	 mpz_add (tmp, tmp, tmp1);
	 mpz_add (fuv[1], f[1], tmp);

	 /* v*g0 */
	 mpz_mul_si (tmp, g[0], v);
	 mpz_add (fuv[0], f[0], tmp);

	 mpz_clear (tmp);
	 mpz_clear (tmp1);
}

/*
  Compute fuv = f+(u*x+v)*g,
  The inputs for f and g are unsigne long.
  Note, u, v must be unsigned long.
  So they should be reduce (mod p) if necessary.
*/
static inline void
compute_fuv_ul ( unsigned long *fuv_ul,
				 unsigned long *f_ul,
				 unsigned long *g_ul,
				 int d,
				 unsigned long u,
				 unsigned long v,
				 unsigned long p )
{
	 int i;
	 modulusul_t mod;
	 residueul_t tmp, tmp1, tmp2;
	 modul_initmod_ul (mod, p);
	 modul_init (tmp, mod);
	 modul_init (tmp1, mod);
	 modul_init (tmp2, mod);

	 for (i = 3; i <= d; i ++)
		  fuv_ul[i] = f_ul[i];

	 /* f + u*g1*x^2
		+ (g0*u* + v*g1)*x
		+ v*g0 */

	 /* u*g1*x^2 */
	 modul_set_ul (tmp, g_ul[1], mod);
	 modul_set_ul (tmp2, u, mod);
	 modul_mul (tmp, tmp, tmp2, mod);
	 modul_set_ul (tmp1, f_ul[2], mod);
	 modul_add (tmp, tmp, tmp1, mod);
	 fuv_ul[2] = modul_get_ul(tmp, mod);

	 /* (g0*u* + v*g1)*x */
	 modul_set_ul (tmp, g_ul[1], mod);
	 modul_set_ul (tmp1, v, mod);
	 modul_mul (tmp, tmp, tmp1, mod);
	 modul_set_ul (tmp1, g_ul[0], mod);
	 // tmp2 = u as set above.
	 modul_mul (tmp1, tmp1, tmp2, mod);
	 modul_add (tmp, tmp, tmp1, mod);
	 modul_set_ul (tmp1, f_ul[1], mod);
	 modul_add (tmp, tmp, tmp1, mod);
	 fuv_ul[1] = modul_get_ul(tmp, mod);

	 /* v*g0 */
	 modul_set_ul (tmp1, v, mod);
	 modul_set_ul (tmp2, g_ul[0], mod);
	 modul_mul (tmp1, tmp1, tmp2, mod);
	 modul_set_ul (tmp, f_ul[0], mod);
	 modul_add (tmp, tmp, tmp1, mod);
	 fuv_ul[0] = modul_get_ul(tmp, mod);
}


/*
  Compute v (mod p) by
  f(r) + u*r*g(r) + v*g(r) = 0 (mod p).
  The inputs for f(r) and g(r) are unsigned long.
*/
static inline unsigned long
compute_v_ul ( unsigned long fx,
			   unsigned long gx,
			   unsigned long r,
			   unsigned long u,
			   unsigned long p)
{
	 modulusul_t mod;
	 residueul_t tmp, tmp1;
	 unsigned long v;
	 modul_initmod_ul (mod, p);
	 modul_init (tmp, mod);
	 modul_init (tmp1, mod);

	 /* g(r)*r*u + f(r) */
	 modul_set_ul (tmp, gx, mod);
	 modul_set_ul (tmp1, r, mod);
	 modul_mul (tmp, tmp, tmp1, mod);
	 modul_set_ul (tmp1, u, mod);
	 modul_mul (tmp, tmp, tmp1, mod);
	 modul_set_ul (tmp1, fx, mod);
	 modul_add (tmp, tmp, tmp1, mod);
	 v = modul_get_ul(tmp, mod);

	 /* solve v in tmp2 + v*g(r) = 0 (mod p) */
	 v = solve_lineq(v, gx, 0, p);

	 modul_clear (tmp, mod);
	 modul_clear (tmp1, mod);
	 modul_clearmod (mod);
	 return v;
}

/*
  Evaluation polynomials at many points.
*/
static void
eval_polys ( mpz_t *f,
			 mpz_t *g,
			 mpz_t *fr,
			 mpz_t *gr,
			 mpz_t *numerator,
			 unsigned int *primes,
			 int d )
{
	 unsigned long k;
	 mpz_t tmp;
	 mpz_init (tmp);

	 for ( k = 0; k <= primes[NP-1]; k ++ ) {
		  mpz_set_ui (fr[k], 0);
		  mpz_set_ui (gr[k], 0);
		  eval_poly_ui (fr[k], f, d, k);
		  eval_poly_ui (gr[k], g, 1, k);
		  eval_poly_diff_ui (numerator[k], f, d, k);
		  mpz_mul (numerator[k], numerator[k], gr[k]);
		  mpz_neg (numerator[k], numerator[k]);
		  mpz_mul (tmp, fr[k], g[1]);
		  mpz_add (numerator[k], tmp, numerator[k]);
		  /* if (DEBUG) { */
		  /* 	   gmp_printf ("numerator[%d]: %Zd\n", k, numerator[k]); */
		  /* 	   gmp_printf ("gr[%d]: %Zd", k, gr[k]); */
		  /* } */
	 }

	 mpz_clear (tmp);
}


/*
  Compute v = f(r) (mod pe), where f is of degree d.
  The input f should be unsigned long.
*/
static inline unsigned long
eval_poly_ui_mod ( unsigned long *f,
				   int d,
				   unsigned long r,
				   unsigned long pe )
{


	 int i;
	 modulusul_t mod;
	 residueul_t vtmp, rtmp, tmp;
	 unsigned long v;

	 modul_initmod_ul (mod, pe);
	 modul_init (vtmp, mod);
	 modul_init (rtmp, mod);
	 modul_init (tmp, mod);

	 /* set vtmp = f[d] (mod p) and rtmp = r (mod p) */
	 modul_set_ul (vtmp, f[d], mod);
	 modul_set_ul (rtmp, r, mod);

	 for (i = d - 1; i >= 0; i--) {
		  modul_mul (vtmp, vtmp, rtmp, mod);
		  modul_set_ul (tmp, f[i], mod);
		  modul_add (vtmp, tmp, vtmp, mod);
	 }

	 v = modul_get_ul (vtmp, mod);
	 modul_clear (vtmp, mod);
	 modul_clear (rtmp, mod);
	 modul_clear (tmp, mod);
	 modul_clearmod (mod);

	 return v;
}


/*
  Compute v = f'(r) (mod pe), where f is of degree d.
  The input f should be unsigned long.
*/
static inline unsigned long
eval_poly_diff_ui_mod ( unsigned long *f,
						int d,
						unsigned long r,
						unsigned long pe )
{
	 int i;
	 modulusul_t mod;
	 residueul_t vtmp, rtmp, tmp, itmp;
	 unsigned long v;

	 modul_initmod_ul (mod, pe);
	 modul_init (vtmp, mod);
	 modul_init (rtmp, mod);
	 modul_init (itmp, mod);
	 modul_init (tmp, mod);

	 /* set vtmp = f[d] (mod p) and rtmp = r (mod p) */
	 modul_set_ul (vtmp, f[d], mod);
	 modul_set_ul (tmp, (unsigned long) d, mod);
	 modul_mul (vtmp, vtmp, tmp, mod);
	 modul_set_ul (rtmp, r, mod);

	 for (i = d - 1; i >= 1; i--) {
		  modul_mul (vtmp, vtmp, rtmp, mod);
		  /* vtmp <- vtmp + i*f[i] */
		  modul_set_ul (tmp, f[i], mod);
		  modul_set_ul (itmp, i, mod);
		  modul_mul (tmp, itmp, tmp, mod);
		  modul_add (vtmp, tmp, vtmp, mod);
	 }

	 v = modul_get_ul (vtmp, mod);
	 modul_clear (vtmp, mod);
	 modul_clear (itmp, mod);
	 modul_clear (rtmp, mod);
	 modul_clear (tmp, mod);
	 modul_clearmod (mod);

	 return v;
}


/*
  Test whether r is a root for f(r) + g(r)*(u*x+v) = 0 (mod p^e);
  Note that, the condition to call this function is that, r is a root
  for f(r) + g(r)*(u*x+v) = 0 (mod p^(e-1)).

  (1) If it can be lifted, then it could be single or multiple root
  -- single root return 1, -- note, the lifted root is not r in general.
  Hence, we need to compute the lifted root and save it in r_lifted.
  -- multiple root return 2.
  (2) If it is can not be lifted, return 0;
*/
static inline int
isroot_fuv ( mpz_t *f,
			 mpz_t *g,
			 int d,
			 unsigned long u,
			 unsigned long v,
			 unsigned long r,
			 unsigned long pe,
			 unsigned long p,
			 unsigned long *r_lifted )
{
	 mpz_t tmp, tmp2, fr, gr;
	 int type;
	 mpz_init (fr);
	 mpz_init (gr);
	 mpz_init (tmp);
	 mpz_init (tmp2);

	 /* r is a root over p^(e-1), we need to see whether it is a
		single or multiple;
		- If it is a single, computer r_lifted
		- If it is a multiple, see whether it can be lifted
		-- consider f(r) (mod p^e)	 */

	 /* f'(x)
		+ u*g(x)
		+ g'(x)*(ux + v) */
	 mpz_set_ui (tmp, r);
	 mpz_mul_ui (tmp, tmp, u);
	 mpz_add_ui (tmp, tmp, v);
	 mpz_mul (tmp, tmp, g[1]);
	 eval_poly_ui (gr, g, 1, r);
	 mpz_mul_ui (gr, gr, u);
	 mpz_add (gr, gr, tmp);
	 eval_poly_diff_ui (fr, f, d, r); // keep fr unchanged for further reference.
	 mpz_add (tmp, fr, gr); // tmp is f'_uv(r)

	 /* cmpute f(x) + g(x)*(ux + v) */
	 mpz_set_ui (tmp2, r);
	 mpz_mul_ui (tmp2, tmp2, u);
	 mpz_add_ui (tmp2, tmp2, v);
	 eval_poly_ui (gr, g, 1, r);
	 mpz_mul (gr, gr, tmp2);
	 eval_poly_ui (fr, f, d, r);
	 mpz_add (fr, fr, gr);  // fr is f_uv(r)

	 /* if f_uv'(r) = 0 (mod p), then it is a double root */
	 if ( mpz_fdiv_ui (tmp, p) == 0) {
		  if (mpz_fdiv_ui (fr, pe) == 0)
			   type = 2;
		  else
			   type = 0;
	 }
	 /* r is a single root, computer r_lifted. */
	 else {
		  unsigned long tmp3;
		  modulusul_t mod;
		  residueul_t res1, res2;

		  /* 1/f_uv'(r) (mod p) */
		  tmp3 = mpz_fdiv_r_ui (tmp, tmp, p);
		  modul_initmod_ul (mod, p);
		  modul_init (res1, mod);
		  modul_set_ul (res1, tmp3, mod);
		  modul_inv (res1, res1, mod);
		  /* -f_uv(r)/p^(e-1) (mod p) */
		  mpz_fdiv_q_ui (tmp2, fr, pe/p);
		  mpz_neg (tmp2, tmp2);
		  tmp3 = mpz_fdiv_r_ui (tmp2, tmp2, p);
		  modul_init (res2, mod);
		  modul_set_ul (res2, tmp3, mod);
		  modul_mul (res1, res1, res2, mod);

		  (*r_lifted) = modul_get_ul(res1, mod);
		  (*r_lifted) = (*r_lifted) * (pe/p) + r;
		  //printf ("lifted_r: %lu\n", *r_lifted);
		  type = 1;

		  modul_clear (res1, mod);
		  modul_clear (res2, mod);
		  modul_clearmod (mod);
	 }

	 mpz_clear (fr);
	 mpz_clear (gr);
	 mpz_clear (tmp);
	 mpz_clear (tmp2);
	 return type;
}


/*
  Do the same thing but consider f_uv in ul.
  Given f (mod pe) and (u, v) (mod pe) pair.
*/
static inline int
isroot_fuv_ul ( unsigned long *f,
				unsigned long *g,
				int d,
				unsigned long u,
				unsigned long v,
				unsigned long r,
				unsigned long p,
				unsigned long pe,
				unsigned long *r_lifted )
{
	 int type;
	 modulusul_t mod, modpe;
	 residueul_t tmp, tmp1, tmp2;
	 unsigned long gr, fr;
	 modul_initmod_ul (mod, p);
	 modul_init (tmp, mod);
	 modul_init (tmp1, mod);
	 modul_init (tmp2, mod);

	 /* r is a root over p^(e-1) = pe / p, we need to see
		whether it is a	single or multiple root;
		- If it is a single, computer r_lifted
		- If it is a multiple, see whether it can be lifted
		-- consider f(r) (mod p^e)	 */

	 /* f_uv'(x) = f'(x) + u*g(x)+ g'(x)*(ux + v) */

	 /* g'(x)*(ux + v) */
	 modul_set_ul (tmp, u, mod);
	 modul_set_ul (tmp1, r, mod);
	 modul_mul (tmp, tmp, tmp1, mod);
	 modul_set_ul (tmp1, v, mod);
	 modul_add (tmp, tmp, tmp1, mod);
	 modul_set_ul (tmp1, g[1], mod);
	 modul_mul (tmp, tmp, tmp1, mod);

	 /* u*g(x) + g'(x)*(ux + v) */
	 gr = eval_poly_ui_mod (g, 1, r, pe);
	 modul_set_ul (tmp1, gr, mod);
	 modul_set_ul (tmp2, u, mod);
	 modul_mul (tmp2, tmp1, tmp2, mod);
	 modul_add (tmp, tmp, tmp2, mod);

	 /* f'(x) */
	 fr = eval_poly_diff_ui_mod (f, d, r, p);
	 modul_set_ul (tmp2, fr, mod);
	 modul_add (tmp, tmp, tmp2, mod); // tmp should be kept unchanged from now on.

	 /* f_uv(r) = f(r) + g(r)*(ur + v) (mod pe) */
	 modul_clear (tmp1, mod);
	 modul_clear (tmp2, mod);
	 modul_initmod_ul (modpe, pe);
	 modul_init (tmp1, modpe);
	 modul_init (tmp2, modpe);
	 /* tmp1 = g(r)*(ur + v) */
	 modul_set_ul (tmp1, u, modpe);
	 modul_set_ul (tmp2, r, modpe);
	 modul_mul (tmp1, tmp1, tmp2, modpe);
	 modul_set_ul (tmp2, v, modpe);
	 modul_add (tmp1, tmp1, tmp2, modpe);
	 modul_set_ul (tmp2, gr, modpe);
	 modul_mul (tmp1, tmp1, tmp2, modpe);
	 /* f(r) */
	 fr = eval_poly_ui_mod (f, d, r, pe);
	 modul_set_ul (tmp2, fr, modpe);
	 modul_add (tmp1, tmp1, tmp2, modpe);
	 modul_clear (tmp1, modpe);

	 /* if f_uv(r) = 0 (mod pe), then it is a root. We
		need to check whether it is single or multiple. */
	 if ( modul_get_ul (tmp1, modpe) == 0 ) {

		  if ( modul_get_ul (tmp, mod) == 0 )
			   type = 2;

		  else {

			   modul_init (tmp1, mod);

			   /* 1/f'(r) (mod p) */
			   modul_inv (tmp, tmp, mod);

			   /* -f(r)/p^(e-1) (mod p) */
			   fr = fr / (pe/p);
			   modul_set_ul (tmp1, fr, mod);
			   modul_neg (tmp1, tmp1, mod);
			   modul_mul (tmp, tmp, tmp1, mod);

			   (*r_lifted) = modul_get_ul(tmp, mod);
			   (*r_lifted) = (*r_lifted) * (pe/p) + r;
			   //printf ("lifted_r: %lu\n", *r_lifted);
			   type = 1;

			   modul_clear (tmp1, mod);
		  }
	 }
	 /* r is not a root mod (pe) */
	 else
		  type = 0;

	 modul_clear (tmp2, modpe);
	 modul_clear (tmp, mod);
	 modul_clearmod (mod);
	 modul_clearmod (modpe);
	 return type;
}


/*
  Given f (mod pe) and (u, v) (mod pe) pair.
  Lift a single root r (mod pe/p) to r (mod p)
  Please ensure that f_{u, v}(r) = 0 (mod pe/p)
  when call this function.
*/
static inline void
liftroot_fuv_ul ( unsigned long *f,
				  int d,
				  unsigned long r,
				  unsigned long p,
				  unsigned long pe,
				  unsigned long *r_lifted )
{
	 modulusul_t mod;
	 residueul_t tmp, tmp1;
	 unsigned long fr;
	 modul_initmod_ul (mod, p);
	 modul_init (tmp, mod);
	 modul_init (tmp1, mod);

	 /* f'(x) */
	 fr = eval_poly_diff_ui_mod (f, d, r, p);
	 modul_set_ul (tmp1, fr, mod);
	 modul_add (tmp, tmp, tmp1, mod); // tmp should be kept unchanged from now on.

	 /* 1/f'(r) (mod p) */
	 modul_inv (tmp, tmp, mod);

	 /* -f(r)/p^(e-1) (mod p) */
	 fr = fr / (pe/p);
	 modul_set_ul (tmp1, fr, mod);
	 modul_neg (tmp1, tmp1, mod);
	 modul_mul (tmp, tmp, tmp1, mod);

	 (*r_lifted) = modul_get_ul(tmp, mod);
	 (*r_lifted) = (*r_lifted) * (pe/p) + r;

	 modul_clear (tmp1, mod);
	 modul_clear (tmp, mod);
	 modul_clearmod (mod);
}


/*
  Reduce mpz_t *f to unsigned long *f_mod;
  Given modulus pe, return f (mod pe).
*/
static inline void
reduce_poly_ul ( unsigned long *f_ul,
				 mpz_t *f,
				 int d,
				 unsigned long pe )
{
	 int i;
	 for (i = 0; i <= d; i ++) {
		  f_ul[i] = mpz_fdiv_ui (f[i], pe);
	 }
}


/*-----------------------------*/
/*  @Good (u, v) sublattices   */
/*-----------------------------*/


/*
  Find good sublattice, the lifted cases.
*/
static void
find_sublattice_lift ( node *firstchild,
					   listnode **top,
					   unsigned long * f_ul,
					   unsigned long * g_ul,
					   unsigned long * fuv_ul,
					   int d,
					   unsigned int p,
					   unsigned int e,
					   unsigned int curr_e )
{
	 if (firstchild == NULL || curr_e > e)
		  return;

	 /* some tmp variables */
	 unsigned short nd, ns;
	 unsigned int i, k, nroots;
	 unsigned long pe, pem1, r_lifted[1], *r, uu, vv, fr, gr;
	 node *currnode, *tmpnode, *l1node;
	 uu = vv = r_lifted[0] = 0;

	 /* until now, all the siblings to the left of "firstchild" has
		been considered. Therefore, we only consider the "nextsibling".
		hence we call it "first" sibling. */
	 currnode = firstchild;
	 tmpnode = l1node = NULL;

	 /* compute p^e */
	 pem1 = 1UL;
	 for (i = 0; i < curr_e - 1; i ++)
		  pem1 = pem1 * p;
	 pe = pem1 * p;

	 /* loop until all siblings are checked. */
	 while ( (currnode != NULL) ) {

		  /* find level 1 nodes, those (u, v) (mod p). */
		  l1node = currnode;
		  while (l1node->e != 1)
			   l1node = l1node->parent;

		  if (DEBUG) {
			   printf("-----\n");
			   printf("p: %u, e: %u -> %u\n", p, curr_e - 1, curr_e);
			   printf("u-v pair: (%lu, %lu) has roots: \n",
					  currnode->u, currnode->v);
			   for (i = 0; i < (currnode->nr); i++)
					printf("\tr[%u]: %lu\n", i, currnode->r[i]);
			   //print_tree(l1node->parent, 0);
			   printf("-----\n");
		  }

		  /* save current roots */
		  r = (unsigned long *) malloc ( (currnode->nr) * sizeof (unsigned long) );
		  if (r == NULL) {
			   fprintf (stderr, "Error, cannot allocate memory in find_sublattice_lift(). \n");
			   exit (1);
		  }
		  for (i = 0; i < currnode->nr; i ++)
			   r[i] = 0;

		  /* save all the multiple and single roots for this node. */
		  nd = 0;
		  ns = 0;
		  for (nroots = 0; nroots < (currnode->nr); nroots++) {
			   /* search all roots of the (u, v) (mod p) to find
				  the level 1 ancestor of the current root. */
			   for (i = 0; i < (l1node->nr); i++) {
					if ( l1node->r[i] == (currnode->r[nroots] % p) ) {
						 if (l1node->roottype[i] == 2) {
							  r[nd] = currnode->r[nroots];
							  nd ++;
						 }
						 else if (l1node->roottype[i] == 1) {
							  r[currnode->nr - 1 - ns] = currnode->r[nroots];
							  ns ++;
						 }
						 else {
							  fprintf (stderr, "Error, roottype wrong in find_sublattice_lift(). \n");
							  exit (1);
						 }
					}
			   }
		  }

		  /* Now a pair (u, v) is fixed. Fix this pair,
			 -- r is multiple, solve lifted (u, v) who has this r + i*p^k as
			 multiple roots; Also for these (u, v), lift any possible single r;
			 -- Note if it is single root. We ignore it. r could be a
			 single root for some other pair (u', v') where (u', v') has
			 some other multiple root. However, this situation will be
			 detected by other (u', v') r' pair. */
		  compute_fuv_ul (fuv_ul, f_ul, g_ul, d, currnode->u, currnode->v, pe);

		  for (nroots = 0; nroots < nd; nroots++) {

			   /* compute g(r) */
			   gr = eval_poly_ui_mod (g_ul, 1, r[nroots], pe);
			   if (gr % p == 0)
					continue;

			   /* compute f_uv(x) and then evaluate it at r. */
			   fr = eval_poly_ui_mod (fuv_ul, d, r[nroots], pe);
			   fr = fr / pem1;
			   /* solve on fr + gr*x = 0 (mod p), where x = uu*r + vv. */
			   fr = solve_lineq (fr, gr, 0, p);
			   fr = fr % p;

			   /* we want to solve (uu, vv) in  fr = uu*r + vv (mod p).
				  - if r is not invertible, fix vv and loop all uu.
				  - otherwise, fix vv and solve uu. */
			   if (r[nroots] % p == 0) {

					for (uu = 0; uu < p; uu ++) {
						 if (DEBUG)
							  printf ("fr: %lu, r: %lu,  (uu, vv): (%lu, %lu) -> (%lu, %lu) (non-invertible, multiple) \n",
									  fr, r[nroots], uu, fr, currnode->u + pem1 * uu,
									  currnode->v + pem1 * fr);

						 /* since now r is a multiple root, add r + k * p^{e-1}. */
						 for (k = 0; k < p; k ++) {
							  insert_node (currnode, &tmpnode, currnode->u + pem1 * uu,
										   currnode->v + fr * pem1, r[nroots] + k * pem1, curr_e, pe, 0);
						 }

						 /* we want to find all the lifted single roots for this pair (u, v) */
						 for (k = nd; k < (currnode->nr); k++) {
							  r_lifted[0] = 0;
							  compute_fuv_ul (fuv_ul, f_ul, g_ul, d, currnode->u + pem1 * uu,
											  currnode->v + fr * pem1, pe);
							  liftroot_fuv_ul (fuv_ul, d, r[k], p, pe, r_lifted);
							  insert_node (currnode, &tmpnode, currnode->u + pem1 * uu,
										   currnode->v + fr * pem1, r_lifted[0], curr_e, pe, 0);
						 }
					}
			   }
			   else {
					for (vv = 0; vv < p; vv ++) {
						 uu = solve_lineq (vv, r[nroots], fr, p);
						 if (DEBUG)
							  printf ("fr: %lu, r: %lu, (uu, vv): (%lu, %lu) -> (%lu, %lu) (invertible, multiple) \n",
									  fr, r[nroots], uu, vv, currnode->u + pem1 * uu, currnode->v + pem1 * vv);

						 /* since now r is a multiple root, add r + k * p^{e-1}. */
						 for (k = 0; k < p; k ++) {
							  insert_node (currnode, &tmpnode, currnode->u + pem1 * uu,
										   currnode->v + pem1 * vv, r[nroots] + k * pem1, curr_e, pe, 0);
						 }

						 /* we want to find all the lifted single roots for this pair (u, v) */
						 for (k = nd; k < (currnode->nr); k++) {
							  r_lifted[0] = 0;
							  compute_fuv_ul (fuv_ul, f_ul, g_ul, d, currnode->u + pem1 * uu,
											  currnode->v + pem1 * vv, pe);
							  liftroot_fuv_ul (fuv_ul, d, r[k], p, pe, r_lifted);
							  insert_node (currnode, &tmpnode, currnode->u + pem1 * uu,
										   currnode->v + pem1 * vv, r_lifted[0], curr_e, pe, 0);
						 }
					}
			   }
		  } // consider next root of current (u, v)

		  free (r);
		  find_sublattice_lift (currnode->firstchild, top, f_ul, g_ul, fuv_ul, d, p, e, curr_e + 1);

		  /* If current node is the 2nd bottom leave, add the bottom level leaves
			 with highest valuations to the list and delete them. */
		  if (curr_e == e) {
			   tmpnode = currnode->firstchild;
			   while (tmpnode != NULL) {
					insert_listnode (top, tmpnode->u, tmpnode->v, tmpnode->val);
					l1node = tmpnode;
					tmpnode = tmpnode->nextsibling;
					//printf ("deleting ... (%lu, %lu)\n", l1node->u, l1node->v);
					free_node (&l1node);
			   }
		  }

		  /* delete current node and move to next sibling. */
		  tmpnode = currnode;
		  currnode = currnode->nextsibling;
		  if (currnode != NULL)
			   (currnode->parent)->firstchild = currnode;
		  //printf ("deleting ... (%lu, %lu)\n", tmpnode->u, tmpnode->v);
		  free_node (&tmpnode);
	 }
	 return;
}


/*
  Find sublattices, the base case.

  Only consider those (u, v) which has at least one multiple root;
  Ignore those pairs which have no multiple root. Note, for the
  (u, v) pairs, we consider all the roots (simple + mul) of them.
*/
static void
find_sublattice ( listnode **top,
				  rsstr_t rs,
				  unsigned int p,
				  unsigned int e )
{
	 mpz_t fx_tmp, gx_tmp, numerator_tmp, tmp;
	 unsigned long pe, i, r, v, a, b, u, *f_ul, *g_ul, *fuv_ul, r_lifted[1];
	 node *currnode, *root;
	 int k;

	 mpz_init (fx_tmp);
	 mpz_init (gx_tmp);
	 mpz_init (tmp);
	 mpz_init (numerator_tmp);

	 f_ul = (unsigned long*) malloc ((rs->d + 1) * sizeof (unsigned long));
	 fuv_ul = (unsigned long*) malloc ((rs->d + 1) * sizeof (unsigned long));
	 g_ul = (unsigned long*) malloc ((2) * sizeof (unsigned long));
	 if ((f_ul == NULL) || (g_ul == NULL) || (fuv_ul == NULL)) {
		  fprintf (stderr, "Error, cannot allocate memory in find_sublattice(). \n");
		  exit (1);
	 }

	 /* compute p^e */
	 pe = 1UL;
	 for (i = 0; i < e; i ++)
		  pe = pe * p;

	 /* compute f (mod pe) */
	 reduce_poly_ul (f_ul, rs->f, rs->d, pe);
	 reduce_poly_ul (g_ul, rs->g, 1, pe);

     /* new (u, v, val) tree */
	 new_tree(&root);
	 root = new_node ();

	 /* for each root 0 <= r < p  */
	 for (r = 0; r < p; r ++) {

		  /* set f(r), g(r), numerator(r) */
		  if (r < primes[NP-1]) {
			   mpz_set (fx_tmp, rs->fx[r]);
			   mpz_set (gx_tmp, rs->gx[r]);
			   mpz_set (numerator_tmp, rs->numerator[r]);
		  }
		  else {
			   fprintf (stderr, "Error, something strange in find_sublattice(). \n");
			   exit(1);
		  }

		  /* ignore those non-invertible g[r] */
		  if (mpz_divisible_ui_p(gx_tmp, p) != 0)
			   continue;

		  /* solve u in g(r)^2 * u = numerators_ext(r) (mod p) */
		  b = mpz_fdiv_ui (numerator_tmp, p);
		  b = p - b;
		  mpz_mul (tmp, gx_tmp, gx_tmp);
		  a = mpz_fdiv_ui (tmp, p);

		  /* call solve_lineq function on b + a*su = 0 (mod p) */
		  u = solve_lineq (b, a, 0, p);

		  /* only need to consider the special u pair (u, v), we
			 want to solve v given u and r. */
		  mpz_mul_ui (tmp, gx_tmp, r);
		  mpz_mul_ui (tmp, tmp, u);
		  mpz_add (tmp, tmp, fx_tmp);
		  b = mpz_fdiv_ui (tmp, p);
		  a = mpz_fdiv_ui (gx_tmp, p);
		  v = solve_lineq (b, a, 0, p);

		  /* For this (u, v) pair which is already known to have a
			 multiple root r, we search for any other possible single
			 and  multiple roots (not necessary for multiple since other r's
			 will generate the same (u, v). ). We exhaustively search. */
		  for (i = 0; i < p; i ++) {

			   /* r_lifted is not correct, but not important here. Only k is used. */
			   k  = isroot_fuv_ul (f_ul, g_ul, rs->d, u, v, i, p, p, r_lifted);
			   if (DEBUG)
					printf ("(u, v): %lu, %lu  r: %lu, is_root: %d\n", u, v, i, k);

			   /* if i is some root, either single or multiple. */
			   if (k != 0) {
					/* insert i, which does not equal to r, if exits.*/
					insert_node (root, &currnode, u, v, i, 1, p, k);
			   }
		  }
	 }

	 /* if e == 1, add nodes to list top. Note that, we add all nodes
		since now e == 1 is too small. It might be inaccurate if we only
		add nodes with best valuation. For accuracy considerations, we
		add all nodes to the list. */
	 if (e == 1) {
		  node *tmpnode;
		  currnode = root->firstchild;
		  while (currnode != NULL) {
			   insert_listnode_plain (top, currnode->u, currnode->v, currnode->val);
			   tmpnode = currnode;
			   currnode = currnode->nextsibling;
			   free_node (&tmpnode);
		  }
		  tmpnode = NULL;
	 }
	 /* if e > 1, lift to higher p^e */
	 else {
		  find_sublattice_lift (root->firstchild, top, f_ul, g_ul, fuv_ul, rs->d, p, e, 2);
	 }

	 free_node(&root);

     /* clear */
	 free (f_ul);
	 free (fuv_ul);
	 free (g_ul);
	 mpz_clear (tmp);
	 mpz_clear (fx_tmp);
	 mpz_clear (gx_tmp);
	 mpz_clear (numerator_tmp);
}

#if 0

/*
  Init a 2D arrary of type ul, the
  dimension X by Y.

*/
static void
init_2D_ld ( long ***array,
			 const unsigned long X,
			 const unsigned long Y )
{
	 unsigned long i, j;

	 /* alloc */
	 (*array) = (long **) malloc (X * sizeof (long *));
	 if ( (*array) != NULL) {
		  for (i = 0; i < X; i ++) {
			   (*array)[i] = (long *) malloc ( Y * sizeof(long) );
			   if ( (*array)[i] == NULL) {
					fprintf (stderr, "Error, cannot allocate memory in init_2D_ld(). \n");
					exit (1);
			   }
		  }
	 }
	 else {
		  fprintf (stderr, "Error, cannot allocate memory in main\n");
		  exit (1);
	 }

	 /* init */
	 for (i = 0; i < X; i ++)
		  for (j = 0; j < Y; j ++)
			   (*array)[i][j] = 0UL;
}


/*
  Free a 2D arrary of type ul, the
  dimension X by Y.
*/
static void
free_2D_ld ( long ***array,
			 const unsigned long X )
{
	 unsigned long i;

	 for (i = 0; i < X; i++) {
		  if ((*array)[i])
			   free ((*array)[i]);
	 }

	 if (*array)
		  free (*array);
}


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
							 unsigned long *pe,
							 unsigned long *ind,
							 unsigned long ***individual_sublattices,
							 sublattice_pq *pqueue )
{
	 int i;
	 mpz_t tmpp1, tmpp2, tmpu1, tmpu2, tmpv, re;
	 mpz_init (tmpp1);
	 mpz_init (tmpp2);
	 mpz_init (tmpu1);
	 mpz_init (tmpu2);
	 mpz_init (tmpv);
	 mpz_init (re);

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
	 /* < 0 by construction */
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

		  insert_sublattice_pq ( pqueue, tmpu1, tmpu2 );
	 }

	 mpz_clear (tmpp1);
	 mpz_clear (tmpp2);
	 mpz_clear (tmpu1);
	 mpz_clear (tmpu2);
	 mpz_clear (tmpv);
	 mpz_clear (re);
}


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

	 unsigned short i, fail = 0;
	 unsigned long j, *pe, *size, *tsize, *ind, ***individual_sublattices;
	 listnode *top, *tmp;

	 /* Sublattice[i][length][] save (u, v) for prime[i] */
	 individual_sublattices = (unsigned long ***)
		  malloc ( rsparam->tlen_e_sl * sizeof (unsigned long **) );
	 size = (unsigned long *) malloc ( rsparam->tlen_e_sl * sizeof (unsigned long));
	 tsize = (unsigned long *) malloc ( rsparam->tlen_e_sl * sizeof (unsigned long));
	 ind = (unsigned long *) malloc ( rsparam->tlen_e_sl * sizeof (unsigned long));
	 pe = (unsigned long *) malloc ( rsparam->tlen_e_sl * sizeof (unsigned long));
	 mpz_set_ui (rsparam->modulus, 1UL);

	 if ( ( (individual_sublattices) == NULL) || (size == NULL)
		  || (tsize == NULL) || (ind == NULL) || (pe == NULL) ) {
		  fprintf (stderr,
				   "Error, cannot allocate memory in return_all_sublattices(). \n");
		  exit (1);
	 }

	 /* For each prime[i], lift the polynomial
		and return (u, v) with best score. */
	 for (i = 0; i < rsparam->tlen_e_sl; i ++) {

		  pe[i] = 1;
		  for (j = 0; j < rsparam->e_sl[i]; j ++)
			   pe[i] *= primes[i];

		  new_list (&top);

		  /* find individual sublattices */
		  find_sublattice ( &top, rs,
							primes[i], rsparam->e_sl[i] );

		  tsize[i] = count_list (top);
		  if (tsize[i] == 0) {
			   fail = 1;
			   break;
		  }

		  size[i] = tsize[i];

		  if (rs->d == 6 && rsparam->ncrts_sl > 0) {
			   if (size[i] > rsparam->ncrts_sl)
					size[i] = rsparam->ncrts_sl;
		  }

		  if (verbose)
			   fprintf ( stderr,
						 "# Info: p: %2u, p^e: %6lu, list_size: %6lu, size_cutoff: %6lu\n",
						 primes[i], pe[i], tsize[i], size[i] );

		  /* allocate array for the i-th prime. */
		  (individual_sublattices)[i] = (unsigned long **)
			   malloc ( size[i] * sizeof(unsigned long *) );

		  if ( (individual_sublattices)[i] == NULL ) {
			   fprintf (stderr, "Error, cannot allocate memory in return_all_sublattices(). \n");
			   exit (1);
		  }

		  tmp = top;
		  for (j = 0; j < size[i]; j ++) {
			   (individual_sublattices)[i][j] =
					(unsigned long *) malloc ( 2 * sizeof(unsigned long));
			   if ( individual_sublattices[i][j] == NULL ) {
					fprintf (stderr, "Error, cannot allocate memory in return_all_sublattices(). \n");
					exit (1);
			   }
			   individual_sublattices[i][j][0] = tmp->u;
			   individual_sublattices[i][j][1] = tmp->v;
			   tmp = tmp->next;
			   //printf ("#%lu, (%lu, %lu)\n", j, individual_sublattices[i][j][0], individual_sublattices[i][j][1]);
		  }
		  free_list (&top);
	 }

	 /* If individual list has 0 len, this set of parameters fails */
	 if (fail == 1) {
		  for (i = 0; i < rsparam->tlen_e_sl; i ++) {
			   for (j = 0; j < size[i]; j ++) {
					free (individual_sublattices[i][j]);
			   }
			   free (individual_sublattices[i]);
		  }
		  free (individual_sublattices);
		  free (size);
		  free (ind);
		  free (pe);
		  return 0;
	 }

	 /* Compute rsparam->modulus */
	 j = 1;
	 for (i = 0; i < rsparam->tlen_e_sl; i ++) {
		  mpz_mul_ui (rsparam->modulus, rsparam->modulus, pe[i]);
		  j *= size[i];
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
														   pe,
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
																pe,
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
																	 pe,
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
																		  pe,
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
																			   pe,
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
																					pe,
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
													  pe,
													  ind,
													  individual_sublattices,
													  pqueue );
	 }
	 /* 2, 3 */
	 else if (rsparam->tlen_e_sl == 2) {
	 	  for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
	 		   for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
					return_all_sublattices_crt ( rsparam,
												 pe,
												 ind,
												 individual_sublattices,
												 pqueue );
	 }
	 /* too aggressive */
	 else {
		  fprintf (stderr, "Error, only len_e_sl: 2, 3, 4, 5, 6, 7, 8 is supported at the moment. \n");
		  exit (1);
	 }

	 /* info */
	 if (verbose)
		  fprintf (stderr, "# Info: computed %lu CRTs\n", j);

	 /* clearence */
	 for (i = 0; i < rsparam->tlen_e_sl; i ++) {
		  for (j = 0; j < size[i]; j ++) {
			   free (individual_sublattices[i][j]);
		  }
		  free (individual_sublattices[i]);
	 }
	 free (individual_sublattices);
	 free (size);
	 free (tsize);
	 free (ind);
	 free (pe);

	 return 1;
}


/*
  Call return_all_sublattice() to return good sublattices.

  Choose the "rsbound->nbest_sl" best ones among them if there
  are more quantities than it. The "best" property is ranked
  by the size of u.
*/
static int
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
			   if (verbose) {
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
	 if (verbose) {
		  fprintf ( stderr,
					"# Info: found %d sublattices with small u, where |u| < %lu (\"rsparam->global_u_bound_rs\")\n",
					pqueue->used - 1, rsparam->global_u_bound_rs);
	 }

	 /* recover, for deg 6 poly */
	 rsparam->global_u_bound_rs = global_u_bound_rs_tmp;

	 return 1;
}


/*-----------------------------*/
/*  @Root sieve                */
/*-----------------------------*/


/*
  root sieve over ARRAY with initial point
  (i, j) (mod pe), where i is fixed, so only
  need to sieve along the line j.
*/
static inline long
rootsieve_run_line ( int16_t *ARRAY,
					 long V,
					 long j,
					 unsigned long pe,
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


#define DEBUG_MULTROOT_LIFT 0
/*
  rootsieve_run_multroot_lift for the multiple root. Since we are
  sure that level 1 node only contains one multiple root
  r, we don't need to (and can't, otherwise, may count
  repeatedly) consider single root at all.
*/
static inline void
rootsieve_run_multroot_lift ( node *currnode,
							  int16_t *ARRAY,
							  unsigned long *f_ul,
							  unsigned long *g_ul,
							  unsigned long *fuv_ul,
							  int d,
							  rsbound_t rsbound,
							  unsigned int p,
							  unsigned int e,
							  unsigned int curr_e,
							  int16_t sub )
{
	 /* recursion end */
	 if (currnode == NULL || curr_e > e || sub == 0)
		  return;

	 /* variables */
	 int16_t subtmp;
	 unsigned int nroots;
	 unsigned long pe, pem1, fr, gr, step;
	 node *tmpnode = NULL, *tmpnode2 = NULL;
	 long j;

	 /* compute p^e */
	 pem1 = 1UL;
	 for (j = 0; j < curr_e - 1; j ++)
		  pem1 = pem1 * p;
	 pe = pem1 * p;

	 /* loop until all siblings are checked. */
	 while ( (currnode != NULL) ) {

		  /* loop all (lifted multiple) roots */
		  for (nroots = 0; nroots < currnode->nr; nroots++) {

			   /* compute g(r) */
			   gr = eval_poly_ui_mod (g_ul, 1, currnode->r[nroots], pe);
			   if (gr % p == 0)
					continue;

			   /* compute f_uv(x) and then evaluate it at r. */
			   compute_fuv_ul (fuv_ul, f_ul, g_ul, d, currnode->u, currnode->v, pe);
			   fr = eval_poly_ui_mod (fuv_ul, d, currnode->r[nroots], pe);
			   fr = fr / pem1;

			   /* solve on fr + gr*x = 0 (mod p), where x = uu*r + vv. */
			   fr = solve_lineq (fr, gr, 0, p); //fr = fr % p;

			   /* If gcd (MOD, p) != 1 && B + MOD*j != v (mod pe), there is no need to
				  record this (u, v). */
			   if (mpz_fdiv_ui (rsbound->MOD, p) == 0) {
					if ( mpz_fdiv_ui (rsbound->B, pe) != (currnode->v + fr * pem1) ) {
						 continue; // no need to insert, skip this loop.
					}
			   }

			   /* insert (u, v), u unchanged, v is fr. to insert roots,
				  r is a multiple root, we need to insert all r + p*j.
				  since if r can be lifted, then all r + p*j are roots */
			   for (j = 0; j < p; j ++) {

#if DEBUG_MULTROOT_LIFT
					fprintf (stderr, "level %u, (%lu, %lu), r: %lu\n",
							 curr_e, currnode->u, currnode->v + fr *pem1, currnode->r[nroots] + j * pem1);
#endif

					insert_node ( currnode, &tmpnode,
								  currnode->u,
								  currnode->v + fr * pem1,
								  currnode->r[nroots] + j * pem1,
								  curr_e, pe, 0 );
			   } // next root of current (u, v)
		  }

		  /* recursieve to next level, curr_e + 1 */
		  rootsieve_run_multroot_lift ( currnode->firstchild,
										ARRAY,
										f_ul,
										g_ul,
										fuv_ul,
										d,
										rsbound,
										p,
										e,
										curr_e + 1,
										sub );

		  /* we are in the second level from bottom, consider all children of this node and sieve */
		  if (curr_e == e) {
			   tmpnode = currnode->firstchild;

			   while (tmpnode != NULL) {

					/* if MOD = 0 (mod p) in B + MOD*j = v (mod pe), no inverse
					   -- if B = v (mod pe), then sieve whole array;
					   -- otherwise, continue; */
					if (mpz_divisible_ui_p (rsbound->MOD, p) != 0) {
						 if ( mpz_fdiv_ui (rsbound->B, pe) == tmpnode->v) {
							  if (mpz_divisible_ui_p (rsbound->MOD, pe ) == 0) {
								   j = 0;
								   step = pe;
							  }
							  /* seems to be redundant since in rootsieve_v(), this
								 situation is discarded */
							  else {
								   tmpnode2 = tmpnode;
								   tmpnode = tmpnode->nextsibling;
								   free_node (&tmpnode2);
								   continue;
							  }
						 }
						 else {
						 	  tmpnode2 = tmpnode;
						 	  tmpnode = tmpnode->nextsibling;
						 	  free_node (&tmpnode2);
						 	  continue;
						 }
					}
					else {
						 j = uv2ab_mod (rsbound->B, rsbound->MOD, tmpnode->v, pe);
						 j = (long) ceil (((double) rsbound->Bmin - (double) j) / (double) pe)
							  * (long) pe + (long) j;
						 j = ab2ij (rsbound->Bmin, j);
						 step = pe;
					}

					/* be careful about this */
					subtmp = (int16_t) ceil ( (double) sub / (double) pem1  * tmpnode->nr);

					rootsieve_run_line ( ARRAY,
										 rsbound->Bmax - rsbound->Bmin,
										 j, step, subtmp );

					tmpnode2 = tmpnode;
					tmpnode = tmpnode->nextsibling;

#if DEBUG_MULTROOT_LIFT
					fprintf (stderr, "deleting bottomnode ... (%lu, %lu) with %u roots in level %u, ",
							 tmpnode2->u, tmpnode2->v, tmpnode2->nr, curr_e);
					fprintf (stderr, "sieving bottomnode ... %d in steps %lu\n", subtmp, step); // pe
#endif

					free_node (&tmpnode2);
			   }
		  }

		  /* delete current node and move to next sibling. */
		  tmpnode = currnode;
		  currnode = currnode->nextsibling;
		  if (currnode != NULL)
			   (currnode->parent)->firstchild = currnode;

		  if (mpz_divisible_ui_p (rsbound->MOD, p) != 0) {
			   if ( mpz_fdiv_ui (rsbound->B, pem1) == tmpnode->v) {
					if (mpz_divisible_ui_p (rsbound->MOD, pem1 ) == 0) {
						 j = 0;
						 step = pem1;
					}
					/* seems to be redundant since in rootsieve_v(), this
					   situation is discarded */
					else {
						 free_node (&tmpnode);
						 continue;
					}
			   }
			   else {
					free_node (&tmpnode);
					continue;
			   }
		  }
		  else {
			   j = uv2ab_mod (rsbound->B, rsbound->MOD, tmpnode->v, pem1);
			   j = (long) ceil (((double) rsbound->Bmin - (double) j) / (double) pem1)
					* (long) pem1 + (long) j;
			   j = ab2ij (rsbound->Bmin, j);
			   step = pem1;
		  }

		  /* be careful about this */
		  subtmp = (int16_t) ceil ( (double) sub / (double) pem1 * (double) p  * tmpnode->nr );
		  rootsieve_run_line ( ARRAY,
							   rsbound->Bmax - rsbound->Bmin,
							   j, step, subtmp );

#if DEBUG_MULTROOT_LIFT
		  fprintf (stderr, "deleting ... (%lu, %lu) with %u roots in level %u, ", tmpnode->u, tmpnode->v, tmpnode->nr, curr_e - 1);
		  fprintf (stderr, "sieving ... %d in steps %lu\n", subtmp, step); //pem1
#endif

		  free_node (&tmpnode);
	 }

	 return;
}


/*
  For multiple root r of the (i, j) (mod p)
*/
static inline  void
rootsieve_run_multroot ( int16_t *ARRAY,
						 rsstr_t rs,
						 rsbound_t rsbound,
						 unsigned int u,
						 long j,
						 unsigned int r,
						 unsigned int p,
						 unsigned int e,
						 int16_t sub )
{

	 /* sieve f(r) + (u*r + (B + MOD * (j+k*p)))*g(r) for k,
		we sieve (mod p) case in the beginning since r would
		be a multiple root for all (j + k*p).
		In the rootsieve_v(), there are two cases in calling
		this function:
		Either flag = 1, e = 1; We can sieve starting from j.
		Note if flag = 2, e > 1. We jump out of this. */
	 if (e <= 1) {
		  rootsieve_run_line ( ARRAY,
							   rsbound->Bmax - rsbound->Bmin,
							   j, p, sub );
#if DEBUG_MULTROOT_LIFT
		  fprintf (stderr, "sieving ... %d in steps %u starting from j: %ld ->, r: %u\n", sub, p, j, r);
#endif
		  return;
	 }

	 /* some variables */
	 unsigned int v;
	 unsigned long pe, *f_ul, *g_ul, *fuv_ul;
	 mpz_t tmpz;

	 mpz_init_set_ui (tmpz, 0UL);

	 /* compute p^e */
	 pe = 1UL;
	 for (v = 0; v < e ; v ++)
		  pe = pe * p;

	 /* use s.p instead of m.p */
	 f_ul = (unsigned long*) malloc ((rs->d + 1) * sizeof (unsigned long));
	 fuv_ul = (unsigned long*) malloc ((rs->d + 1) * sizeof (unsigned long));
	 g_ul = (unsigned long*) malloc ((2) * sizeof (unsigned long));
	 if ((f_ul == NULL) || (g_ul == NULL) || (fuv_ul == NULL)) {
		  fprintf (stderr, "Error, cannot allocate memory in rootsieve_run_multroot(). \n");
		  exit (1);
	 }
	 reduce_poly_ul (f_ul, rs->f, rs->d, pe);
	 reduce_poly_ul (g_ul, rs->g, 1, pe);

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
	 insert_node (root, &tmpnode, u, v, r, 1, p, 2);

	 /* lift to higher p^e */
	 rootsieve_run_multroot_lift ( root->firstchild,
								   ARRAY,
								   f_ul,
								   g_ul,
								   fuv_ul,
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

	 free (f_ul);
	 free (fuv_ul);
	 free (g_ul);
	 mpz_clear (tmpz);
}


#define DEBUG_ROOTSIEVE_V 0
/*
  Root sieve for f + u*x*g(x) + v*g(x) for fixed u, variable v.
*/
static void
rootsieve_v ( int16_t *ARRAY,
			  rsstr_t rs,
			  rsbound_t rsbound,
			  rsparam_t rsparam,
			  const long fixed_i )
{
	 unsigned int np, nb, p, r, u, v, tmp, max_e;
	 int16_t subsgl[rsparam->len_p_rs], submul[rsparam->len_p_rs];
	 unsigned long fx_ul, gx_ul, pe = 1;
	 long  tmp2, start_j_idx[rsparam->len_p_rs], totnb, block_size;
	 float subf;
	 mpz_t tmpz, tmpu;
	 char flag[rsparam->len_p_rs];

	 block_size = L1_SIZE;
	 totnb = (rsbound->Bmax - rsbound->Bmin + 1) / block_size;

#if DEBUG_ROOTSIEVE_V
	 fprintf ( stderr,
			   "# Stat: totnb: %ld, block_size: %ld, total_size: %ld\n",
			   totnb, block_size, totnb*block_size );
#endif

	 mpz_init (tmpz);
	 mpz_init (tmpu);

#if DEBUG_ROOTSIEVE_V
	 int total = 0;
	 int st = 0;
	 int sts1 = 0, sts2 = 0;
	 int stm1 = 0, stm2 = 0;
	 int c = 0, cs = 0, cm = 0;
#endif

	 /* Init subsgl[] for all primes in root sieve */
	 for (np = 0; np < rsparam->len_p_rs; np ++) {
		  p = primes[np];
		  subf = (float) p * log ( (float) p) / ((float) p * (float) p - 1.0);
		  subsgl[np] = (int16_t) ceil (subf * 1000.0);
		  subf = log ( (float) p) / ( (float) p + 1.0);
		  submul[np] = (int16_t) ceil (subf * 1000.0);
		  start_j_idx[np] = 0;
		  flag[np] = 0;
	 }

	 /* Compute u by A + i*MOD = u */
	 ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, fixed_i, tmpu);

	 /* For each r < bound_prime */
	 for (r = 0; r < primes[rsparam->len_p_rs]; r++) {

		  /* u*g(r)^2 - f(r)g'(r) + f'(r)g(r) */
		  mpz_mul (tmpz, rs->gx[r], rs->gx[r]);
		  mpz_mul (tmpz, tmpz, tmpu);
		  mpz_sub (tmpz, tmpz, rs->numerator[r]);

#if DEBUG_ROOTSIEVE_V
		  st = cputime ();
#endif
		  /* For each block */
		  for (nb = 0; nb < totnb + 1; nb ++) {

			   /* For each r < p < Bound_p*/
			   for (np = next_prime_idx[r]; np < rsparam->len_p_rs; np ++) {

					p = primes[np];

					/* skip these */
					if (mpz_divisible_ui_p(rs->gx[r], p) != 0)
						 continue;

					/* e depends on p */
					max_e = log (200.0) / log ((double) p);
					pe = 1UL;
					for (tmp = 0; tmp < max_e; tmp ++)
						 pe = pe * p;

					/* The first block, compute starting points for the sieve */
					if (nb == 0) {

						 /* compute u (mod p) from u */
						 u = (unsigned int) mpz_fdiv_ui (tmpu, p);

						 /* use single precision */
						 fx_ul = mpz_fdiv_ui (rs->fx[r], p);
						 gx_ul = mpz_fdiv_ui (rs->gx[r], p);

						 /* compute v in f(r) + u*r*g(r) + v*g(r) = 0 (mod p) */
						 v = (unsigned int) compute_v_ul (fx_ul, gx_ul, r, u, p);

						 /* reset for the first block */
						 flag[np] = 0;

						 /* If MOD = 0 (mod p) in B + MOD*j = v (mod p), and B = v (mod p) */
						 if (mpz_divisible_ui_p (rsbound->MOD, p) != 0) {

							  if (mpz_fdiv_ui (rsbound->B, p) == v) {

								   /* don't sieve in this case, since all the elements on
									  the array is B + MOD*j = 0 (mod pe), are p-valuation equiv */
								   if (mpz_divisible_ui_p (rsbound->MOD, pe ) != 0)
										continue;
								   /* case where MOD has p-valuation smaller than e.
									  In this case, if r is multiple root, then do the resursion;
									  If r is simple, the all ele on the array are equivalent. So
									  we ignore sieving. */
								   else {
										if (mpz_divisible_ui_p(tmpz, p) == 0)
											 continue;
										/* the multiple root case will be caught by flag == 2 */
										else {
											 start_j_idx[np] = 0;
											 flag[np] = 2;
										}
								   }
							  }
							  else
								   continue;
						 }
						 else {
							  flag[np] = 1;
							  /* v -> j in B + MOD*j = v (mod p) */
							  tmp = (unsigned int) uv2ab_mod (rsbound->B, rsbound->MOD, v, p);
							  /* smallest tmp + p*i such that A0 < tmp + p*i, where A0 < 0,
								 hence i = ceil((A0-tmp)/p). Note, this should be negative. */
							  tmp2 = (long) ceil ( (double) (rsbound->Bmin -tmp) / (double) p)
								   * (long) p + (long) tmp;
							  start_j_idx[np] = ab2ij (rsbound->Bmin, tmp2);
							  /*
								if (r == 1) {
								fprintf (stderr, "tmp: %u, tmp2: %ld, rsbound->Bmin: %ld\n",
								tmp, tmp2, rsbound->Bmin);
								fprintf (stderr, " f(%u) + ( %u * %u + %u ) * g(%u) = 0 (mod %u)\n",
								r, u, r, v, r, p);
								fprintf (stderr, " start_j_idx: %ld\n", start_j_idx[np]);
								}
							  */
						 }
					}



#if DEBUG_ROOTSIEVE_V
					c++;
#endif

					/* cases where we don't want to sieve at all */
					if (flag[np] == 0)
						 continue;

					/* r is a simple root for current (u, v, p). Then r is simple root
					   for (u, v+i*p, p) for all i. */
					if (mpz_divisible_ui_p(tmpz, p) == 0)
					{

#if DEBUG_ROOTSIEVE_V
						 cs ++;
						 sts1 = cputime();
#endif

						 if (nb != totnb)
							  start_j_idx[np] = rootsieve_run_line ( ARRAY,
																	 block_size * (nb + 1), // largest index of the array
																	 start_j_idx[np], p, subsgl[np] );
						 else
							  start_j_idx[np] = rootsieve_run_line ( ARRAY,
																	 rsbound->Bmax - rsbound->Bmin, // largest index of the array
																	 start_j_idx[np], p, subsgl[np] );

#if DEBUG_ROOTSIEVE_V
						 sts2 += cputime() - sts1;
#endif
					}

					/* r is multiple root for current u, p */
					else {

						 /* don't sieve in block, assume cm << cs */
						 if (nb == 0) {
#if DEBUG_ROOTSIEVE_V
							  cm++;
							  stm1 = cputime();
#endif

							  /* Two cases:
								 If flag = 1, then MOD != 0 (mod p), everything is normal;
								 If MOD = 0 (mod p), then at nb = 0, start_j_idx is 0;
								 Also note, must use u (mod pe) instead of u (mod p) */
							  rootsieve_run_multroot ( ARRAY,
													   rs,
													   rsbound,
													   (unsigned int) mpz_fdiv_ui (tmpu, pe),
													   start_j_idx[np],
													   r,
													   p,
													   max_e,
													   submul[np] );
						 }
#if DEBUG_ROOTSIEVE_V
						 stm2 += cputime() - stm1;
#endif
					}
			   }
		  }

#if DEBUG_ROOTSIEVE_V
		  fprintf ( stderr,
					"# Stat: (r = %u,  p >= %u) takes %dms\n",
					r, primes[next_prime_idx[r]], cputime() - st );

		  total += cputime() - st;
#endif
	 }

	 mpz_clear (tmpz);
	 mpz_clear (tmpu);

#if DEBUG_ROOTSIEVE_V
	 fprintf ( stderr,
			   "# Stat: tot_sieve_time: %dms, tot_loops: %d, simple_loops: %d took %dms, mul_loops: %d took %dms\n",
			   total, c, cs, sts2, cm, stm2 );
#endif
}


/*
  init root sieve array with biased alpha.
*/
static void
rootsieve_array_init ( int16_t **A,
					   unsigned long len )
{
	 /* Init array A
		sage: sum ([p/(p^2-1)*log(p) for p in prime_range(200)])
		4.842766583050838  */

	 float tmpf = SUP_ALPHA;
	 int16_t tmpu = (int16_t) ceil (tmpf * 1000.0);
	 unsigned long j;
	 /* allocate matrix A. */
	 *A = (int16_t *) malloc ( len * sizeof (int16_t));

	 if ((*A) != NULL) {
		  for (j = 0; j < len; j++) {
			   (*A)[j] = tmpu;
		  }
	 }
	 else {
		  fprintf (stderr, "Error, cannot allocate memory in rootsieve_array_init().\n");
		  exit (1);
	 }
}


/*
  re-set root sieve array with biased alpha.
*/
static void
rootsieve_array_reset ( int16_t *A,
						unsigned long len )
{
	 float tmpf = SUP_ALPHA;
	 int16_t tmpu = (int16_t) ceil (tmpf * 1000.0);
	 unsigned long j;

	 if (A != NULL) {
		  for (j = 0; j < len; j++) {
			   A[j] = tmpu;
		  }
	 }
	 else {
		  fprintf (stderr, "Error, null memory in rootsieve_array_reset().\n");
		  exit (1);
	 }
}


#if 0
/*
  free array
*/
static void
rootsieve_array_free ( float ***A,
					   unsigned long ibound )
{
	 unsigned long i;
	 for (i = 0; i < ibound; i++)
		  free ((*A)[i]);
	 free (*A);
}
#endif


/*
  Init rsbound.
*/
static void
rsbound_init ( rsbound_t rsbound )
{
	 mpz_init_set_ui (rsbound->Umax, 0UL);
	 mpz_init_set_ui (rsbound->Umin, 0UL);
	 mpz_init_set_ui (rsbound->Vmax, 0UL);
	 mpz_init_set_ui (rsbound->Vmin, 0UL);
	 mpz_init_set_ui (rsbound->A, 0UL);
	 mpz_init_set_ui (rsbound->B, 0UL);
	 mpz_init_set_ui (rsbound->MOD, 0UL);

	 rsbound->Amax = rsbound->Amin = 0;
	 rsbound->Bmax = rsbound->Bmin = 0;
}

/*
  Free rsbound.
*/
static void
rsbound_free ( rsbound_t rsbound )
{
	 mpz_clear (rsbound->Umax);
	 mpz_clear (rsbound->Umin);
	 mpz_clear (rsbound->Vmax);
	 mpz_clear (rsbound->Vmin);
	 mpz_clear (rsbound->A);
	 mpz_clear (rsbound->B);
	 mpz_clear (rsbound->MOD);
}


/*
  Given rsparam->sizebound_ratio_rs and polynomial, set the size for
  sieving region. Set the Amax, Amin, Amax, Amin.
*/
void
rsbound_setup_AB_bound ( rsbound_t rsbound,
						 rsparam_t rsparam,
						 param_t param,
						 mpz_t mod,
						 int verbose )
{
	 /* verbosee == 1 means "tune mode";
		verbose == 0 means "polyselect2 mode"; */
	 if (verbose == 2) {

		  if (param->s2_Amax >= 0 && param->s2_Bmax > 0) {
			   rsbound->Amax = param->s2_Amax;
			   rsbound->Bmax = param->s2_Bmax;
		  }
		  else {
			   unsigned long len;
			   mpz_t q;
			   mpz_init (q);
			   mpz_fdiv_q (q, rsparam->global_v_bound_rs, mod);
			   len =  mpz_get_ui (q);
			   if (len > (SIEVEARRAY_SIZE << 3)) {
					rsbound->Bmax = (SIEVEARRAY_SIZE << 3);
			   }
			   else {
					rsbound->Bmax = (long) len;
			   }
			   mpz_set_ui (q, rsparam->global_u_bound_rs);
			   mpz_fdiv_q (q, q, mod);
			   len =  mpz_get_ui (q);
			   if (len > 8) {
					rsbound->Amax = 8;
			   }
			   else {
					rsbound->Amax = (long) len;
			   }
			   mpz_clear (q);
		  }
	 }
	 else if (verbose == 0) {

		  rsbound->Amax = 0;

		  unsigned long len;
		  mpz_t q;
		  mpz_init (q);
		  mpz_fdiv_q (q, rsparam->global_v_bound_rs, mod);
		  len =  mpz_get_ui (q);

		  if (len > (SIEVEARRAY_SIZE << 3)) {
			   rsbound->Bmax = (SIEVEARRAY_SIZE << 3);
		  }
		  else {
			   rsbound->Bmax = (long) len;
		  }

		  mpz_clear (q);
	 }
	 else {
		  rsbound->Amax = 0;
		  rsbound->Bmax = TUNE_SIEVEARRAY_SIZE;
	 }

	 rsbound->Bmin = -rsbound->Bmax;
	 rsbound->Amin = -rsbound->Amax;
}


/*
  Given Amax, Amin, Bmax, Bmin, set the Umax, Umin, Vmax, Vmin.
  Note that they should have similar size as rsparam->global_v_bound_rs.
*/
void
rsbound_setup_sublattice ( rsbound_t rsbound,
						   mpz_t sl_A,
						   mpz_t sl_B,
						   mpz_t mod )
{
	 mpz_set (rsbound->A, sl_A);
	 mpz_set (rsbound->B, sl_B);
	 /* the rsparam->modulus might be not true since rsparam might be
		changed for different qudratic rotations. Instead, the true mod
		is recorded in the priority queue */
	 mpz_set (rsbound->MOD, mod);

	 ab2uv (rsbound->A, rsbound->MOD, rsbound->Amax, rsbound->Umax);
	 ab2uv (rsbound->A, rsbound->MOD, rsbound->Amin, rsbound->Umin);
	 ab2uv (rsbound->B, rsbound->MOD, rsbound->Bmax, rsbound->Vmax);
	 ab2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, rsbound->Vmin);
}


/*
  Print root sieve bound and sublattice.
*/
static void
rsbound_print ( rsbound_t rsbound )
{
	 gmp_fprintf ( stderr,
				   "# Info: (u, v) = (%Zd + i * %Zd, %Zd + j * %Zd)\n",
				   rsbound->A, rsbound->MOD, rsbound->B, rsbound->MOD );

	 gmp_fprintf ( stderr,
				   "# Info: (Amin: %4ld, Amax: %4ld) -> (Umin: %6Zd, Umax: %6Zd)\n",
				   rsbound->Amin, rsbound->Amax,
				   rsbound->Umin, rsbound->Umax );

	 gmp_fprintf ( stderr,
				   "# Info: (Bmin: %4ld, Bmax: %4ld) -> (Vmin: %6Zd, Vmax: %6Zd)\n",
				   rsbound->Bmin, rsbound->Bmax,
				   rsbound->Vmin, rsbound->Vmax );
}


/*
  Init rsstr_t with.
*/
void
rsstr_init ( rsstr_t rs )
{
	 unsigned int i;
	 /* init */
	 mpz_init (rs->n);
	 mpz_init (rs->m);
	 rs->f = (mpz_t*) malloc ((MAX_DEGREE + 1) * sizeof (mpz_t));
	 rs->g = (mpz_t*) malloc ((MAX_DEGREE + 1) * sizeof (mpz_t));
	 /* pre-computing function values f(r), g(r) for 0 <= r < p*/
	 (rs->fx) = (mpz_t *) malloc ( (primes[NP-1]+1) * sizeof (mpz_t) );
	 (rs->gx) = (mpz_t *) malloc ( (primes[NP-1]+1) * sizeof (mpz_t) );
	 (rs->numerator) = (mpz_t *) malloc ( (primes[NP-1]+1) * sizeof (mpz_t) );

	 if ((rs->f == NULL) || (rs->g == NULL)) {
		  fprintf (stderr, "Error, cannot allocate memory for polynomial coefficients in rsstr_init().\n");
		  exit (1);
	 }

	 if (((rs->fx) == NULL) || ((rs->gx) == NULL) || ((rs->numerator) == NULL)) {
		  fprintf (stderr, "Error, cannot allocate memory for polynomials values in rsstr_init().\n");
		  exit (1);
	 }

	 for (i = 0; i <= MAX_DEGREE; i++)
	 {
		  mpz_init (rs->f[i]);
		  mpz_init (rs->g[i]);
	 }
	 for (i = 0; i <= primes[NP-1]; i++)
	 {
		  mpz_init (rs->fx[i]);
		  mpz_init (rs->gx[i]);
		  mpz_init (rs->numerator[i]);
	 }
}


/*
  Precompute fx, gx and numerator in rsstr_t.
  Note, rs->f, etc must be set already.
*/
static void
rsstr_setup ( rsstr_t rs )
{
	 int i;
	 mpz_t t;
	 mpz_init (t);

	 /* polynomial degree */
	 for ((rs->d) = MAX_DEGREE; (rs->d) > 0 && mpz_cmp_ui ((rs->f[rs->d]), 0) == 0; rs->d --);

	 /* M = -Y0/Y1 mod N */
	 mpz_invert (rs->m, rs->g[1], rs->n);
	 mpz_neg (rs->m, rs->m);
	 mpz_mul (rs->m, rs->m, rs->g[0]);
	 mpz_mod (rs->m, rs->m, rs->n);

	 /* check M ? a root of the algebraic polynomial mod N */
	 mpz_set (t, rs->f[rs->d]);
	 for (i = rs->d - 1; i >= 0; i --) {
		  mpz_mul (t, t, rs->m);
		  mpz_add (t, t, rs->f[i]);
		  mpz_mod (t, t, rs->n);
	 }
	 if (mpz_cmp_ui (t, 0) != 0)
	 {
		  fprintf (stderr, "ERROR: The following polynomial have no common root. \n");
		  print_poly2 (rs->f, rs->g, rs->d, rs->n);
		  exit (1);
	 }

	 /* save f(l) to an array for all l < B
		can we save the differences? not necessary since anyway
		precomputation is only done for once for one polynomial. */
	 eval_polys (rs->f, rs->g, rs->fx, rs->gx, rs->numerator, primes, rs->d);

	 rs->alpha_proj = get_biased_alpha_projective (rs->f, rs->d, 2000);

	 mpz_clear (t);
}


/*
  Free rsstr_t.
*/
static void
rsstr_free ( rsstr_t rs )
{
	 unsigned int i;
	 /* free fl, gl */
	 for (i = 0; i <= primes[NP-1]; i ++)
	 {
		  mpz_clear(rs->fx[i]);
		  mpz_clear(rs->gx[i]);
		  mpz_clear(rs->numerator[i]);
	 }
	 for (i = 0; i <= MAX_DEGREE; i++)
	 {
		  mpz_clear (rs->f[i]);
		  mpz_clear (rs->g[i]);
	 }
	 mpz_clear (rs->n);
	 mpz_clear (rs->m);
	 free (rs->fx);
	 free (rs->gx);
	 free (rs->numerator);
	 free (rs->f);
	 free (rs->g);
}


/*
  Init root sieve parameters. Often, there is no need to change
  them. The real customisation happens in rsparam_setup() function.
*/
static void
rsparam_init ( rsparam_t rsparam )
{
	 /* will be set in rsparam_setup */
	 rsparam->nbest_sl = 0;
	 rsparam->ncrts_sl = 0;
	 rsparam->len_e_sl = 0;
	 rsparam->tlen_e_sl = 0;
	 rsparam->sizebound_ratio_rs = 0;
	 rsparam->exp_min_alpha_rs = 0.0;
	 rsparam->global_w_bound_rs = 0;
	 rsparam->global_u_bound_rs = 0;
	 mpz_init (rsparam->modulus);
	 mpz_init_set_ui (rsparam->global_v_bound_rs, 0UL);

	 /* number of primes beside e_sl[] considered in second stage
		root sieve. Larger takes longer time, but more accurate. */
	 rsparam->len_p_rs = NP - 1;
	 if (rsparam->len_p_rs >= NP)
		  rsparam->len_p_rs = NP - 1;
}


/*
  init the best poly (for the return)
*/
static void
bestpoly_init ( bestpoly_t bestpoly,
				int d )
{
	 int i;
	 bestpoly->f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
	 bestpoly->g = (mpz_t*) malloc (2 * sizeof (mpz_t));
	 if ((bestpoly->f == NULL) || (bestpoly->g == NULL)) {
		  fprintf (stderr, "Error, cannot allocate memory for bestpoly_init()\n");
		  exit (1);
	 }
	 for (i = 0; i <= d; i++)
	 {
		  mpz_init (bestpoly->f[i]);
	 }

	 mpz_init (bestpoly->g[0]);
	 mpz_init (bestpoly->g[1]);
}


/*
  free the best poly (for the return)
*/
static void
bestpoly_free ( bestpoly_t bestpoly,
				int d )
{
	 int i;
	 for (i = 0; i <= d; i++)
	 {
		  mpz_clear (bestpoly->f[i]);
	 }

	 mpz_clear (bestpoly->g[0]);
	 mpz_clear (bestpoly->g[1]);
	 free (bestpoly->f);
	 free (bestpoly->g);
}

/*
  replace f + k0 * x^t * (b*x - m) by f + k * x^t * (b*x - m), and return k to k0
  (modified from auxiliary.c)
*/
static void
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


/*
  Modifed from auxiliary.c using Emmanuel Thome and Paul Zimmermann's ideas.
  Assume lognorm(f + v*g) + E(alpha(f + v*g)) is first decreasing, then
  increasing, then the optimal v corresponds to the minimum of that function.

  Return best V and best E
*/
static double
rotate_bounds_V_mpz ( mpz_t *f,
					  int d,
					  double ratio_margin,
					  mpz_t b,
					  mpz_t m,
					  mpz_t V,
					  int method,
					  const int i_bias_u )
{
	 int i;
	 double lognorm, alpha, E, best_E;
	 mpz_t tmpv, best_V;

	 mpz_init (best_V);
	 mpz_init_set_ui (tmpv, 0);

	 best_E = L2_lognorm (f, d,
						  L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, method),
						  method);

	 /* look for positive k: 2, 4, 8, ... */
	 mpz_set_si (V, 1);
	 for (i = 0; i < 150 - i_bias_u ; i++, mpz_mul_si (V, V, 2) )
	 {
		  rotate_aux_mpz (f, b, m, tmpv, V, 0);
		  lognorm = L2_lognorm (f, d,
								L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, method),
								method);
		  alpha = exp_alpha[i + i_bias_u];
		  E = lognorm + alpha;

		  if (DEBUG)
			   gmp_fprintf (stderr, "# DEBUG  (V_bound): [%d-th] V: %15Zd, E: %.3f, lognorm: %.3f, alpha: %.2f\n",
							i, V, E, lognorm, alpha);

		  if (E < best_E * ratio_margin)
		  {
			   if (E < best_E)
					best_E = E;
		  }
		  else
			   break;
	 }

	 mpz_set (best_V, V);
	 mpz_set_ui (V, 0);
	 /* go back to k=0 */
	 rotate_aux_mpz (f, b, m, tmpv, V, 0);

	 /* look for negative k: -2, -4, -8, ... */
	 mpz_set_si (V, -1);
	 best_E = L2_lognorm (f, d,
						  L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, method),
						  method);
	 for (i = 0; i < 150 - i_bias_u; i++, mpz_mul_si (V, V, 2))
	 {
		  rotate_aux_mpz (f, b, m, tmpv, V, 0);
		  lognorm = L2_lognorm (f, d,
								L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, method),
								method);
		  alpha = exp_alpha[i + i_bias_u];
		  E = lognorm + alpha;

		  if (DEBUG)
			   gmp_fprintf (stderr, "# DEBUG  (V_bound): [%d-th] V: %15Zd, E: %.3f, lognorm: %.3f, alpha: %.2f\n",
							i, V, E, lognorm, alpha);

		  if (E < best_E * ratio_margin)
		  {
			   if (E < best_E) {
					best_E = E;
					mpz_set (best_V, V);
			   }
		  }
		  else
			   break;
	 }

	 mpz_set (V, best_V);
	 mpz_set_ui (best_V, 0);
	 /* go back to k=0 */
	 rotate_aux_mpz (f, b, m, tmpv, best_V, 0);

	 mpz_clear (best_V);
	 mpz_clear (tmpv);
	 return best_E;
}


/* find bound w for qudratic rotation */
static double
rotate_bounds_W ( mpz_t *f,
				  int d,
				  mpz_t b,
				  mpz_t m,
				  unsigned long *W,
				  int rotation_degree,
				  int method,
				  int verbose )
{
	 int i, upper_bound = 12;
	 double skewness, init_lognorm, lognorm;
	 long w0 = 0, w;

	 skewness = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, method);
	 init_lognorm = L2_lognorm (f, d, skewness, method);

	 /* look for positive w: 1, 2, 4, 8, ... */
	 w = 1;

	 /* decide the bound for the loop */
	 if (rotation_degree == 1)
		  upper_bound = 48;

	 for (i = 0; i < upper_bound; i++, w *= 2)
	 {
		  /* rotate by w*x, and fix this polynomial in this loop. */
		  w0 = rotate_aux (f, b, m, w0, w, rotation_degree);

		  skewness = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, method);
		  lognorm = L2_lognorm (f, d, skewness, method);

		  //fprintf (stderr, "# DEBUG --- [%d-th] W: %ld, lognorm: %f, bound_norm: %f\n", i, w, lognorm, init_lognorm * 1.1);

		  if (lognorm > init_lognorm * 1.1)
			   break;
	 }

	 (*W) = (unsigned long) w;

	 if (verbose == 2) {
		  if (rotation_degree == 2)
			   fprintf (stderr, "# Info: w upper bound: %lu, norm bound: %.2f\n", *W, init_lognorm * 1.1); // TBC
		  else if (rotation_degree == 1)
			   fprintf (stderr, "# Info: u upper bound: %lu, norm bound: %.2f\n", *W, init_lognorm * 1.1); // TBC
	 }
	 /* go back to w=0 */
	 rotate_aux (f, b, m, w0, 0, rotation_degree);

	 return init_lognorm * 1.1;
}


/*
  For each U bound, identify the best V bound (such that E is
  smallest for this fixed U). Then return the best (U, V) pair.

  Experimentally, this is better than considering U and U*skew
  but will be slower.
*/
static double
rotate_bounds_UV ( mpz_t *f,
				   int d,
				   double ratio_margin,
				   mpz_t b,
				   mpz_t m,
				   unsigned long *U,
				   mpz_t V,
				   int method,
				   int verbose )
{
	 int i;
	 long k0 = 0, tmpU;
	 double skewness, init_lognorm, E, best_E;
	 mpz_t best_V;
	 unsigned long upper_bound_U = 0;

	 mpz_init (best_V);
	 mpz_set_ui (best_V, 0UL);

	 skewness = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, method);
	 init_lognorm = L2_lognorm (f, d, skewness, method);
	 best_E = init_lognorm;

	 /* First, compute an upper bound for linear rotation U. */
	 rotate_bounds_W ( f,
					   d,
					   b,
					   m,
					   &upper_bound_U,
					   1,
					   DEFAULT_L2_METHOD,
					   verbose );
	 upper_bound_U = (unsigned long) (log ((double) upper_bound_U) / log(2.0));

	 /* Then, look for best (U, V) combinations, where positive k: 2, 4, 8, ... */
	 tmpU = 1;
	 for (i = 0; i < (int) upper_bound_U; i++, tmpU *= 2)
	 {
		  /* rotate by u*x, and fix this polynomial in this loop. */
		  k0 = rotate_aux (f, b, m, k0, tmpU, 1);

		  if (DEBUG)
			   fprintf (stderr, "# DEBUG --- [%d-th] U: %ld ---\n", i, tmpU);

		  /* identify best v in rotating by v */
		  E = rotate_bounds_V_mpz (f, d, ratio_margin,
		  						   b, m, V, DEFAULT_L2_METHOD, i);

		  if (DEBUG)
			   gmp_fprintf (stderr, "# DEBUG: (U: %ld, V: %Zd), best_E: %.3f\n",
							tmpU, V, E);

		  if (E < best_E)
		  {
			   best_E = E;
			   // *U = abs(tmpU); problem for tmpU >= 2^31
			   (*U) = (unsigned long) (tmpU >= 0) ? tmpU : -tmpU;
			   mpz_set (best_V, V);
		  }
	 }

	 /* go back to k=0 */
	 rotate_aux (f, b, m, k0, 0, 1);

	 mpz_set (V, best_V);
	 mpz_clear (best_V);

	 return best_E;
}


/*
  Given rsparam->sizebound_ratio_rs and polynomial information,
  compute rsparam->global_u_bound_rs and rsparam->global_v_bound_rs;
  Then it will set e_sl[];
*/
static void
rsparam_setup ( rsparam_t rsparam,
				rsstr_t rs,
				int verbose )
{

	 /* 1. "rsparam->sizebound_ratio_rs"
		the higher, the more margin in computing the sieving bound
		u and v, hence the larger the sieving bound, and hence
		larger e_sl[] in rsparam_setup(). */
	 rsparam->sizebound_ratio_rs = 1.01;

	 double best_E, lognorm_bound;
	 mpz_t b, m;
	 mpz_init (b);
	 mpz_init (m);
	 mpz_set (b, rs->g[1]);
	 mpz_set (m, rs->g[0]);
	 mpz_neg (m, m);

	 /* 2. "global_w_bound_rs" */
	 if (rs->d == 6) {
		  lognorm_bound = rotate_bounds_W ( rs->f,
											rs->d,
											b,
											m,
											&(rsparam->global_w_bound_rs),
											2,
											DEFAULT_L2_METHOD,
											verbose );
	 }
	 else
		  rsparam->global_w_bound_rs = 0;

	 /* 3. "global_u_bound_rs" and "global_v_bound_rs"
		-- global_u_bound will be used to identify good sublattice and decide e[],
		-- global_v_bound will be used to identify sieving bound */
	 best_E = rotate_bounds_UV ( rs->f,
								 rs->d,
								 rsparam->sizebound_ratio_rs,
								 b,
								 m,
								 &(rsparam->global_u_bound_rs),
								 rsparam->global_v_bound_rs, /* u, v */
								 DEFAULT_L2_METHOD,
								 verbose );
	 mpz_clear (b);
	 mpz_clear (m);

	 /* 4. "rsparam->exp_min_alpha_rs" */
	 int size;
	 size = mpz_sizeinbase (rsparam->global_v_bound_rs, 2);
	 size += (int) (log ( (double) rsparam->global_u_bound_rs ) * 1.442695);
	 rsparam->exp_min_alpha_rs = exp_alpha[size-1];

	 if (verbose == 2)
		  gmp_fprintf ( stderr,
						"# Info: best (u, v) bound (%ld, %Zd) gives exp_best_E: %.3f, exp_min_alpha: %.3f\n",
						rsparam->global_u_bound_rs,
						rsparam->global_v_bound_rs,
						best_E,
						exp_alpha[size-1] );

	 /* 5. "rsparam->nbest_sl" and "rsparam->ncrts_sl" */
	 if (rs->d == 6) {
		  rsparam->nbest_sl = 128;
		  rsparam->ncrts_sl = 64; // for deg 6, there could be too much individual sublattices to do crts. We restrict the num.
	 }
	 else {
		  rsparam->nbest_sl = 128;
	 }

	 /* 6. "rsparam->len_e_sl" and "rsparam->e_sl" */
	 rsparam->len_e_sl = LEN_SUBLATTICE_PRIMES;
	 rsparam->tlen_e_sl = rsparam->len_e_sl;

	 if (rsparam->len_p_rs < rsparam->len_e_sl) {
		  fprintf (stderr, "# Warning: number of primes considered in the root sieve is smaller than that in fin_sublattice(). This might not be accurate. \n");
	 }

	 rsparam->e_sl = (unsigned short*)
		  malloc ( rsparam->len_e_sl * sizeof (unsigned short) );

	 if (rsparam->e_sl == NULL) {
		  fprintf (stderr, "Error, cannot allocate memory in rsparam_setup().\n");
		  exit (1);
	 }

	 for (size = 0; size < rsparam->len_e_sl; size ++) {
		  rsparam->e_sl[size] = 0;
	 }

	 /* Experimental */
	 if (rsparam->global_u_bound_rs <= 1024) {
		  rsparam->e_sl[0] = 2;
		  rsparam->e_sl[1] = 1;
		  rsparam->e_sl[2] = 0;
		  rsparam->e_sl[3] = 0;
	 }
	 else if (rsparam->global_u_bound_rs <= 2048) {
		  rsparam->e_sl[0] = 3;
		  rsparam->e_sl[1] = 1;
		  rsparam->e_sl[2] = 0;
		  rsparam->e_sl[3] = 0;
	 }
	 else if (rsparam->global_u_bound_rs <= 4096) {
		  rsparam->e_sl[0] = 4;
		  rsparam->e_sl[1] = 1;
		  rsparam->e_sl[2] = 0;
		  rsparam->e_sl[3] = 0;
	 }
	 else if (rsparam->global_u_bound_rs <= 8192) {
		  rsparam->e_sl[0] = 4;
		  rsparam->e_sl[1] = 2;
		  rsparam->e_sl[2] = 0;
		  rsparam->e_sl[3] = 0;
	 }
	 else {
		  rsparam->e_sl[0] = 4;
		  rsparam->e_sl[1] = 2;
		  rsparam->e_sl[2] = 1;
		  rsparam->e_sl[3] = 0;
	 }
}


/*
  For existing rsparam and when rs is changed, we update
  the new rotate-bounds.
*/
static void
rsparam_reset ( rsparam_t rsparam,
				rsstr_t rs,
				int verbose )
{
	 double best_E;
	 mpz_t b, m;
	 mpz_init (b);
	 mpz_init (m);
	 mpz_set (b, rs->g[1]);
	 mpz_set (m, rs->g[0]);
	 mpz_neg (m, m);

	 best_E = rotate_bounds_UV ( rs->f,
								 rs->d,
								 rsparam->sizebound_ratio_rs,
								 b,
								 m,
								 &(rsparam->global_u_bound_rs),
								 rsparam->global_v_bound_rs, /* u, v */
								 DEFAULT_L2_METHOD,
								 verbose );
	 mpz_clear (b);
	 mpz_clear (m);

	 int size;
	 size = mpz_sizeinbase (rsparam->global_v_bound_rs, 2);
	 size += (int) (log ( (double) rsparam->global_u_bound_rs ) * 1.442695);
	 rsparam->exp_min_alpha_rs = exp_alpha[size-1];

	 if (verbose == 2)
		  gmp_fprintf ( stderr,
						"# Info: (u, v) bound (%ld, %Zd) gives exp_best_E: %.3f, exp_min_alpha: %.3f\n",
						rsparam->global_u_bound_rs,
						rsparam->global_v_bound_rs,
						best_E,
						exp_alpha[size-1] );

	 /* Reset these, since anyway, there is a tune procedure after rsparam_reset being called. */
	 for (size = 0; size < rsparam->len_e_sl; size ++) {
		  rsparam->e_sl[size] = 0;
	 }

	 if (rsparam->global_u_bound_rs <= 1024) {
		  rsparam->e_sl[0] = 2;
		  rsparam->e_sl[1] = 1;
		  rsparam->e_sl[2] = 0;
		  rsparam->e_sl[3] = 0;
	 }
	 else if (rsparam->global_u_bound_rs <= 2048) {
		  rsparam->e_sl[0] = 3;
		  rsparam->e_sl[1] = 1;
		  rsparam->e_sl[2] = 0;
		  rsparam->e_sl[3] = 0;
	 }
	 else if (rsparam->global_u_bound_rs <= 4096) {
		  rsparam->e_sl[0] = 4;
		  rsparam->e_sl[1] = 1;
		  rsparam->e_sl[2] = 0;
		  rsparam->e_sl[3] = 0;
	 }
	 else if (rsparam->global_u_bound_rs <= 8192) {
		  rsparam->e_sl[0] = 4;
		  rsparam->e_sl[1] = 2;
		  rsparam->e_sl[2] = 0;
		  rsparam->e_sl[3] = 0;
	 }
	 else {
		  rsparam->e_sl[0] = 4;
		  rsparam->e_sl[1] = 2;
		  rsparam->e_sl[2] = 1;
		  rsparam->e_sl[3] = 0;
	 }
}


/*
  Free root sieve parameters.
*/
static void
rsparam_free ( rsparam_t rsparam )
{
	 free(rsparam->e_sl);
	 mpz_clear (rsparam->global_v_bound_rs);
	 mpz_clear (rsparam->modulus);
}


/*
  Init param
*/
void
param_init ( param_t param )
{
	 mpz_init (param->s2_u);
	 mpz_init (param->s2_v);
	 mpz_init (param->s2_mod);
	 mpz_set_ui (param->s2_u, 0);
	 mpz_set_ui (param->s2_v, 0);
	 mpz_set_ui (param->s2_mod, 0);
	 param->w_left_bound = 0;
	 param->w_length = 0;
	 param->s1_num_e_sl = 0;
	 param->s2_Amax = 0;
	 param->s2_Bmax = 0;
	 param->s2_w = 0;
	 param->s1_e_sl = (unsigned short*)
		  malloc ( LEN_SUBLATTICE_PRIMES * sizeof (unsigned short) );
	 int i;
	 for (i = 0; i < LEN_SUBLATTICE_PRIMES; i ++)
		  param->s1_e_sl[i] = 0;
}


void
param_clear ( param_t param )
{
	 mpz_clear (param->s2_u);
	 mpz_clear (param->s2_v);
	 mpz_clear (param->s2_mod);
	 free (param->s1_e_sl);
}


/*
  Rootsieve for f + (u*x +v)*g for each sublattice
*/
static double
rootsieve_uv ( rsstr_t rs,
			   rsbound_t rsbound,
			   rsparam_t rsparam,
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
	 long i, block_size, tmpBmax, tmpBmin;
	 unsigned long j;
	 double MurphyE;
	 mpz_t tmpv, tmpu;

	 mpz_init (tmpv);
	 mpz_init (tmpu);

	 /* test sieve in tune mode ? */
	 if (verbose == 1)
		  block_size = TUNE_SIEVEARRAY_SIZE;
	 else {
		  block_size = SIEVEARRAY_SIZE;

		  /* one block is all ready too long */
		  if ( (rsbound->Bmax - rsbound->Bmin + 1) < block_size ) {
			   block_size = rsbound->Bmax - rsbound->Bmin + 1;
			   }
	 }
#if DEBUG
	 fprintf (stderr, "# Info: totalnb: %lu, %lu, %lu\n", (rsbound->Bmax - rsbound->Bmin + 1), block_size,
			  (rsbound->Bmax - rsbound->Bmin + 1) / block_size);
#endif

	 /* E priority queue for each sublattice */
	 MurphyE_pq *E_pqueue;
	 new_MurphyE_pq (&E_pqueue, TOPE);
	 int st = cputime ();

	 int16_t *MAT;
	 rootscore_pq *alpha_pqueue;

	 /* init the sieve array and priority queue */
	 rootsieve_array_init (&MAT, block_size + 1);
	 new_rootscore_pq (&alpha_pqueue, TOPALPHA);

	 /* for each i -> each u = A + MOD * i */
	 tmpBmax = rsbound->Bmax;
	 tmpBmin = rsbound->Bmin;
	 for (i = 0; i < (rsbound->Amax - rsbound->Amin + 1); i ++)
	 {
		  long k = 0;

		  do {
			   k ++;
			   /* for each block of size = |L2| */
			   rsbound->Bmax = rsbound->Bmin + block_size;
#if DEBUG
			   fprintf (stderr, "[%ld, %ld] ", rsbound->Bmin, rsbound->Bmax);
			   fprintf (stderr, "NEXT i: %lu\n", i);
#endif

			   /* root sieve for v. Note, rsbound with changed Bmin and Bmax will be used in rootsieve_v(). */
			   rootsieve_v (MAT, rs, rsbound, rsparam, i);

			   /* record top n slots */
			   for (j = 0; j < (unsigned long) (rsbound->Bmax - rsbound->Bmin + 1); j++) {
					insert_rootscore_pq ( alpha_pqueue, i, j, MAT[j] );

					/* if (MAT[j] < -2000) { */
					/* 	 ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, i, tmpu); */
					/* 	 ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, j, tmpv); */
					/* 	 gmp_fprintf (stderr, "MAT[]: %d, u: %Zd, v: %Zd, i: %lu, j: %lu\n", MAT[j], tmpu, tmpv, i, j); */
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

			   		MurphyE = print_poly_info (fuv, guv, rs->d, rs->n, rs->m, 0);
					insert_MurphyE_pq ( E_pqueue, w, tmpu, tmpv, MurphyE );
			   }

			   /* next j */
			   rsbound->Bmin = rsbound->Bmax + 1;

			   /* reset sieve array and priority queue */
			   rootsieve_array_reset (MAT, block_size + 1);
			   reset_rootscore_pq ( alpha_pqueue );

		  } while (rsbound->Bmin < tmpBmax);

		  /* next i */
		  rsbound->Bmax = tmpBmax;
		  rsbound->Bmin = tmpBmin;

		  /* reset sieve array and priority queue */
		  rootsieve_array_reset (MAT, block_size + 1);
		  reset_rootscore_pq ( alpha_pqueue );

	 } // FINISHING U-ROTATION

	 /* free priority queue and sieving array */
	 free_rootscore_pq (&alpha_pqueue);
	 free (MAT);

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
		  MurphyE += print_poly_info (fuv, guv, rs->d, rs->n, rs->m, verbose);
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
static double
rootsieve_main_stage2_run ( rsstr_t rs,
							rsbound_t rsbound,
							rsparam_t rsparam,
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
		  fprintf (stderr, "Error, cannot allocate memory in rootsieve_main_stage2_run().\n");
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
static double
rootsieve_main_stage2_prepare ( rsstr_t rs,
								rsparam_t rsparam,
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

	 /* root sieve on this sublattice */
	 rsbound_setup_sublattice ( rsbound,
								u,
								v,
								mod );

	 /* compute exact sieving bounds UV given size AB depending
		on current A, B, MOD */
	 if (verbose == 2)
		  rsbound_print (rsbound);

	 /* root sieve */
	 ave_MurphyE = rootsieve_main_stage2_run ( rs,
											   rsbound,
											   rsparam,
											   global_E_pqueue,
											   w,
											   verbose );
	 /* free rsbound */
	 rsbound_free (rsbound);

	 return ave_MurphyE;
}


/*
  Stage 1: record good sublattices in "alpha_pqueue".
*/
static int
rootsieve_main_stage1 ( rsstr_t rs,
						rsparam_t rsparam,
						sub_alpha_pq *alpha_pqueue,
						int verbose,
						int w )
{
	 int st, i, re;
	 mpz_t *fuv;
	 double alpha_lat, skew, logmu;
	 sublattice_pq *pqueue;

	 /* fuv is f+(u*x+v)*g */
	 fuv = (mpz_t*) malloc ((rs->d + 1) * sizeof (mpz_t));
	 if (fuv == NULL) {
		  fprintf (stderr, "Error, cannot allocate memory in rootsieve_main_stage1_run().\n");
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
		  gmp_fprintf ( stderr, "# Info: find best sublattices over (Mod %Zd) took %dms\n",
						rsparam->modulus, cputime () - st );

	 /* filter all sublattices into another, global queue ranked by the
		sublattice alphas. */
	 for (i = 1; i < pqueue->used; i ++) {

		  compute_fuv_mp (fuv, rs->f, rs->g, rs->d, pqueue->u[i], pqueue->v[i]);
		  //alpha_lat = get_biased_alpha_affine (fuv, rs->d, primes[rsparam->tlen_e_sl - 1]);
		  alpha_lat = get_alpha (fuv, rs->d, 2000);
		  skew = L2_skewness (fuv, rs->d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
		  logmu = L2_lognorm (fuv, rs->d, skew, DEFAULT_L2_METHOD);

#if DEBUG
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
								rsparam->modulus,
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


/*
  auxiliary for tune
*/
static double
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

	 re = rootsieve_main_stage1 ( rs,
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
		  				rsparam->exp_min_alpha_rs );
		  */

		  ave_MurphyE = rootsieve_main_stage2_prepare ( rs,
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
static double
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


/* Record best poly */
static void
rootsieve_main_run_bestpoly ( rsstr_t rs,
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
		  fprintf (stderr, "Error, cannot allocate memory in rootsieve_main_run().\n");
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
		  ave_MurphyE = print_poly_info (fuv, guv, rs->d, rs->n, rs->m, 0); // only output when verbose == 2.

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
static double
rootsieve_main_run ( rsstr_t rs,
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
	 rsparam_init (rsparam);
	 rsparam_setup (rsparam, rs, 0);
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
			   fprintf (stderr,
						"# Info: quadratic rotation by %d*x^2\n", i);
		  old_i = rotate_aux (rs->f, rs->g[1], m, old_i, i, 2);
		  rsstr_setup (rs);

		  /* either use input e_sl[] or tune only once */
		  if (re == 0) {
			   if (param->flag == 1) {
					rsparam_reset (rsparam, rs, verbose);
					for (k = 0; k < LEN_SUBLATTICE_PRIMES; k ++)
						 rsparam->e_sl[k] = param->s1_e_sl[k];
			   }
			   else {
					rsparam_reset (rsparam, rs, verbose);
					rsparam_tune (rs, rsparam, param, 10, i, verbose);
			   }
		  }

		  re ++;
		  k = rootsieve_main_stage1 ( rs,
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
			   gmp_fprintf ( stderr, "# Info: #%4d sublattice (w, u, v): (%d, %Zd, %Zd) (mod %Zd), alpha: %.2f\n",
							 i + 1,
							 w[i],
							 u[i],
							 v[i],
							 mod[i],
							 sub_alpha[i] );
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
		  rsparam_reset (rsparam, rs, 1);

		  //print_poly_info (rs->f, rs->g, rs->d, rs->n, rs->m, 1);

		  ave_MurphyE = rootsieve_main_stage2_prepare ( rs,
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
		  rootsieve_main_run_bestpoly ( rs,
										global_E_pqueue,
										bestpoly );
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
static double
rootsieve_main_run_stage2only ( rsstr_t rs,
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
		  fprintf (stderr, "Error, cannot allocate memory in rootsieve_main_run_stage2only().\n");
		  exit (1);
	 }
	 for (i = 0; i <= rs->d; i++)
		  mpz_init_set (fuv[i], rs->f[i]);
	 for (i = 0; i < 2; i++)
		  mpz_init_set (guv[i], rs->g[i]);

	 rsparam_init (rsparam);
	 rsparam_setup (rsparam, rs, 2);

	 /* rotate polynomial by f + rot*x^2 */
	 old_i = 0;
	 old_i = rotate_aux (rs->f, rs->g[1], m, old_i, param->s2_w, 2);
	 rsstr_setup (rs);
	 rsparam_reset (rsparam, rs, 1);

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

     //print_poly_info (rs->f, rs->g, rs->d, rs->n, rs->m, 2);
	 MurphyE_pq *global_E_pqueue;
	 new_MurphyE_pq (&global_E_pqueue, 4);

	 ave_MurphyE = rootsieve_main_stage2_prepare ( rs,
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
static void
rootsieve_main ( rsstr_t rs,
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
		  rootsieve_main_run_stage2only (rs, param);
	 }
	 else {
		  if (rs->d == 5) {
			   /* not qudratci rot for deg 5 polynomial */
			   param->w_left_bound = 0;
			   param->w_length = 1;
			   rootsieve_main_run (rs, bestpoly, param, verbose);
		  }
		  else if (rs->d == 6) {
			   rootsieve_main_run (rs, bestpoly, param, verbose);
		  }
		  else {
			   fprintf (stderr, "Error: only support deg 5 or 6.\n");
			   exit(1);
		  }
	 }
}


/*
  Do the root sieve on all polynomials in the file.
*/
void
rootsieve_polyselect ( mpz_t *f,
					   int d,
					   mpz_t m,
					   mpz_t l,
					   mpz_t N,
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
		  param->w_left_bound = -32;
		  param->w_length = 64;
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
	 /* start main rootsieve function */
	 rootsieve_main (rs, bestpoly, param, verbose);

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

/*
  Do the root sieve on all polynomials in the file.
*/
void
rootsieve_file ( FILE *file,
				 param_t param )
{
	 unsigned flag = 0UL, count = 0;
	 char str[MAX_LINE_LENGTH];

	 /* rootsieve_struct */
	 rsstr_t rs;
	 rsstr_init (rs);

	 /* for each polynomial, do the root sieve. */
	 while (1) {

		  /* read poly */
		  if (fgets(str, MAX_LINE_LENGTH, file) == NULL)
			   break;

		  if ( str[0] == 'Y' ) {
			   if ( str[1] == '1' ) {
					gmp_sscanf (str, "Y1: %Zd\n", rs->g[1]);
					(flag) ^= (1<<8);
			   }
			   else if ( str[1] == '0' ) {
					gmp_sscanf (str, "Y0: %Zd\n", rs->g[0]);
					(flag) ^= (1<<7);
			   }
			   else {
					fprintf (stderr, "Error in parsing line %s.\n", str);
					exit (1);
			   }
		  }
		  else if ( str[0] == 'c') {
			   if ( str[1] == '0' ) {
					gmp_sscanf (str, "c0: %Zd\n", rs->f[0]);
					(flag) ^= 1;
			   }
			   else if ( str[1] == '1' ) {
					gmp_sscanf (str, "c1: %Zd\n", rs->f[1]);
					(flag) ^= (1<<1);
			   }
			   else if ( str[1] == '2' ) {
					gmp_sscanf (str, "c2: %Zd\n", rs->f[2]);
					(flag) ^= (1<<2);
			   }
			   else if ( str[1] == '3' ) {
					gmp_sscanf (str, "c3: %Zd\n", rs->f[3]);
					(flag) ^= (1<<3);
			   }
			   else if ( str[1] == '4' ) {
					gmp_sscanf (str, "c4: %Zd\n", rs->f[4]);
					(flag) ^= (1<<4);
			   }
			   else if ( str[1] == '5' ) {
					gmp_sscanf (str, "c5: %Zd\n", rs->f[5]);
					(flag) ^= (1<<5);
			   }
			   else if ( str[1] == '6' ) {
					gmp_sscanf (str, "c6: %Zd\n", rs->f[6]);
					(flag) ^= (1<<6);
			   }


			   else
			   {
					fprintf (stderr, "Error in parsing line %s.\n", str);
					exit (1);
			   }
		  }
		  else if ( str[0] == 'n') {
			   gmp_sscanf (str, "n: %Zd\n", rs->n);
			   (flag) ^= (1<<9);
		  }

		  else
			   continue;

		  if (flag == 1023UL) {

			   /* pre-compute and setup rs */
			   rsstr_setup (rs);

			   bestpoly_t bestpoly;
			   bestpoly_init (bestpoly, rs->d);

			   /* set best poly */
			   int i;
			   for (i = 0; i <= rs->d; i++)
			   {
					mpz_set (bestpoly->f[i], rs->f[i]);
			   }
			   mpz_set (bestpoly->g[0], rs->g[0]);
			   mpz_set (bestpoly->g[1], rs->g[1]);

			   fprintf (stderr, "\n# Polynomial (# %5d).\n", count);
			   print_poly_info (rs->f, rs->g, rs->d, rs->n, rs->m, 2);

			   /* start main rootsieve function */
			   rootsieve_main (rs, bestpoly, param, 2);
			   bestpoly_free (bestpoly, rs->d);

			   count += 1;
			   flag = 0UL;
		  }
	 }

	 /* free */
	 rsstr_free (rs);
}


/*
  Do the root sieve for one polynomial from stdin.
  For debugging purpose.
*/
void
rootsieve_stdin ( param_t param )
{
	 /* rootsieve_struct */
	 rsstr_t rs;
	 rsstr_init (rs);

	 /* read poly to rs */
	 read_ggnfs (rs->n, rs->f, rs->g, rs->m);

	 /* pre-compute and setup rs */
	 rsstr_setup (rs);

	 bestpoly_t bestpoly;
	 bestpoly_init (bestpoly, rs->d);

	 /* set best poly */
	 int i;
	 for (i = 0; i <= rs->d; i++)
	 {
		  mpz_set (bestpoly->f[i], rs->f[i]);
	 }
	 mpz_set (bestpoly->g[0], rs->g[0]);
	 mpz_set (bestpoly->g[1], rs->g[1]);

	 fprintf (stderr, "\n# Polynomial (# 0).\n");
	 print_poly_info (rs->f, rs->g, rs->d, rs->n, rs->m, 2);

	 /* start main rootsieve function */
	 rootsieve_main (rs, bestpoly, param, 2);

	 rsstr_free (rs);
	 bestpoly_free (bestpoly, rs->d);
}
