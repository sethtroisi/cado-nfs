/*
  Root sieve for degree 5 polynomials.

  [1. Run and formats]

  For one polynomial in stdin,
  [bash$] ./rootsieve5 < poly.rsa190 > OUTPUT 2> ERRINFO
  where OUTPUT contains results.

  For several polynomials in a file,
  [bash$] ./rootsieve5 -f poly.rsa190 > OUTPUT 2> ERRINFO

  The format of the input polynomial and output polynomials
  are in CADO format.

  For instance, (another RSA190 poly)

$cat OUTPUT | grep "alpha: -8.81" -B 10 -A 1
n: 1907556405060696491061450432646028861081179759533184460647975622318915025587184175754054976155121593293492260464152630093238509246603207417124726121580858185985938946945490481721756401423481
c5: 255190140
c4: -53848884782898166
c3: 12145206120084356715262130
c2: 19191955705982033644872060007757
c1: -6389745371991562680529625140501961075766
c0: 36296980536983889816218353601467302801149433240
Y1: 2642639550249635903
Y0: -1495280603346979841037701044714157977
m: 295069311424266280232575503850212943400547646818258930277687211148503810942166908587111397436966880114029226961961653073130638698706117118357920306697651826715993311914433392591688423380145
# skew: 32808960.00, lognorm: 63.79, alpha: -8.81, E: 54.98, nr: 3
# MurphyE: 1.30e-14 (Bf=10000000, Bg=5000000, area=1.00e+16)

  [2. Parameters]

  Please add parameters in rsparam_init() for different sizes of
  numbers, and this may also depending on polynomial, like skewness.

  [3. Algorithm]

  There are two steps in the root sieve.

  Given a polynomial, the code first tries to find good
  (u, v) (mod (p1^e1 * p2^e2* \cdots *pn^en)) for small prime powers.
  This step is done by hensel-lift-like procedure in a p-ary tree data
  structure (for each such p) and discard any bad/impossible (u, v)
  during the tree-building. The dual roots are updated in a way following
  Emmanuel Thome's idea. Finally, Finally we use CRT to find some good
  (u, v) pairs, those with small u. Note if there are too many pi^ei,
  the CRTs will domimated the running time.

  In the second step, we do the actually root sieve for all primes
  up to some bound. We only consider the (u, v) points lying on the
  sublattice within the sieving range (U, V). For those p appearing
  in the modulus of the sublattice, we have already considred their
  valuations from the previous step and hence we may ignore them in
  the actual sieve. For other primes, we sieve them upt to some bound.
  Currently, the code doesn't consider any p^e for e>1 for these primes.
  This seems to be fine for some polynomials for RSA 190.

  2. The code to deal with special-u (in the root sieve stage) is slow and
  since most of the (u, v) pairs lying on the sublattice have dual roots
  over small primes. Maybe we could do root sieve over small prime powers.

  [4. TODO and history]

  -- reduced memory usage in the crt in find_sublattices()
  -- computed rotate_bounds automatically depending on some information of the
     polynomial. However, this is (too) inaccurate now (TBA)
  -- changed quick_sort to quick_selection
  -- consider sieve in cache blocks (TBA)
  -- put the lengthy file aparts (TBA)
  -- the code is not intelligent. We need to set parameters (TBA). These
     parameter mainly determines the running time between stage 1 and stage 2,
	 and control the trade-offs between speed and accuracy. Until now, it seems
	 that there is no need to consider many primes in stage 2.
  -- beside parameter, rsbound_setup_AB_bound() decides the minimum size for
     sieving, in case the computed rotation bound is too small. Currenly it
	 uses 10^7, and this size can be quickly doen on modern computer.

  [5. Bugs]

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

#define MAX_LINE_LENGTH 4096
#define BLOCK_SIZE 12288
#define TUNE_SIEVE_SIZE BLOCK_SIZE / 2

#define PI 3.14159265358979324
#define MAX_DEGREE 6
#define DEBUG 0
#define SUP_ALPHA 4.843

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
	 double skew, logmu, alpha, e;
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

	 if (verbose) {

		  /* output original poly */
		  gmp_printf ("\nn: %Zd\n", N);
		  for (i = d; i >= 0; i --) {
			   gmp_printf ("c%d: %Zd\n", i, f[i]);
		  }
		  for (i = 1; i >= 0; i --) {
			   gmp_printf ("Y%d: %Zd\n", i, g[i]);
		  }
		  if (verbose == 2) // don't want m in general
			   gmp_printf ("m: %Zd\n", M);
	 }

	 /* compute skew, logmu, nroots */
	 nroots = numberOfRealRoots (f, d, 0, 0);
	 skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
	 logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);
	 alpha = get_alpha (f, d, ALPHA_BOUND);

	 mpz_set (cpoly->n, N);
	 cpoly->degree = d;
	 cpoly->degreeg = 2;
	 cpoly->skew = skew;
	 e = MurphyE (cpoly, BOUND_F, BOUND_G, AREA, MURPHY_K);

	 if (verbose) {
		  printf ("# skew: %.2f, ", skew);
		  printf ("lognorm: %.2f, alpha: %.2f, E: %.2f, nr: %u \n# MurphyE: %1.2e (Bf=%.0f, Bg=%.0f, area=%1.2e)\n",
				  logmu,
				  alpha,
				  logmu + alpha,
				  nroots,
				  e,
				  BOUND_F,
				  BOUND_G,
				  AREA );
	 }

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
ffree_tree ( node *root )
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
  any overlapping. In out case, we call alloc_r_node(ptr, 1)
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
  --- If not, add r and/or "is_dual=k".
  - If (u, v) doesnot exit;
  -- Add a node with (u, v, r, curr_e, val) and/or "is_dual=k".

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
  Insert a listnode to current list, record all.
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
	 /* solve x in a + p1*x = b (mod p) */
	 tmp = solve_lineq (a, p1, b, p2);
	 tmp = tmp * p1 + a;
	 return tmp;
}


#if 0
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
		  /* easy ! */
		  int e;
		  e = poly_roots_ulong(NULL, f, d, p);
		  return (pd * e) / (pd * pd - 1);
	 }
	 else if (pvaluation_disc == 1) {
		  /* special case where p^2 does not divide disc */
		  int e;
		  e = poly_roots_ulong(NULL, f, d, p);
		  /* something special here. */
		  return (pd * e - 1) / (pd * pd - 1);
	 }
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
static double
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
	 alpha =  (- e) * log (2.0);

	 /* FIXME: generate all primes up to B and pass them to get_alpha */
	 for (p = 3; p <= B; p += 2)
		  if (isprime (p)) {
			   e = special_valuation_affine (f, d, p, disc);
			   alpha += (- e) * log ((double) p);
		  }
	 mpz_clear (disc);
	 return alpha;
}


#if 0
/*
  Contribution from a particular dual root r of the polynomial f
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
uv2ab_mod ( long A,
			long MOD,
			long u,
			unsigned long p )
{
	 /* solve_lineq only accept unsigned long */
	 if (A < 0)
		  A = A - p * (long) floor ((double)A / (double)p);

	 /* compute the A + MOD * a = u (mod p) */
	 return solve_lineq(A, MOD, u, p);
}


/*
  Change coordinate from (a, b) to (u, v),
  where A + MOD*a = u.
*/
static inline long
ab2uv ( long A,
		long MOD,
		long a )
{
	 return ( A + MOD * a );
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
static inline long
ij2uv ( long A,
		long MOD,
		long Amin,
		long i )
{
	 return ( ab2uv(A, MOD, ij2ab(Amin, i)) );
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
  Similar to above, but test whether r is a dual root.
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

  (1) If it can be lifted, then it could be single or dual root
  -- single root return 1, -- note, the lifted root is not r in general.
  Hence, we need to compute the lifted root and save it in r_lifted.
  -- dual root return 2.
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
		single or dual;
		- If it is a single, computer r_lifted
		- If it is a dual, see whether it can be lifted
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
		whether it is a	single or dual root;
		- If it is a single, computer r_lifted
		- If it is a dual, see whether it can be lifted
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
		need to check whether it is single or dual. */
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
  Find good sublattice, the lifted cases. For convenience,
  use recursive calls. Note that we could also combine
  this and find_sublattice (base case) together.
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

		  /* save all the dual and single roots for this node. */
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
			 -- r is dual, solve lifted (u, v) who has this r + i*p^k as
			 dual roots; Also for these (u, v), lift any possible single r;
			 -- Note if it is single root. We ignore it. r could be a
			 single root for some other pair (u', v') where (u', v') has
			 some other dual root. However, this situation will be
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
							  printf ("fr: %lu, r: %lu,  (uu, vv): (%lu, %lu) -> (%lu, %lu) (non-invertible, dual) \n",
									  fr, r[nroots], uu, fr, currnode->u + pem1 * uu,
									  currnode->v + pem1 * fr);

						 /* since now r is a dual root, add r + k * p^{e-1}. */
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
							  printf ("fr: %lu, r: %lu, (uu, vv): (%lu, %lu) -> (%lu, %lu) (invertible, dual) \n",
									  fr, r[nroots], uu, vv, currnode->u + pem1 * uu, currnode->v + pem1 * vv);

						 /* since now r is a dual root, add r + k * p^{e-1}. */
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

		  /* If current node is the 2nd bottom leve, add the bottom level leaves
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
  Only consider those (u, v) which has at least one dual root;
  Ignore those pairs which have no dual root. However, for the
  considered (u, v) pairs, we consider all the roots of them.
*/
static void
find_sublattice ( listnode **top,
				  rsstr_t rs,
				  unsigned int p,
				  unsigned int e )
{

	 /* some variables */
	 mpz_t fx_tmp, gx_tmp, numerator_tmp, tmp, *f, *g, *fx, *gx, *numerator;
	 unsigned long pe, i, r, v, a, b, u, *f_ul, *g_ul, *fuv_ul, r_lifted[1];
	 node *currnode, *root;
	 int d, k;

	 new_tree(&root); // freed in an on-line way inside find_sublattice_lift().
	 root = new_node ();

	 mpz_init (fx_tmp);
	 mpz_init (gx_tmp);
	 mpz_init (tmp);
	 mpz_init (numerator_tmp);
	 f = rs->f;
	 g = rs->g;
	 fx = rs->fx;
	 gx = rs->gx;
	 numerator = rs->numerator;
	 d = rs->d;

	 /* compute p^e, shouldn't overflow. */
	 pe = 1UL;
	 for (i = 0; i < e; i ++)
		  pe = pe * p;

	 /* record f, g (mod pe) for each pe. */
	 f_ul = (unsigned long*) malloc ((d + 1) * sizeof (unsigned long));
	 fuv_ul = (unsigned long*) malloc ((d + 1) * sizeof (unsigned long));
	 g_ul = (unsigned long*) malloc ((2) * sizeof (unsigned long));
	 if ((f_ul == NULL) || (g_ul == NULL) || (fuv_ul == NULL)) {
		  fprintf (stderr, "Error, cannot allocate memory in find_sublattice(). \n");
		  exit (1);
	 }

	 /* compute f (mod pe), then we only refer to this instead
		of mpz_t. Note f (mod p) = (f (mod pe)) (mod p). */
	 reduce_poly_ul (f_ul, f, d, pe);
	 reduce_poly_ul (g_ul, g, 1, pe);

	 /* for each root 0 <= r < p  */
	 for (r = 0; r < p; r ++) {

		  /* compute f(r), g(r), numerator(r).
			 if already computed in fx, then copy;
			 other wise compute now.*/
		  if (r < primes[NP-1]) {
			   mpz_set (fx_tmp, fx[r]);
			   mpz_set (gx_tmp, gx[r]);
			   mpz_set (numerator_tmp, numerator[r]);
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
		  /* call solve_lineq function on b + a*v = 0 (mod p) */
		  v = solve_lineq (b, a, 0, p);

		  /* compute fuv (mod p) -- we use isroot_fuv_ul to detect directly.
			 hence don't need to compute f_uv and call isroot_f_ul().
			 compute_fuv_ul (fuv_ul, f_ul, g_ul, d, u, v, p);		  */

		  /* For this (u, v) pair which is already known to have a
			 dual root r, we search for any other possible single
			 and  dual roots (not necessary for dual since other r's
			 will generate the same (u, v). ). We exhaustively search. */
		  for (i = 0; i < p; i ++) {

			   /* r_lifted is not correct, but not important here. Only k is used. */
			   k  = isroot_fuv_ul (f_ul, g_ul, d, u, v, i, p, p, r_lifted);
			   if (DEBUG)
					printf ("(u, v): %lu, %lu  r: %lu, is_root: %d\n", u, v, i, k);

			   /* if i is some root, either single or dual. */
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
		  find_sublattice_lift (root->firstchild, top, f_ul, g_ul, fuv_ul, d, p, e, 2);
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

/*
  Return all sublattices by calling CRT, where for each sublattice,
  the seperate (mod p) valuations are the best.
*/
static int
return_all_sublattice ( rsstr_t rs,
						unsigned short *e,
						unsigned short len_e,
						unsigned long *modulus,
						unsigned long U_bound,
						listnode **all_sublattice,
						int verbose )
{
	 /* at least consider three primes */
	 if (len_e < 2) {
		  fprintf (stderr,
				   "Error: At least consider two primes (2, 3) in return_all_sublattice. \n");
		  exit(1);
	 }

	 unsigned short i, fail = 0;
	 unsigned long j;
	 unsigned long *pe, *size, *ind, ***sublattice_array, tmpu, tmpv, tmpmod = 0UL;
	 listnode *top, *tmp;

	 /* sublattice[i][length][] save (u, v) for prime[i] */
	 sublattice_array = (unsigned long ***) malloc ( len_e * sizeof (unsigned long **) );
	 size = (unsigned long *) malloc ( len_e * sizeof (unsigned long));
	 ind = (unsigned long *) malloc ( len_e * sizeof (unsigned long));
	 pe = (unsigned long *) malloc ( len_e * sizeof (unsigned long));

	 if ( ( (sublattice_array) == NULL) || (size == NULL) || (ind == NULL) || (pe == NULL) ) {
		  fprintf (stderr, "Error, cannot allocate memory in return_all_sublattice(). \n");
		  exit (1);
	 }

	 /* for each prime[i], lift the polynomial
		and return (u, v) with best score. */
	 for (i = 0; i < len_e; i ++) {
		  pe[i] = 1UL;

		  for (j = 0; j < e[i]; j ++)
			   pe[i] = pe[i] * primes[i];

		  /* put best (u, v) into list */
		  new_list (&top);
		  find_sublattice (&top, rs, primes[i], e[i]);
		  size[i] = count_list (top);

		  if (verbose)
			   fprintf ( stderr,
						 "# Info: p: %2u, p^e: %3lu, size_list: %lu\n",
						 primes[i], pe[i], size[i]);

		  /* check whether current list is valid. "len=0" shows
			 the polynomials has no root p, we discard it. */
		  if (size[i] == 0) {
			   fail = 1;
		  }

		  tmp = top;
		  /* allocate array for the i-th prime. */
		  (sublattice_array)[i] = (unsigned long **) malloc ( size[i] * sizeof(unsigned long *));
		  if ( (sublattice_array)[i] == NULL ) {
			   fprintf (stderr, "Error, cannot allocate memory in return_all_sublattice(). \n");
			   exit (1);
		  }
		  for (j = 0; j < size[i]; j ++) {
			   (sublattice_array)[i][j] = (unsigned long *) malloc ( 2 * sizeof(unsigned long));
			   if ( sublattice_array[i][j] == NULL ) {
					fprintf (stderr, "Error, cannot allocate memory in return_all_sublattice(). \n");
					exit (1);
			   }
			   sublattice_array[i][j][0] = tmp->u;
			   sublattice_array[i][j][1] = tmp->v;
			   tmp = tmp->next;
			   //printf ("#%lu, (%lu, %lu)\n", j, sublattice_array[i][j][0], sublattice_array[i][j][1]);
		  }
		  free_list (&top);
	 }

	 /* some len = 0, this set of parameters fails */
	 if (fail == 1) {

		  /* clearence */
		  for (i = 0; i < len_e; i ++) {
			   for (j = 0; j < size[i]; j ++) {
					free (sublattice_array[i][j]);
			   }
			   free (sublattice_array[i]);
		  }
		  free (sublattice_array);
		  free (size);
		  free (ind);
		  free (pe);
		  return 0;
	 }

	 /* Loop over combinations of all arrays. This is awkward.
		We could map 0 ... \prod pe[i] to the indices of the
 		arrays in the price of using more arithmetic. */

	 /* 2, 3, 5, 7 */
	 if (len_e == 4) {
		  for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
			   for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
					for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
						 for (ind[3] = 0; ind[3] < size[3]; ind[3] ++) {

							  tmpmod = pe[0];
							  modulus[0] = size[0];
							  tmpu = sublattice_array[0][ind[0]][0];

							  /* compute u */
							  for (i = 0; i < len_e - 1; i ++) {
								   tmpu = crt_pair(tmpu, tmpmod,
												   sublattice_array[i+1][ind[i+1]][0], pe[i+1]);

								   tmpmod = tmpmod * pe[i+1];
								   modulus[0] = modulus[0] * size[i+1];
							  }

							  /* if u is good, compute v */
							  if ( (tmpu <= U_bound) || (tmpmod - U_bound <= tmpu) ) {

								   tmpmod = pe[0];
								   modulus[0] = size[0];
								   tmpv = sublattice_array[0][ind[0]][1];

								   for (i = 0; i < len_e - 1; i ++) {
										tmpv = crt_pair(tmpv, tmpmod,
														sublattice_array[i+1][ind[i+1]][1], pe[i+1]);

										tmpmod = tmpmod * pe[i+1];
										modulus[0] = modulus[0] * size[i+1];
								   }
								   top = new_listnode (tmpu, tmpv, 0.0);
								   top->next = (*all_sublattice);
								   (*all_sublattice) = top;
							  }
						 }
	 }

	 /* 2, 3, 5, 7, 11 */
	 else if (len_e == 5) {
		  for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
			   for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
					for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
						 for (ind[3] = 0; ind[3] < size[3]; ind[3] ++)
							  for (ind[4] = 0; ind[4] < size[4]; ind[4] ++)	{

								   tmpmod = pe[0];
								   modulus[0] = size[0];
								   tmpu = sublattice_array[0][ind[0]][0];

								   /* compute u */
								   for (i = 0; i < len_e - 1; i ++) {
										tmpu = crt_pair(tmpu, tmpmod,
														sublattice_array[i+1][ind[i+1]][0], pe[i+1]);

										tmpmod = tmpmod * pe[i+1];
										modulus[0] = modulus[0] * size[i+1];
								   }

								   /* if u is good, compute v */
								   if ( (tmpu <= U_bound) || (tmpmod - U_bound <= tmpu) ) {

										tmpmod = pe[0];
										modulus[0] = size[0];
										tmpv = sublattice_array[0][ind[0]][1];

										for (i = 0; i < len_e - 1; i ++) {
											 tmpv = crt_pair(tmpv, tmpmod,
															 sublattice_array[i+1][ind[i+1]][1], pe[i+1]);

											 tmpmod = tmpmod * pe[i+1];
											 modulus[0] = modulus[0] * size[i+1];
										}
										top = new_listnode (tmpu, tmpv, 0.0);
										top->next = (*all_sublattice);
										(*all_sublattice) = top;
								   }
							  }
	 }

	 /* 2, 3, 5, 7, 11, 13 */
	 else if (len_e == 6) {
		  for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
			   for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
					for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
						 for (ind[3] = 0; ind[3] < size[3]; ind[3] ++)
							  for (ind[4] = 0; ind[4] < size[4]; ind[4] ++)
								   for (ind[5] = 0; ind[5] < size[5]; ind[5] ++) {

										tmpmod = pe[0];
										modulus[0] = size[0];
										tmpu = sublattice_array[0][ind[0]][0];

										/* compute u */
										for (i = 0; i < len_e - 1; i ++) {
											 tmpu = crt_pair(tmpu, tmpmod,
															 sublattice_array[i+1][ind[i+1]][0], pe[i+1]);

											 tmpmod = tmpmod * pe[i+1];
											 modulus[0] = modulus[0] * size[i+1];
										}

										/* if u is good, compute v */
										if ( (tmpu <= U_bound) || (tmpmod - U_bound <= tmpu) ) {

											 tmpmod = pe[0];
											 modulus[0] = size[0];
											 tmpv = sublattice_array[0][ind[0]][1];

											 for (i = 0; i < len_e - 1; i ++) {
												  tmpv = crt_pair(tmpv, tmpmod,
																  sublattice_array[i+1][ind[i+1]][1], pe[i+1]);

												  tmpmod = tmpmod * pe[i+1];
												  modulus[0] = modulus[0] * size[i+1];
											 }
											 top = new_listnode (tmpu, tmpv, 0.0);
											 top->next = (*all_sublattice);
											 (*all_sublattice) = top;
										}
								   }
	 }

	 /* 2, 3, 5, 7, 11, 13, 17 */
	 else if (len_e == 7) {
		  for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
			   for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
					for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
						 for (ind[3] = 0; ind[3] < size[3]; ind[3] ++)
							  for (ind[4] = 0; ind[4] < size[4]; ind[4] ++)
								   for (ind[5] = 0; ind[5] < size[5]; ind[5] ++)
										for (ind[6] = 0; ind[6] < size[6]; ind[6] ++) {

											 tmpmod = pe[0];
											 modulus[0] = size[0];
											 tmpu = sublattice_array[0][ind[0]][0];

											 /* compute u */
											 for (i = 0; i < len_e - 1; i ++) {
												  tmpu = crt_pair(tmpu, tmpmod,
																  sublattice_array[i+1][ind[i+1]][0], pe[i+1]);

												  tmpmod = tmpmod * pe[i+1];
												  modulus[0] = modulus[0] * size[i+1];
											 }

											 /* if u is good, compute v */
											 if ( (tmpu <= U_bound) || (tmpmod - U_bound <= tmpu) ) {

												  tmpmod = pe[0];
												  modulus[0] = size[0];
												  tmpv = sublattice_array[0][ind[0]][1];

												  for (i = 0; i < len_e - 1; i ++) {
													   tmpv = crt_pair(tmpv, tmpmod,
																	   sublattice_array[i+1][ind[i+1]][1], pe[i+1]);

													   tmpmod = tmpmod * pe[i+1];
													   modulus[0] = modulus[0] * size[i+1];
												  }
												  top = new_listnode (tmpu, tmpv, 0.0);
												  top->next = (*all_sublattice);
												  (*all_sublattice) = top;
											 }
										}
	 }

	 /* 2, 3, 5, 7, 11, 13, 17, 19 */
	 else if (len_e == 8) {
		  for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
			   for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
					for (ind[2] = 0; ind[2] < size[2]; ind[2] ++)
						 for (ind[3] = 0; ind[3] < size[3]; ind[3] ++)
							  for (ind[4] = 0; ind[4] < size[4]; ind[4] ++)
								   for (ind[5] = 0; ind[5] < size[5]; ind[5] ++)
										for (ind[6] = 0; ind[6] < size[6]; ind[6] ++)
											 for (ind[7] = 0; ind[7] < size[7]; ind[7] ++) {

												  tmpmod = pe[0];
												  modulus[0] = size[0];
												  tmpu = sublattice_array[0][ind[0]][0];

												  /* compute u */
												  for (i = 0; i < len_e - 1; i ++) {
													   tmpu = crt_pair(tmpu, tmpmod,
																	   sublattice_array[i+1][ind[i+1]][0], pe[i+1]);

													   tmpmod = tmpmod * pe[i+1];
													   modulus[0] = modulus[0] * size[i+1];
												  }

												  /* if u is good, compute v */
												  if ( (tmpu <= U_bound) || (tmpmod - U_bound <= tmpu) ) {

													   tmpmod = pe[0];
													   modulus[0] = size[0];
													   tmpv = sublattice_array[0][ind[0]][1];

													   for (i = 0; i < len_e - 1; i ++) {
															tmpv = crt_pair(tmpv, tmpmod,
																			sublattice_array[i+1][ind[i+1]][1], pe[i+1]);

															tmpmod = tmpmod * pe[i+1];
															modulus[0] = modulus[0] * size[i+1];
													   }
													   top = new_listnode (tmpu, tmpv, 0.0);
													   top->next = (*all_sublattice);
													   (*all_sublattice) = top;
												  }
											 }
	 }

	 /* 2, 3, 5 */
	 else if (len_e == 3) {
		  for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
			   for (ind[1] = 0; ind[1] < size[1]; ind[1] ++)
					for (ind[2] = 0; ind[2] < size[2]; ind[2] ++) {

						 tmpmod = pe[0];
						 modulus[0] = size[0];
						 tmpu = sublattice_array[0][ind[0]][0];

						 /* compute u */
						 for (i = 0; i < len_e - 1; i ++) {
							  tmpu = crt_pair(tmpu, tmpmod,
											  sublattice_array[i+1][ind[i+1]][0], pe[i+1]);

							  tmpmod = tmpmod * pe[i+1];
							  modulus[0] = modulus[0] * size[i+1];
						 }

						 /* if u is good, compute v */
						 if ( (tmpu <= U_bound) || (tmpmod - U_bound <= tmpu) ) {

							  tmpmod = pe[0];
							  modulus[0] = size[0];
							  tmpv = sublattice_array[0][ind[0]][1];

							  for (i = 0; i < len_e - 1; i ++) {
								   tmpv = crt_pair(tmpv, tmpmod,
												   sublattice_array[i+1][ind[i+1]][1], pe[i+1]);

								   tmpmod = tmpmod * pe[i+1];
								   modulus[0] = modulus[0] * size[i+1];
							  }
							  top = new_listnode (tmpu, tmpv, 0.0);
							  top->next = (*all_sublattice);
							  (*all_sublattice) = top;
						 }
					}
	 }

	 /* 2, 3 */
	 else if (len_e == 2) {
		  for (ind[0] = 0; ind[0] < size[0]; ind[0] ++)
			   for (ind[1] = 0; ind[1] < size[1]; ind[1] ++) {

					tmpmod = pe[0];
					modulus[0] = size[0];
					tmpu = sublattice_array[0][ind[0]][0];

					/* compute u */
					for (i = 0; i < len_e - 1; i ++) {
						 tmpu = crt_pair(tmpu, tmpmod,
										 sublattice_array[i+1][ind[i+1]][0], pe[i+1]);

						 tmpmod = tmpmod * pe[i+1];
						 modulus[0] = modulus[0] * size[i+1];
					}

					/* if u is good, compute v */
					if ( (tmpu <= U_bound) || (tmpmod - U_bound <= tmpu) ) {

						 tmpmod = pe[0];
						 modulus[0] = size[0];
						 tmpv = sublattice_array[0][ind[0]][1];

						 for (i = 0; i < len_e - 1; i ++) {
							  tmpv = crt_pair(tmpv, tmpmod,
											  sublattice_array[i+1][ind[i+1]][1], pe[i+1]);

							  tmpmod = tmpmod * pe[i+1];
							  modulus[0] = modulus[0] * size[i+1];
						 }
						 top = new_listnode (tmpu, tmpv, 0.0);
						 top->next = (*all_sublattice);
						 (*all_sublattice) = top;
					}
			   }
	 }

	 /* too aggressive */
	 else {
		  fprintf (stderr, "Error, only len_e_sl: 2, 3, 4, 5, 6, 7, 8 is supported at the moment. \n");
		  exit (1);
	 }

	 /* info */
	 if (verbose)
		  fprintf (stderr, "# Info: computed %lu CRTs\n",
				   modulus[0]);

	 /* fix this later. */
	 modulus[0] = tmpmod;

	 /* clearence */
	 for (i = 0; i < len_e; i ++) {
		  for (j = 0; j < size[i]; j ++) {
			   free (sublattice_array[i][j]);
		  }
		  free (sublattice_array[i]);
	 }
	 free (sublattice_array);
	 free (size);
	 free (ind);
	 free (pe);

	 return 1;
}


/*
  Call return_all_sublattice() to return good sublattices.
  Then choose the "rsbound->nbest_sl" ones among them if there
  are more quantities than it. The "best" property is ranked
  by the size of u.

  If "rsbound->nbest_sl" == 0 || >= len; then return all the
  len-best sublattices. This should be the common case if
  len is not too large. "rsbound->nbest_sl" will also be modifed.

  If "rsbound->nbest_sl" != 0 && < len; then this function
  will return the k-best ones by ranking the u. Good ones
  often mean those with small u \in (u, v). We omit the size of v.

  RETURN: *sublattice_array which has k slots.
*/
static int
return_best_sublattice ( rsstr_t rs,
						 rsparam_t rsparam,
						 long ***sublattice_array,
						 int verbose )
{
	 unsigned long i = 0UL, len = 0UL, global_u_bound_rs_tmp;
	 listnode *sublattice_list, *tmp;
	 long ** sublattice_array_all;
	 unsigned short len_e_sl_tmp;

	 /* Repair step: important to get correct length;
		Its assumed that all 0 exponents happens in the tail. */
	 len_e_sl_tmp = rsparam->len_e_sl;
	 global_u_bound_rs_tmp = rsparam->global_u_bound_rs;
	 for (i = 0; i < rsparam->len_e_sl; i ++) {
		  if (rsparam->e_sl[i] == 0)
			   len_e_sl_tmp --;
	 }

	 /* find and return good sublattices to sublattice_list */
	 new_list (&sublattice_list);

	 int re = return_all_sublattice ( rs,
									  rsparam->e_sl,
									  len_e_sl_tmp,
									  &(rsparam->modulus),
									  rsparam->global_u_bound_rs,
									  &sublattice_list,
									  verbose );
	 /* failed, return */
	 if (re == 0) {
		  free_list (&sublattice_list);
		  return 0;
	 }

	 /* decide the length of returnd best sublattices. */
	 len = count_list (sublattice_list);

	 /* if no sublattice is found with u < rsparam->global_u_bound_rs,
		then we try to enlarge u bound. However, it might be better
		to enlarge e_sl[] to allow to check more sublattices. */
	 int count = 0;
	 while (len == 0) {
		  free_list (&sublattice_list);
		  new_list (&sublattice_list);
		  rsparam->global_u_bound_rs *= 2;
		  return_all_sublattice ( rs,
								  rsparam->e_sl,
								  len_e_sl_tmp,
								  &(rsparam->modulus),
								  rsparam->global_u_bound_rs,
								  &sublattice_list,
								  verbose );
		  len = count_list (sublattice_list);
		  if (verbose) {
			   fprintf (stderr,
						"# Warn: Not enough sublattice classes. Reset \"rsparam->global_u_bound_rs = %lu\" (#%d)\n",
						rsparam->global_u_bound_rs, count);

			   /* fprintf (stderr, */
			   /* 			"# Warn: It's better to enlarge 'len_e_sl' and 'e_sl' if want to restrict \"rsparam->global_u_bound_rs\".\n"); */
		  }
		  count ++;
	 }

	 //printf ("rsparam->nbest_sl: %lu, len: %lu\n", rsparam->nbest_sl, len);
	 //if ( (rsparam->nbest_sl == 0) || (rsparam->nbest_sl > len) )

	 if (rs->d == 5)
		  rsparam->nbest_sl = len;
	 else
		  rsparam->nbest_sl = 10;

	 /* info */
	 if (verbose) {
		  fprintf (stderr,
				   "# Info: found %lu sublattices with small u, where |u| < %lu (\"rsparam->global_u_bound_rs = %lu\")\n",
				   len, rsparam->global_u_bound_rs, rsparam->global_u_bound_rs);
		  fprintf (stderr,
				   "# Info: choose %lu best ones among the above all (\"rsparam->nbest_sl = %lu\")\n",
				   rsparam->nbest_sl, rsparam->nbest_sl);
	 }
	 /* now rsbound->nbest_sl <= len. */
	 init_2D_ld (&sublattice_array_all, len, 2);
	 init_2D_ld (sublattice_array, rsparam->nbest_sl, 2);

	 /* copy the list to an 2d array (len*2) and
		free the list in the mean time. */
	  for (i = 0; i < len; i ++) {
		   if (sublattice_list->u < rsparam->modulus - rsparam->global_u_bound_rs) {
				/* sublattice_list->u, v are ul */
				sublattice_array_all [i][0] = (long) sublattice_list->u;
				sublattice_array_all [i][1] = (long) sublattice_list->v;
		   }
		   else {
				sublattice_array_all [i][0] = (long) sublattice_list->u - (long) rsparam->modulus;
				sublattice_array_all [i][1] = (long) sublattice_list->v;
		   }

		  tmp = sublattice_list;
		  sublattice_list = sublattice_list->next;
		  free_listnode (&tmp);
	 }

	 /* do a best-k selection algorithm based on the quick
		sort partition, best means those with small u
		(for size considerations). */
	 quick_sort_2d_ld (sublattice_array_all, 0, 0, len - 1, rsparam->nbest_sl);

	 /* copy to another smaller array -- this is awkward -- */
	 for (i = 0; i < rsparam->nbest_sl; i ++) {
		  if (verbose) {
			   fprintf (stderr, "# Info: sublattice #%4lu, (u, v): (%6ld, %10ld)\n",
						i,
						sublattice_array_all [i][0],
						sublattice_array_all [i][1] );
		  }
		  (*sublattice_array)[i][0] = sublattice_array_all [i][0];
		  (*sublattice_array)[i][1] = sublattice_array_all [i][1];
	 }

	 /* recover */
	 rsparam->global_u_bound_rs = global_u_bound_rs_tmp;

	 free_2D_ld (&sublattice_array_all, len);

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
					 int16_t subu )
{
	 /* Bounding j with B0 < tmp + p*j < B1 */
	 while ( j < V ) {
		  ARRAY[j] = ARRAY[j] - subu;
		  j = j + (long) pe;
	 }
	 return j + (long) pe;
}


/*
  rootsieve_run_lift for the dual root. Since we are
  sure that level 1 node only contains one dual root
  r, we don't need to (and can't, otherwise, may count
  repeatedly) consider single root at all.
*/
static inline void
rootsieve_run_dualroot_lift ( node *currnode,
							  int16_t *ARRAY,
							  unsigned long *f_ul,
							  unsigned long *g_ul,
							  unsigned long *fuv_ul,
							  int d,
							  rsbound_t rsbound,
							  unsigned int p,
							  unsigned int e,
							  unsigned int curr_e,
							  int16_t subu )
{
	 /* recursion end */
	 if (currnode == NULL || curr_e > e)
		  return;

	 /* variables */
	 unsigned int nroots;
	 unsigned long pe, pem1, pep1, fr, gr;
	 node *tmpnode = NULL;
	 long j;

	 /* compute p^e */
	 pem1 = 1UL;
	 for (j = 0; j < curr_e - 1; j ++)
		  pem1 = pem1 * p;
	 pe = pem1 * p;
	 pep1 = pe * p;

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

		  /* insert (u, v), u unchanged, v is fr. to insert roots,
			 r is a multiple root, we need to insert all r + p*j */
		  for (j = 0; j < p; j ++) {
			   /* printf ("insert.. (%lu, %lu) level %lu\n", currnode->u + pem1 * uu, */
			   /* 		  currnode->v + fr *pem1, curr_e); */
			   insert_node ( currnode, &tmpnode,
							 currnode->u,
							 currnode->v + fr * pem1,
							 currnode->r[nroots] + j * pem1,
							 curr_e, pe, 0 );
		  }

		  // next root of current (u, v)
	 }

	 /* recursieve to next level, curr_e + 1 */
	 rootsieve_run_dualroot_lift ( currnode->firstchild,
								   ARRAY,
								   f_ul,
								   g_ul,
								   fuv_ul,
								   d,
								   rsbound,
								   p,
								   e,
								   curr_e + 1,
								   subu );

	 /* delete current node and move to next sibling. */
	 j = uv2ab_mod (rsbound->B, rsbound->MOD, currnode->v, pe);
	 j = (long) ceil (((float) rsbound->Bmin - (float) j) / (float) pe)
		  * (long) pe + (long) j;
	 j = ab2ij (rsbound->Bmin, j);


	 subu = (int16_t) ceil ( (double) subu / (double) pep1 );

	 rootsieve_run_line ( ARRAY,
						  rsbound->Bmax - rsbound->Bmin + 1,
						  j, pe, subu );

	 //rootsieve_run_ij (ARRAY, rsbound, i, j, pe, val / (float) pep1);

	 /* printf ("deleting.. (%lu, %lu) level %lu\n", tmpnode->u, */
	 /* 		  tmpnode->v, curr_e); */

	 free_node (&currnode);
	 tmpnode = NULL;

	 return;
}


/*
  rootsieve_run_lift for the dual root r of the
  congruence class (i, j) (mod p)
*/
static inline  void
rootsieve_run_dualroot ( int16_t *ARRAY,
						 rsstr_t rs,
						 rsbound_t rsbound,
						 unsigned int u,
						 long j,
						 unsigned int r,
						 unsigned int p,
						 unsigned int e,
						 int16_t subu )
{
	 /* some variables */
	 int tmp;
	 unsigned int v;
	 unsigned long pe, *f_ul, *g_ul, *fuv_ul;

	 /* compute p^e */
	 pe = 1UL;
	 for (v = 0; v < e ; v ++)
		  pe = pe * p;

	 /* use s.p instead of m.p */
	 f_ul = (unsigned long*) malloc ((rs->d + 1) * sizeof (unsigned long));
	 fuv_ul = (unsigned long*) malloc ((rs->d + 1) * sizeof (unsigned long));
	 g_ul = (unsigned long*) malloc ((2) * sizeof (unsigned long));
	 if ((f_ul == NULL) || (g_ul == NULL) || (fuv_ul == NULL)) {
		  fprintf (stderr, "Error, cannot allocate memory in rootsieve_run_dualroot(). \n");
		  exit (1);
	 }
	 reduce_poly_ul (f_ul, rs->f, rs->d, pe);
	 reduce_poly_ul (g_ul, rs->g, 1, pe);

	 /* j -> v (mod p) */
	 tmp = ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, j) % p;
	 v = (unsigned int) ( (tmp >= 0) ? tmp : tmp + (int) p );

	 /* we've already known that r is a dual root for f_{u, v}. */
	 node *tmpnode, *root;
	 new_tree (&root);
	 root = new_node ();
	 insert_node (root, &tmpnode, u, v, r, 1, p, 2);

	 /* sieve f(r) + (u*r + (B + MOD * (j+k*p)))*g(r) for k*/
	 rootsieve_run_line ( ARRAY,
						  rsbound->Bmax - rsbound->Bmin + 1,
						  j, p, subu );


	 /* lift to higher p^e */
	 rootsieve_run_dualroot_lift ( root->firstchild,
								   ARRAY,
								   f_ul,
								   g_ul,
								   fuv_ul,
								   rs->d,
								   rsbound,
								   p,
								   e,
								   2,
								   subu );

     /* free, either root itself or root with a level 1 node.
		Note, if e>1, freeing will be done in the lift function;
		However, if e=1, we need to manaully free root->firstchild. */

	 free_node ( &(root->firstchild) );
	 free_node (&root);
	 tmpnode = NULL;

	 free (f_ul);
	 free (fuv_ul);
	 free (g_ul);
}

//#define OLD_ROOTSIEVE_V

#ifdef OLD_ROOTSIEVE_V
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
	 unsigned int np, p, max_e, r, u, v, tmp, tmp1, special_u, fx_ul, gx_ul, nblock;
	 long tmp_ld, i, j;
	 float sub;
	 int16_t subu;

	 // nblock = (rsbound->block rsbound->Bmax - rsbound->Bmin + 1) / 4096;

	 /* For each prime p */
	 for (np = 0; np < rsparam->len_p_rs; np ++) {
		  p = primes[np];

		  /* For each prime not appearing in MOD */
		  if ( (rsbound->MOD) % p != 0 ) {
			   sub = (float) p * log ( (float) p) /
					((float) p * (float) p - 1.0);
			   subu = (int16_t) ceil (sub * 1000.0);

			   /* find u (mod p) which correponds to current i */
			   tmp_ld = ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, fixed_i) % p;
			   u = (unsigned int) ((tmp_ld >= 0) ? tmp_ld : tmp_ld + p);

			   /* approximation level, this should depend on p */
			   /* max_e = (unsigned short) (log((float) rsbound->Bmax - rsbound->Bmin + 1) / logp); */
			   max_e = 1;

			   if (DEBUG)
					fprintf (stderr, "# DEBUG: p: %u, logp_k: %u\n", p, max_e);

			   /* For each possible root < p. */
			   for (r = 0; r < p; r++) {

			   		/* If this happens, f(r) != 0 (mod p). Hence
					   we can through this r. Note f(r), g(r) donot
					   simultaneously = 0 (mod p) in gnfs */
			   		if (mpz_divisible_ui_p(rs->gx[r], p) != 0)
			   			 continue;

					//fprintf (stderr, "# p: %u, r: %u\n", p, r);

			   		/* compute special_u */
					tmp = mpz_fdiv_ui (rs->numerator[r], p);
					tmp1 = mpz_fdiv_ui (rs->gx[r], p);
					special_u = (unsigned int) compute_special_u (tmp, tmp1, p);

					fx_ul = mpz_fdiv_ui (rs->fx[r], p);
					gx_ul = mpz_fdiv_ui (rs->gx[r], p);

					/* find v (mod p) in f(r) + u*r*g(r) + v*g(r) = 0 (mod p) */
					v = (unsigned int) compute_v_ul (fx_ul, gx_ul, r, u, p);
					tmp1 = (unsigned int) uv2ab_mod (rsbound->B, rsbound->MOD, v, p);

					/* smallest tmp + p*i such that A0 < tmp + p*i, where A0 < 0,
					   hence i = ceil((A0-tmp)/p). Note, this should be negative. */
					tmp_ld = (long) ceil ( (float) (rsbound->Bmin -tmp1) / (float) p)
						 * (long) p + (long) tmp1;
					j = ab2ij (rsbound->Bmin, tmp_ld);

					/* sieve array indices (i, j) to real pairs (u, v) */
					if (DEBUG) {
						 gmp_printf ("\nf(%lu): %Zd\n", r, rs->fx[r]);
						 gmp_printf ("g(%lu): %Zd\n",
									 r, rs->gx[r], u, v);
						 gmp_printf ("numerator(%lu): %Zd\n", r, rs->numerator[r]);
						 printf (" -- %ld + [%2u] * %lu = %u (mod %u), ",
								 rsbound->A, tmp, rsbound->MOD, u, p);
						 printf (" %lu + [%2u] * %lu = %u (mod %u), ",
								 rsbound->B, tmp1, rsbound->MOD, v, p);
						 printf (" f(%u) + ( %u * %u + %u ) * g(%u) = 0 (mod %u)\n",
								 r, u, r, v, r, p);

						 printf (" -- %ld + %lu * (%ld + %ld) = %ld = %u (mod %u), ",
								 rsbound->A, rsbound->MOD, i, rsbound->Amin,
								 ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, i), u, p);
						 printf ("%lu + %lu * (%ld + %ld) = %ld = %u (mod %u)\n",
								 rsbound->B, rsbound->MOD, j, rsbound->Bmin,
								 ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, j), v, p);
					}

					/* If r is not a dual root for this f + (u*x+v)*g,
					   then r is not a dual root for
					   f + {[A+MOD*(tmp+p*i)]*x + [J_ST+MOD*(tmp1+p*j)]}*g
					   for	all i, j.  */
					if (u != special_u) {
						 rootsieve_run_line ( ARRAY,
											  rsbound->Bmax - rsbound->Bmin + 1,
											  j, p, subu );
					}
					/* r is multiple root */
					else {
						 rootsieve_run_dualroot ( ARRAY,
												  rs,
												  rsbound,
												  u,
												  j,
												  r,
												  p,
												  max_e,
												  subu );
					}
			   }
		  }
	 }
	 /*
	   print_sievearray (ARRAY, rsbound->Bmin, rsbound->Bmax, rsbound->Amin, rsbound->Amax,
	   rsbound->A, rsbound->B, rsbound->MOD);
	 */
}

#else

/*
  Root sieve for f + u*x*g(x) + v*g(x) for fixed u, variable v.
*/
static void
rootsieve_v ( int16_t *ARRAY,
			  rsstr_t rs,
			  rsbound_t rsbound,
			  rsparam_t rsparam,
			  const long fixed_i,
			  int verbose )
{
	 unsigned int np, nb, p, r, u, v, tmp, max_e;
	 int16_t sub[rsparam->len_p_rs];
	 unsigned long fx_ul, gx_ul;
	 long tmp1, tmp2, start_j_idx[rsparam->len_p_rs], totnb, block_size;
	 float subf;
	 mpz_t tmpz;

	 block_size = BLOCK_SIZE;
	 totnb = (rsbound->Bmax - rsbound->Bmin + 1) / block_size;

	 if (verbose) {
		  fprintf ( stderr,
					"# Stat: totnb: %ld, block_size: %ld, total_size: %ld\n",
					totnb, block_size, totnb*block_size );
	 }

	 mpz_init (tmpz);

	 int total = 0;
	 int st = 0;
	 int c = 0, cs = 0, cm = 0;

	 /* initialize sub[] for all primes in root sieve */
	 for (np = 0; np < rsparam->len_p_rs; np ++) {
		  p = primes[np];
		  subf = (float) p * log ( (float) p) /
			   ((float) p * (float) p - 1.0);
		  sub[np] = (int16_t) ceil (subf * 1000.0);
		  start_j_idx[np] = 0;
	 }

	 /* this is sublattice_u + fixed_i*MOD, */
	 tmp1 = ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, fixed_i);

	 /* For each r < bound_prime */
	 for (r = 0; r < primes[rsparam->len_p_rs]; r++) {

		  /* compute u*g(r)^2 - f(r)g'(r) + f'(r)g(r) */
		  mpz_mul (tmpz, rs->gx[r], rs->gx[r]);
		  mpz_mul_si (tmpz, tmpz, tmp1);
		  mpz_sub (tmpz, tmpz, rs->numerator[r]);

		  st = cputime ();

		  /* For each block */
		  for (nb = 0; nb < totnb; nb ++) {

			   /* For each r < p < Bound_p*/
			   for (np = next_prime_idx[r]; np < rsparam->len_p_rs; np ++) {

					p = primes[np];

					/* skip if need */
					if ((rsbound->MOD) % p == 0)
						 continue;
					if (mpz_divisible_ui_p(rs->gx[r], p) != 0)
						 continue;

					/* compute u (mod p) from u */
					tmp2 = tmp1 % p;
					u = (unsigned int) ((tmp2 >= 0) ? tmp2 : tmp2 + p);

					/* The first block */
					if (nb == 0) {

						 max_e = 1; // should depend on p!!!

						 /* use single precision */
						 fx_ul = mpz_fdiv_ui (rs->fx[r], p);
						 gx_ul = mpz_fdiv_ui (rs->gx[r], p);

						 /* compute v in f(r) + u*r*g(r) + v*g(r) = 0 (mod p) */
						 v = (unsigned int) compute_v_ul (fx_ul, gx_ul, r, u, p);

						 /* v -> j in B + MOD*j = v (mod p) */
						 tmp = (unsigned int) uv2ab_mod (rsbound->B, rsbound->MOD, v, p);

						 /* smallest tmp + p*i such that A0 < tmp + p*i, where A0 < 0,
							hence i = ceil((A0-tmp)/p). Note, this should be negative. */
						 tmp2 = (long) ceil ( (float) (rsbound->Bmin -tmp) / (float) p)
							  * (long) p + (long) tmp;
						 start_j_idx[np] = ab2ij (rsbound->Bmin, tmp2);

						 /* DEBUG */
						 if (DEBUG) {
							  printf (" f(%u) + ( %u * %u + %u ) * g(%u) = 0 (mod %u)\n",
									  r, u, r, v, r, p);

							  printf (" -- %ld + %lu * (%ld + %ld) = %ld = %u (mod %u), ",
									  rsbound->A, rsbound->MOD, fixed_i, rsbound->Amin,
									  ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, fixed_i), u, p);

							  printf ("%lu + %lu * (%ld + %ld) = %ld = %u (mod %u)\n",
									  rsbound->B, rsbound->MOD, start_j_idx[np], rsbound->Bmin,
									  ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, start_j_idx[np]), v, p);
						 }
						 /* !DEBUG */

					}

					c++;

					/* r is simple root for current u, p */
					if (mpz_divisible_ui_p(tmpz, p) == 0)
					{
						 cs ++;
						 if (nb != (totnb - 1)) {
							  start_j_idx[np] = rootsieve_run_line ( ARRAY,
																	 block_size * (nb + 1),
																	 start_j_idx[np], p, sub[np] );
						 }
						 /* boundary */
						 else {
							  start_j_idx[np] = rootsieve_run_line ( ARRAY,
																	 rsbound->Bmax - rsbound->Bmin + 1,
																	 start_j_idx[np], p, sub[np] );
						 }
					}
					/* r is multiple root for current u, p */
					else {
						 cm++;
						 /* don't sieve in block as the max_e may > 1. */
						 if (nb == 0) {
							  rootsieve_run_dualroot ( ARRAY,
													   rs,
													   rsbound,
													   u,
													   start_j_idx[np],
													   r,
													   p,
													   max_e,
													   sub[np] );
						 }
					}
			   }
		  }

		  if (verbose) {
			   fprintf (stderr,
						"# Stat: r = %u,  p >= %u, takes %dms\n",
						r, primes[next_prime_idx[r]], cputime() - st);
		  }

		  total += cputime() - st;
	 }

	 mpz_clear (tmpz);

	 if (verbose)
		  fprintf ( stderr,
					"# Stat: sieve took %dms, total_loops: %d, simple_loops: %d, mul_loops: %d, estimate_mul_loops: %.2f\n",
					total, c, cs, cm,  (double) primes[rsparam->len_p_rs] / log( (double) primes[rsparam->len_p_rs]));

}


#endif


/*
  init root sieve array with biased alpha.
*/
static void
rootsieve_array_init ( int16_t **A,
					   unsigned long jbound,
					   float alpha_bias )
{
	 /* Init array A
		sage: sum ([1/(p-1)*log(p) for p in prime_range(200)])
		4.842766583050838  */

	 float tmpf = SUP_ALPHA + alpha_bias;

	 int16_t tmpu = (int16_t) ceil (tmpf * 1000.0);

	 unsigned long j;
	 /* allocate matrix A. */
	 *A = (int16_t *) malloc ( jbound * sizeof (int16_t));

	 if ((*A) != NULL) {
		  for (j = 0; j < jbound; j++) {
			   (*A)[j] = tmpu;
		  }
	 }
	 else {
		  fprintf (stderr, "Error, cannot allocate memory in rootsieve_array_init().\n");
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
  Init root sieve bound and sublattice.
*/
static void
rsbound_init ( rsbound_t rsbound )
{
	 rsbound->Umax = rsbound->Umin =
		  rsbound->Vmax = rsbound->Vmin =
		  rsbound->Amax = rsbound->Amin = 0;
		  rsbound->Bmax = rsbound->Bmin = 0;

	 rsbound->A = 0;
	 rsbound->B = 0;
	 rsbound->MOD = 0;
}


/*
  Given rsparam->sizebound_ratio_rs and polynomial, set the size for
  sieving region. Set the Amax, Amin, Amax, Amin.
*/
void
rsbound_setup_AB_bound ( rsbound_t rsbound,
						 int verbose )
{
	 /* this is all most always 0 since MOD is much larger than u,
		due to our parameters in rsparam_setup(). Hence we only
		want to sieve over matrix of 1xV. */
	 rsbound->Amax = 0;
	 rsbound->Amin = 0;

	 /* tune mode ? */
	 if (verbose)
		  rsbound->Bmax = 13107200;
	 else
		  rsbound->Bmax = TUNE_SIEVE_SIZE;
	 rsbound->Bmin = -rsbound->Bmax;
}


/*
  Given Amax, Amin, Bmax, Bmin, set the Umax, Umin, Vmax, Vmin.
  Note that they should have similar size as rsparam->global_v_bound_rs.
*/
void
rsbound_setup_sublattice ( rsbound_t rsbound,
						   long sl_A,
						   unsigned long sl_B,
						   rsparam_t rsparam )
{
	 rsbound->A = sl_A;
	 rsbound->B = sl_B;
	 rsbound->MOD = rsparam->modulus;
	 rsbound->Umax = ab2uv (rsbound->A, rsbound->MOD, rsbound->Amax);
	 rsbound->Umin = ab2uv (rsbound->A, rsbound->MOD, rsbound->Amin);
	 rsbound->Vmax = ab2uv (rsbound->A, rsbound->MOD, rsbound->Bmax);
	 rsbound->Vmin = ab2uv (rsbound->A, rsbound->MOD, rsbound->Bmin);
}


/*
  Print root sieve bound and sublattice.
*/
static void
rsbound_print ( rsbound_t rsbound )
{
	 fprintf (stderr, "# Info: (u, v) = (%ld + i * %lu, %lu + j * %lu)\n",
			  rsbound->A, rsbound->MOD, rsbound->B, rsbound->MOD);
	 //fprintf (stderr, "[Umin: %ld, Umax: %ld], ", rsbound->Umin, rsbound->Umax);
	 fprintf (stderr, "# Info: (Amin: %4ld, Amax: %4ld) -> (Umin: %6ld, Umax: %6ld)\n",
			  rsbound->Amin, rsbound->Amax, rsbound->Amin*rsbound->MOD
			  + rsbound->A, rsbound->Amax*rsbound->MOD + rsbound->A);
	 //printf ("[Vmin: %ld, Vmax: %ld], ", rsbound->Vmin, rsbound->Vmax);
	 fprintf (stderr, "# Info: (Bmin: %4ld, Bmax: %4ld) -> (Vmin: %6ld, Vmax: %6ld)\n",
			  rsbound-> Bmin, rsbound->Bmax,  rsbound->Bmin*rsbound->MOD
			  + rsbound->B, rsbound->Bmax*rsbound->MOD + rsbound->B);
}


/*
  Init rsstr_t with.
*/
static void
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

	 /* often set this to 0 unless there are too much sublattices
		found from the find_sublattice() and you want to use the
		first nbest_sl ones. */
	 rsparam->nbest_sl = 0UL;

	 /* the higher, the more margin in computing the sieving bound
		u and v, hence the larger the sieving bound, and hence
		larger e_sl[] in rsparam_setup(). */
	 rsparam->sizebound_ratio_rs = 1.01;

	 /* They will be computed in rsparam_setup */
	 rsparam->global_u_bound_rs = 0UL;
	 mpz_init_set_ui (rsparam->global_v_bound_rs, 0UL);

	 /* number of primes beside e_sl[] considered in second stage
		root sieve. Larger takes longer time, but more accurate. */
	 rsparam->len_p_rs = 20;

	 if (rsparam->len_p_rs >= NP)
		  rsparam->len_p_rs = NP - 1;

	 /* only (further) consider those f_{u, v} which has smaller
		alpha than this which will be set in rsparam_setup(). */
	 rsparam->alpha_bound_rs = 0.0;
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


/*
  For each U bound, identify the best V bound (such that E is
  smallest for this fixed U). Then return the best (U, V) pair.

  Experimentally, this is better than considering U and U*skew
  but will be slower. Anyway, this will only be done once for
  each polynomial (sorry if there are thousands).
*/
static double
rotate_bounds_UV ( mpz_t *f,
				   int d,
				   double ratio_margin,
				   mpz_t b,
				   mpz_t m,
				   unsigned long *U,
				   mpz_t V,
				   int method )
{
	 int i;
	 long k0 = 0, tmpU;
	 double skewness, init_lognorm, E, best_E;
	 mpz_t best_V;

	 mpz_init (best_V);
	 mpz_set_ui (best_V, 0UL);

	 skewness = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, method);
	 init_lognorm = L2_lognorm (f, d, skewness, method);
	 best_E = init_lognorm;

	 /* look for positive k: 2, 4, 8, ... */
	 tmpU = 1;
	 for (i = 0; i < 48; i++, tmpU *= 2)
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

	 /* allow some margin */
	 (*U) *= 2;
	 mpz_add (V, best_V, best_V);
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
				rsstr_t rs )
{
	 double best_E;
	 int size;
	 mpz_t b, m;
	 mpz_init (b);
	 mpz_init (m);
	 mpz_set (b, rs->g[1]);
	 mpz_set (m, rs->g[0]);
	 mpz_neg (m, m);

	 /* compute global bound_v, bound_u.
		-- global_u_bound will be used to identify good sublattice and decide e[],
		-- global_v_bound will be used to identify sieving bound */
	 best_E = rotate_bounds_UV ( rs->f,
								 rs->d,
								 rsparam->sizebound_ratio_rs,
								 b,
								 m,
								 &(rsparam->global_u_bound_rs),
								 rsparam->global_v_bound_rs, /* u, v */
								 DEFAULT_L2_METHOD );

	 /* compute expected alpha */
	 size = mpz_sizeinbase (rsparam->global_v_bound_rs, 2);
	 size += (int) (log ( (double) rsparam->global_u_bound_rs ) * 1.442695);

	 /* sublattice u will often < than this u */
	 gmp_fprintf ( stderr,
				   "# Info: global (u, v) bound (%ld, %Zd) gives best_E: %.3f, exp_min_alpha: %.3f\n",
				   rsparam->global_u_bound_rs,
				   rsparam->global_v_bound_rs,
				   best_E,
				   exp_alpha[size-1] );

	 /* for degree five polynomial */
	 if (rs->d == 5 || rs->d == 6) {

		  /* find_sublattice() consider the first four primes
			 2, 3, 5, 7 with exponents in e[]. */
		  rsparam->len_e_sl = 4;

		  rsparam->e_sl = (unsigned short*) malloc ( rsparam->len_e_sl * sizeof (unsigned short) );
		  if (rsparam->e_sl == NULL) {
			   fprintf (stderr, "Error, cannot allocate memory in rsparam_setup().\n");
			   exit (1);
		  }

		  /* pure experimental, they are the lower bounds in rsparam_tune */
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

	 /* this is the expected min alpha */
	 rsparam->alpha_bound_rs = exp_alpha[size-1];

	 if (rsparam->len_p_rs < rsparam->len_e_sl) {
		  fprintf (stderr, "# Warning: number of primes considered in the root sieve is smaller than that in fin_sublattice(). This might not be accurate. \n");
	 }

	 mpz_clear (b);
	 mpz_clear (m);

	 return;
}


/*
  Free root sieve parameters.
*/
static void
rsparam_free ( rsparam_t rsparam )
{
	 free(rsparam->e_sl);
	 mpz_clear (rsparam->global_v_bound_rs);
}


/*
  Rootsieve for f + (u*x +v)*g for each sublattice
*/
static double
rootsieve_uv ( rsstr_t rs,
			   rsbound_t rsbound,
			   rsparam_t rsparam,
			   mpz_t *fuv,
			   mpz_t *guv,
			   float alpha_bias,
			   int verbose )
{
	 int16_t *MAT, tmp;
	 long i, u, v;
	 unsigned long j;
	 double MurphyE = 0.0, ave_MurphyE = 0.0;
	 int found = 0, WANT;

	 WANT = 64;
	 /* if tune, we want fewer polynomials for speed */
	 if (!verbose)
		  WANT /= 4;

	 /* for each i -> each u = A + MOD * i */
	 for (i = 0; i < (rsbound->Amax - rsbound->Amin + 1); i ++)
	 {
		  /* init the sieve array */
		  rootsieve_array_init (&MAT, rsbound->Bmax - rsbound->Bmin + 1, alpha_bias);

		  int st = cputime ();

		  /* root sieve for v */
		  rootsieve_v (MAT, rs, rsbound, rsparam, i, 0);

		  if (verbose)
			   fprintf ( stderr, "# Info: root sieve took %dms\n",
						 cputime () - st );

		  /* find suitable tmp2 */
		  found = 0;
		  tmp = (int16_t) ceil (rsparam->alpha_bound_rs * 1000.0);
		  while (1) {
			   found = 0;
			   for (j = 0; j < (unsigned long) (rsbound->Bmax - rsbound->Bmin + 1); j++) {

					if (MAT[j] <= tmp)
						 found ++;
			   }

			   if (found < WANT) {
					tmp += 100; // 100 / 1000 = 0.1 in actual alpha

					if (verbose)
						 fprintf ( stderr,
								   "# Warn: Reset \"rsparam->alpha_bound_rs = %.3f\"\n",
								   (float) tmp / 1000.0 );
			   }
			   else
					break;
		  }

		  /* output polynomial */
		  found = 0;
		  for (j = 0; j < (unsigned long) (rsbound->Bmax - rsbound->Bmin + 1); j++) {
			   if (MAT[j] < tmp) {

					found ++;
					u = ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, i);
					v = ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, j);

					compute_fuv (fuv, rs->f, rs->g, rs->d, u, v);
					mpz_set (guv[0], rs->g[0]);
					mpz_set (guv[1], rs->g[1]);

					//optimize (fuv, rs->d, guv, 0, 0); not correct for deg 6 polynomial
					optimize_aux (fuv, rs->d, guv, 0, 0, CIRCULAR);

					if (verbose)
						 fprintf ( stderr,
								   "\n# Found #%2d of (u=%ld, v=%ld) on sublattice #%lu",
								   found, u, v, i );

					MurphyE += print_poly_info (fuv, guv, rs->d, rs->n, rs->m, verbose);
			   }
		  }

		  ave_MurphyE += MurphyE / (double) (found);

		  if (verbose) {
			   fprintf (stderr, "\n");
			   fprintf (stderr, "# Stat: average_MurphyE on %d polynomials: %1.2e (on sublattice %ld, %ld)\n",
						found, MurphyE / (double) (found), rsbound->A, rsbound->B);
		  }

		  /* free sieving array. */
		  free (MAT);
	 }

	 ave_MurphyE = ave_MurphyE / (double) (rsbound->Amax - rsbound->Amin + 1);

	 return ave_MurphyE;
}


/*
  Stage 2: root sieve on the sublattice points;
  For each sublattice, call rootsieve_uv().
*/
static double
rootsieve_main_stage2_run ( rsstr_t rs,
							rsbound_t rsbound,
							rsparam_t rsparam,
							float alpha_proj,
							int verbose )
{
	 int i;
	 mpz_t *fuv, *guv;
	 float alpha_lat, tmp;
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
	 compute_fuv (fuv, rs->f, rs->g, rs->d, rsbound->A, rsbound->B);
	 alpha_lat = (float) get_biased_alpha_affine (fuv, rs->d, primes[rsparam->len_e_sl - 1]);
	 tmp = rsparam->alpha_bound_rs;

	 if (verbose)
		  fprintf (stderr,
				   "# Info: affine_alpha (sublattice): %.3f, projective_alpha: %.3f, total: %.3f, exp_min_alpha: %.3f \n",
			  alpha_lat, alpha_proj, alpha_lat + alpha_proj, tmp);

	 /* On this sublattice, sieve (i, j) */
	 ave_MurphyE = rootsieve_uv ( rs,
								  rsbound,
								  rsparam,
								  fuv,
								  guv,
								  alpha_lat + alpha_proj,
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
								long U,
								long V,
								float alpha_proj,
								int verbose )
{
	 double ave_MurphyE = 0.0;
	 /* init rsbound */
	 rsbound_t rsbound;
	 rsbound_init (rsbound);

	 /* set root sieve bounds, depending on U, V */
	 rsbound_setup_AB_bound (rsbound, verbose);

	 if (verbose)
		  fprintf (stderr, "# Info: sieving matrix size: 1 x [%ld, %ld]\n",
				   rsbound->Bmin, rsbound->Bmax);

	 /* root sieve on this sublattice */
	 rsbound_setup_sublattice ( rsbound,
								U,
								V,
								rsparam );

	 /* compute exact sieving bounds UV given size AB depending
		on current A, B, MOD */
	 if (verbose)
		  rsbound_print (rsbound);

	 /* root sieve */
	 ave_MurphyE = rootsieve_main_stage2_run ( rs,
											   rsbound,
											   rsparam,
											   alpha_proj,
											   verbose );
	 return ave_MurphyE;
}


/*
  Stage 1: find good sublattices.
*/
static double
rootsieve_main_stage1 ( rsstr_t rs,
						rsparam_t rsparam,
						int verbose,
						unsigned long num_sublattices )
{
	 unsigned long i;
	 int st;
	 long ** sublattice_array; /* array contains good sulattices */
	 float alpha_proj;
	 double ave_MurphyE = 0.0, ave2_MurphyE = 0.0;

	 /* alpha projective + alpha contributions from sublattices primes.
		This will be pre-added to the sieve array. */
	 alpha_proj = (float) get_biased_alpha_projective (rs->f, rs->d, 2000);

	 if (verbose)
		  fprintf (stderr, "# Info: projective alpha: %.3f\n", alpha_proj);

	 /* return the first nbest good sulattices to array */
	 st = cputime ();

	 int re = return_best_sublattice (rs, rsparam, &sublattice_array, verbose);

	 /* failed, dont need to free sublattice_array since it has not
		been created in return_best_sublattice(). */
	 if (re == 0) {
		  return -1;
	 }

	 if (verbose)
		  fprintf (stderr, "# Info: find best sublattices over (Mod %lu) took %dms\n",
				   rsparam->modulus, cputime () - st);

	 /* For each sublattice, do the root sieve */
	 for (i = 0; i < rsparam->nbest_sl; i ++) {

		  if (i > num_sublattices)
			   break;

		  if (verbose)
			   fprintf (stderr, "\n# Info: Sieve on sublattice (# %2lu)\n", i);

		  ave_MurphyE = rootsieve_main_stage2_prepare ( rs,
														rsparam,
														sublattice_array [i][0],
														sublattice_array [i][1],
														alpha_proj,
														verbose );
		  ave2_MurphyE += ave_MurphyE;
	 }

	 ave2_MurphyE /= (double) i;

	 /* free */
	 free_2D_ld (&sublattice_array, rsparam->nbest_sl);

	 return ave2_MurphyE;
}


/*
  Tune parameters for find_sublattice().
*/
static double
rsparam_tune ( rsstr_t rs,
			   rsparam_t rsparam,
			   int num_trials )
{
	 unsigned short i, j, best_j = 0, tmp_e_sl[rsparam->len_e_sl], k;
	 unsigned int p, pearr[rsparam->len_e_sl];
	 double ave_MurphyE = 0.0, best_MurphyE = 0.0;

	 /* setup initial */
	 for (i = 0; i < rsparam->len_e_sl; i ++) {
		  p = primes[i];
		  pearr[i] = 1;
		  tmp_e_sl[i] = rsparam->e_sl[i];
		  for (j = 0; j < rsparam->e_sl[i]; j ++) {
			   pearr[i] *= p;
		  }
	 }

	 /* First set of parameters */
	 fprintf (stderr, "# Tune:");
	 for (i = 0; i < rsparam->len_e_sl; i ++) {
		  fprintf (stderr, " %u^%u=%u ", primes[i],
				   rsparam->e_sl[i], pearr[i]);
	 }

	 best_MurphyE = rootsieve_main_stage1 (rs, rsparam, 0, 3UL);

	 if (best_MurphyE == -1)
		  fprintf (stderr, " ave_MurphyE: failed\n");
	 else
		  fprintf (stderr, " ave_MurphyE: %.3e\n", best_MurphyE);

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

		  /* some info and tune */
		  fprintf (stderr, "# Tune:");
		  for (i = 0; i < rsparam->len_e_sl; i ++) {
			   fprintf (stderr, " %u^%u=%u ", primes[i],
						rsparam->e_sl[i], pearr[i]);
		  }

		  ave_MurphyE = rootsieve_main_stage1 (rs, rsparam, 0, 3UL);

		  if (best_MurphyE == -1)
			   fprintf (stderr, " ave_MurphyE: failed\n");
		  else
			   fprintf (stderr, " ave_MurphyE: %.3e\n", ave_MurphyE);

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

	 /* finally, save best parameters back */
	 for (k = 0; k < rsparam->len_e_sl; k ++) {
		   rsparam->e_sl[k] = tmp_e_sl[k];
		   //printf ("e: %u\n", rsparam->e_sl[k]);
	 }

	 return best_MurphyE;
}


/*
  Deg 6
*/
static void
rootsieve_main6 ( rsstr_t rs )
{
	 rsparam_t rsparam, rsparam2;
	 int i, old_i = 0, j;
	 double bestE = 0.0;
	 mpz_t b, m;

	 mpz_init_set (b, rs->g[1]);
	 mpz_init_set (m, rs->g[0]);
	 mpz_neg (m, m);

	 int LEFT = -512, RIGHT = 512;

	 /* Tune parameters:
		fixed a set of parameter e_sl[] for ALL
		(rotated) polynomials: f(x) + i*x^2*g(x). */
	 for (i = 0; i <= RIGHT; i++)
	 {

		  /* rotate polynomial by f + w*x^2 for various w */
		  old_i = rotate_aux (rs->f, b, m, old_i, i, 2);
		  fprintf (stderr, "\n# Info: rotated by %d*x^2\n", i);
		  rsstr_setup (rs);

		  //print_poly_info (rs->f, rs->g, rs->d, rs->n, rs->m, 1);

		  /* set initial parameters */
		  rsparam_init (rsparam);
		  rsparam_setup (rsparam, rs);

		  /* tune parameters in rsparam (try 8 parameter-setups) */
		  bestE = rsparam_tune (rs, rsparam, 8);

		  if (bestE != -1) {
			   break;
		  }
		  rsparam_free (rsparam);
	 }

	 /* rotate back */
	 rotate_aux (rs->f, b, m, old_i, 0, 2);
	 old_i = 0;

	 /* best parameters */
	 fprintf (stderr, "\n# Info: best parameters ");
	 for (j = 0; j < rsparam->len_e_sl; j ++) {
		  fprintf (stderr, "%u:%u ", primes[j], rsparam->e_sl[j]);
	 }
	 fprintf (stderr, "\n");

	 /* Start sieve with above fixed pamaters */
	 for (i = LEFT; i <= RIGHT; i++) {

		  /* rotate polynomial by f + w*x^2 for various w */
		  old_i = rotate_aux (rs->f, b, m, old_i, i, 2);
		  /* reset fx, gx, numerator since f is different */
		  rsstr_setup (rs);

		  fprintf (stderr, "\n# Info: rotated by %d*x^2\n", i);
		  //print_poly_info (rs->f, rs->g, rs->d, rs->n, rs->m, 1);

		  rsparam_init (rsparam2);
		  rsparam_setup (rsparam2, rs);

		  /* save parameters. Note they should have
			 the save length "len_e_sl" */
		  for (j = 0; j < rsparam->len_e_sl; j ++) {
			   rsparam2->e_sl[j] = rsparam->e_sl[j];
		  }

		  /* Call stage 1 & 2. Note, rsparam need to be fixed
			 now since it will be used to find sublattices. */
		  rootsieve_main_stage1 (rs, rsparam2, 1, ULONG_MAX);

		  /* "Tune mode"
			 bestE = rootsieve_main_stage1 (rs, rsparam2, 0, ULONG_MAX);
			 printf ("# TUNE: ave_E: %.3g\n", bestE); */

		  rsparam_free (rsparam2);
	 }

	 rsparam_free (rsparam);

	 /* rotate back */
	 rotate_aux (rs->f, b, m, old_i, 0, 2);

	 mpz_clear (b);
	 mpz_clear (m);
}


/*
  Deg 5
*/
static void
rootsieve_main5 ( rsstr_t rs )
{
	 rsparam_t rsparam;

	 /* set initial parameters */
	 rsparam_init (rsparam);
	 rsparam_setup (rsparam, rs);

	 /* tune parameters in rsparam (try 6 parameter-setups) */
	 rsparam_tune (rs, rsparam, 6);

	 /* Call stage 1 & 2. Note, rsparam need to be fixed
		now since it will be used to find sublattices. */
	 rootsieve_main_stage1 (rs, rsparam, 1, ULONG_MAX);

	 rsparam_free (rsparam);
}


/*
  For the current polynomial (rs), start two-stage root sieve.
*/
static void
rootsieve_main ( rsstr_t rs )
{

	 if (rs->d == 5) {
		  rootsieve_main5 (rs);
	 }
	 else if (rs->d == 6) {
		  rootsieve_main6 (rs);
	 }
	 else {
		  fprintf (stderr, "Error: only support deg 5 or 6.\n");
	 }
}


/*
  Do the root sieve on all polynomials in the file.
*/
static void
rootsieve_file ( FILE *file )
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

			   fprintf (stderr, "\n# Polynomial (# %5d).\n", count);
			   print_poly_info (rs->f, rs->g, rs->d, rs->n, rs->m, 1);

			   /* start main rootsieve function */
			   rootsieve_main (rs);

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
static void
rootsieve_stdin ( void )
{

	 /* rootsieve_struct */
	 rsstr_t rs;
	 rsstr_init (rs);

	 /* read poly to rs */
	 read_ggnfs (rs->n, rs->f, rs->g, rs->m);

	 /* pre-compute and setup rs */
	 rsstr_setup (rs);

	 fprintf (stderr, "\n# Polynomial (# 0).\n");
	 print_poly_info (rs->f, rs->g, rs->d, rs->n, rs->m, 1);

	 /* start main rootsieve function */
	 rootsieve_main (rs);

	 rsstr_free (rs);
}


/*
  Usage
*/
static void
usage (char **argv)
{
	 fprintf (stderr, "# Error: Unexpected argument: %s\n", argv[1]);
	 fprintf (stderr, "# Usage: \n# %s < stdin, where stdin only has one polynomial.\n", argv[0]);
	 fprintf (stderr, "# %s -f filename, where file may contain many polynomials.\n", argv[0]);
	 exit(1);
}


/*
  Main
*/
int
main (int argc, char **argv)
{

	 int i;
	 /* print command-line arguments */
	 fprintf (stderr, "# %s.r%s", *argv, CADO_REV);
	 for (i = 1; i < argc; i++)
		  fprintf (stderr, " %s", *(argv+i));
	 fprintf (stderr, "\n");

	 /* read polynomials from "-f file" */
	 if (argc == 3) {

		  FILE *file = NULL;
		  char *filename = NULL;

		  if (argv[1][0] == '-') {
			   if (strcmp(argv[1], "-f") == 0) {
					filename = argv[2];
					argv += 2;
					argc -= 2;
			   }
			   else {
					usage (argv);
			   }
		  }
		  else
			   usage (argv);

		  // read
		  file = fopen(filename, "r");
		  if (file == NULL) {
			   fprintf(stderr, "# Error in reading file\n");
			   exit (1);
		  }

		  // optimize the raw polys in the file
		  rootsieve_file (file);
		  fclose (file);
		  return 0;
	 }

	 /* read polynomials from stdin*/
	 else if (argc == 1)
	 {
		  rootsieve_stdin ();
		  return 0;
	 }
	 else {
		  usage (argv);
	 }
}
