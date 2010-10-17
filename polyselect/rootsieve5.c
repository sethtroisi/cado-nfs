/*
  Root sieve for degree 5 polynomials.

== USAGE ==

./rootsieve5 < polyfile

== Example ==

$ cat poly.rsa190
Y1: 2642639550249635903
Y0: -1495280603346979841037701044714157977
c5: 255190140
c4: -53848884782898166
c3: 12145206120084356715262130
c2: 19191955694671536369803618342917
c1: 10055610397656329903908039923878920631
c0: 1750978141146217618952056794190360073311557
n: 1907556405060696491061450432646028861081179759533184460647975622318915025587184175754054976155121593293492260464152630093238509246603207417124726121580858185985938946945490481721756401423481

$ ./rootsieve5 < poly.rsa190

<... output omitted ...>

# Totally found (  1120) sublattices, choose the (  10) best ones. (change parameter "nbest" if you want more.)
# sublattice (# 0), (u, v): 500, 654981
# sublattice (# 1), (u, v): 500, 919581
# sublattice (# 2), (u, v): 500, 125781
# sublattice (# 3), (u, v): 500, 390381
# sublattice (# 4), (u, v): 4280, 420621
# sublattice (# 5), (u, v): 4280, 685221
# sublattice (# 6), (u, v): 4280, 156021
# sublattice (# 7), (u, v): 4280, 949821
# sublattice (# 8), (u, v): 8060, 715461
# sublattice (# 9), (u, v): 8060, 186261
# Find sublattice took 0ms

# -- sieving on sublattice (#  0) --
# (u, v) = (500 + i * 1058400, 654981 + j * 1058400)
# (Amin:    0, Amax:    0) -> (Umin:    500, Umax:    500)
# (Bmin: -10000, Bmax: 10000) -> (Vmin: -10583345019, Vmax: 10584654981)

# Found (         500,  -6183576219), alpha: -7.08
# Found (         500,  -2716257819), alpha: -7.03
# Found (         500,    729892581), alpha: -7.27
# Found (         500,   3502900581), alpha: -6.94
# Found (         500,   4355970981), alpha: -6.92
# Found (         500,   8648841381), alpha: -6.92

<... output omitted ...>


== CUSTOM PARAMETERS ==

1. At the moment all parameters are HARDWIRED. Please change them in "CUSTOM
PARAMETERS" section in the main().

2. The most important custom parameters are

"nbest": use the first-nbest sublattices in the root sieve. The root sieve
is done on each of the sublattice. Which (u, v) sublattice is better
is decided by the size of u. The smaller, the better.

"alpha_bound": this refer to the affine alpha value bound. Only those (u, v)
pairs which gives smaller affine_alpha values are reported.

"sieve_size_u" and "sieve_size_v": this is the actually size of the sieve array.
The "sieve_size_u = 0" means only consider the sublattice u itself, and this
is sufficient in practice (since we don't want u to be too large.) Hence the
"sieve_size_v" is the dominant parameter. The larger you set, the larger the
sieving region it runs.


== PROCEDURE ==

1. Given a polynomial, the code first tries to find good
    (u, v) (mod (p1^e1 * p2^e2* \cdots *pn^en))
for small prime powers. This step is done by hensel-lift-like procedure in a
p-ary tree data structure (for each such p) and discard any bad/impossible
(u, v) during the tree-building. Then we need to use CRT to find some good
(u, v) pairs.

2. We do the actually root sieve for all primes up to B. We only consider the
points lying on the sublattice within the sieving range (U, V). For those p
appearing in the points on the sublattice, they normally will have dual
roots and we need special treatments to them.

Currently, we spent more time on part 2 . It might be better to spend more
time on part 1 (for better size property) and less time on part 2.

== TODO ==

1. We are loose with the size of u in the current code. We could only save
those (u, v) with small u and dicard all those not, since any way, they are
less possible to have good size. Also, the code is not intelligent. It is
ideal that, for all sublattice, we do some small rootsieve_run() expriments
to compare them, and do a heavier rootsieve_run() on those good ones.

2. The code to deal with special-u (in the root sieve stage) is slow and
since most of the (u, v) pairs lying on the sublattice have dual roots
over small primes. Maybe we could do root sieve over small prime powers.
For other primes (not on the sublattice) this situation becomes less frequent.

3. The choice of sublattice is the key. Currently it uses a greedy-like way.

4. The file is lengthy. However, it is sort of organized into the following
sets of functions. Put them into new files should be straightforward. Need
to consult with people for a matter of better organisation of files.
- Input-output functions.
- N-ary tree related functions -- change all to non-recursive way.
- listnode functions.

*/


#include "rootsieve5.h"
#define MAX_BUF 4096
#define MAX_DEGREE 6
#define PI 3.14159265358979324
#define DEBUG 0


/*-----------------------------*/
/*   @Input-output functions.  */
/*-----------------------------*/


/*
  Read routines
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
	 char s[MAX_BUF]; /* input buffer */

	 while (feof (stdin) == 0) {
		  ret = readline (s);
		  if (ret == EOF)
			   break;
		  if (strlen (s) + 1 >= MAX_BUF) {
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


/*
   Print polynomial and info: lognorm, skew, alpha.
*/
static void
poly_info ( mpz_t *f,
		   mpz_t *g,
		   int d,
		   mpz_t N,
		   mpz_t M )
{

	 /* print info about the polynomial */
	 unsigned int nroots = 0;
	 double skew, logmu, alpha;
	 int i;

	 /* output original poly */
	 gmp_printf ("\nn: %Zd\n", N);
	 for (i = d; i >= 0; i --) {
		  gmp_printf ("c%d: %Zd\n", i, f[i]);
	 }
	 for (i = 1; i >= 0; i --) {
		  gmp_printf ("Y%d: %Zd\n", i, g[i]);
	 }
	 gmp_printf ("m: %Zd\n", M);

	 /* compute skew, logmu, nroots */
	 nroots = numberOfRealRoots (f, d, 0, 0);
	 skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
	 logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);
	 alpha = get_alpha (f, d, ALPHA_BOUND);
	 printf ("# skew: %.2f, ", skew);
	 printf ("lognorm: %.2f, alpha: %.2f, E: %.2f, nr: %u\n",
			 logmu, alpha, logmu + alpha, nroots);
}


/*-----------------------------*/
/*   @Tree-related functions.  */
/*-----------------------------*/

#if DEBUG
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
	 pnode = (node *)malloc(sizeof(node));
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
  Insert a listnode to current list <- top.
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
}


/*
   Scan a tree, return the best (u, v).
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
  To be changed: find alpha_projective for a poly f.
  Note that, the returned alpha is the  -val * log(p)
  part in the alpha. Hence, we can just add this to
  our affine part.
*/
static double
get_alpha_projective ( mpz_t *f,
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
  Find coordinate i such that
  A + MOD*a = u (mod p).
*/
static inline long
uv2ab_mod ( long A,
			long MOD,
			long u,
			unsigned long p )
{
	 /* compute the (K_ST + MOD * tmp = u) */
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
  Compute fuv = f+(u*x+v)*g, where
  f(r) + u*r*g(r) + v*g(r) = 0 (mod p).
  The inputs for f and g are mpz.
*/
static inline void
compute_fuv ( mpz_t *fuv,
			  mpz_t *f,
			  mpz_t *g,
			  long u,
			  long v)
{
	 mpz_t tmp, tmp1;
	 mpz_init (tmp);
	 mpz_init (tmp1);

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
		  mpz_init (fr[k]);
		  mpz_init (gr[k]);
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

	 /* if (u == 45 && v == 33) */
	 /* 	  gmp_printf("r: %lu, f'(r): %Zd\n", r, tmp); */

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
  Reduce mpz_t *f to unsigned long *f_mod;
  Given modulus pe, return f (mod pe).
*/
static inline void
reduce_f_ul ( unsigned long *f_ul,
			  mpz_t *f,
			  int d,
			  unsigned long pe)
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
   For convenience, use recursive calls.
*/
static void
find_sublattice_lift ( node *firstchild,
					   unsigned long * f_ul,
					   unsigned long * g_ul,
					   unsigned long * fuv_ul,
					   int d,
					   unsigned long p,
					   unsigned int e,
					   unsigned int curr_e )
{
	 if (firstchild == NULL)
		  return;

	 /* some tmp variables */
	 unsigned short roottype = 0;
	 unsigned int i, k, nroots;
	 unsigned long pe, pem1, r_lifted[1], r, uu, vv, fr, gr;
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
	 while ((currnode != NULL)) {

		  if (curr_e > e)
			   return;

		  if (DEBUG) {
			   printf("-----\n");
			   printf("e: %u\n", curr_e);
			   printf("u-v pair: (%lu, %lu) has roots: \n",
					  currnode->u, currnode->v);
			   for (i = 0; i < (currnode->nr); i++)
					printf("\tr[%u]: %lu\n", i, currnode->r[i]);
		  }

		  /*
			Now a pair (u, v) is fixed. Within this pair,
			 - For each single root r;
		       -- For each r + p*i for 0 <= i < p;
			      --- Solve (u, v);
			 - For each dual root r;
			   -- Solve (u, v);
		  */

		  /* find level 1 nodes, those (u, v) (mod p). */
		  l1node = currnode;
		  while (l1node->e != 1)
			   l1node = l1node->parent;

		  for (nroots = 0; nroots < (currnode->nr); nroots++) {
			   /* search all roots of the (u, v) (mod p) to find
				  the ancestor of the current root. */
			   for (i = 0; i < (l1node->nr); i++) {
					if ( l1node->r[i] == (currnode->r[nroots] % p) ) {
						 roottype = l1node->roottype[i];
						 break;
					}
			   }

			   /* check if currnode->r[nroots] is a lifted dual root. if it is,
				  we only need to look at all those (u, v) pairs which have it
				  as a root; otherwise, we have to consider all the lifted cases
				  of currnode->r[nroots], and consider all pairs (u, v) for all
				  of these cases. */

			   /* CASE 1: current r is a dual root which is often in our case. */
			   if (roottype == 2) {
					r = currnode->r[nroots];

					/* compute g(r) */
					gr = eval_poly_ui_mod (g_ul, 1, r, pe);
					if (gr % p == 0)
						 continue;

					/* compute f_uv(x) and then evaluate it at r. */
					compute_fuv_ul (fuv_ul, f_ul, g_ul, d, currnode->u, currnode->v, pe);
					fr = eval_poly_ui_mod (fuv_ul, d, r, pe);
					fr = fr / pem1;
					/* solve on fr + gr*x = 0 (mod p), where x = uu*r + vv. */
					fr = solve_lineq (fr, gr, 0, p);
					fr = fr % p;

					/* we want to solve (uu, vv) in  fr = uu*r + vv (mod p).
					   - if r is not invertible, fix vv and loop all uu.
					   - otherwise, fix vv and solve uu. */
					if (r % p == 0) {
						 for (uu = 0; uu < p; uu ++) {
							  if (DEBUG)
								   printf ("fr: %lu, r: %lu,  (uu, vv): (%lu, %lu) -> (%lu, %lu) (non-invertible, dual) \n",
										   fr, r, uu, fr, currnode->u + pem1 * uu,
										   currnode->v + pem1 * fr);

							  for (k = 0; k < p; k ++) {
								   /* since now r is a single root, consider only r. */
								   insert_node (currnode, &tmpnode, currnode->u + pem1 * uu,
												currnode->v + fr * pem1, r + k * pem1, curr_e, pe, 0);
							  }
						 }
					}

					else {
						 for (vv = 0; vv < p; vv ++) {
							  uu = solve_lineq (vv, r, fr, p);
							  if (DEBUG)
								   printf ("fr: %lu, r: %lu, (uu, vv): (%lu, %lu) -> (%lu, %lu) (invertible, dual) \n",
										   fr, r, uu, vv, currnode->u + pem1 * uu,
										   currnode->v + pem1 * vv);
							  for (k = 0; k < p; k ++) {
								   insert_node (currnode, &tmpnode, currnode->u + pem1 * uu,
												currnode->v + pem1 * vv, r + k * pem1, curr_e, pe, 0);
							  }
						 }
					}
			   }
			   /* CASE2: current r is a single root, it can be lifted forever. but
				we need to find the lifted root for the next recursion. */
			   else if (roottype == 1) {
					/* consider all r + pem1*i */
					for (i = 0; i < p; i++) {

						 r = currnode->r[nroots] + i * pem1;

						 /* compute g(r) */
						 gr = eval_poly_ui_mod (g_ul, 1, r, pe);
						 if (gr % p == 0)
							  continue;
						 /* compute f_uv(x) and then evaluate it at r. */
						 compute_fuv_ul (fuv_ul, f_ul, g_ul, d, currnode->u, currnode->v, pe);
						 fr = eval_poly_ui_mod (fuv_ul, d, r, pe);
						 fr = fr / pem1;
						 /* solve on fr + gr*x = 0 (mod p), where x = uu*r + vv. */
						 fr = solve_lineq (fr, gr, 0, p);

						 /* we want to solve (uu, vv) in  fr = uu*r + vv (mod p).
							- if r is not invertible, fix vv and loop all uu.
							- otherwise, fix vv and solve uu. */
						 if (r % p == 0) {
							  for (uu = 0; uu < p; uu ++) {
								   if (DEBUG)
										printf ("fr: %lu, r: %lu, (uu, vv): (%lu, %lu) -> (%lu, %lu) (non-invertible, single)\n",
												fr, r, uu, fr, currnode->u + pem1 * uu,
												currnode->v + pem1 * fr);

								   /* since now r is a single root, consider insert r. */
								   insert_node (currnode, &tmpnode, currnode->u + pem1 * uu,
												currnode->v + pem1 * fr, r, curr_e, pe, 0);
							  }
						 }
						 else {
							  for (vv = 0; vv < p; vv ++) {
								   uu = solve_lineq (vv, r, fr, p);
								   if (DEBUG)
										printf ("fr: %lu, r: %lu,  (uu, vv): (%lu, %lu) -> (%lu, %lu) (invertible, single)\n",
												fr, r, uu, vv, currnode->u + pem1 * uu,
												currnode->v + pem1 * vv);
								   insert_node (currnode, &tmpnode, currnode->u + pem1 * uu,
												currnode->v + pem1 * vv, r, curr_e, pe, 0);
							  }
						 }
					}
			   }
			   else {
					printf ("roottype: %u, r: %lu, (%lu, %lu) \n", roottype, currnode->r[nroots], currnode->u, currnode->v);
					fprintf (stderr, "Error, something strange in find_sublattice_lift(). \n");
					exit (1);
			   }
		  } // consider next root of current (u, v)

		  find_sublattice_lift (currnode->firstchild, f_ul, g_ul, fuv_ul, d, p, e, curr_e + 1);
		  currnode = currnode->nextsibling;
	 }
	 return;
}


/*
   Find sublattices, the base case.
   Only consider those (u, v) which has at least one dual root;
   Ignore those pairs which have no dual root. However, for the
   considred (u, v) pairs, we considre all the roots of them.
*/
static void
find_sublattice ( node *root,
				  rs_t rs,
				  unsigned long p,
				  unsigned int e )
{
	 /* restriction */
	 assert (e <= 20);
	 if (p >= 5)
		  assert (e <= 8);
	 else if (p >= 11)
		  assert (e <= 3);
	 if (p >= 37) {
		  fprintf (stderr, "Error, only for primes < 37 in find_sublattice(). \n");
		  exit(1);
	 }

	 /* some variables */
	 mpz_t fx_tmp, gx_tmp, numerator_tmp, tmp, *f, *g, *fx, *gx, *numerator;
	 unsigned long pe, i, r, v, a, b, u, *f_ul, *g_ul, *fuv_ul, r_lifted[1];
	 node *currnode;
	 int d, k;

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
	 reduce_f_ul (f_ul, f, d, pe);
	 reduce_f_ul (g_ul, g, 1, pe);

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

	 //print_tree(root, 0);

	 if ((root->firstchild )!= NULL)
		  /* lift to higher p^e */
		  find_sublattice_lift (root->firstchild, f_ul, g_ul, fuv_ul, d, p, e, 2);

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
init_2D_ul ( unsigned long ***array,
			 const unsigned long X,
			 const unsigned long Y )
{
	 unsigned long i, j;

	 /* alloc */
	 (*array) = (unsigned long **) malloc (X * sizeof (unsigned long *));
	 if ( (*array) != NULL) {
		  for (i = 0; i < X; i ++) {
			   (*array)[i] = (unsigned long *) malloc ( Y * sizeof(unsigned long) );
			   if ( (*array)[i] == NULL) {
					fprintf (stderr, "Error, cannot allocate memory in init_2D_ul(). \n");
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
free_2D_ul ( unsigned long ***array,
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
  Init a 2D arrary of type long, the
  dimension X by Y.

*/
static void
init_2D ( long ***array,
		  const unsigned long X,
		  const unsigned long Y )
{
	 unsigned long i, j;

	 /* alloc */
	 (*array) = (long **) malloc (X * sizeof (long *));
	 if ( (*array) != NULL) {
		  for (i = 0; i < X; i ++) {
			   (*array)[i] = ( long *) malloc ( Y * sizeof(long) );
			   if ( (*array)[i] == NULL) {
					fprintf (stderr, "Error, cannot allocate memory in init_2D_ul(). \n");
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
			   (*array)[i][j] = 0;
}


/*
  Free a 2D arrary of type long, the
  dimension X by Y.
*/
static void
free_2D ( long ***array,
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

#if 0
/*
  THE CODE IS NOT CORRECT !
  Do a quick selection of the best k on the
  "dim"-dimension of a 2d array. Note that
  the 2d array must has dimension X*2. The
  parameter "dim" used to indicate sort either
  the array[X][0] or array[X][1]. Any "dim" > 1
  will break the code.
*/
static void
quick_selection_2d_ul ( unsigned long **array,
						const unsigned short dim,
						unsigned long l,
						unsigned long h,
						const unsigned long k )
{
	 unsigned long pivot = 0UL, tmp = 0UL, pivot_i = 0UL;
	 unsigned long i = l, j = h;

	 if (j > i) {
		  /* the middle one, as usual */
		  pivot_i = (i + j) / 2;
		  pivot = array [pivot_i][dim];

		  /* partition */
		  do
		  {

			   while ( array [i][dim] <= pivot )
					i++;
			   while ( array [j][dim] > pivot )
					j--;
			   if (i <= j)
			   {
					tmp = array [i][dim];
					array [i][dim] = array [j][dim];
					array [j][dim] = tmp;

					tmp = array [i][1 - dim];
					array [i][1 - dim] = array [j][1 - dim];
					array [j][1 - dim] = tmp;

					if (i == pivot_i)
						 pivot_i ++;
					if (j == pivot_i)
						 pivot_i --;

					i++;
					j--;
			   }
		  }
		  while (i <= j);

		  /* recursion */
		  if (pivot_i > k)
			   quick_selection_2d_ul ( array, dim, l, pivot_i - 1, k );
		  if (pivot_i < k)
			   quick_selection_2d_ul ( array, dim, pivot_i + 1, h, k );
	 }
}
#endif


/*
  partition.
*/
static unsigned long
quick_sort_2d_ul_partition  ( unsigned long **array,
							  const unsigned short dim,
							  const long l,
							  const long h )
{
	 unsigned long pivot = 0, tmp = 0;
	 long i = l - 1, j;
	 pivot = array [h][dim];
	 for (j = l; j < h; j ++) {
		  if ( array [j][dim] < pivot ) {
			   i ++;
			   tmp = array [i][dim];
			   array [i][dim] = array [j][dim];
			   array [j][dim] = tmp;
			   tmp = array [i][1 - dim];
			   array [i][1 - dim] = array [j][1 - dim];
			   array [j][1 - dim] = tmp;
		  }
	 }
	 tmp = array [i+1][dim];
	 array [i+1][dim] = array [h][dim];
	 array [h][dim] = tmp;
	 tmp = array [i+1][1 - dim];
	 array [i+1][1 - dim] = array [h][1 - dim];
	 array [h][1 - dim] = tmp;
	 return i+1;
}


/*
  Do a quicksort -- change this to quick_selection since
  we only want the best k.
*/
static void
quick_sort_2d_ul ( unsigned long **array,
				   const unsigned short dim,
				   const long l,
				   const long h )
{
	 long pivot = 0;
	 if (l < h) {
		  pivot = quick_sort_2d_ul_partition  (array, dim, l, h);
		  quick_sort_2d_ul (array, dim, l, pivot - 1);
		  quick_sort_2d_ul (array, dim, pivot + 1, h);
	 }
}


/*
  Return sublattices by calling crt.

  "" sidenotes: (u, v): (4280,156021) for rsa190 poly
  could be found through this.
  h = f + (4280*x - 24273189579)*g
  sage: Mod(-24273189579, 2^4*3^3*5^2*7^2)
  156021 ""
*/
static void
return_sublattice ( rs_t rs,
					unsigned short *e,
					unsigned short len_e,
					listnode **sublattice,
					unsigned long *modulus )
{

	 /* at least consider two primes */
	 if (len_e < 2) {
		  fprintf (stderr, "Error: At least consider two primes (2, 3) in find_sublattice. \n");
		  exit(1);
	 }

	 /* consider prime 2 only */
	 unsigned short i, j;
	 unsigned long pe1 = 1UL, pe2 = 1UL;
	 node *root;
	 listnode *top1, *top2, *top12;

	 for (j = 0; j < e[0]; j ++)
		  pe1 = pe1 * primes[0];

	 /* sublattice over 2^e[0] */
	 new_tree (&root);
	 root = new_node ();
	 find_sublattice (root, rs, primes[0], e[0]);
	 new_list (&top1);
	 scan_tree (root, e[0], &top1);
	 free_tree (root);

	 /* consider primes 3, 5, 7 ... */
	 for (i = 1; i < len_e; i ++) {

		  if (DEBUG) {
			   printf ("-- list 1 (mod %lu)--\n", pe1);
			   print_list (top1);
		  }

		  /* compute p_{i}^e_{i} */
		  pe2 = 1UL;
		  for (j = 0; j < e[i]; j ++)
			   pe2 = pe2 * primes[i];

		  /* sublattice over p_{i}^e_{i} */
	 	  new_tree (&root);
	 	  root = new_node ();
		  find_sublattice (root, rs, primes[i], e[i]);
		  new_list (&top2);
		  scan_tree (root, e[i], &top2);
		  free_tree (root);

		  if (DEBUG) {
			   printf ("-- list 2 (mod %lu)--\n", pe2);
			   print_list (top2);
		  }

		  /* crt for each elements in top1 and top2, and
			 save results in top12. !!! MEMORY */
		  new_list (&top12);
		  crt_list (top1, pe1, top2, pe2, &top12);
		  free_list (&top2);
		  free_list (&top1);

		  /* next list */
		  top1 = top12;
		  pe1 = pe1 * pe2;
	 }
	 //print_list (top12);
	 (*sublattice) = top12;
	 modulus[0] = pe1;
}


/*
  Call find_sublattice() and then return the k-best ones.
  Good ones often mean those with small u \in (u, v).
  We omit the size of v.

  !Note that, return to *sublattice_array which has k slots.
*/
static void
return_good_sublattice ( rs_t rs,
						 unsigned short *e,
						 unsigned short len_e,
						 unsigned long *modulus,
						 unsigned long **sublattice_array,
						 const unsigned long k )
{
	 unsigned long i = 0UL, len = 0UL;
	 listnode *sublattice, *tmp;
	 unsigned long ** sublattice_array_all;

	 new_list (&sublattice);

	 /* find and return good sublattices. */
	 return_sublattice (rs, e, len_e, &sublattice, modulus);

	 len = count_list (sublattice);

	 printf ("# Totally found (%6lu) sublattices, choose the (%4lu) best ones. (change parameter \"nbest\" if you want more.)\n", len, k);

	 if (len <= k) {
		  fprintf (stderr, "Error, not enough sublattice classes. Please enlarge 'len_e' and 'e' in CUSTOM PARAMETERS.\n");
		  exit (1);
	 }

	 init_2D_ul (&sublattice_array_all, len, 2);

	 /* copy the list to an 2d array (len*2) and
		free the list in the mean time. */
	 for (i = 0; i < (len - 1); i ++) {

		  /* copy (u, v) in the list to the 2d array. */
		  sublattice_array_all [i][0] = sublattice->u;
		  sublattice_array_all [i][1] = sublattice->v;
		  tmp = sublattice;
		  sublattice = sublattice->next;
		  free_listnode (&tmp);
	 }
	 sublattice_array_all [len - 1][0] = sublattice->u;
	 sublattice_array_all [len - 1][1] = sublattice->v;
	 free_listnode (&sublattice);

	 /* do a best-k selection algorithm based on the quick
		sort partition, best means those with small u
		(for size considerations). */
	 quick_sort_2d_ul (sublattice_array_all, 0, 0, len - 1);

	 /* copy to another smaller array -- this is awkward -- */
	 for (i = 0; i < k; i ++) {
		  ///*
		  printf ("# sublattice (# %lu), (u, v): %lu, %lu\n", i, sublattice_array_all [i][0],
				  sublattice_array_all [i][1]);
		  //*/
		  sublattice_array[i][0] = sublattice_array_all [i][0];
		  sublattice_array[i][1] = sublattice_array_all [i][1];
	 }

	 free_2D_ul (&sublattice_array_all, len);
}



/*-----------------------------*/
/*  @Root sieve                */
/*-----------------------------*/


/*
  Main function for root sieve.
*/
static void
rootsieve_run ( double **ARRAY,
				rs_t rs,
				rsbound_t rsbound )
{
	 unsigned long r, p, u, v, tmp, tmp1, special_u, fx_ul, gx_ul, *r_uv;
	 long tmp2, i, j;
	 int np, d, k, l;
	 mpz_t *f, *g, *fx, *gx, *numerator, *fuv;
	 double sub[NP], val, logp;
	 f = rs->f;
	 g = rs->g;
	 fx = rs->fx;
	 gx = rs->gx;
	 numerator = rs->numerator;
	 d = rs->d;

	 /* fuv is f+(u*x+v)*g */
	 fuv = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
	 if (fuv == NULL) {
		  fprintf (stderr, "Error, cannot allocate memory in rootsieve_run().\n");
		  exit (1);
	 }
	 for (k = 0; k <= d; k++)
		  mpz_init_set (fuv[k], f[k]);
	 /* roots of fuv used for those p | MOD. */
	 r_uv = (unsigned long*) malloc (d * sizeof (unsigned long));
	 if (r_uv == NULL) {
		  fprintf (stderr, "Error, cannot allocate memory in rootsieve_run(). \n");
		  exit (1);
	 }

	 /* For each small prime p. */
	 for (np = 0; np < NP; np ++) {
		  p = primes[np];
		  logp = log ((double) p);
		  sub[np] = (double) p * logp /
			   ((double) p * (double) p - 1.0);

		  /* we consider primes not appearing in MOD */
		  if ((rsbound->MOD) % p != 0) {

			   for (r = 0; r < p; r++) {

			   		/* f(r), g(r) cannot simultaneously = 0 (mod p),
					   hence no such pair (U, V) exists. */
			   		if (mpz_divisible_ui_p(gx[r], p) != 0)
			   			 continue;

			   		/* compute special_u */
					tmp = mpz_fdiv_ui (numerator[r], p);
					tmp1 = mpz_fdiv_ui (gx[r], p);
					special_u = compute_special_u (tmp, tmp1, p);
					fx_ul = mpz_fdiv_ui (fx[r], p);
					gx_ul = mpz_fdiv_ui (gx[r], p);

			   		for (u = 0; u < p; u ++) {

						 /* find v in f(r) + u*r*g(r) + v*g(r) = 0 (mod p) */
						 v = compute_v_ul (fx_ul, gx_ul, r, u, p);

						 /* find i, j (mod p) which correponds to
							the current u, v */
						 tmp = uv2ab_mod (rsbound->A, rsbound->MOD, u, p);
						 tmp1 = uv2ab_mod (rsbound->B, rsbound->MOD, v, p);

						 /* print f_(u, v) and its root r. */
						 if (DEBUG) {
							  gmp_printf ("\nf(%lu): %Zd\n", r, fx[r]);
							  gmp_printf ("g(%lu): %Zd\n",
										  r, gx[r], u, v);
							  gmp_printf ("numerator(%lu): %Zd\n", r, numerator[r]);
							  printf (" -- %lu + [%2lu] * %lu = %lu (mod %lu), ",
									  rsbound->A, tmp, rsbound->MOD, u, p);
							  printf (" %lu + [%2lu] * %lu = %lu (mod %lu), ",
									  rsbound->B, tmp1, rsbound->MOD, v, p);
							  printf (" f(%lu) + ( %lu * %lu + %lu ) * g(%lu) = 0 (mod %lu)\n",
									  r, u, r, v, r, p);
						 }

						 /* smallest tmp + p*i such that A0 < tmp + p*i,
							hence i = ceil((A0-tmp)/p). Note, this should be negative. */
						 tmp2 = ((rsbound->Amin - (long) tmp) / (long) p) * (long) p
							  + (long) tmp;
						 i = ab2ij (rsbound->Amin, tmp2);

						 tmp2 = ((rsbound->Bmin - (long) tmp1) / (long) p) * (long) p
							  + (long) tmp1;
						 j = ab2ij (rsbound->Bmin, tmp2);

						 /* sieve array indices (i, j) to real pairs (u, v) */
						 if (DEBUG) {
							  printf (" -- %lu + %lu * (%ld + %ld) = %ld = %lu (mod %lu), ",
									  rsbound->A, rsbound->MOD, i, rsbound->Amin,
									  ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, i), u, p);

							  printf ("%lu + %lu * (%ld + %ld) = %ld = %lu (mod %lu)\n",
									  rsbound->B, rsbound->MOD, j, rsbound->Bmin,
									  ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, j), v, p);
						 }

						 /* If r is not a dual root for this f + (u*x+v)*g,
							then r is not a dual root for
							f + {[A+MOD*(tmp+p*i)]*x + [J_ST+MOD*(tmp1+p*j)]}*g
							for	all i, j.  */
						 if (u != special_u) {
							  /* Bounding i with B0 < tmp + p*i < B1 */
							  while ( i < (rsbound->Amax - rsbound->Amin + 1) ) {
								   if (DEBUG)
										printf ("[i: %3lu] ", i);
								   while ( j < (rsbound->Bmax - rsbound->Bmin + 1) ) {

										if (DEBUG)
											 printf ("%3lu ", j);
										ARRAY[i][j] = ARRAY[i][j] - sub[np];

										if (DEBUG) {
											 if (i == 27 && j == 220)
												  printf ("case 1: (i=%lu, j=%lu), alpha: %f, r: %lu, p: %lu\n",
														  i, j, sub[np], r, p);
										}
										j = j + (long) p;

								   }
								   if (DEBUG)
										printf ("\n");
								   i = i + (long) p;
								   j = ab2ij (rsbound->Bmin, tmp2);
							  }
						 }

						 /* u is special and r is a dual root of f_uv */
						 else {

							  //printf ("*** p: %lu, r: %lu, potential_special_u: %lu, val: %f\n", p, r, u, val);
							  while ( i < (rsbound->Amax - rsbound->Amin + 1) ) {

								   if (DEBUG)
										printf ("[i: %3lu] ", i);
								   while ( j < (rsbound->Bmax - rsbound->Bmin + 1) ) {

										if (DEBUG)
											 printf ("%3lu ", j);

										/* compute new polynomial f_uv = f + (u*x+v)*g, (u, v) is computed
										 from the current (i, j) using ij2uv(). */
										compute_fuv (fuv, f, g,
													 ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, i),
													 ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, j));
										/* update the contribution of this root r to f_uv for all uv.
										   Note that, we need to look recursively for each u, v pair. */

										val = average_valuation_affine_root (fuv, d, p, r);
										ARRAY[i][j] = ARRAY[i][j] - val * logp;

										j = j + (long) p;
								   }
								   if (DEBUG)
										printf ("\n");
								   i = i + (long) p;
								   j = ab2ij (rsbound->Bmin, tmp2);
							  }
						 }
			   		}
			   }
		  }
		  /* For those p which divides MOD. */
		  else {
			   //printf ("primes[%d]: %u divides %lu\n",np, primes[np], rsbound->MOD);
			   /* only these (u, v) are possible such that A + MOD*i = u (mod p) */
			   u = rsbound->A % p;
			   v = rsbound->B % p;
			   compute_fuv (fuv, f, g, u, v);
			   k = poly_roots_ulong (r_uv, fuv, d, p);

			   /* for all the roots fo this f_uv. */
			   for (l = 0; l < k; l++) {

					r = r_uv[l];
					/* test whether r is dual, we could call roottype(), but its heavier. */
					tmp = mpz_fdiv_ui (numerator[r], p);
					tmp1 = mpz_fdiv_ui (gx[r], p);

					/* test u * tmp1^2 == tmp, whether r is a dual root for this u. */
					if (test_special_u (tmp, tmp1, u, p) != 1) {

						 for (i = 0; i < (rsbound->Amax - rsbound->Amin + 1); i ++) {

								   if (DEBUG)
										printf ("[i: %3lu] ", i);
								   for (j = 0; j < (rsbound->Bmax - rsbound->Bmin + 1); j ++) {

										if (DEBUG)
											 printf ("%3lu ", j);

										if (DEBUG) {
											 if (i == 27 && j == 220)
												  printf ("case 3:(i=%lu, j=%lu), alpha: %f, r: %lu, p: %lu\n",
														  i, j, sub[np], r, p);
										}
										ARRAY[i][j] = ARRAY[i][j] - sub[np];
								   }
								   if (DEBUG)
										printf ("\n");
								   j = 0;
						 }
					}
					/* r is dual root for f_uv */
					else {

						 //printf ("*** p: %lu, r: %lu, potential_special_u: %lu, val: %f\n", p, r, u, val);

						 for (i = 0; i < (rsbound->Amax - rsbound->Amin + 1); i ++) {

							  if (DEBUG)
								   printf ("[i: %3lu] ", i);
							  for (j = 0; j < (rsbound->Bmax - rsbound->Bmin + 1); j ++) {

								   if (DEBUG)
										printf ("%3lu ", j);

								   /* compute new polynomial f_uv = f + (u*x+v)*g, (u, v) is computed
									  from the current (i, j) using ij2uv(). */
								   compute_fuv (fuv, f, g,
												ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, i),
												ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, j));
								   /* update the contribution of this root r to f_uv for all uv.
									  Note that, we need to look recursively for each u, v pair. */

								   val = average_valuation_affine_root (fuv, d, p, r);

								   if (DEBUG) {
										if (i == 27 && j == 220) {
											 printf ("case 4: (i=%ld, j=%ld), alpha: %f, r: %lu, p: %lu\n", i, j, val * logp, r, p);
											 printf ("(u, v): (%ld, %ld), %lu, %lu, %ld\n", ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, i),
													 ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, j), rsbound->A, rsbound->MOD, rsbound->Amin);
										}
								   }

								   ARRAY[i][j] = ARRAY[i][j] - val * logp;

							  }
							  if (DEBUG)
								   printf ("\n");
							  j = 0;
						 }
					}
			   }
		  }
	 }

	 /*
	 print_sievearray (ARRAY, rsbound->Bmin, rsbound->Bmax, rsbound->Amin, rsbound->Amax,
	 rsbound->A, rsbound->B, rsbound->MOD);
	 */
	 for (k = 0; k <= d; k++) {
		  mpz_clear (fuv[k]);
	 }
	 free (fuv);
	 free (r_uv);
}



/*
  init rs with polynomials information.
*/
static void
rootsieve_array_init ( double ***A,
					   unsigned long ibound,
					   unsigned long jbound )
{
	 unsigned long i, j;

	 /* allocate matrix A. */
	 (*A) = (double **) malloc ( ibound * sizeof (double *) );
	 if ((*A) != NULL) {
		  for (i = 0; i < ibound; i ++) {
			   (*A)[i] = (double *) malloc ( jbound * sizeof(double) );
			   if ((*A)[i] == NULL) {
					fprintf (stderr, "Error, cannot allocate memory in main\n");
					exit (1);
			   }
		  }
	 }
	 else {
		  fprintf (stderr, "Error, cannot allocate memory in main\n");
		  exit (1);
	 }

     /* Init array A
		sage: sum ([1/(p-1)*log(p) for p in prime_range(200)])
		4.842766583050838  */
	 for (i = 0; i < ibound; i++) {
		  for (j = 0; j < jbound; j++) {
			   (*A)[i][j] = 4.842767;
		  }
	 }
}

#if 0
/*
  free array
*/
static void
rootsieve_array_free ( double ***A,
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
rootsieve_bound_init ( rsbound_t rsbound )
{
	 rsbound->Umax = rsbound->Umin =
		  rsbound->Vmax = rsbound->Vmin =
		  rsbound->Amax = rsbound->Amin =
		  rsbound->Bmax = rsbound->Bmin = 0;
	 rsbound->A = 0;
	 rsbound->B = 0;
	 rsbound->MOD = 0;
}

/*
  Print root sieve bound and sublattice.
*/
static void
rootsieve_bound_print ( rsbound_t rsbound )
{
	 printf ("# (u, v) = (%lu + i * %lu, %lu + j * %lu)\n", rsbound->A, rsbound->MOD, rsbound->B, rsbound->MOD);
	 //printf ("[Umin: %ld, Umax: %ld], ", rsbound->Umin, rsbound->Umax);
	 printf ("# (Amin: %4ld, Amax: %4ld) -> (Umin: %6ld, Umax: %6ld)\n",
			 rsbound->Amin, rsbound->Amax, rsbound->Amin*rsbound->MOD
			 + rsbound->A, rsbound->Amax*rsbound->MOD + rsbound->A);
	 //printf ("[Vmin: %ld, Vmax: %ld], ", rsbound->Vmin, rsbound->Vmax);
	 printf ("# (Bmin: %4ld, Bmax: %4ld) -> (Vmin: %6ld, Vmax: %6ld)\n", rsbound->Bmin,
			 rsbound->Bmax,  rsbound->Bmin*rsbound->MOD
			 + rsbound->B, rsbound->Bmax*rsbound->MOD + rsbound->B);
}


/*
  Init rs_t with polynomials information.
*/
static void
rootsieve_init ( rs_t rs )
{
	 int i;
	 /* init */
	 mpz_init (rs->n);
	 mpz_init (rs->m);
	 rs->f = (mpz_t*) malloc ((MAX_DEGREE + 1) * sizeof (mpz_t));
	 rs->g = (mpz_t*) malloc ((MAX_DEGREE + 1) * sizeof (mpz_t));
	 if ((rs->f == NULL) || (rs->g == NULL)) {
		  fprintf (stderr, "Error, cannot allocate memory in rootsieve_init().\n");
		  exit (1);
	 }
	 for (i = 0; i <= MAX_DEGREE; i++)
	 {
		  mpz_init (rs->f[i]);
		  mpz_init (rs->g[i]);
	 }

	 /* read poly */
	 read_ggnfs (rs->n, rs->f, rs->g, rs->m);
	 for ((rs->d) = MAX_DEGREE; (rs->d) > 0 && mpz_cmp_ui ((rs->f[rs->d]), 0) == 0; rs->d --);
	 poly_info (rs->f, rs->g, rs->d, rs->n, rs->m);

	 /* pre-computing function values f(r), g(r) for 0 <= r < p*/
	 (rs->fx) = (mpz_t *) malloc ( (primes[NP-1]+1) * sizeof (mpz_t) );
	 (rs->gx) = (mpz_t *) malloc ( (primes[NP-1]+1) * sizeof (mpz_t) );
	 (rs->numerator) = (mpz_t *) malloc ( (primes[NP-1]+1) * sizeof (mpz_t) );
	 if (((rs->fx) == NULL) || ((rs->gx) == NULL) || ((rs->numerator) == NULL)) {
		  fprintf (stderr, "Error, cannot allocate memory for polynomials values.\n");
		  exit (1);
	 }
	 /* save f(l) to an array for all l < B
		can we save the differences? not necessary since anyway
		precomputation is only done for once for one polynomial. */
	 eval_polys (rs->f, rs->g, rs->fx, rs->gx, rs->numerator, primes, rs->d);
}


/*
  Free rs_t.
*/
static void
rootsieve_free ( rs_t rs )
{
	 unsigned long i;
	 /* free fl, gl */
	 for (i = 0; i < primes[NP-1] + 1; i ++)
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
  Usage
*/
static void
usage (char *argv1)
{
  fprintf (stderr, "Unexpected argument: %s\n", argv1);
  fprintf (stderr, "Usage: ./rootsieve5 < POLYFILE \n");
  exit(1);
}


/*
  Main
*/
int
main (int argc, char *argv[])
{
	 /* only use stdin */
	 if (argc > 1)
	 {
           usage (argv[1]);
           exit (1);
	 }

	 int st;
	 unsigned long i, j, k, nbest, modulus[1];
	 rs_t rs;
	 rsbound_t rsbound;
	 double **MAT;

	 /********** CUSTOM PARAMETERS ***********/
	 unsigned short *e;
 	 /* consider the first four primes 2, 3, 5, 7 with exponents in e[] */
	 unsigned short len_e = 4;
	 e = (unsigned short *) malloc ( len_e * sizeof (unsigned short) );
	 e[0] = 5;
	 e[1] = 3;
	 e[2] = 2;
	 e[3] = 2;
	 /* use the nbest good sublattices. */
	 nbest = 10UL;
	 /* only want those with affine alpha < -5 */
	 double alpha_bound = -4;
	 /* only use u in the sublattice itself. This should be fine. */
	 unsigned long sieve_size_u = 0;
	 unsigned long sieve_size_v = 100000;
	 /************* END CUSTOM ***************/

	 /* init and read polynomials */
 	 st = cputime ();
	 rootsieve_init (rs);
	 printf ("# Pre-computing f(r), g(r) took %dms\n", cputime () - st);
 	 st = cputime ();

	 /* contain good sulattices */
	 unsigned long ** sublattice_array;
	 init_2D_ul (&sublattice_array, nbest, 2);

	 /* given sublattice, return the first k good sulattices to	a 2d
		array; the modulus is also set to be p1^e1 * p2^e2 ... */
	 return_good_sublattice (rs, e, len_e, modulus, sublattice_array, nbest);
	 free (e);
	 printf ("# Find sublattice took %dms\n", cputime () - st);

	 /* init sieve bound array and keep it for the loop. */
	 rootsieve_bound_init (rsbound);
	 rsbound->Amax = sieve_size_u;
	 rsbound->Amin = -(rsbound->Amax);
	 rsbound->Bmax = sieve_size_v;
	 rsbound->Bmin = -(rsbound->Bmax);

	 /* this takes too much memory, TBC. */
	 double *alpha_array;
	 long **uv_array;
	 unsigned long size =  (rsbound->Amax - rsbound->Amin + 1) * (rsbound->Bmax - rsbound->Bmin + 1);
	 alpha_array = (double *) malloc ( size * sizeof (double) );
	 if (!alpha_array) {
		  fprintf (stderr, "Error, cannot allocate memory for alpha_array[].\n");
		  exit (1);
	 }
	 init_2D (&uv_array, size, 2);

	 /* alpha valuaes, need to add in the final. */
	 double alpha_p = 0.0;
	 alpha_p = get_alpha_projective (rs->f, rs->d, 200);

	 /* For each sublattice, do the sieve. */
	 for (i = 0; i < nbest; i ++) {

		  /* fix memory usage by fixing A*B space, U*V space should
			 be similar among all sublattices. */
		  rsbound->A = sublattice_array [i][0];
		  rsbound->B = sublattice_array [i][1];
		  rsbound->MOD = modulus[0];
		  rsbound->Umax = ab2uv (rsbound->A, rsbound->MOD, rsbound->Amax);
		  rsbound->Umin = ab2uv (rsbound->A, rsbound->MOD, rsbound->Amin);
		  rsbound->Vmax = ab2uv (rsbound->A, rsbound->MOD, rsbound->Bmax);
		  rsbound->Vmin = ab2uv (rsbound->A, rsbound->MOD, rsbound->Bmin);

		  printf ("\n# -- sieving on sublattice (# %2lu) --\n", i);
		  rootsieve_bound_print (rsbound);

		  /*  keep this for future reference. FIX UV first and transform to AB.
			  given bounds for U, V, compute bounds for A, B. Note that
			  actually areas is reduced slightly, error bounded by MOD
			  rsbound->Amax = uv2ab (rsbound->A, rsbound->MOD, rsbound->Umax);
			  rsbound->Amin = uv2ab (rsbound->A, rsbound->MOD, rsbound->Umin);
			  rsbound->Bmax = uv2ab (rsbound->A, rsbound->MOD, rsbound->Vmax);
			  rsbound->Bmin = uv2ab (rsbound->A, rsbound->MOD, rsbound->Vmin); */

		  /* init the sieve array and do the sieve on the sublattice. */
		  rootsieve_array_init (&MAT, rsbound->Amax - rsbound->Amin + 1,
								rsbound->Bmax - rsbound->Bmin + 1);
		  rootsieve_run (MAT, rs, rsbound);

		  /* change the 2d array to 1d array. ~awkful. Need to change
			 to use a 1d sieve array.  */
		  unsigned long len = 0UL;
		  for (k = 0; k < (unsigned long) (rsbound->Amax - rsbound->Amin + 1); k++)
		  {
			   for (j = 0; j < (unsigned long) (rsbound->Bmax - rsbound->Bmin + 1); j++)
			   {
					if (MAT[k][j] < alpha_bound) {
						 alpha_array[len] = MAT[k][j];
						 uv_array[len][0] = ij2uv (rsbound->A, rsbound->MOD, rsbound->Amin, k);
						 uv_array[len][1] = ij2uv (rsbound->B, rsbound->MOD, rsbound->Bmin, j);
						 len ++;
					}
			   }
			   free (MAT[k]);
		  }

		  printf ("\n");
		  for (k = 0; k < len; k ++) {
			   printf ("# Found (%16ld, %16ld), alpha: %.2f\n", uv_array[k][0], uv_array[k][1], alpha_array[k] + alpha_p);

			   /* NEED optimize, output polynomial here. To be added. */
		  }

		  /* free sieving array. */
		  free (MAT);
	 }

	 /* free */
	 free (alpha_array);
	 free_2D_ul (&sublattice_array, nbest);
	 free_2D (&uv_array, size);
	 rootsieve_free (rs);
	 return 0;
}
