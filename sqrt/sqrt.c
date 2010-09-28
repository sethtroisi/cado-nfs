#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <inttypes.h>
#include <math.h> /* for log */
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>

#include "cado.h"
#include "utils.h"

#define DEBUG 0
#define MAPLE 0

/* Although the functions in plain_poly are not readily available in the
 * publicized interface of utils.h, it's ok to use them if we explicitly
 * include the corresponding header.
 */
#include "plain_poly.h"

//#include "poly.h"

// #define VERBOSE



/********** DEP **********/
int
checkVector(int *vec, int ncols)
{
    int ok, i;

    ok = 1;
    for(i = 0; i < ncols; i++)
	if(vec[i] & 1){
	    ok = 0;
	    break;
	}
#if 0
    if(ok)
	printf(" -> y\n");
    else
	printf(" -> n (%d)\n", i);
#endif
    return ok;
}

void
computefree(relation_t *rel)
{
    int i;

    for(i = 0; i < rel->nb_ap; i++){
        rel->ap[i].r = rel->ap[i].p;
        rel->ap[i].p = rel->a;
        rel->ap[i].e = 1;
    }
}

void
treatRelation(FILE *depfile, hashtable_t *H, relation_t rel)
{
    int j, h;

    if(rel.b > 0)
	computeroots(&rel);
    else
	computefree(&rel);
    fprintf(depfile, "%ld %lu\n", rel.a, rel.b);
#if MAPLE >= 1
    fprintf(stderr, "NORM:=NORM*mynorm(f, %ld, %lu):\n", rel.a, rel.b);
    fprintf(stderr, "AX:=AX*(%ld-%lu*X):\n", rel.a, rel.b);
#endif
    for(j = 0; j < rel.nb_ap; j++){
#if MAPLE >= 2
	fprintf(stderr, "A2:=A2*[%ld, %ld]^%d:\n",
		rel.ap[j].p, rel.ap[j].r, rel.ap[j].e);
#endif
	h = getHashAddr(H, rel.ap[j].p, rel.ap[j].r);
	if(H->hashcount[h] == 0){
	    // new empty place
            SET_HASH_P(H,h,rel.ap[j].p);
            SET_HASH_R(H,h,rel.ap[j].r);
	}
	H->hashcount[h] += rel.ap[j].e;
    }
}

// returns the sign of m1*a+m2*b
int
treatSign(relation_t rel, cado_poly pol)
{
    mpz_t tmp1, tmp2;
    int s;

    /* first perform a quick check */
    s = (rel.a > 0) ? mpz_sgn (pol->g[1]) : -mpz_sgn (pol->g[1]);
    if (mpz_sgn (pol->g[0]) == s)
      return s;

    mpz_init(tmp1);
    mpz_mul_si(tmp1, pol->g[1], rel.a);
    mpz_init(tmp2);
    mpz_mul_ui(tmp2, pol->g[0], rel.b);
    mpz_add(tmp1, tmp1, tmp2);
    s = (mpz_sgn(tmp1) >= 0 ? 1 : -1);
    mpz_clear(tmp1);
    mpz_clear(tmp2);

    return s;
}

// RETURN VALUE: -1 if file exhausted.
//                0 if product is negative...
//                1 if product is positive...
int
treatDep(char *prefix, int numdep, cado_poly pol, char *relname, char *purgedname, char *indexname, char *kername, int verbose)
{
    FILE *depfile, *relfile, *indexfile, *kerfile, *purgedfile = NULL;
    char depname[200], str[1024];
    hashtable_t H;
    relation_t rel;
    uint64_t w;
    int i, j, ret, nlimbs, nrows, ncols, small_nrows, small_ncols;
    int nrel, r, irel, sg, ind;
    char *small_row_used, *rel_used;
    int *vec; // useful to check dependency relation in the purged matrix
    int need64 = (pol->lpba > 32) || (pol->lpbr > 32);
    int rc;

    purgedfile = gzip_open(purgedname, "r");
    rc = fscanf(purgedfile, "%d %d", &nrows, &ncols);
    ASSERT_ALWAYS(rc == 2);

    indexfile = gzip_open(indexname, "r");
    rc = fscanf(indexfile, "%d %d", &small_nrows, &small_ncols);
    ASSERT_ALWAYS(rc == 2);

    nlimbs = ((small_nrows - 1) / 64) + 1;
    // first read used rows in the small matrix
    small_row_used = (char *)malloc(small_nrows * sizeof(char));
    rel_used = (char *)malloc(nrows * sizeof(char));

    // skip first numdep-1 relations
    kerfile = fopen(kername, "r");
    for(j = 0; j < numdep; j++)
	for(i = 0; i < nlimbs; ++i)
	    ret = fscanf (kerfile, "%" SCNx64, &w);

    // use a hash table to rebuild P2
    hashInit(&H, nrows, 1, need64);
	sprintf(depname, "%s.%03d", prefix, numdep);

    memset(small_row_used, 0, small_nrows * sizeof(char));
    // now use this dep
    for(i = 0; i < nlimbs; ++i){
		ret = fscanf (kerfile, "%" SCNx64, &w);
		if(ret == -1)
	    	return ret;
		ASSERT (ret == 1);
		if (verbose)
	    	fprintf(stderr, "w=%" PRIx64 "\n", w);
		for(j = 0; j < 64; ++j){
	    	if(w & 1UL){
				ind = (i * 64) + j;
			if(verbose)
		    	fprintf(stderr, "+R_%d\n", ind);
			small_row_used[ind] = 1;
	    	}
	    w >>= 1;
		}
    }
	fclose(kerfile);
    // now map to the rels of the purged matrix
    memset(rel_used, 0, nrows * sizeof(char));
    for(i = 0; i < small_nrows; i++){
	rc = fscanf(indexfile, "%d", &nrel);
        ASSERT_ALWAYS(rc == 1);
	for(j = 0; j < nrel; j++){
	    rc = fscanf(indexfile, PURGE_INT_FORMAT, &r);
            ASSERT_ALWAYS(rc == 1);
	    if(small_row_used[i]){
#if DEBUG >= 1
		fprintf(stderr, "# Small[%d] -> %d\n", i, r);
#endif
#if DEBUG >= 1
		if(rel_used[r])
		    fprintf(stderr, "WARNING: flipping rel_used[%d]\n", r);
#endif
		rel_used[r] ^= 1;
	    }
	}
    }
    gzip_close(indexfile, indexname);
#if MAPLE >= 1
	fprintf(stderr, "A2:=1;\n");
	fprintf(stderr, "AX:=1;\n");
	fprintf(stderr, "NORM:=1;\n");
#endif
    // FIXME: sg should be removed...
    sg = 1;
    depfile = fopen(depname, "w");
    // now really read the purged matrix in
    vec = (int *)malloc(nrows * sizeof(int));
    char * rp;
    rp = fgets(str, 1024, purgedfile); // get rid of end of first line
    ASSERT_ALWAYS(rp);


    // we assume purgedfile is stored in increasing order of the indices
    // of the real relations, so that one pass in the rels file is needed...!
    relfile = gzip_open(relname, "r");
    irel = 0; // we are ready to read relation irel
    for(i = 0; i < nrows; i++){
        rp = fgets(str, 1024, purgedfile);
        ASSERT_ALWAYS(rp);
        sscanf(str, "%d", vec+i);
        if(rel_used[i]){
#if DEBUG >= 1
            fprintf(stderr, "Reading in rel %d of index %d\n", i, vec[i]);
#endif
            skip_relations_in_file(relfile, vec[i] - irel);
            fread_relation(relfile, &rel);
            irel = vec[i]+1;
            treatRelation(depfile, &H, rel);
            sg *= treatSign(rel, pol);
            clear_relation(&rel);
        }
    }
    gzip_close(relfile, relname);
    //ASSERT(checkVector(vec, ncols));
    if(sg == -1)
	fprintf(stderr, "prod(a-b*m) < 0\n");
    fclose(depfile);
	gzip_close(purgedfile, purgedname);

	fprintf(stderr, "# Treated dependency #%d at %2.2lf\n",
		numdep, seconds());
	hashClear(&H);
    hashFree(&H);
    free(small_row_used);
    free(rel_used);
    free(vec);
    return (sg == -1 ? 0 : 1);
}



/********** RATSQRT **********/

/* Returns memory usage, in KB 
 * This is the VmSize field in the status file of /proc/pid/ dir
 * This is highly non portable.
 * Return -1 in case of failure.
 */
static long Memusage() {
  pid_t pid = getpid();

  char str[1024];
  char *truc;
  snprintf(str, 1024, "/proc/%d/status", pid);

  FILE *file;
  file = fopen(str, "r");
  if (file == NULL)
    return -1;

  long mem;
  for(;;) {
    truc = fgets(str, 1023, file);
    if (truc == NULL) {
      fclose(file);
      return -1;
    }
    int ret = sscanf(str, "VmSize: %ld", &mem);
    if (ret == 1) {
      fclose(file);
      return mem;
    }
  }
}

/* same as above, for resident memory (column RES of top) */
static long Memusage2() {
  pid_t pid = getpid();

  char str[1024];
  char *truc;
  snprintf(str, 1024, "/proc/%d/status", pid);

  FILE *file;
  file = fopen(str, "r");
  if (file == NULL)
    return -1;

  long mem;
  for(;;) {
    truc = fgets(str, 1023, file);
    if (truc == NULL) {
      fclose(file);
      return -1;
    }
    int ret = sscanf(str, "VmRSS: %ld", &mem);
    if (ret == 1) {
      fclose(file);
      return mem;
    }
  }
}

/* Returns peak memory usage, in KB 
 * This is the VmPeak field in the status file of /proc/pid/ dir
 * This is highly non portable.
 * Return -1 in case of failure.
 */
static long PeakMemusage() {
  pid_t pid = getpid();

  char str[1024];
  char *truc;
  snprintf(str, 1024, "/proc/%d/status", pid);

  FILE *file;
  file = fopen(str, "r");
  if (file == NULL)
    return -1;

  long mem;
  for(;;) {
    truc = fgets(str, 1023, file);
    if (truc == NULL) {
      fclose(file);
      return -1;
    }
    int ret = sscanf(str, "VmPeak: %ld", &mem);
    if (ret == 1) {
      fclose(file);
      return mem;
    }
  }
}

static void
my_mpz_mul (mpz_t a, mpz_t b, mpz_t c)
{
  int large, st;

  large = mpz_size (b) + mpz_size (c) >= 5000000;
  if (large)
    {
      fprintf (stderr, "[multiplying %lu*%lu limbs: ",
               mpz_size (b), mpz_size (c));
      fflush (stderr);
      st = cputime ();
    }
  mpz_mul (a, b, c);
  mpz_realloc2 (c, 0);
  if (large)
    {
      fprintf (stderr, "%dms]\n", cputime () - st);
      fflush (stderr);
    }
}

#define THRESHOLD 2 /* must be >= 2 */

/* accumulate up to THRESHOLD products in prd[0], 2^i*THRESHOLD in prd[i].
   nprd is the number of already accumulated values: if nprd = n0 + 
   n1 * THRESHOLD + n2 * THRESHOLD^2 + ..., then prd[0] has n0 entries,
   prd[1] has n1*THRESHOLD entries, and so on.
*/
static mpz_t*
accumulate_fast (mpz_t *prd, mpz_t a, unsigned long *lprd, unsigned long nprd)
{
  unsigned long i;

  my_mpz_mul (prd[0], prd[0], a);
  nprd ++;

  for (i = 0; nprd % THRESHOLD == 0; i++, nprd /= THRESHOLD)
    {
      /* need to access prd[i + 1], thus i+2 entries */
      if (i + 2 > *lprd)
        {
          lprd[0] ++;
          prd = (mpz_t*) realloc (prd, *lprd * sizeof (mpz_t));
          mpz_init_set_ui (prd[i + 1], 1);
        }
      my_mpz_mul (prd[i + 1], prd[i + 1], prd[i]);
      mpz_set_ui (prd[i], 1);
    }

  return prd;
}

/* prd[0] <- prd[0] * prd[1] * ... * prd[lprd-1] */
static void
accumulate_fast_end (mpz_t *prd, unsigned long lprd)
{
  unsigned long i;

  for (i = 1; i < lprd; i++)
    my_mpz_mul (prd[0], prd[0], prd[i]);
}

static size_t
stats (mpz_t *prd, unsigned long lprd)
{
  unsigned long i;
  size_t s = 0;

  for (i = 0; i < lprd; i++)
    s += mpz_size (prd[i]);
  return s;
}

int 
calculateSqrtRat (char *prefix, int numdep, cado_poly pol)
{
  char depname[200];
  char ratname[200];
  sprintf(depname, "%s.%03d", prefix, numdep);
  sprintf(ratname, "%s.rat.%03d", prefix, numdep);
  FILE *depfile = NULL;
  FILE *ratfile;
  //int sign;
  long a, b;
  int ret;
  unsigned long ab_pairs = 0, line_number, freerels = 0;
  mpz_t v, *prd;
  unsigned long lprd; /* number of elements in prd[] */
  unsigned long nprd; /* number of accumulated products in prd[] */
  unsigned long res, peakres = 0;

  if (pol->degreeg != 1)
    {
      fprintf (stderr, "Error, calculateSqrtRat called with non-linear polynomial\n");
      exit (EXIT_FAILURE);
    }

#ifdef __MPIR_VERSION
  fprintf (stderr, "Using MPIR %s\n", mpir_version);
#else
  fprintf (stderr, "Using GMP %s\n", gmp_version);
#endif

  mpz_init (v);

  lprd = 1;
  nprd = 0;
  prd = (mpz_t*) malloc (lprd * sizeof (mpz_t));
  mpz_init_set_ui (prd[0], 1);

  depfile = fopen (depname, "r");
  /*if(!depfile) {
	fprintf (stderr, "Error: file %s not exist\n", depname);
	exit(1);
  }*/
  ASSERT_ALWAYS(depfile != NULL);
	
  line_number = 2;
  for (;;)
    {
      ret = fscanf (depfile, "%ld %ld\n", &a, &b);

      if (ret != 2)
        {
          fprintf (stderr, "Invalid line %lu\n", line_number);
          break;
        }

      ab_pairs ++;
      line_number ++;

      if (ab_pairs % 1000000 == 0)
        {
          res = Memusage2 ();
          if (res > peakres)
            peakres = res;
            fprintf (stderr, "%lu pairs: size %zuMb, %dms, VIRT %luM (peak %luM), RES %luM (peak %luM)\n",
                       ab_pairs, stats (prd, lprd) >> 17, cputime (),
                       Memusage () >> 10, PeakMemusage () >> 10,
                       res >> 10, peakres >> 10);
        }

        if (b == 0)
          freerels ++;

        /* accumulate g1*a+g0*b */
        mpz_mul_si (v, pol->g[1], a);
        mpz_addmul_si (v, pol->g[0], b);

        prd = accumulate_fast (prd, v, &lprd, nprd++);
          
        if (feof (depfile))
          break;
      }
    fprintf (stderr, "%lu (a,b) pairs\n", line_number);

    fclose (depfile);

  fprintf (stderr, "Read %lu (a,b) pairs, including %lu free\n", ab_pairs,
           freerels);

  accumulate_fast_end (prd, lprd);

  /* we must divide by g1^ab_pairs: if the number of (a,b) pairs is odd, we
     multiply by g1, and divide by g1^(ab_pairs+1) */
  if (ab_pairs & 1)
    mpz_mul (prd[0], prd[0], pol->g[1]);

  fprintf (stderr, "Size of product = %zu bits\n", mpz_sizeinbase (prd[0], 2));

  if (mpz_sgn (prd[0]) < 0)
    {
      fprintf (stderr, "Error, product is negative: try another dependency\n");
      exit (1);
    }

  fprintf (stderr, "Starting square root at %dms\n", cputime ());

  /* since we know we have a square, take the square root */
  mpz_sqrtrem (prd[0], v, prd[0]);
  
  fprintf (stderr, "Computed square root at %dms\n", cputime ());

  if (mpz_cmp_ui (v, 0) != 0)
    {
      unsigned long p = 2, e;
      mpz_t pp;

      mpz_init (pp);
      fprintf (stderr, "Error, square root remainder is not zero\n");
      /* reconstruct the initial value of prd[0] to debug */
      mpz_mul (prd[0], prd[0], prd[0]);
      mpz_add (prd[0], prd[0], v);
      while (mpz_cmp_ui (prd[0], 1) > 0)
        {
          e = 0;
          printf ("Removing p=%lu:", p);
          mpz_set_ui (pp, p);
          e = mpz_remove (prd[0], prd[0], pp);
          printf (" exponent=%lu, remaining %lu bits\n", e,
                  mpz_sizeinbase (prd[0], 2));
          if ((e % 2) != 0)
            {
              fprintf (stderr, "Prime %lu appears to odd power %lu\n", p, e);
              break;
            }
          p = getprime (p);
        }
      mpz_clear (pp);
      p = getprime (0);
      exit (1);
    }

  mpz_mod (prd[0], prd[0], pol->n);

  fprintf (stderr, "Reduced mod n at %dms\n", cputime ());

  /* now divide by g1^(ab_pairs/2) if ab_pairs is even, and g1^((ab_pairs+1)/2)
     if ab_pairs is odd */
  
  mpz_powm_ui (v, pol->g[1], (ab_pairs + 1) / 2, pol->n);
  fprintf (stderr, "Computed g1^(nab/2) mod n at %dms\n", cputime ());

  mpz_invert (v, v, pol->n);
  mpz_mul (prd[0], prd[0], v);
  mpz_mod (prd[0], prd[0], pol->n);
  
  ratfile = fopen (ratname, "w");
  gmp_fprintf (ratfile, "%Zd\n", prd[0]);
  fclose(ratfile);

  gmp_fprintf (stderr, "rational square root is %Zd\n", prd[0]);

  fprintf (stderr, "Rational square root time at %dms\n", cputime ());

  mpz_clear (prd[0]);

  mpz_clear (v);
  return 0;
}



/********** ALGSQRT **********/
static void
polymodF_from_ab(polymodF_t tmp, long a, unsigned long b) {
  tmp->v = 0;
  tmp->p->deg = (b != 0) ? 1 : 0;
  mpz_set_si (tmp->p->coeff[1], - (long) b);
  mpz_set_si (tmp->p->coeff[0], a);
}

/* Reduce the coefficients of R in [-m/2, m/2], which are assumed in [0, m[ */
static void
poly_mod_center (poly_t R, const mpz_t m)
{
  int i;
  mpz_t m_over_2;

  mpz_init (m_over_2);
  mpz_div_2exp (m_over_2, m, 2);
  for (i=0; i <= R->deg; i++)
    {
      ASSERT_ALWAYS(mpz_cmp_ui (R->coeff[i], 0) >= 0);
      ASSERT_ALWAYS(mpz_cmp (R->coeff[i], m) < 0);
      if (mpz_cmp (R->coeff[i], m_over_2) > 0)
        mpz_sub (R->coeff[i], R->coeff[i], m);
    }
  mpz_clear (m_over_2);
}

#if 0
/* Check whether the coefficients of R (that are given modulo m) are in
   fact genuine integers. We assume that x mod m is a genuine integer if
   x or |x-m| is less than m/10^6, i.e., the bit size of x or |x-m| is
   less than that of m minus 20.
   Assumes the coefficients x satisfy 0 <= x < m.
*/
static int
poly_integer_reconstruction (poly_t R, const mpz_t m)
{
  int i;
  size_t sizem = mpz_sizeinbase (m, 2), sizer;

  for (i=0; i <= R->deg; ++i)
    {
      sizer = mpz_sizeinbase (R->coeff[i], 2);
      if (sizer + 20 > sizem)
        {
          mpz_sub (R->coeff[i], R->coeff[i], m);
          sizer = mpz_sizeinbase (R->coeff[i], 2);
          if (sizer + 20 > sizem)
            return 0;
        }
    }
  return 1;
}
#endif

// compute res := sqrt(a) in Fp[x]/f(x)
static void
TonelliShanks (poly_t res, const poly_t a, const poly_t F, unsigned long p)
{
  int d = F->deg;
  mpz_t q;
  poly_t delta;  // a non quadratic residue
  poly_t auxpol;
  mpz_t aux;
  mpz_t t;
  int s;
  mpz_t myp;

  mpz_init_set_ui(myp, p);

  mpz_init(aux);
  mpz_init(q);
  poly_alloc(auxpol, d);
  mpz_ui_pow_ui(q, p, (unsigned long)d);

  // compute aux = (q-1)/2
  // and (s,t) s.t.  q-1 = 2^s*t
  mpz_sub_ui(aux, q, 1);
  mpz_divexact_ui(aux, aux, 2);
  mpz_init_set(t, aux);
  s = 1;
  while (mpz_divisible_2exp_p(t, 1)) {
    s++;
    mpz_divexact_ui(t, t, 2);
  }
  // find a non quadratic residue delta
  {
    poly_alloc(delta, d);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    do {
      int i;
      // pick a random delta
      for (i = 0; i < d; ++i)
	mpz_urandomm(delta->coeff[i], state, myp);
      cleandeg(delta, d-1);
      // raise it to power (q-1)/2
      poly_power_mod_f_mod_ui(auxpol, delta, F, aux, p);
    } while ((auxpol->deg != 0) || (mpz_cmp_ui(auxpol->coeff[0], p-1)!= 0));
    gmp_randclear (state);
  }

  // follow the description of Crandall-Pomerance, page 94
  {
    poly_t A, D;
    mpz_t m;
    int i;
    poly_alloc(A, d);
    poly_alloc(D, d);
    mpz_init_set_ui(m, 0);
    poly_power_mod_f_mod_ui(A, a, F, t, p);
    poly_power_mod_f_mod_ui(D, delta, F, t, p);
    for (i = 0; i <= s-1; ++i) {
      poly_power_mod_f_mod_ui(auxpol, D, F, m, p);
      poly_mul_mod_f_mod_mpz(auxpol, auxpol, A, F, myp, NULL);
      mpz_ui_pow_ui(aux, 2, (s-1-i));
      poly_power_mod_f_mod_ui(auxpol, auxpol, F, aux, p);
      if ((auxpol->deg == 0) && (mpz_cmp_ui(auxpol->coeff[0], p-1)== 0))
	mpz_add_ui(m, m, 1UL<<i);
    }
    mpz_add_ui(t, t, 1);
    mpz_divexact_ui(t, t, 2);
    poly_power_mod_f_mod_ui(res, a, F, t, p);
    mpz_divexact_ui(m, m, 2);
    poly_power_mod_f_mod_ui(auxpol, D, F, m, p);

    poly_mul_mod_f_mod_mpz(res, res, auxpol, F, myp, NULL);
    poly_free(D);
    poly_free(A);
    mpz_clear(m);
  }

  poly_free(auxpol);
  poly_free(delta);
  mpz_clear(q);
  mpz_clear(aux);
  mpz_clear(myp);
  mpz_clear (t);
}

// res <- Sqrt(AA) mod F, using p-adic lifting, at prime p.
void
polymodF_sqrt (polymodF_t res, polymodF_t AA, poly_t F, unsigned long p)
{
  poly_t A, *P;
  int v;
  int d = F->deg;
  int k, lk, target_k, logk, K[32];
  size_t target_size; /* target bit size for Hensel lifting */

  /* The size of the coefficients of the square root of A should be about half
     the size of the coefficients of A. Here is an heuristic argument: let
     K = Q[x]/(f(x)) where f(x) is the algebraic polynomial. The square root
     r(x) might be considered as a random element of K: it is smooth, not far
     from an integer, but except that has no relationship with the coefficients
     of f(x). When we square r(x), we obtain a polynomial with coefficients
     twice as large, before reduction by f(x). The reduction modulo f(x)
     produces A(x), however that reduction should not decrease the size of
     the coefficients. */
  target_size = poly_sizeinbase (AA->p, AA->p->deg, 2);
  target_size = target_size / 2;
  target_size += target_size / 10;
  fprintf (stderr, "target_size=%lu\n", (unsigned long int) target_size);

  poly_alloc(A, d-1);
  // Clean up the mess with denominator: if it is an odd power of fd,
  // then multiply num and denom by fd to make it even.
  if (((AA->v)&1) == 0) {
    v = AA->v / 2;
    poly_copy(A, AA->p);
  } else {
    v = (1+AA->v) / 2;
    poly_mul_mpz(A, AA->p, F->coeff[d]);
  }

  // Now, we just have to take the square root of A (without denom) and
  // divide by fd^v.

  // Variables for the lifted values
  poly_t invsqrtA;
  // variables for A and F modulo pk
  poly_t a;
  poly_alloc(invsqrtA, d-1);
  poly_alloc(a, d-1);
  // variable for the current pk
  mpz_t pk, invpk;
  mpz_init (pk);
  mpz_init (invpk);

  /* Jason Papadopoulos's trick: since we will lift the square root of A to at
     most target_size bits, we can reduce A accordingly */
  double st = seconds ();
  target_k = (int) ((double) target_size * log ((double) 2) / log((double) p));
  mpz_ui_pow_ui (pk, p, target_k);
  while (mpz_sizeinbase (pk, 2) <= target_size)
    {
      mpz_mul_ui (pk, pk, p);
      target_k ++;
    }
  poly_reduce_mod_mpz (A, A, pk);
  for (k = target_k, logk = 0; k > 1; k = (k + 1) / 2, logk ++)
    K[logk] = k;
  K[logk] = 1;
  fprintf (stderr, "Reducing A mod p^%d took %2.2lf\n", target_k,
           seconds () - st);

  // Initialize things modulo p:
  mpz_set_ui (pk, p);
  k = 1; /* invariant: pk = p^k */
  lk = 0; /* k = 2^lk */
  st = seconds ();
  P = poly_base_modp_init (A, p, K, logk);
  fprintf (stderr, "poly_base_modp_init took %2.2lf\n", seconds () - st);

  poly_copy (a, P[0]);

  // First compute the inverse square root modulo p
  {
    mpz_t q, aux;
    mpz_init(q);
    mpz_init(aux);
    mpz_ui_pow_ui(q, p, (unsigned long)d);

#if 0
    // compute (q-2)(q+1)/4   (assume q == 3 mod 4, here !!!!!)
    // where q = p^d, the size of the finite field extension.
    // since we work mod q-1, this gives (3*q-5)/4
    mpz_mul_ui(aux, q, 3);
    mpz_sub_ui(aux, aux, 5);
    mpz_divexact_ui(aux, aux, 4);		        // aux := (3q-5)/4
    poly_power_mod_f_mod_ui(invsqrtA, a, F, aux, p);
#else
    TonelliShanks(invsqrtA, a, F, p);
    mpz_sub_ui(aux, q, 2);
    poly_power_mod_f_mod_ui(invsqrtA, invsqrtA, F, aux, p);
#endif

    mpz_clear(aux);
    mpz_clear(q);
  }

  // Now, the lift begins
  // When entering the loop, invsqrtA contains the inverse square root
  // of A computed modulo pk.

  poly_t tmp, tmp2;
  poly_alloc(tmp, 2*d-1);
  poly_alloc(tmp2, 2*d-1);
  do {
    double st;

    if (mpz_sizeinbase (pk, 2) > target_size)
      {
        fprintf (stderr, "Failed to reconstruct an integer polynomial\n");
        printf ("Failed\n");
        exit (1);
      }

    /* invariant: invsqrtA = 1/sqrt(A) bmod p^k */

    st = seconds ();
    poly_base_modp_lift (a, P, ++lk, pk);
    st = seconds () - st;

    /* invariant: k = K[logk] */
    ASSERT_ALWAYS(k == K[logk]);

    mpz_set (invpk, pk);
    mpz_mul (pk, pk, pk);	// double the current precision
    logk --;
    if (K[logk] & 1)
      {
        mpz_div_ui (pk, pk, p);
        k --;
      }
    k = K[logk];
    barrett_init (invpk, pk); /* FIXME: we could lift 1/p^k also */
    fprintf (stderr, "start lifting mod p^%d (%lu bits) at %2.2lf\n",
             k, (unsigned long int) mpz_sizeinbase (pk, 2), seconds ());
#ifdef VERBOSE
    fprintf (stderr, "   poly_base_modp_lift took %2.2lf\n", st);
#endif

    // now, do the Newton operation x <- 1/2(3*x-a*x^3)
    st = seconds ();
    poly_sqr_mod_f_mod_mpz(tmp, invsqrtA, F, pk, invpk); /* tmp = invsqrtA^2 */
#ifdef VERBOSE
    fprintf (stderr, "   poly_sqr_mod_f_mod_mpz took %2.2lf\n", seconds () - st);
#endif

    /* Faster version which computes x <- x + x/2*(1-a*x^2).
       However I don't see how to use the fact that the coefficients
       if 1-a*x^2 are divisible by p^(k/2). */
    st = seconds ();
    poly_mul_mod_f_mod_mpz (tmp, tmp, a, F, pk, invpk); /* tmp=a*invsqrtA^2 */
#ifdef VERBOSE
    fprintf (stderr, "   poly_mul_mod_f_mod_mpz took %2.2lf\n", seconds () - st);
#endif
    poly_sub_ui (tmp, 1); /* a*invsqrtA^2-1 */
    poly_div_2_mod_mpz (tmp, tmp, pk); /* (a*invsqrtA^2-1)/2 */
    st = seconds ();
    poly_mul_mod_f_mod_mpz (tmp, tmp, invsqrtA, F, pk, invpk);
#ifdef VERBOSE
    fprintf (stderr, "   poly_mul_mod_f_mod_mpz took %2.2lf\n", seconds () - st);
#endif
    /* tmp = invsqrtA/2 * (a*invsqrtA^2-1) */
    poly_sub_mod_mpz (invsqrtA, invsqrtA, tmp, pk);

  } while (k < target_k);

  /* multiply by a to get an approximation of the square root */
  poly_mul_mod_f_mod_mpz (tmp, invsqrtA, a, F, pk, invpk);
  poly_mod_center (tmp, pk);

  poly_base_modp_clear (P);

  poly_copy(res->p, tmp);
  res->v = v;

  mpz_clear (pk);
  mpz_clear (invpk);
  poly_free(tmp);
  poly_free(tmp2);
  poly_free (A);
  poly_free (invsqrtA);
  poly_free (a);

  size_t sqrt_size = poly_sizeinbase (res->p, F->deg - 1, 2);
  fprintf (stderr, "maximal sqrt bit-size = %zu (%.0f%% of target size)\n",
          sqrt_size, 100.0 * (double) sqrt_size / target_size);
}

static unsigned long
FindSuitableModP (poly_t F)
{
  unsigned long p = 2;
  int dF = F->deg;

  plain_poly_t fp;

  plain_poly_init (fp, dF);
  while (1)
    {
    int d;

    p = getprime (p);
    if (! plain_poly_fits (dF, p))
      {
        fprintf (stderr, "You are in trouble. Please contact the CADO support team at cado-nfs-commits@lists.gforge.inria.fr.\n");
        exit (1);
      }
    d = plain_poly_set_mod (fp, F->coeff, dF, p);
    if (d != dF)
      continue;
    if (plain_poly_is_irreducible (fp, p))
      break;
    }
  plain_poly_clear(fp);
  getprime (0);

  return p;
}

// Products are computed modulo the polynomial F.
polymodF_t*
accumulate_fast_F (polymodF_t *prd, const polymodF_t a, const poly_t F,
                 unsigned long *lprd, unsigned long nprd)
{
  unsigned long i;

  polymodF_mul (prd[0], prd[0], a, F);
  nprd ++;

  for (i = 0; nprd % THRESHOLD == 0; i++, nprd /= THRESHOLD)
    {
      /* need to access prd[i + 1], thus i+2 entries */
      if (i + 2 > *lprd)
        {
          lprd[0] ++;
          prd = (polymodF_t*) realloc (prd, *lprd * sizeof (polymodF_t));
	  poly_alloc(prd[i+1]->p, F->deg);
          mpz_set_ui(prd[i + 1]->p->coeff[0], 1);
	  prd[i+1]->p->deg = 0;
	  prd[i+1]->v = 0;
        }
      polymodF_mul (prd[i+1], prd[i+1], prd[i], F);
      mpz_set_ui(prd[i]->p->coeff[0], 1);
      prd[i]->p->deg = 0;
      prd[i]->v = 0;
    }

  return prd;
}

/* prd[0] <- prd[0] * prd[1] * ... * prd[lprd-1] */
void
accumulate_fast_F_end (polymodF_t *prd, const poly_t F, unsigned long lprd)
{
  unsigned long i;

  for (i = 1; i < lprd; i++)
    polymodF_mul (prd[0], prd[0], prd[i], F);
}

/* side=0: consider the polynomial f
   side=1: consider the polynomial g
*/
int
calculateSqrtAlg (char *prefix, int numdep, cado_poly pol, int side)
{
  char depname[200];
  char algname[200];
  FILE *depfile = NULL;
  FILE *algfile;
  poly_t F;
  polymodF_t prd, tmp;
  long a;
  unsigned long b;
  unsigned long p;
  double t0 = seconds ();
  mpz_t algsqrt, aux;
  int i, deg;
  mpz_t *f;
  int nab = 0, nfree = 0;

  ASSERT_ALWAYS(side == 0 || side == 1);

  sprintf (depname, "%s.%03d", prefix, numdep);
  if (side == 0)
    sprintf (algname, "%s.alg.%03d", prefix, numdep);
  else
    sprintf (algname, "%s.rat.%03d", prefix, numdep);
  depfile = fopen (depname, "r");
  ASSERT_ALWAYS(depfile != NULL);

  deg = (side == 0) ? pol->degree : pol->degreeg;
  f = (side == 0) ? pol->f : pol->g;

  ASSERT_ALWAYS(deg > 1);
  
  // Init F to be the corresponding polynomial
  poly_alloc (F, deg);
  for (i = deg; i >= 0; --i)
    poly_setcoeff (F, i, f[i]);
  
  // Init prd to 1.
  poly_alloc (prd->p, deg);
  mpz_set_ui (prd->p->coeff[0], 1);
  prd->p->deg = 0;
  prd->v = 0;
  
  // Allocate tmp
  poly_alloc (tmp->p, 1);
  
  // Accumulate product
  #if 0
    // Naive version, without subproduct tree
    while(fscanf(depfile, "%ld %lu", &a, &b) != EOF){
      if(!(nab % 100000))
        fprintf(stderr, "# Reading ab pair #%d at %2.2lf\n",nab,seconds());
      if((a == 0) && (b == 0))
        break;
      polymodF_from_ab (tmp, a, b);
      polymodF_mul (prd, prd, tmp, F);
      nab++;
    }
  #else
    // With a subproduct tree
    {
      polymodF_t *prd_tab;
      unsigned long lprd = 1; /* number of elements in prd_tab[] */
      unsigned long nprd = 0; /* number of accumulated products in prd_tab[] */
      prd_tab = (polymodF_t*) malloc (lprd * sizeof (polymodF_t));
      poly_alloc (prd_tab[0]->p, F->deg);
      mpz_set_ui (prd_tab[0]->p->coeff[0], 1);
      prd_tab[0]->p->deg = 0;
      prd_tab[0]->v = 0;
      while(fscanf(depfile, "%ld %lu", &a, &b) != EOF){
        if(!(nab % 100000))
  	fprintf(stderr, "# Reading ab pair #%d at %2.2lf\n", nab, seconds ());
        if((a == 0) && (b == 0))
  	break;
        polymodF_from_ab(tmp, a, b);
        prd_tab = accumulate_fast_F (prd_tab, tmp, F, &lprd, nprd++);
        nab++;
        if(b == 0)
  	  nfree++;
      }
      fprintf (stderr, "# Read %d including %d free relations\n", nab, nfree);
      ASSERT_ALWAYS ((nab & 1) == 0);
      ASSERT_ALWAYS ((nfree & 1) == 0);
      /* nfree being even is forced by a specific character column added
       * by character.c. The correspond assert should not fail.
       *
       * nab being even is kind of a mystery. There is a character column
       * that gives the sign of the rational side. It could well be that
       * with the parameters we usually use, it is negative for all
       * relations; that would force the overall number of relations to
       * be even. Another possibility is when f_d contains a big prime
       * number that does not occur anywhere else, so that the power of
       * this prime is controlled by the usual characters, and since f_d
       * is always present... 
       *
       * But: I wouldn't be surprised if the assert(even(nab)) fails.
       * Then, a patch would be:
       *    - remove the assert (!)
       *    - in the numerator of the big product, eliminate powers of
       *       f_d that divides all coefficients.
       *    - this should finally give an even power of f_d in the
       *      denominator, and the algorithm can continue.
       */
      accumulate_fast_F_end (prd_tab, F, lprd);
      fclose(depfile);
  
      poly_copy(prd->p, prd_tab[0]->p);
      prd->v = prd_tab[0]->v;
      for (i = 0; i < (long)lprd; ++i)
        poly_free(prd_tab[i]->p);
      free(prd_tab);
    }
  #endif
  
    fprintf(stderr, "Finished accumulating the product at %2.2lf\n", seconds());
    fprintf(stderr, "nab = %d, nfree = %d, v = %d\n", nab, nfree, prd->v);
    fprintf (stderr, "maximal polynomial bit-size = %lu\n",
             (unsigned long) poly_sizeinbase (prd->p, deg - 1, 2));
  
    p = FindSuitableModP(F);
    fprintf(stderr, "Using p=%lu for lifting\n", p);
  
    double tm = seconds();
    polymodF_sqrt (prd, prd, F, p);
    fprintf (stderr, "Square root lifted in %2.2lf\n", seconds()-tm);
  
    mpz_init(algsqrt);
    mpz_init(aux);
    poly_eval_mod_mpz(algsqrt, prd->p, pol->m, pol->n);
    mpz_invert(aux, F->coeff[F->deg], pol->n);  // 1/fd mod n
    mpz_powm_ui(aux, aux, prd->v, pol->n);      // 1/fd^v mod n
    mpz_mul(algsqrt, algsqrt, aux);
    mpz_mod(algsqrt, algsqrt, pol->n);
  
	algfile = fopen (algname, "w");
	gmp_fprintf (algfile, "%Zd\n", algsqrt);
	fclose(algfile);

    gmp_fprintf(stderr, "algebraic square root is: %Zd\n", algsqrt);
    fprintf (stderr, "Algebraic square root time is %2.2lf\n", seconds() - t0);
    mpz_clear(aux);
    mpz_clear(algsqrt);
    poly_free(prd->p);
    poly_free(tmp->p);
    poly_free(F);
	return 0;
}



/********** GCD **********/
int
calculateGcd(char *prefix, int numdep, cado_poly pol)
{
	char ratname[200];
	char algname[200];
    sprintf(ratname, "%s.rat.%03d", prefix, numdep);
    sprintf(algname, "%s.alg.%03d", prefix, numdep);
	FILE *ratfile = NULL;
	FILE *algfile = NULL;
    int found = 0;
    mpz_t ratsqrt, algsqrt, g1, g2;

    mpz_init(ratsqrt);
    mpz_init(algsqrt);
    mpz_init(g1);
    mpz_init(g2);
  
	ratfile = fopen(ratname, "r");
	algfile = fopen(algname, "r");
	ASSERT_ALWAYS(ratfile != NULL);
	ASSERT_ALWAYS(algfile != NULL);
	gmp_fscanf(ratfile, "%Zd", ratsqrt);
	gmp_fscanf(algfile, "%Zd", algsqrt);
	fclose(ratfile);
	fclose(algfile);

    // First check that the squares agree
    mpz_mul(g1, ratsqrt, ratsqrt);
    mpz_mod(g1, g1, pol->n);

    mpz_mul(g2, algsqrt, algsqrt);
    mpz_mod(g2, g2, pol->n);

    if (mpz_cmp(g1, g2)!=0) {
      fprintf(stderr, "Bug: the squares do not agree modulo n!\n");
	  exit(1);
      //      gmp_printf("g1:=%Zd;\ng2:=%Zd;\n", g1, g2);
    }

    mpz_sub(g1, ratsqrt, algsqrt);
    mpz_gcd(g1, g1, pol->n);
    if (mpz_cmp(g1,pol->n)) {
      if (mpz_cmp_ui(g1,1)) {
        found = 1;
        gmp_printf("%Zd\n", g1);
      }
    }

    mpz_add(g2, ratsqrt, algsqrt);
    mpz_gcd(g2, g2, pol->n);
    if (mpz_cmp(g2,pol->n)) {
      if (mpz_cmp_ui(g2,1)) {
        found = 1;
        gmp_printf("%Zd\n", g2);
      }
    }
    mpz_clear(g1);
    mpz_clear(g2);

    if (!found)
      printf("Failed\n");
  
    mpz_clear(ratsqrt);
    mpz_clear(algsqrt);
  
    return 0;
}


int main(int argc, char *argv[])
{
	char *polyname = NULL, *prefix = NULL;
    char *relname = NULL, *purgedname = NULL, *indexname = NULL, *kername = NULL;
    cado_poly pol;
    int verbose = 0, numdep = -1, opt = 0, ret, i;

    /* print the command line */
    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

	if (argc == 1 || (argc > 1 && strcmp(argv[1], "--help") == 0 )) {
		fprintf(stderr, "Usage: %s [-ab || -rat || -alg || -gcd] -poly polyname -dep prefix numdep", argv[0]);
		fprintf(stderr, " -rel relname -purged purgedname -index indexname -ker kername\n");
		fprintf(stderr, "or %s (-rat || -alg || -gcd) -poly polyname -dep prefix numdep\n\n", argv[0]);
		fprintf(stderr, "(a,b) pairs of dependency relation 'numdep' will be r/w in file 'prefix.numdep',");
		fprintf(stderr, " rational sqrt in 'prefix.rat.numdep' ...\n");
		exit(1);
	}
		
	/* args */
    while (argc > 1 && argv[1][0] == '-') {
		if(argc > 1 && strcmp(argv[1], "-ab") == 0) {
			opt += 1;
			argc--;
			argv++;
		}
		else if(argc > 1 && strcmp(argv[1], "-rat") == 0) {
			opt += 2;
			argc--;
			argv++;
		}
		else if(argc > 1 && strcmp(argv[1], "-alg") == 0) {
			opt += 4;
			argc--;
			argv++;
		}
		else if(argc > 1 && strcmp(argv[1], "-gcd") == 0) {
			opt += 8;
			argc--;
			argv++;
		}
        else if (argc > 2 && strcmp(argv[1], "-poly") == 0) {
	    	polyname = argv[2];
            argc -= 2;
            argv += 2;
        }
        else if (argc > 3 && strcmp(argv[1], "-dep") == 0) {
	    	prefix = argv[2];
	    	numdep = atoi(argv[3]);
            argc -= 3;
            argv += 3;
        }
        else if (argc > 2 && strcmp(argv[1], "-rel") == 0) {
	    	relname = argv[2];
            argc -= 2;
            argv += 2;
        }
        else if (argc > 2 && strcmp(argv[1], "-purged") == 0) {
	    	purgedname = argv[2];
            argc -= 2;
            argv += 2;
        }
        else if (argc > 2 && strcmp(argv[1], "-index") == 0) {
	    	indexname = argv[2];
            argc -= 2;
            argv += 2;
        }
        else if (argc > 2 && strcmp(argv[1], "-ker") == 0) {
	    	kername = argv[2];
            argc -= 2;
            argv += 2;
        }
		else {
			fprintf(stderr, "Error: %s parameter unknown or bad arguments in command line\n", argv[1]);
			exit(1);
		}
	}

	/* if no options then -ab -rat -alg -gcd */
	if (!opt)
		opt = 15;

    ASSERT_ALWAYS(polyname != NULL);
    ASSERT_ALWAYS(prefix != NULL);
    ASSERT_ALWAYS(numdep != -1);

    cado_poly_init(pol);
    ret = cado_poly_read(pol, polyname);
    ASSERT (ret);

    /* if bit 0 of opt is set: compute the (a,b) pairs
       if bit 1 of opt is set: compute the rational square root
       if bit 2 of opt is set: compute the algebraic square root
       if bit 3 of opt is set: compute the gcd */
           
    int count = 0;
    while (opt)
      {
        if (opt%2 == 1)
          {
            if (count == 0)
              { /* compute the (a,b) pairs */
                ASSERT_ALWAYS(relname != NULL);
                ASSERT_ALWAYS(purgedname != NULL);
                ASSERT_ALWAYS(indexname != NULL);
                ASSERT_ALWAYS(kername != NULL);
                treatDep (prefix, numdep, pol, relname, purgedname, indexname,
                          kername, verbose);
              } 
            else if (count == 1) /* compute the square root on the g-side */
              {
                if (pol->degreeg == 1)
                  calculateSqrtRat (prefix, numdep, pol);
                else
                  calculateSqrtAlg (prefix, numdep, pol, 1);
              }
            else if (count == 2) /* compute the square root on the f-side */
              calculateSqrtAlg (prefix, numdep, pol, 0);
            else 
              calculateGcd (prefix, numdep, pol);
            opt--;
          }
        opt = opt / 2;
        count++;
      }
	
    cado_poly_clear (pol);
    return 0;
}
  
