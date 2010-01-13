#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <inttypes.h>
#include <math.h> /* for log */

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

// str is a row coming from purgedfile
void
str2Vec(int *vec, char *str)
{
    char *t = str;
    int j, k = 0, nc;

    // skip first integer which is the index
    for(; *t != ' '; t++);
#if 1
    sscanf(t+1, "%d", &nc);
    for(j = 0; j < nc; j++){
	// t is on a "space"
	t++;
	// t begins on a number
	sscanf(t, PURGE_INT_FORMAT, &k);
	vec[k]++;
	if(j == (nc-1))
	    // no need to read further
	    break;
	// go to end of number
	for(; *t != ' '; t++);
    }
#else
    // skip second integer which is the number of primes
    for(++t; *t != ' '; t++);
    t++;
    while(1){
	if((*t == '\n') || (*t == ' ')){
	    // new integer read
	    vec[k]++;
	    k = 0;
	    if(*t == '\n')
		break;
	}
	else
	    k = k*10+(*t - '0');
	t++;
    }
#endif
}

void
treatRationalRelation(hashtable_t *H, relation_t rel)
{
    int j, h;
    uint64_t minus2 = (H->need64) ? (uint64_t) (-2) : (uint32_t) (-2);

#if MAPLE >= 1
    fprintf(stderr, "P:=P * (%ld-m*%lu):\n", rel.a, rel.b);
#endif
    for(j = 0; j < rel.nb_rp; j++){
#if MAPLE >= 1
	fprintf(stderr, "R2:=R2*%ld^%d:\n", rel.rp[j].p, rel.rp[j].e);
#endif
	h = getHashAddr(H, rel.rp[j].p, minus2);
	if(H->hashcount[h] == 0){
	    // new empty place
            SET_HASH_P(H,h,rel.rp[j].p);
            SET_HASH_R(H,h,minus2);
	}
	H->hashcount[h] += rel.rp[j].e;
    }
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
treatAlgebraicRelation(FILE *algfile, hashtable_t *H, relation_t rel)
{
    int j, h;

    if(rel.b > 0)
	computeroots(&rel);
    else
	computefree(&rel);
    fprintf(algfile, "%ld %lu\n", rel.a, rel.b);
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

void
finishRationalSqrt(FILE *ratfile, hashtable_t *H, cado_poly pol)
{
    mpz_t prod;
    int j;
    unsigned int i;
    uint64_t minus2 = (H->need64) ? (uint64_t) (-2) : (uint32_t) (-2);

    mpz_init_set_ui(prod, 1);
#if MAPLE >= 1
    fprintf(stderr, "R:=1;\n");
#endif
    for(i = 0; i < H->hashmod; i++){
	if(H->hashcount[i] > 0){
          if (GET_HASH_R(H,i) != minus2)
		continue;
	    if ((H->hashcount[i] & 1)) {
	        fprintf(stderr, "  Odd valuation! At rational prime %"PRIi64"\n",
                        GET_HASH_P(H,i));
		exit(1);
	    }
#if MAPLE >= 1
	    fprintf(stderr, "R:=R*%ld^%d:\n",
		    H->hashtab_p[i], H->hashcount[i]>>1);
#endif
	    // TODO: do better
	    for(j = 0; j < (H->hashcount[i]>>1); j++)
              mpz_mul_ui(prod, prod, GET_HASH_P(H,i));
	    mpz_mod(prod, prod, pol->n);
	}
    }
#if DEBUG >= 1
    gmp_fprintf(stderr, "prod:=%Zd;\n", prod);
    fprintf(stderr, "# We print the squareroot of the rational side...\n");
#endif
    // TODO: humf, I should say...!
    gmp_fprintf(ratfile, "%Zd\n", prod);
    mpz_clear(prod);
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
//
// If check == 0, then we don't need to read purgedfile again and again
// so we store the interesing values in the vec array (ohhhhhh!).
int
treatDep(char *ratname, char *algname, char *relname, char *purgedname, char *indexname, FILE *kerfile, cado_poly pol, int nrows, int ncols, char *small_row_used, int small_nrows, hashtable_t *H, int nlimbs, char *rel_used, int *vec, int rora, int verbose, int check)
{
    FILE *ratfile, *algfile, *relfile, *indexfile, *purgedfile = NULL;
    relation_t rel;
    uint64_t w;
    int ret, i, j, nrel, r, irel, nr, sg, ind;
    char str[1024];

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
    // now map to the rels of the purged matrix
    memset(rel_used, 0, nrows * sizeof(char));
    indexfile = gzip_open(indexname, "r");
    fscanf(indexfile, "%d %d", &i, &j); // skip first line
    for(i = 0; i < small_nrows; i++){
	fscanf(indexfile, "%d", &nrel);
	for(j = 0; j < nrel; j++){
	    fscanf(indexfile, PURGE_INT_FORMAT, &r);
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
    if((rora == 1) || (rora == 3))
	fprintf(stderr, "R2:=1; P:=1;\n");
    if((rora == 2) || (rora == 3)){
	fprintf(stderr, "A2:=1;\n");
	fprintf(stderr, "AX:=1;\n");
	fprintf(stderr, "NORM:=1;\n");
    }
#endif
    if(check)
	memset(vec, 0, ncols * sizeof(int));
    // FIXME: sg should be removed...
    sg = 1;
    ratfile = fopen(ratname, "w");
    algfile = fopen(algname, "w");
    // now really read the purged matrix in
    if(check){
	purgedfile = gzip_open(purgedname, "r");
	fgets(str, 1024, purgedfile); // skip first line
    }
    // we assume purgedfile is stored in increasing order of the indices
    // of the real relations, so that one pass in the rels file is needed...!
    relfile = gzip_open(relname, "r");
    irel = 0; // we are ready to read relation irel
    for(i = 0; i < nrows; i++){
	if(check)
	    fgets(str, 1024, purgedfile);
	if(rel_used[i]){
	    if(check)
		sscanf(str, "%d", &nr);
	    else
		nr = vec[i];
#if DEBUG >= 1
	    fprintf(stderr, "Reading in rel %d of index %d\n", i, nr);
#endif
	    if(check)
		str2Vec(vec, str);
            skip_relations_in_file(relfile, nr - irel);
            fread_relation(relfile, &rel);
	    irel = nr+1;
	    if((rora == 1) || (rora == 3))
		treatRationalRelation(H, rel);
	    if((rora == 2) || (rora == 3))
		treatAlgebraicRelation(algfile, H, rel);
	    sg *= treatSign(rel, pol);
	    clear_relation(&rel);
	}
    }
    gzip_close(relfile, relname);
    ASSERT(!check || checkVector(vec, ncols));
    if(sg == -1)
	fprintf(stderr, "prod(a-b*m) < 0\n");
    else{
	if((rora == 1) || (rora == 3))
	    finishRationalSqrt(ratfile, H, pol);
    }
    fclose(ratfile);
    fclose(algfile);
    if(check)
	gzip_close(purgedfile, purgedname);
    return (sg == -1 ? 0 : 1);
}

void
SqrtWithIndexAll(char *prefix, char *relname, char *purgedname, char *indexname, FILE *kerfile, cado_poly pol, int rora, int numdep, int verbose, int check)
{
    FILE *indexfile, *purgedfile;
    char ratname[200], algname[200], str[1024];
    hashtable_t H;
    uint64_t w;
    int i, j, ret, nlimbs, nrows, ncols, small_nrows, small_ncols;
    char *small_row_used, *rel_used;
    int *vec; // useful to check dependency relation in the purged matrix
    int need64 = (pol->lpba > 32) || (pol->lpbr > 32);

    purgedfile = gzip_open(purgedname, "r");
    fscanf(purgedfile, "%d %d", &nrows, &ncols);

    indexfile = gzip_open(indexname, "r");
    fscanf(indexfile, "%d %d", &small_nrows, &small_ncols);
    gzip_close(indexfile, indexname);

    nlimbs = ((small_nrows - 1) / 64) + 1;
    // first read used rows in the small matrix
    small_row_used = (char *)malloc(small_nrows * sizeof(char));
    rel_used = (char *)malloc(nrows * sizeof(char));
    if(check)
	vec = (int *)malloc(ncols * sizeof(int));
    else{
	// ohhhhhhhhhhhhh!
	vec = (int *)malloc(nrows * sizeof(int));
	// get rid of end of first line
	fgets(str, 1024, purgedfile);
	for(i = 0; i < nrows; i++){
	    fgets(str, 1024, purgedfile);
	    sscanf(str, "%d", vec+i);
	}
    }
    gzip_close(purgedfile, purgedname);

    // skip first numdep-1 relations
    for(j = 0; j < numdep; j++)
	for(i = 0; i < nlimbs; ++i)
	    ret = fscanf (kerfile, "%" SCNx64, &w);

    // use a hash table to rebuild P2
    hashInit(&H, nrows, 1, need64);
	sprintf(ratname, "%s.rat.%03d", prefix, numdep);
	sprintf(algname, "%s.alg.%03d", prefix, numdep);
	ret = treatDep(ratname, algname, relname, purgedname, indexname, kerfile, pol, nrows, ncols, small_row_used, small_nrows, &H, nlimbs, rel_used, vec, rora, verbose, check);
	fprintf(stderr, "# Treated dependency #%d at %2.2lf\n",
		numdep, seconds());
	hashClear(&H);
    hashFree(&H);
    free(small_row_used);
    free(rel_used);
    free(vec);
}

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
}

unsigned long FindSuitableModP(poly_t F) {
  unsigned long p = 2;
  int dF = F->deg;

  plain_poly_t fp;
  plain_poly_init(fp, dF);
  while (1) {
    int d;
    // select prime congruent to 3 mod 4 to have easy sqrt.
  //  do {
      p = getprime(p);
  //  } while ((p%4)!= 3);
    d = plain_poly_set_mod (fp, F->coeff, dF, p);
    if (d!=dF)
      continue;
    if (plain_poly_is_irreducible(fp, p))
      break;
  }
  plain_poly_clear(fp);
  getprime (0);

  return p;
}

#define THRESHOLD 2 /* must be >= 2 */

/* accumulate up to THRESHOLD products in prd[0], 2^i*THRESHOLD in prd[i].
   nprd is the number of already accumulated values: if nprd = n0 +
   n1 * THRESHOLD + n2 * THRESHOLD^2 + ..., then prd[0] has n0 entries,
   prd[1] has n1*THRESHOLD entries, and so on.
*/
// Products are computed modulo the polynomial F.
polymodF_t*
accumulate_fast (polymodF_t *prd, const polymodF_t a, const poly_t F,
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
accumulate_fast_end (polymodF_t *prd, const poly_t F, unsigned long lprd)
{
  unsigned long i;

  for (i = 1; i < lprd; i++)
    polymodF_mul (prd[0], prd[0], prd[i], F);
}

int main(int argc, char *argv[])
{
    char *relname, *purgedname, *indexname, *kername, *polyname, *prefix;
    FILE *kerfile;
    cado_poly pol;
    int verbose = 1, numdep, rora, ret, i, check;

	char algname[200];
	char ratname[200];
	FILE * algfile;
	FILE * ratfile;
	poly_t F;
	polymodF_t prd, tmp;
	long a;
    unsigned long b;
    unsigned long p;
    double t0 = seconds ();

    if(argc != 9){
	fprintf(stderr, "Usage: %s relname purgedname indexname", argv[0]);
	fprintf(stderr, " kername polyname numdep r|a prefix\n");
	fprintf(stderr, "Dependency relation i will be put in files prefix.i\n");
	return 0;
    }

    /* print the command line */
    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

    relname = argv[1];
    purgedname = argv[2];
    indexname = argv[3];
    kername = argv[4];
    polyname = argv[5];
    numdep = atoi(argv[6]);
	prefix = argv[8];

    kerfile = fopen(kername, "r");

    cado_poly_init(pol);
    ret = cado_poly_read(pol, polyname);
    ASSERT (ret);

    if(!strcmp(argv[7], "r"))
	rora = 1;
    else if(!strcmp(argv[7], "a"))
	rora = 2;
    else if(!strcmp(argv[7], "ar") || !strcmp(argv[7], "ra"))
	rora = 3;
    else
      {
        fprintf (stderr, "Error, 8th argument must be r, a, ar or ra\n");
        exit (1);
      }
    verbose = 0;
    check = 0;
    SqrtWithIndexAll(prefix, relname, purgedname, indexname, kerfile, pol, rora, numdep, verbose, check);

    fclose(kerfile);

	/********** ALGSQRT **********/
  
    sprintf(algname, "%s.alg.%03d", prefix, numdep);
    algfile = fopen(algname, "r");
    ASSERT_ALWAYS(algfile != NULL);
  
    // Init F to be the algebraic polynomial
    poly_alloc(F, pol->degree);
    for (i = pol->degree; i >= 0; --i)
      poly_setcoeff(F, i, pol->f[i]);
  
    // Init prd to 1.
    poly_alloc(prd->p, pol->degree);
    mpz_set_ui(prd->p->coeff[0], 1);
    prd->p->deg = 0;
    prd->v = 0;
  
    // Allocate tmp
    poly_alloc(tmp->p, 1);
  
    // Accumulate product
    int nab = 0, nfree = 0;
  #if 0
    // Naive version, without subproduct tree
    while(fscanf(algfile, "%ld %lu", &a, &b) != EOF){
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
      poly_alloc(prd_tab[0]->p, F->deg);
      mpz_set_ui(prd_tab[0]->p->coeff[0], 1);
      prd_tab[0]->p->deg = 0;
      prd_tab[0]->v = 0;
      while(fscanf(algfile, "%ld %lu", &a, &b) != EOF){
        if(!(nab % 100000))
  	fprintf(stderr, "# Reading ab pair #%d at %2.2lf\n", nab, seconds ());
        if((a == 0) && (b == 0))
  	break;
        polymodF_from_ab(tmp, a, b);
        prd_tab = accumulate_fast (prd_tab, tmp, F, &lprd, nprd++);
        nab++;
        if(b == 0)
  	  nfree++;
      }
      fprintf (stderr, "# Read %d including %d free relations\n", nab, nfree);
      ASSERT_ALWAYS ((nab & 1) == 0);
      ASSERT_ALWAYS ((nfree & 1) == 0);
      accumulate_fast_end (prd_tab, F, lprd);
      fclose(algfile);
  
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
             (unsigned long) poly_sizeinbase (prd->p, pol->degree - 1, 2));
  
    p = FindSuitableModP(F);
    fprintf(stderr, "Using p=%lu for lifting\n", p);
  
    double tm = seconds();
    polymodF_sqrt (prd, prd, F, p);
    fprintf (stderr, "Square root lifted in %2.2lf\n", seconds()-tm);
  
    mpz_t algsqrt, aux;
    mpz_init(algsqrt);
    mpz_init(aux);
    poly_eval_mod_mpz(algsqrt, prd->p, pol->m, pol->n);
    mpz_invert(aux, F->coeff[F->deg], pol->n);  // 1/fd mod n
    mpz_powm_ui(aux, aux, prd->v, pol->n);      // 1/fd^v mod n
    mpz_mul(algsqrt, algsqrt, aux);
    mpz_mod(algsqrt, algsqrt, pol->n);
  
    gmp_fprintf(stderr, "Algebraic square root is: %Zd\n", algsqrt);
  
    sprintf(ratname, "%s.rat.%03d", prefix, numdep);
    ratfile = fopen(ratname, "r");
    ASSERT_ALWAYS(ratfile != NULL);
    gmp_fscanf(ratfile, "%Zd", aux);
    fclose(ratfile);
  
    gmp_fprintf(stderr, "Rational square root is: %Zd\n", aux);
    fprintf (stderr, "Total square root time is %2.2lf\n", seconds() - t0);
  
    int found = 0;
    {
      mpz_t g1, g2;
      mpz_init(g1);
      mpz_init(g2);
  
      // First check that the squares agree
      mpz_mul(g1, aux, aux);
      mpz_mod(g1, g1, pol->n);
      if(mpz_cmp_ui(pol->g[1], 1) != 0){
  	// case g(X)=m1*X+m2 with m1 != 1
  	// we should have prod (a+b*m2/m1) = A^2 = R^2/m1^(nab-nfree)
  	// and therefore nab should be even
  	if(nab & 1){
  	    fprintf(stderr, "Sorry, but #(a, b) is odd\n");
  	    fprintf(stderr, "Bug: this should be patched! Please report your buggy input\n");
  	    printf("Failed\n");
  	    return 0;
  	}
  	mpz_powm_ui(g2, pol->g[1], (nab-nfree)>>1, pol->n);
  	mpz_mul(algsqrt, algsqrt, g2);
  	mpz_mod(algsqrt, algsqrt, pol->n);
      }
      mpz_mul(g2, algsqrt, algsqrt);
      mpz_mod(g2, g2, pol->n);
      if (mpz_cmp(g1, g2)!=0) {
        fprintf(stderr, "Bug: the squares do not agree modulo n!\n");
        //      gmp_printf("g1:=%Zd;\ng2:=%Zd;\n", g1, g2);
      }
      mpz_sub(g1, aux, algsqrt);
      mpz_gcd(g1, g1, pol->n);
  
      if (mpz_cmp(g1,pol->n)) {
        if (mpz_cmp_ui(g1,1)) {
          found = 1;
          gmp_printf("%Zd\n", g1);
        }
      }
      mpz_add(g2, aux, algsqrt);
      mpz_gcd(g2, g2, pol->n);
      if (mpz_cmp(g2,pol->n)) {
        if (mpz_cmp_ui(g2,1)) {
          found = 1;
          gmp_printf("%Zd\n", g2);
        }
      }
      mpz_clear(g1);
      mpz_clear(g2);
    }
    if (!found)
      printf("Failed\n");
  
    cado_poly_clear (pol);
    mpz_clear(aux);
    mpz_clear(algsqrt);
    poly_free(F);
    poly_free(prd->p);
    poly_free(tmp->p);
  
    return 0;
}
  
