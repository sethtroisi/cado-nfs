#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "mod_ul.c"

#define WANT_ASSERT

#include "cado.h"
#include "utils/utils.h"

#include "files.h"

#define DEBUG 0

int kernel(mp_limb_t* mat, mp_limb_t** ker, int nrows, int ncols,
    int limbs_per_row, int limbs_per_col);

// Data structure for one algebraic prime and for a table of those.
typedef struct {
  unsigned long prime;
  unsigned long root;
} rootprime_t;

static int eval_char(long a, unsigned long b, rootprime_t ch) {
  unsigned long ua, aux, saveb = b;
  int res;
#if DEBUG >= 1
  printf("a := %ld; b := %lu; p := %lu; r := %lu;", a, b, ch.prime, ch.root);
#endif
  if (a < 0) {
    ua = ((unsigned long) (-a)) % ch.prime;
    b = b % ch.prime;
    modul_mul(&aux, &b, &ch.root, &ch.prime);
    modul_add(&aux, &ua, &aux, &ch.prime);
    modul_neg(&aux, &aux, &ch.prime);
    res = modul_jacobi(&aux, &ch.prime);
  } else {
    ua = ((unsigned long) (a)) % ch.prime;
    b = b % ch.prime;
    modul_mul(&aux, &b, &ch.root, &ch.prime);
    modul_sub(&aux, &ua, &aux, &ch.prime);
    res = modul_jacobi(&aux, &ch.prime);
  }
  if (res == 0) {
    fprintf(stderr, "Strange: have a jacobi symbol that is zero:\n");
    fprintf(stderr, "  a = %ld, b = %lu, p := %lu; r := %lu;", a, saveb, ch.prime, ch.root);
   }
#if DEBUG >= 1
  printf("res := %d; assert (JacobiSymbol(a-b*r, p) eq res);\n", res);
#endif
  return res;
}

static void
create_characters (rootprime_t * tabchar, int k, cado_poly pol)
{
  unsigned long p;
  int ret;
  int i = 0;
  mpz_t pp;
  LONG *roots;

  mpz_init_set_ui (pp, 1 << (pol->lpba + 1));
  roots = malloc (pol->degree * sizeof (LONG));
  
  do
    {
      mpz_nextprime (pp, pp);
      p = mpz_get_ui (pp);
      ret = roots_mod_long (roots, pol->f, pol->degree, p);
      if (ret == 0)
        continue;
      tabchar[i].prime = p;
      if (roots[0] > 0 )
	tabchar[i].root = roots[0];
      else 
	tabchar[i].root = p - (-roots[0]);
      i++;
    }
  while (i < k);

  free (roots);
  mpz_clear(pp);
}

static void
readOneKer(mp_limb_t *vec, FILE *file, int nlimbs) {
  unsigned long w;
  int ret, i;
  for (i = 0; i < nlimbs; ++i) {
    ret = fscanf(file, "%lx", &w);
    ASSERT (ret == 1);
    *vec = w;
    vec++;
  }
}

typedef struct {
  unsigned int nrows;
  unsigned int ncols;
  mp_limb_t *data;
  unsigned int limbs_per_row;
  unsigned int limbs_per_col;
} dense_mat_t;

// charmat is small_nrows x k
static void
computeAllCharacters(char **charmat, int i, int k, rootprime_t * tabchar, long a, unsigned long b)
{
    int j;
    
    for(j = 0; j < k; j++)
	charmat[i][j] = (char)eval_char(a, b, tabchar[j]);
}

// charmat is small_nrows x k
// relfile is the big rough relation files;
// purgedfile contains crunched rows, with their label referring to relfile;
// indexfile contains the coding "row[i] uses rows i_0...i_r in purgedfile".
static void
buildCharacterMatrix(char **charmat, int k, rootprime_t * tabchar, FILE *purgedfile, FILE *indexfile, FILE *relfile)
{
    relation_t rel;
    int i, j, r, small_nrows, nr, nrows, ncols, irel, kk;
    char str[1024];
    char **charbig;

    // let's dump purgedfile which is a nrows x ncols matrix
    rewind(purgedfile);
    fgets(str, 1024, purgedfile);
    sscanf(str, "%d %d", &nrows, &ncols);
    fprintf(stderr, "Reading indices in purgedfile\n");

    fprintf(stderr, "Reading all (a, b)'s just once\n");
    // charbig is nrows x k
    charbig = (char **)malloc(nrows * sizeof(char *));
    for(i = 0; i < nrows; i++)
	charbig[i] = (char *)malloc(k * sizeof(char));
    irel = 0;
    //    rewind(relfile); // useless?
    for(i = 0; i < nrows; i++){
	fgets(str, 1024, purgedfile);
	sscanf(str, "%d", &nr);
	jumpToRelation(&rel, relfile, irel, nr);
	irel = nr+1;
	computeAllCharacters(charbig, i, k, tabchar, rel.a, rel.b);
	clear_relation(&rel);
    }
    fprintf(stderr, "Reading index file to reconstruct the characters\n");
    rewind(indexfile);
    // skip small_nrows, small_ncols
    fscanf(indexfile, "%d %d", &small_nrows, &r);
    // read all relation-sets and update charmat accordingly
    for(i = 0; i < small_nrows; i++){
	if(!(i % 10000))
	    fprintf(stderr, "Treating relation #%d / %d at %2.2lf\n",
		    i, small_nrows, seconds());
	fscanf(indexfile, "%d", &nr);
	for(j = 0; j < nr; j++){
	    fscanf(indexfile, "%d", &r);
	    for(kk = 0; kk < k; kk++)
		charmat[i][kk] *= charbig[r][kk];
	}
    }
    for(i = 0; i < nrows; i++)
	free(charbig[i]);
    free(charbig);
}

#if DEBUG >= 1
static void
printCompactMatrix(mp_limb_t **A, int nrows, int ncols)
{
    int i, j;

    fprintf(stderr, "array([\n");
    for(i = 0; i < nrows; i++){
        fprintf(stderr, "[");
	for(j = 0; j < ncols; j++){
	    int j0 = j/GMP_NUMB_BITS;
            int j1 = j - j0*GMP_NUMB_BITS;
	    if((A[i][j0]>>j1) & 1UL)
		fprintf(stderr, "1");
	    else
		fprintf(stderr, "0");
	    if(j < ncols-1)
		fprintf(stderr, ", ");
	}
	fprintf(stderr, "]");
	if(i < nrows-1)
	    fprintf(stderr, ", ");
	fprintf(stderr, "\n");
    }
    fprintf(stderr, "])");
}

static void
printTabMatrix(dense_mat_t *mat, int nrows, int ncols)
{
    int i, j;

    fprintf(stderr, "array([\n");
    for(i = 0; i < nrows; i++){
        fprintf(stderr, "[");
	for(j = 0; j < ncols; j++){
	    int j0 = j/GMP_NUMB_BITS;
            int j1 = j - j0*GMP_NUMB_BITS;
	    if((mat->data[i*mat->limbs_per_row+j0]>>j1) & 1UL)
		fprintf(stderr, "1");
	    else
		fprintf(stderr, "0");
	    if(j < ncols-1)
		fprintf(stderr, ", ");
	}
	fprintf(stderr, "]");
	if(i < nrows-1)
	    fprintf(stderr, ", ");
	fprintf(stderr, "\n");
    }
    fprintf(stderr, "])");
}
#endif

static void
handleKer(dense_mat_t *mat, rootprime_t * tabchar, FILE * purgedfile,
          mp_limb_t ** ker, /*int nlimbs, cado_poly pol, */
          FILE *indexfile, FILE *relfile,
	  int small_nrows)
{  
    int i, n;
    unsigned int j, k;
    char **charmat;

    n = mat->nrows;
    k = mat->ncols;
    charmat = (char **) malloc(small_nrows * sizeof(char *));
    ASSERT (charmat != NULL);
    for (i = 0; i < small_nrows; ++i) {
	charmat[i] = (char *) malloc(k*sizeof(char));
	ASSERT (charmat[i] != NULL);
	for (j = 0; j < k; ++j)
	    charmat[i][j] = 1;
    }
    buildCharacterMatrix(charmat, k, tabchar, purgedfile, indexfile, relfile);
#ifdef WANT_ASSERT
    for (i = 0; i < small_nrows; ++i)
	for(j = 0; j < k; j++)
	    ASSERT((charmat[i][j] == 1) || (charmat[i][j] == -1));
#endif
#if DEBUG >= 1
    fprintf(stderr, "charmat:=array([");
    for (i = 0; i < small_nrows; ++i){
	fprintf(stderr, "[");
	for(j = 0; j < k; j++){
	    fprintf(stderr, "%d", (charmat[i][j] == 1 ? 0 : 1));
	    if(j < k-1)
		fprintf(stderr, ", ");
	}
	fprintf(stderr, "]");
	if(i < small_nrows-1)
	    fprintf(stderr, ", ");
	fprintf(stderr, "\n");
    }
    fprintf(stderr, "]);\n");
    fprintf(stderr, "ker:=");
    printCompactMatrix(ker, n, small_nrows);
    fprintf(stderr, ";\n");
#endif

    // now multiply: ker * charmat 
    for (i = 0; i < n; ++i)
	for (j = 0; j < mat->limbs_per_row; ++j)
	    mat->data[i*mat->limbs_per_row+j] = 0UL;

    // mat[i, j] = sum ker[i][u] * charmat[u][j]
    for(i = 0; i < n; ++i){
	for(j = 0; j < k; ++j){
	    int j0 = j/GMP_NUMB_BITS;
	    int j1 = j - j0*GMP_NUMB_BITS;
	    int u;
	    for(u = 0; u < small_nrows; u++){
		int u0 = u/GMP_NUMB_BITS;
		int u1 = u - u0*GMP_NUMB_BITS;
		if((ker[i][u0]>>u1) & 1UL)
		    if(charmat[u][j] == -1)
			mat->data[i*mat->limbs_per_row+j0] ^= (1UL<<j1);
	    }
	}
    }
#if DEBUG >= 1
    fprintf(stderr, "prod:=");
    printTabMatrix(mat, n, k);
    fprintf(stderr, ";\n");
#endif
    // TODO: clean data
}

// matrix M_purged is nrows x ncols
// matrix M_small is small_nrows x small_ncols, the kernel of which
// is contained in ker and is n x small_nrows; small_[nrows,ncols] are
// accessible in file indexfile.
// charmat is small_nrows x k; ker * charmat will be n x k, yielding M_tiny.
int main(int argc, char **argv) {
  FILE *purgedfile, *kerfile, *indexfile, *relfile;
  int ret;
  int k, isz;
  unsigned int i, j, n, nlimbs;
  rootprime_t *tabchar;
  int *charval;
  cado_poly pol;
  mp_limb_t **ker;
  mp_limb_t *newker;
  dense_mat_t mymat;
  mp_limb_t **myker;
  unsigned int dim;

  if (argc != 8) {
    fprintf(stderr, "usage: %s purgedfile kerfile polyfile", argv[0]);
    fprintf(stderr, "indexfile relfile n k\n");
    fprintf(stderr, "  where n is the number of kernel vector to deal with\n");
    fprintf(stderr, "    and k is the number of characters you want to use\n");
    exit(1);
  }

  /* print the command line */
  fprintf (stderr, "%s.r%s", argv[0], REV);
  for (k = 1; k < argc; k++)
    fprintf (stderr, " %s", argv[k]);
  fprintf (stderr, "\n");

  purgedfile = fopen(argv[1], "r");
  ASSERT (purgedfile != NULL);
  kerfile = fopen(argv[2], "r");
  ASSERT (kerfile != NULL);

  read_polynomial(pol, argv[3]);

  indexfile = fopen(argv[4], "r");
  ASSERT (indexfile != NULL);
  relfile = fopen(argv[5], "r");
  ASSERT (relfile != NULL);

  ret = gmp_sscanf(argv[6], "%d", &n);
  ASSERT (ret == 1);
  ret = gmp_sscanf(argv[7], "%d", &k);
  ASSERT (ret == 1);

  tabchar = (rootprime_t *)malloc(k*sizeof(rootprime_t));
  ASSERT (tabchar != NULL);
  charval = (int *)malloc(k*sizeof(int));
  ASSERT (charval != NULL);

  create_characters(tabchar, k, pol);

#if DEBUG >= 2
  fprintf(stderr, "using characters (p,r):\n");
  for (i = 0; i < k; ++i)
      fprintf(stderr, "\t%lu %lu\n", tabchar[i].prime, tabchar[i].root);
#endif

  int small_nrows, small_ncols;
  {
    ret = fscanf(indexfile, "%d %d", &small_nrows, &small_ncols);
    ASSERT (ret == 2);
    nlimbs = (small_nrows / GMP_NUMB_BITS) + 1;
  }

  ker = (mp_limb_t **)malloc(n*sizeof(mp_limb_t *));
  ASSERT (ker != NULL);
  for (j = 0; j < n; ++j) {
    fprintf (stderr, "Reading dependency %u\n", j);
    ker[j] = (mp_limb_t *) malloc (nlimbs*sizeof(mp_limb_t));
    ASSERT (ker[j] != NULL);
    readOneKer(ker[j], kerfile, nlimbs);
  }
  fprintf(stderr, "finished reading kernel file\n");

  mymat.nrows = n;
  mymat.ncols = k;
  mymat.limbs_per_row =  ((mymat.ncols-1) / GMP_NUMB_BITS) + 1;
  mymat.limbs_per_col =  ((mymat.nrows-1) / GMP_NUMB_BITS) + 1;

  mymat.data = (mp_limb_t *) malloc(mymat.limbs_per_row*mymat.nrows*sizeof(mp_limb_t));
  ASSERT (mymat.data != NULL);

  fprintf(stderr, "start computing characters...\n");

  handleKer(&mymat, tabchar, purgedfile, ker, /*nlimbs, *pol, */
            indexfile, relfile, small_nrows);

  myker = (mp_limb_t **)malloc(mymat.nrows*sizeof(mp_limb_t *));
  ASSERT (myker != NULL);
  for (i = 0; i < mymat.nrows; ++i) {
    myker[i] = (mp_limb_t *)malloc(mymat.limbs_per_col*sizeof(mp_limb_t));
    ASSERT (myker[i] != NULL);
    for (j = 0; j < mymat.limbs_per_col; ++j) 
      myker[i][j] = 0UL;
  }
  fprintf(stderr, "Computing tiny kernel\n");
  dim = kernel(mymat.data, myker, mymat.nrows, mymat.ncols, mymat.limbs_per_row, mymat.limbs_per_col);
  fprintf(stderr, "dim of ker = %d\n", dim);

#if DEBUG >= 1
  fprintf(stderr, "newker:=");
  printCompactMatrix(myker, dim, k);
  fprintf(stderr, ";\n");
#endif
    
  newker = (mp_limb_t *) malloc (nlimbs*sizeof(mp_limb_t));
  ASSERT (newker != NULL);
  for (i = 0; i < dim; ++i) {
    for (j = 0; j < nlimbs; ++j)
      newker[j] = 0;
    for (j = 0; j < mymat.limbs_per_col; ++j) {
      int jj;
      unsigned int kk;
      unsigned long w = myker[i][j];
      for (jj = 0; jj < GMP_NUMB_BITS; ++jj) {
	if (w & 1UL) {
	  for (kk = 0; kk < nlimbs; ++kk)
	    newker[kk] ^= ker[j*GMP_NUMB_BITS+jj][kk];
	}
	w >>= 1;
      }
    }
    // do not print zero vector...!
    isz = 1;
    for(j = 0; j < nlimbs; ++j)
	if(newker[j]){
	    isz = 0;
	    break;
	}
    if(isz)
	fprintf(stderr, "Sorry %d-th vector is 0\n", i);
    else{
	for(j = 0; j < nlimbs; ++j)
	    printf("%lx ", newker[j]);
	printf("\n");
    }
  }

  return 0;
}
