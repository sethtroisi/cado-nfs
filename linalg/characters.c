/* Characters

Copyright 2008 Pierrick Gaudry, Fran\c{c}ois Morain, Emmanuel Thom\'e, Paul Zimmermann

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <ctype.h>
#include "mod_ul.c"

#include "cado.h"
#include "utils/utils.h"

#include "files.h"
#include "gzip.h"

#define DEBUG 0

int kernel(mp_limb_t * mat, mp_limb_t ** ker, int nrows, int ncols,
           int limbs_per_row, int limbs_per_col);

// Data structure for one algebraic prime and for a table of those.
typedef struct {
    unsigned long prime;
    unsigned long root;
} rootprime_t;

static int eval_char(long a, unsigned long b, rootprime_t ch)
{
    unsigned long ua, aux, saveb = b;
    int res;
#if DEBUG >= 1
    printf("a := %ld; b := %lu; p := %lu; r := %lu;", a, b, ch.prime,
           ch.root);
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
        fprintf(stderr, "  a = %ld, b = %lu, p := %lu; r := %lu;", a, saveb,
                ch.prime, ch.root);
    }
#if DEBUG >= 1
    printf("res := %d; assert (JacobiSymbol(a-b*r, p) eq res);\n", res);
#endif
    return res;
}


// returns the sign of m1*a+m2*b
static int eval_rat_char(long a, unsigned long b, cado_poly pol)
{
    mpz_t tmp1, tmp2;
    int s;

    /* first perform a quick check */
    s = (a > 0) ? mpz_sgn(pol->g[1]) : -mpz_sgn(pol->g[1]);
    if (mpz_sgn(pol->g[0]) == s)
        return s;

    mpz_init(tmp1);
    mpz_mul_si(tmp1, pol->g[1], a);
    mpz_init(tmp2);
    mpz_mul_ui(tmp2, pol->g[0], b);
    mpz_add(tmp1, tmp1, tmp2);
    s = (mpz_sgn(tmp1) >= 0 ? 1 : -1);
    mpz_clear(tmp1);
    mpz_clear(tmp2);

    return s;
}

static void create_characters(rootprime_t * tabchar, int k, cado_poly pol)
{
    unsigned long p;
    int ret;
    int i = 0;
    mpz_t pp;
    unsigned long *roots;

    /* we want some prime beyond the (algebraic) large prime bound */
    mpz_init_set_ui (pp, 1UL << pol->lpba);
    roots = malloc(pol->degree * sizeof(unsigned long));

    do {
        mpz_nextprime(pp, pp);
        p = mpz_get_ui(pp);
        ret = poly_roots_ulong(roots, pol->f, pol->degree, p);
        if (ret == 0)
            continue;
        tabchar[i].prime = p;
        /* FIXME: Why do we take only one root per prime ??? */
        tabchar[i].root = roots[0];
        i++;
    }
    while (i < k);

    free(roots);
    mpz_clear(pp);
}

/* Read one dependency from file 'file', and stores it into 'vec'.
   The dependency in 'file' is stored on one line, with 64-bit words
   written in hexadecimal.
 */
static void
readOneKer (mp_limb_t * vec, FILE * file, int nlimbs)
{
    uint64_t w;
    int ret, i, j;

    /* the code below assumes the number of bits per limb divides 64 */
    ASSERT_ALWAYS ((64 % GMP_NUMB_BITS) == 0);

    for (i = 0; i < nlimbs;)
      {
        ret = fscanf (file, "%" SCNx64, &w);
	//        ret = fscanf (file, "%lx", &w);
        ASSERT(ret == 1);
	for (j = 0; j < 64; j += GMP_NUMB_BITS)
	  {
	    *vec = w; /* automatic mask of low GMP_NUMB_BITS bits */
	    vec ++;
	    i ++;
	    w >>= GMP_NUMB_BITS;
	  }
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
computeAllCharacters(char **charbig, int i, int k, rootprime_t * tabchar,
                     long a, unsigned long b, cado_poly pol)
{
    int j;

    for (j = 0; j < k - 2; j++)
        charbig[i][j] = (char) eval_char(a, b, tabchar[j]);
    charbig[i][k - 2] = (char) eval_rat_char(a, b, pol);
    // last column is for the free relations...
    charbig[i][k - 1] = (b == 0 ? (char) (-1) : 1);        // ohhhhhhhhhhh!
}

// charmat is small_nrows x k
// relfile is the big rough relation files;
// purgedfile contains crunched rows, with their label referring to relfile;
// indexfile contains the coding "row[i] uses rows i_0...i_r in purgedfile".
static void
buildCharacterMatrix(char **charmat, int k, rootprime_t * tabchar,
                     FILE * purgedfile, FILE * indexfile, FILE * relfile,
                     cado_poly pol, int small_nrows)
{
    relation_t rel;
    int i, j, r, nr, nrows, ncols, irel, kk;
    char str[1024];
    char **charbig;

    // let's dump purgedfile which is a nrows x ncols matrix
    //    rewind(purgedfile);
    fgets(str, 1024, purgedfile);
    sscanf(str, "%d %d", &nrows, &ncols);
    fprintf(stderr, "Reading indices in purgedfile\n");

    fprintf(stderr, "Reading all (a, b)'s just once\n");
    // charbig is nrows x k
    charbig = (char **) malloc(nrows * sizeof(char *));
    for (i = 0; i < nrows; i++)
        charbig[i] = (char *) malloc(k * sizeof(char));
    irel = 0;
    //    rewind(relfile); // useless?
    for (i = 0; i < nrows; i++) {
        fgets(str, 1024, purgedfile);
        sscanf(str, "%d", &nr);
        jumpToRelation(&rel, relfile, irel, nr);
        irel = nr + 1;
        computeAllCharacters(charbig, i, k, tabchar, rel.a, rel.b, pol);
        clear_relation(&rel);
        if ((i + 1) % 100000 == 0) {
            fprintf (stderr, "   read %d/%d (a,b) pairs\r", i + 1, nrows);
        }
    }
    fprintf(stderr, "%d rows done     \n",nrows);

    fprintf(stderr, "Reading index file to reconstruct the characters\n");
    // read all relation-sets and update charmat accordingly
    for (i = 0; i < small_nrows; i++) {
        if (!(i % 10000))
            fprintf(stderr, "Treating relation #%d / %d at %2.2lf\n",
                    i, small_nrows, seconds());
        fscanf(indexfile, "%d", &nr);
        for (j = 0; j < nr; j++) {
            fscanf(indexfile, "%d", &r);
            for (kk = 0; kk < k; kk++)
                charmat[i][kk] *= charbig[r][kk];
        }
    }
    for (i = 0; i < nrows; i++)
        free(charbig[i]);
    free(charbig);
}

#if DEBUG >= 1
static void printCompactMatrix(mp_limb_t ** A, int nrows, int ncols)
{
    int i, j;

    fprintf(stderr, "array([\n");
    for (i = 0; i < nrows; i++) {
        fprintf(stderr, "[");
        for (j = 0; j < ncols; j++) {
            int j0 = j / GMP_NUMB_BITS;
            int j1 = j - j0 * GMP_NUMB_BITS;
            if ((A[i][j0] >> j1) & 1UL)
                fprintf(stderr, "1");
            else
                fprintf(stderr, "0");
            if (j < ncols - 1)
                fprintf(stderr, ", ");
        }
        fprintf(stderr, "]");
        if (i < nrows - 1)
            fprintf(stderr, ", ");
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "])");
}

static void printTabMatrix(dense_mat_t * mat, int nrows, int ncols)
{
    int i, j;

    fprintf(stderr, "array([\n");
    for (i = 0; i < nrows; i++) {
        fprintf(stderr, "[");
        for (j = 0; j < ncols; j++) {
            int j0 = j / GMP_NUMB_BITS;
            int j1 = j - j0 * GMP_NUMB_BITS;
            if ((mat->data[i * mat->limbs_per_row + j0] >> j1) & 1UL)
                fprintf(stderr, "1");
            else
                fprintf(stderr, "0");
            if (j < ncols - 1)
                fprintf(stderr, ", ");
        }
        fprintf(stderr, "]");
        if (i < nrows - 1)
            fprintf(stderr, ", ");
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "])");
}
#endif

// We add the heaviest columns of M_small to the last [kmin..kmin+skip[
// columns of charmat.
static void
addHeavyBlock(char **charmat, FILE * smallfile, int kmin, int skip)
{
    int i, nc, j, u;

    fprintf(stderr, "Adding heavy block of width %d\n", skip);
    fscanf(smallfile, "%d %d", &i, &nc);
    i = 0;
    while (fscanf(smallfile, "%d", &nc) != EOF) {
        for (u = 0; u < nc; u++) {
            fscanf(smallfile, "%d", &j);
            if (j < skip)
                charmat[i][kmin + j] = (char) (-1);        // humf...!
        }
        i++;
    }
    fprintf (stderr, "   done at %2.2lf\n", seconds ());
}

static void
handleKer(dense_mat_t * mat, rootprime_t * tabchar, FILE * purgedfile,
          mp_limb_t ** ker, cado_poly pol,
          FILE * indexfile, FILE * relfile,
          int small_nrows, FILE * smallfile, int skip)
{
    int i, n;
    unsigned int j, k;
    char **charmat;

    n = mat->nrows;
    k = mat->ncols;                // [(nchar-1)+1]+skip+1
    charmat = (char **) malloc(small_nrows * sizeof(char *));
    ASSERT(charmat != NULL);
    for (i = 0; i < small_nrows; ++i) {
        charmat[i] = (char *) malloc(k * sizeof(char));
        ASSERT(charmat[i] != NULL);
        for (j = 0; j < k; ++j)
            charmat[i][j] = 1;
    }
    buildCharacterMatrix(charmat, k - skip, tabchar,
                         purgedfile, indexfile, relfile, pol, small_nrows);
    if (skip > 0)
        addHeavyBlock(charmat, smallfile, k - skip, skip);
#ifndef NDEBUG
    fprintf (stderr, "Checking character matrix\n");
    for (i = 0; i < small_nrows; ++i)
        for (j = 0; j < k; j++)
            ASSERT((charmat[i][j] == 1) || (charmat[i][j] == -1));
    fprintf (stderr, "   done at %2.2lf\n", seconds ());
#endif
#if DEBUG >= 1
    fprintf(stderr, "charmat:=array([");
    for (i = 0; i < small_nrows; ++i) {
        fprintf(stderr, "[");
        for (j = 0; j < k; j++) {
            fprintf(stderr, "%d", (charmat[i][j] == 1 ? 0 : 1));
            if (j < k - 1)
                fprintf(stderr, ", ");
        }
        fprintf(stderr, "]");
        if (i < small_nrows - 1)
            fprintf(stderr, ", ");
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "]);\n");
    fprintf(stderr, "ker:=");
    printCompactMatrix(ker, n, small_nrows);
    fprintf(stderr, ";\n");
#endif

    // now multiply: ker * charmat
    fprintf (stderr, "Multiply ker and character matrix\n");
    for (i = 0; i < n; ++i)
        for (j = 0; j < mat->limbs_per_row; ++j)
            mat->data[i * mat->limbs_per_row + j] = 0UL;

    // mat[i, j] = sum ker[i][u] * charmat[u][j]
    {
      int u;
      for (u = 0; u < small_nrows; u++)
        {
          int u0 = u / GMP_NUMB_BITS;
          int u1 = u - u0 * GMP_NUMB_BITS;
          for (i = 0; i < n; ++i)
            if ((ker[i][u0] >> u1) & 1UL)
              for (j = 0; j < k; ++j)
                {
                  /* GCC is smart enough to do the following with masks and
                     shifts when GMP_NUMB_BITS is a power of two */
                  int j0 = j / GMP_NUMB_BITS;
                  int j1 = j - j0 * GMP_NUMB_BITS;
                  mat->data[i * mat->limbs_per_row + j0] ^=
                    (charmat[u][j] == -1) << j1;
                }
          }
      }
#if DEBUG >= 1
    fprintf(stderr, "prod:=");
    printTabMatrix(mat, n, k);
    fprintf(stderr, ";\n");
#endif

    for (i = 0; i < small_nrows; ++i)
        free(charmat[i]);
    free(charmat);
}

// matrix M_purged is nrows x ncols
// matrix M_small is small_nrows x small_ncols, the kernel of which
// is contained in ker and is n x small_nrows; small_[nrows,ncols] are
// accessible in file indexfile.
// charmat is small_nrows x k; ker * charmat will be n x k, yielding M_tiny.
int main(int argc, char **argv)
{
    FILE *purgedfile = NULL;
    FILE *kerfile = NULL;
    FILE *indexfile = NULL;
    FILE *relfile = NULL;
    FILE *smallfile = NULL;
    int ret;
    int k, isz, skip = 0;
    unsigned int i, j, n, nlimbs;
    rootprime_t *tabchar;
    cado_poly pol;
    mp_limb_t **ker;
    mp_limb_t *newker;
    dense_mat_t mymat;
    mp_limb_t **myker;
    unsigned int dim;

#if 0
    // FIXME...
    if (argc != 8) {
        fprintf(stderr, "usage: %s purgedfile kerfile polyfile", argv[0]);
        fprintf(stderr, "indexfile relfile n k\n");
        fprintf(stderr,
                "  where n is the number of kernel vector to deal with\n");
        fprintf(stderr,
                "    and k is the number of characters you want to use\n");
        exit(1);
    }
#endif

    /* print the command line */
    fprintf(stderr, "%s.r%s", argv[0], REV);
    for (k = 1; k < argc; k++)
        fprintf(stderr, " %s", argv[k]);
    fprintf(stderr, "\n");

    n = UINT_MAX;

    cado_poly_init(pol);
    while (argc > 1 && argv[1][0] == '-') {
        if (argc > 2 && strcmp(argv[1], "-purged") == 0) {
            purgedfile = gzip_open(argv[2], "r");
            argc -= 2;
            argv += 2;
        }
        if (argc > 2 && strcmp(argv[1], "-ker") == 0) {
            kerfile = fopen(argv[2], "r");
            argc -= 2;
            argv += 2;
        }
        if (argc > 2 && strcmp(argv[1], "-poly") == 0) {
            cado_poly_read(pol, argv[2]);
            argc -= 2;
            argv += 2;
        }
        if (argc > 2 && strcmp(argv[1], "-index") == 0) {
            indexfile = fopen(argv[2], "r");
            argc -= 2;
            argv += 2;
        }
        if (argc > 2 && strcmp(argv[1], "-rel") == 0) {
            relfile = gzip_open (argv[2], "r");
            argc -= 2;
            argv += 2;
        }
        if (argc > 2 && strcmp(argv[1], "-nker") == 0) {
            n = atoi(argv[2]);
            argc -= 2;
            argv += 2;
        }
        if (argc > 2 && strcmp(argv[1], "-nchar") == 0) {
            k = atoi(argv[2]);
            argc -= 2;
            argv += 2;
        }
        if (argc > 2 && strcmp(argv[1], "-small") == 0) {
            smallfile = fopen(argv[2], "r");
            argc -= 2;
            argv += 2;
        }
        if (argc > 2 && strcmp(argv[1], "-skip") == 0) {
            skip = atoi(argv[2]);
            argc -= 2;
            argv += 2;
        }
    }

    ASSERT_ALWAYS(n != UINT_MAX);
    ASSERT_ALWAYS(purgedfile != NULL);
    ASSERT_ALWAYS(kerfile != NULL);
    ASSERT_ALWAYS(indexfile != NULL);
    ASSERT_ALWAYS(relfile != NULL);
    ASSERT_ALWAYS(smallfile != NULL);

    // only k-1, since the k-th character is sign on rational side
    tabchar = (rootprime_t *) malloc((k - 1) * sizeof(rootprime_t));
    ASSERT(tabchar != NULL);

    create_characters(tabchar, k - 1, pol);

#if DEBUG >= 2
    fprintf(stderr, "using characters (p,r):\n");
    for (i = 0; i < k - 1; ++i)
        fprintf(stderr, "\t%lu %lu\n", tabchar[i].prime, tabchar[i].root);
#endif

    int small_nrows, small_ncols;
    {
        ret = fscanf(indexfile, "%d %d", &small_nrows, &small_ncols);
        ASSERT(ret == 2);
        nlimbs = ((small_nrows - 1) / GMP_NUMB_BITS) + 1;
    }

    ker = (mp_limb_t **) malloc(n * sizeof(mp_limb_t *));
    ASSERT(ker != NULL);
    fprintf(stderr, "Reading dependency");
    for (j = 0; j < n; ++j) {
        fprintf(stderr, " %u", j);
        fflush (stderr);
        ker[j] = (mp_limb_t *) malloc(nlimbs * sizeof(mp_limb_t));
        ASSERT(ker[j] != NULL);
        readOneKer(ker[j], kerfile, nlimbs);
    }
    fprintf(stderr, "\nfinished reading kernel file\n");

    mymat.nrows = n;
    mymat.ncols = k + skip + 1;        // 1 (rat sign)+(k-1) characters+skip+1 (freerels)
    mymat.limbs_per_row = ((mymat.ncols - 1) / GMP_NUMB_BITS) + 1;
    mymat.limbs_per_col = ((mymat.nrows - 1) / GMP_NUMB_BITS) + 1;

    mymat.data =
        (mp_limb_t *) malloc(mymat.limbs_per_row * mymat.nrows *
                             sizeof(mp_limb_t));
    ASSERT(mymat.data != NULL);

    fprintf(stderr, "start computing characters...\n");

    handleKer(&mymat, tabchar, purgedfile, ker, pol,
              indexfile, relfile, small_nrows, smallfile, skip);

    /* this is a cheap loop, since mymat.nrows is small, typically 64 */
    myker = (mp_limb_t **) malloc(mymat.nrows * sizeof(mp_limb_t *));
    ASSERT(myker != NULL);
    for (i = 0; i < mymat.nrows; ++i) {
        myker[i] =
            (mp_limb_t *) malloc(mymat.limbs_per_col * sizeof(mp_limb_t));
        ASSERT(myker[i] != NULL);
        for (j = 0; j < mymat.limbs_per_col; ++j)
            myker[i][j] = 0UL;
    }
    fprintf(stderr, "%d rows done      \n", mymat.nrows);
    fprintf(stderr, "Computing tiny kernel\n");
    dim =
        kernel(mymat.data, myker, mymat.nrows, mymat.ncols,
               mymat.limbs_per_row, mymat.limbs_per_col);
    fprintf(stderr, "dim of ker = %d\n", dim);

#if DEBUG >= 1
    fprintf(stderr, "newker:=");
    printCompactMatrix(myker, dim, k);
    fprintf(stderr, ";\n");
#endif

    newker = (mp_limb_t *) malloc(nlimbs * sizeof(mp_limb_t));
    ASSERT(newker != NULL);
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
                        newker[kk] ^= ker[j * GMP_NUMB_BITS + jj][kk];
                }
                w >>= 1;
            }
        }
        // do not print zero vector...!
        isz = 1;
        for (j = 0; j < nlimbs; ++j)
            if (newker[j]) {
                isz = 0;
                break;
            }
        if (isz)
            fprintf(stderr, "Sorry %d-th vector is 0\n", i);
        else {
            for (j = 0; j < nlimbs; ++j)
                printf("%lx ", newker[j]);
            printf("\n");
        }
    }

    free(tabchar);
    for (j = 0; j < n; ++j)
        free(ker[j]);
    free(ker);
    free(newker);
    free(mymat.data);
    for (i = 0; i < mymat.nrows; ++i)
        free(myker[i]);
    free(myker);

    return 0;
}
