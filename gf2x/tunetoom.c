/* Program to tune Toom-Cook multiplication over GF(2).
   
  Copyright 2007 Richard Brent and Paul Zimmermann.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02111-1307, USA.
*/

/* How to use this program:

   1) Back up the files besttoom.h and mparam.h as these will
      be overwritten by tunetoom and tunefft respectively.
   2) compile this program (tunetoom) in the NTL source directory.
      You will need to copy the mulfft-bit.c and HalfGCD.c files to this
      directory.  (And some other files?)
   3) run tunetoom, giving as argument the maximum word size, for example
      ./tunetoom 2000 will tune multiplication of polynomials up to
      degree 128000 on a 64-bit machine. For higher degrees the FFT
      is probably faster - see tunefft.
      tunetoom prints detailed logs on stdout, and writes the tuning
      information in the file besttoom.h.  There are two phases - first
      the "balanced" routines with equal-sized inputs are tuned, then the 
      "unbalanced" routines where one input is about twice as large as the 
      other are tuned.
   4) compile and run tunefft to tune FFT multiplication
      (see instructions in tunefft.c).
*/

#define _BSD_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>		/* for LONG_MAX */
#include <assert.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <sys/utsname.h>        /* for uname */
#include "gf2x.h"
#include "timing.h"
#include "replace.h"


struct hash_define replacements[100];
unsigned int nrepl=0;

/* Peek into mul-toom.c */
extern short best_tab[TOOM_TUNING_LIMIT];
extern short best_utab[TOOM_TUNING_LIMIT];

#define MINI_MUL_TOOM_THRESHOLD 	17
#define MINI_MUL_TOOMW_THRESHOLD	8
#define MINI_MUL_TOOM4_THRESHOLD	30
#define MINI_MUL_TOOMU_THRESHOLD	33

#define BESTMIN (MUL_KARA_THRESHOLD-1)
#define BESTMINU (MUL_TOOMU_THRESHOLD-1)

#define MINTIME 0.5		/* time resolution */

FILE *fp;

#define TIME(x, i)				\
  { long j, k = 1;				\
    double s0 = seconds ();			\
    do {					\
      for (j = 0; j < k; j++) (i);		\
      k = 2 * k;				\
      x = seconds () - s0;			\
    } while (x < MINTIME);			\
    (x) = (x) / (k - 1);			\
  }

void random_wordstring(unsigned long *a, long n)
{
    long i;
    for (i = 0; i < n; i++)
	a[i] = random();
}

void dump(const unsigned long *a, long m, const unsigned long *b, long n,
	  const unsigned long *c, const unsigned long *d)
{
    printf("failed:=[");
    printf("[");
    int j;
    for (j = 0; j < m; j++) {
	if (j)
	    printf(", ");
	printf("%lu", a[j]);
    }
    printf("],");
    printf("[");
    for (j = 0; j < n; j++) {
	if (j)
	    printf(", ");
	printf("%lu", b[j]);
    }
    printf("],");
    printf("[");
    for (j = 0; j < m + n; j++) {
	if (j)
	    printf(", ");
	printf("%lu", c[j]);
    }
    printf("],");
    printf("[");
    for (j = 0; j < m + n; j++) {
	if (j)
	    printf(", ");
	printf("%lu", d[j]);
    }
    printf("]");
    printf("];\n");
}
void check(const unsigned long *a, long m,
	   const unsigned long *b, long n,
	   const unsigned long *c, const unsigned long *d, long flag)
{
    if (memcmp(c, d, (m + n) * sizeof(unsigned long)) != 0) {
	if (flag == 2)
	    printf
		("Error detected: Toom3Mul and KarMul give different results\n");
	if (flag == 3)
	    printf
		("Error detected: Toom3WMul and KarMul give different results\n");
	if (flag == 4)
	    printf
		("Error detected: Toom4Mul and KarMul give different results\n");
	if ((flag < 2) || (flag > 4))
	    printf("Error in call to check, illegal flag %ld\n", flag);
	dump(a, m, b, n, c, d);
	abort();
    }
}


void tunetoom(long tablesz)
{
    long high, n, k;
    double T3[1], TK[1], TW[1], T4[1];
    double mint;
    unsigned long *a, *b, *c, *d, *t;

    high = tablesz;
    if (high < BESTMIN)
	high = BESTMIN;

    if (high > TOOM_TUNING_LIMIT) {
	printf
	    ("Increase constant TOOM_TUNING_LIMIT in thresholds.h to %ld\n",
	     high);
	exit(1);
    }

    for (n = 1; n <= TOOM_TUNING_LIMIT; n++)
	best_tab[n - 1] = 0;

    printf("Generating entries best_tab[%ld..%ld]\n", 1L, high);
    fflush(stdout);

    a = (unsigned long *) malloc(high * sizeof(unsigned long));
    b = (unsigned long *) malloc(high * sizeof(unsigned long));
    c = (unsigned long *) malloc(2 * high * sizeof(unsigned long));
    d = (unsigned long *) malloc(2 * high * sizeof(unsigned long));
    t = (unsigned long *) malloc(toomspace(high) * sizeof(unsigned long));

    for (n = BESTMIN + 1; n <= high; n++) {
	srandom(1);
	TK[0] = T3[0] = TW[0] = T4[0] = 0.0;
	printf("%ld ", n);
	fflush(stdout);
	random_wordstring(a, n);
	random_wordstring(b, n);
	if (n >= MUL_KARA_THRESHOLD)
	    TIME(TK[0], mul_kara(c, a, b, n, t));
	if (n >= MINI_MUL_TOOM_THRESHOLD) {
	    TIME(T3[0], mul_tc3(d, a, b, n, t));
	    check(a, n, b, n, c, d, 2);
	}
	if (n >= MINI_MUL_TOOMW_THRESHOLD) {
	    TIME(TW[0], mul_tc3w(d, a, b, n, t));
	    check(a, n, b, n, c, d, 3);
	}
	if (n >= MINI_MUL_TOOM4_THRESHOLD) {
	    TIME(T4[0], mul_tc4(d, a, b, n, t));
	    check(a, n, b, n, c, d, 3);
	}
	printf("TC2:%1.2e TC3:%1.2e TC3W:%1.2e TC4:%1.2e ",
	       TK[0], T3[0], TW[0], T4[0]);
	mint = TK[0];
	k = GF2X_SELECT_KARA;
	if ((T3[0] < mint) && (n >= MINI_MUL_TOOM_THRESHOLD)) {
	    mint = T3[0];
	    k = GF2X_SELECT_TC3;
	}
	if ((TW[0] < mint) && (n >= MINI_MUL_TOOMW_THRESHOLD)) {
	    mint = TW[0];
	    k = GF2X_SELECT_TC3W;
	}
	if ((T4[0] < mint) && (n >= MINI_MUL_TOOM4_THRESHOLD)) {
	    mint = T4[0];
	    k = GF2X_SELECT_TC4;
	}
	best_tab[n - 1] = (short) k;	// Final value of best_tab[n]  
	printf("best:%1.2e ", mint);
	if (k == GF2X_SELECT_KARA)
	    printf("TC2");
	if (k == GF2X_SELECT_TC3)
	    printf("TC3");
	if (k == GF2X_SELECT_TC3W)
	    printf("TC3W");
	if (k == GF2X_SELECT_TC4)
	    printf("TC4");
	printf("\n");
	fflush(stdout);
    }

    free(a);
    free(b);
    free(c);
    free(d);
    free(t);

    return;
}

/* Forms c := a*b where b has size sb, a has size sa = (sb+1)/2 (words),
   representing polynomials over GF(2).  c needs space for sa+sb words. 
   Needs space 2*sa + toomspace(sa) words in stk[0] ... 
   
   The code is essentially the same as in HalfGCD.c   */

static void mul21(unsigned long *c, const unsigned long *b, long sb,
		  const unsigned long *a, unsigned long *stk)
{
    long i, j;
    long sa = (sb + 1) / 2;
    long sc = sa + sb;
    unsigned long *v;
    v = stk;
    stk += 2 * sa;
    for (i = 0; i < sc; i++)
	c[i] = 0;
    do {
	if (sa == 0)
	    break;

	if (sa == 1) {
	    c[sb] ^= addmul_1_n(c, c, b, sb, a[0]);
	    break;
	}

	for (i = 0; i + sa <= sb; i += sa) {
	    mul_toom(v, a, b + i, sa, stk);	// Generic Toom-Cook mult.
	    for (j = 0; j < 2 * sa; j++)
		c[i + j] ^= v[j];
	}

	{
	    const unsigned long *t;
	    t = a;
	    a = b + i;
	    b = t;
	}
	{
	    long t;
	    t = sa;
	    sa = sb - i;
	    sb = t;
	}
	c = c + i;
    }
    while (1);
}

void checku(const unsigned long *a, const unsigned long *b, long n)
{
    long i;
    for (i = 0; i < n; i++) {
	if (a[i] != b[i]) {
	    printf
		("Error detected: mul_toom3u and default give different results\n");
	    printf("index %ld\n", i);
	    exit(1);
	}
    }
}

void tuneutoom(long tabsz)
{
    long high, n, k;
    double T3[1], TK[1];
    double mint;
    unsigned long *a, *b, *c, *d, *t;

    high = tabsz;
    if (high < BESTMINU)
	high = BESTMINU;

    if (high > TOOM_TUNING_LIMIT) {
	printf
	    ("Increase constant TOOM_TUNING_LIMIT in thresholds.c to %ld\n",
	     high);
	exit(1);
    }

    for (n = 1; n <= TOOM_TUNING_LIMIT; n++)
	best_utab[n - 1] = 0;

    printf("Generating entries best_utab[%ld..%ld]\n", 1L, high);
    fflush(stdout);

    long sa = high;
    long sb = (sa + 1) / 2;

    long sp1 = toomuspace(sa);	// space for mul_toom3u
    long sp2 = toomspace(sb) + 2 * sb;	// space for mul21
    long sp = (sp1 > sp2) ? sp1 : sp2;

    a = (unsigned long *) malloc(sa * sizeof(unsigned long));
    b = (unsigned long *) malloc(sb * sizeof(unsigned long));
    c = (unsigned long *) malloc(3 * sb * sizeof(unsigned long));
    d = (unsigned long *) malloc(3 * sb * sizeof(unsigned long));
    t = (unsigned long *) malloc(sp * sizeof(unsigned long));


    for (sa = BESTMINU + 1; sa <= high; sa++) {
	sb = (sa + 1) / 2;
	random_wordstring(a, sa);
	random_wordstring(b, sb);
	TK[0] = T3[0] = 0.0;
	printf("%ld ", sa);
	fflush(stdout);
	TIME(TK[0], mul21(c, a, sa, b, t));
	if (sa >= MINI_MUL_TOOMU_THRESHOLD) {
	    TIME(T3[0], mul_tc3u(d, a, sa, b, t));
	    checku(c, d, sa + sb);
	}
	printf("default:%1.2e TC3U:%1.2e ", TK[0], T3[0]);
	mint = TK[0];
	k = GF2X_SELECT_UNB_DFLT;
	if ((T3[0] < mint) && (sa >= MINI_MUL_TOOMU_THRESHOLD)) {
	    mint = T3[0];
	    k = GF2X_SELECT_UNB_TC3U;
	}
	best_utab[sa - 1] = (short) k;	// Final value
	printf("best:%1.2e ", mint);
	if (k == GF2X_SELECT_UNB_DFLT)
	    printf("default");
	if (k == GF2X_SELECT_UNB_TC3U)
	    printf("TC3U");
	printf("\n");
	fflush(stdout);
    }

    free(a);
    free(b);
    free(c);
    free(d);
    free(t);

    return;
}

void prepare_and_push_hash_define_tbl(const char * name, short * tbl, size_t n)
{
    size_t sz = 3 + 3*n+(n+19)/20*4 + 10;
    size_t h = 0;
    char * table_string = malloc(sz);
    size_t i;

    h += snprintf(table_string + h, sz-h, "{");
    for (i = 1 ; i <= n ; i++) {
        if ((i-1) % 20 == 0) {
            h += snprintf(table_string + h, sz-h, "\t\\\n\t");
        }
        h += snprintf(table_string + h, sz-h, "%d, ", tbl[i-1]);
    }
    h += snprintf(table_string + h, sz-h, "}\n");
    if (h >= sz) abort();

    set_hash_define(replacements + nrepl++, name, table_string);

    free(table_string);
}


void print_tuning_results(unsigned int blim, unsigned int ulim)
{
    /* Get the values for the balanced thresholds:
     * MUL_TOOM_THRESHOLD
     * MUL_TOOMW_THRESHOLD
     * MUL_TOOM4_THRESHOLD
     * MUL_TOOM4_ALWAYS_THRESHOLD
     */
    unsigned int t3 = UINT_MAX;
    unsigned int tw = UINT_MAX;
    unsigned int t4 = UINT_MAX;
    unsigned int t4a = UINT_MAX;

    unsigned int i;

    for (i = BESTMIN + 1; i <= blim; i++) {
	if (best_tab[i - 1] == GF2X_SELECT_TC3) {
	    t3 = i;
	    break;
	}
    }
    for (i = BESTMIN + 1; i <= blim; i++) {
	if (best_tab[i - 1] == GF2X_SELECT_TC3W) {
	    tw = i;
	    break;
	}
    }
    for (i = BESTMIN + 1; i <= blim; i++) {
	if (best_tab[i - 1] == GF2X_SELECT_TC4) {
	    t4 = i;
	    break;
	}
    }
    if (t4 < blim) {
	for (t4a = t4; t4a >= MINI_MUL_TOOM4_THRESHOLD; t4a--) {
	    if (best_tab[t4a - 1] != GF2X_SELECT_TC4) {
		t4a++;
		break;
	    }
	}
    } else {
        t4 = blim;
	t4a = blim;
    }
    /* Now do some sanity checks */

    int err = 0;

    if (!(tw <= t3)) {
	fprintf(stderr,
		"MUL_TOOMW_THRESHOLD(%d) must be below MUL_TOOM_THRESHOLD(%d)\n",
		tw, t3);
	err = 1;
    }

    set_hash_define_int(replacements + nrepl++, "MUL_TOOM_THRESHOLD", t3);
    set_hash_define_int(replacements + nrepl++, "MUL_TOOMW_THRESHOLD", tw);
    set_hash_define_int(replacements + nrepl++, "MUL_TOOM4_THRESHOLD", t4);
    set_hash_define_int(replacements + nrepl++,
            "MUL_TOOM4_ALWAYS_THRESHOLD", t4a);
    prepare_and_push_hash_define_tbl("BEST_TOOM_TABLE", best_tab, t4a);
    /* Get the values for the unbalanced thresholds:
     * MUL_TOOM3U_THRESHOLD
     * MUL_TOOM3U_ALWAYS_THRESHOLD
     */
    unsigned int tu = UINT_MAX;
    unsigned int tua = UINT_MAX;

    for (i = BESTMINU + 1; i <= ulim; i++) {
	if (best_utab[i - 1] == GF2X_SELECT_UNB_TC3U) {
	    tu = i;
	    break;
	}
    }
    if (tu < ulim) {
	for (tua = tu; tua >= MINI_MUL_TOOMU_THRESHOLD; tua--) {
	    if (best_utab[tua - 1] != GF2X_SELECT_UNB_TC3U) {
		tua++;
		break;
	    }
	}
    } else {
        tu = ulim;
	tua = ulim;
    }
    set_hash_define_int(replacements + nrepl++, "MUL_TOOMU_THRESHOLD", tu);
    set_hash_define_int(replacements + nrepl++, "MUL_TOOMU_ALWAYS_THRESHOLD", tua);
    prepare_and_push_hash_define_tbl("BEST_UTOOM_TABLE",best_utab,tua);

    if (err) {
        printf("/* Warning: Something fishy happened with this tuning */\n");
    }

    {
        char id[80];
        char date[40];
        time_t t;
        size_t u;
        struct utsname buf;
        time(&t);
        ctime_r(&t, date);
        u = strlen(date);
        for(;u && isspace(date[u-1]);date[--u]='\0');
        uname(&buf);
        snprintf(id,sizeof(id),"\"%s (%d %d) run on %s on %s\"",
                __FILE__,blim,ulim,buf.nodename,date);
        set_hash_define(replacements + nrepl++, "TOOM_TUNING_INFO", id);
    }

    replace(replacements, nrepl, "thresholds.h");
    for(i = 0 ; i < nrepl; i++ ) {
        free(replacements[i].identifier);
        free(replacements[i].string);
    }
}

int main(int argc, char *argv[])
{
    long tabsz1, tabsz2;
    if ((argc != 2) && (argc != 3)) {
	fprintf(stderr, "Usage: tunetoom table-size1 [table-size2]\n");
	fprintf(stderr, " where %d <= table-size1 <= %d\n",
		BESTMIN, TOOM_TUNING_LIMIT);
	fprintf(stderr, " and  %d <= table-size2 <= %d\n",
		BESTMINU, TOOM_TUNING_LIMIT);
	exit(1);
    }


    tabsz1 = tabsz2 = atoi(argv[1]);
    if (argc == 3)
	tabsz2 = atoi(argv[2]);

    tunetoom(tabsz1);		// Tune balanced routines

    tuneutoom(tabsz2);		// Tune unbalanced routines

    print_tuning_results(tabsz1, tabsz2);

    fflush(stdout);

    return 0;
}

/* vim: set sw=4 sta et: */
