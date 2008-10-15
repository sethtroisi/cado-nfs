/* Program to tune the FFT multiplication over GF(2).
   
  Copyright 2007 Paul Zimmermann.

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

   0) First run tunetoom to tune Toom-Cook multiplication (see the
      instructions in tunetoom.c).
   1) compile this program (tunefft)
   2) run tunefft. See the usage() function below for the arguments and
      options. Be aware that for large sizes, you might need a lot of memory !
   3) compile factor.c.
*/

#define _BSD_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>		/* for LONG_MAX */
#include <float.h>		/* for DBL_MAX */
#include <time.h>
#include <sys/utsname.h>	/* for uname */
#include "gf2x.h"
#include "timing.h"
#include "tuning-common.h"

/* This version of tunefft uses the midpoint of each stair */

/* Must be at least >=MUL_FFT_THRESHOLD (28), but it saves time to set larger */
#define MUL_FFT_BEGIN_TUNE 1000

#define STEPMAX 50

/* Return the largest m >= n such that a product of m*m words with an FFT of
   length K leads to the same size of pointwise products as with n words.
   More precisely if f(n) = K/3 * ceil(6*n*WORDSIZE/K^2) then m is the
   largest integer such that
   ceil(f(n)/WORDSIZE) = ceil(f(m)/WORDSIZE).
   In fact we take for m the largest integer such that g(n) = g(m) with
   g(n) = ceil(6*n*WORDSIZE/K^2).
*/

long end_of_stair(long n, long K)
{
    long g = 6 * n * WORDSIZE;

    g = 1 + (g - 1) / K;	/* ceil(g/K) */
    g = 1 + (g - 1) / K;	/* ceil(g/K^2) */
    g = g * K * K;
    return g / (6 * WORDSIZE);
}

long next_step(long n, long K)
{
    // Return the minimum so steps are not too large
    return MIN(end_of_stair(n, K), MAX(n+1, n*mulstep));

}

void usage(int rc)
{
    FILE * f = rc ? stderr : stdout;

    fprintf(f, "Usage: tunefft [options] [min_word_size] max_word_size\n");
    fprintf(f, "Allowed options:\n");
    fprintf(f, "\t-s <multiplicative step>\n");
    fprintf(f, "\t-o <output file> (default is stdout or file desc. 3)\n");
    exit(rc);
}

int main(int argc, char *argv[])
{
    long minn, maxn, mid, n, n2, ns, i;
    long besti;			/* 0 for TC, 1, 2, ... for FFT(K0*3^(bestK-1)) */
    long bestK;
    long K, K0 = 3;		/* try K0, 3*K0, 9*K0 */
    double T[4];		/* T[0] is for TC, T[1] for K0, T[2] for 3*K0, T[3] for 9*K0 */
    double t1[4], t2[4];
    ulong *a, *b, *c, *t, *u, *v;
    int nsz = 0;
    int tc_takes_too_long = 0;
    const char * reference = "TC";

    maxn = 1000000;		// default
    minn = MUL_FFT_BEGIN_TUNE / 2 + 1;

    argc--,argv++;
    for( ; argc ; argc--,argv++) {
        int r;

        if (strcmp(argv[0], "--help") == 0) {
            usage(0);
        }

        if (strcmp(argv[0], "--no-toom") == 0) {
            tc_takes_too_long = 1;
            reference = "F1(K0)";
            continue;
        }

        r = handle_tuning_mulstep(&argc, &argv);
        if (r < 0) usage(1); else if (r) continue;
        r = handle_tuning_outfile(&argc, &argv);
        if (r < 0) usage(1); else if (r) continue;

        if (strcmp(argv[0], "-k0") == 0) {
            argc--,argv++;
            if (! argc) usage(1);
            K0 = atoi(argv[0]);
            continue;
        }

        if (nsz == 0) {
            maxn = atoi(argv[0]);
            nsz++;
            continue;
        }
        if (nsz == 1) {
            minn = maxn;
            maxn = atoi(argv[0]);
            nsz++;
            continue;
        }
        usage(1);
    }

    if (nsz == 0)
        usage(1);

    set_tuning_output();

    {
	char date[40];
	time_t t;
	size_t u;
	struct utsname buf;
	time(&t);
	ctime_r(&t, date);
	u = strlen(date);
	for (; u && isspace(date[u - 1]); date[--u] = '\0');
	uname(&buf);
        fprintf(rp, "info-fft \"%s -s %.2f %ld run on %s on %s ; based on %s\"\n",
                __FILE__,mulstep,maxn,buf.nodename,date,TOOM_TUNING_INFO);
    }

    printf("Tuning FFT multiplication to wordsize %ld\n\n", maxn);

    a = (ulong *) malloc(maxn * sizeof(ulong));
    b = (ulong *) malloc(maxn * sizeof(ulong));
    c = (ulong *) malloc(2 * maxn * sizeof(ulong));
    u = (ulong *) malloc(2 * maxn * sizeof(ulong));
    v = (ulong *) malloc(2 * maxn * sizeof(ulong));
    t = (ulong *) malloc(toomspace(maxn) * sizeof(ulong));

    random_wordstring(a, maxn);
    random_wordstring(b, maxn);

/* Skip n if (2*n < MUL_FFT_BEGIN_TUNE) as this is too small for the FFT */


    for (n = minn; n <= maxn;) {
	n2 = next_step(n, 3 * K0);	// End of interval
	if (n2 > maxn)		// Only go as far
	    n2 = maxn;		// as maxn.
	mid = (n + n2) / 2;	// Mid-point
	printf("%ld..%ld ", n, n2);
	fflush(stdout);

        if (tc_takes_too_long) {
            T[0] = DBL_MAX;
        } else {
            TIME(T[0], mul_toom(u, a, b, mid, t));	// Time Toom-Cook
            printf("TC:%1.1e ", T[0]);
        }
	fflush(stdout);
	besti = 0;
	bestK = 1;
        K = K0;
        i = 1;
ugly_label:
	for ( ; i <= 3; i++, K *= 3) {
	    TIME(t1[i], FFTMul(c, a, mid, b, mid, K));
            if (tc_takes_too_long) {
                memcpy(u, c, 2 * maxn * sizeof(ulong));
            }
            check(a, mid, b, mid, reference, u, "F1", c);
	    TIME(t2[i], FFTMul2(v, a, mid, b, mid, K));
	    check(a, mid, b, mid, "F1", c, "F2", v);
	    if (t1[i] < t2[i]) {
		T[i] = t1[i];
		printf("F1(%ld):%1.1e ", K, T[i]);
	    } else {
		T[i] = t2[i];
		printf("F2(%ld):%1.1e ", K, T[i]);
	    }
	    fflush(stdout);
	    if (T[i] < T[besti]) {
		besti = i;
		bestK = (t2[i] > t1[i]) ? K : -K;	/* -K for FFT2(|K|) */
	    }
	}

	if (T[3] < T[1] && T[3] < T[2]) {
            if (besti) {
                if (besti == 1)
                    abort();
                besti--;
            }
	    K0 *= 3;
            /* K just stays as it was */
            i = 3;
            T[1] = T[2];
            T[2] = T[3];
            goto ugly_label;
            /* Notice that we can't loop forever here. If we have T[3] <
             * T[2], this will ensure T[2] < T[1] at the next turn,
             * thereby forcing the other case not to happen */
        } else if (T[1] < T[2] && T[1] < T[3] && K0 > 3) {
	    K0 /= 3;
        }

        /* OK, this stair is done */

	if (bestK == 1)
	    printf("TC");
	else {
	    if (bestK > 0)
		printf("F1(%ld)", bestK);
	    else
		printf("F2(%ld)", -bestK);
	}
	printf("\n");
	fflush(stdout);

        if (T[0] >= 4 * T[besti] && !tc_takes_too_long) {
            printf("TC is taking too long, disabling for next sizes\n");
            tc_takes_too_long = 1;
            reference = "F1(K0)";
        }

	/* go to next size */
	ns = n;
	n = next_step(n, 3 * K0);	/* middle value of K */
	if (n > n2)
	    n = n2;		/* end of last stair if K0 increased */
	n++;
	if (n < mid) {		/* redo the last stair if K0 decreased */
	    n = ns;
        } else {
            fprintf(rp, "fft %ld %ld\n",
                    ns == minn ? 1 : ns, ns == minn ? 1 : bestK);
        }
    }

    free(a);
    free(b);
    free(c);
    free(t);
    free(u);
    free(v);

    return 0;
}
