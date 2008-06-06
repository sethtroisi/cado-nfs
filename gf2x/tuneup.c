/* Tuning program for GF(2)[x] schoolbook/Karatsuba/Toom-Cook multiplication.

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

#define	_SVID_SOURCE	/* for lrand48 and friends */	
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <string.h>

#include "gf2x.h"

#ifdef TUNE_MUL1
#include INLINES_FILE
#endif

#define N 2000000

int cputime()
{
    struct rusage rus;

    getrusage(0, &rus);
    return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

int main(int argc, char *argv[])
{
    unsigned long n, m;
    unsigned long i, *c, *a, *b;
#ifndef TUNE_MUL1
    unsigned long *stk;
#endif
    int st;
    int nokara = 0;
    double tries;

    if (strcmp(argv[1], "-nokara") == 0) {
	nokara = 1;
	argc--;
	argv++;
    }

    n = atoi(argv[1]);

    m = N / n;
    tries = (double) N / (double) n / 1e6;

    a = (unsigned long *) malloc((n + m) * sizeof(unsigned long));
    b = (unsigned long *) malloc((n + m) * sizeof(unsigned long));
    c = (unsigned long *) malloc(2 * n * sizeof(unsigned long));

    /* suggestion from Pierrick Gaudry: call the routines on different 
       operands to avoid branch prediction optimizations (in particular
       for n=1) */
    for (i = 0; i < n + m; i++) {
	a[i] = (unsigned long) lrand48();
	b[i] = (unsigned long) lrand48();
    }

    if (!nokara)
	printf("n=%lu: ", n);

    st = cputime();
    for (i = 0; i < m; i++)
	mul_basecase(c, a + i, n, b + i, n);
    printf("mul took %1.0f ns", (cputime() - st) / tries);

#ifndef TUNE_MUL1
    stk = (unsigned long *) malloc(toomspace(n) * sizeof(unsigned long));
    if (nokara == 0 && n >= 2) {
	st = cputime();
	for (i = 0; i < m; i++)
	    KarMul(c, a + i, b + i, n, stk);
	printf(", KarMul took %1.0f ns", (cputime() - st) / tries);

	st = cputime();
	for (i = 0; i < m; i++)
	    Toom(c, a + i, b + i, n, stk);
	printf(", Toom took %1.0f ns", (cputime() - st) / tries);
    }
#endif

    printf("\n");

    return 0;
}
