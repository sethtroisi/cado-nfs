/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

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

/* Tuning program for GF(2)[x] low-level multiplication. */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <string.h>

/* gf2x is used for checking correctness -- not otherwise */
#include "gf2x.h"

#include "timing.h"

/* Enough to be measurable. */
#define N 2000000

#ifndef TSIZE
#error "Define TSIZE to be the size of the operands of the test function"
#endif

#define PAD_(X,Y)       X ## Y
#define PAD(X,Y)        PAD_(X,Y)
#define TFUNCTION       PAD(tuning_gf2x_mul,TSIZE)

#if TSIZE == 1
extern void TFUNCTION(unsigned long *c, unsigned long a, unsigned long b);
#else
extern void TFUNCTION(unsigned long *c, const unsigned long *a, const unsigned long *b);
#endif

int main(int argc, char *argv[])
{
    unsigned long i, *c0, *c, *a, *b;
    uint64_t st;

    unsigned int m = N / TSIZE;

    char * progname = "me";
    if (argc >= 1) {
        progname = argv[0];
    }

    if (m < 10) m = 10;

    a = (unsigned long *) malloc((m + TSIZE) * sizeof(unsigned long));
    b = (unsigned long *) malloc((m + TSIZE) * sizeof(unsigned long));
    c = (unsigned long *) malloc(2 * TSIZE * sizeof(unsigned long));
    c0 = (unsigned long *) malloc(2 * TSIZE * sizeof(unsigned long));

    for (i = 0; i < m; i++) {
        a[i] = (unsigned long) rand();
        b[i] = (unsigned long) rand();
    }

    for (i = 0; i < 10 && i < m; i++) {
        /* Use this one as a reference implementation */
        gf2x_mul(c0, a + i, TSIZE, b + i, TSIZE);
#if TSIZE == 1
        TFUNCTION (c, a[i], b[i]);
#else
        TFUNCTION (c, a + i, b + i);
#endif
        if (memcmp(c, c0, 2 * TSIZE * sizeof(unsigned long)) != 0) {
            fprintf(stderr, "Error, computed test values do not match\n");
            exit(255);
        }
    }

    st = microseconds();
    for (i = 0; i < m; i++) {
#if TSIZE == 1
        TFUNCTION (c, a[i], b[i]);
#else
        TFUNCTION (c, a + i, b + i);
#endif
    }

    printf("%s : %1.0f ns\n", progname, (microseconds() - st) / (double) m * 1.0e3);


    free(c0);
    free(c);
    free(b);
    free(a);
    return 0;
}
