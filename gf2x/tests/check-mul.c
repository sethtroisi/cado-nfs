/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at
   your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with CADO-NFS; see the file COPYING.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#include "gf2x.h"
#include "gf2x/gf2x-impl.h"

/* Fill the ulong* array a by n1 32-bit values in geometric progression.
 * Returns a value which is appropriate for chaining.
 */
uint32_t fill(unsigned long * a, int n1, uint32_t v, uint32_t start)
{
    int i;
#if GF2X_WORDSIZE == 32
    for(i = 0 ; i < n1 ; i++) { a[i] = v; v *= start; }
#elif GF2X_WORDSIZE == 64
    int N1 = n1 / 2 + (n1 & 1);
    unsigned long w = 1;
    for(i = 0 ; i < N1 ; i++) {
        w = (unsigned long) v * (unsigned long) start;
        a[i] = w << 32 | (unsigned long) v;
        v = w * start;
    }
    if (n1 & 1) { v = w; a[N1-1] &= (1UL << 32) - 1; }
#else
#error "Config problem"
#endif
    return v;
}

#define RED(l,h) do {							\
        /* Compute crc mod x^32 + x^7 + x^6 + x^2 + 1 */		\
        l ^= h ^ h << 2 ^ h << 6 ^ h << 7;				\
        h  = h >> 30 ^ h >> 26 ^ h >> 25;				\
        /* h is at most 7 bits now. */					\
        l ^= h ^ h << 2 ^ h << 6 ^ h << 7;				\
        h = 0;								\
} while (0)

uint32_t crc32(unsigned long * c, int n3)
{
    int i;
    uint32_t v = 0UL;
#if GF2X_WORDSIZE == 32
    for(i = 0 ; i < n3 ; i++) {
        uint32_t l,h;
        l = c[i];
        h = v;
        RED(l,h);
        v = l;
    }
#elif GF2X_WORDSIZE == 64
    int N3 = n3 / 2 + (n3 & 3);
    i = 0;
    for(int j = 0 ; j < N3 ; j++) {
        uint32_t l, h;
        unsigned long cj = c[j];

        h = v;
        l = cj;
        RED(l,h);
        i++;
        v = l;
        if (i == n3) break;

        cj >>= 32;

        h = v;
        l = cj;
        RED(l,h);
        i++;
        v = l;
        if (i == n3) break;
    }
#endif
    return v;
}


int main(int argc, char * argv[])
{
    if (argc != 2 && argc != 3) {
        fprintf(stderr, "Usage: ./check <size1> [<size2>]\n");
        fprintf(stderr, "Output of this program must be"
                " checked against precomputed tables\n");
        exit(1);
    }

    int n1 = atoi(argv[1]);
    int n2 = n1;
    if (argc >= 3) {
        n2 = atoi(argv[2]);
    }
    int n3 = n1 + n2;

    /* n1 and n2 must be understood as a number of 32-bit words */

    int N1 = n1;
    int N2 = n2;

#if GF2X_WORDSIZE == 64
    N1 = N1 / 2 + (N1 & 1);
    N2 = N2 / 2 + (N2 & 1);
#endif

    int N3 = N1 + N2;


    unsigned long * a = malloc(N1 * sizeof(unsigned long));
    unsigned long * b = malloc(N2 * sizeof(unsigned long));
    unsigned long * c = malloc(N3 * sizeof(unsigned long));

    uint32_t start=0x76540123UL;

    start |= 1UL;

    uint32_t v;
    uint32_t check0, check1;

    v = start;
    v = fill(a, n1, v, start);
    v = fill(b, n2, v, start);
    gf2x_mul(c,a,N1,b,N2);
    check0 = crc32(c, n3);

    v = start;
    v = fill(c, n1, v, start);
    gf2x_mul(c,c,N1,b,N2);
    check1 = crc32(c, n3);

    if (check0 != check1) {
        printf("aliasing test failed\n");
        exit(1);
    }

    v = start;
    v = fill(a, n1, v, start);
    v = fill(c, n2, v, start);
    gf2x_mul(c,a,N1,c,N2);
    check1 = crc32(c, n3);

    if (check0 != check1) {
        printf("aliasing test failed\n");
        exit(1);
    }

    printf("%d %d %08" PRIx32 " %08" PRIx32 "\n", n1, n2, check0, start);

    free(a);
    free(b);
    free(c);

    return 0;
}

