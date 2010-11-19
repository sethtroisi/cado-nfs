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
#define _BSD_SOURCE
#define _POSIX_C_SOURCE 200112L
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "tuning-common.h"

double mulstep = 1;
FILE * rp;
const char * outfile = NULL;

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
	   const char * cname, const unsigned long *c,
           const char * dname, const unsigned long *d)
{
    long i = 0;
    for(i = 0 ; i < m + n ; i++) {
        if (c[i] != d[i]) {
            fprintf(stderr,
                    "Error: %s and %s differ for %ldx%ld at word %ld\n",
                    cname, dname, m, n, i);
            if (m + n < 1000) {
                dump(a, m, b, n, c, d);
            }
            abort();
        }
    }
}


void set_tuning_output()
{
    if (outfile) {
        if ((rp = fopen(outfile, "w")) == NULL) {
            fprintf(stderr, "fopen(%s): %s\n", outfile, strerror(errno));
            exit(1);
        }
    } else {
#if 0
        /* Not really useful since there is -o anyway */
        /* If file descriptor # 3 is open for writing, use it. */
        if ((rp = fdopen(3, "w")) == NULL) {
            rp = stdout;
        }
#else
        rp = stdout;
#endif
    }

    setbuf(rp, NULL);
    setbuf(stdout, NULL);
}

int handle_tuning_mulstep(int * p_argc, char *** p_argv)
{
    if (strcmp((*p_argv)[0], "--step") == 0 || strcmp((*p_argv)[0], "-s") == 0) {
        (*p_argc)--,(*p_argv)++;
        if (! (*p_argc)) {
            return -1;
        }
        mulstep = atof((*p_argv)[0]);
        return 1;
    }
    return 0;
}

int handle_tuning_outfile(int * p_argc, char *** p_argv)
{
    if (strcmp((*p_argv)[0], "--output") == 0 || strcmp((*p_argv)[0], "-o") == 0) {
        (*p_argc)--,(*p_argv)++;
        if (! (*p_argc)) {
            return -1;
        }
        outfile = (*p_argv)[0];
        return 1;
    }
    return 0;
}
