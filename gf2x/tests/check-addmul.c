/* This file is part of the gf2x library.

   Copyright 2010 Paul Zimmermann

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
#include "gf2x/gf2x_mul1.h"

unsigned long longrand()
{
#if RAND_MAX < LONG_MAX
    return rand() * ((long)RAND_MAX+1) + rand();
#else
    return rand();
#endif
}

/* check gf2x_mul_1_n against gf2x_addmul_1_n */
static void
check_gf2x_mul_1_n (long sb)
{
  unsigned long *cp = malloc ((sb + 1) * sizeof (unsigned long));
  unsigned long *refcp = malloc ((sb + 1) * sizeof (unsigned long));
  unsigned long *bp = malloc (sb * sizeof (unsigned long));
  long i;
  unsigned long a;

  for (i = 0; i < sb; i++)
    bp[i] = longrand ();
  a = longrand ();

  /* using gf2x_mul_1_n */
  /* put initial noise into cp */
  for (i = 0; i < sb; i++)
    cp[i] = longrand ();
  cp[sb] = gf2x_mul_1_n (cp, bp, sb, a);

  /* using gf2x_addmul_1_n */
  for (i = 0; i < sb + 1; i++)
    refcp[i] = (unsigned long) 0;
  cp[sb] ^= gf2x_addmul_1_n (refcp, refcp, bp, sb, a);

  for (i = 0; i < sb + 1; i++)
    if (cp[i] != refcp[i])
      {
        fprintf (stderr, "Error, gf2x_mul_1_n and gf2x_addmul_1_n differ\n");
        fprintf (stderr, "bp=");
        for (i = 0; i < sb; i++)
          fprintf (stderr, "%lu ", bp[i]);
        fprintf (stderr, "\na=%lu\n", a);
        fprintf (stderr, "gf2x_mul_1_n gives");
        for (i = 0; i < sb + 1; i++)
          fprintf (stderr, " %lu", cp[i]);
        fprintf (stderr, "\ngf2x_addmul_1_n gives");
        for (i = 0; i < sb + 1; i++)
          fprintf (stderr, " %lu", refcp[i]);
        fprintf (stderr, "\n");
        printf(" FAILED\n");
        exit (1);
      }
}


int main(int argc, char * argv[])
{
    if (argc != 1) {
        fprintf(stderr, "Usage: ./check-addmul\n");
        exit(1);
    }

    int sb;
    printf ("Checking gf2x_mul_1_n against gf2x_addmul_1_n...");
    fflush(stdout);
    for (sb = 1; sb <= 1000; sb++)
      check_gf2x_mul_1_n (sb);
    printf ("ok\n");

    return 0;
}

