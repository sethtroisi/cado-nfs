/*
  Copyright 2007 Pierrick Gaudry.

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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>    
#include <sys/resource.h> 
#include <gmp.h>  // for random

// cputime in millisec.
static int
cputime() {
  struct rusage rus;
  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}



void usage_and_die(char** argv) {
  fprintf(stderr, "usage: %s nmin nmax [-s step|-f factor]\n", argv[0]);
  fprintf(stderr, "  where nmin and nmax are numbers of 64-words\n");
  exit(1);
}

/* usage:
 *   ./bench nmin nmax [step]
 */

int main(int argc, char** argv) {
  int nmin, nmax, n;
  int step = 1;
  double factor = 1.0;
  int additive = 1;
  uint64_t *f, *g, *h;

  if ((argc < 3) || (argc > 5)) {
    usage_and_die(argv);
  }

  nmin = atoi(argv[1]);
  nmax = atoi(argv[2]);
  if (nmin < 33) {
    fprintf(stderr, "Warning: nmin adjusted to the minimal value allowed: 33.\n");
    nmin = 33;
  }
  if (nmin > nmax) {
    fprintf(stderr, "Error: nmin should be less or equal to nmax\n");
    usage_and_die(argv);
  }
  if (argc == 4) {
    usage_and_die(argv);
  }
  if (argc == 5) {
    if (argv[3][0] != '-')
      usage_and_die(argv);
    if (argv[3][1] == 's') {
      additive = 1;
      step = atoi(argv[4]);
    } else {
      additive = 0;
      factor = atof(argv[4]);
    }
  }

  f = (uint64_t *)malloc(nmax*sizeof(uint64_t));
  g = (uint64_t *)malloc(nmax*sizeof(uint64_t));
  h = (uint64_t *)malloc(2*nmax*sizeof(uint64_t));

  mpn_random((mp_limb_t *)f, (sizeof(uint64_t)/sizeof(mp_limb_t))*nmax);
  mpn_random((mp_limb_t *)g, (sizeof(uint64_t)/sizeof(mp_limb_t))*nmax);

 
#ifndef MULCANTOR
#define MULCANTOR mulCantor128
#endif

  n = nmin;
  while (n <= nmax) {
    int tm;
    int i, rep;
    tm = cputime();
    MULCANTOR(h, f, n, g, n);
    tm = cputime()-tm;
    if (tm == 0) 
      rep = 1000;
    else {
      rep = 200 / tm; // to make about 0.2 sec.
      if (rep < 2)
	rep = 2;
    }

    //printf("benching n = %d, repeating %d\n", n, rep);
    tm = cputime();
    for (i = 0; i < rep; ++i)
      MULCANTOR(h, f, n, g, n);
    tm = cputime()-tm;
    printf("%d\t%f\n", n, ((double)tm) / ((double)rep));
    fprintf(stderr, "%d\t%f\n", n, ((double)tm) / ((double)rep));
    if (additive)
      n += step;
    else {
      int new_n = (int)(((double)n)*factor);
      if (new_n <= n)
	n++;
      else
	n = new_n;
    }
  }

  return 0;
}
