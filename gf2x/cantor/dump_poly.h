#ifndef DUMP_POLY_H_
#define DUMP_POLY_H_
/*
  This file is part of the MPFQ library

  Copyright 2007 Pierrick Gaudry and Emmanuel Thomé

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


#include <gmp.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

static void dump_poly(char * name, int n, int s, mp_limb_t * x)
	__attribute__((unused));

static void dump_poly(char * name, int n, int s, mp_limb_t * x)
{
	int i;
	int j = 0;
	int some = 0;
	printf("%s:=", name);
        for(i = 0 ; i < n ; i++) {
		if ((i-j*s) >= s) {
			j++;
		}
		if (((x[j] >> (i%s)) & 1UL) == 0)
			continue;
                if (some) printf(" + ");
		switch(i) {
			case 0 : printf("1"); break;
			case 1 : printf("z"); break;
			default: printf("z^%d", i); break;
		}
		some = 1;
        }
	if (some == 0) {
		printf("0");
	}
        printf(";\n");
}


#ifdef __cplusplus
}
#endif

#endif	/* DUMP_POLY_H_ */
