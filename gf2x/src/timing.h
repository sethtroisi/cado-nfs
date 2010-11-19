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
#ifndef TIMING_H_
#define TIMING_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

extern uint64_t microseconds();

/* cputime */
static inline int cputime(void) { return (int) microseconds() / 1000; }
static inline double seconds(void) { return (double) microseconds() /1.0e6; }


#ifdef __cplusplus
}
#endif

#endif	/* TIMING_H_ */
