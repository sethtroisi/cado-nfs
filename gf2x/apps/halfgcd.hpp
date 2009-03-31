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

#ifndef HALFGCD_HPP_
#define HALFGCD_HPP_

#include <NTL/GF2X.h>

extern void mul_gen(NTL::GF2X& c, const NTL::GF2X& a, const NTL::GF2X& b);
extern void FastGCD(NTL::GF2X& d, const NTL::GF2X& u, const NTL::GF2X& v);

#endif	/* HALFGCD_HPP_ */
