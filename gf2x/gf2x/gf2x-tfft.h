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

#ifndef GF2X_TFFT_H_
#define GF2X_TFFT_H_

#ifdef __cplusplus
extern "C" {
#endif

struct gf2x_tfft_info_s {
    size_t bits_a;  // number of bits of operand1
    size_t bits_b;  // number of bits of operand2
    size_t K;       // 0 indicates fallback.
    size_t M;
    unsigned long * tmp;
    size_t * perm;
    int split;  // boolean
};

typedef struct gf2x_tfft_info_s gf2x_tfft_info_t[1];
typedef struct gf2x_tfft_info_s * gf2x_tfft_info_ptr;
typedef const struct gf2x_tfft_info_s * gf2x_tfft_info_srcptr;

typedef unsigned long gf2x_tfft_t;
typedef gf2x_tfft_t * gf2x_tfft_ptr;
typedef const gf2x_tfft_t * gf2x_tfft_srcptr;

extern size_t gf2x_tfft_size(gf2x_tfft_info_srcptr o);
extern gf2x_tfft_ptr gf2x_tfft_get(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr ptr, size_t k);
extern void gf2x_tfft_zero(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr ptr, size_t n);
extern void gf2x_tfft_cpy(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr y, gf2x_tfft_srcptr x);
extern gf2x_tfft_ptr gf2x_tfft_alloc(gf2x_tfft_info_srcptr o, size_t n);
extern void gf2x_tfft_free(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr ptr, size_t n);
extern void gf2x_tfft_dft(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tr, const unsigned long * a, size_t bits_a);
extern void gf2x_tfft_compose(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tc, gf2x_tfft_srcptr ta, gf2x_tfft_srcptr tb);
extern void gf2x_tfft_addcompose(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tc, gf2x_tfft_srcptr ta, gf2x_tfft_srcptr tb);
extern void gf2x_tfft_add(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tc, gf2x_tfft_srcptr ta, gf2x_tfft_srcptr tb);
extern void gf2x_tfft_ift(gf2x_tfft_info_srcptr o, unsigned long * c, size_t bits_c, gf2x_tfft_ptr tr);
extern void gf2x_tfft_init(gf2x_tfft_info_ptr o, size_t bits_a, size_t bits_b, ...);
extern void gf2x_tfft_init_similar(gf2x_tfft_info_ptr o, size_t bits_a, size_t bits_b, gf2x_tfft_info_srcptr other);
extern int gf2x_tfft_compatible(gf2x_tfft_info_srcptr o1, gf2x_tfft_info_srcptr o2);
extern void gf2x_tfft_clear(gf2x_tfft_info_ptr o);


#ifdef __cplusplus
}
#endif

#endif	/* GF2X_TFFT_H_ */
