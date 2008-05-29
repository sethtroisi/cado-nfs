#ifndef MPFQ_GF2N_COMMON_H_
#define MPFQ_GF2N_COMMON_H_

/*
  This file is part of the MPFQ library

  Copyright 2007 Pierrick Gaudry and Emmanuel Thome'

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

#ifdef __cplusplus
extern "C" {
#endif

// type for small field char 2
typedef struct {
    int io_type;
} mpfq_2_field_struct;

typedef mpfq_2_field_struct mpfq_2_field[1];
typedef mpfq_2_field_struct * mpfq_2_dst_field;


#ifdef __cplusplus
}
#endif

#endif	/* MPFQ_GF2N_COMMON_H_ */
