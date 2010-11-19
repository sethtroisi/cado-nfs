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
#ifndef THRESHOLDS_H_
#define THRESHOLDS_H_

/* This file contains sensible defaults to start with, which are updated
 * by the tuning program */

/* If you read "placeholder" here, it means that the tuning program has
 * not been run (or has not completed) */
#define	GF2X_TOOM_TUNING_INFO		"placeholder"
#define GF2X_FFT_TUNING_INFO		"placeholder"

/* we rely on the fallback 8*sizeof(unsigned long) as provided by
 * gf2x-impl.h */
/* #undef GF2X_WORDSIZE   */

#define GF2X_STORAGE_CLASS_mul1        static
#define GF2X_STORAGE_CLASS_mul_1_n     static
#define GF2X_STORAGE_CLASS_addmul_1_n  static
#define GF2X_STORAGE_CLASS_mul2        static
#define GF2X_STORAGE_CLASS_mul3        static
#define GF2X_STORAGE_CLASS_mul4        static
#define GF2X_STORAGE_CLASS_mul5        static
#define GF2X_STORAGE_CLASS_mul6        static
#define GF2X_STORAGE_CLASS_mul7        static
#define GF2X_STORAGE_CLASS_mul8        static
#define GF2X_STORAGE_CLASS_mul9        static

/* First size for which KARA is used. Essentially hard-coded, since the
 * sizes up to 9 words are already karatsuba, unrolled. The unrolled
 * routines handle their own temporary storage on the stack.
 */
#define GF2X_MUL_KARA_THRESHOLD	10

/* First size for which TC3W is used. It is assumed that TC3W is used
 * before TC3 kicks in. */
/* must be >= 8 */
#define GF2X_MUL_TOOMW_THRESHOLD		18

/* First size for which TC3 is used. This threshold is informative, it is
 * not used in the code */
/* must be >= 17 */
#define GF2X_MUL_TOOM_THRESHOLD		21

/* First size for which TC4 is used. This threshold is informative, it is
 * not used in the code */
/* must be >= 30 */
#define GF2X_MUL_TOOM4_THRESHOLD		235

/* Size above which TC4 is the only TC variant used. */
#define GF2X_MUL_TOOM4_ALWAYS_THRESHOLD		1600

/* First size for which TC3U is used */
/* must be >= 33 */
#define GF2X_MUL_TOOMU_THRESHOLD		49

/* Size above which TC3U is the only TCU variant used. */
#define GF2X_MUL_TOOMU_ALWAYS_THRESHOLD		2049

/* The default values here are appropriate for the tuning program.
 * Appropriate values are substituted later on. Note that the tuning
 * table is always created with this size, to that re-tuning is possible
 * at any time. */
#define	GF2X_TOOM_TUNING_LIMIT	2048
#define GF2X_BEST_TOOM_TABLE {}

#define	GF2X_BEST_UTOOM_TABLE {}

/* This macro is not what you think, and does not completely belong here.
 * It merely recalls that the FFT code *DOES NOT WORK* below this size.
 * So GF2X_MUL_FFT_TABLE should not wander in this bleak range.
 */
#define GF2X_MUL_FFT_THRESHOLD 28

/* {n, K} means use FFT(|K|) up from n words, */
/* where |K|<3 stands for Toom-Cook 3, K < 0 means use FFT2 */
#undef GF2X_MUL_FFT_TABLE

#endif	/* THRESHOLDS_H_ */

/* vim: set ft=c: */
