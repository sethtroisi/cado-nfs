#ifndef BW_SCALAR_H_
#define BW_SCALAR_H_

#include <gmp.h>

/* This files contains the mandatory definitions for handling basic
 * scalar types (elements of Z/NZ).
 *
 * Any multiprecision package can be used in place of the one used here,
 * provided you are able to replace the corresponding lines in this file
 * and in bw_scalar.c */

#ifdef	__cplusplus
extern "C" {
#endif

/* These constants indicate whether the memory allocated for the element
 * should be tight (BW_SCALAR_SHORT), or somewhat larger (BW_SCALAR_LONG)
 * to handle very lazy reduction scheme (p-q).
 *
 * BW_SCALAR_MEDIUM is currently one limb over the small_size */

/* scalars */
#define BW_SCALAR_SHORT		0
#define BW_SCALAR_LONG		1
#define BW_SCALAR_MEDIUM	2

typedef mp_limb_t * bw_scalar;

extern mp_size_t _bw_scalar_alloc(bw_scalar *, int);
#define bw_scalar_alloc(x,t) _bw_scalar_alloc(&(x),t)
#define bw_scalar_free(x) free(x)
#define bw_scalar_set_zero(x) memset(x,0,bw_allocsize*sizeof(mp_limb_t));
#define bw_long_scalar_set_zero(x) memset(x,0,bw_longsize*sizeof(mp_limb_t));
#define bw_scalar_is_zero(x) (mpn_cmp(x,zero,bw_allocsize)==0)
extern int bw_scalar_fits_word(bw_scalar, type32 *, int *);
extern int bw_scalar_read(bw_scalar, FILE *);
extern int bw_scalar_write(FILE *, bw_scalar);
extern void bw_scalar_set_random(bw_scalar);
extern void bw_scalar_set_one(bw_scalar);

/* vectors */
typedef mp_limb_t * bw_vector;
extern mp_size_t _bw_vector_alloc(bw_vector *, coord_t, int);
#define bw_vector_alloc(x,n,t) _bw_vector_alloc(&(x),n,t)
#define bw_vector_free(x,n) free(x)
#define bw_vector_set_zero(x,n) memset(x,0,bw_allocsize*n*sizeof(mp_limb_t));
#define bw_long_vector_set_zero(x,n) memset(x,0,bw_longsize*n*sizeof(mp_limb_t));
extern int bw_vector_read(bw_vector, int, FILE *);
extern int bw_vector_write(FILE *, bw_vector, int);
/*int bw_long_vector_read(bw_vector, int, FILE *);*/
/*int bw_long_vector_write(FILE *, bw_vector, int);*/
extern void bw_vector_set_random(bw_vector, int);

/* lemmings (tiny animals who step through vectors!) */
typedef mp_limb_t * bw_lemming;
#define bw_lemming_reset(x,v) x=v
#define bw_lemming_step(x,n) x+=n*bw_allocsize
#define bw_lemming_clear(x)

/* long lemmings do the same through long vectors */
/*typedef mp_limb_t * bw_long_lemming;
#define bw_long_lemming_reset(x,v) x=v
#define bw_long_lemming_step(x,n) x+=n*bw_longsize
#define bw_long_lemming_clear(x)*/

extern void bw_shorten_long_scalar(bw_scalar);
extern void bw_reduce_short_scalar(bw_scalar);

extern int bw_read_modulus_info(const char *, int);

/* XXX README : This #ifdef is intended to have some files reference
 * bw_allocsize and bw_longsize the way they should. Those files belong
 * to the ``slave'' branch. For bw_allocsize and bw_longsize to be
 * interpreted as macros, the sources must include some header file first
 * with the HARDCODE_PARAMS symbol defined (if wanted). PLEASE NOTE that
 * c/common/bw_scalar.c DOES NOT do this ! The procedures in this piece
 * of code are not speed-critical and thus always refer to the variables.
 * (all occurences of bw_allocsize and bw_longsize should read as
 * computed_bw_allocsize and computed_bw_longsize in this source file).
 *
 * The consistency checking between the variables and the macros has to
 * be done at runtime by some part of the code (e.g. c/slave/slave.c)
 *
 * No part of the code should refer explicitly to either
 * computed_bw_allocsize or HARD_bw_allocsize (same with longsize). Only
 * the name bw_allocsize should be used. (With an exception when setting
 * these variable, or doing the consistency check). XXX*/
#ifdef HARDCODE_PARAMS
#define bw_allocsize HARD_bw_allocsize
#define bw_longsize HARD_bw_longsize
#else
#define bw_allocsize computed_bw_allocsize
#define bw_longsize computed_bw_longsize
#endif
extern mp_size_t computed_bw_allocsize;
extern mp_size_t computed_bw_longsize;

extern mp_size_t bw_filesize;
extern bw_scalar zero;
extern mpz_t modulus;
extern bw_scalar xxmodulus_hbs;
extern mp_size_t shift_amount;
#define modulus_plain	(modulus->_mp_d)
#define xmodulus_hbs	(xxmodulus_hbs+1)
#define modulus_hbs	(xxmodulus_hbs+2)

extern mp_limb_t * two_n_plus_2_mod_p;

#ifdef	__cplusplus
}
#endif

#endif /* BW_SCALAR_H_ */
