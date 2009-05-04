#ifndef FIELD_USAGE_H
#define FIELD_USAGE_H

#include "field_def.h"

#ifdef	__cplusplus
extern "C" {
#endif

/* There are two fields actually concerned by this module:
 *
 * The base field k.
 * The extension field l, in which lie the required 2^i-th roots of unity.
 *
 * k and l might be equal.
 * More often, l is a finite extension of k.
 *
 * Separate modules can exist, one for each type of definition field.
 *
 * The current implementation, although trivially extendable, is limited
 * to the case where l is a quadratic extension of k.
 */

/* Used structures in the program */
extern struct field * field_k,* field_l;

extern void restrict_scalars(mp_limb_t *, mp_limb_t *);
extern void extend_scalars(mp_limb_t *, mp_limb_t *);
extern void prepare_fields_for_fft(int,int);
extern mp_limb_t * fetch_primitive_root(int);

#define k_add(a,b,c)		field_k->add(a,b,c,field_k)
#define k_sub(a,b,c)		field_k->sub(a,b,c,field_k)
#define k_neg(a,b)		field_k->neg(a,b,field_k)
#define k_addto(a,b)		field_k->addto(a,b,field_k)
#define k_mul(a,b,c)		field_k->mul(a,b,c,field_k)
#define k_addmul(a,b,c)		field_k->addmul(a,b,c,field_k)
#define k_inv(a,b)		field_k->inv(a,b,field_k)
#define k_set_x(a,x)		field_k->set_x(a,x,field_k)
#define k_set_mpn(a,b,n)	field_k->set_mpn(a,b,n,field_k)
#define k_set_mpz(a,x)		field_k->set_mpz(a,x,field_k)
#define k_is_x(a,x)		field_k->is_x(a,x,field_k)
#define k_exp(a,b,c,n)		field_k->exp(a,b,c,n,field_k)
#define k_eq(a,b)		field_k->eq(a,b,field_k)
#define k_reduce(a)		field_k->reduce(a,field_k)
#define k_conjugate(a,b)	field_k->conjugate(a,b,field_k)
#define k_print(a)		field_k->print(a,field_k)
#define	k_size			field_k->size
#define k_set_zero(x)		field_k->set_int(x,0,field_k)
#define k_set_one(x)		field_k->set_int(x,1,field_k)
#define k_set_int(x,a)		field_k->set_int(x,a,field_k)
#define k_is_zero(x)		field_k->is_int(x,0,field_k)
#define k_is_one(x)		field_k->is_int(x,1,field_k)
#define k_set(a,b)		field_k->set(a,b,field_k)

#define l_add(a,b,c)		field_l->add(a,b,c,field_l)
#define l_sub(a,b,c)		field_l->sub(a,b,c,field_l)
#define l_neg(a,b)		field_l->neg(a,b,field_l)
#define l_addto(a,b)		field_l->addto(a,b,field_l)
#define l_mul(a,b,c)		field_l->mul(a,b,c,field_l)
#define l_addmul(a,b,c)		field_l->addmul(a,b,c,field_l)
#define l_inv(a,b)		field_l->inv(a,b,field_l)
#define l_set_x(a,x)		field_l->set_x(a,x,field_l)
#define l_set_mpn(a,b,n)	field_l->set_mpn(a,b,n,field_l)
#define l_set_mpz(a,x)		field_l->set_mpz(a,x,field_l)
#define l_is_x(a,x)		field_l->is_x(a,x,field_l)
#define l_exp(a,b,c,n)		field_l->exp(a,b,c,n,field_l)
#define l_eq(a,b)		field_l->eq(a,b,field_l)
#define l_reduce(a)		field_l->reduce(a,field_l)
#define l_conjugate(a,b)	field_l->conjugate(a,b,field_l)
#define l_print(a)		field_l->print(a,field_l)
#define	l_size			field_l->size
#define l_set_zero(x)		field_l->set_int(x,0,field_l)
#define l_set_one(x)		field_l->set_int(x,1,field_l)
#define l_set_int(x,a)		field_l->set_int(x,a,field_l)
#define l_is_zero(x)		field_l->is_int(x,0,field_l)
#define l_is_one(x)		field_l->is_int(x,1,field_l)
#define l_set(a,b)		field_l->set(a,b,field_l)

#ifdef	__cplusplus
}
#endif

#endif	/* FIELD_USAGE_H */
