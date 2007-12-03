#ifndef FIELD_DEF_H_
#define FIELD_DEF_H_

#ifdef	__cplusplus
extern "C" {
#endif

/* Lazy reduction is allowed all around in the field manipulation
 * routines.
 *
 * Call ->reduce to be certain of the uniqueness of your results.
 *
 * Equality tests always do it.
 */
struct field {
	void (*add)	(mp_limb_t *, mp_limb_t *, mp_limb_t *,	struct field *);
	void (*sub)	(mp_limb_t *, mp_limb_t *, mp_limb_t *,	struct field *);
	void (*neg)	(mp_limb_t *, mp_limb_t *,		struct field *);
	void (*addto)	(mp_limb_t *, mp_limb_t *,		struct field *);
	void (*mul)	(mp_limb_t *, mp_limb_t *, mp_limb_t *,	struct field *);
	void (*addmul)	(mp_limb_t *, mp_limb_t *, mp_limb_t *,	struct field *);
	void (*inv)	(mp_limb_t *, mp_limb_t *,		struct field *);
	void (*set_int)	(mp_limb_t *, int,			struct field *);
	void (*set_mpn)	(mp_limb_t *, mp_limb_t *, mp_size_t,	struct field *);
	void (*set_mpz)	(mp_limb_t *, mpz_t,			struct field *);
	int  (*is_int)	(mp_limb_t *, int,			struct field *);
	void (*exp)	(mp_limb_t *, mp_limb_t *,
				mp_limb_t *, mp_size_t, struct field *);
	int  (*eq)	(mp_limb_t *, mp_limb_t *,		struct field *);
	void (*set)	(mp_limb_t *, mp_limb_t *,		struct field *);
	void (*reduce)	(mp_limb_t *,				struct field *);
	void (*conjugate)(mp_limb_t *, mp_limb_t *,		struct field *);
	void (*print)	(mp_limb_t *,				struct field *);
	void (*destroy)	(struct field *);
	struct field *	base_field;
	unsigned int degree;
	mp_limb_t *	modulus;	/* only meaningful for prime fields */
	mp_limb_t *	minpoly;	/* only meaningful for normal exts. */
	mp_size_t size;	
};

#ifdef	__cplusplus
}
#endif

#endif	/* FIELD_DEF_H_ */
