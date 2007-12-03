#ifndef FIELD_COMMON_H_
#define FIELD_COMMON_H_

#ifdef	__cplusplus
extern "C" {
#endif

/* Internals needed for the inner field layers. */

#define UNIMPLEMENTED(x) if (1) { fprintf(stderr,"Function " x " unimplemented\n"); exit(100); }

extern void	generic_field_set(mp_limb_t *, mp_limb_t *, struct field *);
extern void	generic_field_destroy(struct field *);
extern int	generic_field_eq(mp_limb_t *, mp_limb_t *, struct field *);
extern void	generic_field_exp(mp_limb_t *, mp_limb_t *,
				mp_limb_t *, mp_size_t, struct field *);

#ifdef	__cplusplus
}
#endif

#endif	/* FIELD_COMMON_H_ */
