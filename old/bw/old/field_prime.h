#ifndef FIELD_PRIME_H_
#define FIELD_PRIME_H_

#ifdef	__cplusplus
extern "C" {
#endif

struct field * new_prime_field(mp_limb_t *, mp_size_t);

#ifdef	__cplusplus
}
#endif

#endif	/* FIELD_PRIME_H_ */
