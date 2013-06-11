#ifndef	CADO_TYPEDEFS_H_
#define	CADO_TYPEDEFS_H_

#include <stdint.h>

/* data type to store the (p,r) values */
#ifndef p_r_values_size
#define p_r_values_size 32
#endif

/* data type to store the renumber table */
#ifndef index_size
#define index_size 32
#endif

#if p_r_values_size == 32
#define p_r_values_t uint32_t
#else
#define p_r_values_t uint64_t
#endif

#if index_size == 32
#define index_t uint32_t
#else
#define index_t uint64_t
#endif 

#define weight_t uint8_t
#define exponent_t uint8_t

#endif	/* CADO_TYPEDEFS_H_ */
