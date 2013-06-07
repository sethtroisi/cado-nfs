#ifndef	CADO_TYPEDEFS_H_
#define	CADO_TYPEDEFS_H_

#include <inttypes.h>

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
#define PRpr PRIu32
#else
#define p_r_values_t uint64_t
#define PRpr PRIu64
#endif

#if index_size == 32
#define index_t uint32_t
#define PRid PRIu32
#else
#define index_t uint64_t
#define PRid PRIu64
#endif 

/* The weight of ideals saturates at 255 */
/* For relations, we hope that there will never be more */
/* than 255 ideals per relation */
#define weight_t uint8_t
#define exponent_t uint8_t

#endif	/* CADO_TYPEDEFS_H_ */
