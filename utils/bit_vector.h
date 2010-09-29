#ifndef BIT_VECTOR_H_
#define BIT_VECTOR_H_

#include <stdint.h>
#include <stddef.h>

struct bit_vector_s {
    uint64_t * p;
    size_t n;
};
typedef struct bit_vector_s bit_vector[1];
typedef struct bit_vector_s * bit_vector_ptr;
typedef const struct bit_vector_s * bit_vector_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

extern void bit_vector_init(bit_vector_ptr b, size_t n);
extern void bit_vector_init_set(bit_vector_ptr b, size_t n, int s);
extern void bit_vector_set(bit_vector_ptr b, int s);
extern void bit_vector_clear(bit_vector_ptr b);

extern int bit_vector_getbit(bit_vector_srcptr b, size_t pos);
/* The value returned by the two following functions reflect the _old_
 * value of the flag, which is overwritten by the function */
extern int bit_vector_setbit(bit_vector_ptr b, size_t pos);
extern int bit_vector_clearbit(bit_vector_ptr b, size_t pos);

/* In contrast, this returns the new value */
extern int bit_vector_flipbit(bit_vector_ptr b, size_t pos);

extern void bit_vector_write_to_file(bit_vector_srcptr b, const char * fname);
extern void bit_vector_read_from_file(bit_vector_ptr b, const char * fname);

#ifdef __cplusplus
}
#endif

#endif	/* BIT_VECTOR_H_ */
