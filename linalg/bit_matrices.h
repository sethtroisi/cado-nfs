#ifndef BIT_MATRICES_H_
#define BIT_MATRICES_H_

#include <stdint.h>
typedef uint64_t mat64[64];
typedef uint64_t * mat64_ptr;
typedef const uint64_t * mat64_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

/* This takes, in row major order, an Nx64 matrix A (transpose of a 64xN
 * matrix), together with another Nx64 matrix B, and xors the output
 * matrix with the product transpose(A)*B -- this may as well be seen as
 * the block dot product of A and B.
 */
extern void addmul_TN64_N64(mat64 b, uint64_t * A, uint64_t * x, unsigned int ncol);

/* This transposes a 64x64 matrix. Relatively slow */
extern void transp_6464(mat64 dst, mat64 src, int mask);

/* This multiplies A by B. The chosen function is optimal for N about
 * 20000. At N=2000000, a twice faster version can be obtained. However,
 * it's not critical for cado, so we stick with the slower version.
 */
extern void addmul_N64_6464(uint64_t *C,
                   const uint64_t *A,
                   mat64 B, unsigned long m);


/* this one is simple enough to figure out */
extern void mul_6464_6464(mat64_ptr C, mat64_srcptr A, mat64_srcptr B);

/* No other interface is exposed for now, although some more exist in the
 * development source */

#ifdef __cplusplus
}
#endif

#endif	/* BIT_MATRICES_H_ */
