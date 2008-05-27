#ifndef  BL_RANDOM_H_
#define  BL_RANDOM_H_

#include <stdint.h>
#include "struct.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void init_random();
extern void close_random();
extern unsigned long random_coeff();
extern unsigned long random_coeff2();
extern void RandomDMatrix(unsigned long m, unsigned long n,
			  unsigned long **a);
extern void RandomSMatrix(unsigned long m, unsigned long n,
			  unsigned long **a);
extern void RandomSymmetricSMatrix(unsigned long m, unsigned long n,
				   unsigned long k, unsigned long *a);
extern void RandomSymetricDMatrix(unsigned long n, unsigned long **a);
extern void RandomDMatrixBit(unsigned long N, unsigned long m,
			     unsigned long n, unsigned long *a);
extern void CreateRandomMatrixFile(char *f, unsigned long N, unsigned long m,
				   unsigned long n, unsigned long k);
extern void CreateRandomMatrix(unsigned long m, unsigned long n,
			       unsigned long k, unsigned long *a,
			       unsigned long dev);
extern void CreateRandomSymmetricMatrix(unsigned long n, unsigned long k,
					unsigned long *a, unsigned long dev);
extern unsigned long *CreateRandomSymmetricMatrixTest(unsigned long m,
						      unsigned long n,
						      unsigned long k);
extern unsigned long *RandomDMatrixBitTest(unsigned long m, unsigned long n,
					   unsigned long *a);

extern void RandomDMatrixBitTest_new(unsigned long m,
				    unsigned long n,DenseMatrix a);



#ifdef __cplusplus
}
#endif

#endif  /* BL_RANDOM_H_ */
