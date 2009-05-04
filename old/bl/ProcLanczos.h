#ifndef BL_PROCLANCZOS_H_
#define BL_PROCLANCZOS_H_




#ifdef __cplusplus
extern "C" {
#endif


extern unsigned long *LanczosIterations(unsigned long *a, unsigned long *Y,
					unsigned long m, unsigned long n,
					unsigned long Block,
					unsigned long *Resultado);


extern void LanczosIterations_new(DenseMatrix Result, SparseMatrix M, DenseMatrix Y);

extern void KernelSparse(unsigned long *a, unsigned long *R,
			 unsigned long m, unsigned long n,
			 unsigned long Block, DenseMatrix Ker);


extern void KernelSparse_new(DenseMatrix Ker, SparseMatrix M, DenseMatrix R);

extern void Lanczos(SparseMatrix M, unsigned long Block,DenseMatrix Kernel);

extern void Lanczos_new(DenseMatrix Kernel, SparseMatrix M,
	     unsigned long Block);


#ifdef __cplusplus
}
#endif

#endif	/* BL_PROCLANCZOS_H_ */
