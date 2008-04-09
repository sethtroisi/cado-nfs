

extern unsigned long *LanczosIterations(unsigned long *a, unsigned long *Y,
					unsigned long m, unsigned long n,
					unsigned long Block,
					unsigned long *Resultado);



extern void KernelSparse(unsigned long *a, unsigned long *R,
			 unsigned long m, unsigned long n,
			 unsigned long Block, struct DenseMatrix Ker);


extern void Lanczos(struct DenseMatrix Kernel, struct SparseMatrix M,
		    unsigned long Block);
