

extern unsigned long *LanczosIterations(unsigned long *a, unsigned long *Y,
					unsigned long m, unsigned long n,
					unsigned long Block,
					unsigned long *Resultado);
extern unsigned long *KernelSparse(unsigned long *a, unsigned long *R,
				   unsigned long m, unsigned long n,
				   unsigned long Block, unsigned long *Ker,
				   unsigned long *Ind);
//extern unsigned long Lanczos(struct SparseMatrix M, unsigned long Block, struct DenseMatrix Kernel);

extern struct DenseMatrix Lanczos(struct SparseMatrix M, unsigned long Block);
