#include <stdint.h>

extern void DMatrixSumBit(unsigned long m, unsigned long n,
			  const unsigned long *A, const unsigned long *B,
			  unsigned long *C);
extern void SMultDmatrixBit(unsigned long m, unsigned long n, unsigned long p,
			    unsigned long *a, unsigned long *b,
			    unsigned long *c, unsigned long start,
			    unsigned long size);
extern void TSMultDmatrixBit(unsigned long m, unsigned long n,
			     unsigned long p, unsigned long *a,
			     unsigned long *b, unsigned long *c,
			     unsigned long start, unsigned long size);
extern void TransposeBit(unsigned long m, unsigned long n,
			 const unsigned long *a, unsigned long *b);
extern void SelectColumn(unsigned long m, unsigned long n,
			 const unsigned long *V, unsigned long k,
			 unsigned long *C);
extern unsigned long AddRowBit(unsigned long m, unsigned long n,
			       unsigned long *a, unsigned long u,
			       unsigned long i, unsigned long j);
extern unsigned long SwapRowsBit(unsigned long m, unsigned long n,
				 unsigned long *a, unsigned long i,
				 unsigned long j);
extern unsigned long MultiplyRowBit(unsigned long m, unsigned long n,
				    unsigned long *a, unsigned long u,
				    unsigned long j);
extern unsigned long SelectLinesBit(unsigned long n, const unsigned long *a,
				    unsigned long m, unsigned long *b);
extern void DSumBit(unsigned long m, unsigned long n, unsigned long *a,
		    unsigned long *b, unsigned long *c);
extern void DMultBit(unsigned long m, unsigned long n, unsigned long p,
		     const unsigned long *A, const unsigned long *B,
		     unsigned long *C);
extern void SMultDmatrixBitNew(unsigned long NrLines, unsigned long NrCols_b,
			       unsigned long *a, unsigned long *b,
			       unsigned long *c);
extern void StoDBit(unsigned long m, unsigned long n, unsigned long *a,
		    unsigned long *b);
extern void TSSMultDmatrix(unsigned long m, unsigned long n, unsigned long p,
			   unsigned long *a, unsigned long *b,
			   unsigned long *c);
extern unsigned long GaussElimBit(unsigned long m, unsigned long n,
				  unsigned long *a, unsigned long *b,
				  unsigned long *c, unsigned long *d,
				  unsigned long *Lend);
extern unsigned long SelectLinesListBit(unsigned long n,
					const unsigned long *a,
					unsigned long *m, unsigned long Lengm,
					unsigned long *b);
extern unsigned long InverseList(unsigned long N, unsigned long Lengm,
				 const unsigned long *m, unsigned long *b);






extern void DMultBitOld(unsigned long N, unsigned long m, unsigned long n,
			const unsigned long *A, const unsigned long *B,
			unsigned long *C);
extern void DMultBitOld(unsigned long N, unsigned long m, unsigned long n,
			const unsigned long *A, const unsigned long *B,
			unsigned long *C);
extern unsigned long NumbRows(unsigned long **a);
extern void SMultDmatrix(unsigned long m, unsigned long **a,
			 unsigned long **b, unsigned long **c);
extern void StoD(unsigned long m, unsigned long n, unsigned long **a,
		 unsigned long **b);
extern void DtoS(unsigned long m, unsigned long n, unsigned long **b,
		 unsigned long **a);
extern void DSum(unsigned long m, unsigned long n, unsigned long **a,
		 unsigned long **b, unsigned long **c);
extern void DMult(unsigned long m, unsigned long n, unsigned long p,
		  unsigned long **a, unsigned long **b, unsigned long **c);
extern void Transpose(unsigned long m, unsigned long n, unsigned long **a,
		      unsigned long **b);
extern unsigned long SelectLines(unsigned long n, unsigned long **a,
				 unsigned long m, unsigned long **b);
extern unsigned long AddRow(unsigned long m, unsigned long n,
			    unsigned long a[m][n], unsigned long u,
			    unsigned long i, unsigned long j);
extern unsigned long SwapRows(unsigned long m, unsigned long n,
			      unsigned long a[m][n], unsigned long i,
			      unsigned long j);
extern unsigned long MultiplyRow(unsigned long m, unsigned long n,
				 unsigned long a[m][n], unsigned long u,
				 unsigned long j);


extern unsigned long NumbCoeffSMatrix(unsigned long *a, unsigned long m);
extern void TVUBit(unsigned long m, unsigned long n, const unsigned long *A,
		   const unsigned long *B, unsigned long *C);
extern void VUBit(unsigned long m, unsigned long n, const unsigned long *A,
		  const unsigned long *B, unsigned long *C);
extern void TVUBit_v2(unsigned long m, unsigned long n, uint64_t * A,
		      uint64_t * B, uint64_t * C);
extern inline void addmul(uint64_t a, uint64_t b, uint64_t * w);
