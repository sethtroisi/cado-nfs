
extern unsigned long ***Alloc3Test(unsigned long ***a, unsigned long nrows,
				   unsigned long ncolumns);
extern unsigned long **Alloc2Test(unsigned long **a, unsigned long nrows,
				  unsigned long ncolumns);
extern void FreeAlloc2(unsigned long **a, unsigned long nrows);
extern void FreeAlloc3(unsigned long ***a, unsigned long nrows);
extern unsigned long ***Alloc3(unsigned long nrows, unsigned long ncolumns);
extern unsigned long **Alloc2(unsigned long nrows, unsigned long ncolumns);
extern unsigned long **Allocmn3(unsigned long m, unsigned long n);
extern unsigned long *Allocmn(unsigned long m, unsigned long n);
extern void FreeAllocmn3(unsigned long **a);
