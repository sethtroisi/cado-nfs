#ifndef BL_BASIC_H_
#define BL_BASIC_H_

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned long Refresh_Array(unsigned long m, unsigned long n,
				   unsigned long ***A, unsigned long **a);
extern unsigned long Mult_TV_A_V(unsigned long n, unsigned long N,
				 unsigned long **v, unsigned long **a,
				 unsigned long **Res2);
extern unsigned long bigdeal(unsigned long n, unsigned long N,
			     unsigned long **a, unsigned long **v,
			     unsigned long **X);
extern unsigned long bigdeal_Fast(unsigned long n, unsigned long N,
				  unsigned long **a, unsigned long **v,
				  unsigned long **X);
extern unsigned long LanczosOld(unsigned long m, unsigned long n,
				unsigned long N, unsigned long *a,
				unsigned long *v, unsigned long *X);
extern unsigned long Refresh_ArrayBit(unsigned long m, unsigned long n,
				      unsigned long **A, unsigned long *a);
extern unsigned long TestZero(unsigned long m, unsigned long n,
			      unsigned long *a);

extern char * TimeConvert(float a);

#ifdef __cplusplus
}
#endif

#endif	/* BL_BASIC_H_ */
