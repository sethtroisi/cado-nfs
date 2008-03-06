#ifndef FFT_ADAPTER_HPP_
#define FFT_ADAPTER_HPP_


#define CAT(X,Y) X ## Y

#define DEFINE_FFT_ADAPTER(visible,Z)					\
struct visible {							\
	CAT(Z,info_t) o;						\
	typedef CAT(Z,t) t;						\
	typedef CAT(Z,src_t) src_t;					\
									\
	visible() {}							\
	visible(int d1, int d2) { CAT(Z,setup)(o, d1, d2); }		\
	visible(int d1, int d2,					\
			int d3 MAYBE_UNUSED,				\
			int acc MAYBE_UNUSED)			\
	{								\
		CAT(Z,setup)(o, d1, d2);				\
	}								\
									\
	inline t alloc(int n) const { return CAT(Z,alloc)(o, n); }	\
	inline void free(t x, int n) const { CAT(Z,free)(o,x,n); }	\
	inline void zero(t x, int n) const { CAT(Z,zero)(o,x,n); }	\
	inline t get(t x, int n) const { return CAT(Z,get)(o,x,n); }	\
	inline void dft(t x, unsigned long * F, int d) const {		\
		return CAT(Z,dft)(o,x,F,d);				\
	}								\
	inline void ift(unsigned long * F, int d, src_t x) const {	\
		return CAT(Z,ift)(o,F,d,x);				\
	}								\
	inline void compose(t y, src_t x1, src_t x2) const {		\
		return CAT(Z,compose)(o,y,x1,x2);			\
	}								\
	inline void add(t y, src_t x1, src_t x2) const {		\
		return CAT(Z,add)(o,y,x1,x2);				\
	}								\
};

#include "fake_fft.h"
DEFINE_FFT_ADAPTER(fake_fft, fake_)

#include "cantor128.h"
DEFINE_FFT_ADAPTER(cantor_fft, c128_)


#endif	/* FFT_ADAPTER_HPP_ */
