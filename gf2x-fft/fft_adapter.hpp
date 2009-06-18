#ifndef FFT_ADAPTER_HPP_
#define FFT_ADAPTER_HPP_


#define CAT(X,Y) X ## Y

/* _setup and _dft, _ift integer arguments correspond to lengths in
 * number of bits of polynomials. For _setup, we want the maximum number
 * of bits of the expected result.
 *
 * _alloc(), _free(), _get(), take number of complete polynomials.
 *
 * the size() result is the number of type-t objects which represent
 * _one_ polynomial. alloc(1) is thus the same as malloc(size())
 */

#define DEFINE_FFT_ADAPTER(visible,Z)					\
struct visible {							\
	CAT(Z,info_t) o;						\
	typedef CAT(Z,t) t;						\
	typedef CAT(Z,ptr) ptr;						\
	typedef CAT(Z,srcptr) srcptr;					\
									\
	visible() {}							\
	visible(int n1, int n2) { CAT(Z,init)(o, n1, n2); }		\
	~visible() { CAT(Z,clear)(o); }                 		\
									\
        /* The extra <acc> argument gives the number of times the */    \
        /* fft is going to be reused. It might be replaced by some */   \
        /* other mechanism at some point. */                            \
									\
	visible(int n1, int n2,				        	\
			int n3 MAYBE_UNUSED,				\
			int acc MAYBE_UNUSED)		        	\
	{								\
		CAT(Z,init)(o, n1, n2, acc);				\
	}								\
	inline ptr alloc(int n) const { return CAT(Z,alloc)(o, n); }	\
	inline void free(ptr x, int n) const { CAT(Z,free)(o,x,n); }	\
	inline void zero(ptr x, int n) const { CAT(Z,zero)(o,x,n); }	\
	inline ptr get(ptr x, int n) const { return CAT(Z,get)(o,x,n); }\
	inline void dft(ptr x, unsigned long * F, int n) const {	\
		return CAT(Z,dft)(o,x,F,n);				\
	}								\
	inline void ift(unsigned long * F, int n, ptr x) const {	\
		return CAT(Z,ift)(o,F,n,x);				\
	}								\
	inline void compose(ptr y, srcptr x1, srcptr x2) const {	\
		return CAT(Z,compose)(o,y,x1,x2);			\
	}								\
	inline void add(ptr y, srcptr x1, srcptr x2) const {		\
		return CAT(Z,add)(o,y,x1,x2);				\
	}								\
	inline void cpy(ptr y, srcptr x) const {		        \
		return CAT(Z,cpy)(o,y,x);				\
	}								\
	inline int size() const {	        	                \
		return CAT(Z,size)(o);	        			\
	}								\
};


#endif	/* FFT_ADAPTER_HPP_ */
