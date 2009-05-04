#ifndef OPS_POLY_FFT_H_
#define OPS_POLY_FFT_H_

#include <gmp.h>
#include <list>
#include "structure.h"
#include "auxfuncs.h"

class ops_poly_fft {
	int s;
	static int coeff_stride;
public:
	typedef ops_poly_fft self;
	typedef bw_mbdft mat_mb;
	typedef bw_bbdft mat_bb;
	typedef mp_limb_t * transform_t;

	void zero(transform_t) const;
	void one(transform_t) const;
	void convolution(transform_t, transform_t, transform_t) const;
#define	HAS_CONVOLUTION_SPECIAL
	void convolution_special(transform_t, transform_t, transform_t,
			unsigned int dg_kill) const;
	void transform(transform_t, mp_limb_t *, ptrdiff_t, int) const;
	void itransform(mp_limb_t *, ptrdiff_t, int, transform_t) const;
	operator int() const;
	void set(int);
	static void init(unsigned int nmax, std::list<char *> const& args);
	static void cleanup();

	bool operator==(self const&) const;
	bool fits(int) const;
	self& operator=(self const&);

	void * mat_mb_alloc(mat_mb& x) const;
	void mat_mb_free(mat_mb x) const;
	transform_t mat_mb_get(mat_mb x, int i, int j) const;

	void * mat_bb_alloc(mat_bb& x) const;
	void mat_bb_free(mat_bb x) const;
	transform_t mat_bb_get(mat_bb x, int i, int j) const;
};

#define	HAS_CONVOLUTION_SPECIAL

#endif	/* OPS_POLY_FFT_H_ */
