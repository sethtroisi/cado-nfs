#ifndef OPS_POLY_INTEGERFFT_HPP_
#define OPS_POLY_INTEGERFFT_HPP_

/* needs a patched gmp */
#include <gmp.h>
#include <list>

class ops_poly_ifft {
	int ncoeffs;
	int n;
	int k;
	static mp_size_t limbs_per_coeff;
public:
        typedef ops_poly_ifft self;
	typedef mp_fft_t * mat_bb;
	typedef mp_fft_t * mat_mb;
	typedef mp_fft_t * transform_t;

        void zero(transform_t) const;
        void one(transform_t) const;
        void convolution(transform_t, transform_t, transform_t) const;
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
#define	xxxHAS_CONVOLUTION_SPECIAL

#endif	/* OPS_POLY_INTEGERFFT_HPP_ */
