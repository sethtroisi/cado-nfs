#ifndef OPS_POLY_WRAP_H_
#define OPS_POLY_WRAP_H_

#include "ops_poly_fft.hpp"
#include "ops_poly_schoolbook.hpp"
#include "ops_poly_integerfft.hpp"

class ops_poly_wrap {
	ops_poly_fft fft;
	ops_poly_scbk scbk;
	ops_poly_ifft ifft;
public:
        typedef ops_poly_wrap self;
	struct mat_mb {
		ops_poly_fft::mat_mb fft;
		ops_poly_scbk::mat_mb scbk;
		ops_poly_ifft::mat_mb ifft;
	};
	struct mat_bb {
		ops_poly_fft::mat_bb fft;
		ops_poly_scbk::mat_bb scbk;
		ops_poly_ifft::mat_bb ifft;
	};
	struct transform_t {
		ops_poly_fft::transform_t fft;
		ops_poly_scbk::transform_t scbk;
		ops_poly_ifft::transform_t ifft;
	};
        inline void zero(transform_t x) const {
		fft.zero(x.fft);
		scbk.zero(x.scbk);
		ifft.zero(x.ifft);
	}
        inline void one(transform_t x) const {
		fft.one(x.fft);
		scbk.one(x.scbk);
		ifft.one(x.ifft);
	}
        inline void convolution(transform_t r,
			transform_t p, transform_t q) const
	{
		fft.convolution(r.fft, p.fft, q.fft);
		scbk.convolution(r.scbk, p.scbk, q.scbk);
		ifft.convolution(r.ifft, p.ifft, q.ifft);
	}
#if 0
        inline void convolution_special(transform_t r,
			transform_t p, transform_t q,
                        unsigned int dg_kill) const
	{
		fft.convolution_special(r.fft, p.fft, q.fft, dg_kill);
		scbk.convolution_special(r.scbk, p.scbk, q.scbk, dg_kill);
		ifft.convolution_special(r.ifft, p.ifft, q.ifft, dg_kill);
	}
#endif
        inline void transform(transform_t dst, mp_limb_t * src,
			ptrdiff_t stride, int deg) const
	{
		fft.transform(dst.fft, src, stride, deg);
		scbk.transform(dst.scbk, src, stride, deg);
		ifft.transform(dst.ifft, src, stride, deg);
	}
        void itransform(mp_limb_t *, ptrdiff_t, int, transform_t) const;
        operator int() const { return (int) scbk; }
        void set(int n) {
		fft.set(n);
		scbk.set(n);
		ifft.set(n);
	}
        static void init(unsigned int nmax, std::list<char *> const& args) {
		ops_poly_fft::init(nmax, args);
		ops_poly_scbk::init(nmax, args);
		ops_poly_ifft::init(nmax, args);
	}
        static void cleanup() {
		ops_poly_fft::cleanup();
		ops_poly_scbk::cleanup();
		ops_poly_ifft::cleanup();
	}
        inline bool operator==(self const& o) const {
		return fft == o.fft
			&& scbk == o.scbk
			&& ifft == o.ifft;
	}
        bool fits(int n) const {
		return fft.fits(n) && scbk.fits(n) && ifft.fits(n);
	}
        self& operator=(self const& o) {
		fft = o.fft;
		scbk = o.scbk;
		ifft = o.ifft;
		return *this;
	}

        inline void * mat_mb_alloc(mat_mb& x) const {
		void * res;
		fft.mat_mb_alloc(x.fft);
		res =
		scbk.mat_mb_alloc(x.scbk);
		ifft.mat_mb_alloc(x.ifft);
		return res;
	}
        inline void mat_mb_free(mat_mb x) const {
		fft.mat_mb_free(x.fft);
		scbk.mat_mb_free(x.scbk);
		ifft.mat_mb_free(x.ifft);
	}
        inline transform_t mat_mb_get(mat_mb x, int i, int j) const {
		transform_t res;
		res.fft = fft.mat_mb_get(x.fft, i, j);
		res.scbk = scbk.mat_mb_get(x.scbk, i, j);
		res.ifft = ifft.mat_mb_get(x.ifft, i, j);
		return res;
	}
        inline void * mat_bb_alloc(mat_bb& x) const {
		void * res;
		fft.mat_bb_alloc(x.fft);
		res =
		scbk.mat_bb_alloc(x.scbk);
		ifft.mat_bb_alloc(x.ifft);
		return res;
	}
        inline void mat_bb_free(mat_bb x) const {
		fft.mat_bb_free(x.fft);
		scbk.mat_bb_free(x.scbk);
		ifft.mat_bb_free(x.ifft);
	}
        inline transform_t mat_bb_get(mat_bb x, int i, int j) const {
		transform_t res;
		res.fft = fft.mat_bb_get(x.fft, i, j);
		res.scbk = scbk.mat_bb_get(x.scbk, i, j);
		res.ifft = ifft.mat_bb_get(x.ifft, i, j);
		return res;
	}
};
#endif	/* OPS_POLY_WRAP_H_ */
