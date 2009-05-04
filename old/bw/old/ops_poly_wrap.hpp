#ifndef OPS_POLY_WRAP_H_
#define OPS_POLY_WRAP_H_

#include "ops_poly_fft.hpp"
#include "ops_poly_schoolbook.hpp"
#ifdef	PATCHED_GMP
#include "ops_poly_integerfft.hpp"
#endif

class ops_poly_wrap {
	ops_poly_fft fft;
	ops_poly_scbk scbk;
#ifdef	PATCHED_GMP
	ops_poly_ifft ifft;
#endif
public:
        typedef ops_poly_wrap self;
	struct mat_mb {
		ops_poly_fft::mat_mb fft;
		ops_poly_scbk::mat_mb scbk;
#ifdef	PATCHED_GMP
		ops_poly_ifft::mat_mb ifft;
#endif
	};
	struct mat_bb {
		ops_poly_fft::mat_bb fft;
		ops_poly_scbk::mat_bb scbk;
#ifdef	PATCHED_GMP
		ops_poly_ifft::mat_bb ifft;
#endif
	};
	struct transform_t {
		ops_poly_fft::transform_t fft;
		ops_poly_scbk::transform_t scbk;
#ifdef	PATCHED_GMP
		ops_poly_ifft::transform_t ifft;
#endif
	};
        inline void zero(transform_t x) const {
		fft.zero(x.fft);
		scbk.zero(x.scbk);
#ifdef	PATCHED_GMP
		ifft.zero(x.ifft);
#endif
	}
        inline void one(transform_t x) const {
		fft.one(x.fft);
		scbk.one(x.scbk);
#ifdef	PATCHED_GMP
		ifft.one(x.ifft);
#endif
	}
        inline void convolution(transform_t r,
			transform_t p, transform_t q) const
	{
		fft.convolution(r.fft, p.fft, q.fft);
		scbk.convolution(r.scbk, p.scbk, q.scbk);
#ifdef	PATCHED_GMP
		ifft.convolution(r.ifft, p.ifft, q.ifft);
#endif
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
#ifdef	PATCHED_GMP
		ifft.transform(dst.ifft, src, stride, deg);
#endif
	}
        void itransform(mp_limb_t *, ptrdiff_t, int, transform_t) const;
        operator int() const { return (int) scbk; }
        void set(int n) {
		fft.set(n);
		scbk.set(n);
#ifdef	PATCHED_GMP
		ifft.set(n);
#endif
	}
        static void init(unsigned int nmax, std::list<char *> const& args) {
		ops_poly_fft::init(nmax, args);
		ops_poly_scbk::init(nmax, args);
#ifdef	PATCHED_GMP
		ops_poly_ifft::init(nmax, args);
#endif
	}
        static void cleanup() {
		ops_poly_fft::cleanup();
		ops_poly_scbk::cleanup();
#ifdef	PATCHED_GMP
		ops_poly_ifft::cleanup();
#endif
	}
        inline bool operator==(self const& o) const {
		return fft == o.fft
			&& scbk == o.scbk
#ifdef	PATCHED_GMP
			&& ifft == o.ifft
#endif
			;
	}
        bool fits(int n) const {
		return fft.fits(n) && scbk.fits(n)
#ifdef	PATCHED_GMP
			&& ifft.fits(n)
#endif
			;
	}
        self& operator=(self const& o) {
		fft = o.fft;
		scbk = o.scbk;
#ifdef	PATCHED_GMP
		ifft = o.ifft;
#endif
		return *this;
	}

        inline void * mat_mb_alloc(mat_mb& x) const {
		void * res;
		fft.mat_mb_alloc(x.fft);
		res =
		scbk.mat_mb_alloc(x.scbk);
#ifdef	PATCHED_GMP
		ifft.mat_mb_alloc(x.ifft);
#endif
		return res;
	}
        inline void mat_mb_free(mat_mb x) const {
		fft.mat_mb_free(x.fft);
		scbk.mat_mb_free(x.scbk);
#ifdef	PATCHED_GMP
		ifft.mat_mb_free(x.ifft);
#endif
	}
        inline transform_t mat_mb_get(mat_mb x, int i, int j) const {
		transform_t res;
		res.fft = fft.mat_mb_get(x.fft, i, j);
		res.scbk = scbk.mat_mb_get(x.scbk, i, j);
#ifdef	PATCHED_GMP
		res.ifft = ifft.mat_mb_get(x.ifft, i, j);
#endif
		return res;
	}
        inline void * mat_bb_alloc(mat_bb& x) const {
		void * res;
		fft.mat_bb_alloc(x.fft);
		res =
		scbk.mat_bb_alloc(x.scbk);
#ifdef	PATCHED_GMP
		ifft.mat_bb_alloc(x.ifft);
#endif
		return res;
	}
        inline void mat_bb_free(mat_bb x) const {
		fft.mat_bb_free(x.fft);
		scbk.mat_bb_free(x.scbk);
#ifdef	PATCHED_GMP
		ifft.mat_bb_free(x.ifft);
#endif
	}
        inline transform_t mat_bb_get(mat_bb x, int i, int j) const {
		transform_t res;
		res.fft = fft.mat_bb_get(x.fft, i, j);
		res.scbk = scbk.mat_bb_get(x.scbk, i, j);
#ifdef	PATCHED_GMP
		res.ifft = ifft.mat_bb_get(x.ifft, i, j);
#endif
		return res;
	}
};
#endif	/* OPS_POLY_WRAP_H_ */
