#ifndef BINARY_ULONG_TRAITS_HPP_
#define BINARY_ULONG_TRAITS_HPP_

#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include "traits_globals.hpp"
#include <ios>
#include <iomanip>
#include "matrix_repr_binary.hpp"
#include <string>
#include <sstream>

template<typename T>
struct binary_pod_traits {
	typedef matrix_repr_binary representation;
	static const int max_accumulate = UINT_MAX;
	static const int max_accumulate_wide = UINT_MAX;

	static const int nbits = sizeof(T) * CHAR_BIT;
	typedef binary_field coeff_field;

	struct scalar_t { T p; };
	typedef scalar_t wide_scalar_t;

	struct name {
		char x[40];
		name() {
			snprintf(x, sizeof(x), "binary POD code [%d bits]",
					nbits);
		}
		const char * operator()() const {
			return x;
		}
	};

	static int can() {
		return globals::nbys == nbits && globals::modulus_ulong == 2;
	}

	/* FIXME -- this used to return an mpz_class */
	static inline int get_y(scalar_t const & x, int i) {
		BUG_ON(i >= nbits || i < 0);
		int bit = (x.p >> i) & 1UL;
		return bit;
	}
	static inline void reduce(scalar_t & d, wide_scalar_t & s) {
		d.p = s.p;
	}
	static inline void reduce(scalar_t * dst, wide_scalar_t * src,
		unsigned int i0, unsigned int i1)
	{
		for(unsigned int i = i0 ; i < i1 ; i++) {
			reduce(dst[i], src[i - i0]);
		}
	}
	static inline void addmul(wide_scalar_t & dst, scalar_t const & src) {
		dst.p ^= src.p;
	}
	/* TODO: It's quite unfortunate to call this function. It does
	 * *NOT* get called in the most critical section, that is, when
	 * doing multiplication. However, the check loop uses it. Since
	 * the check is essentially useless in characteristic two, that's
	 * bad.
	 */
	static inline void addmul(wide_scalar_t & dst, scalar_t const & src, int64_t x)
	{
		ulong mask = -x;
		dst.p ^= mask & src.p;
	}

	static inline void zero(scalar_t * p, unsigned int n) {
		memset(p, 0, n * sizeof(scalar_t));
	}
	static inline void zero(scalar_t & x) { zero(&x, 1); }
	static inline void copy(scalar_t * q, const scalar_t * p, unsigned int n) {
		memcpy(q, p, n * sizeof(scalar_t));
	}
        /*      [duplicate]
	static inline void copy(wide_scalar_t * q, const wide_scalar_t * p, unsigned int n) {
		memcpy(q, p, n * sizeof(wide_scalar_t));
	}
        */


	static inline bool is_zero(scalar_t const& x) {
		return !x.p;
	}
	/*
	static inline void assign(mpz_class& z, scalar_t const & x) {
		MPZ_SET_MPN(z.get_mpz_t(), x.p, width);
	}
	*/
	static inline void addmul_wide(scalar_t & dst,
			scalar_t const & a,
			scalar_t const & b)
	{
		dst.p ^= a.p & b.p;
	}
	static inline void
	assign(scalar_t & x, std::vector<mpz_class> const& z,  unsigned int i)
	{
		/* This one affects from a vector. We provide the
		 * position, in case we're doing SIMD.
		 */
		// WARNING("slow");
		/* FIXME -- should we go up to 128 here, or restrict to
		 * nbys ??? */
		x.p = 0;
		for(int j = 0 ; j < globals::nbys ; j++)
			x.p ^= (T) (z[i+j] != 0) << j;
	}
	static inline void assign(std::vector<mpz_class>& z, scalar_t const & x) {
		BUG_ON(z.size() != ULONG_BITS);
		for(unsigned int i = 0 ; i < ULONG_BITS ; i++) {
			z[i] = (x.p >> i) & (T) 1;
		}
	}

	static std::ostream& print(std::ostream& o, scalar_t const& x) {
		/*
		 * TODO: allow some sort of compressed I/O. The problem
		 * is that it interfers a lot with other stuff. a 4x, or
		 * even 16x reduction of the output size would be nice...
		 */
#if 0
		std::ios_base::fmtflags f(o.flags());
		o << std::setw(16) << std::hex;
		for(int i = 0 ; i < 2 ; i++) {
			if (i) { o << " "; }
			o << foo[i];
		}
		o.flags(f);
#endif
		for(int i = 0 ; i < globals::nbys ; i++) {
			if (i) { o << " "; }
			o << ((x.p >> i) & (T) 1);
		}
		return o;
	}
};

typedef binary_pod_traits<ulong> binary_ulong_traits;
typedef binary_pod_traits<uint64_t> binary_uint64_traits;
typedef binary_pod_traits<uint32_t> binary_uint32_traits;
typedef binary_pod_traits<uint16_t> binary_uint16_traits;
typedef binary_pod_traits<uint8_t> binary_uint8_traits;

#endif	/* BINARY_ULONG_TRAITS_HPP_ */
