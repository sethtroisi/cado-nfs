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
#include "manu.h"

template<typename T>
struct binary_pod_traits {
	typedef matrix_repr_binary representation;
	static const unsigned int max_accumulate = UINT_MAX;
	static const unsigned int max_accumulate_wide = UINT_MAX;

	/* This is unfortunately *NOT* what it should */
	static const unsigned int nbits = sizeof(T) * CHAR_BIT;
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
	static inline int get_y(scalar_t const & x, unsigned int i) {
		BUG_ON(i >= nbits);
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
	static inline void addmul_wide(uint64_t & dst,
			scalar_t const & a,
			scalar_t const & b)
	{
		uint64_t x = a.p & b.p;
#if (GMP_LIMB_BITS == 64)
                dst = __builtin_parityl(x);
#elif (GMP_LIMB_BITS == 32)
                dst = __builtin_parityll(x);
#else
#error "ouch"
#endif
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
		for(unsigned int j = 0 ; j < nbits ; j++)
			x.p ^= (T) (z[i+j] != 0) << j;
	}
	static inline void assign(std::vector<mpz_class>& z, scalar_t const & x) {
		BUG_ON(z.size() != nbits);
		for(unsigned int i = 0 ; i < nbits ; i++) {
			z[i] = (x.p >> i) & (T) 1;
		}
	}

        static std::ostream& print(std::ostream& o, scalar_t const& x)
        {
            T mask = 1;
            for(unsigned int i = 0 ; i < nbits ; i++) {
                if (i) o << " ";
                o << ((x.p & mask) != 0);
                mask <<=1;
            }
            return o;
        }

        static std::istream& get(std::istream& o, scalar_t & x)
        {
            T z = 0;
            unsigned long v;
            for(unsigned int i = 0 ; i < nbits ; i++) {
                o >> v;
                z |= (v << i);
            }
            x.p = z;
            return o;
        }
};

typedef binary_pod_traits<ulong> binary_ulong_traits;
typedef binary_pod_traits<uint64_t> binary_uint64_traits;
typedef binary_pod_traits<uint32_t> binary_uint32_traits;
typedef binary_pod_traits<uint16_t> binary_uint16_traits;
typedef binary_pod_traits<uint8_t> binary_uint8_traits;

#endif	/* BINARY_ULONG_TRAITS_HPP_ */
