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
#include "bitstring.hpp"

template<typename T>
struct binary_pod_traits {
	typedef matrix_repr_binary representation;
	static const unsigned int max_accumulate = UINT_MAX;
	static const unsigned int max_accumulate_wide = UINT_MAX;

	/* This is unfortunately *NOT* what it should */
	static const unsigned int nbits = sizeof(T) * CHAR_BIT;
	typedef binary_field coeff_field;
        typedef coeff_field::elt Kelt;
        typedef coeff_field::vec_t Kvec_t;

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

	static inline Kelt get_y(scalar_t const & x, unsigned int i) {
		BUG_ON(i >= nbits);
		Kelt bit = (x.p >> i) & 1UL;
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
	assign(scalar_t & x, Kvec_t const& z,  unsigned int i)
	{
            const unsigned long * ptr = z + i / ULONG_BITS;
            i %= ULONG_BITS;
            x.p = 0;
            for(unsigned int stuffed = 0; stuffed < nbits ; ) {
                x.p ^= *ptr++ >> i;
                stuffed += ULONG_BITS-i;
                if (stuffed >= nbits) break;
                x.p ^= *ptr << (ULONG_BITS-i);
                stuffed += i;
            }
        }
        static inline void assign(Kvec_t& z, scalar_t const & x) {
            unsigned int i = 0;
            for(unsigned int stuffed = 0; stuffed < nbits ; ) {
                z[i] = x.p >> stuffed;
                stuffed += ULONG_BITS;
                i++;
            }
        }

        static std::ostream& print(std::ostream& o, scalar_t const& x)
        {
            /* This seems silly, but we're disturbed by uint64_t on
             * 32-bit machines.
             */
            unsigned long v[nbits / ULONG_BITS];
            for(unsigned int i = 0 ; i < nbits ; i += ULONG_BITS) {
                v[i / ULONG_BITS] = x.p >> i;
            }
            field::vec_write(o, v, nbits);
            return o;
        }

        static std::istream& get(std::istream& is, scalar_t & x)
        {
            unsigned long v[nbits / ULONG_BITS];
            field::vec_read(is, v, nbits);
            x.p = 0;
            for(unsigned int i = 0 ; i < nbits ; i += ULONG_BITS) {
                x.p |= ((T) v[i / ULONG_BITS]) << i;
            }
            return is;
        }
};

typedef binary_pod_traits<ulong> binary_ulong_traits;
typedef binary_pod_traits<uint64_t> binary_uint64_traits;
typedef binary_pod_traits<uint32_t> binary_uint32_traits;
typedef binary_pod_traits<uint16_t> binary_uint16_traits;
typedef binary_pod_traits<uint8_t> binary_uint8_traits;

#endif	/* BINARY_ULONG_TRAITS_HPP_ */
