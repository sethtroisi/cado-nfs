#ifndef BINARY_SSE2_TRAITS_HPP_
#define BINARY_SSE2_TRAITS_HPP_

#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include "traits_globals.hpp"
#include <ios>
#include <iomanip>
#include "matrix_repr_binary.hpp"
#include "bitstring.hpp"

struct binary_sse2_traits {
	typedef matrix_repr_binary representation;
	static const unsigned int max_accumulate = UINT_MAX;
	static const unsigned int max_accumulate_wide = UINT_MAX;
#if defined(__OpenBSD__) && defined(__x86_64)
	/* I know it's outright stupid */
	typedef long inner_type;
#else
	typedef int64_t inner_type;
#endif
	static const int sse2_vectoring = 2;
	typedef inner_type sse2_scalars __attribute__((vector_size(16)));
	typedef inner_type vec_t[sse2_vectoring] __attribute__((aligned(16)));

	typedef binary_field coeff_field;
        typedef coeff_field::elt Kelt;
        typedef coeff_field::vec_t Kvec_t;


	struct scalar_t { sse2_scalars p; };
	typedef scalar_t wide_scalar_t;

	struct name {
		const char * operator()() const {
			return "binary SSE-2 code [ 128 bits ]";
		}
	};

	static int can() {
		return globals::nbys == 128 && globals::modulus_ulong == 2;
	}

	static inline Kelt get_y(scalar_t const & x, int i) {
		vec_t foo;
		memcpy(foo, &x, sizeof(vec_t));
		BUG_ON(i >= 128 || i < 0);
		Kelt bit = (foo[i >> 6] >> (i & 63)) & 1;
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
	static inline void addmul(wide_scalar_t & dst, scalar_t const & src, inner_type x)
	{
		sse2_scalars mask = (sse2_scalars) { -x, -x, };
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
		vec_t foo;
		memcpy(foo, &x, sizeof(vec_t));
		for(unsigned int i = 0 ; i < 2 ; i++) {
			if (foo[i])
				return false;
		}
		return true;
	}
	/*
	static inline void assign(mpz_class& z, scalar_t const & x) {
		MPZ_SET_MPN(z.get_mpz_t(), x, width);
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
            vec_t foo;
            const unsigned long * ptr = z + i / ULONG_BITS;
            i %= ULONG_BITS;
            foo[0] = 0;
            for(unsigned int stuffed = 0; stuffed < 64 ; ) {
                foo[0] ^= *ptr++ >> i;
                stuffed += ULONG_BITS-i;
                if (stuffed >= 64) break;
                foo[0] ^= *ptr << (ULONG_BITS-i);
                stuffed += i;
            }
            foo[1] = 0;
            for(unsigned int stuffed = 0; stuffed < 64 ; ) {
                foo[1] ^= *ptr++ >> i;
                stuffed += ULONG_BITS-i;
                if (stuffed >= 64) break;
                foo[1] ^= *ptr << (ULONG_BITS-i);
                stuffed += i;
            }
		x.p = (sse2_scalars) { foo[0], foo[1], };
	}
	static inline void assign(std::vector<mpz_class>& z, scalar_t const & x) {
		vec_t foo;
		memcpy(foo, &x, sizeof(vec_t));
		BUG_ON(z.size() != 128);
		for(unsigned int i = 0 ; i < 128 ; i++) {
			z[i] = (foo[i>>6] >> (i & 63)) & 1UL;
		}
	}

        static std::ostream& print(std::ostream& o, scalar_t const& x)
        {
            vec_t foo;
            memcpy(foo, &x, sizeof(vec_t));
            unsigned long v[128 / ULONG_BITS];
            for(unsigned int i = 0 ; i < 64 ; i += ULONG_BITS) {
                v[i / ULONG_BITS] = foo[0] >> i;
            }
            for(unsigned int i = 0 ; i < 64 ; i += ULONG_BITS) {
                v[(i + 64) / ULONG_BITS] = foo[1] >> i;
            }
            write_hexstring(o, v, 128);
            return o;
        }

        static std::istream& get(std::istream& is, scalar_t & x)
        {
            unsigned long v[128 / ULONG_BITS];
            vec_t foo;
            read_hexstring(is, v, 128);
            foo[0] = 0;
            for(unsigned int i = 0 ; i < 64 ; i += ULONG_BITS) {
                foo[0] |= ((inner_type) v[i / ULONG_BITS]) << i;
            }
            foo[1] = 0;
            for(unsigned int i = 0 ; i < 64 ; i += ULONG_BITS) {
                foo[1] |= ((inner_type) v[(i + 64) / ULONG_BITS]) << i;
            }
            memcpy(&x, foo, sizeof(vec_t));
            return is;
        }
};

#endif	/* BINARY_SSE2_TRAITS_HPP_ */
