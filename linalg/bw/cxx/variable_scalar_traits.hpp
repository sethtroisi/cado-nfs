#ifndef VARIABLE_SCALAR_TRAITS_HPP_
#define VARIABLE_SCALAR_TRAITS_HPP_

#include <gmp.h>
#include <gmpxx.h>
#include "traits_globals.hpp"
#include "matrix_repr_prime.hpp"

struct variable_scalar_traits {
	typedef matrix_repr_prime representation;

	typedef mpz_class scalar_t;
	typedef mpz_class wide_scalar_t;
	static const unsigned int max_accumulate = INT_MAX;
	static const unsigned int max_accumulate_wide = INT_MAX;

	typedef prime_field_any coeff_field;
	struct name {
		const char * operator()() const {
			return "variable-size code (SLOW SLOW SLOW !)";
		}
	};

	static int can() {
		return globals::nbys == 1;
	}

	/*
	static inline mpz_class reduce(wide_scalar_t & t) {
		mpz_class r;
		mpz_fdiv_r(r.get_mpz_t(), t.get_mpz_t(), globals::modulus.get_mpz_t());
		return r;
	}
	*/
	static inline mpz_class get_y(scalar_t const & x, unsigned int i) {
		BUG_ON(i != 0);
		return x;
	}
	static inline void reduce(scalar_t & d, wide_scalar_t & s) {
		mpz_fdiv_r(d.get_mpz_t(), s.get_mpz_t(), globals::modulus.get_mpz_t());
	}
	static inline void reduce(scalar_t * dst, wide_scalar_t * src,
		unsigned int i0, unsigned int i1)
	{
		for(unsigned int i = i0 ; i < i1 ; i++) {
			reduce(dst[i], src[i - i0]);
		}
	}
	static inline void addmul(wide_scalar_t & dst, scalar_t const & src, int32_t x) {
		dst += src * x;
	}

	static inline void zero(mpz_class & x) {
		x = 0;
	}
	static inline void zero(scalar_t * p, unsigned int n) {
		for(unsigned int i = 0 ; i < n ; i++) {
			p[i] = 0;
		}
	}
	static inline void copy(scalar_t * q, const scalar_t * p, unsigned int n) {
		for(unsigned int i = 0 ; i < n ; i++) {
			q[i] = p[i];
		}
	}
        /*       [duplicate]
	static inline void copy(wide_scalar_t * q, const wide_scalar_t * p, unsigned int n) {
		for(unsigned int i = 0 ; i < n ; i++) {
			q[i] = p[i];
		}
	}
        */
	static inline bool is_zero(mpz_class const& x) {
		return x == 0;
	}
	/*
	static inline void wzero(wide_scalar_t & x) {
		core_ops::zero<width+2>(x.p);
	}
	*/
	static inline void assign(std::vector<mpz_class>& z, scalar_t const & x) {
		BUG_ON(z.size() != 1);
		z[0] = x;
	}
	static inline void addmul_wide(scalar_t & dst,
			scalar_t const & a,
			scalar_t const & b)
	{
		mpz_class x = dst + a * b;
		reduce(dst, x);
	}
	static inline void
	assign(scalar_t & x, std::vector<mpz_class> const& z,  unsigned int i)
	{
		/* This one affects from a vector. We provide the
		 * position, in case we're doing SIMD.
		 */
		x = z[i];
	}
        static inline std::ostream& print(std::ostream& o, scalar_t const& x)
        {
            return o << x;
        }

        static inline std::istream& get(std::istream& i, scalar_t & x)
        {
            return i >> x;
        }
};




#endif	/* VARIABLE_SCALAR_TRAITS_HPP_ */
