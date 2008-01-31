#ifndef VECTORING_STACK_HPP_
#define VECTORING_STACK_HPP_

#include <gmp.h>
#include <gmpxx.h>

/* This vectoring code does nothing apart allowing any number of y
 * vectors which is a multiple of the underlying vectoring constraint */


template<typename traits>
struct vectoring_stack : private traits {
    typedef typename traits::scalar_t * scalar_t;
    typedef typename traits::wide_scalar_t * wide_scalar_t;
    using traits::max_accumulate;
    using traits::max_accumulate_wide;
    using traits::min_y;

    using traits::field;
    using traits::elt;
    using traits::elt_ur;

    std::string name() const {
        std::ostringstream os;
        os << "multiplexed [" << traits::name() << "]";
    }

    unsigned int stride;

    vectoring_stack() : stride(globals::nbys / min_y) {}

    void init(scalar_t ** p, unsigned int n) const
        { *p = new scalar_t[n]; }
    void clear(scalar_t ** p, unsigned int n) const
        { delete[] *p; *p = NULL; }
    void init_wide(wide_scalar_t ** p, unsigned int n) const
        { *p = new wide_scalar_t[n]; }
    void clear_wide(wide_scalar_t ** p, unsigned int n) const
        { delete[] *p; *p = NULL; }
    scalar_t & get(scalar_t * p, unsigned int i) const
        { return p[i]; }
    scalar_t const & get(const scalar_t * p, unsigned int i) const
        { return p[i]; }
    // not duplicating for wide get */
    int can() const { return globals::nbys == 1; }
    inline void get_y(elt& r, scalar_t const & x, unsigned int i) const
        { BUG_ON(i != 0); r = x; }
    inline void reduce(scalar_t & d, wide_scalar_t & s) const {
        mpz_fdiv_r(d.get_mpz_t(), s.get_mpz_t(), globals::modulus.get_mpz_t());
    }
    inline void reduce(scalar_t * dst, wide_scalar_t * src,
            unsigned int i0, unsigned int i1) const
    {
        for(unsigned int i = i0 ; i < i1 ; i++) {
            reduce(get(dst,i), get(src,i - i0));
        }
    }
    inline void addmul(wide_scalar_t & dst, scalar_t const & src, int32_t x) const {
        dst += src * x;
    }

    inline void zero(scalar_t & x) const { x = 0; }
    // inline void zero(scalar_t & x) const { x = 0; }
    inline void zero(scalar_t * p, unsigned int n) const
        { for(unsigned int i = 0 ; i < n ; i++) { zero(get(p,i)); } }
    inline void copy(scalar_t * q, const scalar_t * p, unsigned int n) const
        { for(unsigned int i = 0 ; i < n ; i++) { get(q,i) = get(p,i); }
    // not duplicating for wide copy
    inline bool is_zero(scalar_t const& x) const
        { return x == 0; }
    inline void assign(std::vector<elt>& z, scalar_t const & x) const
        { BUG_ON(z.size() != 1); z[0] = x; }
    inline void addmul_wide_fold(elt & dst,
            scalar_t const & a,
            scalar_t const & b) const
        { elt x = dst + a * b; reduce(dst, x); }
            /* This one affects from a vector. We provide the
             * position, in case we're doing SIMD.
             */
    inline void assign(scalar_t & x,
            std::vector<elt> const& z, unsigned int i) const
        { x = z[i]; }
    inline std::ostream& print(std::ostream& o, scalar_t const& x) const
        { return o << x; }
    inline std::istream& get(std::istream& i, scalar_t & x) const
        { return i >> x; }
};

#endif	/* VECTORING_STACK_HPP_ */
