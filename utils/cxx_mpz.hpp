#ifndef CXX_MPZ_HPP_
#define CXX_MPZ_HPP_

struct cxx_mpz {
    mpz_t x;
    cxx_mpz() { mpz_init(x); }
    ~cxx_mpz() { mpz_clear(x); }
    cxx_mpz(cxx_mpz const & o) {
        mpz_init(x);
        mpz_set(x, o.x);
    }
    cxx_mpz & operator=(cxx_mpz const & o) {
        mpz_set(x, o.x);
        return *this;
    }
#if __cplusplus >= 201103L
    cxx_mpz(cxx_mpz && o) {
        mpz_init(x);
        mpz_swap(x, o.x);
    }
    cxx_mpz& operator=(cxx_mpz && o) {
        mpz_swap(x, o.x);
        return *this;
    }
#endif
    operator mpz_ptr() { return x; }
    operator mpz_srcptr() const { return x; }
    mpz_ptr operator->() { return x; }
    mpz_srcptr operator->() const { return x; }
};

#endif	/* CXX_MPZ_HPP_ */
