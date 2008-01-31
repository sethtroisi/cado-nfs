#ifndef BASEFIELD_HPP_
#define BASEFIELD_HPP_

#include <gmp.h>
#include <gmpxx.h>
#include <istream>
#include <ostream>
#include "manu.h"       /* For BITS_TO_WORDS */

struct binary_field {
    typedef unsigned long elt;
    typedef unsigned long elt_ur;
    typedef unsigned long * vec_t;
    static inline void vec_init(vec_t * v, unsigned int n)
    { *v = new unsigned long[BITS_TO_WORDS(n, ULONG_BITS)]; }
    static inline void vec_clear(vec_t * v, unsigned int n) { delete[] *v; }
    static void vec_reverse(vec_t a, vec_t const b, unsigned int n);
    static unsigned int vec_read(std::istream& o, vec_t, unsigned int);
    static unsigned int vec_write(std::ostream& o, vec_t const, unsigned int);
};

struct prime_field_any {
    typedef mpz_class elt;
    typedef mpz_class elt_ur;
    typedef mpz_class * vec_t;
    static inline void vec_init(vec_t * v, unsigned int n) { *v = new mpz_class[n]; }
    static inline void vec_clear(vec_t * v, unsigned int n) { delete[] *v; }
    static void vec_reverse(vec_t a, vec_t const b, unsigned int n);
    static unsigned int vec_read(std::istream& o, vec_t, unsigned int);
    static unsigned int vec_write(std::ostream& o, vec_t const, unsigned int);
};

#endif	/* BASEFIELD_HPP_ */
