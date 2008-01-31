
#include "basefield.hpp"
#include "manu.h"
#include "bitstring.hpp"

#include <gmp.h>

#include <istream>
#include <ostream>

using namespace std;


unsigned long bitrev(unsigned long a)
{
    unsigned long m;
#if GMP_LIMB_BITS == 64
    a = (a >> 32) | (a << 32);
    m = 0xffff0000ffff0000UL;
    a = (a & m) >> 16 | (a << 16 & m);
    m ^= m >> 8;
    /* m = 0xff00ff00ff00ff00UL */
#elif GMP_LIMB_BITS == 32
    a = a >> 16 | a << 16;
    m = 0xff00ff00ff00ff00UL
#else
#error "Bad ULONG_BITS"
#endif
    a = (a & m) >> 8 | (a << 8 & m);
    m ^= m >> 4;
    /* m = 0xf0f0f0f0f0f0f0f0UL */
    a = (a & m) >> 4 | (a << 4 & m);
    m ^= m >> 2;
    a = (a & m) >> 2 | (a << 2 & m);
    m ^= m >> 1;
    a = (a & m) >> 1 | (a << 1 & m);
    return a;
}

/* This is going to be enabled someday when the vector type for
 * bitstrings specifies 128-bit values.
 */
#if 0
typedef uint64_t two64 __attribute__((vector_size(16)));
typedef uint32_t four32 __attribute__((vector_size(16)));
typedef uint16_t eight16 __attribute__((vector_size(16)));
void bitrev_sse2(two64 * x)
{
    *x =__builtin_ia32_pslldqi128(*x, 8) ^
        __builtin_ia32_psrldqi128(*x, 8);
    *x =__builtin_ia32_psllqi128(*x, 32) ^
        __builtin_ia32_psrlqi128(*x, 32);
    *x =(two64) (
        __builtin_ia32_pslldi128((four32) *x, 16) ^
        __builtin_ia32_psrldi128((four32) *x, 16));
    *x =(two64) (
        __builtin_ia32_psllwi128((eight16) *x, 8) ^
        __builtin_ia32_psrlwi128((eight16) *x, 8));
    eight16 m = { ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, };
    m ^= __builtin_ia32_psllwi128(m, 8);
    m ^= __builtin_ia32_psllwi128(m, 4);
    *x =(two64) (
        __builtin_ia32_psllwi128(m & (eight16) *x, 4) ^
        __builtin_ia32_psrlwi128((eight16) *x, 4) & m);
    m ^= __builtin_ia32_psllwi128(m, 2);
    *x =(two64) (
        __builtin_ia32_psllwi128(m & (eight16) *x, 2) ^
        __builtin_ia32_psrlwi128((eight16) *x, 2) & m);
    m ^= __builtin_ia32_psllwi128(m, 1);
    *x =(two64) (
        __builtin_ia32_psllwi128(m & (eight16) *x, 1) ^
        __builtin_ia32_psrlwi128((eight16) *x, 1) & m);
}
#endif


void binary_field::vec_reverse(binary_field::vec_t b, binary_field::vec_t const a, unsigned int n)
{
    unsigned int nw = BITS_TO_WORDS(n, ULONG_BITS);
    if (a == b) {
        unsigned long t;
        unsigned int i
        for(i = 0 ; i < nw - 1 - i ; i++) {
            t = bitrev(a[i]);
            b[i] = bitrev(a[nw-1-i]);
            b[nw-1-i] = t;
        }
        if (i == nw-1-i) {
            b[i] = bitrev(b[i]);
        }
    } else {
        for(unsigned int i = 0 ; i < nw; i++) {
            b[nw-i] = bitrev(a[i]);
        }
    }
    unsigned int rs = nw * ULONG_BITS - n;
    if (rs) {
        mpn_rshift(b, b, rs, nw);
    }
}

void prime_field_any::vec_reverse(prime_field_any::vec_t b, prime_field_any::vec_t const a, unsigned int n)
{
    if (a == b) {
        for(unsigned int i = 0 ; i < n - 1 - i ; i++) {
            swap(b[i], b[n-1-i]);
        }
    } else {
        for(unsigned int i = 0 ; i < n; i++) {
            b[n-1-i] = a[i];
        }
    }
}

unsigned int binary_field::vec_read(std::istream& is, binary_field::vec_t a, unsigned int n)
{
    return read_hexstring(is, (unsigned long *) a, n);
}

unsigned int binary_field::vec_write(std::istream& is, binary_field::vec_t const a, unsigned int n)
{
    return write_hexstring(is, (unsigned long const *) a, n);
}

unsigned int prime_field_any::vec_read(std::istream& is, prime_field_any::vec_t a, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        if (!(is >> a[i]))
            return i;
    }
    return n;
}

unsigned int prime_field_any::vec_write(std::ostream& o, prime_field_any::vec_t const a, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        if (!(o >> a[i]))
            return i;
    }
    return n;
}


