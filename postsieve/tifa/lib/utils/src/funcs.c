//
// Copyright (C) 2006, 2007 INRIA (French National Institute for Research in
// Computer Science and Control)
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
//

/**
 * \file    funcs.c
 * \author  Jerome Milan
 * \date    Tue Sep 4 2007
 * \version 1.1
 */

 /*
  *  History:
  *    1.1: Tue Sep 4 2007 by JM:
  *          - Added modinv_ui function (modular inverse).
  *          - Added sqrtm_p2 function (modular square root).
  *    1.0: Mon Mar 13 2006 by JM:
  *          - Initial version.
  */

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "first_primes.h"

#include "funcs.h"
#include "macros.h"

/*
 *-----------------------------------------------------------------------------
 *                      Miscellaneous functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static const uint8_t log2table[] = {
    //
    // From public domain code from the Bit Twiddling Hacks web page:
    // http://graphics.stanford.edu/~seander/bithacks.html
    //
  0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
};
//-----------------------------------------------------------------------------
uint32_t most_significant_bit(uint32_t n) {
    //
    // Returns the most significant bit of an integer in (more or less)
    // constant time. Bit 0 is the least significant bit, bit 31 is the
    // most significant.
    //
    // From public domain code from the Bit Twiddling Hacks web page:
    // http://graphics.stanford.edu/~seander/bithacks.html
    //
    uint32_t tmp1;
    uint32_t tmp2;

    tmp2 = (n >> 16);

    if (tmp2 != 0) {
        tmp1 = (n >> 24);
        if (tmp1 != 0) {
            return (24 + log2table[tmp1]);
        } else {
            return (16 + log2table[tmp2 & 0xFF]);
        }
    } else {
        tmp1 = (n >> 8);
        if (tmp1 != 0) {
            return (8 + log2table[tmp1]);
        } else {
            return log2table[n];
        }
    }
}
//-----------------------------------------------------------------------------
inline uint32_t floor_log2(uint32_t n) {
    return most_significant_bit(n);
}
//-----------------------------------------------------------------------------
inline uint32_t ceil_log2(uint32_t n) {
    uint32_t e = most_significant_bit(n);
    if (! IS_POWER_OF_2_UI(n)) {
        e++;
    }
    return e;
}
//-----------------------------------------------------------------------------
void augment_coprime_base(mpz_t f, mpz_array_t* base) {
    //
    // Try to augment the coprime base with 'f':
    //
    //     - If 'f' is coprime with all other integers in the 'base' array:
    //          - add 'f' in the base
    //
    //     - Else if y is the first element in the base for which
    //       gcd(y, 'f') != 1:
    //          - If gcd(y, 'f') == 'f':
    //                - add 'f' in the base
    //                - remove y from the base
    //                - call augment_coprime_base for y/gcd(y, 'f')
    //          - If gcd(y, 'f') == y:
    //                - keep y in the array
    //                - call augment_coprime_base for 'f'/gcd(y, 'f')
    //          - Otherwise:
    //                - add gcd in the array
    //                - call augment_coprime_base for f/gcd(y, 'f')
    //                - call augment_coprime_base for y/gcd(y, 'f')
    //
    if (is_in_mpz_array(f, base)) {
        return;
    }
    if (mpz_cmp_ui(f, 1) == 0) {
        return;
    }

    uint32_t size_f = mpz_sizeinbase(f, 2);
    uint32_t len    = base->length;

    mpz_t gcd;
    mpz_init2(gcd, size_f);

    bool coprime_with_all_others = true;

    //
    // Add the integer f in the base if it is coprime with all other
    // integers in the base.
    //
    for (uint32_t i = 0; i < len; i++) {

        mpz_gcd(gcd, base->data[i], f);

        if ((mpz_cmp_ui(gcd, 1) != 0) ) {
            coprime_with_all_others = false;
            break;
        }
    }
    if (coprime_with_all_others) {
        if (! is_in_mpz_array(f, base)) {
            append_mpz_to_array(base, f);
        }
        goto clear_gcd_and_return;
    }

    mpz_t cofactor_1;
    mpz_t cofactor_2;

    mpz_init2(cofactor_1, size_f);
    mpz_init2(cofactor_2, size_f);

    for (uint32_t i = 0; i < len; i++) {

        mpz_gcd(gcd, base->data[i], f);

        if ((mpz_cmp_ui(gcd, 1) == 0) ) {
            continue;
        }

        if ((mpz_cmp(gcd, f) == 0) ) {
            //
            // Keep the gcd in the base and call recursively
            // augment_coprime_base with base->data[i]'s cofactor.
            //
            mpz_divexact(cofactor_2, base->data[i], gcd);
            mpz_set(base->data[i], gcd);
            augment_coprime_base(cofactor_2, base);

            break;
        }
        if ((mpz_cmp(gcd, base->data[i]) == 0) ) {
            //
            // Keep the gcd in the base and call recursively
            // augment_coprime_base with f's cofactor.
            //
            mpz_divexact(cofactor_2, f, gcd);
            augment_coprime_base(cofactor_2, base);

            break;
        }
        //
        // Here, we have the following inequalities:
        //     1) gcd != 1
        //     2) gcd != f
        //     3) gcd != base->data[i]
        //
        // Now, keep gcd in the base array and call recursively
        // augment_coprime_base for the found cofactors
        //
        mpz_divexact(cofactor_1, f, gcd);
        mpz_divexact(cofactor_2, base->data[i], gcd);

        mpz_set(base->data[i], gcd);

        augment_coprime_base(cofactor_1, base);
        augment_coprime_base(cofactor_2, base);

        break;
    }

    mpz_clear(cofactor_1);
    mpz_clear(cofactor_2);

  clear_gcd_and_return:

    mpz_clear(gcd);
}
//-----------------------------------------------------------------------------
void find_coprime_base(mpz_array_t* const base, const mpz_t n,
                       const mpz_array_t* const factors) {

    mpz_t cofactor;
    mpz_init2(cofactor, mpz_sizeinbase(n, 2));

    uint32_t len = factors->length;

    for (uint32_t i = 0 ; i < len; i++) {
        augment_coprime_base(factors->data[i], base);
        mpz_divexact(cofactor, n, factors->data[i]);
        augment_coprime_base(cofactor, base);
    }
    mpz_clear(cofactor);
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                      Number theoretical functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
const int8_t kronecker_data[8] = {0, 1, 0, -1, 0, -1, 0, 1};
//-----------------------------------------------------------------------------
int8_t  kronecker_ui(uint32_t a, uint32_t b) {
    //
    // Algorithm 1.4.10 from the book "A Course in Computational Algebraic
    // Number Theory" by Henri Cohen, Springer-Verlag 1993.
    //
    if (b == 0) {
        if (a != 1) {
            return 0;
        } else {
            return 1;
        }
    }
    if (ARE_EVEN(a, b)) {
        return 0;
    }
    uint32_t v = 0;
    while (IS_EVEN(b)) {
        v++;
        b = b >> 1;
    }
    int32_t k = 1;
    if (! IS_EVEN(v)) {
        k = kronecker_data[a & 7];      // k = (-1)^((a^2-1)/8)
    }
    uint32_t r = 0;
    while (1) {
        if (a == 0) {
            if (b == 1) {
                return k;
            } else {
                return 0;
            }
        }
        v = 0;
        while (IS_EVEN(a)) {
            v++;
            a = a >> 1;
        }
        if (! IS_EVEN(v)) {
            k *= kronecker_data[b & 7]; // k = (-1)^((b^2-1)/8) * k
        }
        if (a & b & 2) {                // k = (-1)^((a-1)(b-1)/4) * k
            k = -k;
        }
        r = a;
        a = b % r;
        b = r;
    }
}
//-----------------------------------------------------------------------------
uint32_t powm(uint32_t base, uint32_t power, uint32_t modulus) {
    //
    // Left-right binary algorithm to compute (base^power) mod modulus.
    //
    // See for example algorithm 1.2.3 from the book "A Course in Computational
    // Algebraic Number Theory" by Henri Cohen, Springer-Verlag 1993.
    //
    if (power == 0) {
        return 1;
    }
    uint32_t n = power;
    uint32_t z = base % modulus;

    uint64_t y = z;
    uint32_t f = most_significant_bit(power);

    while (f != 0) {
        f--;
        y *= y;
        y %= modulus;
        if (BIT(n, f) == 1) {
            y *= z;
            y %= modulus;
        }
    }
    return (uint32_t)y;
}
//-----------------------------------------------------------------------------
uint32_t sqrtm(uint32_t a, uint32_t p) {
    //
    // Shanks' algorithm for modular square roots.
    //
    // See for example algorithm 1.5.1 from the book "A Course in Computational
    // Algebraic Number Theory" by Henri Cohen, Springer-Verlag 1993.
    //
    if (a == 0) {
        //
        // a does not have a modular square root!
        //
        return NO_SQRT_MOD_P;
    }
    //
    // Find n such that p = q.2^n + 1 with q odd
    //
    uint32_t n = ffs(p - 1) - 1;
    uint32_t q = (p - 1) >> n;
    uint32_t k = n;
    //
    // Find u, a quadratic non-residue mod p.
    //
    // _WARNING_: Here p should be an odd prime!
    //
    // _NOTE_: In his book "A Course in Computational Algebraic Number Theory",
    //         Henri Cohen suggests to choose u at random until a non-residue
    //         is found instead of a sequential search. In practice however, it
    //         should not make any significant difference so we follow here the
    //         path of least resistance.
    //
    uint32_t u = 1;
    while (kronecker_ui(u, p) != -1) {
        u++;
    }
    uint32_t z   = powm(u, q, p);
    uint32_t x   = powm(a, (q + 1) / 2, p);
    uint32_t b   = powm(a, q, p);
    uint32_t m   = 0;
    uint32_t t   = 0;
    uint32_t b2m = 0;

    while (1) {
        b2m = b;
        m   = 0;
        //
        // Find the least integer m such that b^(2^m) = 1 (mod p).
        //
        while (b2m != 1) {
            m++;
            b2m = (uint32_t)( ((uint64_t)b2m * (uint64_t)b2m) % p);
        }
        if (m == k) {
            //
            // a does not have a modular square root!
            //
            return NO_SQRT_MOD_P;
        }
        t = powm(z, 1 << (k - m - 1), p);
        z = (uint32_t)( ((uint64_t)t * (uint64_t)t) % p);
        b = (uint32_t)( ((uint64_t)b * (uint64_t)z) % p);
        x = (uint32_t)( ((uint64_t)x * (uint64_t)t) % p);

        if (b == 1) {
            return x;
        }
        k = m;
    }
}
//-----------------------------------------------------------------------------
static const unsigned short qres_mod_315[315] = {
    1, 1, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 1, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0
};
static const unsigned short qres_mod_256[256] = {
    1, 1, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 1, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0
};
static const unsigned short qres_mod_221[221] = {
    1, 1, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 1, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
    1, 0, 0, 0, 0, 1, 1, 0, 1, 0,
    0, 0, 1, 1, 0, 0, 0, 0, 0, 1,
    0, 1, 1, 1, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 1, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    1, 1, 0, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
    1, 1, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 0, 1, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 1, 1,
    1, 0, 1, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 1, 0, 1, 1, 0, 0, 0,
    0, 1, 0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 1, 0,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
    1
};
//-----------------------------------------------------------------------------
unsigned long int is_square(unsigned long int x) {
    //
    // Basic perfect square detection test.
    //
    // See for example algorithm 1.7.3 from the book "A Course in Computational
    // Algebraic Number Theory" by Henri Cohen, Springer-Verlag 1993.
    //
    // The description given in this book has been adapted to use larger
    // tables in exchange of a slight performance boost (about 30% on
    // Opteron 250).
    //
    if (qres_mod_256[x & 255] == 0) {
        //
        // Get rid of about 82.8% of non squares
        //
        return 0;
    }
    if (qres_mod_315[x % 315] == 0) {
        //
        // Get rid of about 84.8% of non squares
        //
        return 0;
    }
    if (qres_mod_221[x % 221] == 0) {
        //
        // Get rid of about 71.5% of non squares
        //
        return 0;
    }
    unsigned long int root = (unsigned long int)sqrt(x);
    if ((root * root) == x) {
        return root;
    }
    return 0;
}
//-----------------------------------------------------------------------------
unsigned long int gcd_ulint(unsigned long int a, unsigned long int b) {
    //
    // Standard "right-shift" binary gcd algorithm.
    //
    // See for example algorithm 1.3.5 from the book "A Course in Computational
    // Algebraic Number Theory" by Henri Cohen, Springer-Verlag 1993.
    //
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    unsigned long int tmp = 0;
    if (a < b) {
        tmp = a;
        a   = b;
        b   = tmp;
    }
    tmp = a % b;
    a   = b;
    b   = tmp;

    if (b == 0) {
        return a;
    }
    unsigned long int shift = 0;

    while (ARE_EVEN(a, b)) {
        shift++;
        a >>= 1;
        b >>= 1;
    }
    while (IS_EVEN(a)) {
        a >>= 1;
    }
    while (IS_EVEN(b)) {
        b >>= 1;
    }
    unsigned long int t = 0;

    do {
        if (a > b) {
            t   = a;
            t  -= b;
            t >>= 1;
            if (t == 0) {
                return (a << shift);
            }
            while (IS_EVEN(t)) {
                t >>= 1;
            }
            a = t;
        } else {
            t   = b;
            t  -= a;
            t >>= 1;
            if (t == 0) {
                return (a << shift);
            }
            while (IS_EVEN(t)) {
                t >>= 1;
            }
            b = t;
        }
    } while (1);

    return (a << shift);
}
//-----------------------------------------------------------------------------
//
// The following macros are used by the modinv_ui function, where the
// *_bit_mask variables are declared, and are thus not intended to be reused
// elsewhere.
//
#define BIT_N(x)           ( (x) & (nth_bit_mask))
#define BIT_N_MINUS_ONE(x) ( (x) & (nth_minus_one_bit_mask) )
#define HIGH_BITS(x)       ( (x) & (high_bit_mask) )
#define LOW_BITS(x)        ( (x) & (low_bit_mask) )
//-----------------------------------------------------------------------------
unsigned long int modinv_ui(unsigned long int num, unsigned long int p) {
    //
    // Reference:
    //      "New Algorithm for Classical Modular Inverse",
    //      Robert Lorencz,
    //      Lecture Notes In Computer Science, Volume 2523/2003,
    //      Revised Papers from the 4th International Workshop on
    //      Cryptographic Hardware and Embedded Systems
    //
    unsigned long int a = num % p;
    long int u = p;
    long int v = a;
    long int r = 0;
    long int s = 1;
    unsigned long int cu = 0;
    unsigned long int cv = 0;
    long int exp_cu = 1;
    long int exp_cv = 1;
    unsigned long int n  = most_significant_bit(p) + 1;
    //
    // If n is the number of bits of p, computes masks to be used by the
    // previously defined macros, to:
    //   - keep the nth bit (i.e. the most significant one);
    //   - keep the nth-1 bit (i.e. the next most significant one);
    //   - keep the 2 most significant bits;
    //   - keep the nth-2 least significant bits.
    //
    const unsigned long int high_bit_mask = 3 << (n-1);
    const unsigned long int low_bit_mask  = (1 << (n-1)) - 1;
    const unsigned long int nth_bit_mask  = 1 << n;
    const unsigned long int nth_minus_one_bit_mask  = 1 << (n-1);

    while ((u != exp_cu) && (u != -exp_cu) && (v != exp_cv) && (v != -exp_cv)) {

        if (!HIGH_BITS(u) || (BIT_N(u) && BIT_N_MINUS_ONE(u) && LOW_BITS(u))) {
            if (cu >= cv) {
                u <<= 1;
                r <<= 1;
                cu++;
                exp_cu <<= 1;
            } else {
                u <<= 1;
                s >>= 1;
                cu++;
                exp_cu <<= 1;
            }
        } else {
            if (   !HIGH_BITS(v)
                || (BIT_N(v) && BIT_N_MINUS_ONE(v) && LOW_BITS(v)) ) {

               if (cv >= cu) {
                   v <<= 1;
                   s <<= 1;
                   cv++;
                   exp_cv <<= 1;
               } else {
                   v <<= 1;
                   r >>= 1;
                   cv++;
                   exp_cv <<= 1;
               }
            } else {
                if (BIT_N(v) == BIT_N(u)) {
                    if (cu <= cv) {
                        u -= v;
                        r -= s;
                    } else {
                        v -= u;
                        s -= r;
                    }
                } else {
                    if (cu <= cv) {
                        u += v;
                        r += s;
                    } else {
                        v += u;
                        s += r;
                    }
                }
            }
        }
    }
    if ((v == exp_cv) || (v == -exp_cv)) {
        r = s;
        u = v;
    }
    if (BIT_N(u) != 0) {
        if (r < 0) {
            r = -r;
        } else {
            r = p - r;
            if (r < 0) {
                r += p;
            }
        }
    } else {
        if (r < 0) {
            r += p;
        }
    }
    return r;
}
//-----------------------------------------------------------------------------
unsigned long int sqrtm_p2(uint32_t n, uint32_t p) {
    //
    // _NOTE_: Now this is getting a bit messy. Indeed, the code now mix
    //         apparently randomly (unsigned) long ints with C99 integer types
    //         like uint32_t. This can be seen as bad pratice and indeed, it is
    //         not very pretty. The reason for this is to convey an idea of
    //         the size of the numbers involved. Maybe only using long ints and
    //         clearly documenting the expected range of the variables would
    //         have been better.
    //
    long int s = 0;

    uint32_t np = n % p;
    uint32_t sp = sqrtm(np, p);

    if (sp == NO_SQRT_MOD_P) {
        //
        // No solution
        //
        return NO_SQRT_MOD_P2;
    }

    s = sp * sp;
    s = n - s;

    uint32_t inv = modinv_ui(sp << 1, p);

    s /= (long int)p;
    s *= inv;

    if (s < 0) {
        s = p - ((-s) % p);
    } else {
        s = s % p;
    }

    s *= p;
    s += sp;

    return (unsigned long int)s;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                          Hash functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
uint32_t hash_rj_32(const void* const keyptr) {

    uint32_t hash = *((uint32_t*)keyptr);
    //
    // Robert Jenkins' 32 bit mix function.
    // http://www.concentric.net/~Ttwang/tech/inthash.htm
    //
    hash += (hash << 12);
    hash ^= (hash >> 22);
    hash += (hash << 4);
    hash ^= (hash >> 9);
    hash += (hash << 10);
    hash ^= (hash >> 2);
    hash += (hash << 7);
    hash ^= (hash >> 12);
    return hash;
}
//-----------------------------------------------------------------------------
uint32_t hash_pjw(const void* const keyptr) {
    //
    // Found in the Dragon book:
    // "Compilers: Principles, Techniques and Tools", Aho, Sethi, & Ullman.
    //
    // Attributed to P. J. Weinberger.
    // Supposed to perform well on strings.
    //
    char *p;
    uint32_t hash = 0;
    uint32_t g = 0;

    for (p = (char*)keyptr; *p != '\0'; p++) {
        hash = (hash<<4) + (*p);
        g = (hash & 0xf0000000);
        if (0 != g) {
            hash ^= (g>>24);
            hash ^= g;
        }
    }
    return hash;
}
//-----------------------------------------------------------------------------
#define GET_16_BITS(d) (*((const uint16_t *) (d)))
//-----------------------------------------------------------------------------
uint32_t hash_sfh_ph(const void* const keyptr) {
    //
    // The so-called "SuperFastHash" function By Paul Hsieh.
    // http://www.azillionmonkeys.com/qed/hash.html
    //
    // The original code has been slightly altered, but the logic is
    // unchanged. Refer to the aforementioned URL for P. Hsieh's original code.
    //
    if ((keyptr == NULL)) {
        return 0;
    }
    const char* data = (const char*)keyptr;

    uint32_t length = strlen(data)*sizeof(char);

    uint32_t hash = length;
    uint32_t tmp;

    uint32_t rem = length & 3;
    length >>= 2;

    while (length > 0) {
        //
        // Use the next 32 bits of *data
        //
        length--;
        hash  += GET_16_BITS(data);
        tmp    = (GET_16_BITS(data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2 * sizeof(uint16_t);
        hash  += hash >> 11;
    }
    switch (rem) {
        //
        // If the number of bits of *data is not a multiple of 32
        // i.e. (length % 4) != 0
        //
        case 3:
            hash += GET_16_BITS(data);
            hash ^= hash << 16;
            hash ^= data[sizeof(uint16_t)] << 18;
            hash += hash >> 11;
            break;
        case 2:
            hash += GET_16_BITS(data);
            hash ^= hash << 11;
            hash += hash >> 17;
            break;
        case 1:
            hash += *data;
            hash ^= hash << 10;
            hash += hash >> 1;
    }
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 2;
    hash += hash >> 15;
    hash ^= hash << 10;

    return hash;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                          Comparison functions
 *-----------------------------------------------------------------------------
 */

//
// The following functions return an int instead of an int32_t only
// for compatibility with pre-C99 code such as the standard library
// qsort function...
//

//-----------------------------------------------------------------------------
int mpz_cmp_func(const void* const mpza, const void* const mpzb) {

    int cmp = mpz_cmp(*((mpz_t*)mpza), *((mpz_t*)mpzb));
    if (cmp == 0) return 0;
    if (cmp > 0) {
        return 1;
    } else {
        return -1;
    }
}
//-----------------------------------------------------------------------------
int uint32_cmp_func(const void* const uinta, const void* const uintb) {

    if (*((uint32_t*)uinta) == *((uint32_t*)uintb)) return 0;
    if (*((uint32_t*)uinta) > *((uint32_t*)uintb)) {
        return 1;
    } else {
        return -1;
    }
}
//-----------------------------------------------------------------------------
int string_cmp_func(const void* const stra, const void* const strb) {

    int cmp = strcmp((char*)stra, (char*)strb);
    if (cmp == 0) return 0;
    if (cmp > 0) {
        return 1;
    } else {
        return -1;
    }
}
//-----------------------------------------------------------------------------
int cmp_mult_data(const void* mda, const void* mdb) {

    const mult_data_t* const a = (mult_data_t*)mda;
    const mult_data_t* const b = (mult_data_t*)mdb;

    if (a->count > b->count) {
        return 1;
    }
    if (a->count < b->count) {
        return -1;
    }
    //
    // Here, a->count == b->count
    //
    if (a->sum_inv_pi > b->sum_inv_pi) {
        return 1;
    }
    if (a->sum_inv_pi < b->sum_inv_pi) {
        return -1;
    }
    //
    // Here, a->count      == b->count
    // and   a->multiplier == b->multiplier
    //
    if (a->multiplier < b->multiplier) {
        //
        // a's multiplier is smaller, which actually means that
        // it is a better candidate multiplier than b's multiplier.
        // Hence return 1 (and not -1).
        //
        return 1;
    }
    if (a->multiplier > b->multiplier) {
        //
        // a's multiplier is larger, which actually means that
        // it is a worse candidate multiplier than b's multiplier.
        // Hence return -1 (and not 1).
        //
        return -1;
    }
    return 0;
}
//-----------------------------------------------------------------------------

