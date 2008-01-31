#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_

#include <gmp.h>
#include <gmpxx.h>
#include "config.hpp"

namespace globals {
    extern unsigned int m, n;
    extern unsigned int nr, nc;
    extern unsigned int nbys;

    /* This is a placeholder in most cases */
    extern field K;

    extern mpz_class modulus;

    /* These will probably go sooner or later. */
    extern uint8_t modulus_u8;
    extern uint16_t modulus_u16;
    extern uint32_t modulus_u32;
    extern unsigned long modulus_ulong;

    void set_modulus(mpz_class const& modulus);
}

#endif	/* GLOBALS_HPP_ */
