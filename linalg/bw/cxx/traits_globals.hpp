#ifndef TRAITS_GLOBALS_HPP_
#define TRAITS_GLOBALS_HPP_

#include <boost/cstdint.hpp>
#include <gmp.h>
#include <gmpxx.h>
#include <istream>
#include <ostream>
#include "gmp-hacks.h"

namespace globals {
	extern uint8_t modulus_u8;
	extern uint16_t modulus_u16;
	extern uint32_t modulus_u32;
	extern mpz_class modulus;
	extern unsigned long modulus_ulong;
	extern unsigned int nbys;
}

#include "coeff_field_tags.hpp"

#define	MODBITS()	mpz_sizeinbase(globals::modulus.get_mpz_t(),2)
#define	MODLIMBS()	SIZ(globals::modulus.get_mpz_t())

#endif	/* TRAITS_GLOBALS_HPP_ */
