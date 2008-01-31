#ifndef TRAITS_GLOBALS_HPP_
#define TRAITS_GLOBALS_HPP_

#include <boost/cstdint.hpp>
#include <gmp.h>
#include <gmpxx.h>
#include <istream>
#include <ostream>
#include "gmp-hacks.h"
#include "basefield.hpp"
#include "globals.hpp"

#define	MODBITS()	mpz_sizeinbase(globals::modulus.get_mpz_t(),2)
#define	MODLIMBS()	SIZ(globals::modulus.get_mpz_t())

#endif	/* TRAITS_GLOBALS_HPP_ */
