#ifndef STATE_HPP_
#define STATE_HPP_

#include "constants.hpp"
#include <gmp.h>
#include <ostream>
#include <fstream>

extern int recoverable_iteration(int, int, int);
extern int recover_vector(int, int, int, mp_limb_t (*)[MODULUS_SIZE]);
extern int recover_iteration(int, int, int);

#endif	/* STATE_HPP_ */
