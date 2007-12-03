#ifndef STATE_HPP_
#define STATE_HPP_

#include "constants.hpp"
#include <gmp.h>
#include <ostream>
#include <fstream>

extern int recoverable_iteration(int, int, int, int);
extern int recover_iteration(int, int, int, int);

template<typename traits>
int recover_vector(int nr, int col, int nbys, int r, typename traits::scalar_t * w);

#include "state.tcc"

#endif	/* STATE_HPP_ */
