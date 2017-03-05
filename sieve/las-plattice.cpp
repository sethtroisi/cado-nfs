#include "cado.h"

#include "fb.hpp"
#include "utils.h"           /* lots of stuff */
#include "las-types.hpp"
#include "las-plattice.hpp"

uint32_t plattice_enumerate_t::maskI;
plattice_x_t plattice_enumerate_t::even_mask;

template <>
plattice_x_t
plattice_enumerate_area<1>::value = 0;

template <>
plattice_x_t
plattice_enumerate_area<2>::value = 0;

template <>
plattice_x_t
plattice_enumerate_area<3>::value = 0;


void plattice_enumerate_t::advance_to_next_area(int level) {
    switch (level) {
        case 1:
            x -= plattice_enumerate_area<1>::value; break;
        case 2:
            x -= plattice_enumerate_area<2>::value; break;
        case 3:
            x -= plattice_enumerate_area<3>::value; break;
        default:
            ASSERT_ALWAYS(0);
    }
}
