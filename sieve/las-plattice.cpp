#include "cado.h"

#include "fb.h"                                                                 
#include "utils.h"           /* lots of stuff */                                
#include "las-types.h"
#include "las-plattice.h"

uint32_t plattice_enumerate_t::maskI;
plattice_x_t plattice_enumerate_t::even_mask;
//plattice_x_t plattice_enumerate_t::area;

template <>
plattice_x_t
plattice_enumerate_area<1>::value = 0;

template <>
plattice_x_t
plattice_enumerate_area<2>::value = 0;

template <>
plattice_x_t
plattice_enumerate_area<3>::value = 0;
