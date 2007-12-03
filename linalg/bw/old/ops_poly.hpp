#ifndef OPS_POLY_HPP_
#define OPS_POLY_HPP_

#include "ops_poly_fft.hpp"
#include "ops_poly_schoolbook.hpp"
#include "ops_poly_integerfft.hpp"
#include "ops_poly_wrap.hpp"

/*
template<typename T>
struct ops_poly : public T
{ };
*/

// typedef ops_poly_fft ft_order_t;
// typedef ops_poly_scbk ft_order_t;
#undef	HAS_CONVOLUTION_SPECIAL
typedef ops_poly_ifft ft_order_t;
// typedef ops_poly_wrap ft_order_t;

#endif	/* OPS_POLY_HPP_ */
