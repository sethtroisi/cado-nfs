#ifndef TRAITS_HPP_
#define TRAITS_HPP_

#include "typical_scalar_traits.hpp"
#include "variable_scalar_traits.hpp"
#include "sse2_8words_traits.hpp"
#include "sse2_2words_traits.hpp"
#include "ulong_traits.hpp"
#include "binary_sse2_traits.hpp"
#include "binary_ulong_traits.hpp"

#define TRY_ONE_BINARY_CODE(tmpl,impl)					\
	do {								\
		std::cerr << impl::name()();				\
		if (impl::can()) {					\
			std::cerr << " --> ok, go\n";			\
			return tmpl<impl>();				\
		} else {						\
			std::cerr << " --> does not apply\n";		\
		}							\
	} while (0)

/* This list defines the different code variants that are included in the
 * binary. For convenience, everything you might think of is allowed to
 * be included here, although obviously reducing the list does save
 * trees.
 */
#define	TRY_ALL_BINARY_CODES(tmpl) do {				\
	TRY_ONE_BINARY_CODE(tmpl,binary_ulong_traits);		\
	TRY_ONE_BINARY_CODE(tmpl,binary_sse2_traits);		\
	TRY_ONE_BINARY_CODE(tmpl,sse2_8words_traits);		\
	TRY_ONE_BINARY_CODE(tmpl,sse2_2words_traits);		\
	TRY_ONE_BINARY_CODE(tmpl,ulong_traits);			\
	TRY_ONE_BINARY_CODE(tmpl,typical_scalar_traits<1>);	\
	TRY_ONE_BINARY_CODE(tmpl,typical_scalar_traits<2>);	\
	TRY_ONE_BINARY_CODE(tmpl,typical_scalar_traits<3>);	\
	TRY_ONE_BINARY_CODE(tmpl,typical_scalar_traits<4>);	\
	TRY_ONE_BINARY_CODE(tmpl,typical_scalar_traits<5>);	\
	TRY_ONE_BINARY_CODE(tmpl,typical_scalar_traits<6>);	\
	TRY_ONE_BINARY_CODE(tmpl,typical_scalar_traits<7>);	\
	TRY_ONE_BINARY_CODE(tmpl,typical_scalar_traits<8>);	\
	TRY_ONE_BINARY_CODE(tmpl,typical_scalar_traits<9>);	\
	TRY_ONE_BINARY_CODE(tmpl,typical_scalar_traits<10>);	\
	TRY_ONE_BINARY_CODE(tmpl,variable_scalar_traits);	\
} while (0)


#endif	/* TRAITS_HPP_ */
