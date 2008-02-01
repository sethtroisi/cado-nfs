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
	TRY_ONE_BINARY_CODE(tmpl,binary_sse2_traits);		\
	TRY_ONE_BINARY_CODE(tmpl,binary_uint8_traits);		\
	TRY_ONE_BINARY_CODE(tmpl,binary_uint16_traits);		\
	TRY_ONE_BINARY_CODE(tmpl,binary_uint32_traits);		\
	TRY_ONE_BINARY_CODE(tmpl,binary_uint64_traits);         \
} while (0)

/* removed lots of traits types used by the prime field version */



#endif	/* TRAITS_HPP_ */
