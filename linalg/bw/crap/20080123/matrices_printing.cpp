#include "matrices_printing.hpp"

namespace hidden {
	int const magma_index = std::ios_base::xalloc();
}

namespace std {
	std::ios_base& set_magma(std::ios_base& o) {
		o.iword(::hidden::magma_index) = 1;
		return o;
	}
	std::ios_base& unset_magma(std::ios_base& o) {
		o.iword(::hidden::magma_index) = 0;
		return o;
	}
	bool test_magma(std::ios_base& o) {
		return o.iword(::hidden::magma_index);
	}
}
