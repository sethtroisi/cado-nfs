// #include <cstdio>
// #include <cstdlib>
// #include "auxfuncs.h"
// #include "manu.h"
// 
// #include <sys/types.h>
// #include <unistd.h>
#include <gmp.h>
#include <gmpxx.h>
#include "traits.hpp"
// #include "gmp-hacks.h"
// 
// #include "arguments.hpp"
// #include "common_arguments.hpp"
// #include "constants.hpp"
// #include "slave_arguments.hpp"
#include "files.hpp"
// #include "matrix_header.hpp"
// #include "matrix_line.hpp"
#include "must_open.hpp"
// #include "parsing_tools.hpp"
#include <iterator>
#include "state.hpp"

namespace state_details {
template<typename T> inline std::istream_iterator<T> beginof(std::ifstream& i)
{ return std::istream_iterator<T>(i); }
template<typename T> inline std::istream_iterator<T> endof(std::ifstream&)
{ return std::istream_iterator<T>(); }
}

template<typename traits>
int recover_vector(int nr, int col, int r, typename traits::scalar_t * w)
{
	std::ifstream v;
	using namespace state_details;
	if (r) {
		must_open(v, files::v % col % r);
	} else {
		must_open(v, files::y % col);
	}
	std::istream_iterator<mpz_class> it(v);
	for( ; nr && it != endof<mpz_class>(v) ; it++, nr--, w++) {
		traits::assign(*w, *it);
	}
	BUG_ON(nr != 0);
	return 0;
}
