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
#include <vector>
#include "state.hpp"

namespace state_details {
template<typename T> inline std::istream_iterator<T> beginof(std::ifstream& i)
{ return std::istream_iterator<T>(i); }
template<typename T> inline std::istream_iterator<T> endof(std::ifstream&)
{ return std::istream_iterator<T>(); }
}

#if 0
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
	for( ; nr && it != endof<mpz_class>(v) ; nr--, w++) {
		std::vector<mpz_class> foo(globals::nbys);
		for(int i = 0 ; i < globals::nbys ; i++) {
			BUG_ON(it == endof<mpz_class>(v));
			foo[i] = *it++;
		}
		traits::assign(*w, foo, 0);
	}
	BUG_ON(nr != 0);
	return 0;
}
#endif

template<typename traits>
int recover_vector(int nr, int col, int nbys, int r, typename traits::scalar_t * w)
{
	std::ifstream vs[nbys];

	using namespace state_details;
	for(int i = 0 ; i < nbys ; i++) {
		if (r) {
			must_open(vs[i], files::v % (col+i) % r);
		} else {
			must_open(vs[i], files::y % (col+i));
		}
	}
	for( ; nr ; nr--, w++) {
		std::vector<mpz_class> foo;
		int nfail = 0;
		for(int i = 0 ; i < nbys ; i++) {
			mpz_class x;
			if (vs[i] >> x) {
				foo.push_back(x);
			} else {
				nfail++;
			}
		}
		// nr should be such that everything is ok.
		BUG_ON(nfail);
		traits::assign(*w, foo, 0);
	}
	BUG_ON(nr != 0);
	return 0;
}
