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
// #include <iterator>
#include <iomanip>
#include "state.hpp"

/* 
namespace state_details {
template<typename T> inline std::istream_iterator<T> beginof(std::ifstream& i)
{ return std::istream_iterator<T>(i); }
template<typename T> inline std::istream_iterator<T> endof(std::ifstream&)
{ return std::istream_iterator<T>(); }
}
*/

template<typename traits>
int recover_vector(int nr, int col, int nbys, int r, typename traits::scalar_t * w)
{
        using namespace std;
	std::ifstream vs;

        int nfail = 0;
        int offset;
        int expected;

        if (r) {
            must_open(vs, files::v % col % r);
            offset = 0;
            expected = nbys;
        } else {
            BUG_ON(col % nbys != 0);
            must_open(vs, files::y);
            offset = col;
            expected = globals::n;
        }

        for( ; nr ; nr--, w++) {
            std::string s;
            getline(vs, s);
            std::istringstream str(s);
            int z;
            typename traits::scalar_t tmp;

            for(z = 0 ; z < offset    ; z += nbys) {
                traits::get(str, tmp);
                str >> ws;
            }
            traits::get(str, *w);
            nfail += !str;
            str >> ws;
            for(      ; z < expected ; z += nbys)  {
                traits::get(str, tmp);
                str >> ws;
            }
            nfail += !str.eof();
        }

        BUG_ON(nfail);
	BUG_ON(nr != 0);
	return 0;
}
