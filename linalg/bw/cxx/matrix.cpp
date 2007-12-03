#include "matrix.hpp"
#include "manu.h"
#include "matrix_line.hpp"
#include "parsing_tools.hpp"

#include <iterator>
#include <string>
#include <sstream>

/* The routines here populate an in-memory structure with the
 * speed-critical version of the matrix data. Everything is permitted
 * as for HOW this data is represented. */

/* the fill() routine is not speed-critical by itself. */

/* typical layout: */
#if 0
struct matrix_repr_XXXXX {
    struct ptr { XXXXX };
    static void fill(std::istream& mtx, std::streampos pos, uint32_t nc,
	    uint i0, uint i1, ptr p);
    static void count(std::istream& mtx,
	    unsigned int nr,
	    vector<std::streampos>& offsets,
	    vector<boost::uint32_t>& ncoeffs);
};
#endif

using namespace std;

void count_matrix_coeffs(istream& mtx,
                unsigned int nr,
                vector<streampos>& offsets,
                vector<boost::uint32_t>& ncoeffs)
{
    BUG_ON(offsets.empty());
    BUG_ON(offsets.size() != ncoeffs.size());
    unsigned int nt = offsets.size();

    offsets.clear();
    ncoeffs.clear();

    /* Note that the std::istream_iterator ctor advances beforehand in the
     * file. So we keep track of the file position in advance. */
    /* TODO: now that the iterator is gone here, perhaps it's
     * sensible to simplify this function. */

    std::streampos pos = mtx.tellg();

    comment_strip cs(mtx, "//");

    // std::istream_iterator<matrix_line> mtxi(mtx);

    for(uint j = 0 ; j < nt ; j++) {
	uint i0, i1;
	boost::uint32_t nc = 0;

	offsets.push_back(pos);

	i0 = (j * nr) / nt;
	i1 = ((j + 1) * nr) / nt;
	for(uint i = i0 ; i < i1 ; i++) {
	    pos = mtx.tellg();
	    std::string s;
	    cs.getline(s);
	    std::istringstream st(s);
	    uint z;
	    if (!(st >> z)) {
		BUG();
	    }
	    nc += z;
	}
	ncoeffs.push_back(nc);
    }
}


/* vim: set sw=4: */
