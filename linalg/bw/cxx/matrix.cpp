#include "matrix.hpp"
#include "manu.h"
#include "matrix_line.hpp"

#include <iterator>

using namespace std;

template<typename T> inline istream_iterator<T> beginof(istream& i)
{ return istream_iterator<T>(i); }
template<typename T> inline istream_iterator<T> endof(istream&)
{ return istream_iterator<T>(); }

void fill_matrix_data(istream& mtx,
		streampos pos, uint32_t nc,
		uint i0, uint i1, 
		uint32_t * idx, int32_t * val)
{
	mtx.seekg(pos);
	istream_iterator<matrix_line> mit(mtx);
	uint i;
	uint32_t * ip = idx;
	int32_t * vp = val;
	for(i = i0 ; i < i1 && mit != endof<matrix_line>(mtx); i++, mit++) {
		matrix_line l;
		l = *mit;
		int v = 0;
		typedef matrix_line::const_iterator lit_t;
		for(lit_t lit = l.begin() ; lit != l.end() ; lit++) {
			*ip++ = lit->first - v;
			*vp++ = lit->second;
			v = lit->first;
		}
		*ip++ = 0;
		*vp++ = 0;
		BUG_ON((uint) (ip - idx) > (i1 - i0 + nc));
		BUG_ON((uint) (vp - val) > (i1 - i0 + nc));
	}
	BUG_ON(i != i1);
}

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

	/* Note that the istream_iterator ctor advances beforehand in the
	 * file. So we keep track of the file position in advance. */

	streampos pos = mtx.tellg();
	istream_iterator<matrix_line> mtxi(mtx);

	for(uint j = 0 ; j < nt ; j++) {
		uint i0, i1;
		boost::uint32_t nc = 0;

		offsets.push_back(pos);

		i0 = (j * nr) / nt;
		i1 = ((j + 1) * nr) / nt;
		for(uint i = i0 ; i < i1 ; i++) {
			pos = mtx.tellg();
			matrix_line z = *mtxi++;
			nc += z.size();
		}
		ncoeffs.push_back(nc);
	}
}
