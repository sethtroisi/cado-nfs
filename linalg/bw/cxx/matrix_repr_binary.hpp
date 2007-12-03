#ifndef MATRIX_REPR_BINARY_HPP_
#define MATRIX_REPR_BINARY_HPP_

#include <iterator>
#include <istream>

#include "matrix_line.hpp"
#include "matrix_repr.hpp"

struct matrix_repr_binary {
    struct ptr { uint32_t * idx; };
    struct const_ptr {
	const uint32_t * idx;
	const_ptr(ptr const& p) :idx(p.idx) {}
    };

    struct matrix_rowset : public ptr {
	uint32_t nr;
	uint32_t nc;
	matrix_rowset() { idx = NULL; }
	/* This thing does not have a copy ctor ! */

	void alloc(uint r, uint c) {
	    nr = r;
	    nc = c;
	    BUG_ON(idx != NULL); idx = new uint32_t[nr + nc];
	}
	~matrix_rowset() {
	    if (idx != NULL) delete[] idx;
	}

	void fill(std::istream& mtx, std::streampos pos, uint i0, uint i1)
	{
	    mtx.seekg(pos);
	    std::istream_iterator<matrix_line> mit(mtx);
	    uint i;
	    ptr q((ptr const&) *this);
	    BUG_ON(i1 - i0 != nr);
	    for(i = i0 ; i < i1 && mit != endof<matrix_line>(mtx); i++) {
		matrix_line l = *mit++;
		int v = 0;
		typedef matrix_line::const_iterator lit_t;
		*q.idx++ = l.size();
		for(lit_t lit = l.begin() ; lit != l.end() ; lit++) {
		    *q.idx++ = lit->first - v;
		    v = lit->first;
		}
		BUG_ON((uint) (q.idx - idx) > (nr + nc));
	    }
	    BUG_ON(i != i1);
	}

	template<typename traits>
	void mul(
		typename traits::wide_scalar_t * dst,
		const typename traits::scalar_t * src) const
	{
	    const_ptr q(*this);
	    for(uint i = 0 ; i < nr ; i++) {
		traits::zero(dst[i]);
		unsigned int nc = *q.idx++;
		unsigned int c = 0;
		for( ; nc-- ; q.idx++) {
		    c += *q.idx;
		    /* if we're here, then the traits:: class
		     * should have a bare addmul() routine that
		     * does NOT take a scalar as third int32_t
		     * argument.
		     */
		    traits::addmul(dst[i], src[c]);
		}
	    }
	}
    };	/* matrix_rowset */
};

/* vim: set sw=4: */

#endif	/* MATRIX_REPR_BINARY_HPP_ */
