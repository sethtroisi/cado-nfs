#ifndef MATRIX_REPR_PRIME_HPP_
#define MATRIX_REPR_PRIME_HPP_

#include <iterator>
#include <istream>

#include "manu.h"
#include "matrix_line.hpp"
#include "matrix_repr.hpp"

struct matrix_repr_prime {
    struct ptr {
	uint32_t * idx;
	int32_t * val;
    };
    struct const_ptr {
        const uint32_t * idx;
        const int32_t * val;
        const_ptr(ptr const& p) : idx(p.idx), val(p.val) {}
    };

    struct matrix_rowset : public ptr {
	uint32_t nr;
	uint32_t nc;
	matrix_rowset() {
	    idx = NULL;
	    val = NULL;
	}
	/* This thing does not have a copy ctor ! */

	void alloc(uint r, uint c) {
	    nr = r;
	    nc = c;
	    BUG_ON(idx != NULL); idx = new uint32_t[nr + nc];
	    BUG_ON(val != NULL); val = new  int32_t[nr + nc];
	}
	~matrix_rowset() {
	    if (idx != NULL) delete[] idx;
	    if (val != NULL) delete[] val;
	}

	void fill(std::istream& mtx, std::streampos pos,
		uint i0, uint i1)
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
		for(lit_t lit = l.begin() ; lit != l.end() ; lit++) {
		    *q.idx++ = lit->first - v;
		    *q.val++ = lit->second;
		    v = lit->first;
		}
		*q.idx++ = 0; BUG_ON((uint) (q.idx - idx) > (nr + nc));
		*q.val++ = 0; BUG_ON((uint) (q.val - val) > (nr + nc));
	    }
	    BUG_ON(i != i1);
	}
	template<typename traits>
	void mul(
		typename traits::wide_scalar_t * dst,
		const typename traits::scalar_t * src) const
	{
	    int acc=0;
	    const_ptr q(*this);
	    for(uint i = 0 ; i < nr ; i++) {
		traits::zero(dst[i]);
		unsigned int c = 0;
		for( ; *q.val != 0 ; q.idx++, q.val++) {
		    c += *q.idx;
		    traits::addmul(dst[i], src[c], *q.val);
		    if (++acc == traits::max_accumulate) {
			traits::reduce(dst[i], dst[i]);
			acc = 1;
		    }
		}
		q.idx++, q.val++;
	    }
	}
    };	/* matrix_rowset */
};

/* vim: set sw=4: */

#endif	/* MATRIX_REPR_PRIME_HPP_ */
