#ifndef MUL_HPP_
#define MUL_HPP_

#include <boost/cstdint.hpp>

template<typename traits>
void multiply_ur(typename traits::wide_scalar_t * dst, const typename traits::scalar_t * src, const boost::uint32_t * ip, const boost::int32_t  * vp, int nrows)
{
	int acc=0;
	for(int i = 0 ; i < nrows ; i++) {
		traits::zero(dst[i]);
		unsigned int c = 0;
		for( ; *vp != 0 ; ip++, vp++) {
			c += *ip;
			traits::addmul(dst[i], src[c], *vp);
				if (++acc == traits::max_accumulate) {
					traits::reduce(dst[i], dst[i]);
					acc = 1;
				}
		}
		ip++;
		vp++;
	}
}

#endif	/* MUL_HPP_ */
