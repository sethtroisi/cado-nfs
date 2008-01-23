#ifndef MATRICES_PRINTING_HPP_
#define MATRICES_PRINTING_HPP_

#include <ostream>
#include "matrices.hpp"

namespace hidden {
	extern int const magma_index;
}

namespace std {
	extern std::ios_base& set_magma(std::ios_base& o);
	extern std::ios_base& unset_magma(std::ios_base& o);
	extern bool test_magma(std::ios_base& o);

	/* Specialize for one and two dimensions */
        template<typename T, typename N1>
        ostream& operator<<(ostream& o, const row<T, cons<N1, void> >& v)
        {
		typedef cons<N1, void> Dims;
                typedef typename row<const T, Dims>::value_type Z;
		if (test_magma(o)) {
			o << "Vector([KP|";
			copy(v.begin(), v.end() - 1,
				std::ostream_iterator<Z>(o, ", "));
			o << v.back() << "])";
		} else {
			copy(v.begin(), v.end() - 1,
				std::ostream_iterator<Z>(o, " "));
			o << v.back();
		}
                return o;
        }

        template<typename T, typename N1, typename N2>
        ostream& operator<<(ostream& o,
			const row<T, cons<N1, cons<N2, void> > >& v)
        {
		typedef cons<N1, cons<N2, void> > Dims;
                typedef typename row<const T, Dims>::value_type Z;
		if (test_magma(o)) {
			o << "Matrix([";
			copy(v.begin(), v.end() - 1,
				std::ostream_iterator<Z>(o, ",\n"));
			o << v.back() << "])";
		} else {
			copy(v.begin(), v.end() - 1,
				std::ostream_iterator<Z>(o, "\n"));
			o << v.back();
		}
                return o;
        }

#if 0
        /* print a matrix polynomial */
        template<typename T, typename N1, typename N2, typename N3, typename TAIL>
        std::ostream& operator<<(std::ostream& o,
                        const row<T, cons<N1, cons<N2, cons<N3, TAIL> > > >& v)
        {
                typedef cons<N1, cons<N2, cons<N3, TAIL> > > Dims;
                typedef typename row<const T, Dims>::value_type Z;
                typedef typename row<const T, Dims>::const_iterator cit_t;
                for(cit_t i = v.begin() ; i != v.end() ; i++) {
                        o << "[" << (i - v.begin()) << "]\n" << *i;
                }
                return o;
        }
#endif

        template<typename T>
        std::ostream& operator<<(std::ostream& o, const holder<T>& v)
        {
                return operator<<(o, (const typename holder<T>::super&)v);
        }
}


#endif	/* MATRICES_PRINTING_HPP_ */
