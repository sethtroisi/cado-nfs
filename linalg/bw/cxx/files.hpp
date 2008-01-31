#ifndef FILES_HPP_
#define FILES_HPP_

#include "fmt.hpp"

#ifndef	INIT_WITH
#define	INIT_WITH(x)	/**/
#endif

#ifndef	EXTERN_LINKAGE
#define	EXTERN_LINKAGE	extern
#endif

namespace files {
	struct meta_filename : public fmt {
		meta_filename(const fmt& c) : fmt(c) {}
		meta_filename(const char * c) : fmt(c) {}
		template<class T>
		meta_filename operator%(const T& arg) const {
			return meta_filename(((const fmt&)*this) % arg);
		}
		std::string str() const { return (std::string) *this; }
		// do not mess with this one ! str() is a temporary,
		// therefore you must not assign the result. Assign the
		// string instead, and use .c_str() each time.
		const char * c_str() const { return str().c_str(); }
	};

	EXTERN_LINKAGE	std::string matrix      INIT_WITH("matrix.txt");
	EXTERN_LINKAGE	std::string params      INIT_WITH("bw.cfg");
	EXTERN_LINKAGE	meta_filename a		INIT_WITH("A-%[f0w2]-%[f0w2]");
	EXTERN_LINKAGE	meta_filename v		INIT_WITH("V%[f0w2].%[f0w3]");
	EXTERN_LINKAGE	meta_filename w		INIT_WITH("W%[f0w2]");
	EXTERN_LINKAGE	meta_filename r		INIT_WITH("R%[f0w2]");
	EXTERN_LINKAGE	std::string x		INIT_WITH("X");
	EXTERN_LINKAGE	std::string y		INIT_WITH("Y");
	EXTERN_LINKAGE	std::string precond	INIT_WITH("precond.txt");
	EXTERN_LINKAGE	std::string m0		INIT_WITH("M0-vector");
	EXTERN_LINKAGE	std::string x0		INIT_WITH("X0-vector");
	EXTERN_LINKAGE	meta_filename f		INIT_WITH("F%[f0w2]");
	EXTERN_LINKAGE	std::string f_init	INIT_WITH("F_INIT");
	EXTERN_LINKAGE	std::string f_initq	INIT_WITH("F_INIT_QUICK");
	EXTERN_LINKAGE	meta_filename fxy
				INIT_WITH("F%[f0w2]-Y%[f0w2].%[f0w3]");
	EXTERN_LINKAGE	meta_filename pi	INIT_WITH("pi-%-%");

	/* TODO: maybe bump up the format width a wee bit, or drop the f0
	 * thing.  */
}

#undef	INIT_WITH
#undef	EXTERN_LINKAGE

#endif	/* FILES_HPP_ */
