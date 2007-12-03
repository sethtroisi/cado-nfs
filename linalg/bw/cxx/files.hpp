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

	/*
	* std::string wip_tag	INIT_WITH("WORKING");
	* std::string matrix	INIT_WITH("matrix");
	* std::string indexes	INIT_WITH("indexes");
	* std::string values	INIT_WITH("values");
	std::string w		INIT_WITH("W");
	* std::string pre_stats	INIT_WITH("FIRST.INFO");
	meta_filename a_sub	INIT_WITH("sA-%[f0w2]-%[f0w2].%[f0w3]-%[f0w3]");
	meta_filename s		INIT_WITH("S-%[f0w2].%[f0w3]");
	meta_filename ds	INIT_WITH("ndS-%[f0w2].%[f0w3]-%[f0w3]");
	meta_filename v		INIT_WITH("V-%[f0w2].%[f0w3]");
	meta_filename f		INIT_WITH("rF-%[f0w2]-%[f0w2].%[f0w4]");
	meta_filename f_sub	INIT_WITH("srF-%[f0w2]-%[f0w2].%[f0w4].%-%");
	meta_filename f_base	INIT_WITH("F_INIT");
	meta_filename pi	INIT_WITH("P-%-%");
	meta_filename h		INIT_WITH("H.%[f0w2]");
	* meta_filename l		INIT_WITH("LOCAL_INFO-%[f0w2]");
	meta_filename x		INIT_WITH("X%[f0w2]");
	meta_filename y		INIT_WITH("Y%[f0w2]");
	meta_filename z		INIT_WITH("Z%[f0w2]");
	meta_filename valu	INIT_WITH("VALUATION-%[f0w2].%[f0w4]");
	* meta_filename i_matrix	INIT_WITH("matrix.%d");
	* meta_filename i_vector	INIT_WITH("vector.%d");
	* meta_filename certif	INIT_WITH("CERTIF-%-%-%");
	*/
	EXTERN_LINKAGE	std::string matrix      INIT_WITH("matrix.txt");
	EXTERN_LINKAGE	std::string params      INIT_WITH("bw.cfg");
	EXTERN_LINKAGE	meta_filename a		INIT_WITH("A-%[f0w2]-%[f0w2]");
	EXTERN_LINKAGE	meta_filename v		INIT_WITH("V%[f0w2].%[f0w3]");
	EXTERN_LINKAGE	meta_filename w		INIT_WITH("W%[f0w2]");
	EXTERN_LINKAGE	meta_filename x		INIT_WITH("X%[f0w2]");
	EXTERN_LINKAGE	meta_filename y		INIT_WITH("Y%[f0w2]");
	EXTERN_LINKAGE	meta_filename z		INIT_WITH("Z%[f0w2]");
	EXTERN_LINKAGE	std::string precond	INIT_WITH("precond.txt");
	EXTERN_LINKAGE	std::string m0		INIT_WITH("M0-vector");
	EXTERN_LINKAGE	std::string x0		INIT_WITH("X0-vector");
	EXTERN_LINKAGE	meta_filename f		INIT_WITH("F%[f0w2]");
	EXTERN_LINKAGE	meta_filename f_init	INIT_WITH("F_INIT");
	EXTERN_LINKAGE	meta_filename fxy
				INIT_WITH("F%[f0w2]-Y%[f0w2].%[f0w3]");

	/* TODO: maybe bump up the format width a wee bit, or drop the f0
	 * thing.
	 */
}

#undef	INIT_WITH
#undef	EXTERN_LINKAGE

#endif	/* FILES_HPP_ */
