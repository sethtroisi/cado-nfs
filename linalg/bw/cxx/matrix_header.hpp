#ifndef MATRIX_HEADER_HPP_
#define MATRIX_HEADER_HPP_

#include <istream>
#include <string>

extern void get_matrix_header(std::istream & mtx,
		unsigned int &nr,
		std::string & mstr);
extern void get_matrix_header(std::istream & mtx,
		unsigned int &nr,
		unsigned int &nc,
		std::string & mstr);
extern void put_matrix_header(std::ostream & mtx,
		unsigned int nr,
		const std::string & mstr);

#endif	/* MATRIX_HEADER_HPP_ */
