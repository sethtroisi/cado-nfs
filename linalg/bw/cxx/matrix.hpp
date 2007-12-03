#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <istream>
#include <vector>
#include <boost/cstdint.hpp>

extern void fill_matrix_data(std::istream&,
		std::streampos, uint32_t,
		uint, uint, 
		uint32_t *, int32_t *);
extern void count_matrix_coeffs(std::istream&,
		unsigned int,
		std::vector<std::streampos>&,
		std::vector<boost::uint32_t>&);

#endif	/* MATRIX_HPP_ */
