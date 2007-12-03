#ifndef MATRIX_LINE_HPP_
#define MATRIX_LINE_HPP_

#include <istream>
#include <ostream>
#include <utility>
#include <vector>
#include <boost/cstdint.hpp>

/* These structures have to do with I/O on matrices -- they are NOT
 * critical
 */
struct matrix_line : public std::vector<std::pair<boost::uint32_t, boost::int32_t> > {
};

namespace std {
	std::ostream& operator<<(std::ostream&, const matrix_line&);
	std::istream& operator>>(std::istream&, matrix_line&);
}

std::ostream& print_line_without_ones(std::ostream&, const matrix_line&);


#endif	/* MATRIX_LINE_HPP_ */
