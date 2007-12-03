#ifndef CORNERS_HPP_
#define CORNERS_HPP_

#include <utility>
#include <istream>
#include <ostream>

struct corner : public std::pair<int, int> {
	corner() : std::pair<int, int>(-1, -1) {}
	bool unset() const { return first < 0 && second < 0; }
};

namespace std {
std::istream& operator>>(std::istream&, corner&);
std::ostream& operator<<(std::ostream&, const corner&);
}


#endif	/* CORNERS_HPP_ */
