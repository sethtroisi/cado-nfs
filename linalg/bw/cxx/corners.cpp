#include "corners.hpp"
#include "parsing_tools.hpp"

std::istream& std::operator>>(std::istream& i, corner& c)
{
	i >> c.first >> wanted(',') >> std::noskipws >> c.second;
	if (i.fail()) {
		c.first = c.second = -1;
	}
	return i;
}
	
std::ostream& std::operator<<(std::ostream& o, const corner& c)
{
	return o << c.first << ',' << c.second;
}
