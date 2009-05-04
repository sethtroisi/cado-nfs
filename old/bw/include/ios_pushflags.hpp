#ifndef IOS_PUSHFLAGS_HPP_
#define IOS_PUSHFLAGS_HPP_

#include <ios>

class ios_pushflags {
	std::ios_base * parent;
	std::ios_base::fmtflags f;
public:
	ios_pushflags(std::ios_base& io) : parent(&io), f(io.flags()) {}
	~ios_pushflags() { parent->flags(f); }
};

	
#endif	/* IOS_PUSHFLAGS_HPP_ */
