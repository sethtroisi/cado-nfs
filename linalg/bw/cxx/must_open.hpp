#ifndef MUST_OPEN_HPP_
#define MUST_OPEN_HPP_

#include <cstdio>
#include <cstdlib>
#include <cerrno>

#include "auxfuncs.h"

#include <fstream>
#include <istream>
#include <ostream>
#include <string>

#include "files.hpp"

inline void must_open(std::ifstream& o, const char * s)
{
	o.open(s);
	if (!o.is_open()) die("%s : %s", 1, s, strerror(errno));
}
inline void must_open(std::ofstream& o, const char * s)
{
	o.open(s);
	if (!o.is_open()) die("%s : %s", 1, s, strerror(errno));
}
inline void must_open(std::ofstream& o, const std::string& s)
{ must_open(o, s.c_str()); }
inline void must_open(std::ifstream& i, const std::string& s)
{ must_open(i, s.c_str()); }
/*
inline void must_open(std::ofstream& o, const files::meta_filename& s)
{ must_open(o, (std::string) s); }
inline void must_open(std::ifstream& i, const files::meta_filename& s)
{ must_open(i, (std::string) s); }
*/
inline bool open(std::ofstream& o, const char * s)
{
	o.open(s);
	return o.is_open();
}
inline bool open(std::ifstream& i, const char * s)
{
	i.open(s);
	return i.is_open();
}
inline bool open(std::ofstream& o, const std::string& s)
{ return open(o, s.c_str()); }
inline bool open(std::ifstream& i, const std::string& s)
{ return open(i, s.c_str()); }
/*
inline bool open(std::ofstream& o, const files::meta_filename& s)
{ return open(o, (std::string) s); }
inline bool open(std::ifstream& i, const files::meta_filename& s)
{ return open(i, (std::string) s); }
*/

#endif	/* MUST_OPEN_HPP_ */
