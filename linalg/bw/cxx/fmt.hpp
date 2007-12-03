#ifndef FORMATTER_HPP_
#define FORMATTER_HPP_

#include <iomanip>
#include <string>
#include <sstream>
#include <ostream>
#include <stdexcept>

/* Provides type-safe printf-like formatting for C++.
 *
 * Use as : cout << fmt("Got % shells in % buckets") % nsh % nb << "\n";
 *
 * %% is a literal %
 * %[...] is the extended form for the formatting specifier, which has
 * the following form:
 * (<letter><argument>)*
 *
 * <letter> <argument> <effect>
 *     f      char     setfill
 *     w      uint     setw
 *
 * Extensions can come easily, of course.
 */

/* boost::format is much, much better than this, but also much,
 * much bigger.
 */
class fmt {
	struct error : public std::runtime_error {
		inline error(const char s[], const std::string& v)
			: std::runtime_error(v + ": " + s) {}
	};
	std::string s;
	int pos;
	public:
	fmt(const std::string& c, int p = 0) : s(c), pos(p) {}
	fmt(const char * c, int p = 0) : s(c), pos(p) {}
	template<class T>
	fmt operator%(const T& arg) const {
		std::string::size_type s0, s1, s2;
		std::ostringstream m;
		s0 = s.find('%', (std::string::size_type) pos);
		m << s.substr(0, s0);
		for(; s0 != s.npos && s0 + 1 < s.size() && s[s0 + 1] == '%'; ) {
			m << '%';
			s0 += 2;
			if (s0 == s.npos)
				throw error("too many arguments", s);
			s1 = s.find('%', s0);
			m << s.substr(s0, s1 - s0);
			s0 = s1;
		}
		if (s0 == s.npos)
			throw error("too many arguments", s);
		s1 = s0 + 1;
		s2 = s0;
		std::string fm;
		if (s1 < s.size() && s[s1] == '[') {
			s2 = s.find(']', s1+1);
			if (s2 == s.npos) {
				throw error("bad pattern", s);
			}
			fm = s.substr(s1+1, s2-s1-1);
		}

		std::istringstream f(fm);
		for(;!f.eof();) {
			char o;
			int w;
			char fl;
			std::streampos oldpos;
			f >> o;
			if (f.eof()) break;
			switch(o) {
				case 'f':
					f >> fl;
					m << std::setfill(fl);
					break;
				case 'w':
					oldpos = f.tellg();
					f >> w;
					if (f.tellg() > oldpos)
						f.clear();
					m << std::setw(w);
					break;
				case  'h':
					m << std::hex;
					break;
				case  'F':
					m << std::fixed;
					break;
				case '.':
					oldpos = f.tellg();
					f >> w;
					if (f.tellg() > oldpos)
						f.clear();
					m << std::setprecision(w);
					break;
				default:
					throw error("bad format", s);
			}
			if (f.fail()) {
				throw error("bad format", s);
			}
		}
		m << arg;
		m.flush();
		int newpos = m.str().size();
		m << s.substr(s2 + 1);
		return fmt(m.str(), newpos);
	}
	operator std::string() const {
		// do not forget to transform %%'s in the string tail
		//
		std::string v = s;
		std::string::size_type s0;
		for(s0 = pos ; (s0 = v.find("%", s0)) < v.size() ; ) {
			if (v[s0 + 1] != '%')
				throw error("too few arguments", s);
			v = v.erase(++s0, 1);
		}
		return v;
	}
};


namespace std {
	inline std::ostream& operator<<(std::ostream& o, const fmt& f)
	{
		return o << (std::string) f;
	}
}

#endif	/* FORMATTER_HPP_ */
