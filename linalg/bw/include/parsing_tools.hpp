#ifndef PARSER_BASE_HPP_
#define PARSER_BASE_HPP_

#include <istream>	/* FIXME when code goes to .cpp instead. */
#include <string>

template<typename T> struct wanted_item;

template<>
struct wanted_item<char> {
	char c;
	wanted_item(char c) : c(c) {}
	std::istream& operator>>(std::istream& is) const {
		if (is.peek() != c) {
			is.setstate(std::ios::failbit);
		} else {
			is.ignore(1);
		}
		return is;
	}
};

template<unsigned int sz>
struct wanted_item<char [sz]> {
	const char * c;
	wanted_item(const char c[]) : c(c) {}
	std::istream& operator>>(std::istream& is) const {
		is >> std::ws;
		const char * ptr = c;
		for( ; *ptr && isspace(*ptr) ; ptr++);
		for( ; *ptr ; ptr++) {
			if (is.peek() != *ptr) {
				is.setstate(std::ios::failbit);
				return is;
			}
			is.ignore(1);
		}
		return is;
	}
};

inline bool try_match(std::istream& is, char c)
{
	if (is.peek() == c) {
		is.ignore(1);
		return true;
	}
	return false;
}

template<typename T>
inline wanted_item<T> wanted(const T& c) { return wanted_item<T>(c); }

namespace std {
	template<typename T>
	inline std::istream&
	operator>>(std::istream& is, const wanted_item<T>& w) {
		return w.operator>>(is);
	}
}

inline bool startswith(const std::string& msg, const std::string& prefix)
{
	return std::string(msg, 0, prefix.length()) == prefix;
}

class comment_strip {
	std::istream& is;
	std::string clead;
public:
	comment_strip(std::istream& is, const char * c = "//")
		: is(is), clead(c) {}
	void getline(std::string& buf) {
		if (is.eof()) return;
		do {
			buf.erase();
			std::getline(is, buf);
			if (is.eof()) {
				/*
				if (!buf.empty())
					cerr << "Missing newline at end of file"
						", last line dropped\n";
						*/
				return;
			}
			if (is.fail()) {
				/* We do not want the pain of throwing an
				 * exception here */
				/* cerr << "out of string space\n"; */
				return;
			}
		} while (startswith(buf, clead));
	}
};

namespace std {
	inline void getline(comment_strip& c, std::string& b) {
		return c.getline(b);
	}
}

// void getline_nocomment(std::istream&, std::string&, const char * = "//");

#endif	/* PARSER_BASE_HPP_ */
