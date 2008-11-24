#ifndef ARGUMENTS_HPP_
#define ARGUMENTS_HPP_

/* This is a tentative port of the curves/cab argument framework. It's
 * got its rough edges, and could be improved.
 */

#include <string>
#include <sstream>

#include <iostream>	/* for cerr in process_arguments */
#include <cstdlib>	/* for exit() in process_arguments */
#include <cstring>	/* for strlen() */

#include <vector>
#include <cstring>

#include "manu.h"	/* for MAYBE_UNUSED */

namespace argparser {
class situation {
	public:
	int argc;
	char ** p;
	private:
	bool check(const char * opt, int k) {
                if (argc == 0) return false;
		if (opt == NULL)
			return argc >= k;
                size_t olen = strlen(opt);
                if (opt[olen-1] == '=' && std::string(p[0],olen) == opt) {
                    /* Special case for arguments of the form x=123 */
                    p[0] += olen;
                    // for( ; isspace(p[0][0]) && p[0][0] ; p[0]++) ;
                } else {
                    if (operator != (opt)) return false;
                    operator++();
                }
		if (argc < k) {
			std::ostringstream m;
			m << opt;
			m << " needs ";
			if (k == 1) m << "an";
			else m << k;
			m << " argument";
			if (k > 1) m << "s";
			m.flush();
			throw m.str().c_str();
		}
		return true;
	}
	public:
	situation& operator++() { argc--, p++; return *this; }
	situation operator++(int) { situation t = *this; ++(*this); return t; }
	operator const char *() const {
		if (argc == 0) return NULL; else return p[0];
	}
	bool operator==(const char *x) const { return std::string(x) == p[0]; }
	bool operator==(std::string & x) const { return x == p[0]; }
	bool operator!=(const char *x) const { return !operator==(x); }
	bool operator!=(std::string & x) const { return !operator==(x); }
	template<typename T>
	bool operator()(const char * opt, std::vector<T> & var) {
		if (!check(opt, 1)) return false;
		std::istringstream s(p[0]);
                T v;
		s >> v;
                var.push_back(v);
		return s.eof() && !s.fail();
	}
	template<typename T>
	bool operator()(const char * opt, T & var) {
		if (!check(opt, 1)) return false;
		std::istringstream s(p[0]);
		s >> var;
		return !s.fail()
                    && s.peek() == std::istringstream::traits_type::eof();
	}
	bool operator()(const char * opt, char const * & var) {
		if (!check(opt, 1)) return false;
		var = p[0];
		return true;
	}
	bool operator()(const char * opt, std::string & var) {
		if (!check(opt, 1)) return false;
		var = p[0];
		return true;
	}
	bool operator()(const char * opt, bool & var) {
		if (operator==(opt)) { var = true; return true; }
		/* check whether we can form a --no-blah argument */
		if (strncmp(opt,"--", 2) != 0) return false;
		std::string v(opt + 2);
		if ((std::string("--no-") + v) == p[0]) {
			var = false; return true;
		}
		return false;
	}
};

}	/* namespace argparser */

struct no_arguments {
	bool parse(argparser::situation& s MAYBE_UNUSED) { return false; }
	void doc(std::ostream& o) const {
		o
		<< "Accepted options:\n"
		<< "--help\t\tshow this help\n";
	}
	bool check(std::ostream& o MAYBE_UNUSED) { return true; }
	void trigger() const {}
};

template<class common, class special>
void process_arguments(int argc, char * argv [], common& com, special& spe,
        bool can_print = true)
{
    using namespace argparser;
    situation s;

    s.p = argv;
    s.argc = argc;
    s++;
    for(;s;s++) {
        if (s == "--help") {
            if (can_print) {
                com.doc(std::cerr);
                spe.doc(std::cerr);
            }
            exit(0);
        }
        try {
            if (com.parse(s)) continue;
            if (spe.parse(s)) continue;
        } catch (const char * err) {
            if (can_print) {
                std::cerr << "ERROR: " << err << "\n";
                com.doc(std::cerr);
                spe.doc(std::cerr);
            }
            exit(1);
        }
        if (can_print) {
            std::cerr << s.p[0] << " : unexpected option\n";
            com.doc(std::cerr);
            spe.doc(std::cerr);
        }
        exit(1);
    }
    bool ok = true;
    ok = ok & com.check(std::cerr, can_print);
    ok = ok & spe.check(std::cerr, can_print);
    if (!ok) {
        if (can_print) {
            std::cerr << "\n";
            com.doc(std::cerr);
            spe.doc(std::cerr);
        }
        exit(1);
    }
    com.trigger(can_print);
    spe.trigger(can_print);
}
#endif	/* ARGUMENTS_HPP_ */
