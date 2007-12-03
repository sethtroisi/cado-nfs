#ifndef PREP_ARGUMENTS_HPP_
#define PREP_ARGUMENTS_HPP_

#include "arguments.hpp"

struct prep_arguments {
	int m, n;
	prep_arguments() :
		m(-1), n(-1)
	{}
	bool parse(argparser::situation& s) {
		if (m == -1 && s(NULL, m)) return true;
		if (n == -1 && s(NULL, n)) return true;
		return false;
	}
	void doc(std::ostream& o) {
		o
		<< "<m>\tnumber or blocking rows\n"
		<< "<n>\tnumber or blocking columns\n";
	}
	bool check(std::ostream& o) {
		bool b = true;
		if (m == -1 || n == -1) {
			o << "<m> and <n> are mandatory\n";
			b = false;
		}
		return b;
	}
	void trigger() const {}
};

#endif	/* PREP_ARGUMENTS_HPP_ */
