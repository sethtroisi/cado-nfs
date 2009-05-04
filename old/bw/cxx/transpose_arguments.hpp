#ifndef TRANSPOSE_ARGUMENTS_HPP_
#define TRANSPOSE_ARGUMENTS_HPP_

#include "arguments.hpp"
#include "fmt.hpp"
#include <gmp.h>
#include <gmpxx.h>

struct transpose_arguments {
	std::string out;
	transpose_arguments() {}
	bool parse(argparser::situation& s) {
		if (s("--out", out)) return true;
		return false;
	}
	void doc(std::ostream& o) {
		o << "--out <filename>\t output matrix\n";
	}
	bool check(std::ostream& o) { return true; }
	void trigger() {
	}
};

#endif	/* TRANSPOSE_ARGUMENTS_HPP_ */
