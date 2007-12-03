#ifndef MASTER_ARGUMENTS_HPP_
#define MASTER_ARGUMENTS_HPP_

#include "corners.hpp"
#include "arguments.hpp"
#include "fmt.hpp"

struct master_arguments {
	unsigned int threshold;
	master_arguments() {
		threshold = 900;
	}

	bool parse(argparser::situation& s) {
		if (s("--threshold", threshold)) return true;
		return false;
	}
	void doc(std::ostream& o) {
		o << "--threshold\tFFT threshold\n";
	}
	bool check(std::ostream& o) const {
		return true;
	}
	void trigger() const {}
};

#endif	/* MASTER_ARGUMENTS_HPP_ */
