#ifndef LINGEN_ARGUMENTS_HPP_
#define LINGEN_ARGUMENTS_HPP_

#include "corners.hpp"
#include "arguments.hpp"
#include "fmt.hpp"

struct lingen_arguments {
	unsigned int threshold;
	lingen_arguments() {
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

#endif	/* LINGEN_ARGUMENTS_HPP_ */
