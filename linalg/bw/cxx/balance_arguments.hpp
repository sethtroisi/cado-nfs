#ifndef BALANCE_ARGUMENTS_HPP_
#define BALANCE_ARGUMENTS_HPP_


#include "arguments.hpp"

#include <cstdlib>
#include <ctime>

struct balance_arguments {
    unsigned int nbuckets;
	balance_arguments() : nbuckets(24) {}
	bool parse(argparser::situation& s) {
		if (s("--nbuckets", nbuckets)) return true;
		return false;
	}
	void doc(std::ostream& o) {
		o
		<< "[--nbuckets <s>]\tBalance for s buckets instead of default 24.\n";
	}
	bool check(std::ostream& o) {
		return true;
	}
	void trigger() const {
	}
};

#endif	/* BALANCE_ARGUMENTS_HPP_ */
