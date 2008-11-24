#ifndef PREP_ARGUMENTS_HPP_
#define PREP_ARGUMENTS_HPP_

#include "arguments.hpp"

#include <cstdlib>
#include <ctime>
#include <gmpxx.h>

#include "defaults.h"

struct prep_arguments {
	int m, n;
	unsigned int seed;
        mpz_class p;
        int check_interval;
	prep_arguments() :
		m(0), n(0), seed(0), p(2),
                check_interval(CHECK_INTERVAL_DEFAULT)
	{}
	bool parse(argparser::situation& s) {
		if (s("-m", m) || s("m=", m)) return true;
		if (s("-n", n) || s("n=", n)) return true;
		if (s("--mn", n) || s("-mn", n) || s("mn=", n)) {
                    m = n;
                    return true;
                }
		if (s("-p", p) || s("--modulus", p) || s("p=",p)) return true;
		if (s("--seed", seed)) return true;
		if (s("--check-interval", check_interval)) return true;
		return false;
	}
	void doc(std::ostream& o) {
		o
		<< "[--seed <s>]\tUse s as a seed value\n"
		<< "[-m <m>]\tnumber or blocking rows\n"
		<< "[-n <n>]\tnumber or blocking columns\n"
		<< "[-p <p> | --modulus <p>]\tnumber or blocking columns\n"
		<< "[--check-interval <n>]"
                    << "\tDo dot-product check every <n> iterations\n"
                ;
	}
	bool check(std::ostream& o MAYBE_UNUSED, bool can_print MAYBE_UNUSED)
        { return true; }
	void trigger(bool can_print MAYBE_UNUSED) const { }
};

#endif	/* PREP_ARGUMENTS_HPP_ */
