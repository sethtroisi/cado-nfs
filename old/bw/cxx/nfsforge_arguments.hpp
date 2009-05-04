#ifndef NFSFORGE_ARGUMENTS_HPP_
#define NFSFORGE_ARGUMENTS_HPP_

#include "arguments.hpp"
#include "fmt.hpp"
#include <gmp.h>
#include <gmpxx.h>

struct nfsforge_arguments {
	mpz_class p;
	std::string pstring;
	std::string in;
	std::string out;
	std::string clink;
	std::string rlink;
	int nchar;
	nfsforge_arguments() : nchar(0) {
		in = "complete";
		out = "matrix.txt";
		clink = "clink.txt";
		rlink = "rlink.txt";
	}
	bool parse(argparser::situation& s) {
		if (s("--modulus", pstring)) return true;
		if (s("--in", in)) return true;
		if (s("--out", out)) return true;
		if (s("--clink", clink)) return true;
		if (s("--rlink", rlink)) return true;
		if (s("--nchar", nchar)) return true;
		return false;
	}
	void doc(std::ostream& o) {
		o << "--modulus <p>\t modulus\n";
		o << "--in <filename>\t input row list\n";
		o << "--out <filename>\t output matrix\n";
		o << "--clink <filename>\t output clinks\n";
		o << "--rlink <filename>\t output rlinks\n";
		o << "--nchar <integer>\t # of ideals to treat like characters\n";
	}
	bool check(std::ostream& o) {
		bool ok = true;
		if (in.empty() || out.empty() || clink.empty() || rlink.empty()) {
			o << "--in, --out, --clink, -rlink are mandatory\n";
			ok = false;
		}
		if (pstring.empty()) {
			o << "--modulus is mandatory\n";
			ok = false;
		}
		return ok;
	}
	void trigger() {
		p = mpz_class(pstring);
	}
};

#endif	/* NFSFORGE_ARGUMENTS_HPP_ */
