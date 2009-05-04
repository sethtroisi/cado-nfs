#ifndef MUL_ARGUMENTS_HPP_
#define MUL_ARGUMENTS_HPP_

#include "arguments.hpp"
#include "files.hpp"

#include <string>
#include <vector>

struct mul_arguments {
	std::vector<std::string> file_names;
	bool integer;
	mul_arguments() {
		integer = false;
	}
	bool parse(argparser::situation& s) {
		if (s("--matrix", files::matrix)) return true;
		if (s("--integer", integer)) return true;
		file_names.push_back((const char*)s);
		return true;
	}
	void doc(std::ostream& o) {
		o << "--matrix <filename>\tUse filename instead of the default matrix.txt\n";
		o << "--integer\tCompute over Z\n";
		o << "[<fname> ... ]\tvector files\n";
	}
	bool check(std::ostream& o) {
		bool ok = true;
		if (file_names.empty()) {
			o << "file info is necessary\n";
			ok = false;
		}
		return ok;
	}
	void trigger() const { }
};

#endif	/* MUL_ARGUMENTS_HPP_ */
