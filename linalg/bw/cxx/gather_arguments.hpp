#ifndef GATHER_ARGUMENTS_HPP_
#define GATHER_ARGUMENTS_HPP_

#include "arguments.hpp"
#include "fmt.hpp"

#include <string>
#include <vector>

struct gather_arguments {
	unsigned int scol;
	std::vector<std::string> file_names;
	gather_arguments() {
		scol = (unsigned int) -1;
	}

	bool parse(argparser::situation& s) {
		if (file_names.empty() && s(NULL, scol)) return true;
		file_names.push_back((const char*)s);
		return true;
	}
	void doc(std::ostream& o) {
		o << "<col>\tbase of solution\n";
		o << "[<fname> ... ]\tfragment files\n";
	}
	bool check(std::ostream& o) {
		bool ok = true;
		if (scol == (unsigned int) -1 && file_names.empty()) {
			o << "file info is necessary\n";
			ok = false;
		}
		return ok;
	}
	void trigger() const { }
};

#endif	/* GATHER_ARGUMENTS_HPP_ */
