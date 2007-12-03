#ifndef COMMON_ARGUMENTS_HPP_
#define COMMON_ARGUMENTS_HPP_

#include "auxfuncs.h"

#include <iostream>
#include <string>
#include "arguments.hpp"

#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <cstdlib>

struct common_arguments {
	std::string subdir;
	bool core_ok;
	bool skip_lim_check;
	common_arguments() {
		core_ok = skip_lim_check = false;
	}
	bool parse(argparser::situation& s) {
		if (s("--subdir", subdir)) return true;
		if (s("--core-ok", core_ok)) return true;
		if (s("--skip-limits-check", skip_lim_check)) return true;
		return false;
	}
	void doc(std::ostream& o) {
		o << "Accepted options:\n";
		o << "--help\t\tshow this help\n";
		o << "--subdir <dir>\tdirectory of work files\n";
		o << "--core-ok\t\tallow core dumps\n";
		o << "--skip-limits-check\tdo not check OS limits\n";
	}
	bool check(std::ostream& o) {
		bool b = true;
		if (subdir.empty()) {
			o << "--subdir is mandatory\n";
			b = false;
		}
		return b;
	}
	void trigger() const {
		using namespace std;
		struct stat sbuf[1];
		int rc;

		rc = stat(subdir.c_str(), sbuf);
		
		if (rc < 0 && errno == ENOENT) {
			cerr << subdir << " not found, creating\n";
			rc = mkdir(subdir.c_str(), 0777);
			if (rc == 0)
				rc = stat(subdir.c_str(), sbuf);
		}
		if (rc < 0) {
			cerr << subdir << " can not be created :"
				<< strerror(errno) << "\n";
			exit(1);
		}
		cout << "// changing to directory " << subdir.c_str() << "\n";
		rc = chdir(subdir.c_str());
		if (rc < 0) {
			cerr << "can not chdir to " << subdir << ": "
				<< strerror(errno) << "\n";
			exit(1);
		}
		coredump_limit(core_ok);
		if (!skip_lim_check) {
			check_limit_requirements();
		}
	}
};

#endif	/* COMMON_ARGUMENTS_HPP_ */
