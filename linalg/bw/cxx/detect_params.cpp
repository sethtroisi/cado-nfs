#include "files.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"

#include <string>
#include <istream>
#include <fstream>
#include <sstream>

using namespace std;

void detect_mn(unsigned int& m, unsigned int& n)
{
	ifstream f;

	m = n = 0;

	/* TODO : Support a proper config file instead of this shit
	 */
	if (open(f, files::params)) {
		for(;;) {
			string buf;
			getline(f, buf);
			if (f.eof()) break;
			istringstream b(buf);
			b >> wanted("n") >> wanted("=") >> n;
			b.seekg(0); b.clear();
			b >> wanted("m") >> wanted("=") >> m;
		}
	}

	if (n == 0) for(n = 0 ; ; n++) {
		if (!open(f, files::y % n)) {
			break;
		}
		f.close();
	}

	if (m == 0) for(m = 0 ; ; m++) {
		if (!open(f, files::x % m)) {
			break;
		}
		f.close();
	}
}
