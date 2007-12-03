#include "constants.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"

#include <cstdlib>
#include <iostream>

#include <string>
#include <sstream>

using namespace std;

int random_coeff()
{
	int b = random() % 256;
	int v = (1 - (b & 1) * 2) * ((b / 2) + 1);
	return v;
}

int main(int argc, char * argv[])
{
	ios_base::sync_with_stdio(false);
	cerr.tie(&cout);
	cout.rdbuf()->pubsetbuf(0,0);
	cerr.rdbuf()->pubsetbuf(0,0);

	if (argc != 3 && argc != 4) {
		cerr << "Usage : ./random <msize> <modulus> [<density>]\n";
		exit(1);
	}

	unsigned int nr = atoi(argv[1]);
	string mstr(argv[2]);
	double d = 0.5;
	if (argc > 3) {
		istringstream foo(argv[3]);
		foo >> d;
	}

	put_matrix_header(cout, nr, mstr);
	for(unsigned int i = 0 ; i < nr ; i++) {
		matrix_line l;
		if (i) for(unsigned int j = 0 ; j < nr ; j++) {
			double t = random() / (double) RAND_MAX;
			if (t > d) continue;
			l.push_back(make_pair(j, random_coeff()));
		}
		cout << l << "\n";
	}
	return 0;
}

