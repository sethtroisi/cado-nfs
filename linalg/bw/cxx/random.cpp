#include "constants.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"

#include <cstdlib>
#include <iostream>

#include <string>
#include <sstream>

#include <gmp.h>
#include <gmpxx.h>

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

	unsigned int seed = time(NULL);

	if (argc >= 3 && strcmp(argv[1], "--seed") == 0) {
		seed = atoi(argv[2]);
		argv++,argc--;
		argv++,argc--;
	}
	if (argc != 3 && argc != 4) {
		cerr << "Usage : ./random [--seed <seed>] <msize> <modulus> [<density>]\n";
		exit(1);
	}
	srandom(seed);

	unsigned int nr = atoi(argv[1]);
	string mstr(argv[2]);
	double d = 0.5;
	if (argc > 3) {
		istringstream foo(argv[3]);
		foo >> d;
	}
	if (d > 1) { d = d/nr; }

	mpz_class px(mstr);

	put_matrix_header(cout, nr, mstr);
	if (px > 128) {
		for(unsigned int i = 0 ; i < nr ; i++) {
			matrix_line l;
			if (i == 0) {
				cout << l << "\n";
				continue;
			}
			for(unsigned int j = 0 ; j < nr ; j++) {
				double t = random() / (double) RAND_MAX;
				if (t > d) continue;
				l.push_back(make_pair(j, random_coeff()));
			}
			cout << l << "\n";
		}
	} else if (px != 2) {
		int p = px.get_si();
		for(unsigned int i = 0 ; i < nr ; i++) {
			matrix_line l;
			if (i == 0) {
				cout << l << "\n";
				continue;
			}
			for(unsigned int j = 0 ; j < nr ; j++) {
				double t = random() / (double) RAND_MAX;
				if (t > d) continue;
				int v = random_coeff() % p;
				if (v < (-p/2)) {
					v += p;
				} else if (v > p/2) {
					v -= p;
				}
				if (v == 0) {
					continue;
				}
				l.push_back(make_pair(j, v));
			}
			cout << l << "\n";
		}
	} else {
		for(unsigned int i = 0 ; i < nr ; i++) {
			matrix_line l;
			if (i == 0) {
				cout << 0 << "\n";
				continue;
			}
			for(unsigned int j = 0 ; j < nr ; j++) {
				double t = random() / (double) RAND_MAX;
				if (t > d) continue;
				if (random_coeff() & 1 == 0)
					continue;
				l.push_back(make_pair(j, 1));
			}
			print_line_without_ones(cout,l) << "\n";
		}
	}
	return 0;
}

