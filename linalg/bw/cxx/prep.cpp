#include <cstdio>
#include <cstdlib>
#include "auxfuncs.h"

#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "files.hpp"
#include "matrix_header.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"
#include "prep_arguments.hpp"
#include <gmp.h>
#include <gmpxx.h>

#include <iostream>
#include <fstream>

/* bw-prep
 *
 * This program is responsible of the following
 *
 * - Choose X vectors
 * - Choose random Y vectors
 */

using namespace std;

void setup_vectors(unsigned int nrows, int m, int n, const mpz_class& p)
{
	ofstream o;

	mpz_class t;

	for (int j = 0; j < n; j++) {
		must_open(o, files::y % j);
		for (unsigned int c = 0; c < nrows; c++) {
			t = random() | 0x01UL;
			t %= p;
			o << t << "\n";
		}
		o.close();
	}

	for (int i = 0; i < m ; i++) {
		must_open(o, files::x % i);
#if 1
		o << fmt("e%") % i << "\n";
#else
		o << fmt("e%") % ((rand() + time(NULL)) % nrows) << "\n";
#endif
		o.close();
	}
}

int main(int argc, char *argv[])
{
	ios_base::sync_with_stdio(false);
	cerr.tie(&cout);

	common_arguments common;
	prep_arguments mine;

	process_arguments(argc, argv, common, mine);

	unsigned int nr;
	string mstr;
	ifstream mtx;
	
	must_open(mtx, files::matrix);

	get_matrix_header(mtx, nr, mstr);

	mtx.close();

	setup_vectors(nr, mine.m, mine.n, mpz_class(mstr));
}


/* vim:set sw=8: */
