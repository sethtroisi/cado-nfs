#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <set>
#include "auxfuncs.h"

#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "files.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
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

set<unsigned int> zrows;

void setup_vectors(unsigned int nrows, int m, int n, const mpz_class& p)
{
	ofstream o;

	/* as ridiculous as it seems, but we're far from critical here,
	 * and want to be able to do %p simply ! */
	mpz_class t;

	WARNING("Beware: algorithm for setting up Y has changed\n");

	for (int j = 0; j < n; j++) {
		must_open(o, files::y % j);
#if 0
		for (unsigned int c = 0; c < nrows; c++) {
			t = random() | 0x01UL;
			t %= p;
			o << t << "\n";
		}
#endif
		typedef set<unsigned int>::const_iterator suci_t;
		suci_t sc = zrows.begin();
		for(unsigned int c = 0; c < nrows; ) {
			uint next_c;
			if (sc == zrows.end()) {
				next_c = nrows;
			} else {
				next_c = *sc;
				sc++;
			}
			/* aside zero rows, there's no constraint */
			for( ; c < next_c ; c++) {
				t = random() % p;
				o << t << "\n";
			}
			/* For zero rows, we'd better stay OUT of the
			 * subspace of the image ! */
			if (c < nrows) {
				t = random() % p;
				t += (t == 0);
				o << t << "\n";
				c++;
			}
		}
		o.close();
	}

	for (int i = 0; i < m ; i++) {
		must_open(o, files::x % i);
		/* If we choose X vectors such that the Krylov subspace
		 * respective to M^T and X has low dimension, then the
		 * process will fail : the expected degree necessary for
		 * producing a dependency will be higher than N/m. This
		 * means that the quadratic bw-master will fail, and that
		 * the sub-quadratic one will crash, because the
		 * assumptions on the degree of pi will be wrong.
		 *
		 * so idx = i is clearly wrong here.
		 */
		unsigned int idx;
		for(;;) {
			idx= ((random()) % nrows);
			if (zrows.find(idx) == zrows.end()) {
				break;
			}

			cout << fmt("// Avoiding zero row %\n") % idx;
		}
		o << fmt("e%") % idx << "\n";
		o.close();
	}
}

void build_zero_table(set<unsigned int> & zrows, ifstream& is)
{
	cout << "// looking for zero rows..." << flush;

	/* clear the fail bit */
	is.clear();

	comment_strip cs(is, "//");
	string buffer;

	for(unsigned int p = 0 ; ; p++) {
		getline(cs, buffer);
		istringstream st (buffer);
		unsigned int n;
	
		if (!(st >> n)) {
			is.setstate(ios::failbit);
			break;
		} else if (n == 0) {
			zrows.insert(p);
			cout << " [" << p << "]";
		}
	}
	cout << " ok" << endl;
}

int main(int argc, char *argv[])
{
	ios_base::sync_with_stdio(false);
	cerr.tie(&cout);

	common_arguments common;
	prep_arguments mine;

	process_arguments(argc, argv, common, mine);

	unsigned int nr, nc;
	string mstr;
	ifstream mtx;
	
	/* Random state seeding is done from within the argument parser
	 */

	must_open(mtx, files::matrix);

	get_matrix_header(mtx, nr, nc, mstr);
	build_zero_table(zrows,mtx);

	mtx.close();

	setup_vectors(nc, mine.m, mine.n, mpz_class(mstr));
}


/* vim:set sw=8: */
