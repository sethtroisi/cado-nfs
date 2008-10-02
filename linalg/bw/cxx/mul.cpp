#define __STDC_LIMIT_MACROS
#include <cstdio>
#include <cstdlib>
#include "auxfuncs.h"
#include "manu.h"
#include "gmp-hacks.h"

#include "addmul.hpp"
#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "mul_arguments.hpp"
#include "files.hpp"
#include "matrix.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"
#include "state.hpp"
#include "threads.hpp"
#include "ticks.hpp"
#include "traits.hpp"

#include <stdint.h>
#include <iterator>

/* This is a standalone program that performs an example matrix
 * multiplication */

mul_arguments mine;

using namespace std;
// using namespace boost;
using namespace core_ops;

typedef unsigned int uint;

typedef variable_scalar_traits traits;

namespace globals {
	uint nr, nc;
	mpz_class modulus;

	uint nb_coeffs;

	traits::representation::matrix_rowset mat;
}

typedef traits::scalar_t scalar_t;
typedef traits::wide_scalar_t wide_scalar_t;

scalar_t *v;
wide_scalar_t *w;

void init_data()
{
	using namespace globals;
	v = new scalar_t[nc];
	w = new wide_scalar_t[nr];

	traits::zero(v, nc);
	traits::zero(w, nr);
}

void clear_data()
{
	using namespace globals;
	delete[] v;
	delete[] w;

}

void read(const std::string& fn)
{
	ifstream f;
	must_open(f, fn);
	using namespace globals;

	std::vector<mpz_class> foo(1);
	unsigned int i;
	for(i = 0 ; i < nc ; i++) {
		if (!(f >> foo[0])) {
			break;
		}
		traits::assign(v[i], foo, 0);
	}
	if (i < nc) {
		cerr << "Vector " << fn << " has only " << i << " coordinates\n";
	}
	foo[0]=0;
	for( ; i < nc ; i++) {
		traits::assign(v[i], foo, 0);
	}

	f >> ws;
	BUG_ON(!f.eof());
}

/*
bool is_zero(const scalar_t * vec) {
	for(uint i = 0 ; i < globals::nr ; i++) {
		if (!traits::is_zero(vec[i])) {
			return false;
		}
	}
	return true;
}
*/

void one(std::string const& s)
{
	std::string outfile;
	using namespace globals;

	read(s);
	mat.mul<traits>(w,v);
	if (!mine.integer) {
		traits::reduce(w, w, 0, nr);
		// we prefer to have something which does not care on the
		// input matrix shape.
		// traits::zero(w + nr, nc - nr);
	}

	outfile = s;
	outfile.append(".mul");

	ofstream wf;
	must_open(wf, outfile);
	for(uint i = 0 ; i < globals::nr ; i++) {
		traits::print(wf,w[i]) << "\n";
	}
	wf.close();
}

void go()
{
	for(unsigned int i = 0 ; i < mine.file_names.size() ; i++) {
		one(mine.file_names[i]);
	}
}

int main(int argc, char *argv[])
{
	ios_base::sync_with_stdio(false);
	cerr.tie(&cout);
	cout.rdbuf()->pubsetbuf(0,0);
	cerr.rdbuf()->pubsetbuf(0,0);

	common_arguments common;

	process_arguments(argc, argv, common, mine);

	string mstr;
	ifstream mtx;
	
	using namespace globals;

	must_open(mtx, files::matrix);

        matrix_stats stat;
        std::vector<matrix_slice> slices(1);

        stat(mtx);
	
	cout.flush();
	cerr.flush();

        mat.fill(mtx, slices[0]);
	mtx.close();

	cout << "// Using generic code\n";
	init_data();
	go();
	clear_data();
}

/* vim:set sw=8: */
