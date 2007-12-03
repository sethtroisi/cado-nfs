#include <cstdio>
#include <cstdlib>
#include "auxfuncs.h"
#include "manu.h"
#include "gmp-hacks.h"

#include "addmul.hpp"
#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "detect_params.hpp"
#include "gather_arguments.hpp"
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
#include "mul.hpp"

#include <boost/cstdint.hpp>
#include <iterator>

/*
 * Note on vectorization. Unlike slave and mksol, this program is NOT
 * vectorized. Mostly because it does not make much sense, since this one
 * is so simple. We even use generic code inconditionnally.
 *
 * One way to use this program in a vectorized manner would be to gather
 * several solutions at once, but again it would not make much sense.
 *
 * This being said, a --nbys option exists, and is mandatory in order to
 * read from the output of a vectorizaed mksol.
 */

gather_arguments mine;

using namespace std;
using namespace boost;
using namespace core_ops;

typedef unsigned int uint;

namespace globals {
	uint m, n;
	uint col;
	uint nr;
	int nbys;
	mpz_class modulus;
	uint8_t modulus_u8;
	uint16_t modulus_u16;
	uint32_t modulus_u32;
	unsigned long modulus_ulong;

	uint nb_coeffs;

	uint32_t * idx;
	int32_t  * val;
}

template<typename traits>
struct program {
	typedef typename traits::scalar_t scalar_t;
	typedef typename traits::wide_scalar_t wide_scalar_t;

	scalar_t *v;
	scalar_t *w;
	wide_scalar_t *scrap;

void init_data()
{
	using namespace globals;
	v = new scalar_t[nr];
	w = new scalar_t[nr];
	scrap = new wide_scalar_t[nr];

	traits::zero(v, nr);
	traits::zero(w, nr);
}

void clear_data()
{
	using namespace globals;
	delete[] v;
	delete[] w;
	delete[] scrap;

}

void add_file(const std::string& fn)
{
	ifstream f;
	must_open(f, fn);
	using namespace globals;

	for(unsigned int i = 0 ; i < nr ; i++) {
		/* This really works only with one member, as the
		 * operation is not too much focused on vectoring things
		 * so far.
		 */
		std::vector<mpz_class> foo(1);
		traits::assign(foo, v[i]);
		std::string line;
		getline(f, line);
		istringstream ss(line);
		int ny;
		mpz_class blah;
		for(ny = 0 ; ss >> blah ; foo[0] += blah, ny++);
		BUG_ON(ny != nbys);
		traits::assign(v[i], foo, 0);
	}
	f >> ws;
	BUG_ON(!f.eof());
}

void do_sum()
{
	using namespace globals;
	if (mine.file_names.empty()) {
		for(uint j = 0 ; j < m ; j++) {
			for(uint r = 1 ; ; r++) {
				ifstream fxy;
				string nm = files::fxy % mine.scol % j % r;
				if (!open(fxy, nm))
					break;
				mine.file_names.push_back(nm);
				cout << "// selecting " << nm << endl;
			}
		}
	}

	BUG_ON(mine.file_names.empty());

	typedef vector<string>::const_iterator vsci_t;
	for(vsci_t it = mine.file_names.begin() ; 
			it != mine.file_names.end() ;
			it++)
	{
		add_file(*it);
	}
}

bool is_zero(const scalar_t * vec) {
	for(uint i = 0 ; i < globals::nr ; i++) {
		if (!traits::is_zero(vec[i])) {
			return false;
		}
	}
	return true;
}

void multiply()
{
	multiply_ur<traits>(scrap, v, globals::idx, globals::val, globals::nr);
	traits::reduce(w, scrap, 0, globals::nr);
}

program() { init_data(); }
~program() { clear_data(); }

void go()
{
	using namespace globals;

	init_data();
	do_sum();

	traits::reduce(v, v, 0, nr);

	if (is_zero(v)) {
		cerr << "// trivial solution encountered !\n";
		BUG();
	}
	traits::copy(w, v, nr);

	uint i;

	for(i = 0 ; i < 100 && !is_zero(w) ; i++) {
		traits::copy(v, w, nr);
		cout << fmt("// trying B^%*y") % i << endl;
		multiply();
	}

	if (i == 100 && !is_zero(w)) {
		cerr << "// no solution found\n";
		BUG();
	}

	cout << fmt("// B^%*y is a solution\n") % (i-1);

	mpz_class middle = globals::modulus / 2;

	ofstream wf;
	must_open(wf, files::w % mine.scol);
	for(uint i = 0 ; i < nr ; i++) {
		if (v[i] > middle)
			v[i] -= globals::modulus;
		traits::print(wf,v[i]) << "\n";
	}
	wf.close();
}

};

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
	get_matrix_header(mtx, nr, mstr);

	globals::modulus = mpz_class(mstr);
	globals::modulus_u8	= globals::modulus.get_ui();
	globals::modulus_u16	= globals::modulus.get_ui();
	globals::modulus_u32	= globals::modulus.get_ui();
	globals::modulus_ulong	= globals::modulus.get_ui();

	/*
	if (SIZ(globals::modulus.get_mpz_t()) != MODULUS_SIZE) {
		cerr << fmt("ERROR: RECOMPILE WITH"
				" ``#define MODULUS_SIZE %''\n")
			% SIZ(globals::modulus.get_mpz_t());
		exit(1);
	}
	*/

	detect_mn(m, n);
	cout << fmt("// detected m = %\n") % m;
	cout << fmt("// detected n = %\n") % n;

	cout << "// counting coefficients\n" << flush;
	vector<std::streampos> foo(1);
	vector<uint> bar(1);
	count_matrix_coeffs(mtx, nr, foo, bar);
	
	nb_coeffs = bar[0];

	cout << fmt("// matrix has % coeffs\n") % nb_coeffs;
	mtx.close();

	cout.flush();
	cerr.flush();

	idx = new uint32_t[nr + nb_coeffs];
	val = new  int32_t[nr + nb_coeffs];

	{
		ifstream mtx;
		must_open(mtx, files::matrix);
		fill_matrix_data(mtx, 0, nb_coeffs, 0, nr, idx, val);
	}

	globals::nbys = mine.nbys;

#if 0
	if (MODULUS_BITS < 8 && nbys == 8) {
		cout << "// Using SSE-2 code\n";
		program<sse2_8words_traits > x;
		x.go();
	} else if (nbys == 1 && SIZ(globals::modulus.get_mpz_t()) == MODULUS_SIZE) {
		cout << "// Using code for " << MODULUS_SIZE << " words\n";
		program<typical_scalar_traits<MODULUS_SIZE> > x;
		x.go();
	} else if (nbys == 1) {
		/* This amounts to at least a 2x slowdown */
		cout << "// Using generic code\n";
		program<variable_scalar_traits> x;
		x.go();
	} else {
		cerr << "no available code\n";
		exit(1);
	}
#endif
	cout << "// Using generic code\n";
	program<variable_scalar_traits> x;
	x.go();

	delete[] idx;
	delete[] val;
}

/* vim:set sw=8: */
