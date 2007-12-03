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

#include <boost/cstdint.hpp>
#include <iterator>

gather_arguments mine;

using namespace std;
using namespace boost;
using namespace core_ops;

typedef unsigned int uint;

namespace globals {
	uint m, n;
	uint col;
	uint nr;
	mpz_class modulus;

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

	memset(v, 0, sizeof(scalar_t[nr]));
	memset(w, 0, sizeof(scalar_t[nr]));
}

void clear_data()
{
	using namespace globals;
	delete[] v;
	delete[] w;
	delete[] scrap;
	delete[] idx;
	delete[] val;

}

void add_file(const std::string& fn)
{
	ifstream f;
	must_open(f, fn);
	using namespace globals;

	istream_iterator<mpz_class> it(f);

	for(unsigned int i = 0 ; i < nr ; i++) {
		mpz_class foo;
		traits::assign(foo, v[i]);
		foo += *it++;
		traits::assign(v[i], foo);
	}
	BUG_ON(it != istream_iterator<mpz_class>());
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
		if (!traits::is_zero(v[i])) {
			return false;
		}
	}
	return true;
}

void multiply()
{
	using namespace globals;
	const uint32_t * ip = idx;
	const int32_t  * vp = val;
	for(uint i = 0 ; i < nr ; i++) {
		traits::zero(scrap[i]);
		uint c = 0;
		for( ; *vp != 0 ; ip++, vp++) {
			c += *ip;
			traits::addmul(scrap[i], v[c], *vp);
		}
		ip++;
		vp++;
	}
	traits::reduce(w, scrap, 0, nr);
}

program() { init_data(); }
~program() { clear_data(); }

void go()
{
	using namespace globals;

	init_data();
	do_sum();

	if (is_zero(v)) {
		cerr << "// trivial solution encountered !\n";
		BUG();
	}
	memcpy(w, v, nr * sizeof(scalar_t));

	uint i;

	for(i = 0 ; i < 100 && !is_zero(w) ; i++) {
		memcpy(v, w, nr * sizeof(scalar_t));
		cout << fmt("// trying B^%*y") % i << endl;
		multiply();
	}

	if (i == 100 && !is_zero(w)) {
		cerr << "// no solution found\n";
		BUG();
	}

	cout << fmt("// B^%*y is a solution\n") % (i-1);

	ofstream w;
	must_open(w, files::w % mine.scol);
	for(uint i = 0 ; i < nr ; i++) {
		traits::print(w,v[i]) << "\n";
	}

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

	if (SIZ(globals::modulus.get_mpz_t()) != MODULUS_SIZE) {
		cerr << fmt("ERROR: RECOMPILE WITH"
				" ``#define MODULUS_SIZE %''\n")
			% SIZ(globals::modulus.get_mpz_t());
		exit(1);
	}

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

	program<typical_scalar_traits<MODULUS_SIZE> > x;
	x.go();
}

/* vim:set sw=8: */
