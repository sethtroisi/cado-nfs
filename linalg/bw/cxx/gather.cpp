/* The norm requires that UINT16_MAX be defined only when this is on. It
 * is used in matrix_repr_binary_sliced.hpp */
#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <cstdio>
#include <cstdlib>
#include "auxfuncs.h"
#include "manu.h"
#include "gmp-hacks.h"

#include "addmul.hpp"
#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "config_file.hpp"
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
#include "preconditioner.hpp"
#include "matmul_toy.hpp"
#include "bitstring.hpp"

#include <stdint.h>
#include <iterator>

no_arguments mine;

using namespace std;
// using namespace boost;
using namespace core_ops;

typedef unsigned int uint;

namespace globals {
	unsigned int m, n;
	unsigned int col;
	unsigned int nr, nc;
	unsigned int nsols;
	mpz_class modulus;
	uint8_t modulus_u8;
	uint16_t modulus_u16;
	uint32_t modulus_u32;
	unsigned long modulus_ulong;
}

template<typename traits>
struct program_common: public simple_matmul_toy<traits>
{
	typedef typename traits::scalar_t scalar_t;
	typedef typename traits::wide_scalar_t wide_scalar_t;
	typename traits::representation::matrix_rowset mat;
        typedef simple_matmul_toy<traits> super;
	using super::v;
	using super::w;

        program_common() : super(files::matrix) {}
void add_file(const std::string& fn)
{
	ifstream f;
	must_open(f, fn);
	using namespace globals;

	for(unsigned int i = 0 ; i < nc ; i++) {
#if 0
		/* This really works only with one member, as the
		 * operation is not too much focused on vectoring things
		 * so far.
		 */
		std::vector<mpz_class> foo(1);
		traits::assign(foo, v[i]);
		std::string line;
		getline(f, line);
		istringstream ss(line);
		unsigned int ny;
		mpz_class blah;
		for(ny = 0 ; ss >> blah ; foo[0] += blah, ny++);
		if (nbys == 0 || ny != 0) {
			// cerr << "// auto-detected nbys=" << ny << "\n";
			nbys = ny;
		}
		if (ny == 0 && f.eof()) {
			cerr << fmt("// reached EOF after reading % coords from %\n")
				% i % fn;
			abort();
		}
		BUG_ON(ny != nbys);
		traits::assign(v[i], foo, 0);
#endif
                uint64_t x;
                read_hexstring(f, &x, 64);
                v[i].p ^= x;
	}
	f >> ws;
	BUG_ON(!f.eof());
}

void do_sum()
{
	using namespace globals;
        vector<string> file_names;
	if (file_names.empty()) {
		for(uint j = 0 ; j < m ; j++) {
			for(uint r = 1 ; ; r++) {
				ifstream fxy;
				string nm = files::fxy % j % r;
				if (!open(fxy, nm))
					break;
				file_names.push_back(nm);
				cout << "// selecting " << nm << endl;
			}
		}
	}

	BUG_ON(file_names.empty());

	typedef vector<string>::const_iterator vsci_t;
	for(vsci_t it = file_names.begin() ; 
			it != file_names.end() ;
			it++)
	{
		add_file(*it);
	}
}

void shipout(scalar_t * v, string const & s, bool check=true)
{
	using namespace globals;

	mpz_class middle = globals::modulus / 2;
	ofstream wf;
	must_open(wf, s);
	/* Last but not least, apply the preconditioner in order to get a
	 * true solution !!! */
	precond(v);
	if (check && is_zero(v)) {
		cerr << "// trivial solution found after preconditioning !\n";
		BUG();
	}
	for(uint i = 0 ; i < nr ; i++) {
        /*
		if (v[i] > middle)
			v[i] -= globals::modulus;
		traits::print(wf,v[i]) << "\n";
        */
                write_hexstring(wf, &(v[i].p), 64);
                wf << "\n";
	}
	wf.close();
}
};

template<typename traits>
struct program : public program_common<traits> {
	typedef program_common<traits> super;
	using super::v;
	using super::w;
public:
void go()
{
	using namespace globals;


	super::do_sum();

	traits::reduce(v, v, 0, nr);

	shipout(v, files::r, false);

	if (is_zero(v)) {
		cerr << "// trivial solution encountered !\n";
		BUG();
	}
	traits::copy(w, v, nc);

	int i;

	for(i = 0 ; i < 20 ; i++) {
		traits::copy(v, w, nc);
		cout << fmt("// trying B^%*y") % i;
		super::multiply();
		int r = super::nb_nonzero(w);
		if (r) {
			cout << fmt(" -- % non-zero") % r << endl;
		} else {
			cout << " -- zero !" << endl;
			break;
		}
	}

	if (i == 20 && !super::is_zero(w)) {
		cerr << "// no solution found\n";
		cerr << "// this may either be a bug, or an indication that\n";
		cerr << "// the span of x^T*B^i is too small. In the latter\n";
		cerr << "// case, recombining the vectors might work.\n";
		/* doing so automatically is bound to be a pain. */
		exit(1);
		// BUG();
	}

	cout << fmt("// B^%*y is a solution\n") % i;

	shipout(v, files::w);
}
};

template<typename traits>
class program_secondtry : public program_common<traits> {
public:
	typedef typename traits::scalar_t scalar_t;
	typedef typename traits::wide_scalar_t wide_scalar_t;
	typedef program_common<traits> super;
	scalar_t *zv;
	scalar_t *zw;
	using super::v;
	using super::w;
private:
void init_data()
{
	using namespace globals;
	zv = new scalar_t[nr];
	zw = new scalar_t[nr];

	traits::zero(zv, nr);
	traits::zero(zw, nr);
}

void clear_data()
{
	delete[] zv;
	delete[] zw;

}
public:
program_secondtry() { init_data(); }
~program_secondtry() { clear_data(); }
void go()
{
	using namespace globals;

	super::do_sum();

	traits::reduce(v, v, 0, nr);
	super::multiply();

	traits::copy(zw, w, nr);

	int i;

	for(i = 0 ; i < 20 ; i++) {
		traits::copy(v, w, nr);
		cout << fmt("// trying B^%*y") % i;
		super::multiply();

		/* find some point for pivoting */
		int r = super::nb_nonzero(w);
		if (r) {
			cout << fmt(" -- % non-zero") % r << endl;
		} else {
			cout << " -- zero !" << endl;
			break;
		}
	}

	if (i == 20 && !super::is_zero(w)) {
		cerr << "// no solution found\n";
		cerr << "// this may either be a bug, or an indication that\n";
		cerr << "// the span of x^T*B^i is too small. In the latter\n";
		cerr << "// case, recombining the vectors might work.\n";
		/* doing so automatically is bound to be a pain. */
		exit(1);
		// BUG();
	}

	cout << fmt("// B^%*y is a solution\n") % i;

	mpz_class middle = globals::modulus / 2;

	ofstream wf;
	must_open(wf, files::w);
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

    std::cout << "// This is bw-gather, version " VERSION << std::endl;

	common_arguments common;

	process_arguments(argc, argv, common, mine);

	string mstr;
	ifstream mtx;
	
	using namespace globals;

	must_open(mtx, files::matrix);
	get_matrix_header(mtx, nr, nc, mstr);
	mtx.close();

	BUG_ON_MSG(nr > nc, "Matrix has too many rows\n");

	globals::modulus = mpz_class(mstr);
        BUG_ON(globals::modulus != 2);

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

        config_file_t cf;
        read_config_file(cf, files::params);

        bool rc = true;
        rc &= get_config_value(cf, "m",m); cout << fmt("// m = %\n") % m;
        rc &= get_config_value(cf, "n",n); cout << fmt("// n = %\n") % n;

        if (!rc) {
                cerr << "Missing config file " << files::params << "\n";
                exit(1);
        }


	cout.flush();
	cerr.flush();

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
	program<binary_pod_traits<uint64_t> > x;
	x.go();
}

/* vim:set sw=8: */
