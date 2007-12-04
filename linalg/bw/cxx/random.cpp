#include "constants.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"

#include <cstdlib>
#include <cmath>
#include <iostream>

#include <string>
#include <sstream>

#include <set>

#include <gmp.h>
#include <gmpxx.h>

using namespace std;

int random_coeff()
{
	int b = random() % 256;
	int v = (1 - (b & 1) * 2) * ((b / 2) + 1);
	return v;
}

struct rcoeff_large {
	int32_t operator()() const {
		return random_coeff();
	}
};

struct rcoeff_small {
	int p;
	rcoeff_small(mpz_class const& x) : p(x.get_si()) {}
	int32_t operator()() const {
		int v = random_coeff() % p;
		if (v < (-p/2)) {
			v += p;
		} else if (v > p/2) {
			v -= p;
		}
		return v;
	}
};

struct rcoeff_binary { int32_t operator()() const { return 1; } };

class fill_row_binomial {
	unsigned int nc;
	double p;
public:
	fill_row_binomial(unsigned int nc, double p) :
		nc(nc), p(p) {}
template<class T>
void operator()(matrix_line & l, T const& op) const
{
	l.clear();
	for(unsigned int j = 0 ; j < nc ; j++) {
		double t = random() / (double) RAND_MAX;
		if (t > p) continue;
		int32_t v = op();
		if (v == 0) continue;
		l.push_back(make_pair(j, v));
	}
}
};

class fill_row_normal {
	unsigned int nc;
	double mean, sigma;
	double pick() const
	{
		double a = random() / (double) RAND_MAX;
		double b = random() / (double) RAND_MAX;
		double z = sqrt(-2*log(a)) * cos(2*M_PI*b);
		return abs(mean + sigma * z);
	}
public:
	fill_row_normal(unsigned int nc, double p) :
		nc(nc),
		mean(nc * p),
		sigma(sqrt(nc*p*(1-p))) {}

template<class T>
void operator()(matrix_line & l, T const& op) const
{
	l.clear();
	unsigned int n = (unsigned int) pick();
	set<uint32_t> indices;
	for(unsigned int i = 0 ; i < n ; i++) {
		double t = random() / (double) RAND_MAX;
		indices.insert((uint32_t) (t*nc));
	}
	for(set<uint32_t>::const_iterator j = indices.begin() ; j != indices.end() ; j++) {
		int32_t v = op();
		if (v == 0) continue;
		l.push_back(make_pair(*j, v));
	}
}
};

struct fillmat {
	template<typename filler>
	void operator()(unsigned int nr, unsigned int nz,
			mpz_class const& px, filler const& f) const
	{
	if (px > 128) {
		for(unsigned int i = nz ; i < nr ; i++) {
			matrix_line l;
			f(l, rcoeff_large());
			cout << l << "\n";
		}
	} else if (px != 2) {
		rcoeff_small op(px);
		for(unsigned int i = nz ; i < nr ; i++) {
			matrix_line l;
			f(l, op);
			cout << l << "\n";
		}
	} else {
		for(unsigned int i = nz ; i < nr ; i++) {
			matrix_line l;
			f(l, rcoeff_binary());
			print_line_without_ones(cout,l) << "\n";
		}
	}
	}
};

int main(int argc, char * argv[])
{
	ios_base::sync_with_stdio(false);
	cerr.tie(&cout);
	cout.rdbuf()->pubsetbuf(0,0);
	cerr.rdbuf()->pubsetbuf(0,0);

	unsigned int seed = time(NULL);

	if (argc >= 4 && strcmp(argv[1], "--seed") == 0) {
		seed = atoi(argv[2]);
		argv++,argc--;
		argv++,argc--;
	}
	if (argc != 4 && argc != 5) {
		cerr << "Usage : ./random [--seed <seed>] <nrows> <ncols> <modulus> [<density>]\n";
		exit(1);
	}
	srandom(seed);

	unsigned int nr = atoi(argv[1]);
	unsigned int nc = atoi(argv[2]);
	string mstr(argv[3]);
	double p = 0.5;
	if (argc > 3) {
		istringstream foo(argv[4]);
		foo >> p;
	}
	if (p > 1) { p = p/nc; }

	mpz_class px(mstr);

	put_matrix_header(cout, nr, nc, mstr);

	/* Put a zero for the first row */
	cout << 0 << "\n";
	unsigned int nz = 1;

	if (nr > 100) {
		/* We should check nc*p and nc*(1-p) as well if we were
		 * concerned by good approximation, but it's not really
		 * the case */
		cerr << "Using normal approximation\n";
		fillmat()(nr, nz, px, fill_row_normal(nc, p));
	} else {
		fillmat()(nr, nz, px, fill_row_binomial(nc, p));
	}

	return 0;
}

