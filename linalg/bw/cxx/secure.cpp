#include <cstdio>
#include "manu.h"

#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "files.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"

#include <boost/cstdint.hpp>
#include <iterator>
#include <vector>
#include <deque>
#include <algorithm>

#include <gmp.h>
#include <gmpxx.h>

common_arguments common;

using namespace std;

/* The purpose of thie program is to obtain a pair of small row vectors (with
 * small NON-ZERO coefficients) x and m such that m = x M.
 *
 * The algorithm relies on the fact that the rows of M are quite sparse.
 *
 * we start with two zero vectors.
 *
 * for each coefficient of x in order, we substitute a random small
 * integer in place of the existing 0. We verify that the added
 * contribution to the vector m entails no cancellation. We pick random
 * coefficients until everything is well.
 *
 * As randomness is really not a concern, we are happy with the sequence :
 * +1, -1, +2, -2, ...
 *
 *
 * very obviously this check is not very useful in characteristic two.
 */

void contribute_check_overflow(vector<int>& m,
		int v,
		const matrix_line& l)
{
	for(unsigned int k = 0 ; k < l.size() ; k++) {
		mpz_class z = m[l[k].first];
		z += v * l[k].second;
		if (!z.fits_sint_p()) {
			cerr << "reached overflow, sorry\n";
			exit(1);
		}
		m[l[k].first] = z.get_si();
	}
}

int modred(int x, int p)
{
	x += p/2;
	x %= p;
	x -= p/2;
	return x;
}

int main(int argc, char * argv[])
{
	ios_base::sync_with_stdio(false);
	cerr.tie(&cout);
	cout.rdbuf()->pubsetbuf(0,0);
	cerr.rdbuf()->pubsetbuf(0,0);

	no_arguments nothing;
	process_arguments(argc, argv, common, nothing);
	unsigned int nr, nc;
	string mstr;

	ifstream mtx;

	must_open(mtx, files::matrix);
	get_matrix_header(mtx, nr, nc, mstr);
	BUG_ON_MSG(nr != nc, "Matrix is not square\n");
	/* It's a pain to handle the non-square case. And hardly useful
	 */

	vector<int> x(nr, 0);
	vector<int> m(nc, 0);

	cout << "// reading matrix rows" << endl;
	istream_iterator<matrix_line> mit(mtx);
	unsigned int maxb = 2;
	for(unsigned int p = 0 ; mit != istream_iterator<matrix_line>() ; p++) {
		vector<bool> accept(2048, true);
		matrix_line l = *mit++;
		unsigned int b;
		for(unsigned int k = 0 ; k < l.size() ; k++) {
			unsigned int j = l[k].first;
			/* Which coefficient will cancel m[j] ? */
			if (m[j] == 0) continue;
			if (m[j] % l[k].second) continue;
			int forbidden = - m[j] / l[k].second;
			b = 2 * (ABS(forbidden) - 1) + (forbidden < 0);
			if (b >= accept.size())
				continue;
			accept[b] = false;
		}
		for(b = 0 ; b < accept.size() && !accept[b] ; b++);
		if (b == accept.size()) {
			/* perhaps increase accept.size() above in this
			 * case */
			cerr << fmt("// row % : can not obtain"
					" checking vectors\n") % p;
			exit(1);
		}
		int v;
		v = (1 - (b & 1) * 2) * ((b / 2) + 1);
		if (b > maxb) {
			cerr << fmt("// row % : chosen %\n") % p % v;
			maxb = b;
		}
		contribute_check_overflow(m, v, l);
		x[p] = v;
	}
	cout << "// ok" << endl;

	ofstream o;

	mpz_class px(mstr);
	if (px <= 32768) {
		int p = px.get_si();
		for(unsigned int i = 0 ; i < nr ; i++) {
			x[i] = modred(x[i], p);
		}
		for(unsigned int i = 0 ; i < nc ; i++) {
			m[i] = modred(m[i], p);
		}
	}

	must_open(o, files::x0);
	copy(x.begin(), x.end(), ostream_iterator<int>(o, "\n"));
	o.close();

	must_open(o, files::m0);
	copy(m.begin(), m.end(), ostream_iterator<int>(o, "\n"));
	o.close();

}
