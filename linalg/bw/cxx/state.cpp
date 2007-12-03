#include <cstdio>
#include <cstdlib>
#include "auxfuncs.h"
#include "manu.h"

#include <sys/types.h>
#include <unistd.h>
#include <gmp.h>
#include <gmpxx.h>
#include "gmp-hacks.h"

#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "slave_arguments.hpp"
#include "files.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"

#include <iterator>

using namespace std;

int recoverable_iteration(int m, int col, int cp_lag)
{
	mpz_class t;

	vector<int> per_row;
	for(int i = 0 ; i < m ; i++) {
		ifstream a((files::a % i % col).c_str());
		if (!a.is_open()) {
			per_row.push_back(0);
			continue;
		}
		int v;
		for(v = 0 ; a >> t ; v++);
		per_row.push_back(v);
	}
	BUG_ON(per_row.empty());

	int mi = *min_element(per_row.begin(), per_row.end());
	int mx = *max_element(per_row.begin(), per_row.end());
	
	/* Each checkpoint is (cp_lag) integers long ; when k integers
	 * are present in the A files, the last is (x, B^(k-1) y). We
	 * know that the A files are flushed once B^k y has been
	 * computed, for k a multiple of cp_lag. Hence we expect to
	 * discard at least one integer in all cases (except when there
	 * is no data anyway).
	 */
	mi -= (mi != 0);
	mx -= (mx != 0);

	mi -= mi % cp_lag;
	mx -= mx % cp_lag;

	if (mi != mx) {
		die("Found inconsistent amount of data (%d ... %d)\n",1,mi,mx);
	}

	return mi / cp_lag;
}

template<typename T> inline istream_iterator<T> beginof(ifstream& i)
{ return istream_iterator<T>(i); }
template<typename T> inline istream_iterator<T> endof(ifstream&)
{ return istream_iterator<T>(); }

int recover_vector(int nr, int col, int r, mp_limb_t (*w)[MODULUS_SIZE])
{
	ifstream v;
	if (r) {
		must_open(v, files::v % col % r);
	} else {
		must_open(v, files::y % col);
	}
	mpz_class z;
	istream_iterator<mpz_class> it(v);
	for( ; nr && it != endof<mpz_class>(v) ; it++, nr--, w++) {
		z = *it;
		BUG_ON(z < 0);
		BUG_ON(ABS(SIZ(z.get_mpz_t())) > MODULUS_SIZE);
		MPN_SET_MPZ(*w, MODULUS_SIZE, z.get_mpz_t());
	}
	BUG_ON(nr != 0);
	return 0;
}

/* We read the first r == (k * cp_lag) integers in the A files.  */
int recover_iteration(int m, int col, int r)
{
	for(int i = 0 ; i < m ; i++) {
		string a_string = files::a % i % col;
		const char * a_nm = a_string.c_str();
		std::streampos a_pos;
		if (r == 0) {
			a_pos = 0;
		} else {
			/* We must truncate the a files to the right size. */
			ifstream a;
			mpz_class dummy;
			must_open(a, a_nm);
			int j;
			for(j = 0 ; j < r && (a >> dummy); j++);
			BUG_ON(j != r);
			/* Important : get the newline as well ! */
			a >> ws;
			a_pos = a.tellg();
			a.close();
		}
		cout << fmt("// truncating % to size %\n") % a_nm % a_pos;
		truncate(a_nm, a_pos);
	}
	return 0;
}
