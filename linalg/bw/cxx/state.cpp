#include <cstdio>
#include <cstdlib>
#include "auxfuncs.h"
#include "manu.h"

#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <gmp.h>
#include <gmpxx.h>
#include "gmp-hacks.h"

#include "arguments.hpp"
#include "common_arguments.hpp"
#include "constants.hpp"
#include "krylov_arguments.hpp"
#include "files.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"

#include <iterator>
#include <string>
#include <sstream>

using namespace std;

int recoverable_iteration(int m, int col, int nbys, int cp_lag)
{
	mpz_class t;

	vector<int> per_row;
	for(int i = 0 ; i < m ; i++) {
		std::string a_nm = files::a % i % col;
		ifstream a(a_nm.c_str());
		if (!a.is_open()) {
			per_row.push_back(0);
			continue;
		}
		int v = 0;
		std::string line;
		for( ; getline(a, line) ; ) {
			istringstream ss(line);
			int ny;
			for(ny = 0 ; ss >> t ; ny++);
			if (ny != nbys) {
				/* The data file might not have been
				 * flushed */
				BUG_ON(!a.eof());
				break;
			}
			v++;
		}
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

/* We read the first r == (k * cp_lag) integers in the A files.  */
int recover_iteration(int m, int col, int nbys, int r)
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
			for(j = 0 ; j < r * nbys && (a >> dummy); j++);
			BUG_ON(j != r * nbys);
			/* Important : get the newline as well ! */
			a >> ws;
			a_pos = a.tellg();
			a.close();
		}
		struct stat s[1];
		if (stat(a_nm, s) == 0) {
			cout << fmt("// truncating % to size %\n")
				% a_nm % a_pos;
			truncate(a_nm, a_pos);
		} else if (errno != ENOENT) {
			cout << fmt("// % ->> %\n")
				% a_nm % strerror(errno);
		}
	}
	return 0;
}
