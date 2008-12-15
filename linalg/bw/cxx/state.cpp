#include <cstdio>
#include <cstdlib>
#include "auxfuncs.h"
#include "manu.h"

#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <ctype.h>
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
#include <algorithm>

using namespace std;

/* There used to be a fair amount of mess with checkpoints, earlier on.
 * It seems actually that the situation is much easier than what one
 * would be led to think by re-reading the old commit logs.
 *
 * Let K denote the checkpoint periodicity (variable cp_lag, also
 * reported on stdout). Let V_0 denote the vector Y, and V_i = M^(K*i)Y.
 *
 * The files A-*-* must be interpreted as sets of K lines. The i-th such
 * set (starting with i=0) contains X^T V_i, X^T M V_i, ..., X^T M^(K-1) V_i.
 * If this set is complete, and if V_{i+1} has been computed already,
 * then it is possible to continue.
 */
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
                    size_t len = line.size();
                    for( ; len && isspace(line[len-1]) ; ) ;
                    int ny = 4 * len;
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
	
	mi -= mi % cp_lag;
	mx -= mx % cp_lag;

	if (mi != mx) {
		die("Found inconsistent amount of data (%d ... %d)\n",1,mi,mx);
	}

	return mi / cp_lag;
}

/* We read the first r == (k * cp_lag) lines in the A files.  */
int recover_iteration(int m, int col, int nbys MAYBE_UNUSED, int r)
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
			int j = 0;
                        std::string line;
                        for(std::string line ; j < r && getline(a, line) ; j++);
			BUG_ON(j != r);
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
