#include <cstdio>
#include "manu.h"

#include "arguments.hpp"
#include "common_arguments.hpp"
#include "nfsforge_arguments.hpp"
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
#include <string>
#include <sstream>
#include <deque>

/* This program is specially intended for use with the nfsforge code.
 *
 * It takes a file with matrix lines. Each matrix line is prepended with
 * some column link information, in the form of three integers (the first
 * two are nonnegative).
 *
 * A modulus is also expected in argument.
 *
 * The output of the program consists of two files. First, a matrix.txt
 * file with the *transposed* matrix, and coefficients reduced modulo the
 * given modulus.
 *
 * Also, a column_links.txt file is produced with the informations that
 * were attached to lines (because of transposing, they're attached to
 * columns now).
 *
 * The matrix is reduced so that heavy columns are discarded.
 *
 * The returned matrix is square, and is guaranteed to have a non-zero
 * right nullspace. If needed, zero rows are added.
 */
common_arguments common;
nfsforge_arguments mine;

using namespace std;

void reduce_mline(matrix_line& v, boost::int32_t p)
{
	matrix_line::iterator dst;
	matrix_line::const_iterator src;
	dst = v.begin();
	for(src = v.begin() ; src != v.end() ; src++) {
		boost::int32_t pm;
		pm = src->second % mine.p;
		if (pm < -(p/2)) {
			pm+=p;
		} else if (pm > p/2) {
			pm-=p;
		}
		if (pm) {
			*dst++ = make_pair(src->first, pm);
		}
	}
	v.erase(dst, v.end());
}

int main(int argc, char * argv[])
{
	ios_base::sync_with_stdio(false);
	cerr.tie(&cout);
	cout.rdbuf()->pubsetbuf(0,0);
	cerr.rdbuf()->pubsetbuf(0,0);

	process_arguments(argc, argv, common, mine);

	vector<matrix_line> rows;
	vector<unsigned int> sizes;
	vector<string> clink_data;
	ifstream mtx;
	
	must_open(mtx, mine.in);

	unsigned int nindexes = 0UL;

	for( ;; ) {
		matrix_line li;
		string q,r,x;
		ostringstream sstr;

		if (!(mtx >> q >> r >> x >> li)) {
			break;
		}

		matrix_line::const_iterator src;
		matrix_line::iterator dst;
		dst = li.begin();
		for(src = li.begin() ; src != li.end() ; src++) {
			*dst++ = *src;
		}

		li.erase(dst, li.end());

		sstr << q << ' ' << r << ' ' << x << flush;

		clink_data.push_back(sstr.str());
		rows.push_back(li);
		sizes.push_back(li.size());

		if (li.size()) {
			unsigned int nindexes_here = li.back().first + 1;
			if (nindexes_here > nindexes) {
				nindexes = nindexes_here;
			}
		}
	}
	mtx.close();

	cout << "read " << rows.size() << " rows\n";
	cout << "seen " << nindexes << " columns\n";

	/* We have nindexes columns. A dependency can be obtained as soon
	 * as nindexes < the number of rows selected
	 */

	if (rows.size() <= nindexes) {
		cerr << "Cannot ensure a dependency\n";
		exit(0);
	}

	std::sort(sizes.begin(), sizes.end());

	ofstream clink;
	must_open(clink, mine.clink);
	unsigned int discard_weight = sizes[nindexes + 1] + 1;
	{
		vector<matrix_line>::iterator dst;
		vector<matrix_line>::const_iterator src;
		vector<string>::const_iterator c_iter;
		c_iter = clink_data.begin();
		dst = rows.begin();
		unsigned int nselected = 0;
		for(src = rows.begin() ; src != rows.end() ; src++) {
			if (src->size() < discard_weight) {
				*dst++ = *src;
				clink << *c_iter << "\n";
				if (++nselected >= nindexes + 1)
					break;
			}
			c_iter++;
		}
		cout << (dst - rows.begin()) << " rows selected\n";
		rows.erase(dst, rows.end());
	}
	clink.close();

	vector<matrix_line> columns;
	columns.assign(rows.size(), matrix_line());
	vector<matrix_line>::const_iterator src;
	unsigned int total_coeffs = 0;
	for(src = rows.begin() ; src != rows.end() ; src++) {
		for(unsigned int j = 0 ; j < src->size() ; j++) {
			unsigned int jj = (*src)[j].first;
			columns[jj].push_back(make_pair(src - rows.begin(), (*src)[j].second));
		}
		total_coeffs += src->size();
	}

	cout << "total: "  << total_coeffs << " coefficients\n";
	if (mine.nchar) {
		/* Look for good character row candidates */
		/* (yes, these are rows by now) */

		deque<unsigned int> heaviest_rows(mine.nchar, 0);
		for(unsigned int j = 0 ; j < columns.size() ; j++) {
			unsigned int w = columns[j].size();
			if (w > heaviest_rows.front()) {
				heaviest_rows.push_back(w);
				heaviest_rows.pop_front();
				sort(heaviest_rows.begin(), heaviest_rows.end());
			}
		}
		cout << "weight of " << mine.nchar << " heaviest rows:";
		for(int j = 0 ; j < mine.nchar ; j++) {
			cout << " " << heaviest_rows[j];
		}
		cout << "\n";

		for(unsigned int j = 0 ; j < columns.size() && mine.nchar; j++) {
			unsigned int w = columns[j].size();
			if (w >= heaviest_rows.front()) {
				cout << "Discarding row " << j << "\n";
				total_coeffs -= columns[j].size();
				columns[j] = matrix_line();
				mine.nchar--;
			}
		}
		cout << "now: "  << total_coeffs << " coefficients\n";
	}

	ostringstream mstr;
	mstr << mine.p;

	ofstream cfile;

	/* First write the unreduced version of the matrix */
	must_open(cfile, mine.out + ".ur");
	put_matrix_header(cfile, rows.size(), mstr.str());
	for(unsigned int j = 0 ; j < columns.size() ; j++) {
		cfile << columns[j] << "\n";
	}
	cfile.close();

	must_open(cfile, mine.out);
	put_matrix_header(cfile, rows.size(), mstr.str());
	for(unsigned int j = 0 ; j < columns.size() ; j++) {
		reduce_mline(columns[j], mine.p);
		cfile << columns[j] << "\n";
	}
	cfile.close();

	return 0;
}
