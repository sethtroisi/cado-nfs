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

const int extra = 64;

void reduce_mline(matrix_line& v, boost::int32_t p)
{
	matrix_line::iterator dst;
	matrix_line::const_iterator src;
	dst = v.begin();
	for(src = v.begin() ; src != v.end() ; src++) {
		boost::int32_t pm;
		pm = src->second % p;
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

/* Note that it is important to do the clamping in the most trivial
 * manner here. If we do it in a smarter way, then we'll have to face
 * true life later on, and abandon the quirky matrix.txt.ur trick.
 */

/* columns[i]'s are ideal appearance vectors. */
void clamp(vector<matrix_line>& columns)
{
	vector<matrix_line>::iterator dst = columns.begin();
	for(dst = columns.begin() ; dst != columns.end() ; dst++) {
		/* We should really do it more carefully */
		vector<matrix_line::value_type>::iterator q = dst->begin();
		vector<matrix_line::value_type>::const_iterator p = dst->begin();
		for( ; p != dst->end() ; p++) {
			if (p->first < columns.size())
				*q++ = *p;
		}
		dst->erase(q, dst->end());
	}
}

int prune_empty_rows(vector<matrix_line>& columns, vector<int>& ideal_idx)
{
	vector<matrix_line>::iterator dst = columns.begin();
	vector<matrix_line>::const_iterator src = columns.begin();
	vector<int>::iterator idst = ideal_idx.begin();
	vector<int>::const_iterator isrc = ideal_idx.begin();
	assert(columns.size() == ideal_idx.size());
	for( ; src != columns.end() ; src++, isrc++) {
		if (!src->empty()) {
			*dst++=*src;
			*idst++=*isrc;
		}
	}
	int r = (columns.end()-dst);
	cout << "Erasing " << r << " zero rows\n";

	columns.erase(dst, columns.end());
	ideal_idx.erase(idst, ideal_idx.end());

	columns.insert(columns.end(),extra,matrix_line());
	ideal_idx.insert(ideal_idx.end(),extra,-1);
	return r;
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
		string s;
		ostringstream sstr;
		getline(mtx, s);
		if (mtx.eof())
			break;
		unsigned int x = s.find('|');
		istringstream mt(s.substr(x+2));
		
		if (!(mt >> li)) {
			break;
		}

		sstr << s.substr(0,x) << endl;

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
	 * as nindexes < the number of rows selected. */

	if (rows.size() < nindexes + extra) {
		cerr << "Cannot ensure a dependency\n";
		exit(0);
	}

	std::sort(sizes.begin(), sizes.end());

	/* We need to have as many as nindexes + extra ideal relations,
	 * so that it exceeds the # of FB ideals (nindexes) by at least
	 * extra
	 */
	unsigned int discard_weight = sizes[nindexes + extra] + 1;
	{
		vector<matrix_line>::iterator dst;
		vector<matrix_line>::const_iterator src;
		vector<string>::iterator c_dst;
		vector<string>::const_iterator c_src;
		c_src = clink_data.begin();
		c_dst = clink_data.begin();
		dst = rows.begin();
		unsigned int nselected = 0;
		for(src = rows.begin() ; src != rows.end() ; src++) {
			if (src->size() < discard_weight) {
				*dst++ = *src;
				*c_dst++ = *c_src;
				// clink << *c_iter << "\n";
				if (++nselected >= nindexes + extra)
					break;
			}
			c_src++;
		}
		cout << (dst - rows.begin()) << " rows selected\n";
		rows.erase(dst, rows.end());
		clink_data.erase(c_dst, clink_data.end());
	}

	vector<matrix_line> columns;
	vector<int> ideal_idx(rows.size());
	columns.assign(rows.size(), matrix_line());
	for(unsigned int j = 0 ; j < rows.size() ; j++) {
		ideal_idx[j]=j;
	}
	unsigned int total_coeffs = 0;
	{
		vector<matrix_line>::const_iterator src;
		for(src = rows.begin() ; src != rows.end() ; src++) {
			int i = src - rows.begin();
			for(unsigned int j = 0 ; j < src->size() ; j++) {
				unsigned int jj = (*src)[j].first;
				int v = (*src)[j].second;
				columns[jj].push_back(make_pair(i, v));
			}
			total_coeffs += src->size();
		}
	}

	cout << "total: "  << total_coeffs << " coefficients\n";

	/* We might choose some very heavy ideal rows (those which were
	 * columns before -- we've now translated the thing) to become
	 * character rows. It means that we'll no longer have them in the
	 * matrix.
	 */
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
	/* First write the unreduced version of the matrix. We do so
	 * BEFORE pruning part of the matrix, because otherwise the magma
	 * code will get nuts.
	 */ 
	must_open(cfile, mine.out + ".ur");
	put_matrix_header(cfile, columns.size(), columns.size(), mstr.str());
	for(unsigned int j = 0 ; j < columns.size() ; j++) {
		cfile << columns[j] << "\n";

	}
	cfile.close();

#if 1
	for(;;) {
		int r=prune_empty_rows(columns, ideal_idx);
		if (r == extra)
			break;
		clamp(columns);
	}
#endif

	/* The clink info should be linked to pruning/clampind better
	 * than it is now */
	ofstream clink;
	must_open(clink, mine.clink);
	for(unsigned int i = 0 ; i < columns.size() ; i++) {
		clink << clink_data[i] << "\n";
	}
	clink.close();

	ofstream rlink;
	must_open(rlink, mine.rlink);
	for(unsigned int i = 0 ; i < ideal_idx.size() ; i++) {
		rlink << ideal_idx[i] << "\n";
	}
	rlink.close();


	must_open(cfile, mine.out);
	put_matrix_header(cfile, columns.size(), columns.size(), mstr.str());

	if (mine.p < 1000000) {
		for(unsigned int j = 0 ; j < columns.size() ; j++) {
			reduce_mline(columns[j], mine.p.get_si());
			cfile << columns[j] << "\n";
		}
	} else {
		for(unsigned int j = 0 ; j < columns.size() ; j++) {
			cfile << columns[j] << "\n";
		}
	}
	cfile.close();

	return 0;
}
