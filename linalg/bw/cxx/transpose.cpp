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
#include "transpose_arguments.hpp"

#include <boost/cstdint.hpp>
#include <iterator>
#include <vector>
#include <deque>
#include <algorithm>

common_arguments common;
transpose_arguments mine;

using namespace std;

void doit(std::ostream& o)
{
	unsigned int nr;
	unsigned int nc;
	string mstr;

	/* The transpose matrix */
	vector<matrix_line> columns;

	ifstream mtx;
	
	must_open(mtx, files::matrix);
	get_matrix_header(mtx, nr, nc, mstr);
	mtx.close();

	/*
	vector<streampos>  pos(1, 0);
	vector<uint32_t>  ncoeffs(1, 0);

	count_matrix_coeffs(mtx, nr, mtxfile_pos, nb_coeffs);

	uint32_t * idx = new uint32_t[i1 - i0 + nb_coeffs[t]];
	int32_t *  val = new int32_t[i1 - i0 + nb_coeffs[t]];

	fill_matrix_data(mtx, 0, nb_coeffs[t],
				0, nr, idx, val);

				*/

	columns.assign(nc, matrix_line());

	must_open(mtx, files::matrix);
	// istream_iterator<matrix_line> mit(mtx);

	unsigned int max_rows_kept = min(nc, nr);

	for(unsigned int p = 0 ; p < max_rows_kept ; p++) {
		matrix_line li;
		if (!(mtx >> li)) {
			break;
		}

		unsigned int j0 = 0;
		unsigned int j = 0;
#if 0
		if (li[j].first == 0) {
			li[j].second +=  3;
			if (li[j].second) li[j0++] = li[j];
			j++;
		}
		if (li[j].first == 1) {
			li[j].second +=  4;
			if (li[j].second) li[j0++] = li[j];
			j++;
		}
		if (li[j].first == 2) {
			li[j].second +=  1;
			if (li[j].second) li[j0++] = li[j];
			j++;
		}
		if (li[j].first == 3) {
			li[j].second +=  1;
			if (li[j].second) li[j0++] = li[j];
			j++;
		}
		if (li[j].first == 4) {
			li[j].second +=  -1;
			if (li[j].second) li[j0++] = li[j];
			j++;
		}
		if (li[j].first == 5) {
			li[j].second +=  1;
			if (li[j].second) li[j0++] = li[j];
			j++;
		}
#endif
		for( ; j < li.size() ; li[j0++] = li[j++]);

		for(j = 0 ; j < j0 ; j++) {
			unsigned int jj = li[j].first;
			li[j].first = p;
			columns[jj].push_back(li[j]);
		}
	}
	put_matrix_header(o, nc, mstr);
	if (mstr != string("2")) {
		for(unsigned int j = 0 ; j < nc ; j++) {
			o << columns[j] << "\n";
		}
	} else {
		for(unsigned int j = 0 ; j < nc ; j++) {
			print_line_without_ones(o, columns[j]) << "\n";
		}
	}
}

int main(int argc, char * argv[])
{
	ios_base::sync_with_stdio(false);
	cerr.tie(&cout);
	cout.rdbuf()->pubsetbuf(0,0);
	cerr.rdbuf()->pubsetbuf(0,0);

	process_arguments(argc, argv, common, mine);

	if (!mine.out.empty()) {
		ofstream ofile;
		must_open(ofile, mine.out);
		doit(ofile);
	} else {
		doit(cout);
	}

	return 0;
}
