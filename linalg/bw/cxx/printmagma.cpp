#include "constants.hpp"
#include "matrix_header.hpp"
#include "matrix_line.hpp"
#include "must_open.hpp"
#include "parsing_tools.hpp"
#include "arguments.hpp"
#include "common_arguments.hpp"

#include <iostream>
#include <iterator>

using namespace std;
using namespace boost;

void pvec(const std::string& n, const vector<int>& co)
{
	cout << fmt("%:=Vector([GF(p) | \n") % n;
	copy(co.begin(), co.end() - 1,
			ostream_iterator<int>(cout, ", "));
	cout << co.back() << "]);\n";
}

bool rvec(const std::string& n, const std::string& f, vector<int>& co)
{
	ifstream y;
	if (!open(y, f)) return false;
	copy(istream_iterator<int>(y),
			istream_iterator<int>(),
			back_insert_iterator<vector<int> >(co));
	return true;
}

int main(int argc, char * argv[])
{
	ios_base::sync_with_stdio(false);
	cerr.tie(&cout);
	cout.rdbuf()->pubsetbuf(0,0);
	cerr.rdbuf()->pubsetbuf(0,0);

	common_arguments co;
	no_arguments nothing;

	process_arguments(argc, argv, co, nothing);

	unsigned int nr;
	unsigned int nc;

	do {
		ifstream mtx;
		string mstr;
		if (!open(mtx, files::matrix)) break;
		get_matrix_header(mtx, nr, nc, mstr);

		cout << fmt("p:=%;\n") % mstr;
		cout << fmt("M:=SparseMatrix(GF(p), %, %, [\n") % nr % nc;
		istream_iterator<matrix_line> mit(mtx);
		unsigned int p;
		vector<pair<pair<uint32_t, uint32_t>, int32_t> > coeffs;

		for(p = 0 ; mit != istream_iterator<matrix_line>() ; p++) {
			matrix_line l = *mit++;
			for(unsigned int k = 0 ; k < l.size() ; k++) {
				coeffs.push_back(
						make_pair(
							make_pair(
								(1+p),
								(1+l[k].first)),
							l[k].second));
			}
		}

		for(unsigned int k = 0 ; k < coeffs.size() - 1; k++) {
			cout << fmt("<%, %, %>, ")
					% coeffs[k].first.first
					% coeffs[k].first.second
					% coeffs[k].second;
		}
		cout << fmt("<%, %, %>")
			% coeffs.back().first.first
			% coeffs.back().first.second
			% coeffs.back().second;

		cout << "]);\n";
		cout << "M:=Matrix(M);\n";
	} while (0);

	unsigned int m = 0;
	for(int j = 0 ; ; j++, m++) {
		ifstream x;
		unsigned int c;
		if (!open(x, files::x % j)) break;
		vector<int> co(nr, 0);
		x >> wanted('e') >> c;
		co[c] = 1;
		pvec(fmt("X%")%j, co);
	}
	cout << "XX:=Matrix([";
	for(unsigned int j = 0 ; j < m-1 ; j++) cout << fmt("X%, ")%j;
	cout << fmt("X%]);\n")%(m-1);

	unsigned int n = 0;
	for(int j = 0 ; ; j++, n++) {
		vector<int> co;
		if (!rvec(fmt("Y%")%j, files::y % j, co)) break;
		pvec(fmt("Y%")%j, co);
	}
	cout << "YY:=Matrix([";
	for(unsigned int j = 0 ; j < n-1 ; j++) cout << fmt("Y%, ")%j;
	cout << fmt("Y%]);\n")%(n-1);

	do {
		vector<int> co;
		if (!rvec("CM0", files::m0, co)) break;
		pvec("CM0", co);
	} while (0);
	do {
		vector<int> co;
		if (!rvec("CX0", files::x0, co)) break;
		pvec("CX0", co);
	} while (0);
}
