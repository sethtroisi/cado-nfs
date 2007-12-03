#include <cstdlib>
#include <cstdio>
#include "auxfuncs.h"

#include "matrix_header.hpp"
#include "parsing_tools.hpp"
#include "fmt.hpp"

#include <iostream>
#include <sstream>

using namespace std;

void put_matrix_header(ostream & mtx, unsigned int nr, const string & mstr)
{
	mtx << fmt("// % ROWS % COLUMNS ; MODULUS %\n") % nr % nr % mstr;
	mtx << flush;
}

void get_matrix_header(istream & mtx,
		unsigned int &nr,
		unsigned int &nc,
		string & mstr)
{
	for ( ; !mtx.eof() ; ) {
		string buffer;

		getline(mtx, buffer);

		if (!startswith(buffer, "//")) {
			/* Too late, the header should have come
			 * before that */
			break;
		}

		istringstream st(buffer);
		st >> wanted("//")
			>> nr >> wanted("ROWS")
			>> nc >> wanted("COLUMNS")
			>> ws >> optional(';')
			>> wanted("MODULUS") >> mstr;

		if (st) {
			cout << "// Matrix:"
				<< " " << nr << " rows"
				<< " " << nc << " columns"
				<< " modulus " << mstr << endl;
			return;
		}
	}

	die("unable to get matrix header ; file is corrupt\n", 1);
}

void get_matrix_header(istream & mtx, unsigned int &nr, string & mstr)
{
	unsigned int nc;
	get_matrix_header(mtx, nr, nc, mstr);
	if (nr != nc)
		die("matrix is not square\n", 1);
}
