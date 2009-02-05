#include <cstdlib>
#include <cstdio>
#include "auxfuncs.h"

#include "matrix_header.hpp"
#include "parsing_tools.hpp"
#include "fmt.hpp"

#include <iostream>
#include <sstream>

using namespace std;

void put_matrix_header(ostream & mtx, unsigned int nr, unsigned int nc, const string & mstr)
{
    if (mstr == "2") {
        mtx << fmt("% %\n") % nr % nc;
    } else {
        mtx << fmt("// % ROWS % COLUMNS ; MODULUS %\n") % nr % nc % mstr;
    }
    mtx << flush;
}

void get_matrix_header(istream & mtx,
        unsigned int &nr,
        unsigned int &nc,
        string & mstr)
{
    bool ok = false;
    for (uint line=0 ; !mtx.eof() ; line++) {
        string buffer;

        getline(mtx, buffer);

        istringstream st(buffer);

        if (!startswith(buffer, "//")) {
            if (line == 0 && st >> nr >> nc) {
                mstr="2";
                cerr << "// Auto-detected cado format\n";
                ok = true;
            }
            /* Too late, the header should have come
             * before that */
            break;
        }

        st >> wanted("//")
            >> nr >> wanted("ROWS")
            >> nc >> wanted("COLUMNS")
            >> ws >> optional(';')
            >> wanted("MODULUS") >> mstr;

        if (st) {
            ok = true ; break;
        }
    }
    if (ok) {
        cout << "// Matrix:"
            << " " << nr << " rows"
            << " " << nc << " columns"
            << " modulus " << mstr << endl;
        return;
    }

    die("unable to get matrix header ; file is corrupted\n", 1);
}

/* vim: set sw=4 sta et: */
