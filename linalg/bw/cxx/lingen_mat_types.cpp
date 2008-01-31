#include "lingen_mat_types.hpp"
#include "bitstring.hpp"
#include <ostream>
#include <fstream>

using namespace std;

ostream& operator<<(ostream& o, polmat const& P)
{
    for(unsigned int j = 0 ; j < P.ncols ; j++) {
        for(unsigned int i = 0 ; i < P.nrows ; i++) {
            if (i) o << " ";
            write_hexstring(o, P.poly(i,j), P.deg(i,j) + 1);
        }
        o << "\n";
    }
    return o;
}

istream& operator>>(istream& is, polmat & P)
{
    for(unsigned int j = 0 ; j < P.ncols ; j++) {
        for(unsigned int i = 0 ; i < P.nrows ; i++) {
            is >> ws;
            read_hexstring(is, P.poly(i,j), P.ncoef);
        }
        P.setdeg(j);
    }
    return is;
}

