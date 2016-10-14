#include "cado.h"
#include <gmp.h>
#include <istream>
#include <ostream>
#include <string>
#include <stdlib.h>
#include "gmpxx.hpp"

using namespace std;

/* In case we don't have gmpxx, we have to provide the C++ I/O functions
 * by ourselves...
 *
 * We make no effort to do this very accurately, though.
 */

ostream& operator<<(ostream& os, mpz_srcptr x)
{
    char * str;
    os << (str = mpz_get_str(NULL, 10, x));
    free(str);
    return os;
}

ostream& operator<<(ostream& os, mpq_srcptr x)
{
    char * str;
    os << (str = mpq_get_str(NULL, 10, x));
    free(str);
    return os;
}

istream& operator>>(istream& is, mpz_ptr x)
{
    string s;
    is >> s;
    mpz_set_str(x, s.c_str(), 0);
    return is;
}


istream& operator>>(istream& is, mpq_ptr x)
{
    string s;
    is >> s;
    mpq_set_str(x, s.c_str(), 0);
    return is;
}

