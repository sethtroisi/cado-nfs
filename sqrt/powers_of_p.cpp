
/* This file uses the standard C++ functions to provide a look-up table
 * for primes of a power p. It must be given an opaque structure, which
 * in fact is of type power_lookup_table
 */
#include <cstdlib>
#include <cstdio>
#include <map>
#include <gmp.h>
#include "powers_of_p.h"

using namespace std;

struct power_lookup_table {
    unsigned long p;
    typedef map<int, int> m_t;
    mpz_t * z;
    int alloc;
    int nz;
    m_t m;
    int extra_power_swapstore(mpz_t w) {
        if (nz == alloc) {
            alloc += alloc ? alloc / 4 : 16;
            z = (mpz_t *) realloc(z, alloc * sizeof(mpz_t));
        }
        mpz_init(z[nz]);
        mpz_swap(z[nz], w);
        nz++;
        return nz-1;
    }
    power_lookup_table(unsigned long p) : p(p), z(NULL), alloc(0), nz(0) {}
    ~power_lookup_table() {
        for(int i = 0 ; i < nz ; i++) {
            mpz_clear(z[i]);
        }
        free(z);
        z = NULL;
    }
    mpz_srcptr operator()(int i);
};

mpz_srcptr power_lookup_table::operator()(int i)
{
    m_t::const_iterator px = m.find(i);
    if (px != m.end()) return z[px->second];
    mpz_t q;
    // XXX valgrind says that this sometimes leaks. I don't understand
    // why.  Perhaps it's obvious.
    mpz_init(q);
    if (i == 0) {
        mpz_set_ui(q, 1);
    } else if (i == 1) {
        mpz_init_set_ui(q, p);
    } else if (i & 1) {
        mpz_srcptr pj = (*this)(i-1);
        mpz_mul_ui(q,pj,p);
    } else {
        mpz_srcptr pj = (*this)(i/2);
        mpz_mul(q,pj,pj);
    }
    int r = extra_power_swapstore(q);
    m[i]=r;
    mpz_clear(q);  // has been swapped with the other one.
    return z[r];
}

// C entry functions

void * power_lookup_table_init(unsigned long p)
{
    return new power_lookup_table(p);
}

void power_lookup_table_clear(void * t)
{
    power_lookup_table * pt = static_cast<power_lookup_table *>(t);
    return delete pt;
}

mpz_srcptr power_lookup(void * t, int i)
{
    power_lookup_table * pt = static_cast<power_lookup_table *>(t);
    return (*pt)(i);
}
