#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <set>
#include <vector>
#include <fstream>
#include <iosfwd>
#include <gmp.h>
#include <gmpxx.h>
#include <stdint.h>

struct matrix_slice {
    unsigned int i0;
    unsigned int i1;
    uint32_t ncoeffs;
    std::streampos pos;
};

class matrix_stats {
    public:
    unsigned int nr;
    unsigned int nc;
    mpz_class modulus;
    private:
    std::vector<matrix_slice> * slices;
    std::set<uint32_t> * zrows;
    std::set<uint32_t> * zcols;
    public:
    matrix_stats() : slices(NULL), zrows(NULL), zcols(NULL) {}
    void need_slices(std::vector<matrix_slice> * ptr, unsigned int nt = 0);
    void need_zcols(std::set<uint32_t> * ptr);
    void need_zrows(std::set<uint32_t> * ptr);
    void operator()(std::ifstream& mtx);
    void operator()(std::string const & name);
};

#endif	/* MATRIX_HPP_ */
