#include <iostream>
#include <sstream>
#include <string>
#include "must_open.hpp"
#include "matrix.hpp"
#include "manu.h"
#include "matrix_line.hpp"
#include "matrix_header.hpp"
#include "parsing_tools.hpp"
#include "manu.h"

using namespace std;


void matrix_stats::need_slices(vector<matrix_slice> * ptr, unsigned int nt) {
    slices = ptr;
    BUG_ON(nt && slices->size());
    BUG_ON(!nt && slices->empty());
    if (nt) {
        slices->assign(nt, matrix_slice());
    } else {
        nt = slices->size();
    }
}
/* Obtaining the information about zero columns seems pointless at
 * first. However it's useful for tests, and prevents silly mistakes
 * that one might make. And it's kinda cheap anyway.
 */
void matrix_stats::need_zcols(set<uint32_t> * ptr) {
    zcols = ptr;
}
void matrix_stats::need_zrows(set<uint32_t> * ptr) {
    zrows = ptr;
}
void matrix_stats::operator()(std::ifstream& mtx) {
    std::string mstr;
    get_matrix_header(mtx, nr, nc, mstr);
    std::istringstream ss(mstr);
    if (!(ss >> modulus >> ws) || !ss.eof()) {
        die("Garbage while reading modulus string: <<%s>>", 1, mstr.c_str());
    }

    if (!slices && !zrows && !zcols)
        return;

    unsigned int nhslices = slices ? slices->size() : 1;

    std::streampos pos = mtx.tellg();
    comment_strip cs(mtx, "//");

    vector<bool> colmap;
    if (zcols) {
        colmap.assign(nc,false);
    }
    // We'll constitute an horizontal slice info anyway, even if
    // there's only one such.
    for(uint j = 0 ; j < nhslices ; j++) {
        matrix_slice slice;
        slice.i0 = (j * nr) / nhslices;
        slice.i1 = ((j + 1) * nr) / nhslices;
        slice.ncoeffs = 0;
        slice.pos = mtx.tellg();
        for(uint i = slice.i0 ; i < slice.i1 ; i++) {
            std::string s;
            cs.getline(s);
            std::istringstream st(s);
            uint z;
            matrix_line mi;
            if (zcols) {
                if (!(st >> mi)) { BUG(); }
                z = mi.size();
                for(uint k = 0 ; k < mi.size() ; k++) {
                    colmap[mi[k].first] = true;
                }
            } else {
                if (!(st >> z)) { BUG(); }
            }
            if (zrows && !z) {
                zrows->insert(i);
            }
            slice.ncoeffs += z;
        }
        if (slices) {
            (*slices)[j] = slice;
            std::cout << fmt("// slice % : % rows (%..%) % coeffs from pos %\n")
                % j
                % (slice.i1 - slice.i0)
                % slice.i0
                % slice.i1
                % slice.ncoeffs
                % slice.pos;
            std::cout << std::flush;
        }
    }
    typedef set<uint32_t>::const_iterator suci_t;
    if (zrows) {
        std::cout << "// " << zrows->size() << " zero rows:";
        for(suci_t xi = zrows->begin() ; xi != zrows->end() ; xi++) {
            std::cout << " " << *xi;
        }
        std::cout << "\n";
        std::cout << std::flush;
    }
    if (zcols) {
        typedef std::vector<bool>::const_iterator vbci_t;
        uint j=0;
        for(vbci_t xj = colmap.begin() ; xj != colmap.end() ; xj++,j++) {
            if (!*xj) zcols->insert(j);
        }
        std::cout << "// " << zcols->size() << " zero cols:";
        for(suci_t xi = zcols->begin() ; xi != zcols->end() ; xi++) {
            std::cout << " " << *xi;
        }
        std::cout << "\n";
        std::cout << std::flush;
        if (zcols->size() != 0) {
            std::cerr << "// ZERO COLS FOUND -- MOST CERTAINLY A BUG !!\n";
        }
    }
}

void matrix_stats::operator()(std::string const & name) {
    /* It's probably forbidden to close the input stream before
     * reusing the offsets, so most probably it's forbidden to use
     * this entry point when slices are needed.
     */
    ifstream mtx;
    must_open(mtx, name);
    (*this)(mtx);
    mtx.close();
}

/* vim: set sw=4: */
