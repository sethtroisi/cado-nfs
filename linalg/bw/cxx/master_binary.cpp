#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <list>
#include <string.h>
#include <stdio.h>
#include <gmp.h>
#include <errno.h>
#include <assert.h>
#include <unistd.h>
#include <assert.h>
#include <boost/shared_array.hpp>
#include <utility>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <functional>

#define Lmacro(N, m, n) (iceildiv((N)+2*(n),(m))+iceildiv((N)+2*(n),(n))+10)

/* avoid this one ! */
// #define  VARIABLES_H_
// #define  STRUCTURE_H_

#include "constants.hpp"
// #include "masters.h"
#include "params.h"
#include "types.h"
#include "macros.h"
#include "auxfuncs.h"
#include "timer.h"
// #include "e_polynomial.h"
// #include "twisting_polynomials.h"
// #include "fft_on_matrices.hpp"
// #include "ops_poly.hpp"

// #include "master_common.h"
#include "auxfuncs.h"

#include "fmt.hpp"


/* output: F0x is the candidate number x. All coordinates of the
 * candidate are grouped in a row, and coefficients come in order (least
 * significant first
 */


/* Needed operations on critical path.
 *
 * compute_ctaf: compute the coefficient of degree t of the product of E
 * (m*b matrix) by PI (a b*b matrix), to be stored into small-e (m*b
 * matrix)
 *
 * zeroing out a column of small-e (m*b)
 *
 * do a column gaussian elimination on small-e:
 *   - column after column,
 *     - find the first non-zero bit in the column.
 *     - add this column to other column.
 *
 * Then polynomial operations.
 *
 * => small-e is better stored column-wise.
 * => big-E as well, and PI also.
 *
 * Extra un-critical stuff includes
 * - initialization -> computation of F0
 *                  -> computation of initial E = A * F0
 * - termination -> computation of final F = F0 * PI
 *
 * all this stuff has linear complexity.
 */

struct bcol;
struct bmat;
struct polmat;

struct bcol {/*{{{*/
    friend class bmat;
    uint nrows;
    inline uint stride() const { return BITS_TO_WORDS(nrows, ULONG_BITS); }
    private:
    boost::shared_array<ulong> x;
    public:
    bcol(uint nrows)
        : nrows(nrows), x(new ulong[stride()])
    {}
    bcol() : nrows(0) {}
    void clear() { nrows = 0; x.reset(); }
    bcol(bcol const& a) : nrows(a.nrows), x(a.x) {}
    bcol& operator=(bcol const& a) {
        nrows = a.nrows;
        x = a.x;
        return *this;
    }
    void zero() {
        memset(x.get(), 0, stride() * sizeof(ulong));
    }
    bool is_zero_likely() const {
        for(uint i = 0 ; i < stride() ; i++) {
            if (x[i] != 0) return false;
        }
        return true;
    }
    inline ulong coeff(uint i) const
    {
        uint offset = i / ULONG_BITS;
        ulong shift = i % ULONG_BITS;
        return (x[offset] >> shift) & 1UL;
    }
    uint ffs() const {
        uint z = 0;
        const ulong * src = x.get();
        ulong w;
        uint k;
        for(k = stride() ; k && !(w=*src++) ; k--) z += ULONG_BITS;
        if (k == 0) return UINT_MAX;
        return z + __builtin_ctzl(w);
    }
    void add(bcol const& a, unsigned long mask = 1UL) {
        for(uint l = 0 ; l < stride() ; l++) {
            x[l] ^= a.x[l] & -mask;
        }
    }
};/*}}}*/
struct bmat { /* This is for example small-e, an m times b matrix *//*{{{*/
    friend class polmat;
    uint nrows;
    uint ncols;
    private:
    boost::shared_array<ulong> x;
    boost::shared_array<uint> order;
    inline uint stride() const { return BITS_TO_WORDS(nrows, ULONG_BITS); }
    public:
    bmat(uint nrows, uint ncols)
        : nrows(nrows), ncols(ncols), x(new ulong[ncols*stride()]),
        order(new uint[ncols])
    {
        memset(x.get(), 0, ncols*stride()*sizeof(ulong));
        for(uint j = 0 ; j < ncols ; j++) order[j]=j;
    }
    bmat() : nrows(0), ncols(0) {}
    void clear() { nrows = ncols = 0; x.reset(); order.reset(); }
    bmat(bmat const& a) : nrows(a.nrows), ncols(a.ncols), x(a.x), order(a.order) {}
    bmat& operator=(bmat const& a) {
        nrows = a.nrows;
        ncols = a.ncols;
        x = a.x;
        order = a.order;
        return *this;
    }
    bmat clone() const {
        bmat dst(nrows, ncols);
        memcpy(dst.x.get(), x.get(), ncols*stride()*sizeof(ulong));
        memcpy(dst.order.get(), order.get(), ncols*sizeof(uint));
        return dst;
    }
    /* zero column j */
    void zcol(uint j) {
        memset(x.get() + order[j] * stride(), 0, stride() * sizeof(ulong));
    }
    bool col_is_zero(uint j) const {
        const ulong * src = x.get() + order[j] * stride();
        for(uint i = 0 ; i < stride() ; i++) {
            if (src[i] != 0) return false;
        }
        return true;
    }
    bool is_zero() const {
        for(uint i = 0 ; i < ncols * stride() ; i++) {
            if (x[i] != 0) return false;
        }
        return true;
    }
    inline ulong coeff(uint i, uint j) const
    {
        uint offset = i / ULONG_BITS;
        ulong shift = i % ULONG_BITS;
        return (x[order[j] * stride() + offset] >> shift) & 1UL;
    }
    /* add column j0 to column j. Optionally, b indicates a number of rows
     * which are known to be zero in both columns.
     * The mask corresponds to a multiplicating coefficient.
     */
    void acol(uint j, uint j0, uint b = 0, ulong mask = 1UL) {
        ulong * dst = x.get() + order[j] * stride();
        const ulong * src = x.get() + order[j0] * stride();
        for(uint l = b / ULONG_BITS ; l < stride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
    }
    void acol(uint j, bcol& a, uint b = 0, ulong mask = 1UL) {
        ulong * dst = x.get() + order[j] * stride();
        const ulong * src = a.x.get();
        for(uint l = b / ULONG_BITS ; l < stride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
    }
    uint ffs(uint j) const {
        uint z = 0;
        const ulong * src = x.get() + order[j] * stride();
        ulong w;
        uint k;
        for(k = stride() ; k && !(w=*src++) ; k--) z += ULONG_BITS;
        if (k == 0) return UINT_MAX;
        return z + __builtin_ctzl(w);
    }
    bcol extract_col(uint j) const
    {
        bcol a(nrows);
        memcpy(a.x.get(), x.get() + order[j] * stride(), stride() * sizeof(ulong));
        return a;
    }
    /* permute columns: column presently referred to with index i will
     * now be accessed at index p[i]
     */
    void perm(uint * p) {
        boost::shared_array<uint> norder(new uint[ncols]);
        for(uint i = 0 ; i < ncols ; i++) {
            norder[p[i]]=order[i];
        }
        order.swap(norder);
    }
    void addcoeff(uint i, uint j, ulong z)
    {
        /* cow check ??? */
        uint offset = i / ULONG_BITS;
        x[order[j] * stride() + offset] ^= z << (i % ULONG_BITS);
    }
};/*}}}*/
struct polmat { /* {{{ */
    /* polmat may represent:
     * - big-E, an m times b matrix of polynomials having nc max coeffs
     * - PI,     a b times b matrix of polynomials having nc max coeffs
     *
     * coeffs of big-E will eventually be stored in reversed bit order for
     * efficiency, but this is of little concern here.
     *
     * polmat may also represent other things, but for which speed is not
     * important.
     */
    uint nrows;
    uint ncols;
    uint ncoef;
    static bool critical;
    private:
    boost::shared_array<ulong> x;
    boost::shared_array<uint> order;
    boost::shared_array<int> _deg;
    inline uint stride() const { return BITS_TO_WORDS(ncoef, ULONG_BITS); }
    inline uint colstride() const { return nrows * stride(); }
    static void brev_warning();
    public:
    int& deg(uint j) { return _deg[order[j]]; }
    int const & deg(uint j) const { return _deg[order[j]]; }
    polmat(uint nrows, uint ncols, uint ncoef)
        : nrows(nrows), ncols(ncols), ncoef(ncoef),
        x(new ulong[ncols*colstride()]),
        order(new uint[ncols]),
        _deg(new int[ncols])
    {
        memset(x.get(), 0, ncols*colstride()*sizeof(ulong));
        for(uint j = 0 ; j < ncols ; j++) order[j]=j;
        for(uint j = 0 ; j < ncols ; j++) _deg[j]=-1;
    }
    void clear() {
        nrows = ncols = ncoef = 0;
        x.reset();
        order.reset();
        _deg.reset();
    }
    polmat() : nrows(0), ncols(0), ncoef(0) {}
    polmat(polmat const& a) :
        nrows(a.nrows),
        ncols(a.ncols),
        ncoef(a.ncoef),
        x(a.x),
        order(a.order),
        _deg(a._deg)
    {
        /* There's a catch here. Perhaps we should do copy on write on
         * the non-const methods ??? This could be done easily by calling
         * the private method cow(). But such a thing is quite bad
         * because it only works in restricted cases, and certainly not
         * in threaded situations.
         *
         * The article there: http://accu.org/index.php/journals/1413
         * raises a good point that one is much better off taking
         * advantage of copy elision. Heck, it's hard however to make
         * this assumption ?
         */
    }
    private:
    void cow() {
        if (x.get() != NULL && !x.unique()) {
            *this = clone();
        }
    }
    public:
    polmat& operator=(polmat const& a) {
        nrows = a.nrows;
        ncols = a.ncols;
        ncoef = a.ncoef;
        x = a.x;
        order = a.order;
        _deg = a._deg;
        return *this;
    }
    polmat clone() {
        polmat dst(nrows, ncols, ncoef);
        memcpy(dst.x.get(), x.get(), ncols*colstride()*sizeof(ulong));
        memcpy(dst.order.get(), order.get(), ncols*sizeof(uint));
        memcpy(dst._deg.get(), _deg.get(), ncols*sizeof(uint));
        return dst;
    }
    /* this handles expansion and shrinking */
    void resize(uint ncoef2) {
        /* cow check ??? */
        uint newstride = BITS_TO_WORDS(ncoef2, ULONG_BITS);
        if (newstride == stride()) {
            ncoef = ncoef2;
            return;
        }
        uint minstride = std::min(stride(),newstride);
        /* take the opportunity of reallocation for reordering columns. */
        polmat n(ncols,nrows,ncoef2);
        ulong * dst = n.x.get();
        for(uint j = 0 ; j < ncols ; j++) {
            const ulong * src = x.get() + order[j] * colstride();
            for(uint i = 0 ; i < nrows ; i++) {
                memcpy(dst, src, minstride * sizeof(ulong));
                dst += newstride;
                src += stride();
            }
        }
        x.swap(n.x);
        order.swap(n.order);
        _deg.swap(n._deg);
    }
    /* zero column j */
    void zcol(uint j) {
        /* cow check ??? */
        memset(x.get() + order[j] * colstride(), 0, colstride() * sizeof(ulong));
        _deg[order[j]] = -1;
    }
    void xmul_col(uint j, uint s=1) {
        /* cow check ??? */
        mp_limb_t * dst = x.get() + order[j] * colstride();
        mpn_lshift(dst, dst, colstride(), s);
        /* We may have garbage for the low bits. It may matter, or maybe
         * not. If it does, call xclean0_col */
        _deg[order[j]] += _deg[order[j]] >= 0;
    }
    void xmul_poly(uint i, uint j, uint s=1) {
        /* cow check ??? */
        ulong * dst = x.get() + (order[j] * nrows + i) * stride();
        mpn_lshift(dst, dst, stride(), s);
        /* We may have garbage for the low bits. It may matter, or maybe
         * not. If it does, call xclean0_col */
        // _deg[order[j]] += _deg[order[j]] >= 0;
    }
    void addpoly(uint i, uint j, const ulong * src) {
        ulong * dst = x.get() + (order[j] * nrows + i) * stride();
        for(uint k = 0 ; k < stride() ; k++)
            dst[k] ^= src[k];
    }
    void xclean0_col(uint j) {
        /* cow check ??? */
        ulong * dst = x.get() + order[j] * colstride();
        for(uint i = 0 ; i < nrows ; i++, dst += stride())
            dst[0] &= ~1UL;
        _deg[order[j]] -= _deg[order[j]] == 0;
    }
    /* add column j times mask to column i. */
    void acol(uint i, uint j, ulong mask = 1UL) {
        /* cow check ??? */
        ulong * dst = x.get() + order[i] * colstride();
        const ulong * src = x.get() + order[j] * colstride();
        for(uint l = 0 ; l < colstride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
        // ASSERT(_deg[order[i]] >= _deg[order[j]]);
        _deg[order[i]] = std::max(_deg[order[j]], _deg[order[i]]);
    }
    /* permute columns: column presently referred to with index i will
     * now be accessed at index p[i]
     */
    void perm(uint * p) {
        /* cow check ??? */
        boost::shared_array<uint> norder(new uint[ncols]);
        for(uint i = 0 ; i < ncols ; i++) {
            norder[p[i]]=order[i];
        }
        order.swap(norder);
    }
    uint ffs(uint j, uint k) const {
        uint offset = k / ULONG_BITS;
        ulong mask = 1UL << (k % ULONG_BITS);
        const ulong * src = x.get() + j * colstride() + offset;
        for(uint z = 0 ; z < nrows ; z++) {
            if (*src & mask) return z;
            src += stride();
        }
        return UINT_MAX;
    }
    /* cow check ??? */
    /* XXX BEWARE ! This kind of stuff really opens a can of worms !
     * It's possible to do polmat a(b); ulong * ptr = a.poly(i,j), and,
     * by changing ptr, changing b as well. Perhaps a const analogue to
     * this method would be desirable
     */
    ulong * poly(uint i, uint j)
    {
        return x.get() + (order[j] * nrows + i) * stride();
    }
    ulong const * poly(uint i, uint j) const
    {
        return x.get() + (order[j] * nrows + i) * stride();
    }
    ulong * col(uint j) { return poly(0, j); }
    ulong const * col(uint j) const { return poly(0, j); }
    /* shift is understood ``shift left'' (multiply by X) */
    void import_col_shift(uint k, polmat const& a, uint j, int s)
    {
        BUG_ON(a.nrows != nrows);
        ulong const * src = a.col(j);
        ulong * dst = col(k);
        if (s > 0) {
            mpn_lshift(dst, src, colstride(), s);
        } else {
            mpn_rshift(dst, src, colstride(), -s);
        }
    }

    /* these accessors are not the preferred ones for performance */
    inline ulong coeff(uint i, uint j, uint k) const
    {
        BUG_ON(critical);
        brev_warning();
        uint offset = k / ULONG_BITS;
        ulong shift = k % ULONG_BITS;
        return poly(i,j)[offset] >> shift & 1UL;
    }
    bmat extract_coeff(uint k) const
    {
        BUG_ON(critical);
        brev_warning();
        bmat a(nrows, ncols);
        for(uint j = 0 ; j < ncols ; j++) {
            for(uint i = 0 ; i < nrows ; i++) {
                uint boffset = i / ULONG_BITS;
                ulong bmask = 1UL << (i % ULONG_BITS);
                ulong coeff = this->coeff(i,j,k);
                a.x[j*a.stride() + boffset] ^= bmask & -coeff;
            }
        }
        return a;
    }
    void addcoeff(uint i, uint j, uint k, ulong z)
    {
        /* cow check ??? */
        BUG_ON(critical);
        brev_warning();
        uint offset = k / ULONG_BITS;
        poly(i,j)[offset] ^= z << (k % ULONG_BITS);
    }
};

bool polmat::critical = false;
void polmat::brev_warning()
{
    static bool told;
    if (!told) {
        WARNING("attention : bit reverse might come here !");
    }
    told = 1;
}

void dbmat(bmat const * x)
{
    for(uint i = 0 ; i < x->nrows ; i++) {
        for(uint j = 0 ; j < x->ncols ; j++) {
            std::cout << x->coeff(i,j);
        }
        std::cout << "\n";
    }
}
void dpmat(polmat const * pa, uint k)
{
    bmat x = pa->extract_coeff(k);
    dbmat(&x);
}
/*}}}*/

/* Computes a times b, but only coefficients of degree k0 to k1
 * (inclusive). Hence the resulting number of coefficients is at most
 * k1-k0+1
 */
void multiply_slow(polmat& res, polmat const &a, polmat const &b,
        uint k0 = 0, uint k1 = UINT_MAX)
{
    BUG_ON(polmat::critical);
    BUG_ON(a.ncols != b.nrows);
    uint ncmax = std::min(std::max(0u,a.ncoef + b.ncoef - 1 - k0),
            a.ncoef + b.ncoef - 1);
    if (k1 != UINT_MAX) {
        ncmax = std::min(ncmax, k1 - k0 + 1);
    }
    res = polmat(a.nrows, b.ncols, ncmax);
    for (uint j = 0; j < b.ncols; j++) {
        for (uint i = 0; i < a.nrows; i++) {
            uint ncres = 0;
            for (uint k = 0; k < a.ncols; k++) {
                for (uint s = 0; (int) s <= a.deg(k); s++) {
                    for (uint t = 0; (int) t <= b.deg(j); t++) {
                        if (s + t < k0) continue;
                        if (s + t > k1) continue;
                        uint u = s + t - k0;
                        if (u + 1 > ncres)
                            ncres = u + 1;
                        res.addcoeff(i, j, u,
                                a.coeff(i, k, s) & b.coeff(k, j, t));
                    }
                }
            }
            res.deg(j) = ncres - 1;
        }
    }
}

namespace globals {
    uint nrows;
    uint m,n;
    uint t,t0,total_work;
    std::vector<uint> delta;
    std::vector<uint> chance_list;

    polmat E;
    bmat e0;

    // F0 is exactly the n x n identity matrix, plus the X^(s-exponent)e_{cnum}
    // vectors. Here we store the cnum,exponent pairs.
    std::vector<std::pair<uint, uint> > f0_data;

}

#if 0/*{{{*/
// #include "bw_scalar.h"
// #include "variables.h"
// #include "structure.h"
// #include "gmp-hacks.h"
// #include "modulus_hacks.h"

/* This code assumes that m and n are both multiples of the machine word
 * size. */

/* accessors for matrices, matrix polynomials and so on. */
/* mxr stands for matrix accessor. */

struct bxr {    /* (immediate) bit accessor */
    ulong * x;
    uint n;
    bxr(unsigned long * x, uint n) : x(x), n(n) {};
    ulong get(uint j) const {
        return (x[j / ULONG_BITS] >> (j % ULONG_BITS)) & 1UL;
    }
    void set(uint j) { x[j / ULONG_BITS] |= 1UL << (j % ULONG_BITS); }
    void clr(uint j) { x[j / ULONG_BITS] &= ~(1UL << (j % ULONG_BITS)); }
    void add(uint j, ulong a) { x[j / ULONG_BITS] ^= a << (j % ULONG_BITS); }
    void add(bxr const& b) {
        ASSERT(n == b.n);
        for(uint k = 0 ; k < BITS_TO_WORDS(n, ULONG_BITS) ; k ++) {
            x[k] ^= b.x[k];
        }
    }
};

struct strided_bxr {    /* bit accessor */
    ulong * x;
    uint n;
    uint stride;
    uint offset;
    strided_bxr(ulong * x, uint n, uint s = 1, uint o = 0) :
        x(x), n(n), stride(s), offset(o)
    {
        ASSERT(stride % ULONG_BITS == 0);
    }
    ulong get(uint j) const { return bxr(x,n).get(offset + stride * j); }
    void set(uint j) { bxr(x,n).set(x, offset + stride * j); }
    void clr(uint j) { bxr(x,n).clr(x, offset + stride * j); }
    void add(uint j, ulong a) { bxr(x,n).add(x, offset + stride * j, a); }
    void add(strided_bxr const& b) {
        ASSERT(n == b.n);
        ASSERT(stride == b.stride);
        ASSERT(offset == b.offset);
        ASSERT(stride >= ULONG_BITS);   // otherwise it's pointless
        uint j = offset;
        ulong mask = 1UL << (j % ULONG_BITS);
        /* We do insist on having a word-aligned stride here, of course */
        for(uint k = 0 ; k < n ; k++, j+=stride) {
            x[j / ULONG_BITS] ^= mask & b.x[j / ULONG_BITS];
        }
    }
};

struct mxr_byrow {
    ulong * x;
    uint m;
    uint n;
    mxr_byrow(ulong * x, uint m, uint n) : x(x), m(m), n(n) {}
    typedef bxr row_xr;
    typedef strided_bxr col_xr;
    inline uint stride() { return BITS_TO_WORDS(n, ULONG_BITS); }
    row_xr row(uint i) const { return bxr(x + i * stride, n); }
    col_xr col(uint j) const { return bxr_strided(x, m, j, stride); }
    ulong get(uint i, uint j) const {
        return row(i).get(j);
    }
    void set(uint i, uint j) { row(i).set(j); }
    void clr(uint i, uint j) { row(i).clr(j); }
    void add(uint i, uint j, ulong a) { row(i).add(j, a); }
    void add(mxr_byrow const& b) {
        ASSERT(n == b.n);
        ASSERT(m == b.m);
        for(uint k = 0 ; k < m * stride ; k++) {
            x[k] ^= b.x[k];
        }
    }
};

/* could more or less be made as a wrap-up on top of the other one */
struct mxr_bycol {
    ulong * x;
    uint m;
    uint n;
    typedef strided_bxr row_xr;
    typedef bxr col_xr;
    mxr_bycol(ulong * x, uint m, uint n) : x(x), m(m), n(n) {}
    inline uint stride() { return BITS_TO_WORDS(m, ULONG_BITS); }
    row_xr row(uint j) const { return bxr_strided(x, n, i, stride); }
    col_xr col(uint i) const { return bxr(x + j * stride, m); }
    ulong get(uint i, uint j) const {
        return col(j).get(i);
    }
    void set(uint i, uint j) { col(j).set(i); }
    void clr(uint i, uint j) { col(j).clr(i); }
    void add(uint i, uint j, ulong a) { col(j).add(i, a); }
    void add(mxr_byrow const& b) {
        ASSERT(n == b.n);
        ASSERT(m == b.m);
        for(uint k = 0 ; k < n * stride ; k++) {
            x[k] ^= b.x[k];
        }
    }
};

struct prxr_bydeg {
    ulong * x;
    uint nc;    /* max number of coeffs */
    uint n;
    prxr_byrow(ulong * x, uint nc, uint n) : x(x), nc(nc), n(n) {}
    typedef bxr pol_xr;
    typedef strided_bxr row_xr;
    inline uint stride() { return BITS_TO_WORDS(nc, ULONG_BITS); }
    row_xr row(uint i) const { return bxr_strided(x, n, i, tride); }
    pol_xr pol(uint k) const { return bxr(x + k * stride, nc); }
    ulong get(uint i, uint k) const {
        return pol(i).get(j);
    }
    void set(uint i, uint j) { row(i).set(j); }
    void clr(uint i, uint j) { row(i).clr(j); }
    void add(uint i, uint j, ulong a) { row(i).add(j, a); }
    void add(mxr_byrow const& b) {
        ASSERT(n == b.n);
        ASSERT(m == b.m);
        for(uint k = 0 ; k < m * stride ; k++) {
            x[k] ^= b.x[k];
        }
    }
};

#endif/*}}}*/

// To multiply on the right an m x n matrix A by F0, we start by copying
// A into the first n columns. Since we're also dividing out by X^t0, the
// result has to be shifted t0 positions to the right.
// Afterwards, column n+j of the result is column cnum[j] of A, shifted
// exponent[j] positions to the right.
void compute_E_from_A(polmat const &a)
{
    using namespace globals;
    E = polmat(n, m + n, a.ncoef - t0);
    for(uint j = 0 ; j < n ; j++) {
        E.import_col_shift(j, a, j, - (int) t0);
    }
    for(uint j = 0 ; j < m ; j++) {
        uint cnum = f0_data[j].first;
        uint exponent = f0_data[j].second;
        E.import_col_shift(n + j, a, cnum, - (int) exponent);
    }
}

// F is in fact F0 * PI.
// To multiply on the *left* by F, we cannot work directly at the column
// level (we could work at the row level, but it does not fit well with
// the way data is organized). So it's merely a matter of adding
// polynomials.
void compute_final_F_from_PI(polmat& F, polmat const& pi)
{
    using namespace globals;
    // We take t0 rows, so that we can do as few shifts as possible
    polmat tmpmat(t0,1,pi.ncoef);
    F = polmat(n, m + n, globals::t0 + pi.ncoef - 1);
    using namespace std;
    for(uint i = 0 ; i < n ; i++) {
        // What contributes to entries in this row ?
        vector<pair<uint, uint> > l;
        set<uint> sexps;
        // l.push_back(make_pair(i,0));
        for(uint j = 0 ; j < m ; j++) {
            if (f0_data[j].first == i) {
                l.push_back(make_pair(n + j, f0_data[j].second));
                sexps.insert(f0_data[j].second);
            }
        }
        vector<uint> exps(sexps.begin(), sexps.end());
        // So the i-th row of f0 has a 1 at position i, and then
        // X^(t0-f.second) at position n+j whenever f.first == i.
        //
        // Now fill in the row.
        for(uint j = 0 ; j < m + n ; j++) {
            // We zero out the whole column, it's less trouble.
            tmpmat.zcol(0);
            for(uint k = 0 ; k < l.size() ; k++) {
                tmpmat.addpoly(l[k].second, 0, pi.poly(l[k].first, j));
            }
            F.addpoly(i,j,pi.poly(i,j));
            for(uint k = 0 ; k < exps.size() ; k++) {
                tmpmat.xmul_poly(exps[k],0,t0-exps[k]);
                F.addpoly(i,j,tmpmat.poly(exps[k],0));
            }
        }
    }
    for(uint j = 0 ; j < m + n ; j++) {
        F.deg(j) = t0 + pi.deg(j);
    }
}


const char *a_meta_filename = "A-%02d-%02d";
const char *f_meta_filename = "F%02d";
const char *pi_meta_filename = "pi-%d-%d";
const char *f_base_filename = "F_INIT_QUICK";


std::pair<polmat,int> read_data_for_series()/*{{{*/
{
    using namespace globals;

    FILE * files[m][n];
    uint nbys[n];

    polmat a(m,n,total_work+1);

    uint i, j;
    uint k;

    for (j = 0; j < n; j++) { nbys[j] = 0; }

    for (i = 0; i < m; i++) {
        for (j = 0; j < n;) {
            uint y = 0;
            char * ptr;
            int d;

            char filename[FILENAME_LENGTH];
            sprintf(filename, a_meta_filename, i, j);
            files[i][j] = fopen(filename, "r");
            if (files[i][j] == NULL) {
                die("fopen(%s) : %s", errno, filename,
                        strerror(errno));
            }
            printf("Reading file %s", filename);

            /* NOTE : we drop the first coefficient, because
             * we mean to work with the sequence generated
             * on x and By, so that we obtain a generator
             * afterwards.
             */

            char row[1024];
            ulong blah;
            /* We also use this in order to read the number
             * of coefficients per row */
            fgets(row, sizeof(row), files[i][j]);
            ptr = row;
            for(y = 0 ; sscanf(ptr,"%lu%n",&blah,&d) >= 1 ; ptr += d, y++);
            if (y == 0) {
                fprintf(stderr, "\nproblem while reading %s, line = %s\n",
                        filename, row);
                abort();
            }
            if (y > 1) {
                printf(" [ %d values per row ]", y);
            }
            printf("\n");
            if (nbys[j] != 0) {
                BUG_ON(y != nbys[j]);
            }
            nbys[j] = y;
            j += y;
        }
    }

    uint read_coeffs;

    for (i = 0; i < m; i++) {
        ulong * pol[n];
        for (j = 0; j < n; j ++) {
            pol[j] = a.poly(i,j);
        }
        ulong mask = 1UL;
        uint shift = 0;
        for (k = 0; k <= total_work; k++) {
            uint rc = 0;
            for (j = 0; j < n; j += nbys[j]) {
                uint l;
                for(l = 0 ; l < nbys[j] ; l++) {
                    ulong blah;
                    rc += 1 == fscanf(files[i][j], "%lu", &blah);
                    pol[j + l][shift] |= mask & -blah;
                }
            }
            if (rc == 0) {
                break;
            } else if (rc != n) {
                fprintf(stderr, "A files not in sync ; "
                        "read only %u coeffs in X^%d, row %d\n", rc, k, i);
                break;
            }
            mask <<= 1;
            shift += mask == 0;
            mask += mask == 0;
        }
        if (i == 0) {
            read_coeffs = k;
        } else {
            /* This happens when A files have gone live */
            BUG_ON(k != read_coeffs);
        }
    }
    for (j = 0; j < n; j ++) {
        a.deg(j) = read_coeffs - 1;
    }

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j+= nbys[j]) {
            fclose(files[i][j]);
        }
    }

    printf("Stopped after reading %u coefficients (not counting 1st)\n",
            read_coeffs);

    return std::make_pair(a,read_coeffs);
}/*}}}*/
void write_polmat(polmat const& P, const char * fn)/*{{{*/
{
    std::ofstream f(fn);
    if (!f.is_open()) {
        perror(fn);
        exit(1);
    }

    for(uint k = 0 ; k < P.ncoef ; k++) {
        for(uint i = 0 ; i < P.nrows ; i++) {
            for(uint j = 0 ; j < P.ncols ; j++) {
                if (j) f << " ";
                if ((int) k <= P.deg(j)) {
                    f << P.coeff(i,j,k);
                } else {
                    f << 0;
                }
            }
            f << "\n";
        }
        f << "\n";
    }
    printf("Written f to %s\n",fn);
}/*}}}*/

bool recover_f0_data(const char * fn)
{
    std::ifstream f(fn);

    using namespace globals;
    for(uint i = 0 ; i < m ; i++) {
        uint exponent,cnum;
        if (!(f >> cnum >> exponent))
            return false;
        f0_data.push_back(std::make_pair(cnum,exponent));
    }
    std::cout << fmt("recovered init data % on disk with t0=%\n") % fn % t0;
    return true;
}

bool write_f0_data(const char * fn)
{
    std::ofstream f(fn);

    using namespace globals;
    for(uint i = 0 ; i < m ; i++) {
        uint cnum = f0_data[i].first;
        uint exponent = f0_data[i].second;
        if (!(f << " " << cnum << " " << exponent))
            return false;
    }
    f << std::endl;
    std::cout << fmt("written init data % to disk\n") % fn;
    return true;
}

// the computation of F0 says that the column number cnum of the
// coefficient X^exponent of A increases the rank.
//
// These columns are gathered into an invertible matrix at step t0, t0
// being strictly greater than the greatest exponent above. This is so
// that there is no trivial column dependency in F0.
void set_t0_delta_from_F0()
{
    using namespace globals;
    t0 = f0_data.back().second + 1;     // see above for the +1
    for(uint j = 0 ; j < m + n ; j++) {
        delta[j] = t0;
        // Even irrespective of the value of the exponent, we get t0 for
        // delta.
    }
}

void write_F0_from_F0_quick()
{
    using namespace globals;
    uint s = t0;
    /* Now build f */
    polmat F0(n, m + n, s + 1);

    /* First n columns: identity matrix */
    /* rest: X^(s-exponent)cnum's */
    for(uint i=0;i<n;i++) {
        F0.addcoeff(i,i,0,1);
        F0.deg(i) = 0;
    }
    for(uint j=0;j<m;j++) {
        uint cnum = f0_data[j].first;
        uint exponent = f0_data[j].second;
        F0.addcoeff(cnum,n+j,t0-exponent,1);
        F0.deg(n + j) = t0-exponent;
    }
    write_polmat(F0, "F_INIT");
}

#if 0
bool read_f_init(polmat & F0, const char * fn, uint nr, uint nc)/*{{{*/
{
    ulong blah;

    FILE * f=fopen(fn,"r");

    if (f == NULL)
        return false;

    using namespace globals;

    F0 = polmat(nr,nc,1);

    printf("Found F_INIT data on disk. Restoring\n");

    uint k;

    for(k=0;;k++) {
        uint rc = 0;
        for(uint i=0;i<nr;i++) for(uint j=0;j<nc;j++) {
            if (fscanf(f, "%lu", &blah) != 1) 
                continue;
            rc++;
            /* ensures proper allocation. */
            F0.resize(k+1);
            F0.addcoeff(i,j,k,blah);
        }
        if (rc == 0) {
            printf("Stopped after reading %u coeffs\n", k);
            break;
        } else if (rc != nr * nc) {
            die("file not in sync\n", 1);
        }
    }
    fclose(f);

    for(uint j = 0 ; j < nc ; j++) {
        delta[j] = 0;
        for(int lk = k - 1 ; lk >= 0 ; lk--) {
            if (F0.ffs(j, lk) != UINT_MAX) {
                delta[j] = lk + 1;
                F0.deg(j) = lk;
                break;
            }
        }
    }
    return k > 0;
}/*}}}*/
// give_poly_rank_info is no longer useful, compute_f_init does much
// more, much better.
void give_poly_rank_info(polmat& a, int deg)    /* informative *//*{{{*/
{
    uint k;
    uint grank=globals::m;
    uint m = std::min(deg, 200);
    std::cout << fmt("Rank of first % A matrices:") % m << std::flush;
    for(k = 0 ; k < m ; k++) {
        uint r = compute_rank(a.extract_coeff(k));
        std::cout << " " << r;
        if (k < m / 2 || r >= grank)
            grank = r;
        std::cout << std::flush;
    }
    std::cout << std::endl;
    if (grank != globals::m) {
        std::cout << "The program will almost surely fail"
                " because the choice of X-vectors was bad\n";
        sleep(3);
    }
}/*}}}*/
static uint compute_rank(bmat const& a0)/*{{{*/
{
    bmat a = a0.clone();
    uint i,j,k;
    uint rank;

    /* Pay attention here, this is a gaussian elimination on
     * *columns* */

    rank = 0 ;
    for(j = 0 ; j < a.ncols ; j++) {
        /* Find the pivot inside the column. */
        i = a.ffs(j);
        if (i == UINT_MAX)
            continue;
        assert(rank < a0.nrows && rank < a0.ncols);
        rank++;
        /* Cancel this coeff in all other columns. */
        for(k = j + 1 ; k < a0.ncols ; k++) {
            /* TODO: Most probably this branch avoiding trick could be
             * deleted, I doubt it gains anything. */
            a.acol(k,j,i,a.coeff(i,k));
        }

    }
    return rank;
}/*}}}*/
#endif

void bw_commit_f(polmat& F) /*{{{*/
{
    char filename[FILENAME_LENGTH];

    using namespace globals;
    using namespace std;

    for (uint j = 0; j < m + n ; j++) {
        sprintf(filename, f_meta_filename, j);
        ofstream f(filename);
        if (!f.is_open()) {
            perror("writing f");
            eternal_sleep();
        }
        printf("Writing %s\n", filename);
        for (int k = 0; k <= F.deg(j); k++) {
            for (uint i = 0; i < n; i++) {
                f << F.coeff(i,j,k) << (i == (n-1) ? "\n" : " ");
            }
        }
    }
}
/*}}}*/


void compute_f_init(polmat& A)/*{{{*/
{
    using namespace globals;

    /* build F_INIT */
    printf("Computing t0\n");

    /* For each integer i between 0 and m-1, we have a column, picked
     * from column cnum[i] of coeff exponent[i] of A which, once reduced modulo
     * the other ones, has coefficient at row pivots[i] unequal to zero.
     */
    bcol pcols[m];
    uint pivots[m], exponent[m], cnum[m];

    uint r = 0;
    for(uint k=0;r < m && k<A.ncoef;k++) {
        bmat amat = A.extract_coeff(k);
        for(uint j=0;r < m && j<n;j++) {
            /* copy column j, coeff k into acol */
            bcol acol = amat.extract_col(j);
            /* kill as many coeffs as we can */
            for(uint v=0;v<r;v++) {
                uint u=pivots[v];
                /* the v-th column in the matrix reduced_rank is known to
                 * kill coefficient u (more exactly, to have a -1 as u-th
                 * coefficient, and zeroes for the other coefficients
                 * referenced in the pivots[0] to pivots[v-1] indices).
                 */

                acol.add(pcols[v],acol.coeff(u));
            }
            uint u = acol.ffs();
            if (u == UINT_MAX) {
                printf("[X^%d] A, col %d does not increase rank (still %d)\n",
                        k,j,r);
                continue;
            }

            /* Bingo, it's a new independent col. */
            pivots[r]=u;
            cnum[r]=j;
            exponent[r]=k;

            f0_data.push_back(std::make_pair(cnum[r], exponent[r]));

            /* TODO: For non-binary stuff, multiply the column so that
             * acol[u] becomes -1
             */
            pcols[r] = acol;
            r++;

            printf("[X^%d] A, col %d increases rank to %d (head row %d)\n",
                    k,j,r,u);
        }
    }

    t0 = exponent[r-1] + 1;
    printf("Found satisfying init data for t0=%d\n", t0);
                    
    printf("Init e0 matrix\n");
    for(uint i = 0 ; i < m ; i++) {
        for(uint j = 0 ; j < n ; j++) {
            std::cout << A.coeff(i,j,t0);
        }
        for(uint j = 0 ; j < m ; j++) {
            std::cout << A.coeff(i,cnum[j],exponent[j]);
        }
        std::cout << "\n";
    }


    if (r!=m) {
        printf("This amount of data is insufficient. "
                "Cannot find %u independent cols within A\n",m);
        exit(1);
    }

}/*}}}*/

void print_deltas()
{
    using namespace globals;
    std::cout << fmt("[t=%[w4]] delta =") % t;
    for(uint i = 0 ; i < m + n ; i++) {
        std::cout << " " << delta[i];
    }
    std::cout << "\n";
}

/*{{{*//* bw_init
 *
 * fill in the structures with the data available on disk. a(X) is
 * fetched this way. f(X) is chosen in order to make [X^0]e(X)
 * nonsingular. Once this condition is satisfied, e(X) is computed. All
 * further computations are done with e(X).
 *
 * a(X) is thrown away as soon as e(X) is computed.
 *
 * At a given point of the algorithm, one can bet that most of the data
 * for e(X) is useless.
 *
 * XXX Therefore it is wise to shrink e(x) aggressively. The new
 * organisation of the data prevents e(x) from being properly swapped
 * out. XXX TODO !!!
 *
 */
static void bw_init(void)
{
    uint read_coeffs;
    polmat A;

    using namespace globals;

    printf("Using A(X) div X in order to consider Y as starting point\n");
    total_work--;

    {
	printf("Reading scalar data in polynomial ``a''\n");
	std::pair < polmat, uint > a_read = read_data_for_series();
	A = a_read.first;
	read_coeffs = a_read.second;
    }

    /* informational *//*{{{ */
    if (read_coeffs < total_work) {
	printf("Since we do not have the full information yet "
	       "(only %d coefficients while %d were needed)"
	       ", we\n"
	       "merely can tell if the system looks like doable\n",
	       read_coeffs, total_work);
    }





    /*}}} */
    /* Data read stage completed. */
    /* TODO. Prepare the FFT engine for handling polynomial
     * multiplications up to ncmax coeffs */
    // uint ncmax = (total_work<<1)+2;
    // ft_order_t::init(ncmax, ops_poly_args);

    delta.assign(m + n, -1);

    if (!recover_f0_data(f_base_filename)) {
        // This is no longer useful
	// give_poly_rank_info(A, read_coeffs - 1);
	compute_f_init(A);
	write_f0_data(f_base_filename);
    }
    set_t0_delta_from_F0();

    // this is not used. It's only a debugging aid for the helper magma
    // code.
    write_F0_from_F0_quick();

    /*
        if (!read_f_init(F0, f_base_filename, n, m + n)) {
            give_poly_rank_info(A, read_coeffs - 1);
            compute_f_init(A);
            write_f0_data(F0, f_base_filename);
        }
    */
    t = t0;
    printf("t0 = %d\n", t0);


    chance_list.assign(m + n, 0);
    /*
       for(j=0;j< m + n ;j++) {
       global_delta[j]  = t_counter;
       chance_list[j]   = 0;
       }
     */

    total_work = read_coeffs;
    // A must be understood as A_computed + O(X^total_work)
    // therefore it does not make sense to compute coefficients of E at
    // degrees above total_work-1
    std::cout << "Computing value of E(X)=A(X)F(X) "
        << fmt("(degree %) [ +O(X^%) ]\n") % (total_work-1) % total_work;
    
    // multiply_slow(E, A, F0, t0, total_work-1);
    compute_E_from_A(A);

    printf("Throwing out a(X)\n");
    A.clear();
}
/*}}} */

#if 0/*{{{*/
/* Sort the degree table delta, with local ordering obeying the columns
 * of pi. The permutation is applied to delta, but not to pi.
 */
static void column_order(unsigned int * perm, int * delta, struct t_poly * pi)
{
    struct xcol_id * tmp_delta;
    int i;

    tmp_delta   = (xcol_id *) malloc((m+n)  * sizeof(struct xcol_id));

    for(i=0;i< m + n ;i++) {
        tmp_delta[i].pos=i;
        tmp_delta[i].deg=delta[i];
        tmp_delta[i].sdeg=pi->degnom[pi->clist[i]];
    }

    qsort(tmp_delta,(m+n),sizeof(struct xcol_id),(sortfunc_t)&col_cmp);

    for(i=0;i< m + n ;i++) {
        perm[i]=tmp_delta[i].pos;
        delta[i]=tmp_delta[i].deg;
    }

    free(tmp_delta);
}


static void e_transvec(bw_mbmat e,
        unsigned int * clist,
        int o_j1,
        int o_j2,
        int ikill,
        bw_scalar lambda)
{
    int i;
    int j1,j2;

    j1=clist[o_j1];
    j2=clist[o_j2];

    bw_scalar_set_zero(mbmat_scal(e,ikill,j2));

    for(i=ikill+1;i<m;i++) {
        addmul( mbmat_scal(e,i,j2),
                mbmat_scal(e,i,j1),
                lambda);
        k_reduce(mbmat_scal(e,i,j2));
    }
}
#endif/*}}}*/

static void rearrange_ordering(polmat & PI, uint piv[])
{
    /* Sort the columns. It might seem merely cosmetic and useless to
     * sort w.r.t both the global and local nominal degrees. In fact, it
     * is crucial for the corectness of the computations. (Imagine a
     * 2-step increase, starting with uneven global deltas, and hitting
     * an even situation in the middle. One has to sort out the local
     * deltas to prevent trashing the whole picture).
     *
     * The positional sort, however, *is* cosmetic (makes debugging
     * easier).
     */
    using namespace std;
    using namespace globals;
    typedef pair<pair<uint,uint>, int> corresp_t;
    vector<corresp_t> corresp(m+n);
    for(uint i = 0 ; i < m + n ; i++) {
        corresp[i] = make_pair(make_pair(delta[i],PI.deg(i)), i);
    }
    sort(corresp.begin(), corresp.end(), less<corresp_t>());
    uint perm[m+n];
    std::vector<uint> permuted_global_deltas(m+n,-1);
    for(uint i = 0 ; i < m + n ; i++) {
        perm[corresp[i].second] = i;
        permuted_global_deltas[i] = delta[corresp[i].second];
    }
    if (piv) {
        for(uint i = 0 ; i < m ; i++) {
            piv[i] = perm[piv[i]];
        }
    }
    PI.perm(perm);
    e0.perm(perm);
    swap(permuted_global_deltas, delta);
}

static void gauss(uint piv[], polmat& PI)
{
    /* Do one step of (column) gaussian elimination on e0. The columns
     * are assumed to be in the exact order that corresponds to the
     * ordering of the columns of PI (therefore relative to delta)
     */

    /* Note that we do *NOT* modify E. It seems to be the fastest option,
     * although the other deserves being investigated as well (see
     * old/master2.c, there are two versions of the quadratic algorithm).
     */
    uint i,j,k;
    uint rank;

    rank = 0 ;

    using namespace globals;

    /*
    std::cout << "Input matrix\n";
    dbmat(&e0);
    */

    for(j = 0 ; j < e0.ncols ; j++) {
        /* Find the pivot inside the column. */
        i = e0.ffs(j);
        if (i == UINT_MAX)
            continue;
        assert(rank < e0.nrows && rank < e0.ncols);
        // std::cout << fmt("col % is the %-th pivot\n") % j % rank;
        piv[rank++] = j;
        /* Cancel this coeff in all other columns. */
        for(k = j + 1 ; k < e0.ncols ; k++) {
            /* TODO : Over the binary field, this branch avoiding trick
             * could most probably be deleted, I doubt it gains anything. */
            ulong c = e0.coeff(i,k);
            /* add c times column j to column k */
            e0.acol(k,j,i,c);
            // E.acol(k,j,c);
            // This one is tempting, but it's a wrong assert.
            // ASSERT(PI.deg(j) <= PI.deg(k));
            // ASSERT(delta[j] <= delta[k]);
            ASSERT(std::make_pair(delta[j],PI.deg(j)) <= std::make_pair(delta[k],PI.deg(k)));
            PI.acol(k,j,c);
        }
        PI.xmul_col(j);
        // E.xmul_col(j);
        delta[j]++;
    }
    BUG_ON (rank != m);

    /*
    std::cout << "Invertible e0 matrix\n";
    for(uint i = 0 ; i < m ; i++) {
        for(uint j = 0 ; j < m ; j++) {
            std::cout << e0.coeff(i,piv[j]);
        }
        std::cout << "\n";
    }
    std::cout << "[ rows";
    for(uint j = 0 ; j < m ; j++) {
        std::cout << " " << piv[j];
    }
    std::cout << " ]\n";
    */

    /*
    std::cout << "Invertible e0 matrix\n";
    for(uint i = 0 ; i < m ; i++) {
        for(uint j = 0 ; j < m ; j++) {
            std::cout << e0.coeff(i,piv[j]);
        }
        std::cout << "\n";
    }
    std::cout << "[ rows";
    for(uint j = 0 ; j < m ; j++) {
        std::cout << " " << piv[j];
    }
    std::cout << " ]\n";
    */

}

void read_mat_file_header(const char *name)
{
    FILE *f = fopen(name, "r");
    if (f == NULL) {
        die("fopen(%s): %s", errno, name, strerror(errno));
    }
    char modulus[512];
    uint ncols;
    fscanf(f, "// %u ROWS %u COLUMNS ; MODULUS %s",
            &globals::nrows, &ncols, modulus);
    BUG_ON(strcmp(modulus,"2") != 0);
    BUG_ON(globals::nrows != ncols);

    fclose(f);
}

#if 0/*{{{*/
/*
 * ec : A raw polynomial matrix
 * pi : The permutation applied to previous data to obtain ec.
 * e  : [X^t] ec * pi
 * t  : the variable above
 * pivlist: NULL : transport the transvections automatically on ec
 *         !NULL : don't and recompute extensively e for next time if applicable
 */
static void bw_gauss_onestep(bw_mbmat e,
			     struct e_coeff *ec,
			     struct t_poly *pi,
			     int *delta, unsigned int *pivlist, double *tt)
{
    int i, j, k, jr;
    bw_scalar piv;
    int rank;
    mp_limb_t *inv;
    mp_limb_t *lambda;
    unsigned int *pivots_list;
    struct timeval tv;

    timer_r(&tv, TIMER_SET);

    if (pivlist)
	pivots_list = pivlist;
    else
	pivots_list = (unsigned int *) malloc(m * sizeof(int));

    inv = (mp_limb_t *) malloc(k_size * sizeof(mp_limb_t));
    lambda = (mp_limb_t *) malloc(k_size * sizeof(mp_limb_t));

    /* Pay attention here, this is a gaussian elimination on
     * *columns* */

    rank = 0;
    for (j = 0; j < m + n; j++) {
	jr = pi->clist[j];
	/* Find the pivot inside the column. */
	for (i = 0; i < m; i++) {
	    piv = mbmat_scal(e, i, jr);
	    if (!k_is_zero(piv))
		break;
	}
	if (i == m)
	    continue;
	assert(rank < m);
	pivots_list[rank++] = j;
	k_inv(inv, mbmat_scal(e, i, jr));
	/* Cancel this coeff in all other columns. */
	for (k = j + 1; k < m + n; k++) {
	    k_mul(lambda, mbmat_scal(e, i, pi->clist[k]), inv);
	    k_neg(lambda, lambda);
	    assert(delta[j] <= delta[k]);
	    if (!pivlist)
		ec_transvec(ec, j, k, i, lambda);
	    else
		e_transvec(e, pi->clist, j, k, i, lambda);
	    tp_transvec(pi, j, k, i, lambda);
	}
    }

    free(inv);
    free(lambda);

    if (rank != m) {
	fprintf(stderr, "Duh, rank is not == m !\n");
	exit(1);
    }

    for (j = 0; j < m; j++) {
	if (!pivlist)
	    ec_x_multiply(ec, pivots_list[j]);
	tp_x_multiply(pi, pivots_list[j]);
	delta[pivots_list[j]]++;
	global_sum_delta++;
    }

    if (!pivlist)
	free(pivots_list);

    *tt = timer_r(&tv, TIMER_ASK);
}
#endif/*}}}*/

#if 0/*{{{*/
void print_chance_list(unsigned int t, unsigned int *chance_list)
{
    /*
       printf("#RESULT T=%d\n", it->t);
       for(j=0;j< m + n ;j++) {
       if (chance_list[j]) printf("#RESULT J=%d\n", j);
       }
     */

    int j;
    printf("// step %d LOOK [", t);
    for (j = 0; j < m + n ; j++) {
        if (chance_list[j])
            printf(" %d", j);
    }
    printf(" ]\n");
}



bw_nbpoly	f_poly;
int		t_counter;
int	      * global_delta;
int		global_sum_delta;
unsigned int  * chance_list;
static int	recursion_level;	/* static, hence 0 on init */
static int	rec_threshold=1;
static int	print_min=10;
static int	preferred_quadratic_algorithm;	/* Defaults to 0 (old) */
static int	t_init;
static int	check_input = 1;
static int	check_pi = 1;


std::list<char *> ops_poly_args;

#define SAVE_LEVEL_THRESHOLD	4

void reclevel_prolog(void)
{
    int i;
    printf("%2d [ t=%6d ] ",recursion_level, t_counter);
    for(i=0;i<recursion_level;i++) printf("  ");
}

int sum_delta(int * delta)
{
    int i,res=0;
    for(i=0;i< m + n ;i++) res += delta[i];
    return res;
}

int max_delta(int * delta)
{
    int i,res=-1;
    for(i=0;i< m + n ;i++) if (delta[i]>res) res=delta[i];
    return res;
}


static int save_pi(struct t_poly * pi, int t_start, int t_middle, int t_end)
{
    char filename[FILENAME_LENGTH];
    FILE *f;
    int res;

    sprintf(filename, pi_meta_filename, t_start, t_end);
    f=fopen(filename,"w");
    if (f==NULL)
        return -1;
    res=tp_write(f,pi);
    if (fclose(f)<0)
        res=-1;

    if (res==0) {
        printf("Saved file %s\n",filename);
    } else {
        fprintf(stderr,"Failure to save %s: %s\n",
                filename,strerror(errno));
        return -1;
    }

    if (t_middle==-1)
        return res;

#if 0
    /* Unlinking not done, for safety */
    sprintf(filename, pi_meta_filename, t_start, t_middle);

    if (unlink(filename)<0) {
        fprintf(stderr,"Cannot unlink %s: %s\n",
                filename,strerror(errno));
    } else {
        printf("Unlinked %s\n",filename);
    }

    sprintf(filename, pi_meta_filename, t_middle, t_end);

    if (unlink(filename)<0) {
        fprintf(stderr,"Cannot unlink %s: %s\n",
                filename,strerror(errno));
    } else {
        printf("Unlinked %s\n",filename);
    }
#endif



    return res;
}

static void core_if_null(const void * p, const char * v)
{
    if (p==NULL) {
        printf("Could not allocate space for %s, dumping core\n",v);
        abort();
    }
}


static int retrieve_pi_files(struct t_poly ** p_pi, int t_start)
{
    DIR * pi_dir;
    struct dirent * curr;
    int n_pi_files;
    struct couple {int s; int e;} * pi_files;
    const char *pattern;
    int i;
    struct t_poly * left = NULL, * right = NULL;
    ft_order_t order;
    int o_i;
    struct dft_bb * dft_left, * dft_right, * dft_prod;
    double tt;

    *p_pi=NULL;

    if ((pi_dir=opendir("."))==NULL) {
        perror(".");
        return t_start;
    }

    printf("Scanning directory %s for pi files\n", ".");

    pattern = strrchr(pi_meta_filename,'/');
    if (pattern == NULL) {
        pattern = pi_meta_filename;
    } else {
        pattern++;
    }

    for(n_pi_files=0;(curr=readdir(pi_dir))!=NULL;) {
        int s,e;
        if (sscanf(curr->d_name,pattern,&s,&e)==2) {
            printf("Found %s\n", curr->d_name);
            if (s>e) {
                printf("but that's a stupid one\n");
                continue;
            }
            n_pi_files++;
        }
    }

    if (n_pi_files==0) {
        printf("Found no pi files\n");
        return t_start;
    }

    pi_files=(struct couple *) malloc(n_pi_files*sizeof(struct couple));

    rewinddir(pi_dir);

    for(i=0;i<n_pi_files && (curr=readdir(pi_dir))!=NULL;) {
        if (sscanf(curr->d_name,
                    pattern,
                    &(pi_files[i].s),
                    &(pi_files[i].e))==2)
        {
            if (pi_files[i].s>pi_files[i].e)
                continue;
            i++;
        }
    }
    n_pi_files=i;
    closedir(pi_dir);

    /* The rule is: only look at the best candidate. It's not worth
     * bothering about more subtle cases */
    for(;;) {
        int t_max;
        int best=-1;
        FILE *f;
        char filename[FILENAME_LENGTH];

        printf("Scanning for data starting at t=%d\n",t_start);
        t_max=-1;
        for(i=0;i<n_pi_files;i++) {
            if (pi_files[i].s==t_start && pi_files[i].e != -1) {
                printf("candidate : ");
                printf(pattern,pi_files[i].s,pi_files[i].e);
                printf("\n");
                if (pi_files[i].e>t_max) {
                    t_max=pi_files[i].e;
                    best=i;
                }
            }
        }
        if (t_max==-1) {
            printf("Could not find such data\n");
            break;
        }

        sprintf(filename,pi_meta_filename,t_start,t_max);
        printf("trying %s\n", filename);
        f=fopen(filename,"r");
        if (f==NULL) {
            perror(filename);
            pi_files[best].e=-1;
            continue;
        }
        /* Which degree can we expect for t_start..t_max ?
         */

        unsigned int pideg;
        pideg = iceildiv(m * (t_max - t_start), (m+n));
        pideg += 10;
        if (t_max > total_work) {
            pideg += t_max - total_work;
        }

        right=tp_read(f, pideg);
        fclose(f);

        if (right==NULL) {
            printf("%s : bad or nonexistent data\n",filename);
            pi_files[best].e=-1;
            continue;
        }

        if (left==NULL) {
            left=right;
            right=NULL;
            t_start=t_max;
            continue;
        }

        printf("Beginning multiplication\n");
        *p_pi=tp_comp_alloc(left,right);
        core_if_null(*p_pi,"*p_pi");

        order.set((*p_pi)->degree+1);
        o_i = order;

        dft_left=fft_tp_dft(left,order,&tt);
        printf("DFT(pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_left,"dft_left");

        dft_right=fft_tp_dft(right,order,&tt);
        printf("DFT(pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_right,"dft_right");

        dft_prod=fft_bbb_conv(dft_left,dft_right,&tt);
        printf("CONV(pi_left,pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_prod,"dft_prod");

        fft_tp_invdft(*p_pi,dft_prod,&tt);
        printf("IDFT(pi,%d) : %.2fs\n",o_i,tt);

        tp_free(left);
        tp_free(right);
        dft_bb_free(dft_left);
        dft_bb_free(dft_right);
        dft_bb_free(dft_prod);
        left=*p_pi;
        right=NULL;
        t_start=t_max;
    }
    free(pi_files);

    *p_pi=left;
    return t_start;
}

struct xcol_id {
    unsigned int pos;
    int deg;
    int sdeg;
};

static int col_cmp(const struct xcol_id * x, const struct xcol_id * y)
{
    int diff;
    int sdiff;
    diff = x->deg - y->deg;
    sdiff = x->sdeg - y->sdeg;
    return diff?diff:(sdiff?sdiff:(x->pos - y ->pos));
}

static int bw_check_chance(bw_mbmat e, unsigned int * clist)
{
    int i,j;
    unsigned int maxchance;

    maxchance=0;

    for(j=0;j< m + n ;j++) {
        for(i=0;i<m;i++) {
            bw_reduce_short_scalar(mbmat_scal(e,i,j));
        }
    }

    for(j=0;j< m + n ;j++) {
        if (mcol_is_zero(mbmat_col(e,clist?clist[j]:j)))
        {
            if (++chance_list[j] > maxchance)
                maxchance=chance_list[j];
            printf("Column %d happens to be zero ! (%d)\n",
                    j,chance_list[j]);
            if (t_counter < 30) {
                fprintf(stderr, "surely a degenerate case :-(\n");
                BUG();
            }
        } else
            chance_list[j]=0;
    }

    return maxchance;
}

static void banner_traditional(int t, int deg, double inner, double * last)
{
    *last+=inner;
    if ((t+1)%print_min==0 || t==deg) {
        reclevel_prolog();
        printf("avg=%.1f	"
                "step:	%.2fs	last %d : %.2fs\n",
                ((double)global_sum_delta)/(m+n),inner,
                1+(t%print_min),*last);
        *last=0.0;
    }
}
#endif/*}}}*/

#if 0/*{{{*/
static void compute_ctaf(bw_mbmat e,
        struct e_coeff * ec,
        struct t_poly * pi,
        int t,
        unsigned int * known_cols,
        double * tt)
{
    int i,j,k,s;
    int * bounds;
    struct timeval tv;

    timer_r(&tv,TIMER_SET);

    bounds=(int*) malloc((m+n)*sizeof(int));
    memset(bounds,0,(m+n)*sizeof(int));
    if (known_cols) for(j=0;j<m;j++) {
        /* 20070501 -- farewell to this 6+-year-old bug (!) */
        /* bounds[pi->clist[j]]=-1; */
        bounds[pi->clist[known_cols[j]]]=-1;
    }
    for(j=0;j< m + n ;j++) {
        if (bounds[j]==-1)
            continue;
        bounds[j]=pi->degnom[j];
        mcol_zero(mbmat_col(e,j));
    }

    for(j = 0 ; j < m + n ; j++) {
        for(s = 0 ; s <= bounds[j] ; s++) {
            for(i = 0 ; i < m; i++) {
                for(k = 0 ; k < m + n ; k++)	{
                    addmul( mbmat_scal(e,i,j),
                            mbmat_scal(mbpoly_coeff(ec->p,t-s),i,k),
                            bbmat_scal(bbpoly_coeff(pi->p,s),k,j));
                }
            }
        }
    }

    for(i = 0 ; i < m ; i++) {
        for(j = 0 ; j < m + n ; j++) {
            k_reduce(mbmat_scal(e,i,j));
        }
    }
    free(bounds);

    *tt=timer_r(&tv,TIMER_SET);
}
#endif/*}}}*/

/* There are two choices for the quadratic algorithm. The first one
 * (default) seems to be a bit faster. Once I've got a more reasonable
 * implementation, it'll be time to do a serious comparison.
 *
 * The complexity curve is steeper with the first version.
 */

#if 0/*{{{*/
static void bw_traditional_algo_1(struct e_coeff * ec, int * delta,
        struct t_poly * pi, int check_chance)
{
    unsigned int * perm;
    int t;
    bw_mbmat e;
    unsigned int * pivlist;
    double inner_1,inner_2,last;
    int k;

    perm=(unsigned int *) malloc((m+n)*sizeof(unsigned int));
    pivlist=(unsigned int *) malloc(m*sizeof(unsigned int));

    mbmat_alloc(e);
    mbmat_zero(e);

    assert(!ec_is_twisted(ec));
    last=0.0;

    for(t=0;t<=ec->degree;t++) {
        compute_ctaf(e,ec,pi,t,t?pivlist:NULL, &inner_1);
        column_order(perm,delta,pi);
        tp_apply_perm(pi,perm);
        if (check_chance)
            bw_check_chance(e,pi->clist);
        bw_gauss_onestep(e,ec,pi,delta,pivlist,	&inner_2);
        t_counter++;
        banner_traditional(t, ec->degree, inner_1 + inner_2, &last);
    }
    printf("DELTA : ( ");
    for(k=0;k< m + n ;k++) printf("%d ",delta[k]);
    printf(")\n");

    ec_advance(ec,ec->degree+1);	/* cosmetic */
    mbmat_free(e);
    free(pivlist);
    free(perm);
}

static void bw_traditional_algo_2(struct e_coeff * ec, int * delta,
        struct t_poly * pi, int check_chance)
{
    unsigned int * perm;
    int t;
    int deg;
    double inner,last;

    perm=(unsigned int *) malloc((m+n)*sizeof(unsigned int));

    assert(!ec_is_twisted(ec));

    deg=ec->degree;
    last=0.0;

    for(t=0;t<=deg;t++) {
        if (check_chance)
            bw_check_chance(mbpoly_coeff(ec->p,0),ec->clist);
        column_order(perm,delta,pi);
        tp_apply_perm(pi,perm);
        ec_apply_perm(ec,perm);
        bw_gauss_onestep(mbpoly_coeff(ec->p,0),ec,pi,delta,NULL,&inner);
        ec_advance(ec,1);
        t_counter++;
        banner_traditional(t, deg, inner, &last);
    }

    free(perm);
}
#endif/*}}}*/

static ulong extract_coeff_degree_t(uint t, ulong const * a, uint da, ulong const * b, uint db)
{
    ulong c = 0;
    /*
    0 <= s <= t
    0 <= s <= na
    0 <= t-s <= nb
    t-nb <= s <= t
    */
    uint high = std::min(t, da);
    uint low = std::max(0u, t - db);
    for(uint s = low ; s <= high ; s++) {
        uint si = s / ULONG_BITS;
        uint ss = s % ULONG_BITS;
        uint ri = (t-s) / ULONG_BITS;
        uint rs = (t-s) % ULONG_BITS;
        c ^= a[si] >> ss & b[ri] >> rs;
    }
    return c & 1UL;
}


static void extract_coeff_degree_t(uint t, uint piv[], polmat const& PI)
{
    using namespace std;
    using namespace globals;
    vector<bool> known(m+n,false);
    if (piv != NULL) {
        for(uint i = 0 ; i < m ; i++)
            known[piv[i]] = true;
    }

    vector<uint> z;

    for(uint j = 0 ; j < m+n ; j++) {
        if (known[j]) continue;
        e0.zcol(j);
        ulong some_nonzero = 0;
        for(uint i = 0 ; i < m ; i++) {
            ulong c = 0;
            for(uint k = 0 ; k < m + n ; k++) {
                c ^= extract_coeff_degree_t(t,
                        E.poly(i, k), E.ncoef + 1,
                        PI.poly(k, j), PI.deg(j));
            }
            e0.addcoeff(i,j,c);
            some_nonzero |= c;
        }
        if (!some_nonzero) {
            z.push_back(j);
        }
    }
    vector<uint> ncha(m + n, 0);
    if (z.empty()) {
        chance_list.swap(ncha);
    } else {
        for(uint i = 0 ; i < z.size() ; i++) {
            ++ncha[z[i]];
        }
        std::cout << fmt("Step %, zero cols:") % t;
        for(uint i = 0 ; i < m + n ; i++) {
            if (ncha[i] == 0) {
                chance_list[i] = 0;
            } else {
                uint w = chance_list[i] += ncha[i];
                std::cout << " " << i;
                if (w >= 2) {
                    std::cout << fmt("[%]") %w;
                }
            }
        }
        std::cout << "\n";
    }
}

/*
 * Rule : in the following, ec can be trashed at will.
 *
 * Compute pi_left, of degree (ec->degree + 1)*(m/(m+n)), such that ec * pi is
 * divisible by X^(ec->degree + 1) (that is, all coefficients up to
 * degree ec->degree in the product are forced to 0).
 *
 */

/* Forward declaration, it's used by the recursive version */
// static double bw_lingen(struct e_coeff *, int *, struct t_poly **);

static void go_quadratic(polmat & PI)
{
    using namespace globals;

    uint piv[m];
    uint deg = E.ncoef - 1;

    PI = polmat(m + n, m + n, deg * m / (m + n) + 10);
    for(uint i = 0 ; i < m + n ; i++) {
        PI.addcoeff(i,i,0,1UL);
        PI.deg(i) = 0;
    }

    rearrange_ordering(PI, NULL);

    for (uint dt = 0; dt <= deg ; dt++) {
        print_deltas();
	extract_coeff_degree_t(t - t0, dt ? piv : NULL, PI);
        gauss(piv, PI);
        rearrange_ordering(PI, piv);
        t++;
        // write_polmat(PI,std::string(fmt("pi-%-%") % t0 % t).c_str());
    }
    write_polmat(PI,std::string(fmt("pi-%-%") % t0 % t).c_str());
    print_deltas();
}

#if 0/*{{{*/
int check_zero_and_advance(struct e_coeff * ec, unsigned int kill)
{
    unsigned int i;
    for(i = 0 ; i < kill ; i++) {
        int res=1;
        int j;
        for(j=0;res && j< m + n ;j++) {
            if (!mcol_is_zero(mbmat_col(mbpoly_coeff(ec->p,0),j))) {
                die("argh, not zero !\n",1);
                return 1;
            }
        }
        if (!res) {
            return 0;
        }
        ec_advance(ec,1);
    }
    return 1;
}

static double bw_recursive_algorithm(struct e_coeff * ec,
        int * delta,
        struct t_poly ** p_pi)
{
    struct t_poly * pi_left, * pi_right;
    int deg,ldeg,rdeg;
    struct dft_mb *dft_e_left,*dft_e_middle;
    struct dft_bb *dft_pi_left,*dft_pi_right;
    struct dft_bb *dft_pi;
    ft_order_t sub_order;
    int so_i;
    int expected_pi_deg;
    struct timeval tv;
    double	t_dft_e_l,  t_dft_pi_l, t_conv_e, t_idft_e,
                t_dft_pi_r, t_conv_pi,  t_idft_pi, t_ft, t_cv, t_sub;
    int kill;
    int t0 = t_counter;

    timer_r(&tv,TIMER_SET);

    deg=ec->degree;

    /* Repartition of the job:
     *
     *		left		right
     * deg==0	1		0	(never recursive)
     * deg==1	1		1
     * deg==2	2		1
     * deg==n	n/2 + 1		(n+1)/2
     * 
     * The figures are for the number of steps, each one corres-
     * ponding to a m/(m+n) increase of the average degree of pi.
     */

    ldeg=(deg   /2)+1;
    rdeg=(deg+1)/2;

    assert(ldeg && rdeg && ldeg + rdeg == deg + 1);

    /* We aim at computing ec * pi / X^ldeg. The degree of this
     * product will be
     *
     * ec->degree + pi->degree - ldeg
     *
     * (We are actually only interested in the low (ec->degree-ldeg)
     * degree part of the product, but the whole thing is required)
     *
     * The expected value of pi->degree is 
     * 	ceil(ldeg*m/(m+n))
     *
     * The probability that pi exceeds this expected degree
     * depends on the base field, but is actually low.
     * However, by the end of the computations, this does
     * happen because the degrees increase unevenly.
     *
     * The DFTs of e and pi can be computed using only the
     * number of points given above, *even if their actual
     * degree is higher*. The FFT routines need to have
     * provision for this.
     *
     * The number of points will then be the smallest power
     * of 2 above deg+ceil(ldeg*m/(m+n))-ldeg+1
     */

    expected_pi_deg = 10 + iceildiv(ldeg*m, (m+n));
#ifdef	HAS_CONVOLUTION_SPECIAL
    kill=ldeg;
#else
    kill=0;
#endif
    sub_order.set(deg+expected_pi_deg-kill+1);

    so_i = sub_order;

    dft_e_left	= fft_ec_dft(ec,sub_order,&t_dft_e_l);
    reclevel_prolog();
    printf("DFT(e,%d) : %.2fs\n", so_i, t_dft_e_l);
    core_if_null(dft_e_left,"dft_e_left");

    /* at this point E can be shrinked. */
    ec->degree	= ldeg - 1;
    t_sub=bw_lingen(ec,delta,&pi_left);

    if (t_counter < t0 + ldeg) {
        printf("Exceptional situation, small generator ; escaping\n");
        *p_pi = pi_left;
        dft_mb_free(dft_e_left);
        return timer_r(&tv,TIMER_ASK);
    }

    printf("deg(pi_l)=%d, bound is %d\n",pi_left->degree,expected_pi_deg);
    if (!sub_order.fits(deg+pi_left->degree-kill + 1)) {
        printf("Warning : pi grows above its expected degree...\n");
        printf("order %d , while :\n"
                "deg=%d\n"
                "deg(pi_left)=%d\n"
                "ldeg-1=%d\n"
                "hence, %d is too big\n",
                so_i,
                deg,pi_left->degree,ldeg - 1,
                deg+pi_left->degree-kill + 1);

        int n_exceptional = 0;
        for(int i = 0 ; i < m + n ; i++) {
            n_exceptional += chance_list[i] * m;
        }
        if (!n_exceptional) {
            die("This should only happen at the end of the computation\n",1);
        }
        *p_pi = pi_left;
        dft_mb_free(dft_e_left);
        return timer_r(&tv,TIMER_ASK);
    }

    dft_pi_left	= fft_tp_dft(pi_left,sub_order,&t_dft_pi_l);
    reclevel_prolog();
    printf("DFT(pi_l,%d) : %.2fs\n", so_i, t_dft_pi_l);
    core_if_null(dft_pi_left,"dft_pi_left");

#ifdef  HAS_CONVOLUTION_SPECIAL
    dft_e_middle = fft_mbb_conv_sp(dft_e_left,dft_pi_left,ldeg,&t_conv_e);
#else
    dft_e_middle = fft_mbb_conv(dft_e_left,dft_pi_left,&t_conv_e);
#endif
    reclevel_prolog();
    printf("CONV(e*pi_l,%d) : %.2fs\n", so_i, t_conv_e);
    core_if_null(dft_e_middle,"dft_e_middle");

    /* This is a special convolution in the sense that we
     * compute f(w)*g(w) / w^k for k=ldeg, since we are
     * interested in fg div X^k (we know fg mod X^k==0)
     */

    ec_park(ec);
    ec_untwist(ec);

    fft_mb_invdft(ec->p,dft_e_middle,deg - kill,&t_idft_e);
    ec->degree=deg-kill;

    check_zero_and_advance(ec, ldeg - kill);
    reclevel_prolog();
    printf("IDFT(e,%d) : %.2fs\n", (int) (dft_e_middle->order), t_idft_e);

    dft_mb_free(dft_e_middle);
    dft_mb_free(dft_e_left);

    assert(ec->degree==rdeg-1);

    t_sub+=bw_lingen(ec,delta,&pi_right);
    printf("deg(pi_r)=%d, bound is %d\n",pi_right->degree,expected_pi_deg);

    *p_pi		= tp_comp_alloc(pi_left,pi_right);
    core_if_null(*p_pi,"*p_pi");

    printf("deg(pi_prod)=%d (order %d)\n",(*p_pi)->degree, so_i);
    assert(sub_order.fits((*p_pi)->degree + 1));

    dft_pi_right	= fft_tp_dft(pi_right, sub_order, &t_dft_pi_r);
    reclevel_prolog();
    printf("DFT(pi_r,%d) : %.2fs\n", so_i, t_dft_pi_r);
    core_if_null(dft_pi_right,"dft_pi_right");

    dft_pi		= fft_bbb_conv(dft_pi_left,dft_pi_right, &t_conv_pi);
    reclevel_prolog();
    printf("CONV(pi_l*pi_r,%d) : %.2fs\n", so_i, t_conv_pi);
    core_if_null(dft_pi,"dft_pi");


    fft_tp_invdft(*p_pi,dft_pi,&t_idft_pi);
    reclevel_prolog();
    printf("IDFT(pi,%d) : %.2fs\n",(int)(dft_pi->order),t_idft_pi);

    dft_bb_free(dft_pi);
    dft_bb_free(dft_pi_right);
    dft_bb_free(dft_pi_left);
    tp_free(pi_left);
    tp_free(pi_right);

    reclevel_prolog();

    t_ft=t_dft_e_l+t_dft_pi_l+t_idft_e+t_dft_pi_r+t_idft_pi;
    t_cv=t_conv_e+t_conv_pi;

    printf("proper : %.2fs (%.2fs FT + %.2fs CV), sub : %.2fs\n",
            t_ft+t_cv,t_ft,t_cv,t_sub);
    printf("constants : c_ft=%.4e c_cv=%.4e		# %d,%d\n",
            t_ft/(double)(so_i<<so_i),
            t_cv/(double)(1<<so_i),
            deg,so_i);
    printf("Different values for M1:");
    printf("   e_left: M1=%.3e\n",
            t_dft_e_l/ (so_i<<so_i) / (m*(m+n)));
    printf("  pi_left: M1=%.3e\n",
            t_dft_pi_l/(so_i<<so_i) / ((m+n)*(m+n)));
    printf("    e_inv: M1=%.3e\n",
            t_idft_e/  (so_i<<so_i) / (m*(m+n)));
    printf(" pi_right: M1=%.3e\n",
            t_dft_pi_r/(so_i<<so_i) / ((m+n)*(m+n)));
    printf("   pi_inv: M1=%.3e\n",
            t_idft_pi/ (so_i<<so_i) / ((m+n)*(m+n)));
    printf("   e_conv: M1=%.3e\n",
            t_conv_e/ (1<<so_i)  / (m*(m+n)*(m+n)));
    printf("  pi_conv: M1=%.3e\n",
            t_conv_pi/ (1<<so_i) /  ((m+n)*(m+n)*(m+n)));
    return timer_r(&tv,TIMER_ASK);
}

static double bw_lingen(struct e_coeff * ec,
        int * delta,
        struct t_poly ** p_pi)
{
    /* int check_chance; */
    double inner;
    int deg;
    int t_before;
    int did_rec=0;

    t_before=t_counter;
    deg=ec->degree;
    reclevel_prolog();
    printf("Degree %d\n",deg);

    recursion_level++;
    /* check_chance=(t_counter > total_work - 5); */
    if (ec->degree < rec_threshold) {
        /* check_chance=(t_counter > total_work - 5 - ec->degree); */
        inner=bw_traditional_algorithm(ec,delta,p_pi,1);
    } else {
        did_rec=1;
        inner=bw_recursive_algorithm(ec,delta,p_pi);
    }
    recursion_level--;

    reclevel_prolog();
    printf("Degree %d : took %.2fs\n",deg,inner);

    if (recursion_level <= SAVE_LEVEL_THRESHOLD) {
        if (!did_rec || recursion_level == SAVE_LEVEL_THRESHOLD) {
            save_pi(*p_pi,t_before,-1,t_counter);
        } else {
            int t_middle;
            t_middle=t_before+(deg/2)+1;
            save_pi(*p_pi,t_before,t_middle,t_counter);
        }
    }

    return inner;
}

void showuse(void)
{
    die("Usage : bw-master <bank#>\n",1);
}
#endif/*}}}*/

void recycle_old_pi(polmat& pi)
{
#if 0/*{{{*/
    if (pi_left!=NULL && !(!check_input && ec->degree<new_t-t_counter)) {
        int dg_kill;
        struct dft_bb * dft_pi_left;
        struct dft_mb * dft_e_left, * dft_e_middle;
        ft_order_t order;
        int o_i;

        printf("Beginning multiplication (e_left*pi_left)\n");

        /* That's a pity to bring the DFT so far, but we have to
         * check the input is correct */

        dg_kill=new_t-t_counter;

#ifdef	HAS_CONVOLUTION_SPECIAL
        if (check_input) {
            order.set(1+pi_left->degree+ec->degree);
        } else {
            order.set(1+pi_left->degree+ec->degree-dg_kill);
        }
#else
        order.set(1+pi_left->degree+ec->degree);
#endif
        o_i = order;

        dft_e_left=fft_ec_dft(ec,order,&tt);
        printf("DFT(e_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_e_left,"dft_e_left");

        dft_pi_left=fft_tp_dft(pi_left,order,&tt);
        printf("DFT(pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_left,"dft_pi_left");

        ec_untwist(ec);
#ifdef  HAS_CONVOLUTION_SPECIAL
        if (check_input) {
            dft_e_middle=fft_mbb_conv_sp(dft_e_left,dft_pi_left,
                    0,&tt);
        } else {
            dft_e_middle=fft_mbb_conv_sp(dft_e_left,dft_pi_left,
                    dg_kill,&tt);
            ec_advance(ec,dg_kill);
        }
#else
        dft_e_middle=fft_mbb_conv(dft_e_left,dft_pi_left,&tt);
#endif

        printf("CONV(e_left,pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_e_middle,"dft_e_middle");

#ifdef  HAS_CONVOLUTION_SPECIAL
        if (check_input) {
            fft_mb_invdft(ec->p,dft_e_middle,
                    ec->degree,&tt);
            printf("IDFT(e,%d) : %.2fs\n",o_i,tt);
            printf("Verifying the product\n");
            check_zero_and_advance(ec, dg_kill);
            printf("Input OK\n");
        } else {
            fft_mb_invdft(ec->p,dft_e_middle,
                    ec->degree-dg_kill,&tt);
            printf("IDFT(e,%d) : %.2fs\n",o_i,tt);
        }
#else
        fft_mb_invdft(ec->p,dft_e_middle,
                ec->degree,&tt);
        printf("IDFT(e,%d) : %.2fs\n",o_i,tt);
        printf("Verifying the product\n");
        check_zero_and_advance(ec, dg_kill);
        printf("Input OK\n");
#endif
        tp_act_on_delta(pi_left,global_delta);

        dft_bb_free(dft_pi_left);
        dft_mb_free(dft_e_left);
        dft_mb_free(dft_e_middle);
        t_counter=new_t;

        if (check_pi)
            save_pi(pi_left,t_start,-1,t_counter);
    }
    if (pi_left!=NULL && !check_input && ec->degree<new_t-t_counter) {
        printf("We are not interested in the computation of e(X)\n");
        ec_advance(ec,new_t-t_counter);
        tp_act_on_delta(pi_left,global_delta);
        t_counter=new_t;
        if (check_pi)
            save_pi(pi_left,t_start,-1,t_counter);
    }
#endif/*}}}*/
}

void block_wiedemann(void)
{
    bw_init();

    polmat pi_left, pi_right;

#if 0/*{{{*/
    if (check_pi) {
        if (retrieve_pi_files(&pi_left)) {
            recycle_old_pi(pi_left);
        }
    } else {
        printf("Not reading pi files due to --nopi option\n");
    }
#endif/*}}}*/

    using namespace globals;
    using namespace std;

    cout << fmt("E: % coeffs, t=%\n") % E.ncoef % t;

    e0 = bmat(m,m + n);
    for(uint i = 0 ; i < m + n ; i++) {
        E.deg(i) = E.ncoef - 1;
    }

    // E.resize(deg + 1);
    go_quadratic(pi_left);

    polmat F;
    compute_final_F_from_PI(F, pi_left);
    bw_commit_f(F);

	printf("// step %d LOOK [", t);
	for (uint j = 0; j < m + n; j++) {
		if (chance_list[j])
			printf(" %d", j);
	}
	printf(" ]\n");

#if 0/*{{{*/
    if (ec->degree>=0) {
        bw_lingen(ec,global_delta,&pi_right);
    } else {
        pi_right=pi_left;
        pi_left=NULL;
    }

    if (pi_left!=NULL) {
        struct dft_bb * dft_pi_left, * dft_pi_right, * dft_pi_prod;
        ft_order_t order;
        int o_i;

        printf("Beginning multiplication (pi_left*pi_right)\n");
        pi_prod=tp_comp_alloc(pi_left,pi_right);
        core_if_null(pi_prod,"pi_prod");

        order.set(pi_prod->degree+1);
        o_i = order;

        dft_pi_left=fft_tp_dft(pi_left,order,&tt);
        printf("DFT(pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_left,"dft_pi_left");

        dft_pi_right=fft_tp_dft(pi_right,order,&tt);
        printf("DFT(pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_right,"dft_pi_right");

        dft_pi_prod=fft_bbb_conv(dft_pi_left,dft_pi_right,&tt);
        printf("CONV(pi_left,pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_prod,"dft_pi_prod");

        fft_tp_invdft(pi_prod,dft_pi_prod,&tt);
        printf("IDFT(pi,%d) : %.2fs\n",o_i,tt);

        tp_free(pi_left);
        tp_free(pi_right);
        dft_bb_free(dft_pi_left);
        dft_bb_free(dft_pi_right);
        dft_bb_free(dft_pi_prod);
        pi_left=NULL;
        pi_right=NULL;
        if (check_pi) {
            save_pi(pi_prod,t_start,new_t,t_counter);
        }
        /* new_t is the inner value that has to be discarded */
    } else {
        pi_prod=pi_right;
        /* Don't save here, since it's already been saved by
         * lingen (or it comes from the disk anyway).
         */
    }
#endif/*}}}*/


    // print_chance_list(total_work, chance_list);
}

void usage()
{
    die("Usage: ./master2 -t <n> [options]\n",1);
}

int main(int argc, char *argv[])
{
    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    argv++, argc--;

    int pop = 0;

    using namespace globals;

    for( ; argc ; ) {
        if (strcmp(argv[0], "--subdir") == 0) {
            if (argc <= 1) usage();
            int rc = chdir(argv[1]);
            if (rc < 0) {
                perror(argv[1]);
                exit(errno);
            }
            argv+=2;
            argc-=2;
            continue;
        }
#if 0/*{{{*/
        if (strcmp(argv[0], "-t") == 0) {
            if (argc <= 1) usage();
            argv++, argc--;
            rec_threshold = atoi(argv[0]);
            argv++, argc--;
            continue;
        }
        if (strcmp(argv[0], "-p") == 0) {
            if (argc <= 1) usage();
            argv++, argc--;
            print_min = atoi(argv[0]);
            argv++, argc--;
            continue;
        }
        if (strcmp(argv[0], "--enable-complex-field") == 0
                || strcmp(argv[0], "--fermat-prime") == 0)
        {
            ops_poly_args.push_back(argv[0]);
            argv++, argc--;
            continue;
        }
        if (strcmp(argv[0], "--no-check-input") == 0) {
            check_input = 0;
            argv++, argc--;
            continue;
        }
        if (strcmp(argv[0], "--nopi") == 0) {
            check_pi = 0;
            argv++, argc--;
            continue;
        }
        if (strcmp(argv[0], "-q") == 0) {
            if (argc <= 1) usage();
            argv++, argc--;
            preferred_quadratic_algorithm = atoi(argv[0]);
            argv++, argc--;
            continue;
            /* changes the quadratic algorithm employed
               {"coppersmith old uneven",0},
               {"new even",1},
             */
        }
#endif/*}}}*/
        if (pop == 0) {
            read_mat_file_header(argv[0]);
            pop++;
            argv++, argc--;
            continue;
        } else if (pop == 1) {
            m = atoi(argv[0]);
            pop++;
            argv++, argc--;
            continue;
        } else if (pop == 2) {
            n = atoi(argv[0]);
            pop++;
            argv++, argc--;
            continue;
        }

        usage();
    }

    if (m == 0 || n == 0) {
        usage();
    }

    total_work = Lmacro(nrows, m, n);

    // if (rec_threshold == 1) usage();

    coredump_limit(1);

    /********************************************************************/

    block_wiedemann();

    // ft_order_t::cleanup();

    return 0;
}
/* vim: set sw=4 sta et: */
