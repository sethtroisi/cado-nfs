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
#include <utility>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <functional>

#include "cado.h"
#undef  ASSERT
#undef  ASSERT_ALWAYS
#include "utils.h"

#define Lmacro(N, m, n) (iceildiv((N)+2*(n),(m))+iceildiv((N)+2*(n),(n))+10)

unsigned int rec_threshold = 0;

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


// Requirements for the transform interface.
//
// - Must define an info struct[1]. This type may contain info
//   on the FFT size and so on, but NO DATA.
//
// - Must define a data_t typedef. This type is expected
//   to be as bare as possible (typically a pointer), and its
//   interpretation is not expected to be possible without the
//   accompanying fft_order_info_t struct.
//
// - Must define the functions as in the example below.

#include "gf2x.h"
struct fake_fft {
    unsigned int nc;
    unsigned int size;
    public:
    /* Set up the parameters for a multiplications of polynomials of
     * degree d1 and d2, into a polynomial of degree d3. The acc
     * parameter indicates the number of transform that we might occur to
     * add together.
     *
     * d3 might be different from d1 + d2, but not without help 
     */
    fake_fft(unsigned int d1,
            unsigned int d2,
            unsigned int d3,
            unsigned int acc)
    {
        if (d2 == UINT_MAX) d2 = d1;
        if (d3 == UINT_MAX) d3 = d1 + d2;
        nc = std::max(d1,d2) + 1;
        size = BITS_TO_WORDS(nc, ULONG_BITS);
    }

    ulong * alloc(unsigned int n) const {
        return new ulong[n * size];
    }
    void zero(ulong * p, unsigned int n = 1) const {
        memset(p, 0, n * size * sizeof(ulong));
    }
    void clear(ulong * p, unsigned int n) const {
        delete[] p;
    }
    ulong * get(ulong * p, unsigned int i) const {
        return p + size * i;
    }
    ulong const * cget(ulong * const p, unsigned int i) const {
        return p + size * i;
    }

    ulong * alloc_p(unsigned int n) const {
        return new ulong[n * 2 * size];
    }
    void zero_p(ulong * p, unsigned int n = 1) const {
        memset(p, 0, n * 2 * size * sizeof(ulong));
    }
    void clear_p(ulong * p, unsigned int n) const {
        delete[] p;
    }
    ulong * get_p(ulong * p, unsigned int i) const {
        return p + 2 * size * i;
    }
    ulong const * cget_p(ulong * const p, unsigned int i) const {
        return p + 2 * size * i;
    }

    void transform(ulong * dst, ulong const * src) const {
        memcpy(dst, src, size * sizeof(ulong));
    }
    void compose(ulong * dst, ulong const * s1, ulong const * s2) {
        mul_gf2x(dst, s1, size, s2, size);
    }
    void add(ulong * dst, ulong const * s1, ulong const * s2) {
        for(unsigned int i = 0 ; i < 2*size ; i++) {
            dst[i] = s1[i] ^ s2[i];
        }
    }
    void itransform(ulong * dst, ulong const * src) const {
        memcpy(dst, src, BITS_TO_WORDS(2*nc-1, ULONG_BITS) * sizeof(ulong));
    }
};

/* output: F0x is the candidate number x. All coordinates of the
 * candidate are grouped, and coefficients come in order (least
 * significant first), separated by a carriage return.
 */

// applies a permutation on the source indices
template<typename T>
void permute(std::vector<T>& a, uint p[])
{
    std::vector<T> b(a.size());
    for(uint i = 0 ; i < a.size() ; i++) {
        b[p[i]] = a[i];
    }
    a.swap(b);
}


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
 * (while uncritical, the computation of the initial E and the final F
 * may become a problem for debugging situations when the block size is
 * largish).
 *
 * all this stuff has linear complexity.
 */

template<typename POD>
inline void podswap(POD& a, POD& b) { POD c; c = a; a = b; b = c; } 

struct bcol;
struct bmat;
struct polmat;

struct bcol {/*{{{*/
    friend class bmat;
    uint nrows;
    inline uint stride() const { return BITS_TO_WORDS(nrows, ULONG_BITS); }
    private:
    ulong * x;
    public:
    bcol(uint nrows)
        : nrows(nrows), x(new ulong[stride()])
    {}
    bcol() { nrows = 0; x = NULL; }
    void clear() { nrows = 0; delete[] x; x = NULL; }
    private:
    bcol(bcol const& a) : nrows(a.nrows), x(a.x) {}
    bcol& operator=(bcol const& a) {
        nrows = a.nrows;
        x = a.x;
        return *this;
    }
    public:
    ~bcol() { delete[] x; }
    void swap(bcol& a) {
        podswap(nrows, a.nrows);
        podswap(x, a.x);
    }
    void zero() {
        memset(x, 0, stride() * sizeof(ulong));
    }
    bool is_zero_likely() const {
        for(uint i = 0 ; i < stride() ; i++) {
            if (x[i] != 0) return false;
        }
        return true;
    }
    inline ulong coeff(uint i) const
    {
        ASSERT(i < nrows);
        uint offset = i / ULONG_BITS;
        ulong shift = i % ULONG_BITS;
        return (x[offset] >> shift) & 1UL;
    }
    uint ffs() const {
        uint z = 0;
        const ulong * src = x;
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
    ulong * x;
    uint * order;
    inline uint stride() const { return BITS_TO_WORDS(nrows, ULONG_BITS); }
    public:
    bmat(uint nrows, uint ncols)
        : nrows(nrows), ncols(ncols), x(new ulong[ncols*stride()]),
        order(new uint[ncols])
    {
        memset(x, 0, ncols*stride()*sizeof(ulong));
        for(uint j = 0 ; j < ncols ; j++) order[j]=j;
    }
    bmat() {
        nrows = ncols = 0;
        x = NULL;
        order = NULL;
    }
    void clear() {
        nrows = ncols = 0;
        delete[] x; x = NULL;
        delete[] order; order = NULL;
    }
private:
    bmat(bmat const& a) : nrows(a.nrows), ncols(a.ncols), x(a.x), order(a.order) {}
    bmat& operator=(bmat const& a) {
        nrows = a.nrows;
        ncols = a.ncols;
        x = a.x;
        order = a.order;
        return *this;
    }
public:
    void swap(bmat& a) {
        podswap(nrows, a.nrows);
        podswap(ncols, a.ncols);
        podswap(x,a.x);
        podswap(order,a.order);
    }
    ~bmat() {
        delete[] x;
        delete[] order;
    };
#if 0
    bmat clone() const {
        bmat dst(nrows, ncols);
        memcpy(dst.x, x, ncols*stride()*sizeof(ulong));
        memcpy(dst.order, order, ncols*sizeof(uint));
        return dst;
    }
#endif
    /* zero column j */
    void zcol(uint j) {
        ASSERT(j < ncols);
        memset(x + order[j] * stride(), 0, stride() * sizeof(ulong));
    }
    bool col_is_zero(uint j) const {
        ASSERT(j < ncols);
        const ulong * src = x + order[j] * stride();
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
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        uint offset = i / ULONG_BITS;
        ulong shift = i % ULONG_BITS;
        return (x[order[j] * stride() + offset] >> shift) & 1UL;
    }
    /* add column j0 to column j. Optionally, b indicates a number of rows
     * which are known to be zero in both columns.
     * The mask corresponds to a multiplicating coefficient.
     */
    void acol(uint j, uint j0, uint b = 0, ulong mask = 1UL) {
        ASSERT(j < ncols);
        ASSERT(j0 < ncols);
        ulong * dst = x + order[j] * stride();
        const ulong * src = x + order[j0] * stride();
        for(uint l = b / ULONG_BITS ; l < stride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
    }
    void acol(uint j, bcol& a, uint b = 0, ulong mask = 1UL) {
        ASSERT(j < ncols);
        ulong * dst = x + order[j] * stride();
        const ulong * src = a.x;
        for(uint l = b / ULONG_BITS ; l < stride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
    }
    uint ffs(uint j) const {
        ASSERT(j < ncols);
        uint z = 0;
        const ulong * src = x + order[j] * stride();
        ulong w;
        uint k;
        for(k = stride() ; k && !(w=*src++) ; k--) z += ULONG_BITS;
        if (k == 0) return UINT_MAX;
        return z + __builtin_ctzl(w);
    }
    void extract_col(bcol& a, uint j) const
    {
        ASSERT(j < ncols);
        bcol tmp_a(nrows);
        memcpy(tmp_a.x, x + order[j] * stride(), stride() * sizeof(ulong));
        a.swap(tmp_a);
    }
    /* permute columns: column presently referred to with index i will
     * now be accessed at index p[i]
     */
    void perm(uint * p) {
        uint * norder(new uint[ncols]);
        for(uint i = 0 ; i < ncols ; i++) {
            norder[p[i]]=order[i];
        }
        podswap(order, norder);
    }
    void addcoeff(uint i, uint j, ulong z)
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
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
     * polynomials are stored packed, so this is not an array of bmat's.
     *
     * polmat may also represent other things, but for which speed is not
     * important.
     */
    uint nrows;
    uint ncols;
    uint ncoef;
    static bool critical;
    private:
    ulong * x;
    uint  * order;
    int   * _deg;
    inline uint stride() const { return BITS_TO_WORDS(ncoef, ULONG_BITS); }
    inline uint colstride() const { return nrows * stride(); }
    static void brev_warning();
    public:
    int& deg(uint j) { ASSERT(j < ncols); return _deg[order[j]]; }
    int const & deg(uint j) const { ASSERT(j < ncols); return _deg[order[j]]; }
    polmat(uint nrows, uint ncols, uint ncoef)
        : nrows(nrows), ncols(ncols), ncoef(ncoef),
        x(new ulong[ncols*colstride()]),
        order(new uint[ncols]),
        _deg(new int[ncols])
    {
        memset(x, 0, ncols*colstride()*sizeof(ulong));
        for(uint j = 0 ; j < ncols ; j++) order[j]=j;
        for(uint j = 0 ; j < ncols ; j++) _deg[j]=-1;
    }
    void clear() {
        nrows = ncols = ncoef = 0;
        delete[] x; x = NULL;
        delete[] order; order = NULL;
        delete[] _deg; _deg = NULL;
    }
    polmat() {
        nrows = ncols = ncoef = 0;
        x = NULL;
        order = NULL;
        _deg = NULL;
    }
    private:
    polmat(polmat const& a) { }
    polmat& operator=(polmat const&){ return *this;}
    public:
    ~polmat() {
        delete[] x;
        delete[] order;
        delete[] _deg;
    }
    inline void swap(polmat & n) {
        podswap(nrows,n.nrows);
        podswap(ncols,n.ncols);
        podswap(ncoef,n.ncoef);
        podswap(x,n.x);
        podswap(order,n.order);
        podswap(_deg,n._deg);
    }
    public:
#if 0
    polmat clone() {
        polmat dst(nrows, ncols, ncoef);
        memcpy(dst.x, x, ncols*colstride()*sizeof(ulong));
        memcpy(dst.order, order, ncols*sizeof(uint));
        memcpy(dst._deg, _deg, ncols*sizeof(uint));
        return dst;
    }
#endif
    /* this handles expansion and shrinking */
    void resize(uint ncoef2) {
        uint newstride = BITS_TO_WORDS(ncoef2, ULONG_BITS);
        if (newstride == stride()) {
            ncoef = ncoef2;
            return;
        }
        uint minstride = std::min(stride(),newstride);
        /* take the opportunity of reallocation for reordering columns. */
        polmat n(ncols,nrows,ncoef2);
        ulong * dst = n.x;
        for(uint j = 0 ; j < ncols ; j++) {
            const ulong * src = x + order[j] * colstride();
            for(uint i = 0 ; i < nrows ; i++) {
                memcpy(dst, src, minstride * sizeof(ulong));
                dst += newstride;
                src += stride();
            }
        }
        swap(n);
    }
    void xmul_col(uint j, uint s=1) {
        ASSERT(j < ncols);
        mp_limb_t * dst = x + order[j] * colstride();
        ASSERT(1u <= s && s <= GMP_LIMB_BITS-1);
        mpn_lshift(dst, dst, colstride(), s);
        /* We may have garbage for the low bits. It may matter, or maybe
         * not. If it does, call xclean0_col */
        deg(j) += deg(j) >= 0;
        // BUG_ON(deg(j) >= (int) ncoef);
        // we do NOT consider this a bug. Instead, this means that te
        // corresponding column has gone live. 
        // Normally not a problem since it is supposed to occur only
        // when sufficiently many generators are known.
    }
    void xmul_poly(uint i, uint j, uint s=1) {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        ulong * dst = x + (order[j] * nrows + i) * stride();
        ASSERT(1 <= s && s <= GMP_LIMB_BITS-1);
        mpn_lshift(dst, dst, stride(), s);
        /* We may have garbage for the low bits. It may matter, or maybe
         * not. If it does, call xclean0_col */
        // deg(j) += deg(j) >= 0;
    }
    void addpoly(uint i, uint j, const ulong * src) {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        ulong * dst = x + (order[j] * nrows + i) * stride();
        for(uint k = 0 ; k < stride() ; k++)
            dst[k] ^= src[k];
    }
    void xclean0_col(uint j) {
        ASSERT(j < ncols);
        ulong * dst = x + order[j] * colstride();
        for(uint i = 0 ; i < nrows ; i++, dst += stride())
            dst[0] &= ~1UL;
        deg(j) -= deg(j) == 0;
    }
    /* add column j times mask to column i. */
    void acol(uint i, uint j, ulong mask = 1UL) {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        ulong * dst = x + order[i] * colstride();
        const ulong * src = x + order[j] * colstride();
        for(uint l = 0 ; l < colstride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
        // ASSERT(_deg[order[i]] >= _deg[order[j]]);
        deg(i) = std::max(deg(j), deg(i));
    }
    /* permute columns: column presently referred to with index i will
     * now be accessed at index p[i]
     */
    void perm(uint * p) {
        uint * norder(new uint[ncols]);
        for(uint i = 0 ; i < ncols ; i++) {
            norder[p[i]]=order[i];
        }
        podswap(order,norder);
        delete[] norder;
    }
    uint ffs(uint j, uint k) const {
        ASSERT(k < nrows);
        ASSERT(j < ncols);
        uint offset = k / ULONG_BITS;
        ulong mask = 1UL << (k % ULONG_BITS);
        const ulong * src = x + j * colstride() + offset;
        for(uint z = 0 ; z < nrows ; z++) {
            if (*src & mask) return z;
            src += stride();
        }
        return UINT_MAX;
    }
    ulong * poly(uint i, uint j)
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return x + (order[j] * nrows + i) * stride();
    }
    ulong const * poly(uint i, uint j) const
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return x + (order[j] * nrows + i) * stride();
    }
    ulong * col(uint j) { return poly(0, j); }
    ulong const * col(uint j) const { return poly(0, j); }
    /* zero column j */
    void zcol(uint j) {
        ASSERT(j < ncols);
        memset(col(j), 0, colstride() * sizeof(ulong));
        deg(j) = -1;
    }
    /* shift is understood ``shift left'' (multiply by X) */
    void import_col_shift(uint k, polmat const& a, uint j, int s)
    {
        ASSERT(k < ncols);
        ASSERT(j < a.ncols);
        BUG_ON(a.nrows != nrows);
        ulong const * src = a.col(j);
        ulong * dst = col(k);
        /* Beware, trailing bits are lurking here and there */
        if (colstride() == a.colstride()) {
            /* Then this is fast */
            if (s > 0) {
                mp_size_t sw = s / GMP_LIMB_BITS;
                mp_size_t sb = s & (GMP_LIMB_BITS - 1);
                mpn_lshift(dst + sw, src, colstride() - sw, sb);
            } else {
                mp_size_t sw = (-s) / GMP_LIMB_BITS;
                mp_size_t sb = (-s) & (GMP_LIMB_BITS - 1);
                mpn_rshift(dst, src + sw, colstride() - sw, sb);
            }
        } else if (colstride() < a.colstride()) {
            /* Otherwise, we have to resample... */
            ulong tmp[a.colstride()];
            BUG_ON(critical);
            for(uint i = 0 ; i < nrows ; i++) {
                if (s > 0) {
                    mp_size_t sw = s / GMP_LIMB_BITS;
                    mp_size_t sb = s & (GMP_LIMB_BITS - 1);
                    mpn_lshift(tmp + sw, src, a.stride() - sw, sb);
                } else {
                    mp_size_t sw = (-s) / GMP_LIMB_BITS;
                    mp_size_t sb = (-s) & (GMP_LIMB_BITS - 1);
                    mpn_rshift(tmp, src + sw, a.stride() - sw, sb);
                }
                memcpy(dst, tmp, stride() * sizeof(ulong));
                src += a.stride();
                dst += stride();
            }
        } else {
            BUG();
            // I doubt we'll ever end up here, so I do not have a test
            // case. It should be ok though 
            BUG_ON(critical);
            for(uint i = 0 ; i < nrows ; i++) {
                if (s > 0) {
                    mp_size_t sw = s / GMP_LIMB_BITS;
                    mp_size_t sb = s & (GMP_LIMB_BITS - 1);
                    dst[stride()] = mpn_lshift(dst + sw, src, a.stride() -sw, sb);
                } else {
                    mp_size_t sw = (-s) / GMP_LIMB_BITS;
                    mp_size_t sb = (-s) & (GMP_LIMB_BITS - 1);
                    mpn_rshift(dst, src + sw, a.stride() - sw, sb);
                }
                src += a.stride();
                dst += stride();
            }
        }
    }

    /* these accessors are not the preferred ones for performance */
    inline ulong coeff(uint i, uint j, uint k) const
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        ASSERT(k < ncoef);
        BUG_ON(critical);
        brev_warning();
        uint offset = k / ULONG_BITS;
        ulong shift = k % ULONG_BITS;
        return poly(i,j)[offset] >> shift & 1UL;
    }
    void extract_coeff(bmat& a, uint k) const
    {
        ASSERT(k < ncoef);
        BUG_ON(critical);
        brev_warning();
        bmat tmp_a(nrows, ncols);
        for(uint j = 0 ; j < ncols ; j++) {
            for(uint i = 0 ; i < nrows ; i++) {
                uint boffset = i / ULONG_BITS;
                ulong bmask = 1UL << (i % ULONG_BITS);
                ulong coeff = this->coeff(i,j,k);
                tmp_a.x[j*tmp_a.stride() + boffset] ^= bmask & -coeff;
            }
        }
        a.swap(tmp_a);
    }
    void addcoeff(uint i, uint j, uint k, ulong z)
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        ASSERT(k < ncoef);
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
    bmat x;
    pa->extract_coeff(x, k);
    dbmat(&x);
}
/*}}}*/
template<typename fft_type> struct tpolmat /* {{{ */
{
    /* tpolmat is essentially the same as polmat, except that it is
     * expected to hold transform data. Only a few operations are
     * allowed.
     */
    uint nrows;
    uint ncols;
    fft_type o;
    private:
    /* XXX do we force ulong * as a type ? Not necessarily in the long
     * run, but it seems fair enough for the time being.
     */
    ulong * x;
    uint  * order;
    int   * _deg;
    // inline uint stride() const { return BITS_TO_WORDS(ncoef, ULONG_BITS); }
    // inline uint colstride() const { return nrows * stride(); }
    public:
    int& deg(uint j) { ASSERT(j < ncols); return _deg[order[j]]; }
    int const & deg(uint j) const { ASSERT(j < ncols); return _deg[order[j]]; }
    tpolmat(uint nrows, uint ncols, fft_type const& o)
        : nrows(nrows), ncols(ncols), o(o),
        x(o.alloc(nrows * ncols)),
        order(new uint[ncols]),
        _deg(new int[ncols])
    {
        o.zero(x, nrows * ncols);
        for(uint j = 0 ; j < ncols ; j++) order[j]=j;
        for(uint j = 0 ; j < ncols ; j++) _deg[j]=-1;
    }
    void clear() {
        nrows = ncols = 0;
        o.clear(x);
        delete[] order; order = NULL;
        delete[] _deg; _deg = NULL;
    }
    tpolmat() {
        nrows = ncols = 0;
        x = NULL;
        order = NULL;
        _deg = NULL;
    }
    private:
    tpolmat(tpolmat const& a) { }
    tpolmat& operator=(tpolmat const&){ return *this;}
    public:
    ~tpolmat() {
        o.clear(x);
        delete[] order;
        delete[] _deg;
    }
    inline void swap(tpolmat & n) {
        tpodswap(nrows,n.nrows);
        tpodswap(ncols,n.ncols);
        tpodswap(o,n.o);
        tpodswap(x,n.x);
        tpodswap(order,n.order);
        tpodswap(_deg,n._deg);
    }
    public:
#if 0
    tpolmat clone() {
        tpolmat dst(nrows, ncols, o);
        memcpy(dst.x, x, ncols*colstride()*sizeof(ulong));
        memcpy(dst.order, order, ncols*sizeof(uint));
        memcpy(dst._deg, _deg, ncols*sizeof(uint));
        return dst;
    }
#endif
    /* permute columns: column presently referred to with index i will
     * now be accessed at index p[i]
     */
    void perm(uint * p) {
        uint * norder(new uint[ncols]);
        for(uint i = 0 ; i < ncols ; i++) {
            norder[p[i]]=order[i];
        }
        podswap(order,norder);
        delete[] norder;
    }
    ulong * poly(uint i, uint j)
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return o.get(x, order[j] * nrows + i);
    }
    ulong const * poly(uint i, uint j) const
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return o.get(x, order[j] * nrows + i);
    }
    ulong * col(uint j) { return poly(0, j); }
    ulong const * col(uint j) const { return poly(0, j); }
    /* zero column j */
    void zcol(uint j) {
        ASSERT(j < ncols);
        o.zero(col(j));
        deg(j) = -1;
    }
};
/*}}}*/

namespace globals {
    uint nrows;
    uint ncols;
    uint m,n;
    uint t,t0,total_work;
    std::vector<uint> delta;
    std::vector<uint> chance_list;
    std::vector<uint> gone_mad;

    polmat E;
    bmat e0;

    // F0 is exactly the n x n identity matrix, plus the X^(s-exponent)e_{cnum}
    // vectors. Here we store the cnum,exponent pairs.
    std::vector<std::pair<uint, uint> > f0_data;

    double start_time;
}

// To multiply on the right an m x n matrix A by F0, we start by copying
// A into the first n columns. Since we're also dividing out by X^t0, the
// result has to be shifted t0 positions to the right.
// Afterwards, column n+j of the result is column cnum[j] of A, shifted
// exponent[j] positions to the right.
void compute_E_from_A(polmat const &a)
{
    using namespace globals;
    polmat tmp_E(n, m + n, a.ncoef - t0);
    for(uint j = 0 ; j < n ; j++) {
        tmp_E.import_col_shift(j, a, j, - (int) t0);
    }
    for(uint j = 0 ; j < m ; j++) {
        uint cnum = f0_data[j].first;
        uint exponent = f0_data[j].second;
        tmp_E.import_col_shift(n + j, a, cnum, - (int) exponent);
    }
    E.swap(tmp_E);
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
    polmat tmp_F(n, m + n, globals::t0 + pi.ncoef);
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
            if (gone_mad[j])
                continue;
            // We zero out the whole column, it's less trouble.
            tmpmat.zcol(0);
            for(uint k = 0 ; k < l.size() ; k++) {
                tmpmat.addpoly(l[k].second, 0, pi.poly(l[k].first, j));
            }
            tmp_F.addpoly(i,j,pi.poly(i,j));
            for(uint k = 0 ; k < exps.size() ; k++) {
                tmpmat.xmul_poly(exps[k],0,t0-exps[k]);
                tmp_F.addpoly(i,j,tmpmat.poly(exps[k],0));
            }
        }
    }
    for(uint j = 0 ; j < m + n ; j++) {
        if (gone_mad[j])
            continue;
        tmp_F.deg(j) = t0 + pi.deg(j);
    }
    F.swap(tmp_F);
}


const char *a_meta_filename = "A-%02d-%02d";
const char *f_meta_filename = "F%02d";
const char *pi_meta_filename = "pi-%d-%d";
const char *f_base_filename = "F_INIT_QUICK";


void read_data_for_series(polmat& A, int & rc)/*{{{*/
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

    rc = read_coeffs;
    A.swap(a);
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
    std::cout << fmt("recovered init data % on disk\n") % fn;
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

void bw_commit_f(polmat& F) /*{{{*/
{
    char filename[FILENAME_LENGTH];

    using namespace globals;
    using namespace std;

    for (uint j = 0; j < m + n ; j++) {
        sprintf(filename, f_meta_filename, j);
        if (gone_mad[j]) {
            cout << fmt("not writing % -- gone haywire since step %, "
                    "hence % steps ago\n")
                % filename % (t-gone_mad[j]) % gone_mad[j] ;
            continue;
        }
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
        bmat amat;
        A.extract_coeff(amat, k);
        for(uint j=0;r < m && j<n;j++) {
            /* copy column j, coeff k into acol */
            bcol acol;
            amat.extract_col(acol, j);
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
                if (k * n > m + 40) {
                    printf("The choice of starting vectors was bad. "
                            "Cannot find %u independent cols within A\n",m);
                    exit(1);
                }
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
            pcols[r].swap(acol);
            r++;

            if (r == m)
                printf("[X^%d] A, col %d increases rank to %d (head row %d)\n",
                        k,j,r,u);
        }
    }

    t0 = exponent[r-1] + 1;
    printf("Found satisfying init data for t0=%d\n", t0);
                    
    /*
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
    */


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
    uint last = UINT_MAX;
    uint nrep = 0;
    for(uint i = 0 ; i < m + n ; i++) {
        uint d = delta[i];
        if (d == last) {
            nrep++;
            continue;
        }
        // Flush the pending repeats
        if (last != UINT_MAX) {
            std::cout << " " << last;
            if (nrep > 1)
                std::cout << "[" << nrep << "]";
        }
        last = d;
        nrep = 1;
    }
    BUG_ON(last == UINT_MAX);
    std::cout << " " << last;
    if (nrep > 1)
        std::cout << "[" << nrep << "]";
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
    int read_coeffs;
    polmat A;
    using namespace globals;

    printf("Using A(X) div X in order to consider Y as starting point\n");
    total_work--;


    printf("Reading scalar data in polynomial ``a''\n");
    read_data_for_series(A, read_coeffs);


    /* informational *//*{{{ */
    if (read_coeffs < (int) total_work) {
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
    gone_mad.assign(m + n, 0);
    chance_list.assign(m + n, 0);

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
    uint p[m+n];
    for(uint i = 0 ; i < m + n ; i++) {
        p[corresp[i].second] = i;
    }
    permute(delta, p);
    permute(chance_list, p);
    permute(gone_mad, p);
    if (piv) {
        for(uint i = 0 ; i < m ; i++) {
            piv[i] = p[piv[i]];
        }
    }
    PI.perm(p);
    e0.perm(p);
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
        if (PI.deg(j) >= (int) PI.ncoef) {
            if (gone_mad[j] == 0) {
                std::cout << fmt("Column % goes haywire from step %\n")%j%t;
            }
            gone_mad[j]++;
            PI.deg(j) = PI.ncoef-1;
        }
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
    fscanf(f, "// %u ROWS %u COLUMNS ; MODULUS %s",
            &globals::nrows, &globals::ncols, modulus);
    BUG_ON(strcmp(modulus,"2") != 0);
    BUG_ON(globals::nrows > globals::ncols);

    fclose(f);
}

#if 0/* {{{ */
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
                c ^= extract_coeff_degree_t(t - t0,
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
static void compute_lingen(polmat& pi);

static void go_quadratic(polmat& pi)
{
    using namespace globals;

    uint piv[m];
    uint deg = E.ncoef - 1;

    polmat tmp_pi(m + n, m + n, deg * m / (m + n) + 10);
    for(uint i = 0 ; i < m + n ; i++) {
        tmp_pi.addcoeff(i,i,0,1UL);
        tmp_pi.deg(i) = 0;
    }

    rearrange_ordering(tmp_pi, NULL);


    for (uint dt = 0; dt <= deg ; dt++) {
        double delta;
        delta = seconds() - start_time;

        double percent = (double) dt / (deg + 1);
        percent = percent * percent;

        double estim_final = delta / percent;

        percent *= 100.0;
            
        printf("%5.0f / est %-7.0f (%2.0f%%) ",
                delta, estim_final, percent);

        print_deltas();
	extract_coeff_degree_t(t, dt ? piv : NULL, tmp_pi);
        gauss(piv, tmp_pi);
        rearrange_ordering(tmp_pi, piv);
        t++;
        // if (t % 60 < 10 || deg-dt < 30)
            // write_polmat(tmp_pi,std::string(fmt("tmp_pi-%-%") % t0 % t).c_str());
    }
    pi.swap(tmp_pi);
    write_polmat(pi,std::string(fmt("pi-%-%") % t0 % t).c_str());
    print_deltas();
}


typedef fake_fft fft_type;
static void go_recursive(polmat& pi)
{
    using namespace globals;
    uint deg = E.ncoef - 1;
    uint ldeg = deg / 2 + 1;
    uint rdeg = (deg + 1) / 2;
    assert(ldeg && rdeg && ldeg + rdeg == deg + 1);

    uint expected_pi_deg = 10 + iceildiv(ldeg*m, (m+n));
    uint kill;

#ifdef	HAS_CONVOLUTION_SPECIAL
    kill=ldeg;
#else
    kill=0;
#endif

    fft_type o(deg, expected_pi_deg, deg + expected_pi_deg - kill, m + n);
    // tpolmat<fft_type> E_hat(m+n, m+n, o);
    go_quadratic(pi);
}

static void compute_lingen(polmat& pi)
{
    /* reads the data in the global thing, E and delta. ;
     * compute the linear generator from this.
     */
    using namespace globals;
    uint deg = E.ncoef - 1;

    if (deg <= rec_threshold) {
        go_quadratic(pi);
    } else {
        go_recursive(pi);
    }
}

#if 0/*{{{*/
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

    { bmat tmp_e0(m,m + n); e0.swap(tmp_e0); }

    for(uint i = 0 ; i < m + n ; i++) {
        E.deg(i) = E.ncoef - 1;
    }

    start_time = seconds();

    // E.resize(deg + 1);
    polmat pi_left;
    compute_lingen(pi_left);

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
        if (strcmp(argv[0], "-t") == 0) {
            if (argc <= 1) usage();
            argv++, argc--;
            rec_threshold = atoi(argv[0]);
            argv++, argc--;
            continue;
        }
#if 0/*{{{*/
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

    total_work = Lmacro(ncols, m, n);

    // if (rec_threshold == 1) usage();

    coredump_limit(1);

    /********************************************************************/

    block_wiedemann();

    // ft_order_t::cleanup();

    return 0;
}
/* vim: set sw=4 sta et: */
