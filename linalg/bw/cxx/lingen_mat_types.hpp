#ifndef LINGEN_MAT_TYPES_HPP_
#define LINGEN_MAT_TYPES_HPP_

#include "manu.h"

#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include <algorithm>

/* See the discussion in lingen_binary about the pros and cons of data
 * ordering schemes */

template<typename POD>/*{{{*/
inline void podswap(POD& a, POD& b) { POD c; c = a; a = b; b = c; } /*}}}*/

struct bcol;
struct bmat;
struct polmat;

struct bcol {/*{{{*/
    friend class bmat;
    unsigned int nrows;
    inline unsigned int stride() const { return BITS_TO_WORDS(nrows, ULONG_BITS); }
    private:
    unsigned long * x;
    public:
    bcol(unsigned int nrows)
        : nrows(nrows), x(new unsigned long[stride()])
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
        memset(x, 0, stride() * sizeof(unsigned long));
    }
    bool is_zero_likely() const {
        for(unsigned int i = 0 ; i < stride() ; i++) {
            if (x[i] != 0) return false;
        }
        return true;
    }
    inline unsigned long coeff(unsigned int i) const
    {
        ASSERT(i < nrows);
        unsigned int offset = i / ULONG_BITS;
        unsigned long shift = i % ULONG_BITS;
        return (x[offset] >> shift) & 1UL;
    }
    unsigned int ffs() const {
        unsigned int z = 0;
        const unsigned long * src = x;
        unsigned long w;
        unsigned int k;
        for(k = stride() ; k && !(w=*src++) ; k--) z += ULONG_BITS;
        if (k == 0) return UINT_MAX;
        return z + ctzl(w);
    }
    void add(bcol const& a, unsigned long mask = 1UL) {
        for(unsigned int l = 0 ; l < stride() ; l++) {
            x[l] ^= a.x[l] & -mask;
        }
    }
};/*}}}*/
struct bmat { /* This is for example small-e, an m times b matrix *//*{{{*/
    friend class polmat;
    unsigned int nrows;
    unsigned int ncols;
    private:
    unsigned long * x;
    unsigned int * order;
    inline unsigned int stride() const { return BITS_TO_WORDS(nrows, ULONG_BITS); }
    public:
    bmat(unsigned int nrows, unsigned int ncols)
        : nrows(nrows), ncols(ncols), x(new unsigned long[ncols*stride()]),
        order(new unsigned int[ncols])
    {
        memset(x, 0, ncols*stride()*sizeof(unsigned long));
        for(unsigned int j = 0 ; j < ncols ; j++) order[j]=j;
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
        memcpy(dst.x, x, ncols*stride()*sizeof(unsigned long));
        memcpy(dst.order, order, ncols*sizeof(unsigned int));
        return dst;
    }
#endif
    /* zero column j */
    void zcol(unsigned int j) {
        ASSERT(j < ncols);
        memset(x + order[j] * stride(), 0, stride() * sizeof(unsigned long));
    }
    bool col_is_zero(unsigned int j) const {
        ASSERT(j < ncols);
        const unsigned long * src = x + order[j] * stride();
        for(unsigned int i = 0 ; i < stride() ; i++) {
            if (src[i] != 0) return false;
        }
        return true;
    }
    bool is_zero() const {
        for(unsigned int i = 0 ; i < ncols * stride() ; i++) {
            if (x[i] != 0) return false;
        }
        return true;
    }
    inline unsigned long coeff(unsigned int i, unsigned int j) const
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        unsigned int offset = i / ULONG_BITS;
        unsigned long shift = i % ULONG_BITS;
        return (x[order[j] * stride() + offset] >> shift) & 1UL;
    }
    /* add column j0 to column j. Optionally, b indicates a number of rows
     * which are known to be zero in both columns.
     * The mask corresponds to a multiplicating coefficient.
     */
    void acol(unsigned int j, unsigned int j0, unsigned int b = 0, unsigned long mask = 1UL) {
        ASSERT(j < ncols);
        ASSERT(j0 < ncols);
        unsigned long * dst = x + order[j] * stride();
        const unsigned long * src = x + order[j0] * stride();
        for(unsigned int l = b / ULONG_BITS ; l < stride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
    }
    void acol(unsigned int j, bcol& a, unsigned int b = 0, unsigned long mask = 1UL) {
        ASSERT(j < ncols);
        unsigned long * dst = x + order[j] * stride();
        const unsigned long * src = a.x;
        for(unsigned int l = b / ULONG_BITS ; l < stride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
    }
    unsigned int ffs(unsigned int j) const {
        ASSERT(j < ncols);
        unsigned int z = 0;
        const unsigned long * src = x + order[j] * stride();
        unsigned long w;
        unsigned int k;
        for(k = stride() ; k && !(w=*src++) ; k--) z += ULONG_BITS;
        if (k == 0) return UINT_MAX;
        return z + ctzl(w);
    }
    void extract_col(bcol& a, unsigned int j) const
    {
        ASSERT(j < ncols);
        bcol tmp_a(nrows);
        memcpy(tmp_a.x, x + order[j] * stride(), stride() * sizeof(unsigned long));
        a.swap(tmp_a);
    }
    /* permute columns: column presently referred to with index i will
     * now be accessed at index p[i]
     */
    void perm(unsigned int * p) {
        unsigned int * norder(new unsigned int[ncols]);
        for(unsigned int i = 0 ; i < ncols ; i++) {
            norder[p[i]]=order[i];
        }
        podswap(order, norder);
    }
    void addcoeff(unsigned int i, unsigned int j, unsigned long z)
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        unsigned int offset = i / ULONG_BITS;
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
    unsigned int nrows;
    unsigned int ncols;
    unsigned int ncoef;
    static bool critical;
    private:
    unsigned long * x;
    unsigned int  * order;
    int   * _deg;
    inline unsigned int stride() const { return BITS_TO_WORDS(ncoef, ULONG_BITS); }/*{{{*/
    inline unsigned int colstride() const { return nrows * stride(); }/*}}}*/
    public:
    int& deg(unsigned int j) { ASSERT(j < ncols); return _deg[order[j]]; }/*{{{*/
    int deg(unsigned int j) const { ASSERT(j < ncols); return _deg[order[j]]; }
    int maxdeg() const {
        int m = -1;
        for(unsigned int j = 0 ; j < ncols ; j++) {
            if (_deg[order[j]] > m) m = _deg[order[j]];
        }
        return m;
    }/*}}}*/
    /* ctors dtors etc {{{ */
    polmat(unsigned int nrows, unsigned int ncols, unsigned int ncoef)
        : nrows(nrows), ncols(ncols), ncoef(ncoef),
        x(new unsigned long[ncols*colstride()]),
        order(new unsigned int[ncols]),
        _deg(new int[ncols])
    {
        memset(x, 0, ncols*colstride()*sizeof(unsigned long));
        for(unsigned int j = 0 ; j < ncols ; j++) order[j]=j;
        for(unsigned int j = 0 ; j < ncols ; j++) _deg[j]=-1;
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
    /* }}} */
    public:
#if 0
    polmat clone() {
        polmat dst(nrows, ncols, ncoef);
        memcpy(dst.x, x, ncols*colstride()*sizeof(unsigned long));
        memcpy(dst.order, order, ncols*sizeof(unsigned int));
        memcpy(dst._deg, _deg, ncols*sizeof(unsigned int));
        return dst;
    }
#endif
    /* this handles expansion and shrinking */
    void resize(unsigned int ncoef2) {/*{{{*/
        unsigned int newstride = BITS_TO_WORDS(ncoef2, ULONG_BITS);
        if (newstride == stride()) {
            ncoef = ncoef2;
            return;
        }
        unsigned int minstride = std::min(stride(),newstride);
        /* take the opportunity of reallocation for reordering columns. */
        polmat n(nrows,ncols,ncoef2);
        unsigned long * dst = n.x;
        for(unsigned int j = 0 ; j < ncols ; j++) {
            const unsigned long * src = x + order[j] * colstride();
            for(unsigned int i = 0 ; i < nrows ; i++) {
                memcpy(dst, src, minstride * sizeof(unsigned long));
                dst += newstride;
                src += stride();
            }
        }
        swap(n);
    }/*}}}*/
    /* Divide by X^k, keep ncoef2 coefficients *//*{{{*/
    void xdiv_resize(unsigned int k, unsigned int ncoef2) {
        ASSERT(k + ncoef2 <= ncoef);
        ASSERT(k);
        unsigned int newstride = BITS_TO_WORDS(ncoef2, ULONG_BITS);
        if (newstride == stride()) {
            /*
            mp_size_t sw = k / GMP_LIMB_BITS;
            mp_size_t sb = k & (GMP_LIMB_BITS - 1);
            mpn_rshift(x, x + sw, ncols * colstride() - sw, sb);
            */
            ASSERT(k && k < GMP_LIMB_BITS);
            mpn_rshift(x, x, ncols * colstride(), k);
            ncoef = ncoef2;
            return;
        }
        /* take the opportunity of reallocation for reordering columns. */
        polmat n(nrows,ncols,ncoef2);

        ASSERT(GMP_LIMB_BITS == ULONG_BITS);
        mp_size_t sw = k / GMP_LIMB_BITS;
        mp_size_t sb = k & (GMP_LIMB_BITS - 1);
        unsigned int input_length = BITS_TO_WORDS(k + ncoef2, ULONG_BITS) - sw;

        /* sb might equal zero here, in which case sw >0. So k is
         * sw * ULONG_BITS, and input_length == newstride
         */

        /* We know that input_length <= stride */
        if (newstride < input_length) {
            ASSERT(sb && sb < GMP_LIMB_BITS);
            /* {{{ In this case we're encumbered by the fact that
             * mpn_rshift might land outside our main area for the last
             * cell. So we'll use an ugly hack */

            unsigned int i = 0 ;
            unsigned int j = 0 ;
            unsigned long * dst = n.poly(0,0);
            for(j = 0 ; j < ncols-1 ; j++) {
                const unsigned long * src = poly(0,j);
                for(i = 0 ; i < nrows ; i++) {
                    mpn_rshift(dst, src + sw, input_length, sb);
                    dst += newstride;
                    src += stride();
                }
            }

            {
                const unsigned long * src = poly(0,j);
                for(i = 0 ; i < nrows - 1 ; i++) {
                    mpn_rshift(dst, src + sw, input_length, sb);
                    dst += newstride;
                    src += stride();
                }

                unsigned long * tmp = new unsigned long[input_length];
                mpn_rshift(tmp, src + sw, input_length, sb);
                memcpy(dst, tmp, newstride * sizeof(unsigned long));
                delete[] tmp;
            }
            /*}}}*/
        } else {
            unsigned long * dst = n.poly(0,0);
            if (sb == 0) {
                ASSERT(input_length == newstride);
                for(unsigned int j = 0 ; j < ncols ; j++) {
                    const unsigned long * src = poly(0,j);
                    for(unsigned int i = 0 ; i < nrows ; i++) {
                        memcpy(dst, src + sw,
                                newstride * sizeof(unsigned long));
                        dst += newstride;
                        src += stride();
                    }
                }
            } else {
                ASSERT(sb && sb < GMP_LIMB_BITS);
                for(unsigned int j = 0 ; j < ncols ; j++) {
                    const unsigned long * src = poly(0,j);
                    for(unsigned int i = 0 ; i < nrows ; i++) {
                        mpn_rshift(dst, src + sw, input_length, sb);
                        dst += newstride;
                        src += stride();
                    }
                }
            }
        }
        swap(n);
        for(unsigned int j = 0 ; j < ncols ; j++) {
            if (deg(j)) deg(j) -= k;
        }
    }
    /*}}}*/
    void xmul_col(unsigned int j, unsigned int s=1) {/*{{{*/
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
    }/*}}}*/
    void xmul_poly(unsigned int i, unsigned int j, unsigned int s=1) {/*{{{*/
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        unsigned long * dst = x + (order[j] * nrows + i) * stride();
        ASSERT(1 <= s && s <= GMP_LIMB_BITS-1);
        mpn_lshift(dst, dst, stride(), s);
        /* We may have garbage for the low bits. It may matter, or maybe
         * not. If it does, call xclean0_col */
        // deg(j) += deg(j) >= 0;
    }/*}}}*/
    void addpoly(unsigned int i, unsigned int j, const unsigned long * src) {/*{{{*/
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        unsigned long * dst = x + (order[j] * nrows + i) * stride();
        for(unsigned int k = 0 ; k < stride() ; k++)
            dst[k] ^= src[k];
    }/*}}}*/
    void xclean0_col(unsigned int j) {/*{{{*/
        ASSERT(j < ncols);
        unsigned long * dst = x + order[j] * colstride();
        for(unsigned int i = 0 ; i < nrows ; i++, dst += stride())
            dst[0] &= ~1UL;
        deg(j) -= deg(j) == 0;
    }/*}}}*/
    /* {{{ add column j times mask to column i. */
    void acol(unsigned int i, unsigned int j, unsigned long mask = 1UL) {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        unsigned long * dst = x + order[i] * colstride();
        const unsigned long * src = x + order[j] * colstride();
        for(unsigned int l = 0 ; l < colstride() ; l++) {
            dst[l] ^= src[l] & -mask;
        }
        // ASSERT(_deg[order[i]] >= _deg[order[j]]);
        deg(i) = std::max(deg(j), deg(i));
    }
    /* }}} */
    /* permute columns: column presently referred to with index i will
     * now be accessed at index p[i]
     */
    void perm(unsigned int * p) {/*{{{*/
        unsigned int * norder(new unsigned int[ncols]);
        for(unsigned int i = 0 ; i < ncols ; i++) {
            norder[p[i]]=order[i];
        }
        podswap(order,norder);
        delete[] norder;
    }/*}}}*/
    unsigned int ffs(unsigned int j, unsigned int k) const {/*{{{*/
        ASSERT(k < nrows);
        ASSERT(j < ncols);
        unsigned int offset = k / ULONG_BITS;
        unsigned long mask = 1UL << (k % ULONG_BITS);
        const unsigned long * src = x + j * colstride() + offset;
        for(unsigned int z = 0 ; z < nrows ; z++) {
            if (*src & mask) return z;
            src += stride();
        }
        return UINT_MAX;
    }/*}}}*/
    unsigned int valuation() const {/*{{{*/
        /* It's fairly ugly, but not critical */
        unsigned long y[stride() + 1];
        for(unsigned int k = 0 ; k < stride() ; k++) {
            y[k] = 0;
        }
        y[stride()] = 1;
        for(unsigned int j = 0 ; j < ncols ; j++) {
            const unsigned long * src = poly(0,j);
            for(unsigned int i = 0 ; i < nrows ; i++) {
                for(unsigned int k = 0 ; k < stride() ; k++) {
                    y[k] |= src[k];
                }
                src += stride();
            }
        }
        return mpn_scan1(y,0);
    }/*}}}*/
    void setdeg(unsigned int j) { /* {{{ */
        unsigned long y[stride()];
        for(unsigned int k = 0 ; k < stride() ; k++) {
            y[k] = 0;
        }
        const unsigned long * src = poly(0,j);
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int k = 0 ; k < stride() ; k++) {
                y[k] |= src[k];
            }
            src += stride();
        }
        int k;
        /* Keep only (ULONG_BITS-ncoef) mod ULONG BITS in the top word. */
        y[stride()-1] &= ~0UL >> ((-ncoef) & (ULONG_BITS-1));
        for(k = stride() - 1 ; k >= 0 && y[k] == 0 ; k--);
        if (k < 0) {
            deg(j) = -1;
        } else {
            deg(j) = (k+1) * ULONG_BITS - clzl(y[k]) - 1;
        }
    } /* }}} */
    unsigned long * poly(unsigned int i, unsigned int j)/*{{{*/
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return x + (order[j] * nrows + i) * stride();
    }
    unsigned long const * poly(unsigned int i, unsigned int j) const
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return x + (order[j] * nrows + i) * stride();
    }/*}}}*/
    unsigned long * col(unsigned int j) { return poly(0, j); }/*{{{*/
    unsigned long const * col(unsigned int j) const { return poly(0, j); }/*}}}*/
    /* zero column j *//*{{{*/
    void zcol(unsigned int j) {
        ASSERT(j < ncols);
        memset(col(j), 0, colstride() * sizeof(unsigned long));
        deg(j) = -1;
    }/*}}}*/
    /* is zero column j *//*{{{*/
    bool is_zcol(unsigned int j) const {
        ASSERT(j < ncols);
        unsigned long const * y = col(j);
        for(unsigned int x = colstride(); x--;y++)
            if (*y) return false;
        return true;
    }/*}}}*/
    /* shift is understood ``shift left'' (multiply by X) */
    void import_col_shift(unsigned int k, polmat const& a, unsigned int j, int s)/*{{{*/
    {
        ASSERT(k < ncols);
        ASSERT(j < a.ncols);
        BUG_ON(a.nrows != nrows);
        unsigned long const * src = a.col(j);
        unsigned long * dst = col(k);
        unsigned int as = s > 0 ? s : -s;
        mp_size_t sw = as / GMP_LIMB_BITS;
        mp_size_t sb = as & (GMP_LIMB_BITS - 1);
        /* Beware, trailing bits are lurking here and there */
        if (colstride() == a.colstride()) {
            /* Then this is fast *//*{{{*/
            if (sb) {
                if (s > 0) {
                    mpn_lshift(dst + sw, src, colstride() - sw, sb);
                } else {
                    mpn_rshift(dst, src + sw, colstride() - sw, sb);
                }
            } else {
                if (s > 0) {
                    memcpy(dst + sw, src, (colstride()-sw)*sizeof(mp_limb_t));
                } else {
                    memcpy(dst, src + sw, (colstride()-sw)*sizeof(mp_limb_t));
                }
            }/*}}}*/
        } else if (colstride() < a.colstride()) {
            /* Otherwise, we have to resample... *//*{{{*/
            unsigned long tmp[a.colstride()];
            BUG_ON(critical);
            for(unsigned int i = 0 ; i < nrows ; i++) {
                if (sb) {
                    if (s > 0) {
                        mpn_lshift(tmp + sw, src, a.stride() - sw, sb);
                    } else {
                        mpn_rshift(tmp, src + sw, a.stride() - sw, sb);
                    }
                } else {
                    if (s > 0) {
                        memcpy(tmp + sw, src,(a.stride()-sw)*sizeof(mp_limb_t));
                    } else {
                        memcpy(tmp, src + sw,(a.stride()-sw)*sizeof(mp_limb_t));
                    }
                }
                memcpy(dst, tmp, stride() * sizeof(unsigned long));
                src += a.stride();
                dst += stride();
            }/*}}}*/
        } else {
            BUG();/*{{{*/
            // I doubt we'll ever end up here, so I do not have a test
            // case. It should be ok though 
            BUG_ON(critical);
            for(unsigned int i = 0 ; i < nrows ; i++) {
                if (sb) {
                    if (s > 0) {
                        dst[stride()] = mpn_lshift(dst + sw, src, a.stride() -sw, sb);
                    } else {
                        mpn_rshift(dst, src + sw, a.stride() - sw, sb);
                    }
                } else {
                    if (s > 0) {
                        memcpy(dst+sw,src,(a.stride()-sw)*sizeof(mp_limb_t));
                    } else {
                        memcpy(dst,src+sw,(a.stride()-sw)*sizeof(mp_limb_t));
                    }
                }
                src += a.stride();
                dst += stride();
            }/*}}}*/
        }
    }/*}}}*/
    /* these accessors are not the preferred ones for performance */
    inline unsigned long coeff(unsigned int i, unsigned int j, unsigned int k) const/*{{{*/
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        ASSERT(k < ncoef);
        BUG_ON(critical);
        // brev_warning();
        unsigned int offset = k / ULONG_BITS;
        unsigned long shift = k % ULONG_BITS;
        return poly(i,j)[offset] >> shift & 1UL;
    }/*}}}*/
    void extract_coeff(bmat& a, unsigned int k) const/*{{{*/
    {
        ASSERT(k < ncoef);
        BUG_ON(critical);
        // brev_warning();
        bmat tmp_a(nrows, ncols);
        for(unsigned int j = 0 ; j < ncols ; j++) {
            for(unsigned int i = 0 ; i < nrows ; i++) {
                unsigned int boffset = i / ULONG_BITS;
                unsigned long bmask = 1UL << (i % ULONG_BITS);
                unsigned long coeff = this->coeff(i,j,k);
                tmp_a.x[j*tmp_a.stride() + boffset] ^= bmask & -coeff;
            }
        }
        a.swap(tmp_a);
    }/*}}}*/
    void addcoeff(unsigned int i, unsigned int j, unsigned int k, unsigned long z)/*{{{*/
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        ASSERT(k < ncoef);
        BUG_ON(critical);
        // brev_warning();
        unsigned int offset = k / ULONG_BITS;
        poly(i,j)[offset] ^= z << (k % ULONG_BITS);
    }/*}}}*/
};

bool polmat::critical = false;

/*}}}*/
template<typename fft_type> struct tpolmat /* {{{ */
{
    /* tpolmat is essentially the same as polmat, except that it is
     * expected to hold transform data. Only a few operations are
     * allowed.
     */
    unsigned int nrows;
    unsigned int ncols;
    fft_type o;
    private:
    typename fft_type::t x;
    unsigned int  * order;
    int   * _deg;
    // inline unsigned int stride() const { return BITS_TO_WORDS(ncoef, ULONG_BITS); }
    // inline unsigned int colstride() const { return nrows * stride(); }
    public:
    int& deg(unsigned int j) { ASSERT(j < ncols); return _deg[order[j]]; }
    int const & deg(unsigned int j) const { ASSERT(j < ncols); return _deg[order[j]]; }
    tpolmat(unsigned int nrows, unsigned int ncols, fft_type const& o)
        : nrows(nrows), ncols(ncols), o(o),
        x(o.alloc(nrows * ncols)),
        order(new unsigned int[ncols]),
        _deg(new int[ncols])
    {
        for(unsigned int j = 0 ; j < ncols ; j++) order[j]=j;
        for(unsigned int j = 0 ; j < ncols ; j++) _deg[j]=-1;
    }
    void zero()
    {
        o.zero(x, nrows * ncols);
    }
    void clear() {
        nrows = ncols = 0;
        o.free(x, nrows * ncols); x = NULL;
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
        o.free(x,  nrows * ncols);
        delete[] order;
        delete[] _deg;
    }
    inline void swap(tpolmat & n) {
        podswap(nrows,n.nrows);
        podswap(ncols,n.ncols);
        podswap(o,n.o);
        podswap(x,n.x);
        podswap(order,n.order);
        podswap(_deg,n._deg);
    }
    public:
#if 0
    tpolmat clone() {
        tpolmat dst(nrows, ncols, o);
        memcpy(dst.x, x, ncols*colstride()*sizeof(unsigned long));
        memcpy(dst.order, order, ncols*sizeof(unsigned int));
        memcpy(dst._deg, _deg, ncols*sizeof(unsigned int));
        return dst;
    }
#endif
    /* permute columns: column presently referred to with index i will
     * now be accessed at index p[i]
     */
    void perm(unsigned int * p) {
        unsigned int * norder(new unsigned int[ncols]);
        for(unsigned int i = 0 ; i < ncols ; i++) {
            norder[p[i]]=order[i];
        }
        podswap(order,norder);
        delete[] norder;
    }
    typename fft_type::t poly(unsigned int i, unsigned int j)
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return o.get(x, order[j] * nrows + i);
    }
    typename fft_type::src_t poly(unsigned int i, unsigned int j) const
    {
        ASSERT(i < nrows);
        ASSERT(j < ncols);
        return o.get(x, order[j] * nrows + i);
    }
    typename fft_type::t col(unsigned int j) { return poly(0, j); }
    typename fft_type::src_t col(unsigned int j) const { return poly(0, j); }
    /* zero column j */
    void zcol(unsigned int j) {
        ASSERT(j < ncols);
        o.zero(col(j), 1);
        deg(j) = -1;
    }
};
/*}}}*/

template<typename fft_type>
void transform(tpolmat<fft_type>& dst, polmat& src, fft_type& o, int d)
{
    tpolmat<fft_type> tmp(src.nrows, src.ncols, o);
    tmp.zero();
    for(unsigned int j = 0 ; j < src.ncols ; j++) {
        for(unsigned int i = 0 ; i < src.nrows ; i++) {
            o.dft(tmp.poly(i,j), src.poly(i,j), d);
        }
    }
    dst.swap(tmp);
}

/* XXX Do Strassen here ! */
template<typename fft_type>
void compose(
        tpolmat<fft_type>& dst,
        tpolmat<fft_type> const & s1,
        tpolmat<fft_type> const & s2,
        fft_type& o)
{
    tpolmat<fft_type> tmp(s1.nrows, s2.ncols, o);
    ASSERT(s1.ncols == s2.nrows);
    tmp.zero();
    typename fft_type::t x = o.alloc(1);
    for(unsigned int j = 0 ; j < s2.ncols ; j++) {
        for(unsigned int k = 0 ; k < s1.ncols ; k++) {
            for(unsigned int i = 0 ; i < s1.nrows ; i++) {
                o.compose(x, s1.poly(i,k), s2.poly(k,j));
                o.add(tmp.poly(i,j), tmp.poly(i,j), x);
            }
        }
    }
    o.free(x,1);
    dst.swap(tmp);
}

template<typename fft_type>
void itransform(polmat& dst, tpolmat<fft_type>& src, fft_type& o, int d)
{
    polmat tmp(src.nrows, src.ncols, d + 1);
    for(unsigned int j = 0 ; j < src.ncols ; j++) {
        for(unsigned int i = 0 ; i < src.nrows ; i++) {
            o.ift(tmp.poly(i,j), d, src.poly(i,j));
        }
        tmp.setdeg(j);
    }
    dst.swap(tmp);
}


#endif	/* LINGEN_MAT_TYPES_HPP_ */
