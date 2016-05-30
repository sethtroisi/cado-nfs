#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <pthread.h>
#include <vector>
#include <utility>
#include <sstream>
#include <algorithm>
#include <map>

#include "bwc_config.h"
#include "matmul.h"
#include "matmul-common.h"
#include "mpfq_layer.h"
#include "portability.h"

#include "matmul_facade.h"

#include "arith-modp.hpp"

typedef arith_modp::gfp<sizeof(abelt)/sizeof(unsigned long)> gfp;

/* This extension is used to distinguish between several possible
 * implementations of the product */
#define MM_EXTENSION   "-zone"
#define MM_MAGIC_FAMILY        0xb002UL
#define MM_MAGIC_VERSION       0x1001UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)

/* This selects the default behaviour as to which is our best code
 * for multiplying. If this flag is 1, then a multiplication matrix times
 * vector (direction==1) performs best if the in-memory structure
 * reflects the non-transposed matrix. Similarly, a vector times matrix
 * multiplication (direction==0) performs best if the in-memory structure
 * reflects the transposed matrix. When the flag is 1, the converse
 * happens.
 * This flag depends on the implementation, and possibly even on the cpu
 * type under certain circumstances.
 */
#define MM_DIR0_PREFERS_TRANSP_MULT   1

using namespace std;

static inline uint64_t cputicks()
{
        uint64_t r;
        __asm__ __volatile__(
                "rdtsc\n\t"
                "shlq $32, %%rdx\n\t"
                "orq %%rdx, %%rax\n\t"
                : "=a"(r)
                :
                : "rdx");
        return r;
}

/* {{{ some boilerplate related to saving cache files */

struct cachefile {
    FILE * f;
    cachefile(FILE * f) : f(f) {}
    template<typename T> struct basic_seq {
        cachefile& out(cachefile& c, T const * t, size_t n) const {
            for(size_t i = 0 ; i < n ; i++) c << t[i];
            return c;
        }
        cachefile& in(cachefile& c, T * t, size_t n) {
            for(size_t i = 0 ; i < n ; i++) c >> t[i];
            return c;
        }
    };
    template<typename T> struct seq : public basic_seq<T> {
    };
};

template<> struct cachefile::seq<uint16_t> : public cachefile::basic_seq<uint16_t> {
    typedef uint16_t T;
    cachefile& in(cachefile& c, T * t, size_t n) {
        MATMUL_COMMON_READ_MANY16(t, n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        MATMUL_COMMON_WRITE_MANY16(t, n, c.f);
        return c;
    }
};

template<> struct cachefile::seq<uint32_t> : public cachefile::basic_seq<uint32_t> {
    typedef uint32_t T;
    cachefile& in(cachefile& c, T * t, size_t n) {
        MATMUL_COMMON_READ_MANY32(t, n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        MATMUL_COMMON_WRITE_MANY32(t, n, c.f);
        return c;
    }
};

template<> struct cachefile::seq<pair<int, uint32_t>> : public cachefile::basic_seq<pair<int, uint32_t>> {
    typedef pair<int, uint32_t> T;
    cachefile& in(cachefile& c, T * t, size_t n) {
        MATMUL_COMMON_READ_MANY32(t, 2*n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        MATMUL_COMMON_WRITE_MANY32(t, 2*n, c.f);
        return c;
    }
};

template<> struct cachefile::seq<pair<uint16_t, int32_t>> : public cachefile::basic_seq<pair<uint16_t, int32_t>> {
    typedef pair<uint16_t, int32_t> T;
    cachefile& in(cachefile& c, T * t, size_t n) {
        MATMUL_COMMON_READ_MANY16(t, 3*n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        MATMUL_COMMON_WRITE_MANY16(t, 3*n, c.f);
        return c;
    }
};

template<> struct cachefile::seq<pair<uint16_t, uint16_t>> : public cachefile::basic_seq<pair<uint16_t, uint16_t>> {
    typedef pair<uint16_t, uint16_t> T;
    cachefile& in(cachefile& c, T * t, size_t n) {
        MATMUL_COMMON_READ_MANY16(t, 2*n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        MATMUL_COMMON_WRITE_MANY16(t, 2*n, c.f);
        return c;
    }
};

template<> struct cachefile::seq<pair<int, pair<uint32_t, int32_t>>> : public cachefile::basic_seq<pair<int, pair<uint32_t, int32_t>>> {
    typedef pair<int, pair<uint32_t, int32_t>> T;
    cachefile& in(cachefile& c, T * t, size_t n) {
        MATMUL_COMMON_READ_MANY32(t, 3*n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        MATMUL_COMMON_WRITE_MANY32(t, 3*n, c.f);
        return c;
    }
};

template<> struct cachefile::seq<pair<uint16_t, pair<uint16_t, int32_t>>> : public cachefile::basic_seq<pair<uint16_t, pair<uint16_t, int32_t>>> {
    typedef pair<uint16_t, pair<uint16_t, int32_t>> T;
    cachefile& in(cachefile& c, T * t, size_t n) {
        MATMUL_COMMON_READ_MANY16(t, 4*n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        MATMUL_COMMON_WRITE_MANY16(t, 4*n, c.f);
        return c;
    }
};


cachefile& operator>>(cachefile & c, uint32_t & x) {
    MATMUL_COMMON_READ_ONE32(x, c.f);
    return c;
}

cachefile& operator<<(cachefile & c, uint32_t const & x) {
    MATMUL_COMMON_WRITE_ONE32(x, c.f);
    return c;
}

cachefile& operator>>(cachefile & c, int32_t & x) {
    MATMUL_COMMON_READ_ONE32(x, c.f);
    return c;
}

cachefile& operator<<(cachefile & c, int32_t const & x) {
    MATMUL_COMMON_WRITE_ONE32(x, c.f);
    return c;
}

cachefile& operator>>(cachefile & c, uint64_t & x) {
    MATMUL_COMMON_READ_ONE64(x, c.f);
    return c;
}

cachefile& operator<<(cachefile & c, uint64_t const & x) {
    MATMUL_COMMON_WRITE_ONE64(x, c.f);
    return c;
}

cachefile& operator>>(cachefile & c, int64_t & x) {
    MATMUL_COMMON_READ_ONE64(x, c.f);
    return c;
}

cachefile& operator<<(cachefile & c, int64_t const & x) {
    MATMUL_COMMON_WRITE_ONE64(x, c.f);
    return c;
}

template<typename T> cachefile& operator>>(cachefile & c, cachefile::seq<T> & s);
template<typename T> cachefile& operator<<(cachefile & c, cachefile::seq<T> const & s);

template<typename T> cachefile& operator>>(cachefile & c, vector<T>&v)
{
    size_t size;
    c >> size;
    v.insert(v.end(), size, T());
    return cachefile::seq<T>().in(c, &(v[0]), v.size());
}
template<typename T> cachefile& operator<<(cachefile & c, vector<T> const &v)
{
    c << v.size();
    return cachefile::seq<T>().out(c, &(v[0]), v.size());
}
/* }}} */

struct zone {/*{{{*/
    typedef vector<pair<uint16_t, uint16_t>> qpm_t;
    typedef vector<pair<uint16_t, pair<uint16_t, int32_t>>> qg_t;
    unsigned int i0, j0;
    qpm_t qp, qm;
    qg_t qg;
    zone() { i0 = j0 = 0; }
    zone(unsigned int i0, unsigned int j0) : i0(i0), j0(j0) {}
    inline bool empty() const { return qp.empty() && qm.empty() && qg.empty(); }
    inline size_t size() const { return qp.size() + qm.size() + qg.size(); }
    void operator()(gfp::elt_ur *, const gfp::elt *) const;

    struct sort_qpm {
        inline bool operator()(qpm_t::value_type const& a, qpm_t::value_type const& b) const {
            return a.second < b.second;
        }
    };

    struct sort_qg {
        inline bool operator()(qg_t::value_type const& a, qg_t::value_type const& b) const {
            return a.second.first < b.second.first;
        }
    };

    void sort() {
        std::sort(qp.begin(), qp.end(), sort_qpm());
        std::sort(qm.begin(), qm.end(), sort_qpm());
        std::sort(qg.begin(), qg.end(), sort_qg());
    }
    struct rowmajor_sorter {/*{{{*/
        inline bool operator()(zone const& a, zone const& b) const {
            return a.i0 < b.i0 || (a.i0 == b.i0 && a.j0 < b.j0);
        }
    };/*}}}*/
    struct colmajor_sorter {/*{{{*/
        inline bool operator()(zone const& a, zone const& b) const {
            return a.j0 < b.j0 || (a.j0 == b.j0 && a.i0 < b.i0);
        }
    };/*}}}*/

    cachefile& cachefile_load(cachefile& c) {/*{{{*/
        return c >> i0 >> j0 >> qp >> qm >> qg;
    }/*}}}*/
    cachefile& cachefile_save(cachefile& c) const {/*{{{*/
        return c << i0 << j0 << qp << qm << qg;
    }/*}}}*/
};

inline cachefile& operator>>(cachefile& c, zone& z) { return z.cachefile_load(c); }
inline cachefile& operator<<(cachefile& c, zone const & z) { return z.cachefile_save(c); }
/*}}}*/

struct matmul_zone_data {/*{{{*/
    /* repeat the fields from the public interface */
    struct matmul_public_s public_[1];
    /* now our private fields */
    abdst_field xab;

    vector<zone> q;

    /* {{{ timing data */
    struct twn {
        uint64_t tt;
        size_t w;
        unsigned int n;
        // twn() : tt(0), w(0), n(0) {}
    };

    typedef map<unsigned int, twn> tmap_t;
    tmap_t tmap;
    /* }}} */

    /* {{{ dispatchers and combiners */

    struct dispatcher_t : public vector<uint16_t> {
        typedef vector<uint16_t> super;
        unsigned int i0;
        unsigned int j0;
        /* a dispatcher contains: (col id)* */
        cachefile& cachefile_load(cachefile& c) {/*{{{*/
            return c >> i0 >> j0 >> (super&)*this;
        }/*}}}*/
        cachefile& cachefile_save(cachefile& c) const {/*{{{*/
            return c << i0 << j0 << (super&)*this;
        }/*}}}*/
    };
    struct combiner_t {
        unsigned int i0;
        unsigned int j0;
        /* which buffer will we read from, and how many values, from
         * where ?  */
        size_t index, offset, count;
        /* a combiner contains: (dest row id)* (duplicated indices for negative) */
        typedef uint16_t main_value_type;
        vector<main_value_type> main;
        vector<pair<uint16_t, int32_t>> aux;
        struct rowmajor_sorter {
            inline bool operator()(combiner_t const& a, combiner_t const& b) const {
                return a.i0 < b.i0 || (a.i0 == b.i0 && a.j0 < b.j0);
            }
        };
        cachefile& cachefile_load(cachefile& c) {/*{{{*/
            return c >> i0 >> j0
                >> index >> offset >> count
                >> main >> aux;
        }/*}}}*/
        cachefile& cachefile_save(cachefile& c) const {/*{{{*/
            return c << i0 << j0
                << index << offset << count
                << main << aux;
        }/*}}}*/
    };

    size_t maxmaxw;
    size_t dispatch_strips;
    vector<dispatcher_t> dispatchers;
    vector<combiner_t> combiners;
    /* }}} */

    /* {{{ front-end */
    ~matmul_zone_data();
    matmul_zone_data(void* xx, param_list pl, int optimized_direction);
    void build_cache(uint32_t * data);
    int reload_cache();
    void save_cache();
    void mul(void * xdst, void const * xsrc, int d);
    void report(double scale MAYBE_UNUSED);
    void auxv(int op, va_list ap);
    void aux(int op, ...);
    /* }}} */

    private:/*{{{*/
    void create_dispatchers_and_combiners(vector<zone> & q1);
    /*}}}*/
};
inline cachefile& operator>>(cachefile& c, matmul_zone_data::dispatcher_t& z) { return z.cachefile_load(c); }
inline cachefile& operator<<(cachefile& c, matmul_zone_data::dispatcher_t const & z) { return z.cachefile_save(c); }
inline cachefile& operator>>(cachefile& c, matmul_zone_data::combiner_t& z) { return z.cachefile_load(c); }
inline cachefile& operator<<(cachefile& c, matmul_zone_data::combiner_t const & z) { return z.cachefile_save(c); }
/*}}}*/
/**************************************************************************/
/*{{{ trampolines for C bindings */
void MATMUL_NAME(clear)(matmul_ptr mm0)
{
    delete (matmul_zone_data *) mm0;
}

matmul_ptr MATMUL_NAME(init)(void* xx, param_list pl, int optimized_direction)
{
    return (matmul_ptr) new matmul_zone_data(xx, pl, optimized_direction);
}

void MATMUL_NAME(build_cache)(matmul_ptr mm0, uint32_t * data)
{
    ((matmul_zone_data*)mm0)->build_cache(data);
}

int MATMUL_NAME(reload_cache)(matmul_ptr mm0)
{
    return ((matmul_zone_data*)mm0)->reload_cache();
}

void MATMUL_NAME(save_cache)(matmul_ptr mm0)
{
    ((matmul_zone_data*)mm0)->save_cache();
}

void MATMUL_NAME(mul)(matmul_ptr mm0, void * xdst, void const * xsrc, int d)
{
    ((matmul_zone_data*)mm0)->mul(xdst, xsrc, d);
}

void MATMUL_NAME(report)(matmul_ptr mm0 MAYBE_UNUSED, double scale MAYBE_UNUSED) {
    ((matmul_zone_data*)mm0)->report(scale);
}

void MATMUL_NAME(auxv)(matmul_ptr mm0 MAYBE_UNUSED, int op MAYBE_UNUSED, va_list ap MAYBE_UNUSED)
{
    ((matmul_zone_data*)mm0)->auxv(op, ap);
}

void MATMUL_NAME(aux)(matmul_ptr mm0, int op, ...)
{
    va_list ap;
    va_start(ap, op);
    ((matmul_zone_data*)mm0)->auxv(op, ap);
    va_end(ap);
}
/*}}}*/

/**************************************************************************/

matmul_zone_data::~matmul_zone_data() {/*{{{*/
    matmul_common_clear(public_);
}
/*}}}*/
matmul_zone_data::matmul_zone_data(void* xab, param_list pl, int optimized_direction) : xab((abdst_field) xab)/*{{{*/
{
    memset(&public_, 0, sizeof(public_));
    int suggest = optimized_direction ^ MM_DIR0_PREFERS_TRANSP_MULT;
    public_->store_transposed = suggest;
    if (pl) {
        param_list_parse_int(pl, "mm_store_transposed",
                &public_->store_transposed);
        if (public_->store_transposed != suggest) {
            fprintf(stderr, "Warning, mm_store_transposed"
                    " overrides suggested matrix storage ordering\n");
        }
    }
}
/*}}}*/

#define CMAX 4
#define ROWBATCH        4096
#define COLBATCH        65536

struct sort_jc {/*{{{*/
    inline bool operator()(pair<uint32_t, int32_t> const& a, pair<uint32_t,int32_t> const& b) const {
        return a.first < b.first;
    }
};
/*}}}*/
void merge_dispatchers(vector<matmul_zone_data::dispatcher_t>& all, size_t maxmaxw)/*{{{*/
{
    typedef matmul_zone_data::dispatcher_t dispatcher_t;
    vector<dispatcher_t> merged;
    for(size_t k0 = 0, k1; k0 < all.size() ; ) {
        /* recompute maxw, since we haven't kept track */
        size_t maxw = 0;
        unsigned int j0 = all[k0].j0;
        for(k1 = k0 ; k1 < all.size() && all[k1].j0 == all[k0].j0 ; k1++) {
            maxw = max(maxw, all[k1].size());
        }
        /* the range [k0..k1[ has dispatchers for the same set of
         * columns. */
        unsigned int group = maxmaxw / maxw;
        unsigned int nm = 0;
        size_t old_k0 = k0;
        for( ; k0 < k1 ; ) {
            dispatcher_t D;
            D.i0 = all[k0].i0;
            D.j0 = all[k0].j0;
            /* concatenate. */
            for(unsigned int k2 = 0 ; k2 < group && k0 < k1 ; k2++, k0++) {
                D.insert(D.end(), all[k0].begin(), all[k0].end());
                /* free this memory as early as we can */
                all[k0].clear();
                /* XXX at this point, combiner must get the info that it
                 * should read in the buffer from index D.size()
                 */
            }
            merged.push_back(D);
            nm++;
        }
        printf("j0=%u: maxw = %zu ; group %u together. Split %zd dispatchers [%zu..%zu[ into %u dispatchers\n", j0, maxw, group, k1 - old_k0, old_k0, k1, nm);
    }
    printf("We now have %zu dispatchers, combined from %zu original\n",
            merged.size(), all.size());
    all.swap(merged);
}
/*}}}*/
void matmul_zone_data::create_dispatchers_and_combiners(vector<zone> & q1)/*{{{*/
{
    /* q1 has vertical strips which are *not* the first strip */
    sort(q1.begin(), q1.end(), zone::colmajor_sorter());

    maxmaxw = 0;
    for(size_t k0 = 0, k1; k0 < q1.size() ; k0 = k1) {
        /* find the max weight of cells for this j0 */
        /* here, weight is the sum of absolute values. Oh, and we merge
         * qp and qm, too. */
        size_t maxw = 0;
        for(k1 = k0 ; k1 < q1.size() && q1[k1].j0 == q1[k0].j0 ; k1++) {
            /* qp and qm already have repeated coeffs. For qg, we'll
             * store coeffs just at the end of the array */
            maxw = max(maxw, q1[k1].size());
        }
        maxmaxw = max(maxmaxw, maxw);
    }

    /* now prepare the pre-reading of source coeffs into the buffers. We
     * need one array which directs data in the first pass, and another
     * one for the second pass.
     */

    for(size_t k = 0 ; k < q1.size() ; k++) {
        dispatcher_t D;
        combiner_t C;
        D.i0 = C.i0 = q1[k].i0;
        D.j0 = C.j0 = q1[k].j0;
        for(size_t i = 0 ; i < q1[k].qp.size() ; i++) {
            uint16_t col_id = q1[k].qp[i].second;
            uint16_t destrow_id = q1[k].qp[i].first;
            D.push_back(col_id);
            C.main.push_back(destrow_id);
        }
        for(size_t i = 0 ; i < q1[k].qm.size() ; i++) {
            uint16_t col_id = q1[k].qm[i].second;
            uint16_t destrow_id = q1[k].qm[i].first;
            D.push_back(col_id);
            /* specific offset for negative coefficients */
            C.main.push_back(destrow_id + ROWBATCH);
        }
        for(size_t i = 0 ; i < q1[k].qg.size() ; i++) {
            uint16_t col_id = q1[k].qg[i].second.first;
            uint16_t destrow_id = q1[k].qg[i].first;
            int32_t coeff = q1[k].qg[i].second.second;
            D.push_back(col_id);
            C.aux.push_back(make_pair(destrow_id, coeff));
        }
        /* XXX I think that D is sorted column-major, right ? */
        dispatchers.push_back(D);
        combiners.push_back(C);
    }
    q1.clear();

    /* We now merge dispatchers */
    merge_dispatchers(dispatchers, maxmaxw);

    /* and sort combiners so that they are organized row-major */
    /* notice that the overall weight of combiners may be a concern in
     * the end.
     */
    sort(combiners.begin(), combiners.end(), combiner_t::rowmajor_sorter());
}
/*}}}*/
void matmul_zone_data::build_cache(uint32_t * data)/*{{{*/
{
    matmul_zone_data * mm = this;

    ASSERT_ALWAYS(data);

    unsigned int nrows_t = mm->public_->dim[ mm->public_->store_transposed];

#if COLBATCH
    unsigned int ncols_t = mm->public_->dim[!mm->public_->store_transposed];
#endif

    uint32_t * ptr = data;

    /* count coefficients */
    mm->public_->ncoeffs = 0;

    uint64_t ccount[2*CMAX + 1] = {0,};
    double zavg = 0;

    vector<zone> q1;

    for(unsigned int i0 = 0 ; i0 < nrows_t ; i0 += ROWBATCH) {
        /* Create the data structures for the horizontal strip starting
         * at row i0, column 0.
         */

        /* Because this horizontal strip will be split in many blocks, we
         * need to have a batch of pointers for reading each row. */
        uint32_t * pp[ROWBATCH + 1];
        uint32_t * cc[ROWBATCH + 1];
        pp[0] = ptr;
        for(unsigned int k = 0 ; k < ROWBATCH ; k++) {
            cc[k] = pp[k] + 1;
            if (i0 + k < nrows_t) {
                uint32_t weight = *pp[k];
                pp[k+1] = pp[k] + 1 + 2*weight;
                mm->public_->ncoeffs += weight;
                /* This is very important. We must sort rows before
                 * processing. */
                pair<uint32_t, int32_t> * cb = (pair<uint32_t, int32_t> *) cc[k];
                pair<uint32_t, int32_t> * ce = cb + weight;
                sort(cb, ce, sort_jc());
            } else {
                pp[k+1] = pp[k];
            }
        }
        ptr = pp[ROWBATCH];
        for(unsigned int j0 = 0 ; j0 < ncols_t ; j0 += COLBATCH) {
            zone z(i0, j0);
            for(unsigned int k = 0 ; k < ROWBATCH ; k++) {
                for( ; cc[k] < pp[k+1] ; cc[k] += 2) {
                    uint32_t j = cc[k][0] - j0;
                    int32_t c = cc[k][1];
                    if (j >= COLBATCH) break;
                    if (c < 0 && c >= -CMAX) {
                        ccount[CMAX + c]++;
                        for( ; c++ ; ) {
                            z.qm.push_back(make_pair(k, j));
                        }
                    } else if (c > 0 && c <= CMAX) {
                        ccount[CMAX + c]++;
                        for( ; c-- ; ) {
                            z.qp.push_back(make_pair(k, j));
                        }
                    } else {
                        ccount[CMAX]++;
                        z.qg.push_back(make_pair(k, make_pair(j, c)));
                    }
                }
            }
            z.sort();
            // printf("Zone %zu: %zu+%zu+%zu\n", q.size(), z.qp.size(), z.qm.size(), z.qg.size());
            if (!z.empty()) {
                if (j0) {
                    q1.push_back(z);
                    q.push_back(z);
                } else {
                    q.push_back(z);
                }
            }
            zavg += z.size();
        }
    }
    ASSERT_ALWAYS(ptr - data == (ptrdiff_t) (nrows_t + 2 * mm->public_->ncoeffs));
    sort(q.begin(), q.end(), zone::rowmajor_sorter());

    /* each zone in q1 becomes a dispatcher and a combiner */
    /* The first pass sets bucket_id to zero */
    create_dispatchers_and_combiners(q1);

    free(data);
    ostringstream os;
    for(int i = -CMAX ; i <= CMAX ; i++) {
        os << " " << i << ":" << (double) ccount[CMAX + i] / nrows_t;
    }
    printf("Stats: [%" PRIu64 "] %s, %zu zones of average weight %.1f\n", mm->public_->ncoeffs, os.str().c_str(), q.size(), zavg / q.size());
}
/*}}}*/
/* cache load and save {{{ */
int matmul_zone_data::reload_cache()
{
    FILE * f = matmul_common_reload_cache_fopen(sizeof(abelt), public_, MM_MAGIC);
    if (f == NULL) { return 0; }
    cachefile c(f);
    c >> q;
    fclose(f);

    return 1;
}

void matmul_zone_data::save_cache()
{
    FILE * f = matmul_common_save_cache_fopen(sizeof(abelt), public_, MM_MAGIC);
    cachefile c(f);
    c << q;
    fclose(f);
}
/* }}} */

void zone::operator()(gfp::elt_ur * tdst, const gfp::elt * tsrc) const
{
    for(qpm_t::const_iterator u = qp.begin() ; u != qp.end() ; u++) {
        uint16_t i = u->first;
        uint16_t j = u->second;
        gfp::add(tdst[i], tsrc[j]);
    }
    for(qpm_t::const_iterator u = qm.begin() ; u != qm.end() ; u++) {
        uint16_t i = u->first;
        uint16_t j = u->second;
        gfp::sub(tdst[i], tsrc[j]);
    }
    for(qg_t::const_iterator u = qg.begin() ; u != qg.end() ; u++) {
        uint16_t i = u->first;
        uint16_t j = u->second.first;
        int32_t c = u->second.second;
        if (c>0) {
            gfp::addmul_ui(tdst[i], tsrc[j], c);
        } else {
            gfp::submul_ui(tdst[i], tsrc[j], -c);
        }
    }
}

void matmul_zone_data::mul(void * xdst, void const * xsrc, int d)
{
    matmul_zone_data * mm = this;
    ASM_COMMENT("multiplication code");
    abdst_field x = mm->xab;
    const gfp::elt * src = (const gfp::elt *) xsrc;
    gfp::elt * dst = (gfp::elt *) xdst;

    gfp::preinv preinverse;
    gfp::elt prime;

    {
        mpz_t p;
        mpz_init(p);
        abfield_characteristic(x, p);
        prime = p;
        mpz_clear(p);
    }

    gfp::compute_preinv(preinverse, prime);

    /* d == 1: matrix times vector product */
    /* d == 0: vector times matrix product */

    /* However the matrix may be stored either row-major
     * (store_transposed == 0) or column-major (store_transposed == 1)
     */

    gfp::elt::zero(dst, mm->public_->dim[!d]);

    /* TODO: missing in mpfq elt_ur_{add,sub}_elt
     */
    if (d == !mm->public_->store_transposed) {
        /* We need to find the position of all dispatchers */
        vector<vector<dispatcher_t>::const_iterator> start_points;
        vector<vector<dispatcher_t>::const_iterator> current_points;
        vector<abvec_ur> current_buffers;
        vector<abvec_ur> current_buffer_pointers;

        {       /* prepare dispatcher heads and buffers {{{ */
            unsigned int last_j0 = UINT_MAX;
            for(vector<dispatcher_t>::const_iterator v = dispatchers.begin() ; v != dispatchers.end() ; v++) {
                if (v->j0 != last_j0) {
                    start_points.push_back(v);
                    current_points.push_back(v);
                    abvec_ur buf;
                    abvec_ur_init(x, &buf, maxmaxw);
                    current_buffers.push_back(buf);
                    current_buffer_pointers.push_back(buf);
                    last_j0 = v->j0;
                }
            }
            start_points.push_back(dispatchers.end());
            // printf("Found %zu different dispatcher strips\n", current_points.size());
        } /* }}} */

        gfp::elt_ur * tdst = new gfp::elt_ur[ROWBATCH];

        ASM_COMMENT("critical loop");

        unsigned int last_i0 = UINT_MAX;
        /*
        unsigned int last_j0 = UINT_MAX;
        twn current;
        */

        /* This loops processes the blocks in their sorting order. Here
         * we assume it's row-wise.
         */
        for(size_t k = 0 ; k < q.size() ; k++) {
            zone const& z = q[k];

            /* loop through all dispatcher strips, and pre-fill buffers if
             * we happen to need it. */

#if 0
            for(size_t i = 0 ; i < current_points.size() ; i++) {
                /* If this dispatcher strip is done, surely we don't want
                 * to check. Otherwise, the check is whether the ``to be
                 * done'' value i0 is ours */
                if (current_points[i] == start_points[i + 1]) continue;
                if (current_points[i]->i0 > z.i0) continue;
                /* right, so we can fill the buffer. Let's go. */
                absrc_vec tsrc = abvec_subvec_const(x, src, current_points[i]->j0);
                /* oh, but we need to find the offsets of each bucket... */



            }
#endif
            if (z.i0 != last_i0) {
                if (last_i0 != UINT_MAX) {
                    for(unsigned int i = 0 ; i < ROWBATCH ; i++) {
                        gfp::reduce(dst[last_i0 + i], tdst[i], prime, preinverse);
                    }
                }
                /* This operation below sets stuff to zero, so it's here
                 * that we assume blocks are sorted row-wise.
                 */
                gfp::elt_ur::zero(tdst, ROWBATCH);
                last_i0 = z.i0;
            }

            /*
            if (z.j0 != last_j0) {
                if (last_j0 != UINT_MAX) {
                    current.tt += cputicks();
                    tmap_t::iterator it = tmap.find(last_j0);
                    if (it == tmap.end()) {
                        tmap.insert(make_pair(last_j0, current));
                    } else {
                        it->second.tt += current.tt;
                        it->second.w += current.w;
                        it->second.n++;
                    }
                }
                last_j0 = z.j0;
                current.tt = -cputicks();
                current.w = z.size();
                current.n = 1;
            }
            */
            const gfp::elt * tsrc = src + z.j0;

            z(tdst, tsrc);

        }
        {
            /* reduce last batch. It could possibly be incomplete */
            unsigned int active = MIN(ROWBATCH, mm->public_->dim[!d] - last_i0);
            for(unsigned int i = 0 ; i < active ; i++) {
                gfp::reduce(dst[last_i0 + i], tdst[i], prime, preinverse);
            }
        }

        /*
        if (last_j0 != UINT_MAX) {
            current.tt += cputicks();
            tmap_t::iterator it = tmap.find(last_j0);
            if (it == tmap.end()) {
                tmap.insert(make_pair(last_j0, current));
            } else {
                it->second.tt += current.tt;
                it->second.w += current.w;
                it->second.n++;
            }
        }
        */
        ASM_COMMENT("end of critical loop");
        delete[] tdst;
    } else {
#if 0
        abvec_ur tdst;
        abvec_ur_init(x, &tdst, mm->public_->dim[!d]);
        abelt_ur tmp;
        abelt_ur_init(x, &tmp);
        if (mm->public_->iteration[d] == 10) {
            fprintf(stderr, "Warning: Doing many iterations with transposed code\n");
        }
        abvec_set_zero(x, dst, mm->public_->dim[!d]);
        abvec_ur_set_zero(x, tdst, mm->public_->dim[!d]);
        ASM_COMMENT("critical loop (transposed mult)");
        uint32_t * q;
        q = &(qp[0]);
        for(unsigned int i = 0 ; i < mm->public_->dim[d] ; i++) {
            uint32_t len = *q++;
            abelt_ur_set_elt(x, tmp, abvec_coeff_ptr_const(x, src, i));
            for(unsigned int j = 0 ; len-- ; ) {
                j = *q++;
                ASSERT(j < mm->public_->dim[!d]);
                abdst_elt_ur yj = abvec_ur_coeff_ptr(x, tdst, j);
                abelt_ur_add(x, yj, yj, tmp);
            }
        }
        q = &(qm[0]);
        for(unsigned int i = 0 ; i < mm->public_->dim[d] ; i++) {
            uint32_t len = *q++;
            abelt_ur_set_elt(x, tmp, abvec_coeff_ptr_const(x, src, i));
            for(unsigned int j = 0 ; len-- ; ) {
                j = *q++;
                ASSERT(j < mm->public_->dim[!d]);
                abdst_elt_ur yj = abvec_ur_coeff_ptr(x, tdst, j);
                abelt_ur_sub(x, yj, yj, tmp);
            }
        }
        q = &(qq[0]);
        for(unsigned int i = 0 ; i < mm->public_->dim[d] ; i++) {
            uint32_t len = *q++ / 2;
            if (!len) continue;
            abelt_ur_set_elt(x, tmp, abvec_coeff_ptr_const(x, src, i));
            for(unsigned int j = 0 ; len-- ; ) {
                j = *q++;
                int32_t c = *(int32_t*)q++;
                ASSERT(j < mm->public_->dim[!d]);
                abdst_elt_ur yj = abvec_ur_coeff_ptr(x, tdst, j);
                abaddmul_si_ur(x, yj, tmp, c);
            }
        }
        for(unsigned int j = 0 ; j < mm->public_->dim[!d] ; j++) {
            abreduce(x, abvec_coeff_ptr(x, dst, j), abvec_ur_coeff_ptr(x, tdst, j));
        }
        ASM_COMMENT("end of critical loop (transposed mult)");
        abelt_ur_clear(x, &tmp);
        abvec_ur_clear(x, &tdst, mm->public_->dim[!d]);
#endif
    }
    ASM_COMMENT("end of multiplication code");

    mm->public_->iteration[d]++;
}

void matmul_zone_data::report(double scale MAYBE_UNUSED)
{
    static pthread_mutex_t lk = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&lk);
    unsigned int niter = public_->iteration[0] + public_->iteration[1];
    for(tmap_t::const_iterator it = tmap.begin() ; it != tmap.end() ; it++) {
        printf("j0=%u [%u zones]: avg %.1f cycles/c [%.1f coeffs avg] - %.1f Mcycles/iter\n",
                it->first, it->second.n / niter, it->second.tt / scale / it->second.w, (double) it->second.w / it->second.n, (double) it->second.tt / niter * 1.0e-6);
    }
    pthread_mutex_unlock(&lk);
}

void matmul_zone_data::auxv(int op MAYBE_UNUSED, va_list ap MAYBE_UNUSED)
{
}

void matmul_zone_data::aux(int op, ...)
{
    va_list ap;
    va_start(ap, op);
    auxv(op, ap);
    va_end(ap);
}

/* vim: set sw=4: */
