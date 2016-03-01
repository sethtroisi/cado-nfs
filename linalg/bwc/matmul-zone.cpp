#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
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
struct zone {
    typedef vector<pair<int, uint32_t>> qpm_t;
    typedef vector<pair<int, pair<uint32_t, int32_t>>> qg_t;
    unsigned int i0, j0;
    qpm_t qp, qm;
    qg_t qg;
    zone(unsigned int i0, unsigned int j0) : i0(i0), j0(j0) {}
    inline bool empty() const { return qp.empty() && qm.empty() && qg.empty(); }
    inline size_t size() const { return qp.size() + qm.size() + qg.size(); }
    void operator()(abdst_field x, abdst_vec_ur tdst, absrc_vec tsrc) const;

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
#if 1
    struct sorter {
        inline bool operator()(zone const& a, zone const& b) const {
            return a.i0 < b.i0 || (a.i0 == b.i0 && a.j0 < b.j0);
        }
    };
#else
    struct sorter {
        inline bool operator()(zone const& a, zone const& b) const {
            return a.j0 < b.j0 || (a.j0 == b.j0 && a.i0 < b.i0);
        }
    };
#endif
};

struct matmul_zone_data {
    /* repeat the fields from the public interface */
    struct matmul_public_s public_[1];
    /* now our private fields */
    abdst_field xab;

    vector<zone> q;

    struct twn {
        uint64_t tt;
        size_t w;
        unsigned int n;
        // twn() : tt(0), w(0), n(0) {}
    };

    typedef map<unsigned int, twn> tmap_t;
    tmap_t tmap;

    ~matmul_zone_data();
    matmul_zone_data(void* xx, param_list pl, int optimized_direction);
    void build_cache(uint32_t * data);
    int reload_cache();
    void save_cache();
    void mul(void * xdst, void const * xsrc, int d);
    void report(double scale MAYBE_UNUSED);
    void auxv(int op, va_list ap);
    void aux(int op, ...);
};

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

matmul_zone_data::~matmul_zone_data() {
    matmul_common_clear(public_);
}

matmul_zone_data::matmul_zone_data(void* xab, param_list pl, int optimized_direction) : xab((abdst_field) xab)
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

#define CMAX 4
#define ROWBATCH        128
#define COLBATCH        16384

struct sort_jc {
    inline bool operator()(pair<uint32_t, int32_t> const& a, pair<uint32_t,int32_t> const& b) const {
        return a.first < b.first;
    }
};

void matmul_zone_data::build_cache(uint32_t * data)
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

    for(unsigned int i0 = 0 ; i0 < nrows_t ; i0 += ROWBATCH) {
        uint32_t * pp[ROWBATCH + 1];
        uint32_t * cc[ROWBATCH + 1];
        pp[0] = ptr;
        for(unsigned int k = 0 ; k < ROWBATCH ; k++) {
            cc[k] = pp[k] + 1;
            if (i0 + k < nrows_t) {
                pp[k+1] = pp[k] + 1 + 2*(*pp[k]);
                mm->public_->ncoeffs += *pp[k];
                /* This is very important. We must sort rows before
                 * processing. */
                pair<uint32_t, int32_t> * cb = (pair<uint32_t, int32_t> *) cc[k];
                pair<uint32_t, int32_t> * ce = cb + *pp[k];
                sort(cb, ce, sort_jc());
            } else {
                pp[k+1] = pp[k];
            }
        }
        ptr = pp[ROWBATCH];
#if COLBATCH
        for(unsigned int j0 = 0 ; j0 < ncols_t ; j0 += COLBATCH) {
#else
        {
            unsigned int j0 = 0;
#endif
            zone z(i0, j0);
            for(unsigned int k = 0 ; k < ROWBATCH ; k++) {
                for( ; cc[k] < pp[k+1] ; cc[k] += 2) {
                    uint32_t j = cc[k][0] - j0;
                    int32_t c = cc[k][1];
#if COLBATCH
                    if (j >= COLBATCH) break;
#endif
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
            if (!z.empty()) q.push_back(z);
            zavg += z.size();
        }
    }
    sort(q.begin(), q.end(), zone::sorter());
    ASSERT_ALWAYS(ptr - data == (ptrdiff_t) (nrows_t + 2 * mm->public_->ncoeffs));
    free(data);
    ostringstream os;
    for(int i = -CMAX ; i <= CMAX ; i++) {
        os << " " << i << ":" << (double) ccount[CMAX + i] / nrows_t;
    }
    printf("Stats: [%" PRIu64 "] %s, %zu zones of average weight %.1f\n", mm->public_->ncoeffs, os.str().c_str(), q.size(), zavg / q.size());
}

int matmul_zone_data::reload_cache()
{
    FILE * f = matmul_common_reload_cache_fopen(sizeof(abelt), public_, MM_MAGIC);
    if (f == NULL) { return 0; }

    size_t qsize;

    MATMUL_COMMON_READ_ONE32(qsize, f);
    q.insert(q.end(), qsize, zone(0,0));
    for(size_t i = 0 ; i < qsize ; i++) {
        size_t qpsize, qmsize, qgsize;
        MATMUL_COMMON_READ_ONE32(q[i].i0, f);
        MATMUL_COMMON_READ_ONE32(q[i].j0, f);
        MATMUL_COMMON_READ_ONE32(qpsize, f);
        MATMUL_COMMON_READ_ONE32(qmsize, f);
        MATMUL_COMMON_READ_ONE32(qgsize, f);
        q[i].qp.insert(q[i].qp.end(), qpsize, zone::qpm_t::value_type());
        q[i].qm.insert(q[i].qm.end(), qmsize, zone::qpm_t::value_type());
        q[i].qg.insert(q[i].qg.end(), qgsize, zone::qg_t::value_type());
        MATMUL_COMMON_READ_MANY32(&(q[i].qp[0]), 2 * q[i].qp.size(), f);
        MATMUL_COMMON_READ_MANY32(&(q[i].qm[0]), 2 * q[i].qm.size(), f);
        MATMUL_COMMON_READ_MANY32(&(q[i].qg[0]), 3 * q[i].qg.size(), f);
    }
    fclose(f);

    return 1;
}

void matmul_zone_data::save_cache()
{
    FILE * f = matmul_common_save_cache_fopen(sizeof(abelt), public_, MM_MAGIC);

    MATMUL_COMMON_WRITE_ONE32(q.size(), f);
    for(size_t i = 0 ; i < q.size() ; i++) {
        MATMUL_COMMON_WRITE_ONE32(q[i].i0, f);
        MATMUL_COMMON_WRITE_ONE32(q[i].j0, f);
        MATMUL_COMMON_WRITE_ONE32(q[i].qp.size(), f);
        MATMUL_COMMON_WRITE_ONE32(q[i].qm.size(), f);
        MATMUL_COMMON_WRITE_ONE32(q[i].qg.size(), f);
        MATMUL_COMMON_WRITE_MANY32(&(q[i].qp[0]), 2 * q[i].qp.size(), f);
        MATMUL_COMMON_WRITE_MANY32(&(q[i].qm[0]), 2 * q[i].qm.size(), f);
        MATMUL_COMMON_WRITE_MANY32(&(q[i].qg[0]), 3 * q[i].qg.size(), f);
    }
    fclose(f);
}

void zone::operator()(abdst_field x, abdst_vec_ur tdst, absrc_vec tsrc) const
{
    abelt_ur tmp;
    abelt_ur_init(x, &tmp);
    for(qpm_t::const_iterator u = qp.begin() ; u != qp.end() ; u++) {
        unsigned int i = u->first;
        unsigned int j = u->second;
        abelt_ur_set_elt(x,
                tmp,
                abvec_coeff_ptr_const(x, tsrc, j));
        abelt_ur_add(x,
                abvec_ur_coeff_ptr(x, tdst, i),
                abvec_ur_coeff_ptr(x, tdst, i),
                tmp);
    }
    for(qpm_t::const_iterator u = qm.begin() ; u != qm.end() ; u++) {
        unsigned int i = u->first;
        unsigned int j = u->second;
        abelt_ur_set_elt(x,
                tmp,
                abvec_coeff_ptr_const(x, tsrc, j));
        abelt_ur_sub(x,
                abvec_ur_coeff_ptr(x, tdst, i),
                abvec_ur_coeff_ptr(x, tdst, i),
                tmp);
    }
    for(qg_t::const_iterator u = qg.begin() ; u != qg.end() ; u++) {
        unsigned int i = u->first;
        unsigned int j = u->second.first;
        int32_t c = u->second.second;
        abaddmul_si_ur(x,
                abvec_ur_coeff_ptr(x, tdst, i),
                abvec_coeff_ptr_const(x, tsrc, j),
                c);
    }
    abelt_ur_clear(x, &tmp);
}

void matmul_zone_data::mul(void * xdst, void const * xsrc, int d)
{
    matmul_zone_data * mm = this;
    ASM_COMMENT("multiplication code");
    abdst_field x = mm->xab;
    absrc_vec src = (absrc_vec) xsrc; // typical C const problem.
    abdst_vec dst = (abdst_vec) xdst;

    /* d == 1: matrix times vector product */
    /* d == 0: vector times matrix product */

    /* However the matrix may be stored either row-major
     * (store_transposed == 0) or column-major (store_transposed == 1)
     */

    /* TODO: missing in mpfq   elt_ur_{add,sub}_elt
     */
    if (d == !mm->public_->store_transposed) {
        abvec_ur tdst ;
        abvec_ur_init(x, &tdst, ROWBATCH);

        abvec_set_zero(x, dst, mm->public_->dim[!d]);
        ASM_COMMENT("critical loop");
        unsigned int last_i0 = UINT_MAX;
        unsigned int active = 0;

        unsigned int last_j0 = UINT_MAX;
        twn current;
        for(size_t k = 0 ; k < q.size() ; k++) {
            zone const& z = q[k];
            if (z.i0 != last_i0) {
                if (last_i0 != UINT_MAX) {
                    for(unsigned int i = 0 ; i < active ; i++) {
                        abreduce(x, 
                                abvec_coeff_ptr(x, dst, last_i0 + i),
                                abvec_ur_coeff_ptr(x, tdst, i));
                    }
                }
                last_i0 = z.i0;
                active = MIN(ROWBATCH, mm->public_->dim[!d] - z.i0);
                for(unsigned int i = 0 ; i < active ; i++) {
                    abelt_ur_set_elt(x,
                            abvec_ur_coeff_ptr(x, tdst, i),
                            abvec_coeff_ptr_const(x, dst, z.i0 + i));
                }
            }
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
            absrc_vec tsrc = abvec_subvec_const(x, src, z.j0);
            z(x, tdst, tsrc);
        }
        if (last_i0 != UINT_MAX) {
            for(unsigned int i = 0 ; i < active ; i++) {
                abreduce(x, 
                        abvec_coeff_ptr(x, dst, last_i0 + i),
                        abvec_ur_coeff_ptr(x, tdst, i));
            }
        }
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
        ASM_COMMENT("end of critical loop");
        abvec_ur_clear(x, &tdst, ROWBATCH);
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
    for(tmap_t::const_iterator it = tmap.begin() ; it != tmap.end() ; it++) {
        printf("j0=%u [%u zones]: avg %.1f cycles/c [%" PRIu64 " coeffs]\n",
                it->first, it->second.n, it->second.tt / scale / it->second.w, it->second.w);
    }
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
