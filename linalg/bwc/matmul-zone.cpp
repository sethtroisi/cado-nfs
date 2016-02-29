#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>
#include <vector>
#include <utility>
#include <sstream>

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

struct matmul_zone_data {
    /* repeat the fields from the public interface */
    struct matmul_public_s public_[1];
    /* now our private fields */
    abdst_field xab;
    vector<uint32_t> qp, qm, qq;

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

void matmul_zone_data::build_cache(uint32_t * data)
{
    matmul_zone_data * mm = this;

    ASSERT_ALWAYS(data);

    unsigned int nrows_t = mm->public_->dim[ mm->public_->store_transposed];
    
    uint32_t * ptr = data;
    unsigned int i = 0;

    /* count coefficients */
    mm->public_->ncoeffs = 0;

#define CMAX 4

    uint64_t ccount[2*CMAX + 1] = {0,};

    for( ; i < nrows_t ; i++) {
        unsigned int weight = *ptr++;
        vector<uint32_t> lqp, lqm, lqq;
        for(unsigned int i = 0 ; i < weight ; i++, ptr += 2) {
            uint32_t j = ptr[0];
            int32_t c = ptr[1];
            if (c < 0 && c >= -CMAX) {
                ccount[CMAX + c]++;
                for( ; c++ ; ) lqm.push_back(j);
            } else if (c > 0 && c <= CMAX) {
                ccount[CMAX + c]++;
                for( ; c-- ; ) lqp.push_back(j);
            } else {
                ccount[CMAX]++;
                lqq.push_back(j);
                lqq.push_back((uint32_t) c);
            }
        } 
        qp.push_back(lqp.size()); qp.insert(qp.end(), lqp.begin(), lqp.end());
        qm.push_back(lqm.size()); qm.insert(qm.end(), lqm.begin(), lqm.end());
        qq.push_back(lqq.size()); qq.insert(qq.end(), lqq.begin(), lqq.end());

        mm->public_->ncoeffs += weight;
    }
    ASSERT_ALWAYS(ptr - data == (ptrdiff_t) (nrows_t + 2 * mm->public_->ncoeffs));
    free(data);
    ostringstream os;
    for(int i = -CMAX ; i <= CMAX ; i++) {
        os << " " << i << ":" << (double) ccount[CMAX + i] / nrows_t;
    }
    printf("Stats: %s\n", os.str().c_str());
}

int matmul_zone_data::reload_cache()
{
    FILE * f = matmul_common_reload_cache_fopen(sizeof(abelt), public_, MM_MAGIC);
    if (f == NULL) { return 0; }

    uint32_t nqp;
    uint32_t nqm;
    uint32_t nqq;
    MATMUL_COMMON_READ_ONE32(nqp, f);
    MATMUL_COMMON_READ_ONE32(nqm, f);
    MATMUL_COMMON_READ_ONE32(nqq, f);

    qp.insert(qp.end(), nqp, 0);
    MATMUL_COMMON_READ_MANY32(&(qp[0]), nqp, f);
    qm.insert(qm.end(), nqm, 0);
    MATMUL_COMMON_READ_MANY32(&(qm[0]), nqm, f);
    qq.insert(qq.end(), nqq, 0);
    MATMUL_COMMON_READ_MANY32(&(qq[0]), nqq, f);

    fclose(f);

    return 1;
}

void matmul_zone_data::save_cache()
{
    FILE * f = matmul_common_save_cache_fopen(sizeof(abelt), public_, MM_MAGIC);

    MATMUL_COMMON_WRITE_ONE32(qp.size(), f);
    MATMUL_COMMON_WRITE_ONE32(qm.size(), f);
    MATMUL_COMMON_WRITE_ONE32(qq.size(), f);

    MATMUL_COMMON_WRITE_MANY32(&(qp[0]), qp.size(), f);
    MATMUL_COMMON_WRITE_MANY32(&(qm[0]), qm.size(), f);
    MATMUL_COMMON_WRITE_MANY32(&(qq[0]), qq.size(), f);

    fclose(f);
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
        abelt_ur rowsum;
        abelt_ur_init(x, &rowsum);
        abelt_ur tmp;
        abelt_ur_init(x, &tmp);

        abvec_set_zero(x, dst, mm->public_->dim[!d]);
        ASM_COMMENT("critical loop");
        uint32_t * zp, * zm;
        zp = &(qp[0]);
        zm = &(qm[0]);
        uint32_t * q = &(qq[0]);					\
        for(unsigned int i = 0 ; i < mm->public_->dim[!d] ; i++) {
            uint32_t plen = *zp++;
            uint32_t mlen = *zm++;
            uint32_t glen = *q++ / 2;
            abelt_ur_set_elt(x, rowsum, abvec_coeff_ptr(x, dst, i));
            for(unsigned int j = 0 ; plen-- ; ) {
                j = *zp++;
                ASSERT(j < mm->public_->dim[d]);
                abelt_ur_set_elt(x, tmp, abvec_coeff_ptr_const(x, src, j));
                abelt_ur_add(x, rowsum, rowsum, tmp);
            }
            for(unsigned int j = 0 ; mlen-- ; ) {
                j = *zm++;
                ASSERT(j < mm->public_->dim[d]);
                abelt_ur_set_elt(x, tmp, abvec_coeff_ptr_const(x, src, j));
                abelt_ur_sub(x, rowsum, rowsum, tmp);
            }
            for(unsigned int j = 0 ; glen-- ; ) {
                j = *q++;
                int32_t c = *(int32_t*)q++;
                ASSERT(j < mm->public_->dim[d]);
                abaddmul_si_ur(x, rowsum, abvec_coeff_ptr_const(x, src, j), c);
            }
            abreduce(x, abvec_coeff_ptr(x, dst, i), rowsum);
        }
        ASM_COMMENT("end of critical loop");
        abelt_ur_clear(x, &rowsum);
        abelt_ur_clear(x, &tmp);
    } else {
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
    }
    ASM_COMMENT("end of multiplication code");

    mm->public_->iteration[d]++;
}

void matmul_zone_data::report(double scale MAYBE_UNUSED)
{
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
