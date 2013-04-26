#include "cado.h"
#include "abase.h"
#include <stdlib.h>
#include <unistd.h>     // getpid for debugging
#include "portability.h"
#include "macros.h"
#include "polymat.h"
#include "bigpolymat.h"

/* {{{  init/zero/clear interface for bigpolymat */
polymat_ptr bigpolymat_my_cell(bigpolymat_ptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    return bigpolymat_cell(p, irank, jrank);
}

void bigpolymat_init_model(bigpolymat_ptr model, MPI_Comm comm, unsigned int m, unsigned int n)
{
    int size;
    int rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    ASSERT_ALWAYS(size == (int) (m * n));
    int irank = rank / n;
    int jrank = rank % n;

    memset(model, 0, sizeof(bigpolymat));
    model->m1 = m;
    model->n1 = n;
    MPI_Comm_split(comm, irank, jrank, &(model->row));
    MPI_Comm_split(comm, jrank, irank, &(model->col));
}

/* This completes the initialization process. This is _not_ a collective
 * operation */
void bigpolymat_finish_init(bigpolymat_ptr p, unsigned int m, unsigned int n, int len)
{
    ASSERT_ALWAYS(bigpolymat_check_pre_init(p));
    p->m = m;
    p->n = n;
    p->m0 = m / p->m1;
    p->n0 = n / p->n1;
    p->size = len;
    polymat_ptr me = bigpolymat_my_cell(p);
    /* For polymat, init or finish_init are identical */
    polymat_init(me, p->m0, p->n0, len);
}
/* If m,n,len are zero, then this is exactly equivalent to duplicating
 * the mpi_model. In which case the initialization may be completed later
 * on with bigpolymat_finish_init
 */
void bigpolymat_init(bigpolymat_ptr p, bigpolymat_srcptr model, unsigned int m, unsigned int n, int len)
{
    ASSERT_ALWAYS(model);
    ASSERT_ALWAYS(m % model->m1 == 0);
    ASSERT_ALWAYS(n % model->n1 == 0);
    int irank;
    int jrank;
    MPI_Comm_rank(model->col, &irank);
    MPI_Comm_rank(model->row, &jrank);
    memset(p, 0, sizeof(bigpolymat));
    p->m1 = model->m1;
    p->n1 = model->n1;
    MPI_Comm_dup(model->row, &(p->row));
    MPI_Comm_dup(model->col, &(p->col));

    p->cells = malloc(p->m1*p->n1*sizeof(polymat));
    memset(p->cells, 0, p->m1*p->n1*sizeof(polymat));
    /* Either none or all must be non-zero */
    ASSERT_ALWAYS((!m||!n||!len) ^ (m&&n&&len));
    if (!len)
        return;

    /* Allocate only our own cell. If we happen to appear as a left
     * multiplicand someday, we'll need recipient cells for our whole
     * row. If we happen to appear as a right multiplicand, we'll need
     * recipient cells for our whole columns. But there is no reason for
     * both situations (let alone any) to happen. Thus this allocation of
     * recipient cells is dynamic */
    bigpolymat_finish_init(p, m, n, len);
}

/* This checks the validity of the m,n,size fields. If the structure
 * corresponds to a lazy allocation state, this function returns 1. If
 * the structure corresponds to something fully initializaed already,
 * this function returns 0. If the fields are filled in an inconsistent
 * manner, abort() is called.
 */
int bigpolymat_check_pre_init(bigpolymat_srcptr p)
{
    if (p->m && p->n && p->size)
        return 0;
    if (!p->m && !p->n && !p->size)
        return 1;
    abort();
    return 0;
}


void bigpolymat_clear_model(bigpolymat_ptr p)
{
    MPI_Comm_free(&(p->row));
    MPI_Comm_free(&(p->col));
    memset(p, 0, sizeof(bigpolymat));
}

void bigpolymat_clear(bigpolymat_ptr p)
{
    for(unsigned int t = 0 ; t < p->m1 * p->n1 ; t++) {
        if (p->cells[t]->x) {
            polymat_clear(p->cells[t]);
        }
    }
    free(p->cells);
    bigpolymat_clear_model(p);
}


/* We are a left multiplicand. This is a no-op if space for our row has
 * already been allocated */
void bigpolymat_provision_row(bigpolymat_ptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    polymat_srcptr me = bigpolymat_my_cell(p);
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank) continue;
        polymat_ptr them = bigpolymat_cell(p, irank, j);
        if (them->x) continue;
        polymat_init(them, p->m0, p->n0, me->alloc);
    }
}
/* We are a right multiplicand. This is a no-op if space for our col has
 * already been allocated */
void bigpolymat_provision_col(bigpolymat_ptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    polymat_srcptr me = bigpolymat_my_cell(p);
    for(unsigned int i = 0 ; i < p->m1 ; i++) {
        if (i == (unsigned int) irank) continue;
        polymat_ptr them = bigpolymat_cell(p, i, jrank);
        if (them->x) continue;
        polymat_init(them, p->m0, p->n0, me->alloc);
    }
}

/* Set size to be large enough to receive the given number of
 * coefficients. If space is provisioned for other cells in row or
 * column, allocation is triggered for them as well.
 */
void bigpolymat_set_size(bigpolymat_ptr p, size_t size)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    polymat_ptr me = bigpolymat_my_cell(p);
    p->size = size;
    me->size = size;
    ASSERT_ALWAYS(me->size <= me->alloc);
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank) continue;
        polymat_ptr them = bigpolymat_cell(p, irank, j);
        if (!them->x) continue;
        them->size = size;
        ASSERT_ALWAYS(size <= them->alloc);
    }
    for(unsigned int i = 0 ; i < p->m1 ; i++) {
        if (i == (unsigned int) irank) continue;
        polymat_ptr them = bigpolymat_cell(p, i, jrank);
        if (!them->x) continue;
        them->size = size;
        ASSERT_ALWAYS(size <= them->alloc);
    }
}

/* If our row or col cells have already been allocated, then reallocate
 * them as well (XXX is it clear or not ?) */
void bigpolymat_realloc(bigpolymat_ptr p, int newalloc)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    polymat_ptr me = bigpolymat_my_cell(p);
    polymat_realloc(me, newalloc);
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank) continue;
        polymat_ptr them = bigpolymat_cell(p, irank, j);
        if (!them->x) continue;
        polymat_realloc(them, newalloc);
    }
    for(unsigned int i = 0 ; i < p->m1 ; i++) {
        if (i == (unsigned int) irank) continue;
        polymat_ptr them = bigpolymat_cell(p, i, jrank);
        if (!them->x) continue;
        polymat_realloc(them, newalloc);
    }
}

/* This zeroes out _our_ cell */
void bigpolymat_zero(bigpolymat_ptr p)
{
    polymat_ptr me = bigpolymat_my_cell(p);
    polymat_zero(me);
}

/* okay, it's ugly */
void bigpolymat_swap(bigpolymat_ptr a, bigpolymat_ptr b)
{
    bigpolymat x;
    memcpy(x, a, sizeof(bigpolymat));
    memcpy(a, b, sizeof(bigpolymat));
    memcpy(b, x, sizeof(bigpolymat));
}

/* }}} */
/* interface for polymat_ur is just query-replace */

/* Add coefficient of degree ka in a, and kb in b, into coeff of degree
 * kc in c.
 *
 * This involves only our local mat, and can be done with no
 * communication.
 */
void bigpolymat_addmat(abdst_field ab,
        bigpolymat c, unsigned int kc,
        bigpolymat a, unsigned int ka,
        bigpolymat b, unsigned int kb)
{
    polymat_ptr la = bigpolymat_my_cell(a);
    polymat_ptr lb = bigpolymat_my_cell(b);
    polymat_ptr lc = bigpolymat_my_cell(c);
    polymat_addmat(ab, lc, kc, la, ka, lb, kb);
}

/* see add */
void bigpolymat_submat(abdst_field ab,
        bigpolymat c, unsigned int kc,
        bigpolymat a, unsigned int ka,
        bigpolymat b, unsigned int kb)
{
    polymat_ptr la = bigpolymat_my_cell(a);
    polymat_ptr lb = bigpolymat_my_cell(b);
    polymat_ptr lc = bigpolymat_my_cell(c);
    polymat_submat(ab, lc, kc, la, ka, lb, kb);
}

#if 0
void bigpolymat_reducemat(abdst_field ab,
        bigpolymat c, unsigned int kc,
        bigpolymat_ur a, unsigned int ka)
{
    polymat_ur_ptr la = bigpolymat_ur_my_cell(a);
    polymat_ptr lc = bigpolymat_my_cell(c);
    polymat_reducemat(ab, lc, kc, la, ka);
}
#endif

/* {{{ allgather operations */
void bigpolymat_allgather_row(bigpolymat a)
{
    int irank;
    int jrank;
    MPI_Comm_rank(a->col, &irank);
    MPI_Comm_rank(a->row, &jrank);
    bigpolymat_provision_row(a);
    for(unsigned int k = 0 ; k < a->n1 ; k++) {
        polymat_ptr data = bigpolymat_cell(a, irank, k);
        /* XXX: Should we ensure earlier that we agree on the size ? */
        unsigned long size = data->size;
        MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, k, a->row);
        data->size = size;
        ASSERT_ALWAYS(data->size <= data->alloc);
        MPI_Bcast(data->x, data->m * data->n * data->size * sizeof(abelt), MPI_BYTE, k, a->row);
    }
}
void bigpolymat_allgather_col(bigpolymat a)
{
    int irank;
    int jrank;
    MPI_Comm_rank(a->col, &irank);
    MPI_Comm_rank(a->row, &jrank);
    bigpolymat_provision_col(a);
    for(unsigned int k = 0 ; k < a->m1 ; k++) {
        polymat_ptr data = bigpolymat_cell(a, k, jrank);
        /* XXX: Should we ensure earlier that we agree on the size ? */
        unsigned long size = data->size;
        MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, k, a->col);
        data->size = size;
        ASSERT_ALWAYS(data->size <= data->alloc);
        MPI_Bcast(data->x, data->m * data->n * data->size * sizeof(abelt), MPI_BYTE, k, a->col);
    }
}
/* scatter from node 0 */
/* TODO: post asynchronous sends ? */
void bigpolymat_scatter_row(bigpolymat a)
{
    int irank;
    int jrank;
    MPI_Comm_rank(a->col, &irank);
    MPI_Comm_rank(a->row, &jrank);
    bigpolymat_provision_row(a);

    for(unsigned int k = 0 ; k < a->n1 ; k++) {
        /* We are currently dispatching the data we have which belongs to
         * row k, to job k in this row */
        if (k == 0) continue;
        /* Send first data size, then payload from node jrank==0 to node jrank == k */
        /* XXX: Should we ensure earlier that we agree on the size ? */
        if (jrank == 0) {
            polymat_ptr data = bigpolymat_cell(a, irank, k);
            unsigned long size = data->size;
            MPI_Send(&size, 1, MPI_UNSIGNED_LONG, k, 0, a->row);
            MPI_Send(data->x, data->m * data->n * data->size * sizeof(abelt), MPI_BYTE, k, 1, a->row);
        } else if (jrank == (int) k) {
            polymat_ptr data = bigpolymat_my_cell(a);
            unsigned long size = 0;
            MPI_Recv(&size, 1, MPI_UNSIGNED_LONG, 0, 0, a->row, MPI_STATUS_IGNORE);
            data->size = size;
            ASSERT_ALWAYS(data->size <= data->alloc);
            MPI_Recv(data->x, data->m * data->n * data->size * sizeof(abelt), MPI_BYTE, 0, 1, a->row, MPI_STATUS_IGNORE);
        }
    }
}
/* }}} */

#if 0   /* I doubt we'll need these... */
/* This one needs communication */
void bigpolymat_addmulmat_ur(abdst_field ab,
        bigpolymat_ur c, unsigned int kc,
        bigpolymat a, unsigned int ka,
        bigpolymat b, unsigned int kb)
{
    int irank;
    int jrank;
    MPI_Comm_rank(c->col, &irank);
    MPI_Comm_rank(c->row, &jrank);
    bigpolymat_allgather_row(a);
    bigpolymat_allgather_col(b);
    polymat_ur_ptr lc = bigpolymat_ur_my_cell(c);
    ASSERT_ALWAYS(a->n == b->m);
    for(unsigned int k = 0 ; k < a->n ; k++) {
        polymat_addmulmat_ur(ab, lc, kc,
                bigpolymat_cell(a, irank, k), ka
                bigpolymat_cell(b, k, jrank), kb);
    }
}

void bigpolymat_mulmat(abdst_field ab,
        bigpolymat c, unsigned int kc,
        bigpolymat a, unsigned int ka,
        bigpolymat b, unsigned int kb)
{
    bigpolymat_ur cx;
    /* FIXME: interface-wise, having to multiply here is ugly. Prefer n1
     * for the mpi-level sizes ? */
    bigpolymat_ur_init(cx, c, c->m * c->m0, c->n * c->n0, 1);
    bigpolymat_addmulmat_ur(ab, cx, 0, a, ka, b, kb);
    bigpolymat_reducemat(ab, c, kc, cx, 0);
    bigpolymat_ur_clear(cx);
}
#endif

void bigpolymat_mul(abdst_field ab, bigpolymat c, bigpolymat a, bigpolymat b)/*{{{*/
{
    ASSERT_ALWAYS(a->n == b->m);
    ASSERT_ALWAYS(a->n1 == b->m1);
    if (!c->m) {
        ASSERT_ALWAYS(!c->m && !c->n && !c->size);
        bigpolymat_finish_init(c, a->m, b->n, a->size + b->size - 1);
    }
    ASSERT_ALWAYS(c->m);
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    int irank;
    int jrank;
    MPI_Comm_rank(c->col, &irank);
    MPI_Comm_rank(c->row, &jrank);
    bigpolymat_allgather_row(a);
    bigpolymat_allgather_col(b);
    c->size = a->size + b->size - 1;
    // ASSERT_ALWAYS(c->alloc >= c->size);

    /* Create a pre-destination, while accumulation will be done on c */
    bigpolymat c0;
    bigpolymat_init(c0, c, c->m, c->n, c->size);
    polymat_ptr lc0 = bigpolymat_my_cell(c0);
    polymat_ptr lc = bigpolymat_my_cell(c);
    bigpolymat_zero(c);
    lc->size = lc0->size = c->size;
    ASSERT_ALWAYS(a->n == b->m);
    for(unsigned int k = 0 ; k < a->n1 ; k++) {
        polymat_mul(ab, lc0, 
                bigpolymat_cell(a, irank, k),
                bigpolymat_cell(b, k, jrank));
        for(size_t i = 0 ; i < lc0->size ; i++) {
            polymat_addmat(ab, lc, i, lc, i, lc0, i);
        }
    }
    bigpolymat_clear(c0);
}/*}}}*/

void bigpolymat_mp_raw(abdst_field ab,/*{{{*/
        bigpolymat b, unsigned int xb,
        bigpolymat a, unsigned int xa, unsigned int na,
        bigpolymat c, unsigned int xc, unsigned int nc,
        int transpose, int add)
{
    unsigned int nb = MAX(na, nc) - MIN(na, nc) + 1;
    /*
    ASSERT_ALWAYS(a->n == c->m);
    ASSERT_ALWAYS(a->n1 == c->m1);
    ASSERT_ALWAYS(c->m);
    ASSERT_ALWAYS(b->m == a->m);
    ASSERT_ALWAYS(b->n == c->n);
    // ASSERT_ALWAYS(b->alloc >= nb);
    */
    int irank;
    int jrank;
    MPI_Comm_rank(b->col, &irank);
    MPI_Comm_rank(b->row, &jrank);

    if (!transpose) {
        bigpolymat_allgather_row(a);
        bigpolymat_allgather_col(c);
    } else {
        bigpolymat_allgather_row(c);
        bigpolymat_allgather_col(a);
    }

    // bigpolymat b0;
    // bigpolymat_init(b0, b, b->m, b->n, b->size);
    // polymat_ptr lb0 = bigpolymat_my_cell(b0);
    polymat_ptr lb = bigpolymat_my_cell(b);
    lb->size = c->size;
    // lb0->size = c->size;
    if (!add)
        bigpolymat_zero(b);
    bigpolymat_set_size(b, nb);
    for(unsigned int k = 0 ; k < a->n1 ; k++) {
        if (!transpose) {
            polymat_mp_raw(ab, lb, xb,
                    bigpolymat_cell(a, irank, k), xa, na,
                    bigpolymat_cell(c, k, jrank), xc, nc,
                    transpose, 1);
        } else {
            polymat_mp_raw(ab, lb, xb,
                    bigpolymat_cell(a, k, jrank), xa, na,
                    bigpolymat_cell(c, irank, k), xc, nc,
                    transpose, 1);
        }
        /*
        for(size_t i = 0 ; i < lb0->size ; i++) {
            polymat_addmat(ab, lb, i, lb, i, lb0, i);
        }
        */
    }
    // bigpolymat_clear(b0);
}/*}}}*/

/* Collect everything into node 0 */
void bigpolymat_gather_mat(abdst_field ab, polymat dst, bigpolymat src)
{
    if (polymat_check_pre_init(dst)) {
        polymat_init(dst, src->m, src->n, src->size);
    }
    polymat_zero(dst);  /* See comment below about the gather operation */
    ASSERT_ALWAYS(dst->m == src->m);
    ASSERT_ALWAYS(dst->n == src->n);
    int irank;
    int jrank;
    MPI_Comm_rank(src->col, &irank);
    MPI_Comm_rank(src->row, &jrank);
    /* Node 0 must receive data from everyone. Whether we seek
     * collaboration from some of the peer nodes to delegate some of the
     * receiving does not matter much, since _we_ will receive the whole
     * thing in the end. However, this does have an impact on the
     * communication pattern.
     */
    bigpolymat_allgather_row(src);
    /* Do this on all nodes. This will aid the coding of the next gather
     * */
    unsigned int ibase = irank * src->m0;
    for(unsigned int k = 0 ; k < src->size ; k++) {
        for(unsigned int i = 0 ; i < src->m0 ; i++) {
            for(unsigned int j1 = 0 ; j1 < src->n1 ; j1++) {
                /* Collect local row i from node j1 in this row */
                polymat_ptr them = bigpolymat_cell(src, irank, j1);
                bwmat_copy_coeffs(ab,
                        polymat_part(dst, ibase + i, j1 * src->n0, k), 1,
                        polymat_part(them, i, 0, k), 1,
                        src->n0);
            }
        }
    }
    /* Now across rows, each row has a copy of the complete row block it
     * belongs to (in the proper place in dst). Do gather() on node 0 in
     * the column, now (we could as well do allgather).
     *
     * Unfortunately, data is strided, as we have all coefficients
     * corresponding to our row block, but obviously these are spaced
     * apart.
     *
     * There are several possibilities:
     *  - post large number of nonblocking transfers. Imagine submatrices
     *    of size maybe 2*2, with 64-byte data, that means very tiny
     *    transfers (256 bytes). It is hard to imagine that reasonable
     *    performance may be obtained like this.
     *  - allocate some transitory buffer to do larger transfers, still
     *    asynchronous ones. It is perhaps possible to achieve fairly
     *    good performance like this.
     *  - handle striding by just ensuring that data is zero where we
     *    expect no coefficients (see polymat_zero above). We may thus do
     *    a simple MPI_BXOR. It's a bit cheating, but it works. We
     *    do at least a fraction 1/r of useless transfers.
     * Note that it is not clear at this point whether gather_mat is a
     * bottleneck or not.
     */
    /* size-wise, which MPI_Reduce is rather big. It's not unexpected,
     * since in the end, we do _gather_ everything here.
     */
    /* MPI_IN_PLACE semantics for Reduce() and Gather() are stupid. */
    MPI_Reduce(
            irank ? polymat_part(dst, 0, 0, 0) : MPI_IN_PLACE,
            irank ? NULL : polymat_part(dst, 0, 0, 0),
            src->m * src->n * src->size * sizeof(abelt),
            MPI_BYTE,
            MPI_BXOR, 0, src->col);
    dst->size = src->size;
}

/* Exactly the converse of the previous function. */
void bigpolymat_scatter_mat(abdst_field ab, bigpolymat_ptr dst, polymat_ptr src)
{
    int irank;
    int jrank;
    MPI_Comm_rank(dst->col, &irank);
    MPI_Comm_rank(dst->row, &jrank);
    /* src at root must be initialized. src may be pre-init at other
     * nodes */
    int pre_init_status[2];     /* local, max-global */
    pre_init_status[0] = polymat_check_pre_init(src);
    pre_init_status[1] = pre_init_status[0];
    MPI_Allreduce(MPI_IN_PLACE, &(pre_init_status[1]), 1, MPI_INT, MPI_MAX, dst->col);
    MPI_Allreduce(MPI_IN_PLACE, &(pre_init_status[1]), 1, MPI_INT, MPI_MAX, dst->row);

    if (irank == 0 && jrank == 0) {
        ASSERT_ALWAYS(pre_init_status[0] == 0);
    }

    if (pre_init_status[1]) {   /* nodes other than root are uninitialized */
        /* So, initialize them... */
        /* XXX. This is quite unelegant */
        if (irank || jrank) {
            ASSERT_ALWAYS(pre_init_status[0] == 1);
        }
        /* share allocation size. Don't share data, as it's not very
         * relevant. */
        polymat shell;
        shell->m = src->m;
        shell->n = src->n;
        shell->size = src->size;
        shell->alloc = src->alloc;
        /* 2-step broadcast */
        MPI_Bcast(shell, sizeof(polymat), MPI_BYTE, 0, dst->col);
        MPI_Bcast(shell, sizeof(polymat), MPI_BYTE, 0, dst->row);
        if (irank || jrank) {
            polymat_init(src, shell->m, shell->n, shell->alloc);
            src->size = shell->size;
        }
    }

    /* dst must be either initialized, or in pre-init mode */
    if (bigpolymat_check_pre_init(dst)) {
        bigpolymat_finish_init(dst, src->m, src->n, src->size);
    }
    ASSERT_ALWAYS(dst->m == src->m);
    ASSERT_ALWAYS(dst->n == src->n);
    bigpolymat_set_size(dst, src->size);
    /* See similar comment in gather_mat. We can't do scatter() here.
     * We're better off with a simple bcast.
     */
    MPI_Bcast(
            polymat_part(src, 0, 0, 0),
            dst->m * dst->n * dst->size * sizeof(abelt),
            MPI_BYTE,
            0,
            dst->col);
    /* Now scatter across rows */
    bigpolymat_provision_row(dst);
    /* Copy our polymat parts to our local bigpolymat cells. Then we'll
     * do scatter() (admittedly, this is one extra copy which we could
     * seek to avoid) */
    unsigned int ibase = irank * dst->m0;
    for(unsigned int k = 0 ; k < src->size ; k++) {
        for(unsigned int i = 0 ; i < dst->m0 ; i++) {
            for(unsigned int j1 = 0 ; j1 < dst->n1 ; j1++) {
                polymat_ptr them = bigpolymat_cell(dst, irank, j1);
                bwmat_copy_coeffs(ab,
                        polymat_part(them, i, 0, k), 1,
                        polymat_part(src, ibase + i, j1 * dst->n0, k), 1,
                        dst->n0);
            }
        }
    }
    /* Make sure all peer cells do have the data size field correctly set
     * up */
    bigpolymat_set_size(dst, dst->size);
    bigpolymat_scatter_row(dst);
}
