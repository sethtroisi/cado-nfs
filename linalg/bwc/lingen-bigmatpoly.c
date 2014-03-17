#include "cado.h"
#include "abase.h"
#include <stdlib.h>
#include "portability.h"
#include "macros.h"
#include "lingen-matpoly.h"
#include "lingen-bigmatpoly.h"

/* {{{  init/zero/clear interface for bigmatpoly */
matpoly_ptr bigmatpoly_my_cell(bigmatpoly_ptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    return bigmatpoly_cell(p, irank, jrank);
}

void bigmatpoly_init_model(bigmatpoly_ptr model, MPI_Comm comm, unsigned int m, unsigned int n)
{
    int size;
    int rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    ASSERT_ALWAYS(size == (int) (m * n));
    int irank = rank / n;
    int jrank = rank % n;

    memset(model, 0, sizeof(bigmatpoly));
    model->m1 = m;
    model->n1 = n;
    MPI_Comm_dup(comm, &(model->comm));
    MPI_Comm_split(comm, irank, jrank, &(model->row));
    MPI_Comm_split(comm, jrank, irank, &(model->col));
}

/* This completes the initialization process. This is _not_ a collective
 * operation */
void bigmatpoly_finish_init(abdst_field ab, bigmatpoly_ptr p, unsigned int m, unsigned int n, int len)
{
    ASSERT_ALWAYS(bigmatpoly_check_pre_init(p));
    p->m = m;
    p->n = n;
    p->m0 = m / p->m1;
    p->n0 = n / p->n1;
    p->size = len;
    matpoly_ptr me = bigmatpoly_my_cell(p);
    /* For matpoly, init or finish_init are identical */
    matpoly_init(ab, me, p->m0, p->n0, len);
}
/* If m,n,len are zero, then this is exactly equivalent to duplicating
 * the mpi_model. In which case the initialization may be completed later
 * on with bigmatpoly_finish_init
 */
void bigmatpoly_init(abdst_field ab, bigmatpoly_ptr p, bigmatpoly_srcptr model, unsigned int m, unsigned int n, int len)
{
    ASSERT_ALWAYS(model);
    ASSERT_ALWAYS(m % model->m1 == 0);
    ASSERT_ALWAYS(n % model->n1 == 0);
    int irank;
    int jrank;
    MPI_Comm_rank(model->col, &irank);
    MPI_Comm_rank(model->row, &jrank);
    memset(p, 0, sizeof(bigmatpoly));
    p->m1 = model->m1;
    p->n1 = model->n1;
    MPI_Comm_dup(model->comm, &(p->comm));
    MPI_Comm_dup(model->row, &(p->row));
    MPI_Comm_dup(model->col, &(p->col));

    p->cells = malloc(p->m1*p->n1*sizeof(matpoly));
    memset(p->cells, 0, p->m1*p->n1*sizeof(matpoly));
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
    bigmatpoly_finish_init(ab, p, m, n, len);
}

/* This checks the validity of the m,n,size fields. If the structure
 * corresponds to a lazy allocation state, this function returns 1. If
 * the structure corresponds to something fully initializaed already,
 * this function returns 0. If the fields are filled in an inconsistent
 * manner, abort() is called.
 */
int bigmatpoly_check_pre_init(bigmatpoly_srcptr p)
{
    if (p->m && p->n && p->size)
        return 0;
    if (!p->m && !p->n && !p->size)
        return 1;
    abort();
    return 0;
}


void bigmatpoly_clear_model(bigmatpoly_ptr p)
{
    MPI_Comm_free(&(p->comm));
    MPI_Comm_free(&(p->row));
    MPI_Comm_free(&(p->col));
    memset(p, 0, sizeof(bigmatpoly));
}

void bigmatpoly_clear(abdst_field ab, bigmatpoly_ptr p)
{
    for(unsigned int t = 0 ; t < p->m1 * p->n1 ; t++) {
        if (p->cells[t]->x) {
            matpoly_clear(ab, p->cells[t]);
        }
    }
    free(p->cells);
    bigmatpoly_clear_model(p);
}

/* Return a bitmask indicating whether bigmatpoly_provision_{row,col} has
 * been called on this matrix before. bit 0 is for row, bit 1 is for col.
 * If the returned value is zero, then we really have only the local part
 * of this matrix for the moment.
 */
int bigmatpoly_provisioned(bigmatpoly_ptr p)
{
    if (bigmatpoly_check_pre_init(p)) return 0;
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    unsigned int np_row = 0;
    unsigned int np_col = 0;
    for(unsigned int j = 0 ; j < p->n1 ; j++)
        np_row += bigmatpoly_cell(p, irank, j)->x != NULL;
    for(unsigned int i = 0 ; i < p->m1 ; i++)
        np_col += bigmatpoly_cell(p, i, jrank)->x != NULL;
    ASSERT_ALWAYS(np_row == 1 || np_row == p->n1);
    ASSERT_ALWAYS(np_col == 1 || np_col == p->m1);
    return (np_row>1)+((np_col>1)<<1);
}


/* We are a left multiplicand. This is a no-op if space for our row has
 * already been allocated */
void bigmatpoly_provision_row(abdst_field ab, bigmatpoly_ptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    matpoly_srcptr me = bigmatpoly_my_cell(p);
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, irank, j);
        if (them->x) continue;
        matpoly_init(ab, them, p->m0, p->n0, me->alloc);
    }
}
/* Rarely useful. We do need it because we resort to a kludgy
 * implementation of scatter_mat, which calls for provisioning on all
 * rows.
 */

void bigmatpoly_unprovision_row(abdst_field ab, bigmatpoly_ptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, irank, j);
        if (!them->x) continue;
        matpoly_clear(ab, them);
    }
}

/* We are a right multiplicand. This is a no-op if space for our col has
 * already been allocated */
void bigmatpoly_provision_col(abdst_field ab, bigmatpoly_ptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    matpoly_srcptr me = bigmatpoly_my_cell(p);
    for(unsigned int i = 0 ; i < p->m1 ; i++) {
        if (i == (unsigned int) irank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, i, jrank);
        if (them->x) continue;
        matpoly_init(ab, them, p->m0, p->n0, me->alloc);
    }
}

/* Set size to be large enough to receive the given number of
 * coefficients. If space is provisioned for other cells in row or
 * column, allocation is triggered for them as well.
 */
void bigmatpoly_set_size(bigmatpoly_ptr p, size_t size)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    matpoly_ptr me = bigmatpoly_my_cell(p);
    ASSERT_ALWAYS(size <= me->alloc);
    p->size = size;
    me->size = size;
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, irank, j);
        if (!them->x) continue;
        them->size = size;
        ASSERT_ALWAYS(size <= them->alloc);
    }
    for(unsigned int i = 0 ; i < p->m1 ; i++) {
        if (i == (unsigned int) irank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, i, jrank);
        if (!them->x) continue;
        them->size = size;
        ASSERT_ALWAYS(size <= them->alloc);
    }
}

/* If our row or col cells have already been allocated, then reallocate
 * them as well (XXX is it clear or not ?) */
#if 0 /* This function has never been used or needed */
void bigmatpoly_realloc(bigmatpoly_ptr p, int newalloc)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->col, &irank);
    MPI_Comm_rank(p->row, &jrank);
    matpoly_ptr me = bigmatpoly_my_cell(p);
    matpoly_realloc(me, newalloc);
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, irank, j);
        if (!them->x) continue;
        matpoly_realloc(them, newalloc);
    }
    for(unsigned int i = 0 ; i < p->m1 ; i++) {
        if (i == (unsigned int) irank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, i, jrank);
        if (!them->x) continue;
        matpoly_realloc(them, newalloc);
    }
}
#endif

/* This zeroes out _our_ cell */
void bigmatpoly_zero(abdst_field ab, bigmatpoly_ptr p)
{
    matpoly_ptr me = bigmatpoly_my_cell(p);
    matpoly_zero(ab, me);
}

/* okay, it's ugly */
void bigmatpoly_swap(bigmatpoly_ptr a, bigmatpoly_ptr b)
{
    bigmatpoly x;
    memcpy(x, a, sizeof(bigmatpoly));
    memcpy(a, b, sizeof(bigmatpoly));
    memcpy(b, x, sizeof(bigmatpoly));
}

/* }}} */

void bigmatpoly_truncate_loc(abdst_field ab, bigmatpoly_ptr dst, bigmatpoly_ptr src, unsigned int size)/*{{{*/
{
    ASSERT_ALWAYS(bigmatpoly_provisioned(src) == 0);
    ASSERT_ALWAYS(bigmatpoly_provisioned(dst) == 0);
    if (bigmatpoly_check_pre_init(dst)) {
        bigmatpoly_finish_init(ab, dst, src->m, src->n, size);
    }
    matpoly_truncate(ab, bigmatpoly_my_cell(dst), bigmatpoly_my_cell(src), size);
    dst->size = bigmatpoly_my_cell(dst)->size;
}
/*}}}*/

void bigmatpoly_rshift(abdst_field ab, bigmatpoly_ptr dst, bigmatpoly_ptr src, unsigned int k)/*{{{*/
{
    if (bigmatpoly_check_pre_init(dst)) {
        bigmatpoly_finish_init(ab, dst, src->m, src->n, src->size - k);
    }
    ASSERT_ALWAYS(bigmatpoly_provisioned(dst) == 0);
    ASSERT_ALWAYS(bigmatpoly_provisioned(src) == 0);
    matpoly_rshift(ab, bigmatpoly_my_cell(dst), bigmatpoly_my_cell(src), k);
    dst->size = bigmatpoly_my_cell(dst)->size;
}
/*}}}*/


/* {{{ allgather operations */
void bigmatpoly_allgather_row(abdst_field ab, bigmatpoly a)
{
    int irank;
    int jrank;
    MPI_Comm_rank(a->col, &irank);
    MPI_Comm_rank(a->row, &jrank);
    bigmatpoly_provision_row(ab, a);
    for(unsigned int k = 0 ; k < a->n1 ; k++) {
        matpoly_ptr data = bigmatpoly_cell(a, irank, k);
        /* XXX: Should we ensure earlier that we agree on the size ? */
        unsigned long size = data->size;
        MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, k, a->row);
        data->size = size;
        ASSERT_ALWAYS(data->size <= data->alloc);
        MPI_Bcast(data->x, data->m * data->n * data->size, abmpi_datatype(ab), k, a->row);
    }
}
void bigmatpoly_allgather_col(abdst_field ab, bigmatpoly a)
{
    int irank;
    int jrank;
    MPI_Comm_rank(a->col, &irank);
    MPI_Comm_rank(a->row, &jrank);
    bigmatpoly_provision_col(ab, a);
    for(unsigned int k = 0 ; k < a->m1 ; k++) {
        matpoly_ptr data = bigmatpoly_cell(a, k, jrank);
        /* XXX: Should we ensure earlier that we agree on the size ? */
        unsigned long size = data->size;
        MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, k, a->col);
        data->size = size;
        ASSERT_ALWAYS(data->size <= data->alloc);
        MPI_Bcast(data->x, data->m * data->n * data->size, abmpi_datatype(ab), k, a->col);
    }
}
/* scatter from node 0. This is not a very interesting function, in fact
 * it's used only within bigmatpoly_scatter_mat. We assume that data for
 * all rows is currently present at node 0 in each row, and we we
 * dispatch this data to the rows which actually do need it. */
/* TODO: post asynchronous sends ? */
static void bigmatpoly_scatter_row(abdst_field ab, bigmatpoly a)
{
    int irank;
    int jrank;
    MPI_Comm_rank(a->col, &irank);
    MPI_Comm_rank(a->row, &jrank);
    bigmatpoly_provision_row(ab, a);

    for(unsigned int k = 0 ; k < a->n1 ; k++) {
        /* We are currently dispatching the data we have which belongs to
         * row k, to job k in this row */
        if (k == 0) continue;
        /* Send first data size, then payload from node jrank==0 to node jrank == k */
        /* XXX: Should we ensure earlier that we agree on the size ? */
        if (jrank == 0) {
            matpoly_ptr data = bigmatpoly_cell(a, irank, k);
            unsigned long size = data->size;
            MPI_Send(&size, 1, MPI_UNSIGNED_LONG, k, 0, a->row);
            MPI_Send(data->x, data->m * data->n * data->size, abmpi_datatype(ab), k, 1, a->row);
        } else if (jrank == (int) k) {
            matpoly_ptr data = bigmatpoly_my_cell(a);
            unsigned long size = 0;
            MPI_Recv(&size, 1, MPI_UNSIGNED_LONG, 0, 0, a->row, MPI_STATUS_IGNORE);
            data->size = size;
            ASSERT_ALWAYS(data->size <= data->alloc);
            MPI_Recv(data->x, data->m * data->n * data->size, abmpi_datatype(ab), 0, 1, a->row, MPI_STATUS_IGNORE);
        }
    }
}
/* }}} */

void bigmatpoly_mul(abdst_field ab, bigmatpoly c, bigmatpoly a, bigmatpoly b)/*{{{*/
{
    ASSERT_ALWAYS(a->n == b->m);
    ASSERT_ALWAYS(a->n1 == b->m1);
    if (bigmatpoly_check_pre_init(c)) {
        bigmatpoly_finish_init(ab, c, a->m, b->n, a->size + b->size - 1);
    }
    ASSERT_ALWAYS(c->m);
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    int irank;
    int jrank;
    MPI_Comm_rank(c->col, &irank);
    MPI_Comm_rank(c->row, &jrank);
    bigmatpoly_allgather_row(ab, a);
    bigmatpoly_allgather_col(ab, b);
    bigmatpoly_set_size(c, a->size + b->size - 1);
    // ASSERT_ALWAYS(c->alloc >= c->size);

    matpoly_ptr lc = bigmatpoly_my_cell(c);
    bigmatpoly_zero(ab, c);
    lc->size = c->size;
    ASSERT_ALWAYS(a->n == b->m);
    for(unsigned int k = 0 ; k < a->n1 ; k++) {
        matpoly_addmul(ab, lc, 
                bigmatpoly_cell(a, irank, k),
                bigmatpoly_cell(b, k, jrank));
    }
}/*}}}*/

void bigmatpoly_mp(abdst_field ab,/*{{{*/
        bigmatpoly b,
        bigmatpoly a,
        bigmatpoly c)
{
    if (bigmatpoly_check_pre_init(b)) {
        bigmatpoly_finish_init(ab, b, a->m, c->n, MAX(a->size, c->size) - MIN(a->size, c->size) + 1);
    }
    ASSERT_ALWAYS(a->n == c->m);
    ASSERT_ALWAYS(a->n1 == c->m1);
    ASSERT_ALWAYS(b->m);
    ASSERT_ALWAYS(b->m == a->m);
    ASSERT_ALWAYS(b->n == c->n);
    int irank;
    int jrank;
    MPI_Comm_rank(b->col, &irank);
    MPI_Comm_rank(b->row, &jrank);

    bigmatpoly_allgather_row(ab, a);
    bigmatpoly_allgather_col(ab, c);

    matpoly_ptr lb = bigmatpoly_my_cell(b);
    matpoly_zero(ab, lb);
    bigmatpoly_set_size(b, b->size);
    ASSERT_ALWAYS(lb->size == b->size);
    ASSERT_ALWAYS(lb->alloc >= lb->size);

    for(unsigned int k = 0 ; k < a->n1 ; k++) {
        matpoly_addmp(ab, lb,
                bigmatpoly_cell(a, irank, k),
                bigmatpoly_cell(c, k, jrank));
    }
}/*}}}*/

/* Collect everything into node 0 */
void bigmatpoly_gather_mat(abdst_field ab, matpoly dst, bigmatpoly src)
{
    if (matpoly_check_pre_init(dst)) {
        matpoly_init(ab, dst, src->m, src->n, src->size);
    }
    matpoly_zero(ab, dst);  /* See comment below about the gather operation */
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
    bigmatpoly_allgather_row(ab, src);
    dst->size = src->size;
    /* Do this on all nodes. This will aid the coding of the next gather */
    unsigned int ibase = irank * src->m0;
    for(unsigned int i = 0 ; i < src->m0 ; i++) {
        for(unsigned int j1 = 0 ; j1 < src->n1 ; j1++) {
            /* Collect local row i from node j1 in this row */
            matpoly_ptr them = bigmatpoly_cell(src, irank, j1);
            matpoly_extract_row_fragment(ab,
                    dst, ibase+i, j1 * src->n0,
                    them, i, 0, src->n0);
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
     *    expect no coefficients (see matpoly_zero above). We may thus do
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
            irank ? matpoly_part(ab, dst, 0, 0, 0) : MPI_IN_PLACE,
            irank ? NULL : matpoly_part(ab, dst, 0, 0, 0),
            src->m * src->n * src->size,
            abmpi_datatype(ab),
            abmpi_addition_op(ab), 0, src->col);
}

/* Exactly the converse of the previous function. */
void bigmatpoly_scatter_mat(abdst_field ab, bigmatpoly_ptr dst, matpoly_ptr src)
{
    int irank;
    int jrank;
    MPI_Comm_rank(dst->col, &irank);
    MPI_Comm_rank(dst->row, &jrank);
    /* src at root must be initialized. src may be pre-init at other
     * nodes */
    int pre_init_status[2];     /* local, max-global */
    pre_init_status[0] = matpoly_check_pre_init(src);
    pre_init_status[1] = pre_init_status[0];
    MPI_Allreduce(MPI_IN_PLACE, &(pre_init_status[1]), 1, MPI_INT, MPI_MAX, dst->col);
    MPI_Allreduce(MPI_IN_PLACE, &(pre_init_status[1]), 1, MPI_INT, MPI_MAX, dst->row);

    if (irank == 0 && jrank == 0) {
        ASSERT_ALWAYS(pre_init_status[0] == 0);
    }

    /* source argument src at nodes other than root are uninitialized */
    if (pre_init_status[1]) {
        /* So, initialize them... */
        /* XXX. This is quite unelegant */
        if (irank || jrank) {
            ASSERT_ALWAYS(pre_init_status[0] == 1);
        }
        /* share allocation size. Don't share data, as it's not very
         * relevant. */
        matpoly shell;
        shell->m = src->m;
        shell->n = src->n;
        shell->size = src->size;
        shell->alloc = src->alloc;
        /* 2-step broadcast */
        MPI_Bcast(shell, sizeof(matpoly), MPI_BYTE, 0, dst->col);
        MPI_Bcast(shell, sizeof(matpoly), MPI_BYTE, 0, dst->row);
        if (irank || jrank) {
            /* make sure we have exactly the same amount of allocated
             * memory everywhere */
            /* XXX In effect, this allocates the root amount
             * *everywhere*, which is completely insane !
             */
            matpoly_init(ab, src, shell->m, shell->n, shell->alloc);
            src->size = shell->size;
        }
    }

    /* dst must be either initialized, or in pre-init mode */
    if (bigmatpoly_check_pre_init(dst)) {
        bigmatpoly_finish_init(ab, dst, src->m, src->n, src->size);
    }
    ASSERT_ALWAYS(dst->m == src->m);
    ASSERT_ALWAYS(dst->n == src->n);
    bigmatpoly_set_size(dst, src->size);
    /* See similar comment in gather_mat. We can't do scatter() here.
     * We're better off with a simple bcast.
     * 
     * Copy the full range of allocated bytes, not only up to size.
     */
    MPI_Bcast(
            matpoly_part(ab, src, 0, 0, 0),
            src->m * src->n * src->alloc,
            abmpi_datatype(ab),
            0,
            dst->col);
    /* Now scatter across rows */
    bigmatpoly_provision_row(ab, dst);
    /* propagate the size info */
    bigmatpoly_set_size(dst, dst->size);
    /* Copy our matpoly parts to our local bigmatpoly cells. Then we'll
     * do scatter() (admittedly, this is one extra copy which we could
     * seek to avoid) */
    unsigned int ibase = irank * dst->m0;
    for(unsigned int i = 0 ; i < dst->m0 ; i++) {
        for(unsigned int j1 = 0 ; j1 < dst->n1 ; j1++) {
            matpoly_ptr them = bigmatpoly_cell(dst, irank, j1);
            matpoly_extract_row_fragment(ab,
                    them, i, 0,
                    src, ibase+i, j1*dst->n0,
                    dst->n0);
        }
    }
    bigmatpoly_scatter_row(ab, dst);
    bigmatpoly_unprovision_row(ab, dst);
}

/* Exactly the converse of the previous function. */
void bigmatpoly_scatter_mat_alt(abdst_field ab, bigmatpoly_ptr dst, matpoly_ptr src)
{
    int rank;
    int irank;
    int jrank;
    MPI_Comm_rank(dst->comm, &rank);
    MPI_Comm_rank(dst->col, &irank);
    MPI_Comm_rank(dst->row, &jrank);
    /* src at root must be initialized.
     * src at other nodes is untouched. */

    /* share allocation size. Don't share data, as it's not very
     * relevant. */
    matpoly shell;
    shell->m = src->m;
    shell->n = src->n;
    shell->size = src->size;
    shell->alloc = src->alloc;

    MPI_Bcast(shell, sizeof(matpoly), MPI_BYTE, 0, dst->comm);

    /* dst must be in pre-init mode */
    ASSERT_ALWAYS(bigmatpoly_check_pre_init(dst));

    /* Important: allocate the same area */
    bigmatpoly_finish_init(ab, dst, shell->m, shell->n, shell->alloc);

    bigmatpoly_set_size(dst, shell->size);

    /* sanity check, because the code below assumes this. */
    ASSERT_ALWAYS(irank * (int) dst->n1 + jrank == rank);

    if (!rank) {
        MPI_Request * reqs = malloc(dst->m1 * dst->n1 * dst->m0 * sizeof(MPI_Request));
        MPI_Request * req = reqs;
        /* the master sends data to everyone ! */
        for(unsigned int i1 = 0 ; i1 < dst->m1 ; i1++) {
            for(unsigned int i0 = 0 ; i0 < dst->m0 ; i0++) {
                for(unsigned int j1 = 0 ; j1 < dst->n1 ; j1++) {
                    unsigned int ii = i1 * dst->m0 + i0;
                    unsigned int jj = j1 * dst->n0;
                    unsigned int count = dst->n0 * shell->alloc;
                    unsigned int peer = i1 * dst->n1 + j1;
                    unsigned int tag = ii * dst->n + jj;
                    abdst_vec from = matpoly_part(ab, src, ii, jj, 0);

                    if (peer == 0) {
                        /* talk to ourself */
                        matpoly_ptr me = bigmatpoly_my_cell(dst);
                        abdst_vec to = matpoly_part(ab, me, i0, 0, 0);
                        abvec_set(ab, to, from, count);
                    } else {
                        MPI_Isend(from, count, abmpi_datatype(ab),
                                peer, tag, dst->comm, req);
                    }

                    req++;
                }
            }
        }

        req = reqs;
        for(unsigned int i1 = 0 ; i1 < dst->m1 ; i1++) {
            for(unsigned int i0 = 0 ; i0 < dst->m0 ; i0++) {
                for(unsigned int j1 = 0 ; j1 < dst->n1 ; j1++) {
                    unsigned int peer = i1 * dst->n1 + j1;
                    if (peer)
                        MPI_Wait(req, MPI_STATUS_IGNORE);
                    req++;
                }
            }
        }
        free(reqs);
    } else {
        MPI_Request * reqs = malloc(dst->m0 * sizeof(MPI_Request));
        MPI_Request * req = reqs;
        /* receive. Each job will receive exactly dst->m0 transfers */
        matpoly_ptr me = bigmatpoly_my_cell(dst);
        for(unsigned int i0 = 0 ; i0 < dst->m0 ; i0++) {
            unsigned int i1 = irank;
            unsigned int j1 = jrank;
            unsigned int ii = i1 * dst->m0 + i0;
            unsigned int jj = j1 * dst->n0;
            unsigned int count = dst->n0 * shell->alloc;
            unsigned int tag = ii * dst->n + jj;
            abdst_vec to = matpoly_part(ab, me, i0, 0, 0);
            MPI_Irecv(to, count, abmpi_datatype(ab),
                    0, tag, dst->comm, req);
            req++;
        }
        req = reqs;
        for(unsigned int i0 = 0 ; i0 < dst->m0 ; i0++) {
            MPI_Wait(req, MPI_STATUS_IGNORE);
            req++;
        }
        free(reqs);
    }
}
