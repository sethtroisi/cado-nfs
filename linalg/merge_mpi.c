/*
  MPI section
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "utils/utils.h"
#include "files.h"
#include "gzip.h"
#include "sparse.h"
#include "merge_mono.h"
#include "mpi.h"
#include "merge_mpi.h"

#define DEBUG 0

#define MPI_ERROR_TAG             -1
#define MPI_DIE_TAG                1
#define MPI_J_TAG                  2
#define MPI_ADD_TAG                3
#define MPI_SEND_W_TAG             4
#define MPI_SEND_ROW_TAG           5
#define MPI_SEND_ROW_BACK_TAG      6
#define MPI_MIN_M_TAG              7
#define MPI_DO_M_TAG               8
#define MPI_M_DONE_TAG             9
#define MPI_INACTIVATE_ROWS_TAG   10
#define MPI_ASK_MST_TAG           11
#define MPI_SEND_MST              12
#define MPI_SEND_HIS_TAG          13

#define MPI_BUF_SIZE 20000

unsigned int mpi_index = 0;

int *mpi_tab_j;

void
mpi_err(char *str)
{
    FILE *out = stderr;
    int mpi_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    fprintf(out, "MPI#%d[%d]# %s", mpi_rank, mpi_index, str);
    fflush(out);
}

void
mpi_err1(char *format, int i)
{
    FILE *out = stderr;
    int mpi_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    fprintf(out, "MPI#%d[%d]# ", mpi_rank, mpi_index);
    fprintf(out, format, i);
    fflush(out);
}

void
mpi_err2(char *format, int i1, int i2)
{
    FILE *out = stderr;
    int mpi_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    fprintf(out, "MPI#%d[%d]# ", mpi_rank, mpi_index);
    fprintf(out, format, i1, i2);
    fflush(out);
}

void
mpi_err3(char *format, int i1, int i2, int i3)
{
    FILE *out = stderr;
    int mpi_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    fprintf(out, "MPI#%d[%d]# ", mpi_rank, mpi_index);
    fprintf(out, format, i1, i2, i3);
    fflush(out);
}

void
mpi_err_tab(char *str, unsigned int *buf, int ibuf)
{
    FILE *out = stderr;
    int mpi_rank, i;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    fprintf(out, "MPI#%d[%d]# %s", mpi_rank, mpi_index, str);
    for(i = 0; i < ibuf; i++)
	fprintf(out, " %u", buf[i]);
    fprintf(out, "\n");
    fflush(out);
}

void
fprint_report_aux(FILE *out, report_t *rep)
{
    int i, k;

#if DEBUG >= 1
    mpi_err1("Report[0..%d]\n", rep->mark);
#endif
    for(i = 0; i <= rep->mark; i++){
	fprintf(out, "%d", mpi_index);
	for(k = 1; k <= rep->history[i][0]; k++)
	    fprintf(out, " %d", rep->history[i][k]);
	fprintf(out, "\n");
    }
    fflush(out);
}

void
fprint_report(report_t *rep)
{
    fprint_report_aux(rep->outfile, rep);
}

int
mpi_get_proc_for_j(int j)
{
    int i;

    for(i = 0; j >= mpi_tab_j[i]; i++);
    // j < mpi_tab_j[i], we hope
    return i;
}

void
mpi_kill_slaves()
{
    unsigned int imsg[1];
    int nprocs, rank;
    
    printf("MPI# Stopping everybody\n");
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    imsg[0] = 0;
    for(rank = 1; rank < nprocs; ++rank)
	MPI_Send(imsg, 1, MPI_UNSIGNED, rank, MPI_DIE_TAG, MPI_COMM_WORLD);
    printf("MPI# I sent everybody the die signal\n");
}

void
mpi_check_rows(sparse_mat_t *mat, INT i, INT *tab, int ntab)
{
    int k1, k2;
#if DEBUG >= 1
    int k;
    fprintf(stderr, "Original:         ");
    for(k = 1; k <= lengthRow(mat, i); k++)
	fprintf(stderr, " %d", cell(mat, i, k));
    fprintf(stderr, "\n");
    fprintf(stderr, "Reconstructed[%d]:", ntab);
    for(k = 0; k < ntab; k++)
	fprintf(stderr, " %d", tab[k]);
    fprintf(stderr, "\n");
    for(k = 0; k < ntab; k++){
	// check that tab[k] is indeed a position in mat->rows[i]
	int ok = 0, r;
	for(r = 1; r <= lengthRow(mat, i); r++)
	    if(cell(mat, i, r) == (int)tab[k]){
		ok = 1;
		break;
	    }
	if(!ok){
	    fprintf(stderr, "PB? for i=%d jj=%u\n", i, tab[k]);
	    exit(0);
	}
    }
#endif
    k1 = 0;
    k2 = 1;
    while(k1 < ntab){
	if(cell(mat, i, k2) == tab[k1]){
	    k2++; k1++;
	}
	else if(cell(mat, i, k2) == -1)
	    k2++;
	else{
	    fprintf(stderr, "PB? for i=%d/k1=%d/k2=%d\n", i, k1, k2);
	    exit(0);
	}
    }
    for( ; k2 <= lengthRow(mat, i); k2++)
	if(cell(mat, i, k2) != -1){
	    fprintf(stderr, "Non -1 cell!!!!\n");
	    exit(0);
	}
}

#define MAX_ROW_LENGTH 1000 // who cares, really?

#if 0 // useful???
void
mpi_load_rows_for_j(sparse_mat_t *mat, int m, INT j)
{
    MPI_Status status;
    unsigned int buf[MPI_BUF_SIZE];
    int mpi_size, ibuf, k, nrecv, cnt, ind;
    INT i, **newrows;
    char *done;

    fprintf(stderr, "Loading m=%d rows for j=%d\n", m, j);
    buf[0] = (unsigned)m;
    buf[1] = (unsigned)j;
    ibuf = 2;
    for(k = 1; k <= mat->R[j][0]; k++)
        if(mat->R[j][k] != -1)
	    buf[ibuf++] = (unsigned)mat->R[j][k];
    // at this point, we should have ibuf-2 == m...!
    if((ibuf-2) != m){
	fprintf(stderr, "#!# ibuf-2 != m in mpi_load_rows_for_j\n");
	exit(0);
    }
#if DEBUG >= 1
    fprintf(stderr, "Ready to send R[%d]=%d // #i=%d\n", j, 
	    mat->R[j][0], ibuf-2);
    fflush(stderr);
#endif
    // asking for the rows referenced by mat->R[j]
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    done = (char *)malloc(mpi_size * sizeof(char));
    for(k = 1; k < mpi_size; k++){
	MPI_Send(buf, ibuf, MPI_UNSIGNED, k, MPI_SEND_ROW_TAG, MPI_COMM_WORLD);
	done[k] = 0;
    }
    // feed newrows
    newrows = (INT **)malloc(m * sizeof(INT *));
    for(k = 0; k < m; k++){
	newrows[k] = (INT *)malloc(MAX_ROW_LENGTH * sizeof(INT));
	newrows[k][0] = (INT)buf[k+2];
	newrows[k][1] = 1; // last index fed
    }
    nrecv = 0;
    // waiting for the rows to come back from the slaves
    while(1){
#if DEBUG >= 1
	fprintf(stderr, "[m_l_r] Master waiting... nrecv=%d\n", nrecv);
	fflush(stderr);
#endif
	MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#if DEBUG >= 1
	fprintf(stderr, "Something received: %d\n", status.MPI_TAG);
	fflush(stderr);
#endif
	if(status.MPI_TAG == MPI_SEND_ROW_BACK_TAG){
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
#if DEBUG >= 1
	    fprintf(stderr, "[m_l_r] source=%u:", status.MPI_SOURCE);
	    for(k = 0; k < cnt; k++)
		fprintf(stderr, " %u", buf[k]);
	    fprintf(stderr, "\n");
#endif
	    // buf = [m, j, i, j_1, ..., j_r] where the indices j_s
	    // are positions in mat->rows[i]
	    i = (INT)buf[2];
	    if(isRowNull(mat, i)){
		fprintf(stderr, "Can it really happen??? i=%d j=%u cnt=%d\n", 
			i, buf[1], cnt);
		fflush(stderr);
	    }
	    // feed the corresponding row
	    for(ind = 0; ind < m; ind++)
		if(newrows[ind][0] == i)
		    break;
	    for(k = 3; k < cnt; k++){
		newrows[ind][1]++;
		if(newrows[ind][1] >= MAX_ROW_LENGTH){
		    fprintf(stderr, "#!# newrows.length exceeded, sorry!\n");
		    exit(0);
		}
		newrows[ind][newrows[ind][1]] = (INT)buf[k];
	    }
	    done[status.MPI_SOURCE] += 1;
	    if(done[status.MPI_SOURCE] == m){
		nrecv++;
		if(nrecv == (mpi_size-1))
		    break;
	    }
	}
    }
    for(k = 0; k < m; k++){
	qsort(newrows[k]+2, newrows[k][1]-1, sizeof(INT), cmp);
	mpi_check_rows(mat, newrows[k][0], newrows[k]+2, newrows[k][1]-1);
    }
    free(done);
    for(k = 0; k < m; k++)
	free(newrows[k]);
    free(newrows);
}
#endif

void
mpi_send_inactive_rows(int i)
{
    unsigned int buf[MPI_BUF_SIZE];
    int mpi_size, mpi_rank, k;

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    buf[0] = mpi_index;
    buf[1] = i;
    // TODO_MPI: include the master?
    for(k = 1; k < mpi_size; k++)
	if(k != mpi_rank)
	    MPI_Send(buf, 2, MPI_UNSIGNED, k, 
		     MPI_INACTIVATE_ROWS_TAG, MPI_COMM_WORLD);
}

int
mpi_get_number_of_active_colums(sparse_mat_t *mat, INT jmin, INT jmax)
{
    INT j;
    int nb = 0;

    for(j = jmin; j < jmax; j++)
	nb += (mat->wt[j] == 0 ? 0 : 1);
    return nb;
}

void
mpi_inactivate_rows(report_t *rep, sparse_mat_t *mat, unsigned int *tab, int ntab)
{
    int k, i;

    for(k = 1; k < ntab; k++){
	i = (int)tab[k];
#if DEBUG >= 1
	mpi_err1("Row[%d] has to be inactivated\n", i);
#endif
	if(isRowNull(mat, i))
	    mpi_err1("Row[%d] already nulled!!!\n", i);
	else
	    removeRowDefinitely(rep, mat, i);
    }
}

// Sending rows to add; get new weight in return; rows are sent to all
// procs except the current one and the master. At the end, the new weights
// are sent to the master for update.
// the submaster *has* to collect the weights itself, since otherwise the
// master would be receiving a lot of suborders in a probably random order.
void
mpi_add_rows(sparse_mat_t *mat, int m, INT j, INT *ind)
{
    MPI_Status status;
    unsigned int buf[MPI_BUF_SIZE];
    int new_weight[MERGE_LEVEL_MAX];
    int mpi_rank, mpi_size, ibuf, cnt, k, nrecv, finished;
    char *done;
    
#if DEBUG >= 1
    fprintf(stderr, "Column %d is for processor %d (m=%d)\n",
	    j, mpi_get_proc_for_j(j), m);
#endif
    ibuf = 0;
    buf[ibuf++] = mpi_index;
    buf[ibuf++] = m;
    buf[ibuf++] = j;
    for(k = 0; k < m; k++){
	buf[ibuf++] = (unsigned)ind[k];
	buf[ibuf++] = -1; // TODO_MPI: some real length??
    }
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    done = (char *)malloc(mpi_size * sizeof(char));
    // get current proc rank
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    done[mpi_rank] = 2;
    for(k = 1; k < mpi_size; k++){
	if(k != mpi_rank){
	    MPI_Send(buf, ibuf, MPI_UNSIGNED, k, MPI_ADD_TAG, MPI_COMM_WORLD);
	    done[k] = 1;
	}
    }
    nrecv = 0;
    for(k = 0; k < m; k++)
	// this is the local weight
	new_weight[k] = (isRowNull(mat, ind[k]) ? 0 : lengthRow(mat, ind[k]));
    finished = 0;
    while(!finished){
	// we loop until we receive all new weights and then we exit
#if DEBUG >= 1
	fprintf(stdout, "[m_a_r] Proc %d waiting... nrecv=%d\n", 
		mpi_rank, nrecv);
#endif
	MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	switch(status.MPI_TAG){
	case MPI_SEND_W_TAG:
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
# if DEBUG >= 1
	    fprintf(stdout, "[m_a_r] source=%u:", status.MPI_SOURCE);
	    for(k = 0; k < cnt; k++)
		fprintf(stdout, " %u", buf[k]);
	    fprintf(stdout, "\n");
# endif
	    for(k = 5; k < cnt; k += 2)
		new_weight[(k>>1)-1] += (int)buf[k+1];
	    if(done[status.MPI_SOURCE] == 1){
		done[status.MPI_SOURCE] = 2;
		nrecv++;
		finished = (nrecv == (mpi_size-2));
	    }
	    else
		fprintf(stderr, "MPI??\n");
	    break;
	default:
	    // should not happen!!!!
	    mpi_err1("tag received: %d\n", status.MPI_TAG);
	}
    }
    new_weight[0] = 0; // this is the destroyed row
    ibuf = 0;
    buf[ibuf++] = mpi_index;
    buf[ibuf++] = m;
    buf[ibuf++] = j;
    for(k = 0; k < m; k++){
	buf[ibuf++] = ind[k];
	buf[ibuf++] = new_weight[k];
    }
    MPI_Send(buf, ibuf, MPI_UNSIGNED, 0, MPI_SEND_W_TAG, MPI_COMM_WORLD);
#if DEBUG >= 1
    mpi_err("Exiting from mpi_add_rows\n");
#endif
    free(done);
}

// broadcast current history to other procs, so that they perform the same
// operations. Let's remind that this history can be rather large, so that
// some care is needed.
int
mpi_broadcast_history(report_t *rep, int mpi_rank)
{
    unsigned int buf[MPI_BUF_SIZE];
    int mpi_size, i, k, ibuf = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    // feed the buffer
    buf[ibuf++] = mpi_index;
    for(i = 0; i <= rep->mark; i++){
	if((ibuf + rep->history[i][0]) >= MPI_BUF_SIZE){
	    mpi_err("Too many lines in history\n");
	}
	for(k = 0; k <= rep->history[i][0]; k++)
	    buf[ibuf++] = rep->history[i][k];
    }
    // now send it
    for(k = 1; k < mpi_size; k++)
	if(k != mpi_rank)
	    MPI_Send(buf,ibuf,MPI_UNSIGNED,k,MPI_SEND_HIS_TAG,MPI_COMM_WORLD);
    return ibuf;
}

#define FULL_MONTY 0

// buf = [index, m, ...] and surely, we have m > 2.
void
mpi_MST(report_t *rep, sparse_mat_t *mat, int *njrem, unsigned int *buf, int mpi_rank, int mpi_size)
{
#if FULL_MONTY
    MPI_Status status;
    int nrecv, cnt;
#endif
    dclist dcl;
    double tMST;
    INT j, ind[MERGE_LEVEL_MAX];
    int k, ni, m, A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], r, s;

    // first, we need the rows
    m = (int)buf[1];
    dcl = mat->S[m]->next;
    j = dcl->j;
    ni = 0;
    for(k = 1; k <= mat->R[j][0]; k++){
	if(mat->R[j][k] != -1){
	    ind[ni++] = (INT)mat->R[j][k];
	    if(ni == m)
		break;
	}
    }
#if FULL_MONTY
    // second, we have to ask all procs for their shares
    for(r = 0; r < m; r++)
	buf[r+2] = (unsigned int)ind[r];
    for(k = 1; k < mpi_size; k++)
	if(k != mpi_rank)
	    MPI_Send(buf, m+2, MPI_UNSIGNED, k, 
		     MPI_ASK_MST_TAG, MPI_COMM_WORLD);
#endif // FULL_MONTY
    // now, compute my own share
    fillRowAddMatrix(A, mat, m, ind);
#if DEBUG >= 0
    mpi_err("Here is my MST submaster share:");
    for(r = 0; r < m; r++)
	for(s = r+1; s < m; s++)
	    fprintf(stderr, " %d", A[r][s]);
    fprintf(stderr, "\n");
#endif
#if FULL_MONTY
    // now, wait for the other shares
    nrecv = 0;
    while(nrecv != (mpi_size-2)){
	MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if(status.MPI_TAG == MPI_ASK_MST_TAG){
	    // buf = [index, m, A1, A2, ...]
	    nrecv++;
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
#if DEBUG >= 0
	    fprintf(stderr, "MST received from %d:", status.MPI_SOURCE);
	    mpi_err_tab("", buf, cnt);
#endif
	    k = 2;
	    for(r = 0; r < m; r++)
		for(s = r+1; s < m; s++)
		    A[r][s] += (int)buf[k++];
	}
    }
#endif // FULL_MONTY
    // then, we have to finish the MST business
    rep->mark = -1;
    MSTWithA(rep, mat, m, ind, &tMST, A);
    mat->rem_nrows--;
    mat->rem_ncols--;
    remove_j_from_SWAR(mat, j);
    *njrem = deleteHeavyColumns(rep, mat);
}

int
mpi_doOneMerge(report_t *rep, sparse_mat_t *mat, unsigned int *buf, int mpi_rank)
{
    MPI_Status status;
    double totopt = 0.0, totfill = 0.0, totMST = 0.0, totdel = 0.0;
    unsigned int *new_weight, index;
    int m, njdel, njrem = 0, verbose = 1, mpi_size, nrecv, finished, k, cnt;
    int nsup, nbrow;
    char *done;

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    // buf = [index, m]
    index = buf[0];
    m = (int)buf[1];
    // the local part
#if DEBUG >= 1
    mpi_err1("doOneMerge(%d)\n", m);
#endif
    if((m <= 2) || 0) // TODO: tmp
	doOneMerge(rep, mat, &njrem, &totopt, &totfill, &totMST,
		   &totdel, m, 10, 0, verbose); // 10 is rather arbitrary...!
    else
	mpi_MST(rep, mat, &njrem, buf, mpi_rank, mpi_size);
    njdel = deleteEmptyColumns(mat);
#if DEBUG >= 1
    if(njdel > 0)
	mpi_err1("I deleted %d empty columns\n", njdel);
#endif
    // the sub-master must store its history file...
    fprint_report(rep);
    // ... and must force other procs to perform the operations
    nsup = mpi_broadcast_history(rep, mpi_rank);
    // now wait for the partial weights to come back
    done = (char *)malloc(mpi_size * sizeof(char));
    memset(done, 0, mpi_size);
    done[mpi_rank] = 1;

    // let's overshoot a little
    nsup <<= 1;
    new_weight = (unsigned int *)malloc(nsup * sizeof(unsigned int));
    memset(new_weight, 0, nsup * sizeof(unsigned int));
    nrecv = 0;
    finished = 0;
    while(!finished){
	// we loop until we receive all new weights and then we exit
#if DEBUG >= 1
	fprintf(stderr, "[m_a_r] Proc %d waiting... nrecv=%d\n", 
		mpi_rank, nrecv);
#endif
	MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	switch(status.MPI_TAG){
	case MPI_SEND_W_TAG:
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
#if DEBUG >= 1
	    fprintf(stderr, "[m_a_r] source=%u:", status.MPI_SOURCE);
	    for(k = 0; k < cnt; k++)
		fprintf(stderr, " %u", buf[k]);
	    fprintf(stderr, "\n");
#endif
	    if(nrecv == 0){
		// first time: feed array
		memcpy(new_weight, buf+1, (cnt-1) * sizeof(unsigned int));
		nbrow = cnt-1;
#if DEBUG >= 1
		fprintf(stderr, "INIT:");
		for(k = 0; k < nbrow; k++)
		    fprintf(stderr, " %u", new_weight[k]);
		fprintf(stderr, "\n");
#endif
		if(nbrow > nsup)
		    mpi_err2("Gasp: nbrow=%d > nsup=%d\n", nbrow, nsup);
	    }
	    else
		// update
		for(k = 1; k < cnt; k += 2){
		    if(buf[k] != new_weight[k-1])
			mpi_err3("Pb[%d]: %d // %d\n",
				 k, (int)buf[k], (int)new_weight[k-1]);
		    new_weight[k] += buf[k+1];
		}
	    if(done[status.MPI_SOURCE] == 0){
		done[status.MPI_SOURCE] = 1;
		nrecv++;
		finished = (nrecv == (mpi_size-2));
	    }
	    else
		fprintf(stderr, "MPI??\n");
	    break;
	default:
	    // should not happen!!!!
	    mpi_err1("tag received: %d\n", status.MPI_TAG);
	}
    }
    for(k = 0; k < nbrow; k += 2)
	// this is the local weight
	new_weight[k+1] += (isRowNull(mat, new_weight[k]) ?
			    0 : lengthRow(mat, new_weight[k]));
    memcpy(buf+1, new_weight, nbrow * sizeof(unsigned int));
    MPI_Send(buf, nbrow+1, MPI_UNSIGNED, 0, MPI_SEND_W_TAG, MPI_COMM_WORLD);
    free(done);
    free(new_weight);
#if DEBUG >= 1
    mpi_err2("nrows=%d ncols=%d\n", mat->rem_nrows, mat->rem_ncols);
#endif
    return njdel;
}

int
mpi_replay_history(report_t *rep, sparse_mat_t *mat, unsigned int *send_buf, unsigned int *buf, int cnt)
{
    int r = 0, ind, k, kmax, i, i0, isb = 0, sgi0;

#if DEBUG >= 1
    fprintf(stderr, "DUMP_BEGIN [cnt=%d]\n", cnt);
    ind = buf[r++];
    while(1){
        kmax = buf[r++];
	fprintf(stderr, "%d", kmax);
	for(k = 0; k < kmax; k++)
	    fprintf(stderr, " %d", (int)buf[r++]);
	fprintf(stderr, "\n");
	if(r >= cnt)
	    break;
    }
    fprintf(stderr, "DUMP_END\n");
#endif
    rep->mark = -1;
    send_buf[isb++] = buf[0]; // index
    r = 0;
    ind = buf[r++];
    while(1){
	kmax = buf[r++];
	i0 = (int)buf[r++];
	// R[i0] should be added to R[1..kmax[
	if(i0 >= 0){
	    // and destroyed
	    sgi0 = 1;
	    send_buf[isb++] = i0;
	    send_buf[isb++] = 0;
	}
	else{
	    // not destroyed
	    sgi0 = -1;
	    i0 = -i0-1;
	}
	if(isRowNull(mat, i0))
	    for(k = 1; k < kmax; k++){
		send_buf[isb++] = buf[r++];
		send_buf[isb++] = (isRowNull(mat, buf[r-1]) ? 0 :
				   lengthRow(mat,  buf[r-1]));
	    }
	else{
	    for(k = 1; k < kmax; k++){
		i = (int)buf[r++];
		send_buf[isb++] = buf[r-1];
#if DEBUG >= 1
		mpi_err2("R[%d] += R[%d]\n", i, i0);
#endif
		if(isRowNull(mat, i)){
		    fprintf(stderr, "Row %d is null -- %p!!!!!\n", 
			    i, mat->rows[i]);
		    // just copy
		    mat->rows[i] = copyRow(mat->rows[3]);
		    addRowSWAR(mat, i);
		}
		else
		    addRowsSWAR(mat, i, i0, -1);
		send_buf[isb++] = lengthRow(mat, i);
	    }
	    if(sgi0 >= 0){
		// i0 is no longer needed...!
		removeRowSWAR(mat, i0);
		destroyRow(mat, i0);
	    }
	}
	if(r >= cnt)
	    break;
    }
    return isb;
}

// buf = [index, m, i1, ..., im]. Compute the partial sum of all possible
// row additions to be used later in the real MST stuff performed by the
// submaster. A priori, we will be sending A[r][s] for r < s only, so that
// we can gain on the communication load.
int
mpi_mst_share(unsigned int *send_buf, sparse_mat_t *mat, unsigned int *buf)
{
    INT ind[MERGE_LEVEL_MAX];
    int m = (int)buf[1], r, s, A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], ibuf;

    for(r = 0; r < m; r++)
	ind[r] = (INT)buf[r+2];
    fillRowAddMatrix(A, mat, m, ind);
    ibuf = 0;
    send_buf[ibuf++] = buf[0];
    send_buf[ibuf++] = buf[1];
    for(r = 0; r < m; r++)
	for(s = r+1; s < m; s++)
	    send_buf[ibuf++] = A[r][s];
    return ibuf;
}

void
mpi_slave(report_t *rep, int mpi_rank, sparse_mat_t *mat, FILE *purgedfile, char *purgedname)
{
    double totwait = 0.0, tt, wctstart = MPI_Wtime();
    unsigned int buf[MPI_BUF_SIZE];
    unsigned int send_buf[MPI_BUF_SIZE];
    MPI_Status status;
    int jmin, jmax, cnt, i, k, m, njdel = 0, i0, submaster, nbminm = 0;

    // loop forever
    while(1){
	tt = seconds();
#if DEBUG >= 1
	mpi_err("Waiting...\n");
#endif
	MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
		 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	totwait += (seconds()-tt);
#if DEBUG >= 1
	mpi_err3("Received %d from %d [index=%d]\n",
		 status.MPI_TAG, status.MPI_SOURCE, (int)buf[0]);
#endif
	// FIXME: sometimes, only 0 has the right to send things
	switch(status.MPI_TAG){
	case MPI_DIE_TAG:
	    mpi_err3("I was asked to die at %d: wct=%d wait=%d\n", 
		     (int)seconds(),(int)(MPI_Wtime()-wctstart),(int)totwait);
	    break;
	case MPI_J_TAG:
	    jmin = (int)buf[0];
	    jmax = (int)buf[1];
	    mpi_err2("A slave does not need to do too much things...\n... but has to treat C[%d..%d[\n", jmin, jmax);
	    initWeightFromFile(mat, purgedfile, jmin, jmax);
	    gzip_close(purgedfile, purgedname);
	    fillSWAR(mat, jmin, jmax);
	    purgedfile = gzip_open(purgedname, "r");
	    readmat(mat, purgedfile, jmin, jmax);
	    gzip_close(purgedfile, purgedname);
	    MPI_Send(buf, 1, MPI_UNSIGNED, 0, MPI_J_TAG, MPI_COMM_WORLD);
	    break;
	case MPI_MIN_M_TAG:
	    // buf = [index]
	    if(buf[0] >= mpi_index){
		mpi_index = buf[0];
#if DEBUG >= 1
		mpi_err1("New mpi_index = %d\n", mpi_index);
#endif
	    }
	    nbminm++;
	    if((nbminm % 10000) == 0)
		mpi_err2("WCT[%d]: %d\n", nbminm, (int)(MPI_Wtime()-wctstart));
	    // buf <- [index, m_min]
	    m = minColWeight(mat);
	    send_buf[0] = buf[0];
	    send_buf[1] = m;
	    send_buf[2] = mpi_get_number_of_active_colums(mat, jmin, jmax);
	    MPI_Send(send_buf, 3, MPI_UNSIGNED, 0, MPI_MIN_M_TAG, MPI_COMM_WORLD);
	    break;
	case MPI_DO_M_TAG:
	    rep->mark = -1;
	    njdel = mpi_doOneMerge(rep, mat, buf, mpi_rank);
	    send_buf[0] = buf[0];
	    send_buf[1] = buf[1];
	    send_buf[2] = njdel;
	    MPI_Send(send_buf, 3, MPI_UNSIGNED, 0, MPI_M_DONE_TAG, MPI_COMM_WORLD);
	    break;
	case MPI_ADD_TAG:
	    // buf = [index, m, j, i1, l1, i2, l2, ..., ir, lr]
	    // TODO_MPI: receive length of r1+r2???? Perhaps pbs with MST
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
#if DEBUG >= 1
	    fprintf(stderr, "From proc=%d, I received %d values:",
		    status.MPI_SOURCE, cnt);
	    fprintf(stderr, " index=%d m=%d j=%d", buf[0], buf[1], buf[2]);
	    for(k = 3; k < cnt; k++)
		fprintf(stderr, " %d", (int)buf[k]);
	    fprintf(stderr, "\n");
#endif
	    // if buf[3] is non null, we must add it to all buf[k > 3]
	    i0 = (int)buf[3];
	    // humffffff
	    if(!isRowNull(mat, i0)){
		for(k = 5; k < cnt; k += 2){
		    i = (int)buf[k];
#if DEBUG >= 1
		    mpi_err2("R[%d] += R[%d]\n", i, i0);
#endif
		    if(isRowNull(mat, i)){
			fprintf(stderr, "Row %d is null -- %p!!!!!\n", 
				i, mat->rows[i]);
			// just copy
			mat->rows[i] = copyRow(mat->rows[3]);
			addRowSWAR(mat, i);
		    }
		    else
			addRowsSWAR(mat, i, i0, (int)buf[k+1]);
		}
		// i0 is no longer needed...!
		removeRowSWAR(mat, i0);
		destroyRow(mat, i0);
	    }
	    // now, send back the new weights of the rows (skipping 2...)
	    memcpy(send_buf, buf, cnt * sizeof(unsigned int));
	    send_buf[4] = 0;
	    for(i = 5; i < cnt; i += 2){
		if(isRowNull(mat, (int)buf[i]))
		    send_buf[i+1] = 0;
		else
		    send_buf[i+1] = (unsigned)lengthRow(mat, (int)buf[i]);
#if DEBUG >= 1
		fprintf(stderr, "index=%u, m=%u j=%u", buf[0], buf[1], buf[2]);
		fprintf(stderr, " i=%u => w=%u\n", send_buf[i], send_buf[i+1]);
#endif
	    }
	    // back to the current sub-master
	    // send_buf = [index, m, j, i1, 0, i2, new_l2, ..., ir, newlr]
	    MPI_Send(send_buf, cnt, MPI_UNSIGNED, status.MPI_SOURCE,
		     MPI_SEND_W_TAG,MPI_COMM_WORLD);
	    break;
	case MPI_INACTIVATE_ROWS_TAG:
	    // buf = [index, i1, i2, ..., ir]
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
	    rep->mark = -1;
	    mpi_inactivate_rows(rep, mat, buf, cnt);
	    if((cnt-1) != (rep->mark+1)){
		fprintf(stderr, "At the end of inactivate\n");
		fprint_report_aux(stderr, rep);
	    }
	    MPI_Send(buf, 1, MPI_UNSIGNED, 0, 
		     MPI_INACTIVATE_ROWS_TAG, MPI_COMM_WORLD);
	    break;
	case MPI_SEND_HIS_TAG:
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
	    submaster = status.MPI_SOURCE;
	    rep->mark = -1;
	    cnt = mpi_replay_history(rep, mat, send_buf, buf, cnt);
	    // send_buf = [index, i1, w(i1), i2, w(i2), ..., ir, w(ir)]
#if DEBUG >= 1
	    fprintf(stderr, "I have to send %d items:", cnt);
	    mpi_err_tab("", send_buf, cnt);
#endif
	    MPI_Send(send_buf, cnt, MPI_UNSIGNED, submaster,
                     MPI_SEND_W_TAG, MPI_COMM_WORLD);
	    break;
	case MPI_ASK_MST_TAG:
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
            submaster = status.MPI_SOURCE;
	    // buf = [index, m, i1, ..., im]
	    cnt = mpi_mst_share(send_buf, mat, buf);
#if DEBUG >= 0
	    mpi_err_tab("MST share:", send_buf, cnt);
#endif
	    MPI_Send(send_buf, cnt, MPI_UNSIGNED, submaster,
                     MPI_ASK_MST_TAG, MPI_COMM_WORLD);
	    break;
	default:
	    mpi_err1("Unknown tag: %d\n", status.MPI_TAG);
	    break;
	}
	if(status.MPI_TAG == MPI_DIE_TAG)
	    break;
    }
}

void
mpi_start_slaves(int mpi_size, sparse_mat_t *mat)
{
    MPI_Status status;
    double wctstart = MPI_Wtime();
    unsigned int buf[MPI_BUF_SIZE];
    int jmin, jmax, kappa;
    int rk;

    mpi_tab_j = (int *)malloc(mpi_size * sizeof(int));
#if 1 // V1
    kappa = mat->ncols/(mpi_size-1);
#endif
    for(rk = 1; rk < mpi_size; rk++){
#if 1 // V1
	// let kappa = ncols/(mpi_size-1)
	// slave_j for j < mpi_size-1 will have to treat [jmin..jmax[
	// jmin = (j-1)*kappa, jmax = j*kappa-1;
	// for j = mpi_size-1, this will be jmax = ncols
	jmin = (rk-1) * kappa;
	jmax = (rk == (mpi_size-1) ? mat->ncols : rk * kappa);
#else // PB: we would need wt[] just for this...!
	jmin = (rk == 1 ? 0 : jmax);
	jmax = jmin;
	kappa = 0;
	for(jmax = jmin; jmax < mat->ncols; jmax++){
	    kappa += abs(mat->wt[jmax]);
# if DEBUG >= 1
	    fprintf(stderr, "wt=%d kappa=%d\n", mat->wt[jmax], kappa);
# endif
	    if(((unsigned long)kappa) > (mat->weight / ((unsigned long)((mpi_size-1)))))
		break;
	}
	if(rk == (mpi_size-1))
	    // anyway...
	    jmax = mat->ncols;
#endif
	mpi_err3("Weight[%d..%d[ = %d\n", jmin, jmax, kappa);
	mpi_tab_j[rk-1] = jmin;
	buf[0] = (unsigned)jmin;
	buf[1] = (unsigned)jmax;
        MPI_Send(buf, 2, MPI_UNSIGNED, rk, MPI_J_TAG, MPI_COMM_WORLD);
	while(1){
	    MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, rk,
		     MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	    if(status.MPI_TAG == MPI_J_TAG)
		break;
	}
	mpi_err3("#T# started processor %d at %d [wct=%d]\n", 
		 rk, (int)seconds(), (int)(MPI_Wtime()-wctstart));
    }
    mpi_tab_j[mpi_size-1] = jmax;
    fprintf(stderr, "# mpi_tab_j:");
    for(rk = 0; rk < mpi_size; rk++)
	fprintf(stderr, " %d", mpi_tab_j[rk]);
    fprintf(stderr, "\n");
}

int
mpi_get_minimal_m_proc(int *m, int *proc, int mpi_size, unsigned int index)
{
    MPI_Status status;
    unsigned int buf[MPI_BUF_SIZE];
    int k, nrecv, curr_ncols = 0;

    *m = MERGE_LEVEL_MAX, nrecv;
    buf[0] = index;
    for(k = 1; k < mpi_size; k++)
	MPI_Send(buf, 1, MPI_UNSIGNED, k, MPI_MIN_M_TAG, MPI_COMM_WORLD);
    *proc = -1;
    nrecv = 0;
    while(1){
	MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#if DEBUG >= 2
	fprintf(stderr, "source=%u: index=%u m_min=%u\n",
		status.MPI_SOURCE, buf[0], buf[1]);
#endif
	if(status.MPI_TAG == MPI_MIN_M_TAG){
	    if(buf[0] != index){
		fprintf(stderr, 
			"Should not happen: bad index %u instead of %u\n", 
			buf[0], index);
	    }
	    else
		nrecv++;
	    if(*m > (int)buf[1]){
		*proc = status.MPI_SOURCE;
		*m = (int)buf[1];
	    }
	    curr_ncols += (int)buf[2];
	    if(nrecv == (mpi_size-1))
		break;
	}
    }
    return curr_ncols;
}

void
mpi_feed(int *row_weight, sparse_mat_t *mat, FILE *purgedfile)
{
    int i, j, nc, ret, x;

    mat->rem_nrows = mat->nrows;
    mat->rem_ncols = mat->ncols;
    for(i = 0; i < mat->nrows; i++){
	ret = fscanf(purgedfile, "%d", &j); // unused index to rels file
	ASSERT_ALWAYS (ret == 1);
	ret = fscanf (purgedfile, "%d", &nc);
	ASSERT_ALWAYS (ret == 1);
	for(j = 0; j < nc; j++)
	    ret = fscanf(purgedfile, PURGE_INT_FORMAT, &x);
	row_weight[i] = nc;
    }
}

void
mpi_delete_superfluous_rows(report_t *rep, sparse_mat_t *mat, int *row_weight, int mpi_size)
{
    MPI_Status status;
    unsigned int send_buf[MPI_BUF_SIZE];
    int k, r, ni2rem = number_of_superfluous_rows(mat), *tmp, ntmp, nrecv;

    if((mat->rem_nrows - mat->rem_ncols) <= mat->delta)
	return;
    if(ni2rem > MPI_BUF_SIZE-1)
	ni2rem = MPI_BUF_SIZE-1;
    rep->mark = -1;
    mpi_index++;
#if DEBUG >= 0
    mpi_err1("We should drop %d rows\n", ni2rem);
#endif
    tmp = (int *)malloc((mat->nrows << 1) * sizeof(int));
    for(k = 0, ntmp = 0; k < mat->nrows; k++)
	if(row_weight[k] > 0){
	    tmp[ntmp++] = row_weight[k];
	    tmp[ntmp++] = k;
	}
    qsort(tmp, ntmp>>1, 2 * sizeof(int), cmp);
    for(k = ntmp-1, r = 1; k >= 0; k -= 2){
	send_buf[r++] = tmp[k];
	row_weight[tmp[k]] = 0;
	report1(rep, tmp[k]);
	mat->rem_nrows--;
	if((r >= ni2rem) || ((mat->rem_nrows - mat->rem_ncols) < mat->delta))
	    break;
    }
    free(tmp);
    fprint_report(rep);
    for(k = 1; k < mpi_size; k++){
	mpi_index++;
	send_buf[0] = mpi_index;
#if DEBUG >= 1
	mpi_err1("Sending order to %d\n", k);
#endif
	MPI_Send(send_buf, r, MPI_UNSIGNED, k, 
		 MPI_INACTIVATE_ROWS_TAG, MPI_COMM_WORLD);
    }
    // this is kinda wait, no?
    nrecv = 0;
    while(nrecv != (mpi_size-1)){
#if DEBUG >= 1
	mpi_err1("Waiting in superfluous, nrecv=%d\n", nrecv);
#endif
	MPI_Recv(send_buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
		 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if(status.MPI_TAG == MPI_INACTIVATE_ROWS_TAG)
	    nrecv++;
    }
}

// actually, mat is rather empty, since it does not use too much fancy things.
// So we just need to init the row weights.
void
mpi_master(report_t *rep, sparse_mat_t *mat, int mpi_size, FILE *purgedfile)
{
    MPI_Status status;
    double totwait = 0.0, tt, wctstart = MPI_Wtime();
    unsigned int buf[MPI_BUF_SIZE];
    unsigned int send_buf[MPI_BUF_SIZE];
    unsigned int threshold = 0;
    int ibuf, m, proc, cnt, k, done, maxm = 0, nrecv;
    int *row_weight, curr_ncols;

#if 0
    mpi_err("WARNING: maxm=2\n");
    mat->mergelevelmax = 2;
#endif
    // let's initialize a new simpler array with the weights of the rows
    row_weight = (int *)malloc(mat->nrows * sizeof(int));
    mpi_feed(row_weight, mat, purgedfile);
    while(1){
	mpi_index++;
	tt = seconds();
	curr_ncols = mpi_get_minimal_m_proc(&m, &proc, mpi_size, mpi_index);
	totwait += (seconds()-tt);
#if DEBUG >= 1
	if(curr_ncols != mat->rem_ncols){
	    fprintf(stderr, "index=%d cur=%d ncols=%d\n", 
		    mpi_index, curr_ncols, mat->rem_ncols);
# if 0
	    sleep(5);
	    mpi_kill_slaves();
	    MPI_Finalize();
	    exit(0);
# endif
	}
#endif
	if(mpi_index >= threshold){
	    threshold += 10000;
	    fprintf(stderr, 
		    "R%d (wct=%2.2lf): maxm=%d nrows=%d ncols=%d[%d] (%d)\n",
		    mpi_index, MPI_Wtime()-wctstart, maxm, 
		    mat->rem_nrows, mat->rem_ncols, curr_ncols,
		    mat->rem_nrows - mat->rem_ncols);
	    mpi_delete_superfluous_rows(rep, mat, row_weight, mpi_size);
	}
	// look for minimal m for index
	//	mat->rem_ncols = curr_ncols;
#if DEBUG >= 1
	mpi_err2("Minimal m is %d for proc=%d\n", m, proc);
#endif
	if(m > mat->mergelevelmax){
	    mpi_err("Breaking\n");
	    break;
	}
	if(m > maxm)
	    maxm = m;
	// ask the proc to execute one step
	send_buf[0] = mpi_index;
	send_buf[1] = m;
	MPI_Send(send_buf, 2, MPI_UNSIGNED, proc, MPI_DO_M_TAG, MPI_COMM_WORLD);
	while(1){
#if DEBUG >= 1
	    mpi_err("Waiting...\n");
#endif
	    done = 0;
	    tt = seconds();
	    MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
		     MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	    totwait += (seconds()-tt);
	    switch(status.MPI_TAG){
	    case MPI_M_DONE_TAG:
		// buf = [index, m, njdel]
#if DEBUG >= 1
		mpi_err2("We stop here: m=%d njdel=%d\n",
			 (int)buf[1], (int)buf[2]);
#endif
		mat->rem_ncols -= (int)buf[2];
		done = 1;
		break;
	    case MPI_SEND_W_TAG:
		// new rows received after additions
		MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
		// buf = [index, i1, new_l1, i2, new_l2, ..., ir, new_lr]
		for(k = 1; k < cnt; k += 2){
#if DEBUG >= 1
		    mpi_err3("New row %d has now weight %d instead of %d\n",
			     buf[k], buf[k+1], row_weight[buf[k]]);
#endif
		    row_weight[buf[k]] = buf[k+1];
		    if(buf[k+1] == 0){
			mat->rem_nrows--;
			mat->rem_ncols--;
		    }
		}
		memcpy(send_buf, buf, cnt);
		ibuf = 1;
#if 1
		// TODO: here, it is not the master who might decide
		// but all procs at the same time by sharing the
		// row_weight info???????????
		if((mat->rem_nrows - mat->rem_ncols) > mat->delta){
		    // look for too heavy rows
		    rep->mark = -1;
		    mpi_index++;
		    for(k = 1; k < cnt; k += 2){
			if(row_weight[buf[k]] > mat->rwmax){
#if DEBUG >= 1
			    mpi_err3("Row %d is too heavy (%d/%d)\n",
				     buf[k], row_weight[buf[k]], mat->rwmax);
#endif
			    send_buf[ibuf++] = buf[k];
			    row_weight[buf[k]] = 0;
			    report1(rep, (int)buf[k]);
			    mat->rem_nrows--;
			}
		    }
		}
#endif
		if(ibuf > 1){
		    fprint_report(rep);
#if DEBUG >= 1
		    mpi_err1("Sending %d heavy row(s) to be deleted\n",ibuf-1);
#endif
		    for(k = 1; k < mpi_size; k++){
			mpi_index++;
			send_buf[0] = mpi_index;
#if DEBUG >= 1
			mpi_err1("Sending order to %d\n", k);
#endif
			MPI_Send(send_buf, ibuf, MPI_UNSIGNED, k, 
				 MPI_INACTIVATE_ROWS_TAG, MPI_COMM_WORLD);
		    }
		    done = nrecv = 0;
		}
		break;
	    case MPI_INACTIVATE_ROWS_TAG:
		nrecv++;
		if(nrecv == (mpi_size-1))
		    done = 1;
		break;
	    default:
		fprintf(stderr, "Received tag was %d\n", status.MPI_TAG);
	    }
	    if(done)
		break;
	}
    }
    //    fprint_report(rep);// what for?
    fprintf(stderr, "Finally (%d -- %2.2lf -- wait=%2.2lf):",
	    mpi_index, seconds(), totwait);
    fprintf(stderr, " nrows=%d ncols=%d[%d]\n",
	    mat->rem_nrows, mat->rem_ncols, curr_ncols);
    fflush(stderr);
    free(row_weight);
}

void
mpi_start_proc(char *outname, int mpi_rank, int mpi_size, sparse_mat_t *mat, FILE *purgedfile, char *purgedname)
{
    report_t rep;
    char *str;
    char mpitmp[1024];

    if(mpi_rank != 0){
	// redirecting stderr; one day, redirect the master?
	sprintf(mpitmp, "%s%03d.err", outname, mpi_rank);
	freopen(mpitmp, "w", stderr);
    }
    // treat the .gz case
    str = (char *)malloc(strlen(outname)+3);
    if(is_gzip(outname)){
	printf("Case of gz file NYI\n");
	exit(0);
    }
    else
	sprintf(str, "%s%03d", outname, mpi_rank);
    fprintf(stderr, "Outfile for proc=%d will be %s\n", mpi_rank, str);
    mpi_index = 0;
    init_rep(&rep, str, mat, 1);
    if(mpi_rank == 0){
	fprintf(rep.outfile, "0 %d %d\n", mat->nrows, mat->ncols);
	mpi_start_slaves(mpi_size, mat);
	mpi_master(&rep, mat, mpi_size, purgedfile);
	mpi_kill_slaves();
    }
    else
	mpi_slave(&rep, mpi_rank, mat, purgedfile, purgedname);
    gzip_close(rep.outfile, str);
    free(str);
}
