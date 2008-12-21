#ifndef	FAKEMPI_H_
#define	FAKEMPI_H_

#include <stdlib.h>
#include <string.h>

typedef int MPI_Status;
typedef int MPI_Datatype;
typedef int MPI_Comm;
typedef int MPI_Op;
typedef int MPI_User_function;

// type keys are sizeof() values.
#define MPI_BYTE        1
#define MPI_UNSIGNED_LONG       sizeof(unsigned long)

#define MPI_COMM_WORLD	0

#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

static inline int MPI_Comm_rank(int s, int  * p) { *p=0; return 0;}
static inline int MPI_Comm_size(int s, int  * p) { *p=1; return 0;}
static inline int MPI_Initialized(int  * p) { *p=1; return 0; }
static inline int MPI_Init(int * argc, char *** argv) { return 0; }
static inline int MPI_Finalize() {return 0;}
static inline int MPI_Op_create( MPI_User_function *function, int commute, MPI_Op *op ){return 0;}
static inline int MPI_Send( void *buf, int count, MPI_Datatype datatype, int dest,int tag, MPI_Comm comm ){return 0;}
static inline int MPI_Recv( void *buf, int count, MPI_Datatype datatype, int source,int tag, MPI_Comm comm, MPI_Status *status ){ abort(); return 0;}
static inline int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){return 0;}
static inline int MPI_Reduce ( void *sendbuf, void *recvbuf, int count,MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm )
{
    memcpy(recvbuf, sendbuf, count * datatype);
    return 0;
}
static inline int MPI_Comm_split (MPI_Comm x, int color, int key, MPI_Comm * y)
{
    return 0;
}
static inline int MPI_Comm_free (MPI_Comm * x) { return 0; }
static inline int MPI_Comm_dup (MPI_Comm y, MPI_Comm * x) { *x = y; return 0; }
static inline int MPI_Scatterv(void * sendbuf, int * sendcounts, int * displs,  MPI_Datatype st, void * recvbuf, int recvcount, MPI_Datatype rt, int root, MPI_Comm x) {
    ASSERT_ALWAYS(sendcounts[0] * st == recvcount * rt);
    memcpy(recvbuf, ((char *)sendbuf) + displs[0] * st, recvcount * rt);
    return 0;
}

static inline int MPI_Barrier (MPI_Comm x) { return 0; }

static inline int MPI_Gatherv(void * sendbuf, int sendcount,  MPI_Datatype st, void * recvbuf, int * recvcounts, int * displs, MPI_Datatype rt, int root, MPI_Comm x) {
    ASSERT_ALWAYS(sendcount * st == recvcounts[0] * rt);
    memcpy(((char *)recvbuf) + displs[0] * rt, sendbuf, sendcount * st);
    return 0;
}

#endif /* FAKEMPI_H_ */
