#ifndef	FAKEMPI_H_
#define	FAKEMPI_H_

#include <stdlib.h>
#include <string.h>
#include "cado_mpi_config.h"
#include "macros.h"

typedef int MPI_Status;
typedef int MPI_Datatype;
typedef int MPI_Comm;
typedef int MPI_Op;
typedef int MPI_User_function;
typedef int MPI_Errhandler;
typedef int MPI_Request;

// type keys are sizeof() values.
#define MPI_BYTE        1
#define MPI_INT         sizeof(int)
#define MPI_DOUBLE      sizeof(double)
#define MPI_UNSIGNED_LONG sizeof(unsigned long)
#define MPI_UNSIGNED_INT  sizeof(unsigned int)
#define MPI_UNSIGNED      sizeof(unsigned int)

#define MPI_COMM_WORLD	0

#define MPI_THREAD_SINGLE       0
#define MPI_THREAD_MULTIPLE     3

/* We define different ops, but since they're collected amongst
 * communicators of size 1 anyway, the operation does not matter much...
 * */
#define MPI_BXOR       0
#define MPI_SUM        1
#define MPI_MAX        2
#define MPI_LAND       3

#define MPI_ERRORS_ARE_FATAl        0
#define MPI_ERRORS_RETURN        1

#define MPI_IN_PLACE    0

#define MPI_MAX_OBJECT_NAME     64

#define MPI_STATUS_IGNORE       0

/* That would be quite neat, but it's too coarse (it applies to the whole
 * translation unit), and anyway it isn't supported by oldish gcc's.
#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
*/

static inline int MPI_Wait(MPI_Request *request MAYBE_UNUSED, MPI_Status *status MAYBE_UNUSED) { return 0; }
static inline int MPI_Abort(MPI_Comm comm MAYBE_UNUSED, int s) { exit(s); }
static inline int MPI_Comm_rank(int s MAYBE_UNUSED, int  * p) { *p=0; return 0;}
static inline int MPI_Comm_size(int s MAYBE_UNUSED, int  * p) { *p=1; return 0;}
static inline int MPI_Initialized(int  * p) { *p=1; return 0; }
static inline int MPI_Init(int * argc MAYBE_UNUSED, char *** argv MAYBE_UNUSED) { return 0; }
static inline int MPI_Init_thread(int * argc MAYBE_UNUSED, char *** argv MAYBE_UNUSED, int req, int * prov) { if (prov) *prov=req; return 0; }
static inline int MPI_Finalize() {return 0;}
static inline int MPI_Op_create( MPI_User_function *function MAYBE_UNUSED, int commute MAYBE_UNUSED, MPI_Op *op MAYBE_UNUSED ){return 0;}
static inline int MPI_Send( void *buf MAYBE_UNUSED, int count MAYBE_UNUSED, MPI_Datatype datatype MAYBE_UNUSED, int dest MAYBE_UNUSED,int tag MAYBE_UNUSED, MPI_Comm comm  MAYBE_UNUSED){return 0;}
static inline int MPI_Isend( void *buf MAYBE_UNUSED, int count MAYBE_UNUSED, MPI_Datatype datatype MAYBE_UNUSED, int dest MAYBE_UNUSED,int tag MAYBE_UNUSED, MPI_Comm comm  MAYBE_UNUSED, MPI_Request * zz MAYBE_UNUSED){return 0;}
static inline int MPI_Recv( void *buf MAYBE_UNUSED, int count MAYBE_UNUSED, MPI_Datatype datatype MAYBE_UNUSED, int source MAYBE_UNUSED,int tag MAYBE_UNUSED, MPI_Comm comm MAYBE_UNUSED, MPI_Status *status MAYBE_UNUSED ){ abort(); return 0;}
static inline int MPI_Irecv( void *buf MAYBE_UNUSED, int count MAYBE_UNUSED, MPI_Datatype datatype MAYBE_UNUSED, int source MAYBE_UNUSED,int tag MAYBE_UNUSED, MPI_Comm comm MAYBE_UNUSED, MPI_Status *status MAYBE_UNUSED, MPI_Request * zz MAYBE_UNUSED){ abort(); return 0;}
static inline int MPI_Bcast( void *buffer MAYBE_UNUSED, int count MAYBE_UNUSED, MPI_Datatype datatype MAYBE_UNUSED, int root MAYBE_UNUSED, MPI_Comm comm MAYBE_UNUSED){return 0;}
static inline int MPI_Reduce ( void *sendbuf, void *recvbuf, int count,MPI_Datatype datatype, MPI_Op op MAYBE_UNUSED, int root MAYBE_UNUSED, MPI_Comm comm  MAYBE_UNUSED)
{
    if (sendbuf) memcpy(recvbuf, sendbuf, count * datatype);
    return 0;
}
static inline int MPI_Allreduce ( void *sendbuf, void *recvbuf, int count,MPI_Datatype datatype, MPI_Op op MAYBE_UNUSED, MPI_Comm comm  MAYBE_UNUSED)
{
    if (sendbuf) memcpy(recvbuf, sendbuf, count * datatype);
    return 0;
}
static inline int MPI_Comm_split (MPI_Comm x MAYBE_UNUSED, int color MAYBE_UNUSED, int key MAYBE_UNUSED, MPI_Comm * y MAYBE_UNUSED)
{
    return 0;
}
static inline int MPI_Comm_set_errhandler (MPI_Comm x MAYBE_UNUSED, MPI_Errhandler e MAYBE_UNUSED)
{
    return 0;
}
static inline int MPI_Comm_free (MPI_Comm * x MAYBE_UNUSED) { return 0; }
static inline int MPI_Comm_dup (MPI_Comm y, MPI_Comm * x) { *x = y; return 0; }
static inline int MPI_Comm_set_name(MPI_Comm comm MAYBE_UNUSED, char *comm_name MAYBE_UNUSED) { return 0;}
static inline int MPI_Comm_get_name(MPI_Comm comm MAYBE_UNUSED, char *comm_name MAYBE_UNUSED, int * rlen) { *comm_name='\0'; *rlen=0; return 0;}
static inline int MPI_Scatterv(void * sendbuf, int * sendcounts, int * displs,  MPI_Datatype st, void * recvbuf, int recvcount, MPI_Datatype rt, int root MAYBE_UNUSED, MPI_Comm x MAYBE_UNUSED) {
    ASSERT_ALWAYS(sendcounts[0] * st == recvcount * rt);
    memcpy(recvbuf, ((char *)sendbuf) + displs[0] * st, recvcount * rt);
    return 0;
}

static inline int MPI_Barrier (MPI_Comm x MAYBE_UNUSED) { return 0; }

static inline int MPI_Gatherv(void * sendbuf, int sendcount,  MPI_Datatype st, void * recvbuf, int * recvcounts, int * displs, MPI_Datatype rt, int root MAYBE_UNUSED, MPI_Comm x MAYBE_UNUSED) {
    ASSERT_ALWAYS(sendcount * st == recvcounts[0] * rt);
    memcpy(((char *)recvbuf) + displs[0] * rt, sendbuf, sendcount * st);
    return 0;
}
static inline int MPI_Allgather(void * sendbuf, int sendcount,  MPI_Datatype st, void * recvbuf, int recvcount, MPI_Datatype rt, MPI_Comm x MAYBE_UNUSED) {
    ASSERT_ALWAYS(sendcount * st == recvcount * rt);
    if (sendbuf) memcpy(recvbuf, sendbuf, sendcount * st);
    return 0;
}

#endif /* FAKEMPI_H_ */
