#ifndef	FAKEMPI_H_
#define	FAKEMPI_H_

#include <stdlib.h>
#include <string.h>

typedef int MPI_Status;
typedef int MPI_Datatype;
typedef int MPI_Comm;
typedef int MPI_Op;
typedef int MPI_User_function;

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



#endif /* FAKEMPI_H_ */
