#ifndef	FAKEMPI_H_
#define	FAKEMPI_H_

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "cado_mpi_config.h"
#include "macros.h"
#include "portability.h"

typedef int MPI_Status;
typedef int MPI_Datatype;
typedef int MPI_Comm;
typedef int MPI_Op;
typedef void MPI_User_function(void *invec, void *inoutvec,
             int *len, MPI_Datatype *datatype);
typedef int MPI_Errhandler;
typedef int MPI_Request;

// type keys are sizeof() values.
#define MPI_DATATYPE_NULL       0
#define MPI_BYTE        1
#define MPI_INT         sizeof(int)
#define MPI_DOUBLE      sizeof(double)
#define MPI_UNSIGNED_LONG sizeof(unsigned long)
#define MPI_LONG        sizeof(long)
/* It seems that MPI_UNSIGNED_INT is in fact unspecified */
// #define MPI_UNSIGNED_INT  sizeof(unsigned int)
#define MPI_UNSIGNED      sizeof(unsigned int)

#define MPI_Type_size(x, s)     *(s)=(x)

#define MPI_COMM_WORLD	0

#define MPI_THREAD_SINGLE       0
#define MPI_THREAD_MULTIPLE     3

/* We define different ops, but since they're collected amongst
 * communicators of size 1 anyway, the operation does not matter much...
 * */
#define MPI_BXOR       0
#define MPI_SUM        1
#define MPI_MAX        2
#define MPI_MIN        3
#define MPI_LAND       4
#define MPI_BAND       5
#define MPI_BOR        6

#define MPI_ERRORS_ARE_FATAl        0
#define MPI_ERRORS_RETURN        1

#define MPI_IN_PLACE    0

#define MPI_MAX_OBJECT_NAME     64

#define MPI_STATUS_IGNORE       0

#define MPI_ANY_TAG     -1

/* That would be quite neat, but it's too coarse (it applies to the whole
 * translation unit), and anyway it isn't supported by oldish gcc's.
*/
#ifdef  __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

static inline int MPI_Wait(MPI_Request *request , MPI_Status *status ) { return 0; }
static inline int MPI_Abort(MPI_Comm comm , int s) { exit(s); }
static inline int MPI_Comm_rank(int s , int  * p) { *p=0; return 0;}
static inline int MPI_Comm_size(int s , int  * p) { *p=1; return 0;}
static inline int MPI_Initialized(int  * p) { *p=1; return 0; }
static inline int MPI_Init(int * argc , char *** argv ) { return 0; }
static inline int MPI_Init_thread(int * argc , char *** argv , int req, int * prov) { if (prov) *prov=req; return 0; }
static inline int MPI_Finalize() {return 0;}
static inline int MPI_Op_create( MPI_User_function *function , int commute , MPI_Op *op  ){return 0;}
static inline int MPI_Op_free(MPI_Op *op  ){return 0;}
static inline int MPI_Send( void *buf , int count , MPI_Datatype datatype , int dest ,int tag , MPI_Comm comm  ){return 0;}
static inline int MPI_Sendrecv( void *sbuf , int scount , MPI_Datatype sdatatype , int sdest ,int stag ,  void *rbuf , int rcount , MPI_Datatype rdatatype , int rdest ,int rtag , MPI_Comm comm  , MPI_Status *status ){return 0;}
static inline int MPI_Isend( void *buf , int count , MPI_Datatype datatype , int dest ,int tag , MPI_Comm comm  , MPI_Request * zz ){return 0;}
static inline int MPI_Recv( void *buf , int count , MPI_Datatype datatype , int source ,int tag , MPI_Comm comm , MPI_Status *status  ){ abort(); return 0;}
static inline int MPI_Irecv( void *buf , int count , MPI_Datatype datatype , int source ,int tag , MPI_Comm comm , MPI_Status *status , MPI_Request * zz ){ abort(); return 0;}
static inline int MPI_Bcast( void *buffer , int count , MPI_Datatype datatype , int root , MPI_Comm comm ){return 0;}
static inline int MPI_Reduce ( void *sendbuf, void *recvbuf, int count,MPI_Datatype datatype, MPI_Op op , int root , MPI_Comm comm  )
{
    if (sendbuf) memcpy(recvbuf, sendbuf, count * datatype);
    return 0;
}
static inline int MPI_Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts,
                    MPI_Datatype datatype, MPI_Op op , MPI_Comm comm )
{
    if (sendbuf) memcpy(recvbuf, sendbuf, recvcounts[0] * datatype);
    return 0;
}
static inline int MPI_Allreduce ( void *sendbuf, void *recvbuf, int count,MPI_Datatype datatype, MPI_Op op , MPI_Comm comm  )
{
    if (sendbuf) memcpy(recvbuf, sendbuf, count * datatype);
    return 0;
}
static inline int MPI_Comm_split (MPI_Comm x , int color , int key , MPI_Comm * y)
{
    *y=0;
    return 0;
}
static inline int MPI_Comm_set_errhandler (MPI_Comm x , MPI_Errhandler e )
{
    return 0;
}
static inline int MPI_Comm_free (MPI_Comm * x ) { return 0; }
static inline int MPI_Comm_dup (MPI_Comm y, MPI_Comm * x) { *x = y; return 0; }
static inline int MPI_Comm_set_name(MPI_Comm comm , char *comm_name ) { return 0;}
static inline int MPI_Comm_get_name(MPI_Comm comm , char *comm_name , int * rlen) { *comm_name='\0'; *rlen=0; return 0;}
static inline int MPI_Scatterv(void * sendbuf, int * sendcounts, int * displs,  MPI_Datatype st, void * recvbuf, int recvcount, MPI_Datatype rt, int root , MPI_Comm x ) {
    ASSERT_ALWAYS(sendcounts[0] * st == recvcount * rt);
    memcpy(recvbuf, ((char *)sendbuf) + displs[0] * st, recvcount * rt);
    return 0;
}

static inline int MPI_Scatter(void * sendbuf, int sendcount, MPI_Datatype st, void * recvbuf, int recvcount, MPI_Datatype rt, int root , MPI_Comm x ) {
    ASSERT_ALWAYS(sendcount * st == recvcount * rt);
    if (recvbuf && sendbuf)
        memcpy(recvbuf, sendbuf, recvcount * rt);
    return 0;
}

static inline int MPI_Barrier (MPI_Comm x ) { return 0; }

static inline int MPI_Gather(void * sendbuf, int sendcount,  MPI_Datatype st, void * recvbuf, int recvcount, MPI_Datatype rt, int root , MPI_Comm x ) {
    if (sendbuf == MPI_IN_PLACE) return 0;
    memcpy(((char *)recvbuf), (char*) sendbuf, sendcount * st);
    return 0;
}
static inline int MPI_Gatherv(void * sendbuf, int sendcount,  MPI_Datatype st, void * recvbuf, int * recvcounts, int * displs, MPI_Datatype rt, int root , MPI_Comm x ) {
    ASSERT_ALWAYS(sendcount * st == recvcounts[0] * rt);
    memcpy(((char *)recvbuf) + displs[0] * rt, sendbuf, sendcount * st);
    return 0;
}
static inline int MPI_Allgather(void * sendbuf , int sendcount ,  MPI_Datatype st , void * recvbuf , int recvcount , MPI_Datatype rt , MPI_Comm x ) {
    ASSERT_ALWAYS(sendbuf == MPI_IN_PLACE || sendcount * st == recvcount * rt);
    if (sendbuf) memcpy(recvbuf, sendbuf, sendcount * st);
    return 0;
}
static inline int MPI_Allgatherv(void *sendbuf, int sendcount ,
            MPI_Datatype sendtype , void *recvbuf , int *recvcount ,
            int *displs , MPI_Datatype recvtype , MPI_Comm comm )
{
    ASSERT_ALWAYS(sendbuf == MPI_IN_PLACE);
    return 0;
}

static inline int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype * newtype)
{
    *newtype = count * oldtype;
    return 0;
}
static inline int MPI_Type_commit(MPI_Datatype * t ) { return 0; }
static inline int MPI_Type_free(MPI_Datatype * t ) { return 0; }
static inline int MPI_Type_set_attr(MPI_Datatype type , int key , void *value )
{
    /* XXX We are *NOT* storing the type attribute here. This is because
     * we expect that the user function that we use, and which exploits
     * this, will actually *never* be called in this context of ``fake''
     * mpi: by assumption, all communicator sizes are equal to 1, and
     * thus mpi reduction just amounts to copying data.
     */
    return 0;
}
static inline int MPI_Type_get_attr(MPI_Datatype type , int key , void *value, int * flag)
{
    /* Same as above */
    // *(void**)value = NULL;
    memset(value, 0, sizeof(void*));
    *flag=1;
    return 0;
}
static inline int MPI_Type_delete_attr(MPI_Datatype type , int key ) { return 0; }

typedef int MPI_Type_copy_attr_function(MPI_Datatype oldtype,
           int type_keyval, void *extra_state, void *attribute_val_in,
           void *attribute_val_out, int *flag);
typedef int MPI_Type_delete_attr_function(MPI_Datatype type, int type_keyval,
            void *attribute_val, void *extra_state);
#define MPI_TYPE_DUP_FN NULL
#define MPI_TYPE_NULL_DELETE_FN NULL
static inline int MPI_Type_create_keyval(
        MPI_Type_copy_attr_function *type_copy_attr_fn ,
        MPI_Type_delete_attr_function *type_delete_attr_fn ,
        int *type_keyval,
        void *extra_state )
{
    /* same rationale as above: we don't care about providing usable
     * function. Only the prototypes are barely right. The rest will
     * actually never be called.
     */
    *type_keyval = 0;
    return 0;
}
static inline int MPI_Type_free_keyval(int * x )
{
    return 0;
}

#define MPI_MAX_ERROR_STRING    64
static inline int MPI_Error_string(int err, char * msg, int * len)
{
    *len = snprintf(msg, MPI_MAX_ERROR_STRING, "%s", strerror(err));
    return *len >= 0 ? 0 : ENOMEM;
}

#ifdef  __GNUC__
#pragma GCC diagnostic pop
#endif

#endif /* FAKEMPI_H_ */
