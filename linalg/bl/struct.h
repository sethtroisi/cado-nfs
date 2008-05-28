#ifndef BL_STRUCT_H_
#define BL_STRUCT_H_

#ifdef __cplusplus
extern "C" { 
#endif
    
struct matrix_slice_s {
    unsigned long i0;
    unsigned long i1;
    unsigned long nbcoeffs;
    long start;
    long end;
};

typedef struct matrix_slice_s matrix_slice[1];
typedef struct matrix_slice_s * matrix_slice_ptr;

struct SparseMatrix_s {
    /* This is the _global_ info */
    unsigned long Nrows;
    unsigned long Ncols;
    // unsigned long Weight;
    /* This gives the info on all the matrix slices */
    matrix_slice * slices;
    unsigned long *Data;
};


typedef struct SparseMatrix_s SparseMatrix[1];


struct DenseMatrix_s {
    unsigned long Nrows;
    unsigned long Ncols;
    unsigned long *Data;
};

typedef struct DenseMatrix_s DenseMatrix[1];




#ifdef __cplusplus
}
#endif

#endif /* BL_STRUCT_H_ */
