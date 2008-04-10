#ifndef BL_STRUCT_H_
#define BL_STRUCT_H_

#ifdef __cplusplus
extern "C" { 
#endif
    
struct SparseMatrix_s {
    unsigned long Nrows;
    unsigned long Ncols;
    unsigned long Weight;
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
