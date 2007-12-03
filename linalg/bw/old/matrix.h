#ifndef MATRIX_H_
#define MATRIX_H_

#ifdef	__cplusplus
extern "C" {
#endif

/*****************************************************************************/

typedef struct {
	type32  * indexes;
	stype32 * values;
} bw_matrix_t;

bw_matrix_t matrix;

extern coord_t nrows;
extern coord_t ncols;

extern size_t indexes_size;
extern size_t values_size;

/*****************************************************************************/

void load_matrix(void);
void multiply(bw_vector_block, bw_vector_block);
int compute_thread_offsets(void);
void reorder_indexes();

/*****************************************************************************/

#ifdef	__cplusplus
}
#endif

/* vim:set sw=8: */
#endif	/* MATRIX_H_ */
