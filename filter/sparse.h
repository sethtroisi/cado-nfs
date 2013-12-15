#ifndef CADO_SPARSE_H_
#define CADO_SPARSE_H_

#include <stdint.h> /* for int32_t */

#ifndef FOR_DL
#define rowLength(rows, i) rows[(i)][0]
#define rowCell(rows, i, k) rows[(i)][(k)]
#else
#define rowLength(rows, i) rows[(i)][0].id
#define rowCell(rows, i, k) rows[(i)][(k)].id
#endif

// Structures for the data to create the index file.
// This is an array of relation-sets.
// One relation-set is an array of pairs (relation,multiplicity).
//
// Note: the reason why it's here is because the addRows function is
// shared by merge and replay, and needs to know about it.

typedef struct {
    uint32_t ind_row;
    int32_t e;
} multirel_t;

typedef struct {
    unsigned int n;
    multirel_t * rels;
} relset_t;

typedef relset_t * index_data_t;


// Protos.

extern void fprintRow(FILE *file, typerow_t *row);
extern int hasCol(int32_t **rows, int i, int32_t j);

extern int parse_hisfile_line (int32_t *ind, char *t, int32_t *j);

extern void addRowsUpdateIndex(typerow_t **rows, index_data_t index_data_t,
        int i1, int i2, int32_t j);

// The following version without index updating, is for merge.
static inline void addRows(typerow_t **rows, int i1, int i2, int32_t j) {
    addRowsUpdateIndex(rows, NULL, i1, i2, j);
}

#endif  /* CADO_SPARSE_H_ */
