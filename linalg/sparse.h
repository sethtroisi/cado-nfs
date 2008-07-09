#ifndef CADO_SPARSE_H_
#define CADO_SPARSE_H_
/* Data type to store column values: a 32-bit integer should be enough for
   most applications!
*/
#include <stdint.h>

/* TODO: rename this to sparse_coeff_t when Francois is done with the mpi
 * merge stuff.  */

// please, do not use an unsigned type!
// XXX [ET] Why ? 
// [FM] because of some trick on negative numbers in MST related stuff!
#define INT int32_t

extern void fprintRow(FILE *file, INT *row);
extern INT * copyRow(INT *row);
extern void removeWeight(INT **rows, int *wt, int i);
extern void addWeight(INT **rows, int *wt, int i);
extern void addRows(INT **rows, int i1, int i2, int len0);

#endif  /* CADO_SPARSE_H_ */
