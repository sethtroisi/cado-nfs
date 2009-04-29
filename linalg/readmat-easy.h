#ifndef READMAT_EASY_H_
#define READMAT_EASY_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Read the matrix in memory as a stupid straight array of 32-bit
 * unsigned integers. Store the matrix either in normal or transposed
 * order (or both) in the pointers whose address is p_direct and
 * p_transposed (either may be NULL, but obviously not both. Having both
 * non-null is valid). In any case, if *p_direct or *p_transposed is set
 * to something by this function, the resulting pointer must be freed.
 *
 * The number of rows and columns is stored in *p_nr and *p_nc
 */
void read_easy(const char * filename,
        uint32_t ** p_direct, uint32_t ** p_transposed,
        unsigned int * p_nr, unsigned int * p_nc);

#ifdef __cplusplus
}
#endif

#endif	/* READMAT_EASY_H_ */
