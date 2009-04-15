#ifndef READMAT_EASY_H_
#define READMAT_EASY_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* The p_nr and p_nc pointers correspond to the dimensions of the
 * possibly transposed matrix ; IOW, this matches the organization of the
 * data returned.
 */
uint32_t * read_easy(const char * filename,
        unsigned int * p_nr, unsigned int * p_nc);

uint32_t * read_easy_transposed(const char * filename,
        unsigned int * p_nr, unsigned int * p_nc);

/* The return value of these functions is a pointer which must be freed
 */

#ifdef __cplusplus
}
#endif

#endif	/* READMAT_EASY_H_ */
