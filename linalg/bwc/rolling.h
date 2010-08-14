#ifndef ROLLING_H_
#define ROLLING_H_

#include "balancing.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Removes all vector files in bwc's wdir (current directory, in fact,
 * since bwc programs always run in the current directory) according to
 * the rules specified by the parameters keep_rolling_checkpoints and
 * checkpoint_precious.
 *
 * If v == 0, consider all checkpoints. Otherwise, restrict to those
 * whose index is <= v
 *
 * The stem is to be understood according to the
 * COMMON_VECTOR_ITERATE_PATTERN constant in filenames.h
 *
 * The bal argument is used only to compose the filename according to the
 * checksum of the current balancing permutation.
 */
extern void keep_rolling_checkpoints(balancing_ptr bal, const char * stem, unsigned int v);

#ifdef __cplusplus
}
#endif

#endif	/* ROLLING_H_ */
