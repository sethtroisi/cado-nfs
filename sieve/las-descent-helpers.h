#ifndef LAS_DESCENT_HELPERS_H_
#define LAS_DESCENT_HELPERS_H_

#include "relation.h"
#include "las-types.h"

/* This interface is actually implemented by C++ code, but all calls
 * are plain C. */

#ifdef __cplusplus
extern "C" {
#endif

void * las_descent_helper_alloc();

// void las_descent_helper_start_forest(void *);
// void las_descent_helper_end_forest(void *);

void las_descent_helper_new_node(void *, siever_config_srcptr label, int level);
void las_descent_helper_found_relation(void *, relation_t *);
void las_descent_helper_done_node(void *);

void las_descent_helper_display_last_tree(void *, FILE *);
void las_descent_helper_display_all_trees(void *, FILE *);
void las_descent_helper_free(void *);

int las_descent_helper_current_depth(void *);

#ifdef __cplusplus
}
#endif

#endif	/* LAS_DESCENT_HELPERS_H_ */
