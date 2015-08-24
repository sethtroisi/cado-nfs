#ifndef LAS_FILL_IN_BUCKETS_H_
#define LAS_FILL_IN_BUCKETS_H_

#include "las-types.h"
#include "las-threads.h"

typedef std::vector<plattices_vector_t *> precomp_plattice_t [2][FB_MAX_PARTS];

template <int LEVEL>
void
downsort_tree(uint32_t bucket_index,
        uint32_t first_region0_index,
        thread_workspaces &ws,
        thread_pool &pool,
        sieve_info_ptr si,
        precomp_plattice_t precomp_plattice);
void fill_in_buckets_both(thread_pool &, thread_workspaces &, sieve_info_srcptr);

#endif
