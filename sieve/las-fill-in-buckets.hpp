#ifndef LAS_FILL_IN_BUCKETS_HPP_
#define LAS_FILL_IN_BUCKETS_HPP_

#include "las-types.hpp"
#include "las-threads.hpp"
#include "tdict.hpp"

typedef std::vector<plattices_vector_t *> precomp_plattice_t [2][FB_MAX_PARTS];

template <int LEVEL>
void
downsort_tree(timetree_t&,
        uint32_t bucket_index,
        uint32_t first_region0_index,
        thread_workspaces &ws,
        thread_pool &pool,
        sieve_info& si,
        precomp_plattice_t precomp_plattice);
void fill_in_buckets_both(timetree_t&, thread_pool &, thread_workspaces &, sieve_info const &);

#endif
