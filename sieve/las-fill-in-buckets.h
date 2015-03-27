#ifndef LAS_FILL_IN_BUCKETS_H_
#define LAS_FILL_IN_BUCKETS_H_

#include "las-types.h"
#include "las-threads.h"

void fill_in_buckets_both(thread_pool &, thread_workspaces &, int, sieve_info_srcptr);

#endif
