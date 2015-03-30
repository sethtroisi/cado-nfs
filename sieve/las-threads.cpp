#include "cado.h"
#include "memory.h"
#include "las-threads.h"
#include "las-types.h"
#include "las-config.h"

thread_side_data_s::thread_side_data_s()
{
  /* Allocate memory for each side's bucket region */
  bucket_region = (unsigned char *) contiguous_malloc(BUCKET_REGION + MEMSET_MIN);
}

thread_side_data_s::~thread_side_data_s()
{
  contiguous_free(bucket_region);
  bucket_region = NULL;
}

void
thread_side_data_s::allocate_bucket_array(const uint32_t n_bucket, const double fill_ratio)
{
  const size_t bucket_region = (size_t) 1 << LOG_BUCKET_REGION;
  const size_t bucket_size = bucket_region * fill_ratio;
  BA.allocate_memory(n_bucket, bucket_size);
}

thread_data::thread_data() : is_initialized(false)
{
  /* Allocate memory for the intermediate sum (only one for both sides) */
  SS = (unsigned char *) contiguous_malloc(BUCKET_REGION);
  las_report_init(rep);
}

thread_data::~thread_data()
{
  ASSERT_ALWAYS(is_initialized);
  ASSERT_ALWAYS(SS != NULL);
  contiguous_free(SS);
  SS = NULL;
  las_report_clear(rep);
}

void thread_data::init(const int _id, las_info_srcptr _las)
{
  ASSERT_ALWAYS(!is_initialized);
  id = _id;
  las = _las;
  is_initialized = true;
}

void thread_data::pickup_si(sieve_info_ptr _si)
{
  si = _si;
  for (int side = 0 ; side < 2 ; side++) {
    sides[side]->set_fb(si->sides[side]->fb->get_part(1));
    /* Always allocate the max number of buckets (i.e., as if we were using the
       max value for J), even if we use a smaller J due to a poor q-lattice
       basis */
    sides[side]->allocate_bucket_array(si->nb_buckets_max, si->sides[side]->max_bucket_fill_ratio);
  }
}

void thread_data::update_checksums()
{
  for(int s = 0 ; s < 2 ; s++)
    sides[s]->update_checksum();
}

thread_workspaces::thread_workspaces(const size_t n, las_info_ptr las)
  : nr_workspaces(n)
{
    pthread_mutex_init(&mutex, NULL);
    thrs = new thread_data[n];
    ASSERT_ALWAYS(thrs != NULL);
    used = new bool[n];
    ASSERT_ALWAYS(used != NULL);

    for(size_t i = 0 ; i < nr_workspaces; i++) {
        thrs[i].init(i, las);
        used[i] = false;
    }
}

thread_workspaces::~thread_workspaces()
{
    delete[] thrs;
    delete[] used;
    pthread_mutex_destroy(&mutex);
}

/* Prepare to work on sieving a special-q as described by _si.
   This implies allocating all the memory we need for bucket arrays,
   sieve regions, etc. */
void
thread_workspaces::pickup_si(sieve_info_ptr _si)
{
    si = _si;
    for (size_t i = 0; i < nr_workspaces; ++i) {
        thrs[i].pickup_si(_si);
    }
}

void
thread_workspaces::thread_do(void * (*f) (thread_data *))
{
    if (nr_workspaces == 1) {
        /* Then don't bother with pthread calls */
        (*f)(&thrs[0]);
        return;
    }
    pthread_t * th = (pthread_t *) malloc(nr_workspaces * sizeof(pthread_t)); 
    ASSERT_ALWAYS(th);

#if 0
    /* As a debug measure, it's possible to activate this branch instead
     * of the latter. In effect, this causes las to run in a
     * non-multithreaded way, albeit strictly following the code path of
     * the multithreaded case.
     */
    for (size_t i = 0; i < nr_workspaces; ++i) {
        (*f)(&thrs[i]);
    }
#else
    for (size_t i = 0; i < nr_workspaces; ++i) {
        int ret = pthread_create(&(th[i]), NULL, 
		(void * (*)(void *)) f,
                (void *)(&thrs[i]));
        ASSERT_ALWAYS(ret == 0);
    }
    for (size_t i = 0; i < nr_workspaces; ++i) {
        int ret = pthread_join(th[i], NULL);
        ASSERT_ALWAYS(ret == 0);
    }
#endif

    free(th);
}

double
thread_workspaces::buckets_max_full()
{
    double mf0 = 0;
    for (size_t i = 0; i < nr_workspaces; ++i) {
      for(unsigned int side = 0 ; side < 2 ; side++) {
        double mf = thrs[i].sides[side]->BA.max_full();
        if (mf > mf0) mf0 = mf;
      }
    }
    return mf0;
}

void
thread_workspaces::accumulate(las_report_ptr rep, sieve_checksum *checksum)
{
    for (size_t i = 0; i < nr_workspaces; ++i) {
        las_report_accumulate(rep, thrs[i].rep);
        for (int side = 0; side < 2; side++)
            checksum[side].update(thrs[i].sides[side]->checksum_post_sieve);
    }
}

thread_data &
thread_workspaces::reserve_workspace()
{
  pthread_mutex_lock(&mutex);
  size_t i;
  for (i = 0; i < nr_workspaces; ++i) {
    if (!used[i])
      break;
  }
  /* Currently no logic to make threads wait until workspace becomes available */
  ASSERT_ALWAYS(i < nr_workspaces);
  used[i] = true;
  pthread_mutex_unlock(&mutex);
  return thrs[i];
}

void 
thread_workspaces::release_workspace(thread_data &ws)
{
  pthread_mutex_lock(&mutex);
  size_t i;
  for (i = 0; i < nr_workspaces; ++i) {
    if (&ws == &thrs[i])
      break;
  }
  ASSERT_ALWAYS(i < nr_workspaces);
  used[i] = false;
  pthread_mutex_unlock(&mutex);
}
