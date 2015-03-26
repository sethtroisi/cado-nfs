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
  for(int s = 0 ; s < 2 ; s++) {
    sides[s]->set_fb(si->sides[s]->fb->get_part(1));
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

/* {{{ thread_buckets_alloc
 * TODO: Allow allocating larger buckets if bucket_fill_ratio ever grows
 * above the initial guess. This can easily be made a permanent choice.
 *
 * Note also that we could consider having bucket_fill_ratio global.
 */
void
thread_workspaces::buckets_alloc()
{
  for (size_t i = 0; i < nr_workspaces; ++i) {
    for(unsigned int side = 0 ; side < 2 ; side++) {
      thread_side_data_ptr ts = thrs[i].sides[side];
      sieve_info_srcptr si = thrs[i].si;
      /* We used to re-allocate whenever the number of buckets changed. Now we
         always allocate memory for the max. number of buckets, so that we
         never have to re-allocate */
      uint32_t nb_buckets;
      /* If shell environment variable LAS_REALLOC_BUCKETS is *not* set,
         always allocate memory for the max. number of buckets, so that we
         never have to re-allocate. If it is set, allocate just enough for
         for the current number of buckets, which will re-allocate memory
         if number of buckets changes. */
      if (getenv("LAS_REALLOC_BUCKETS") == NULL) {
        nb_buckets = si->nb_buckets_max;
      } else {
        nb_buckets = si->nb_buckets;
      }
      init_buckets(&(ts->BA), &(ts->kBA), &(ts->mBA),
                   si->sides[side]->max_bucket_fill_ratio, nb_buckets);
    }
  }
}

void
thread_workspaces::buckets_free()
{
  for (size_t i = 0; i < nr_workspaces; ++i) {
    for(unsigned int side = 0 ; side < 2 ; side++) {
      // fprintf ("# Freeing buckets, thread->id=%d, side=%d\n", thrs[i].id, side);
      thread_side_data_ptr ts = thrs[i].sides[side];
      /* if there is no special-q in the interval, the arrays are not malloced */
      if (ts->BA.bucket_write != NULL)
        clear_buckets(&(ts->BA), &(ts->kBA), &(ts->mBA));
    }
  }
}

double
thread_workspaces::buckets_max_full()
{
    double mf0 = 0;
    for (size_t i = 0; i < nr_workspaces; ++i) {
      for(unsigned int side = 0 ; side < 2 ; side++) {
        double mf = ::buckets_max_full (thrs[i].sides[side]->BA);
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
