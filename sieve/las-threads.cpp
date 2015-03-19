#include "cado.h"
#include <pthread.h>
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


thread_data *thread_data_alloc(las_info_ptr las, int n)
{
    thread_data *thrs = new thread_data[n];
    ASSERT_ALWAYS(thrs != NULL);

    for(int i = 0 ; i < n ; i++)
        thrs[i].init(i, las);
    return thrs;
}

void thread_data_free(thread_data *thrs)
{
    delete[] thrs;
}

void thread_pickup_si(thread_data *thrs, sieve_info_ptr si, int n)
{
    for (int i = 0; i < n ; ++i) {
        thrs[i].pickup_si(si);
    }
}

void thread_do(thread_data *thrs, void * (*f) (thread_data *), int n)/*{{{*/
{
    if (n == 1) {
        /* Then don't bother with pthread calls */
        (*f)(&thrs[0]);
        return;
    }
    pthread_t * th = (pthread_t *) malloc(n * sizeof(pthread_t)); 
    ASSERT_ALWAYS(th);

#if 0
    /* As a debug measure, it's possible to activate this branch instead
     * of the latter. In effect, this causes las to run in a
     * non-multithreaded way, albeit strictly following the code path of
     * the multithreaded case.
     */
    for (int i = 0; i < n ; ++i) {
        (*f)(&thrs[i]);
    }
#else
    for (int i = 0; i < n ; ++i) {
        int ret = pthread_create(&(th[i]), NULL, 
		(void * (*)(void *)) f,
                (void *)(&thrs[i]));
        ASSERT_ALWAYS(ret == 0);
    }
    for (int i = 0; i < n ; ++i) {
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
void thread_buckets_alloc(thread_data *thrs, unsigned int n)/*{{{*/
{
  for (unsigned int i = 0; i < n ; ++i) {
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
}/*}}}*/

void thread_buckets_free(thread_data *thrs, unsigned int n)/*{{{*/
{
  for (unsigned int i = 0; i < n ; ++i) {
    for(unsigned int side = 0 ; side < 2 ; side++) {
      // fprintf ("# Freeing buckets, thread->id=%d, side=%d\n", thrs[i].id, side);
      thread_side_data_ptr ts = thrs[i].sides[side];
      /* if there is no special-q in the interval, the arrays are not malloced */
      if (ts->BA.bucket_write != NULL)
        clear_buckets(&(ts->BA), &(ts->kBA), &(ts->mBA));
    }
  }
}/*}}}*/

double thread_buckets_max_full(thread_data *thrs, int n)/*{{{*/
{
    double mf0 = 0;
    for (int i = 0; i < n ; ++i) {
        double mf = buckets_max_full (thrs[i].sides[RATIONAL_SIDE]->BA);
        if (mf > mf0) mf0 = mf;
        mf = buckets_max_full (thrs[i].sides[ALGEBRAIC_SIDE]->BA);
        if (mf > mf0) mf0 = mf;
    }
    return mf0;
}/*}}}*/

