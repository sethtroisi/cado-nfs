#ifndef LAS_THREADS_H_
#define LAS_THREADS_H_

#include "las-forwardtypes.h"
#include "bucket.h"
#include "fb.h"
#include "las-report-stats.h"
#include "las-base.hpp"

/* {{{ thread-related defines */
/* All of this exists _for each thread_ */
struct thread_side_data_s : private NonCopyable {
  m_bucket_array_t mBA; /* Not used if not fill_in_m_buckets (3 passes sort) */
  k_bucket_array_t kBA; /* Ditto for fill_in_k_buckets (2 passes sort) */
  bucket_array_t BA;    /* Always used */
  const fb_part *fb;

  /* For small sieve */
  int * ssdpos;
  int * rsdpos;

  unsigned char *bucket_region;
  sieve_checksum checksum_post_sieve;

  thread_side_data_s();
  ~thread_side_data_s();
  void set_fb(const fb_part *_fb) {fb = _fb;}
  void update_checksum(){checksum_post_sieve.update(bucket_region, BUCKET_REGION);}
};
typedef struct thread_side_data_s thread_side_data[1];
typedef struct thread_side_data_s * thread_side_data_ptr;
typedef const struct thread_side_data_s * thread_side_data_srcptr;

struct thread_data : private NonCopyable {
  int id;
  thread_side_data sides[2];
  las_info_srcptr las;
  sieve_info_ptr si;
  las_report rep;
  unsigned char *SS;
  bool is_initialized;
  thread_data();
  ~thread_data();
  void init(int id, las_info_srcptr las);
  void pickup_si(sieve_info_ptr si);
  void update_checksums();
};

class thread_workspaces {
  thread_data *thrs;
  const size_t nr_workspaces;
public:
  thread_workspaces(size_t n, las_info_ptr _las);
  ~thread_workspaces();
  void pickup_si(sieve_info_ptr si);
  void thread_do(void * (*) (thread_data *));
  void buckets_alloc();
  void buckets_free();
  double buckets_max_full();
  void accumulate(las_report_ptr, sieve_checksum *);
};

#endif
