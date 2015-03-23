#include <pthread.h>
#include <queue>
#include <vector>
#include <stdarg.h>

#include "macros.h"

class ThreadNonCopyable {
 protected:
   ThreadNonCopyable() {}
   ~ThreadNonCopyable() {}
 private:
   ThreadNonCopyable(const ThreadNonCopyable&);
   ThreadNonCopyable& operator=(const ThreadNonCopyable&);
};

/* Base for classes that hold parameters for worker functions */
class task_parameters {
  public:
  virtual ~task_parameters(){};
};

/* Base for classes that hold results produced by worker functions */
class task_result {
  public:
  virtual ~task_result(){};
};

class monitor {
  pthread_mutex_t mutex;
public:
  monitor() : mutex(PTHREAD_MUTEX_INITIALIZER) {}
  ~monitor(){ASSERT_ALWAYS(pthread_mutex_destroy(&mutex) == 0);}
  void enter() {ASSERT_ALWAYS(pthread_mutex_lock(&mutex) == 0);}
  void leave() {ASSERT_ALWAYS(pthread_mutex_unlock(&mutex) == 0);}
  void signal(pthread_cond_t &cond) {ASSERT_ALWAYS(pthread_cond_signal(&cond) == 0);}
  void broadcast(pthread_cond_t &cond){ASSERT_ALWAYS(pthread_cond_broadcast(&cond) == 0);}
  void wait(pthread_cond_t &cond) {ASSERT_ALWAYS(pthread_cond_wait(&cond, &mutex) == 0);}
};

typedef task_result *(*task_function_t)(const task_parameters *);

class thread_task;
class thread_pool;
class worker_thread;

class thread_pool : private monitor, private ThreadNonCopyable {
  friend class worker_thread;
  typedef worker_thread *worker_thread_ptr;
  worker_thread_ptr *threads;
  std::queue<thread_task *> tasks;
  std::queue<task_result *> results;

  const size_t nr_threads;

  pthread_cond_t tasks_not_empty_cond, results_not_empty_cond;

  bool kill_threads; /* If true, hands out kill tasks once work queue is empty */

  static void * thread_work_on_tasks(void *pool);
  thread_task *get_task();
  void add_result(task_result *result);

public:
  thread_pool(size_t _nr_threads);
  ~thread_pool();
  void add_task(task_function_t func, task_parameters *params, int id);
  task_result *get_result();
};
