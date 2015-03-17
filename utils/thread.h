#include <pthread.h>
#include <queue>
#include <vector>
#include <stdarg.h>

#include "macros.h"

class NonCopyable {
 protected:
   NonCopyable() {}
   ~NonCopyable() {}
 private:
   NonCopyable(const NonCopyable&);
   NonCopyable& operator=(const NonCopyable&);
};

/* Base for classes that hold parameters for worker functions */
class task_parameters {
};

/* Base for classes that hold results produced by worker functions */
class task_result {
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

class thread_task {
public:
  const task_function_t func;
  const int id;
  task_parameters * const parameters;
  const bool please_die;
  task_result *result;

  thread_task(task_function_t _func, int _id, task_parameters *_parameters) :
    func(_func), id(_id), parameters(_parameters), please_die(false), result(NULL) {};
  thread_task(bool _kill)
    : func(NULL), id(0), parameters(NULL), please_die(true), result(NULL) {
    ASSERT_ALWAYS(_kill);
  }
};

class thread_pool;

class worker_thread : private NonCopyable {
  friend class thread_pool;
  thread_pool &pool;
  pthread_t thread;
public:
  worker_thread(thread_pool &_pool);
  ~worker_thread();
};

class thread_pool : private monitor, private NonCopyable {
  friend class worker_thread;
  typedef worker_thread *worker_thread_ptr;
  worker_thread_ptr *threads;
  std::queue<thread_task *> tasks;
  std::queue<task_result *> results;

  const size_t nr_threads;

  pthread_cond_t tasks_not_empty_cond, results_not_empty_cond;

  bool kill_threads; /* If true, hands out kill tasks once queue is empty */

  static void * thread_work_on_tasks(void *pool);

public:
  thread_pool(size_t _nr_threads);
  ~thread_pool();
  void add_task(task_function_t func, task_parameters *params, int id);
  thread_task *get_task();
  void add_result(task_result *result);
  task_result *get_result();
};
