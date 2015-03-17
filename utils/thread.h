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
class task_parameters;

/* Base for classes that hold results produced by worker functions */
class task_result;

class monitor {
  pthread_mutex_t mutex;
public:
  monitor() : mutex(PTHREAD_MUTEX_INITIALIZER) {}
  ~monitor(){ASSERT_ALWAYS(pthread_mutex_destroy(&mutex) == 0);}
  void enter() {ASSERT_ALWAYS(pthread_mutex_lock(&mutex) == 0);}
  void leave() {ASSERT_ALWAYS(pthread_mutex_unlock(&mutex) == 0);}
  void signal(pthread_cond_t &cond) {ASSERT_ALWAYS(pthread_cond_signal(&cond) == 0);}
  void wait(pthread_cond_t &cond) {ASSERT_ALWAYS(pthread_cond_wait(&cond, &mutex) == 0);}
};

typedef task_result *(*task_function_t)(const task_parameters *);

class thread_task {
public:
  const task_function_t func;
  const int id;
  task_parameters * const parameters;
  task_result *result;
  const bool please_die;

  thread_task(task_function_t _func, int _id, task_parameters *_parameters) :
    func(_func), id(_id), parameters(_parameters), result(NULL), please_die(false) {};
  thread_task(bool _kill)
    : func(NULL), id(0), parameters(NULL), result(NULL), please_die(true) {
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
  thread_pool(size_t _nr_threads) : nr_threads(_nr_threads), kill_threads(false) {
    threads = new worker_thread_ptr[nr_threads];
    for (size_t i = 0; i < nr_threads; i++)
      threads[i] = new worker_thread(*this);
  };
  ~thread_pool() {
    enter();
    kill_threads = true;
    leave();
    for (size_t i = 0; i < nr_threads; i++)
      delete threads[i];
    delete threads;
    ASSERT_ALWAYS(tasks.empty());
    ASSERT_ALWAYS(results.empty());
    pthread_cond_destroy(&tasks_not_empty_cond);
    pthread_cond_destroy(&results_not_empty_cond);
  }

  void add_task(task_function_t func, task_parameters *params, int id) {
    enter();
    tasks.push(new thread_task(func, id, params));
    signal(tasks_not_empty_cond);
    leave();
  }
  
  thread_task *get_task() {
    enter();
    thread_task *task;
    if (kill_threads && tasks.empty()) {
      task = new thread_task(true);
    } else {
      while (tasks.empty()) {
        /* No work -> tell this thread to wait until work becomes available.
           The while() protects against spurious wake-ups that can fire even
           if the queue is still empty */
        wait(tasks_not_empty_cond);
      }
      task = tasks.front();
      tasks.pop();
    }
    leave();
    return task;
  }
  
  void add_result(task_result *const result) {
    enter();
    results.push(result);
    signal(results_not_empty_cond);
    leave();
  }

  task_result *get_result() {
    enter();
    while (results.empty())
      wait(results_not_empty_cond);
    task_result *result = results.front();
    results.pop();
    leave();
    return result;
  }
};
