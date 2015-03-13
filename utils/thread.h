#include <pthread.h>
#include <queue>
#include <list>
#include <stdarg.h>

typedef void *(*task_function_t)(void *);

class thread_task {
  task_function_t func;
  void *parameters;
};

class thread {
  pthread_t thread;
  size_t preferred_queue;
  thread(size_t _preferred_queue) : preferred_queue(_preferred_queue){
  }
};

class thread_pool {
  typedef std::queue<thread_task *> task_queue_t;
  typedef std::list<thread *> threads_t;

  task_queue_t *task_queues;
  threads_t busy_threads, waiting_threads;
  size_t nr_queues;
  size_t nr_threads;
  pthread_mutex_t mutex;
  pthread_cond_t cond;

  void enter() {pthread_mutex_lock(mutex);}
  void leave() {pthread_mutex_unlock(mutex);}

  void thread_assign_task();

public:
  thread_pool(const size_t _nr_queues);
  ~thread_pool();
  void add_threads(const size_t nr_threads, size_t preferred_queue);
  void add_task(size_t which_queue, task_function_t, void *);
};
