#include "cado.h"
#include "thread.h"

worker_thread::worker_thread(thread_pool &_pool) : pool(_pool)
{
  int rc = pthread_create(&thread, NULL, pool.thread_work_on_tasks, this);
  ASSERT_ALWAYS(rc == 0);
}

worker_thread::~worker_thread() {
  int rc = pthread_join(thread, NULL);
  ASSERT_ALWAYS(rc == 0);
}

void *
thread_pool::thread_work_on_tasks(void *arg)
{
  worker_thread *I = (worker_thread *) arg;
  thread_pool &pool = I->pool;
  while (1) {
    thread_task *task = pool.get_task();
    if (task->please_die)
      break;
    task_function_t func = task->func;
    task_parameters *params = task->parameters;
    task_result *result = func(params);
    pool.add_result(result);
    delete task;
  }
  return NULL;
}
