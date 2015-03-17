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

thread_pool::thread_pool(size_t _nr_threads)
  : nr_threads(_nr_threads), tasks_not_empty_cond(PTHREAD_COND_INITIALIZER),
    results_not_empty_cond(PTHREAD_COND_INITIALIZER), kill_threads(false)
{
  threads = new worker_thread_ptr[nr_threads];
  for (size_t i = 0; i < nr_threads; i++)
    threads[i] = new worker_thread(*this);
};

thread_pool::~thread_pool() {
  enter();
  kill_threads = true;
  /* Wakey wakey, time to die */
  broadcast(tasks_not_empty_cond);
  leave();
  for (size_t i = 0; i < nr_threads; i++)
    delete threads[i];
  delete[] threads;
  ASSERT_ALWAYS(tasks.empty());
  ASSERT_ALWAYS(results.empty());
  pthread_cond_destroy(&tasks_not_empty_cond);
  pthread_cond_destroy(&results_not_empty_cond);
}

void *
thread_pool::thread_work_on_tasks(void *arg)
{
  worker_thread *I = (worker_thread *) arg;
  thread_pool &pool = I->pool;
  while (1) {
    thread_task *task = pool.get_task();
    if (task->please_die) {
      delete task;
      break;
    }
    task_function_t func = task->func;
    task_parameters *params = task->parameters;
    task_result *result = func(params);
    pool.add_result(result);
    delete task;
  }
  return NULL;
}

void
thread_pool::add_task(task_function_t func, task_parameters *const params,
                      const int id)
{
    enter();
    tasks.push(new thread_task(func, id, params));
    signal(tasks_not_empty_cond);
    leave();
}
  
thread_task *
thread_pool::get_task() {
  enter();
  while (!kill_threads && tasks.empty()) {
    /* No work -> tell this thread to wait until work becomes available.
       We also leave the loop when the thread needs to die.
       The while() protects against spurious wake-ups that can fire even if
       the queue is still empty. */
    wait(tasks_not_empty_cond);
  }
  thread_task *task;
  if (kill_threads && tasks.empty()) {
    task = new thread_task(true);
  } else {
    task = tasks.front();
    tasks.pop();
  }
  leave();
  return task;
}

void
thread_pool::add_result(task_result *const result) {
  enter();
  results.push(result);
  signal(results_not_empty_cond);
  leave();
}

task_result *
thread_pool::get_result() {
  enter();
  while (results.empty())
    wait(results_not_empty_cond);
  task_result *result = results.front();
  results.pop();
  leave();
  return result;
}
