static void *thread_wait_for_task(void * pool)
{
  /* See if there is a task in my preferred queue. If not, see if there is a task in another queue */
  /* If there is a task, remove it from the list, remove the thread from the list of waiting threads, and work on the task */
  /* If there is no task, wait for the condition variable */
}

thread_pool::thread_pool(const size_t _nr_queues)
    : mutex(PTHREAD_MUTEX_INITIALIZER), cond(PTHREAD_COND_INITIALIZER), nr_queues(_nr_queues)
{
  queues = new task_queue_t[nr_queues];
  busy_threads = new thread_queue_t[nr_queues];
  waiting_threads = new thread_queue_t[nr_queues];
}

thread_pool::~thread_pool()
{
  /* Wait until all the tasks queue is empty and all threads are idle */
  /* Then join all threads, delete the queues */
}

void thread_pool::add_threads(const size_t nr_threads, const size_t preferred_queue)
{
  pthread_t *thread = new pthread_t;
  pthread_attr_t attr;
  int rc;
  rc = pthread_attr_init (&attr);
  rc = pthread_create(thread, attr, thread_wait_for_task, this);
  rc = pthread_attr_destroy(&attr);
  
}

void thread_pool::add_task(size_t which_queue, task_function_t, void *)
{
}
