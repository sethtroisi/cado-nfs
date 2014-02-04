/* Iterators that produce the smallest and largest possible values uint64_t,
   and smallest, largest, and closest to 0 for int64_t. */

struct test_iter_s {
  unsigned long next, max;
};
typedef struct test_iter_s test_iter_t[1];
typedef struct test_iter_s * test_iter_ptr;
typedef const struct test_iter_s * test_iter_srcptr;

void test_iter_init (test_iter_t, unsigned long);
int test_iter_isdone (test_iter_t);
int64_t test_iter_int64_next(test_iter_t);
uint64_t test_iter_uint64_next(test_iter_t);
