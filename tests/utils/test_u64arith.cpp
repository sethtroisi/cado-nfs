#include "cado.h"
#include <cstdint>
#include <cinttypes>
#include "macros.h"
#include "tests_common.h"
#include "test_iter.h"
#include "u64arith.h"

void
test_one_u64arith_gt_2_2(const uint64_t a1, const uint64_t a2,
  const uint64_t b1, const uint64_t b2, const int v) {
  const int r = u64arith_gt_2_2(a1, a2, b1, b2);
  if (v != r) {
    printf("%s: Error, got result %d, but expected %d\n", __func__, v, r);
    abort();
  }
}

void test_u64arith_gt_2_2() {
  test_one_u64arith_gt_2_2(0,0,0,0,0);
  test_one_u64arith_gt_2_2(1,0,0,0,1);
  test_one_u64arith_gt_2_2(0,1,0,0,1);
  test_one_u64arith_gt_2_2(1,1,0,0,1);
  test_one_u64arith_gt_2_2(0,0,1,0,0);
  test_one_u64arith_gt_2_2(1,0,1,0,0);
  test_one_u64arith_gt_2_2(0,1,1,0,1);
  test_one_u64arith_gt_2_2(1,1,1,0,1);
  test_one_u64arith_gt_2_2(0,0,0,1,0);
  test_one_u64arith_gt_2_2(1,0,0,1,0);
  test_one_u64arith_gt_2_2(0,1,0,1,0);
  test_one_u64arith_gt_2_2(1,1,0,1,1);
  test_one_u64arith_gt_2_2(0,0,1,1,0);
  test_one_u64arith_gt_2_2(1,0,1,1,0);
  test_one_u64arith_gt_2_2(0,1,1,1,0);
  test_one_u64arith_gt_2_2(1,1,1,1,0);
}

void
test_one_u64arith_add_1_2(uint64_t r1, uint64_t r2,
    const uint64_t a, const uint64_t v1, const uint64_t v2) {
  u64arith_add_1_2(&r1, &r2, a);
  if (r1 != v1 || r2 != v2) {
    printf("%s: Error, got result %"PRIu64":%"PRIu64", but expected %"PRIu64":%"PRIu64"\n",
        __func__, r2, r1, v2, v1);
    exit(EXIT_FAILURE);
  }
}

void
test_u64arith_add_1_2() {
  test_one_u64arith_add_1_2(1, 4, 1, 2, 4);
  test_one_u64arith_add_1_2(UINT64_MAX, 1, 1, 0, 2);
  test_one_u64arith_add_1_2(UINT64_MAX, UINT64_MAX, 1, 0, 0);
}


void
test_one_u64arith_add_2_2(uint64_t r1, uint64_t r2,
    const uint64_t a1, const uint64_t a2, const uint64_t v1, const uint64_t v2) {
  u64arith_add_2_2(&r1, &r2, a1, a2);
  if (r1 != v1 || r2 != v2) {
    printf("%s: Error, got result %"PRIu64":%"PRIu64", but expected %"PRIu64":%"PRIu64"\n",
        __func__, r2, r1, v2, v1);
    exit(EXIT_FAILURE);
  }
}


void
test_u64arith_add_2_2() {
  test_one_u64arith_add_2_2(1, 4, 1, 0, 2, 4);
  test_one_u64arith_add_2_2(UINT64_MAX, 1, 1, 0, 0, 2);
  test_one_u64arith_add_2_2(UINT64_MAX, UINT64_MAX, 1, 0, 0, 0);
  test_one_u64arith_add_2_2(1, 4, 0, 1, 1, 5);
}


void
test_one_u64arith_add_2_2_cy(uint64_t r1, uint64_t r2,
    const uint64_t a1, const uint64_t a2, const uint64_t v1, const uint64_t v2,
    const char vcy) {
  unsigned char cy = u64arith_add_2_2_cy(&r1, &r2, a1, a2);
  if (r1 != v1 || r2 != v2 || vcy != cy) {
    printf("%s: Error, got result %hhu:%"PRIu64":%"PRIu64", but expected %hhu:%"PRIu64":%"PRIu64"\n",
        __func__, cy, r2, r1, vcy, v2, v1);
    exit(EXIT_FAILURE);
  }
}


void
test_u64arith_add_2_2_cy() {
  test_one_u64arith_add_2_2_cy(1, 4, 1, 0, 2, 4, 0);
  test_one_u64arith_add_2_2_cy(UINT64_MAX, 1, 1, 0, 0, 2, 0);
  test_one_u64arith_add_2_2_cy(UINT64_MAX, UINT64_MAX, 1, 0, 0, 0, 1);
  test_one_u64arith_add_2_2_cy(1, 4, 0, 1, 1, 5, 0);
}


void
test_one_u64arith_addmod_1_1(uint64_t *r, const uint64_t a,
                     const uint64_t b, const uint64_t m) {
  uint64_t r;
  u64arith_addmod_1_1 (&r, a, b, m);

}


void
test_one_u64arith_mul_1_1_2(const uint64_t a, const uint64_t b,
    const uint64_t v1, const uint64_t v2) {
  uint64_t r1, r2;
  u64arith_mul_1_1_2 (&r1, &r2, a, b);
  if (r1 != v1 || r2 != v2) {
    printf("%s: Error, got result %"PRIu64":%"PRIu64", but expected %"PRIu64":%"PRIu64"\n",
        __func__, r2, r1, v2, v1);
    exit(EXIT_FAILURE);
  }
}

void
test_u64arith_mul_1_1_2() {
  test_one_u64arith_mul_1_1_2(1, 2, 2, 0);
  uint64_t a = 1, b = 1;
  a <<= 31; b <<= 33;
  test_one_u64arith_mul_1_1_2(a, b, 0, 1);
  test_one_u64arith_mul_1_1_2(UINT64_MAX, UINT64_MAX - 1, 2, UINT64_MAX - 2);
  test_one_u64arith_mul_1_1_2(UINT64_MAX, UINT64_MAX, 1, UINT64_MAX - 1);
}


void
test_one_u64arith_sqr_1_2(const uint64_t a, const uint64_t v1, const uint64_t v2) {
  uint64_t r1, r2;
  u64arith_sqr_1_2 (&r1, &r2, a);
  if (r1 != v1 || r2 != v2) {
    printf("%s: Error, got result %"PRIu64":%"PRIu64", but expected %"PRIu64":%"PRIu64"\n",
        __func__, r2, r1, v2, v1);
    exit(EXIT_FAILURE);
  }
}


void
test_u64arith_sqr_1_2() {
  test_one_u64arith_sqr_1_2(2, 4, 0);
  uint64_t a = 1;
  a <<= 32;
  test_one_u64arith_sqr_1_2(a, 0, 1);
  test_one_u64arith_sqr_1_2(UINT64_MAX, 1, UINT64_MAX - 1);
}


int
main (int argc, const char *argv[])
{
  tests_common_cmdline(&argc, &argv, 0);
  test_u64arith_gt_2_2();
  test_u64arith_add_1_2();
  test_u64arith_add_2_2();
  test_u64arith_add_2_2_cy();
  test_u64arith_mul_1_1_2();
  test_u64arith_sqr_1_2();
  tests_common_clear();
  exit (EXIT_SUCCESS);
}
