#ifndef TESTS_COMMON_H
#define TESTS_COMMON_H

#include <stdint.h>

#define PARSE_SEED 1

int cmp_double(double, double, double);
int64_t random_int64 ();
uint64_t random_uint64 ();
void tests_common_cmdline(int *, const char ***, uint64_t);

#endif
