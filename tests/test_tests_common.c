#include "cado.h"
#include "tests_common.h"

int main(int argc, const char **argv)
{
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
}
