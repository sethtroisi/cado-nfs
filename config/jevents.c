#include "rdpmc.h"

int main()
{
    struct rdpmc_ctx ctx;
    if (rdpmc_open(PERF_COUNT_HW_CPU_CYCLES, &ctx) < 0)
          return 1;
    rdpmc_close (&ctx);
    return 0;
}
