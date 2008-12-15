#define _POSIX_C_SOURCE 200112L
/* ugly but useful for c++ compilers to recognize inttypes.h */
#define __STDC_FORMAT_MACROS

#include "parallelizing_info.h"
#include <stdio.h>
#include <inttypes.h>
#include "xvectors.h"
#include "utils.h"

void setup_x_random(uint32_t * xs,
        unsigned int m, unsigned int nx, unsigned int nr,
        parallelizing_info_ptr pi)
{
    /* Here, everybody has to agree on an array of random values. The xs
     * pointer is on the stack of each calling thread, so threads must
     * converge to a commmon point of view on the data.
     */
    // job 0 thread 0 decides for everybody.
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        for(unsigned int i = 0 ; i < nx * m ; i++) {
            xs[i] = myrand() % nr;
        }
    }
    complete_broadcast(pi->m, xs, nx * m * sizeof(unsigned int), 0, 0);
}

void load_x(uint32_t * xs, unsigned int m, unsigned int nx,
        parallelizing_info_ptr pi)
{
    /* pretty much the same deal as above */
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        FILE * f = fopen("X", "r");
        FATAL_ERROR_CHECK(f == NULL, "Cannot open file X for reading");
        for (unsigned int i = 0 ; i < nx * m; i++) {
            int rc = fscanf(f, "%" SCNu32, &(xs[i]));
            FATAL_ERROR_CHECK(rc != 1, "short read in file X");
        }
        fclose(f);
    }
    complete_broadcast(pi->m, xs, nx * m * sizeof(unsigned int), 0, 0);
}

void save_x(uint32_t * xs, unsigned int m, unsigned int nx, parallelizing_info_ptr pi)
{
    /* Here, we expect that the data is already available to eveybody, so
     * no synchronization is necessary.
     */
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        FILE * f = fopen("X", "w");
        FATAL_ERROR_CHECK(f == NULL, "Cannot open file X for writing");
        for (unsigned int i = 0 ; i < nx * m; i++) {
            int rc = fprintf(f, "%" PRIu32, xs[i]);
            FATAL_ERROR_CHECK(rc != 1, "short write in file X");
        }
        fclose(f);
    }
}
