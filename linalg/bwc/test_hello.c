#include "cado.h"
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "select_mpi.h"
#include "params.h"
#include "portability.h"
#include "macros.h"

int m,n;
int verbose=0;

void * program(parallelizing_info_ptr pi, param_list pl MAYBE_UNUSED, void * arg MAYBE_UNUSED)
{
    // it is here as a cheap sanity check.
    hello(pi);

    return NULL;
}

void usage()
{
    fprintf(stderr, "Usage: test-hello <options>\n"
        "\tmpi=<int>x<int>\tset number of mpi jobs. Must agree with mpirun\n"
        "\tthr=<int>x<int>\tset number of threads.\n"
        );
    exit(1);
}

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    param_list pl;
    pi_outer_info poi;

    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    param_list_init (pl);
    pi_outer_info_init(poi);
    argv++, argc--;
    param_list_configure_switch(pl, "-v", &verbose);
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage();
    }

    pi_interpret_parameters(poi, pl);
    param_list_lookup_string(pl, "cpubinding");

    if (verbose)
        param_list_display (pl, stderr);
    if (param_list_warn_unused(pl)) {
        usage();
    }

    pi_go(poi, program, pl, 0);

    MPI_Finalize();

    pi_outer_info_clear(poi);
    param_list_clear(pl);

    return 0;
}

