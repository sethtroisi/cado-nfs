#define _POSIX_C_SOURCE 200112L

#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>

#include "parallelizing_info.h"
#include "matmul_top.h"
#include "abase.h"
#include "select_mpi.h"
#include "random_generation.h"

#include "params.h"

abobj_t abase;

int m,n;
mpz_t p;
int interval=1000;
int verbose=0;
char matrix_filename[FILENAME_MAX];

void * program(parallelizing_info_ptr pi)
{
    // Doing the ``hello world'' test is a very good way of testing the
    // global mpi/pthreads setup. So despite its apparent irrelevance, I
    // suggest leaving it here as a cheap sanity check.
    hello(pi);
    matmul_top_data mmt;
    matmul_top_init(mmt, abase, pi, matrix_filename);

    return NULL;
}

void usage()
{
    fprintf(stderr, "Usage: bw-prep <options>\n"
        "Allowed options are\n"
        "\tm=<int>\t(*) set m blocking factor\n"
        "\tn=<int>\t(*) set n blocking factor\n"
        "\tmn=<int>\tset both m and n (exclusive with the two above)\n"
        "\tmpi=<int>x<int>\tset number of mpi jobs. Must agree with mpirun\n"
        "\tthr=<int>x<int>\tset number of threads.\n"
        "\tmatrix=<filename>\tset matrix\n"
        "\tseed=<int>\tset random seed\n"
        "\tinterval=<int>\tset checking interval\n"
        );
    exit(1);
}

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    param_list pl;

    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int can_print = rank == 0 || getenv("CAN_PRINT");

    if (can_print) {
        /* print command line */
        fprintf (stderr, "# %s.r%s", argv[0], REV);
        for (int i = 1; i < argc; i++)
            fprintf (stderr, " %s", argv[i]);
        fprintf (stderr, "\n");
    }

    param_list_init (pl);
    argv++, argc--;
    param_list_configure_knob(pl, "-v", &verbose);
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage();
    }

    param_list_parse_string(pl, "matrix", matrix_filename, sizeof(matrix_filename));
    int mpi_split[2] = {1,1,};
    int thr_split[2] = {1,1,};
    int seed=0;

    param_list_parse_intxint(pl, "mpi", mpi_split);
    param_list_parse_intxint(pl, "thr", thr_split);
    param_list_parse_int(pl, "seed", &seed);
    param_list_parse_int(pl, "interval", &interval);

    char subdir[1024] = { 0, };
    param_list_parse_string(pl, "subdir", subdir, sizeof(subdir));
    if (subdir[0]) {
        if (chdir(subdir) < 0) {
            fprintf(stderr, "chdir(%s): %s\n", subdir, strerror(errno));
            exit(1);
        }
    }

    
    int okm=0, okn=0;
    int mn;
    if (param_list_parse_int(pl, "mn", &mn)) {
        m=mn;
        n=mn;
        okm++;
        okn++;
    }
    okm += param_list_parse_int(pl, "m", &m);
    okn += param_list_parse_int(pl, "n", &n);
    if (okm != 1 || okn != 1)
        usage();

    mpz_init_set_ui(p, 2);
    param_list_parse_mpz(pl, "p", p);

    setup_seeding(seed);

    if (verbose)
        param_list_display (pl, stderr);
    if (param_list_warn_unused(pl)) {
        usage();
    }

    param_list_save(pl, "bw-prep.cfg");
    param_list_clear(pl);

    // abase is our arithmetic type.
    abobj_init(abase);

    pi_go(program, mpi_split[0], mpi_split[1], thr_split[0], thr_split[1]);

    mpz_clear(p);
    MPI_Finalize();

    return 0;
}

