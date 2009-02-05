#define _POSIX_C_SOURCE 200112L

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "bw-common.h"
#include "select_mpi.h"
#include "params.h"

int m,n;
mpz_t p;
int interval=1000;
int verbose=0;
char matrix_filename[FILENAME_MAX];
int can_print = 1;
int start = 0;
int ys[2] = {-1, -1};
unsigned int nx = 0;
int dir = 1;
const char * dirtext[] = { "left", "right" };
int mpi_split[2] = {1,1,};
int thr_split[2] = {1,1,};
int seed=0;

param_list pl;

/* Has to be defined by the program */
extern void usage();

static int init_mpi = 0;

int bw_common_init(int argc, char * argv[])
{
    if (can_print) {
        /* print command line */
        fprintf (stderr, "# (%s) %s", REV, argv[0]);
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

    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "wdir")) != NULL) {
        if (chdir(tmp) < 0) {
            fprintf(stderr, "chdir(%s): %s\n", tmp, strerror(errno));
            exit(1);
        }
    }

    const char * cfg;

    if ((cfg = param_list_lookup_string(pl, "cfg"))) {
        param_list_read_file(pl, cfg);
    } else {
        /* Otherwise we check first that the file exists */
        cfg = "bw.cfg";
        if (access(cfg, R_OK) == 0) {
            param_list_read_file(pl, cfg);
        }
    }


    // Here, the matrix filename really ends up in some heap data that will
    // survive after freeing the param_list (if ever we do so early -- as a
    // matter of fact, in prep.c we don't)
    param_list_parse_string(pl, "matrix", matrix_filename, sizeof(matrix_filename));
    param_list_parse_intxint(pl, "mpi", mpi_split);
    param_list_parse_intxint(pl, "thr", thr_split);
    param_list_parse_int(pl, "seed", &seed);
    param_list_parse_int(pl, "interval", &interval);
    param_list_parse_uint(pl, "nx", &nx);
    param_list_parse_int_and_int(pl, "ys", ys, "..");
    param_list_parse_int(pl, "start", &start);

    mpz_init_set_ui(p, 2);
    param_list_parse_mpz(pl, "p", p);

    if ((tmp = param_list_lookup_string(pl, "nullspace")) != NULL) {
        if (strcmp(tmp, dirtext[0]) == 0) {
            dir = 0;
        } else if (strcmp(tmp, dirtext[1]) == 0) {
            dir = 1;
        } else {
            fprintf(stderr, "Parameter nullspace may only be %s|%s\n",
                    dirtext[0], dirtext[1]);
            exit(1);
        }
    } else {
        param_list_add_key(pl, "nullspace", dirtext[dir], PARAMETER_FROM_FILE);
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

    if (verbose)
        param_list_display (pl, stderr);
    if (param_list_warn_unused(pl)) {
        usage();
    }

    return 0;
}

int bw_common_init_mpi(int argc, char * argv[])
{
#ifdef  MPI_LIBRARY_MT_CAPABLE
    int req = MPI_THREAD_MULTIPLE;
    int prov;
    MPI_Init_thread(&argc, &argv, req, &prov);
    if (req != prov) {
        fprintf(stderr, "Cannot init mpi with MPI_THREAD_MULTIPLE ;"
                " got %d != req %d\n",
                prov, req);
        exit(1);
    }
#else
    MPI_Init(&argc, &argv);
#endif
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    can_print = rank == 0 || getenv("CAN_PRINT");

    init_mpi = 1;

    bw_common_init(argc, argv);

    return 0;
}

int bw_common_clear()
{
    param_list_clear(pl);

    mpz_clear(p);

    if (init_mpi) {
        MPI_Finalize();
    }

    return 0;
}

const char * bw_common_usage_string()
{
    static char t[]=
        "Common options:\n"
        "\twdir=<path>\tchdir to <path> beforehand\n"
        "\tcfg=<file>\timport many settings from <file>\n"
        "\tm=<int>\tset m blocking factor\n"
        "\tn=<int>\tset n blocking factor\n"
        "\tmn=<int>\tset both m and n (exclusive with the two above)\n"
        "\tmpi=<int>x<int>\tset number of mpi jobs. Must agree with mpiexec\n"
        "\tthr=<int>x<int>\tset number of threads.\n"
        "\tmatrix=<filename>\tset matrix\n"
        "\tinterval=<int>\tset checking interval\n"
        "\tseed=<int>\tseed value for picking random numbers\n"
        "\tys=<int>..<int>\tcoordinate range for krylov/mksol task\n"
        ;
    return t;
}
