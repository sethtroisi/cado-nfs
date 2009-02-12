#define _POSIX_C_SOURCE 200112L

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "bw-common.h"
#include "select_mpi.h"
#include "params.h"
#include "filenames.h"

const char * dirtext[] = { "left", "right" };

/* Has to be defined by the program */
extern void usage();

static int init_mpi = 0;

int bw_common_init_defaults(struct bw_params * bw)
{
    /*** defaults ***/
    memset(bw, 0, sizeof(bw));
    bw->interval = 1000;
    bw->can_print = 1;
    bw->ys[0] = bw->ys[1] = -1;
    bw->nx = 0;
    bw->dir = 1;
    bw->mpi_split[0] = bw->mpi_split[1] = 1;
    bw->thr_split[0] = bw->thr_split[1] = 1;

    return 0;
}

int bw_common_init_shared(struct bw_params * bw, int argc, char * argv[])
{
    if (bw->can_print) {
        /* print command line */
        fprintf (stderr, "# (%s) %s", REV, argv[0]);
        for (int i = 1; i < argc; i++)
            fprintf (stderr, " %s", argv[i]);
        fprintf (stderr, "\n");
    }

    param_list_init (bw->pl);

    argv++, argc--;
    param_list_configure_knob(bw->pl, "-v", &bw->verbose);
    for( ; argc ; ) {
        if (param_list_update_cmdline(bw->pl, &argc, &argv)) { continue; }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage();
    }

    const char * tmp;

    if ((tmp = param_list_lookup_string(bw->pl, "wdir")) != NULL) {
        if (chdir(tmp) < 0) {
            fprintf(stderr, "chdir(%s): %s\n", tmp, strerror(errno));
            exit(1);
        }
    }

    const char * cfg;

    if ((cfg = param_list_lookup_string(bw->pl, "cfg"))) {
        param_list_read_file(bw->pl, cfg);
    } else {
        /* Otherwise we check first that the file exists */
        cfg = BW_CONFIG_FILE;
        if (access(cfg, R_OK) == 0) {
            param_list_read_file(bw->pl, cfg);
        }
    }


    // Here, the matrix filename really ends up in some heap data that will
    // survive after freeing the param_list (if ever we do so early -- as a
    // matter of fact, in prep.c we don't)
    param_list_parse_string(bw->pl, "matrix", bw->matrix_filename, sizeof(bw->matrix_filename));
    param_list_parse_intxint(bw->pl, "mpi", bw->mpi_split);
    param_list_parse_intxint(bw->pl, "thr", bw->thr_split);
    param_list_parse_int(bw->pl, "seed", &bw->seed);
    param_list_parse_int(bw->pl, "interval", &bw->interval);
    param_list_parse_uint(bw->pl, "nx", &bw->nx);
    param_list_parse_int_and_int(bw->pl, "ys", bw->ys, "..");
    param_list_parse_int(bw->pl, "start", &bw->start);
    param_list_parse_int(bw->pl, "end", &bw->end);

    mpz_init_set_ui(bw->p, 2);
    param_list_parse_mpz(bw->pl, "p", bw->p);

    if ((tmp = param_list_lookup_string(bw->pl, "nullspace")) != NULL) {
        if (strcmp(tmp, dirtext[0]) == 0) {
            bw->dir = 0;
        } else if (strcmp(tmp, dirtext[1]) == 0) {
            bw->dir = 1;
        } else {
            fprintf(stderr, "Parameter nullspace may only be %s|%s\n",
                    dirtext[0], dirtext[1]);
            exit(1);
        }
    } else {
        param_list_add_key(bw->pl, "nullspace", dirtext[bw->dir], PARAMETER_FROM_FILE);
    }

    
    int okm=0, okn=0;
    int mn;
    if (param_list_parse_int(bw->pl, "mn", &mn)) {
        bw->m=mn;
        bw->n=mn;
        okm++;
        okn++;
    }
    okm += param_list_parse_int(bw->pl, "m", &bw->m);
    okn += param_list_parse_int(bw->pl, "n", &bw->n);
    if (okm != 1 || okn != 1)
        usage();

    if (bw->verbose)
        param_list_display (bw->pl, stderr);
    if (param_list_warn_unused(bw->pl)) {
        usage();
    }

    return 0;
}

int bw_common_init(struct bw_params * bw, int argc, char * argv[])
{
    bw_common_init_defaults(bw);
    return bw_common_init_shared(bw, argc, argv);
}

int bw_common_init_mpi(struct bw_params * bw, int argc, char * argv[])
{
    bw_common_init_defaults(bw);

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

    bw->can_print = rank == 0 || getenv("CAN_PRINT");

    init_mpi = 1;

    return bw_common_init_shared(bw, argc, argv);
}

int bw_common_clear(struct bw_params * bw)
{
    param_list_clear(bw->pl);

    mpz_clear(bw->p);

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
        "\tm=<int>\tset bw->m blocking factor\n"
        "\tn=<int>\tset bw->n blocking factor\n"
        "\tmn=<int>\tset both bw->m and bw->n (exclusive with the two above)\n"
        "\tmpi=<int>x<int>\tset number of mpi jobs. Must agree with mpiexec\n"
        "\tthr=<int>x<int>\tset number of threads.\n"
        "\tmatrix=<filename>\tset matrix\n"
        "\tinterval=<int>\tset checking bw->interval\n"
        "\tseed=<int>\tseed value for picking random numbers\n"
        "\tys=<int>..<int>\tcoordinate range for krylov/mksol task\n"
        ;
    return t;
}
