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
#include "xvectors.h"

abobj_t abase;

int m,n;
mpz_t p;
int interval=1000;
int verbose=0;
char matrix_filename[FILENAME_MAX];
int can_print;

unsigned int nx_main;

// dir is a boolean flag equal to 1 if we are looking for the right nullspace
// of the matrix.
int dir = 1;

const char * dirtext[] = { "left", "right" };

void * program(parallelizing_info_ptr pi, void * arg MAYBE_UNUSED)
{
    // Doing the ``hello world'' test is a very good way of testing the
    // global mpi/pthreads setup. So despite its apparent irrelevance, I
    // suggest leaving it here as a cheap sanity check.
    hello(pi);

    // avoid cluttering output too much.
    int tcan_print = can_print && pi->m->trank == 0;

    unsigned int nx = 1;

    matmul_top_data mmt;

    int flags[2] = { 0, THREAD_SHARED_VECTOR };

    // FIXME: we have set no multiplication algorithm here.
    matmul_top_init(mmt, abase, NULL, pi, flags, matrix_filename);

    uint32_t * xvecs = malloc(nx * m * sizeof(uint32_t));

    for (unsigned ntri = 0;; ntri++) {
        if (tcan_print) {
            printf("// Generating new x,y vector pair\n");
        }
        if (ntri >= nx * 10) {
            ++nx;
            if (tcan_print) {
                printf("// Getting bored. Trying %u x vectors\n", nx);
            }
            xvecs = realloc(xvecs, nx * m * sizeof(uint32_t));
            ASSERT_ALWAYS(xvecs != NULL);
        }

        // if we're looking for the right nullspace, then x is on the left.
        // Otherwise, it's on the right
        setup_x_random(xvecs, m, nx, mmt->n[dir], pi);

        break;
    }

    matmul_top_fill_random_source(mmt, dir);

    // we need to save this starting vector for later use if it turns out
    // that we need to save it for real.

    matmul_top_save_vector(mmt, dir, 0, 0);
    serialize(pi->m);

    matmul_top_load_vector(mmt, dir, 0, 0);
    matmul_top_save_vector(mmt, dir, 0, 1);

    // last 
    serialize(pi->m);

    if (pi->m->trank == 0) {
        nx_main = nx;
    }

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

    can_print = rank == 0 || getenv("CAN_PRINT");

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

    // Here, the matrix filename really ends up in some heap data that will
    // survive after freeing the param_list (if ever we do so early -- as a
    // matter of fact, in prep.c we don't)
    param_list_parse_string(pl, "matrix", matrix_filename, sizeof(matrix_filename));
    int mpi_split[2] = {1,1,};
    int thr_split[2] = {1,1,};
    int seed=0;

    param_list_parse_intxint(pl, "mpi", mpi_split);
    param_list_parse_intxint(pl, "thr", thr_split);
    param_list_parse_int(pl, "seed", &seed);
    param_list_parse_int(pl, "interval", &interval);

    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "nullspace")) != NULL) {
        if (strcmp(tmp, dirtext[0])) {
            dir = 0;
        } else if (strcmp(tmp, dirtext[1])) {
            dir = 1;
        } else {
            fprintf(stderr, "Parameter nullspace may only be %s|%s\n",
                    dirtext[0], dirtext[1]);
            exit(1);
        }
    } else {
        param_list_add_key(pl, "nullspace", dirtext[dir], PARAMETER_FROM_FILE);
    }

    if ((tmp = param_list_lookup_string(pl, "tmp")) != NULL) {
        if (chdir(tmp) < 0) {
            fprintf(stderr, "chdir(%s): %s\n", tmp, strerror(errno));
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

    // abase is our arithmetic type.
    abobj_init(abase);

    pi_go(program, mpi_split[0], mpi_split[1], thr_split[0], thr_split[1], 0);

    // now we're not multithread anymore.

    // we save the parameter list once again, because the prep program
    // generates some useful info.
    char nxstring[10];
    snprintf(nxstring, sizeof(nxstring), "%u", nx_main);
    param_list_add_key(pl, "nx", nxstring, PARAMETER_FROM_FILE);
    param_list_save(pl, "bw.cfg");
    param_list_clear(pl);

    mpz_clear(p);
    MPI_Finalize();

    return 0;
}

