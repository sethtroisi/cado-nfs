#define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE

#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>

#include "parallelizing_info.h"
#include "matmul_top.h"
#include "abase.h"
#include "select_mpi.h"
#include "gauss.h"
//#include "hexstring.h"
#include "manu.h"

#include "params.h"
#include "xvectors.h"

abobj_t abase;

int m,n;
mpz_t p;
int interval=1000;
int verbose=0;
char matrix_filename[FILENAME_MAX];
int can_print;

unsigned int nx;

// dir is a boolean flag equal to 1 if we are looking for the right nullspace
// of the matrix.
int dir = 1;

const char * dirtext[] = { "left", "right" };

void * sec_prog(parallelizing_info_ptr pi, void * arg MAYBE_UNUSED)
{
    /* OK. Now we merely have to build up a check vector. We've got a
     * fairly easy candidate for that: the x vector.
     */
    int tcan_print = can_print && pi->m->trank == 0;
    matmul_top_data mmt;

    int flags[2];
    flags[!dir] = THREAD_SHARED_VECTOR;
    flags[dir] = 0;

    matmul_top_init(mmt, abase, pi, flags, matrix_filename);

    abzero(abase, mmt->wr[!dir]->v, mmt->wr[!dir]->i1 - mmt->wr[!dir]->i0);

    uint32_t * gxvecs = malloc(nx * m * sizeof(uint32_t));

    load_x(gxvecs, m, nx, pi);

    ASSERT_ALWAYS(m >= 64);
    for(int j = 0 ; j < 64 ; j++) {
        for(unsigned int k = 0 ; k < nx ; k++) {
            // set bit j of entry gxvecs[j*nx+k] to 1.
            uint32_t i = gxvecs[j*nx+k];
            if (i < mmt->wr[!dir]->i0 || i >= mmt->wr[!dir]->i1)
                continue;
            abt * where;
            where = mmt->wr[!dir]->v + aboffset(abase, i-mmt->wr[!dir]->i0);
            abset_ui(abase, where, j, 1);
        }
    }

    if (tcan_print) {
        printf("Computing trsp(x)*M^%d\n",interval);
    }

    serialize(pi->m);

    for(int k = 0 ; k < interval ; k++) {
        matmul_top_mul(mmt, !dir);
        if (tcan_print) {
            putchar('.');
            fflush(stdout);
        }
    }
    if (tcan_print) {
        printf("\n");
    }
    serialize(pi->m);

    matmul_top_save_vector(mmt, "C", !dir, 0, interval);

    serialize(pi->m);
    matmul_top_clear(mmt, abase);
    serialize(pi->m);

    free(gxvecs);

    return NULL;
}


void usage()
{
    fprintf(stderr, "Usage: ./secure <options>\n"
        "Allowed options are\n"
        "\tm=<int>\t(*) set m blocking factor\n"
        "\tn=<int>\t(*) set n blocking factor\n"
        "\tmn=<int>\tset both m and n (exclusive with the two above)\n"
        "\twdir=<path>\tchdir to <path> beforehand\n"
        "\tmpi=<int>x<int>\tset number of mpi jobs. Must agree with mpiexec\n"
        "\tthr=<int>x<int>\tset number of threads.\n"
        "\tmatrix=<filename>\tset matrix\n"
        "\tinterval=<int>\tset checking interval\n"
        "\tcfg=<file>\timport many settings from <file>\n"
        );
    fprintf(stderr, "Note: data files must be found in wdir !\n");
    exit(1);
}

int main(int argc, char * argv[])
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

    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "wdir")) != NULL) {
        if (chdir(tmp) < 0) {
            fprintf(stderr, "chdir(%s): %s\n", tmp, strerror(errno));
            exit(1);
        }
    }

    const char * cfg;

    if (!(cfg = param_list_lookup_string(pl, "cfg"))) {
        cfg = "bw.cfg";
    }
    param_list_read_file(pl, cfg);

    // Here, the matrix filename really ends up in some heap data that will
    // survive after freeing the param_list (if ever we do so early -- as a
    // matter of fact, in prep.c we don't)
    param_list_parse_string(pl, "matrix", matrix_filename, sizeof(matrix_filename));
    int mpi_split[2] = {1,1,};
    int thr_split[2] = {1,1,};

    param_list_parse_intxint(pl, "mpi", mpi_split);
    param_list_parse_intxint(pl, "thr", thr_split);
    param_list_parse_int(pl, "interval", &interval);

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

    if (!param_list_parse_uint(pl, "nx", &nx)) {
        fprintf(stderr, "no nx value set\n");
        exit(1);
    }

    if (verbose)
        param_list_display (pl, stderr);
    if (param_list_warn_unused(pl)) {
        usage();
    }

    // abase is our arithmetic type.
    abobj_init(abase);

    // now do the computation of the iterates for the check vector.
    // NOTE that we hard-code the program to handle only 64-bit check
    // vectors. That's good enough.
    abobj_set_nbys(abase, 64);
    pi_go(sec_prog, mpi_split[0], mpi_split[1], thr_split[0], thr_split[1], 0);

    // now we're not multithread anymore.

    // we save the parameter list once again, in case we inherited some
    // new parameter. Is there any point in doing so ???
    param_list_save(pl, "bw.cfg");
    param_list_clear(pl);

    mpz_clear(p);
    MPI_Finalize();

    return 0;
}

