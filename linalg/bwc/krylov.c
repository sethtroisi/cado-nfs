#define _POSIX_C_SOURCE 200112L
#define _GNU_SOURCE

#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#include "parallelizing_info.h"
#include "matmul_top.h"
#include "abase.h"
#include "select_mpi.h"
#include "gauss.h"
#include "manu.h"

#include "params.h"
#include "xvectors.h"
#include "xymats.h"

abobj_t abase;

// FIXME: This is a mess. This stuff, and its parsing, is blatantly
// shared between prep, secure, krylov. Should be damn easy to refactor.
int m,n;
mpz_t p;
int interval=1000;
int verbose=0;
char matrix_filename[FILENAME_MAX];
int can_print;

/* This indicates the starting iteration */
int start = 0;

/* relevant to secure.c as well ; really, this is a hard-coded constant.
 * Don't imagine it's tunable. */
#define NCHECKS_CHECK_VECTOR    64

int ys[2];

unsigned int nx;

// dir is a boolean flag equal to 1 if we are looking for the right nullspace
// of the matrix.
int dir = 1;

const char * dirtext[] = { "left", "right" };

void x_dotprod(matmul_top_data_ptr mmt, abt * v, int m, uint32_t * xv, unsigned int nx)
{
    /* We're reading from the shared right vector data -- this area is
     * written to by the other threads in the column. Some of them might
     * be lingering in reduce operations, so we have to wait for them
     */
    serialize_threads(mmt->pi->wr[dir]);

    for(int j = 0 ; j < m ; j++) {
        abt * where = v + aboffset(abase, j);
        for(unsigned int t = 0 ; t < nx ; t++) {
            uint32_t i = xv[j*nx+t];
            if (i < mmt->wr[dir]->i0 || i >= mmt->wr[dir]->i1)
                continue;
            /* We want the data to match our interval on both
             * directions, because otherwise we'll end up
             * computing rubbish -- recall that no broadcast_down
             * has occurred yet.
             */
            if (i < mmt->wr[!dir]->i0 || i >= mmt->wr[!dir]->i1)
                continue;
            abadd(mmt->abase, where,
                    mmt->wr[dir]->v->v + aboffset(abase, i - mmt->wr[dir]->i0));
        }
    }
}

void * krylov_prog(parallelizing_info_ptr pi, void * arg MAYBE_UNUSED)
{
    int tcan_print = can_print && pi->m->trank == 0;
    matmul_top_data mmt;

    int flags[2];
    flags[dir] = THREAD_SHARED_VECTOR;
    flags[!dir] = 0;

    matmul_top_init(mmt, abase, pi, flags, matrix_filename);

    serialize(pi->m);
    
    abzero(abase, mmt->wr[!dir]->v->v, mmt->wr[!dir]->i1 - mmt->wr[!dir]->i0);

    uint32_t * gxvecs = malloc(nx * m * sizeof(uint32_t));

    load_x(gxvecs, m, nx, pi);

    char * v_name;
    asprintf(&v_name, "V%u-%u", ys[0], ys[1]);

    if (tcan_print) { printf("Loading %s...", v_name); fflush(stdout); }
    matmul_top_load_vector(mmt, v_name, dir, start);
    if (tcan_print) { printf("done\n"); }

    // We must also fetch the check vector C.
    abvobj_t placeholder;
    abvobj_set_nbys(placeholder, NCHECKS_CHECK_VECTOR);
    size_t vstride = abvbytes(placeholder, 1);
    size_t stride =  abbytes(abase, 1);

    mmt_generic_vec check_vector;
    matmul_top_vec_init_generic(mmt, abvbytes(placeholder, 1),
            check_vector, !dir, THREAD_SHARED_VECTOR);
    if (tcan_print) { printf("Loading check vector..."); fflush(stdout); }
    matmul_top_load_vector_generic(mmt, vstride,
            check_vector, "C", !dir, interval);
    if (tcan_print) { printf("done\n"); }

    mmt_vec ahead;
    vec_init_generic(pi->m, stride, (mmt_generic_vec_ptr) ahead, 0, NCHECKS_CHECK_VECTOR);

    /* We'll store all xy matrices locally before doing reductions. Given
     * the small footprint of these matrices, it's rather innocuous.
     */
    mmt_vec xymats;
    if (tcan_print) {
        printf("Each thread allocates %zu kb for the A matrices\n",
                abbytes(abase, m*interval) >> 10);
    }
    vec_init_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, 0, m*interval);
    
    for(int s = start ; ; s += interval ) {
        // Plan ahead. The check vector is here to predict the final A matrix.
        // Note that our share of the dot product is determined by the
        // intersections of the i0..i1 intervals on both sides.
        
        unsigned int how_many;
        unsigned int offset_c;
        unsigned int offset_v;
        how_many = intersect_two_intervals(&offset_c, &offset_v,
                mmt->wr[!dir]->i0, mmt->wr[!dir]->i1,
                mmt->wr[dir]->i0, mmt->wr[dir]->i1);

        abzero(abase, ahead->v, NCHECKS_CHECK_VECTOR);
        if (how_many) {
            size_t bytes_c =  abvbytes(placeholder, offset_c);
            abvt * c = abase_generic_ptr_add(check_vector->v, bytes_c);
            abt * v = mmt->wr[dir]->v->v + aboffset(abase, offset_v);
            abvdotprod(abase, placeholder, ahead->v, c, v, how_many);
        }

        /*
        debug_write(abase, ahead->v, NCHECKS_CHECK_VECTOR, "ahead.%u.j%u.t%u",
                s, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        abzero(abase, xymats->v, m*interval);

        serialize(pi->m);
        for(int i = 0 ; i < interval ; i++) {
            /* Compute the product by x */
            x_dotprod(mmt, xymats->v + aboffset(abase, i * m), m, gxvecs, nx);
            matmul_top_mul(mmt, dir);
        }
        serialize(pi->m);

        /*
        debug_write(abase, xymats->v, m * interval, "xy.%u-%u.j%u.t%u",
                s, s + interval, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        /* Last dot product. This must cancel ! */
        x_dotprod(mmt, ahead->v, m, gxvecs, nx);

        /*
        debug_write(abase, ahead->v, NCHECKS_CHECK_VECTOR, "post%u.j%u.t%u",
                s, mmt->pi->m->jrank, mmt->pi->m->trank);
         */

        reduce_generic(abase, ahead, pi->m, NCHECKS_CHECK_VECTOR);
        if (!abis_zero(abase, ahead->v, NCHECKS_CHECK_VECTOR)) {
            printf("Failed check at iteration %d\n", s + interval);
            exit(1);
        }

        /* Now (and only now) collect the xy matrices */
        reduce_generic(abase, xymats, pi->m, m * interval);

        if (pi->m->trank == 0 && pi->m->jrank == 0) {
            char * tmp;
            asprintf(&tmp, "A%u-%u.%u-%u", ys[0], ys[1], s, s+interval);
            FILE * f = fopen(tmp, "w");
            fwrite(xymats->v, stride, aboffset(abase, m*interval), f);
            fclose(f);
        }
        matmul_top_save_vector(mmt, "Y", dir, s + interval);

        if (serialize(pi->m)) {
            printf("Reached iteration %d\n", s + interval);
        }
    }

    if (tcan_print) {
        printf("\n");
    }
    serialize(pi->m);

    matmul_top_vec_clear_generic(mmt, abvbytes(placeholder, 1), check_vector, !dir);
    vec_clear_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, m*interval);
    vec_clear_generic(pi->m, stride, (mmt_generic_vec_ptr) ahead, NCHECKS_CHECK_VECTOR);
    free(gxvecs);

    serialize(pi->m);
    matmul_top_clear(mmt, abase);
    serialize(pi->m);

    return NULL;
}


void usage()
{
    fprintf(stderr, "Usage: ./krylov <options>\n"
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
        "\tys=<int>..<int>\twork on given coordinate in y (right exclusive)\n"
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

    mpz_init_set_ui(p, 2);
    param_list_parse_mpz(pl, "p", p);

    if (!param_list_parse_uint(pl, "nx", &nx)) {
        fprintf(stderr, "no nx value set\n");
        exit(1);
    }

    // Useful to restart from a given checkpoint.
    param_list_parse_int(pl, "start", &start);

    if (!param_list_parse_int_and_int(pl, "ys", ys, "..")) {
        fprintf(stderr, "no ys value set\n");
        exit(1);
    }

    if (verbose)
        param_list_display (pl, stderr);
    if (param_list_warn_unused(pl)) {
        usage();
    }

    // abase is our arithmetic type.
    abobj_init(abase);

    abobj_set_nbys(abase, ys[1]-ys[0]);
    pi_go(krylov_prog, mpi_split[0], mpi_split[1], thr_split[0], thr_split[1], 0);

    // now we're not multithread anymore.

    param_list_clear(pl);

    mpz_clear(p);
    MPI_Finalize();

    return 0;
}

