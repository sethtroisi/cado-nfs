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
#include "random_generation.h"
#include "gauss.h"
#include "manu.h"
#include "gauss.h"

#include "params.h"
#include "xvectors.h"
#include "xymats.h"

/* Number of copies of m by n matrices to use for trying to obtain a
 * full-rank matrix (rank m).
 *
 * Note that it must be at least m/n, otherwise we stand no chance !
 */
#define NBITER  2

abobj_t abase;

int m,n;
mpz_t p;
int verbose=0;
char matrix_filename[FILENAME_MAX];
int can_print;

unsigned int nx_main;

// dir is a boolean flag equal to 1 if we are looking for the right nullspace
// of the matrix.
int dir = 1;

const char * dirtext[] = { "left", "right" };

#if 0
#include "hexstring.h"
void display_matrix(abobj_ptr abase, mmt_vec_ptr xy, pi_wiring_ptr wr, unsigned int nr, unsigned int nvblocks)
{
    for(int z = 0 ; z < wr->ncores ; z++) {
        serialize_threads(wr);
        if (z != wr->trank) continue;
        for(unsigned int i = 0 ; i < m ; i++) {
            printf("(% 2d)", z);
            for(unsigned int k = 0 ; k < nvblocks ; k++) {
                unsigned long * data;
                data = xy->v + aboffset(abase, i * NBITER + k);
                putchar(' ');
                write_hexstring(stdout, (unsigned long *) data, abnbits(abase));
            }
            putchar('\n');
        }
        putchar('\n');
    }
    serialize_threads(wr);
}
#endif

void * prep_prog(parallelizing_info_ptr pi, void * arg MAYBE_UNUSED)
{
    // Doing the ``hello world'' test is a very good way of testing the
    // global mpi/pthreads setup. So despite its apparent irrelevance, I
    // suggest leaving it here as a cheap sanity check.
    hello(pi);

    // avoid cluttering output too much.
    int tcan_print = can_print && pi->m->trank == 0;

    unsigned int nx = 1;

    matmul_top_data mmt;

    int flags[2];
    flags[dir] = THREAD_SHARED_VECTOR;
    flags[!dir] = 0;

    matmul_top_init(mmt, abase, pi, flags, matrix_filename);

    uint32_t * xvecs = malloc(nx * m * sizeof(uint32_t));

    mmt_vec xymats;
    size_t stride =  abbytes(abase, 1);

    /* We're cheating on the generic init routines */
    vec_init_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, 0, m*NBITER);

    for (unsigned ntri = 0;; ntri++) {
        serialize_threads(pi->m);

        if (tcan_print) {
            printf("// Generating new x,y vector pair (trial # %u)\n", ntri);
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
        // Otherwise, it's on the right.
        setup_x_random(xvecs, m, nx, mmt->n[dir], pi);

        // Compute y.
        matmul_top_fill_random_source(mmt, dir);

        // we need to save this starting vector for later use if it turns out
        // that we need to save it for real.
        matmul_top_save_vector(mmt, "Y", dir, 0);

        // We must compute x^T M y, x^T M^2 y, and so on.
        
        // we have indices mmt->wr[1]->i0..i1 available.
        abzero(abase, xymats->v, m * NBITER);


        for(unsigned int k = 0 ; k < NBITER ; k++) {
#if 0
            {
                char * tmp;
                asprintf(&tmp, "/tmp/y%d.%d.%d.try%d", k, pi->m->jrank, pi->m->trank, ntri);
                FILE * f = fopen(tmp, "w");
                fwrite(mmt->wr[dir]->v, sizeof(abt), aboffset(abase, mmt->wr[dir]->i1 - mmt->wr[dir]->i0), f);
                fclose(f);
            }
            serialize_threads(pi->m);
#endif
    
            for(int j = 0 ; j < m ; j++) {
                abt * where = xymats->v + aboffset(abase, j * NBITER + k);
                for(unsigned int t = 0 ; t < nx ; t++) {
                    uint32_t i = xvecs[j*nx+t];
                    if (i < mmt->wr[dir]->i0 || i >= mmt->wr[dir]->i1)
                        continue;
                    /* We want the data to match our interval on both
                     * directions, because otherwise we'll end up
                     * computing rubbish -- recall that no broadcast_down
                     * has occurred yet.
                     */
                    if (i < mmt->wr[!dir]->i0 || i >= mmt->wr[!dir]->i1)
                        continue;
                    abadd(abase, where,
                            mmt->wr[dir]->v->v + aboffset(abase, i - mmt->wr[dir]->i0));
                }
            }
            matmul_top_mul(mmt, dir);
        }

        /* Make sure computation is over for everyone ! */
        serialize_threads(pi->m);

        /* Now all threads and jobs must collectively reduce the zone
         * pointed to by xymats */
        reduce_generic(abase, xymats, pi->m, m * NBITER);

        /* OK -- now everybody has the same data, as can be seen for
         * instance with the following debugging code.  */
        {
            char * tmp;
            asprintf(&tmp, "/tmp/xy.%d.%d.try%d.postred", pi->m->jrank, pi->m->trank, ntri);
            FILE * f = fopen(tmp, "w");
            fwrite(xymats->v, sizeof(abt), aboffset(abase, m * NBITER), f);
            fclose(f);
        }


        int dimk;
        int * pdimk;
        
        /* the kernel() call is not reentrant */
        if (pi->m->trank == 0) {
            dimk = kernel((mp_limb_t *) xymats->v, NULL,
                    m, NBITER * abnbits(abase), 
                abbytes(abase,NBITER) / sizeof(mp_limb_t), 0);
            pdimk = &dimk;
        }
        thread_agreement(pi->m, (void **) &pdimk, 0);
        dimk = * pdimk;

        if (tcan_print)
            printf("// Dimension of kernel: %d\n", dimk);

        if (dimk == 0) {
            if (tcan_print)
                printf("// Found good x,y vector pair after %u trials\n",
                        ntri+1);
            break;
        }
    }

    if (pi->m->trank == 0) {
        nx_main = nx;
    }

    save_x(xvecs, m, nx, pi);

    // last 
    serialize(pi->m);
    matmul_top_clear(mmt, abase);
    serialize(pi->m);

    /* clean up xy mats stuff */
    vec_clear_generic(pi->m, stride, (mmt_generic_vec_ptr) xymats, m*NBITER);

    free(xvecs);
    return NULL;
}

void usage()
{
    fprintf(stderr, "Usage: ./prep <options>\n"
        "Allowed options are\n"
        "\tm=<int>\t(*) set m blocking factor\n"
        "\tn=<int>\t(*) set n blocking factor\n"
        "\tmn=<int>\tset both m and n to same value (exclusive with the two above)\n"
        "\twdir=<path>\tchdir to <path> beforehand\n"
        "\tmpi=<int>x<int>\tset number of mpi jobs. Must agree with mpiexec\n"
        "\tthr=<int>x<int>\tset number of threads.\n"
        "\tmatrix=<filename>\tset matrix\n"
        "\tseed=<int>\tset random seed\n"
        "\tinterval=<int>\tset checking interval (not used by prep)\n"
        );
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

    setup_seeding(seed);

    // we don't use this parameter, but we tolerate it, since it's going
    // to be used by the secure program afterwards. By looking it up, we
    // force the param_list routines to accept it, so that it'll be
    // written to the config file eventually.
    param_list_lookup_string(pl, "interval");

    if (verbose)
        param_list_display (pl, stderr);
    if (param_list_warn_unused(pl)) {
        usage();
    }

    // abase is our arithmetic type.
    abobj_init(abase);

    abobj_set_nbys(abase, n);
    pi_go(prep_prog, mpi_split[0], mpi_split[1], thr_split[0], thr_split[1], 0);

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

