/* Shirokauer maps 

   This is largely copied from sm_simple.cpp ; however the purpose of
   sm_simple was to be kept simple and stupid, so let's keep it like
   this.

   This program is to be used as an mpi accelerator for reconstructlog-dl

   Given a relation file where each line begins with an a,b pair IN
   HEXADECIMAL FORMAT, output the same line, appending the SM values in
   the end.

   SM computation is offloaded to the (many) MPI jobs.

*/

#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>

#include <vector>
#include <string>

#include "macros.h"
#include "utils.h"
#include "select_mpi.h"
#include "relation.h"

#define BATCH_SIZE 128

using namespace std;

struct ab_pair {
    int64_t a;		/* only a is allowed to be negative */
    uint64_t b;
};

typedef vector<ab_pair> ab_pair_batch;

static void sm_append_master(FILE * in, FILE * out, sm_side_info *sm_info, int nb_polys, int size)
{
    char buf[1024];

    /* need to know how many mp_limb_t's we'll get back from each batch */
    size_t limbs_per_ell = 0;
    int nsm_total=0;
    for(int side = 0; side < nb_polys; side++) {
        nsm_total += sm_info[side]->nsm;
        if (sm_info[side]->nsm) limbs_per_ell = mpz_size(sm_info[side]->ell);
    }

    std::vector<MPI_Request> active(size);
    std::vector<ab_pair_batch> ab_pair_batches(size);
    std::vector<string> history;

    bool eof = false;

    for(int turn = 0 ; !eof ; turn++) {
        history.clear();
        /* {{{ Collect batches, and send them */
        for(int peer = 1; peer < size; peer++) {
            ASSERT_ALWAYS(ab_pair_batches[peer].empty());
            while (!eof &&
                    ab_pair_batches[peer].size() < BATCH_SIZE &&
                    fgets(buf, 1024, in))
            {
                int n = strlen(buf);
                if (!n) {
                    fprintf(stderr, "Got 0-sized buffer in fgets, shouldn't happen. Assuming EOF.\n");
                    eof=true;
                    break;
                }
                buf[n-1]='\0';

                if (buf[0] == '#') {
                    fputs(buf, out);
                    fputc('\n', out);
                    continue;
                }

                history.push_back(string(buf));

                char * p = buf;
                ab_pair ab;
                int64_t sign = 1;
                if (*p=='-') {
                    sign=-1;
                    p++;
                }
                if (sscanf(p, "%" SCNx64 ",%" SCNx64 ":", &ab.a, &ab.b) < 2) {
                    fprintf(stderr, "Parse error at line: %s\n", buf);
                    exit(EXIT_FAILURE);
                }
                ab.a *= sign;

                ab_pair_batches[peer].push_back(ab);
            }
            if (!eof && ab_pair_batches[peer].size() < BATCH_SIZE) {
                eof=true;
                if (ferror(stdin)) {
                    fprintf(stderr, "Error on stdin\n");
                }
            }
            /* 0 bsize will be recognized by slaves as a 
             * reason to stop processing */
            unsigned long bsize = ab_pair_batches[peer].size();
            MPI_Send(&bsize, 1, MPI_UNSIGNED_LONG, peer, turn, MPI_COMM_WORLD);
            if (bsize)
                MPI_Isend((char*) &(ab_pair_batches[peer][0]), bsize * sizeof(ab_pair), MPI_BYTE, peer, turn, MPI_COMM_WORLD, &active[peer]);

        }
        /* }}} */

        /* Now receive the batches from the clients and print the results */
        size_t histidx = 0;
        for(int peer = 1; peer < size; peer++) {
            unsigned long bsize = ab_pair_batches[peer].size();

            if (bsize == 0)
                continue;

            MPI_Wait(&active[peer], MPI_STATUS_IGNORE);
            ab_pair_batches[peer].clear();

            mp_limb_t returns[bsize][nsm_total][limbs_per_ell];

            MPI_Recv(returns, bsize * nsm_total * limbs_per_ell * sizeof(mp_limb_t), MPI_BYTE, peer, turn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(unsigned long i = 0 ; i < bsize ; i++) {
                fputs(history[histidx++].c_str(), out);
                bool comma = false;
                for (int j = 0 ; j < nsm_total ; j++) {
                    gmp_fprintf(out, "%c%Nd", comma ? ',' : ':', returns[i][j], limbs_per_ell);
                    comma=true;
                }
                fputc('\n', out);
            }
            if (eof) {
                /* also tell this one to finish */
                bsize = 0;
                MPI_Send(&bsize, 1, MPI_UNSIGNED_LONG, peer, turn+1, MPI_COMM_WORLD);
            }
        }
        history.clear();
    }
}

static void sm_append_slave(sm_side_info *sm_info, int nb_polys)
{
    /* need to know how many mp_limb_t's we'll get back from each batch */
    size_t limbs_per_ell = 0;
    int nsm_total=0;
    int maxdeg = 0;

    for(int side = 0; side < nb_polys; side++) {
        nsm_total += sm_info[side]->nsm;
        maxdeg = MAX(maxdeg, sm_info[side]->f->deg);
        if (sm_info[side]->nsm) limbs_per_ell = mpz_size(sm_info[side]->ell);
    }

    mpz_poly_t smpol;

    mpz_poly_init(smpol, maxdeg);

    for(int turn = 0 ; ; turn++) {
        unsigned long bsize;
        MPI_Recv(&bsize, 1, MPI_UNSIGNED_LONG, 0, turn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (bsize == 0)
            break;
        ab_pair_batch batch(bsize);
        MPI_Recv((char*) &(batch[0]), bsize * sizeof(ab_pair), MPI_BYTE, 0, turn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        mp_limb_t returns[bsize][nsm_total][limbs_per_ell];
        memset(returns, 0, bsize*nsm_total*limbs_per_ell*sizeof(mp_limb_t));

        for(unsigned long i = 0 ; i < bsize ; i++) {
            mpz_poly_t pol;
            mpz_poly_init_set_ab(pol, batch[i].a, batch[i].b);
            int smidx = 0;
            for (int side = 0; side < nb_polys; ++side) {
                compute_sm_piecewise(smpol, pol, sm_info[side]);
                for(int k = 0 ; k < sm_info[side]->nsm ; k++, smidx++) {
                    if (k <= smpol->deg) {
                        for(size_t j = 0 ; j < limbs_per_ell ; j++) {
                            returns[i][smidx][j] = mpz_getlimbn(smpol->coeff[smpol->deg-k], j);
                        }
                    }
                }
            }
            mpz_poly_clear(pol);
        }

        MPI_Send(returns, bsize * nsm_total * limbs_per_ell * sizeof(mp_limb_t), MPI_BYTE, 0, turn, MPI_COMM_WORLD);
    }
    mpz_poly_clear(smpol);
}

static void sm_append_sync(FILE * in, FILE * out, sm_side_info *sm_info, int nb_polys)
{
    char buf[1024];
    mpz_poly_t pol, smpol;
    int maxdeg = sm_info[0]->f->deg;
    for(int side = 1; side < nb_polys; side++)
        maxdeg = MAX(maxdeg, sm_info[side]->f->deg);
    mpz_poly_init(pol, maxdeg);
    mpz_poly_init(smpol, maxdeg);
    while (fgets(buf, 1024, in)) {
        int n = strlen(buf);
        if (!n) break;
        buf[n-1]='\0';

        if (buf[0] == '#') {
            fputs(buf, out);
            fputc('\n', out);
            continue;
        }

        char * p = buf;
        int64_t a;		/* only a is allowed to be negative */
        uint64_t b;
        int64_t sign = 1;
        if (*p=='-') {
            sign=-1;
            p++;
        }
        if (sscanf(p, "%" SCNx64 ",%" SCNx64 ":", &a, &b) < 2) {
            fprintf(stderr, "Parse error at line: %s\n", buf);
            exit(EXIT_FAILURE);
        }

        mpz_poly_init_set_ab(pol, a*sign, b);

        fputs(buf, out);
        fputc(':', out);
        for (int side = 0; side < nb_polys; ++side) {
            compute_sm_piecewise(smpol, pol, sm_info[side]);
            print_sm2(out, smpol, sm_info[side]->nsm, sm_info[side]->f->deg, ",");
            if (side == 0 && sm_info[0]->nsm > 0 && sm_info[1]->nsm > 0)
                fputc(',', out);
        }
        fputc('\n', out);
        mpz_poly_clear(pol);
    }
    mpz_poly_clear(smpol);
}


static void sm_append(FILE * in, FILE * out, sm_side_info *sm_info, int nb_polys)
{
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size > 1) {
        if (rank == 0) {
            sm_append_master(in, out, sm_info, nb_polys, size);
        } else {
            sm_append_slave(sm_info, nb_polys);
        }
    } else {
        sm_append_sync(in, out, sm_info, nb_polys);
    }
}


static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "(required) poly file");
    param_list_decl_usage(pl, "ell", "(required) group order");
    param_list_decl_usage(pl, "nsm", "number of SMs to use per side");
    param_list_decl_usage(pl, "in", "data input (defaults to stdin)");
    param_list_decl_usage(pl, "out", "data output (defaults to stdout)");
    verbose_decl_usage(pl);
}

static void usage (const char *argv, const char * missing, param_list pl)
{
    if (missing) {
        fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
                missing);
    }
    param_list_print_usage(pl, argv, stderr);
    exit (EXIT_FAILURE);
}

/* -------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
    MPI_Init(&argc, & argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char *argv0 = argv[0];

    const char *polyfile = NULL;

    param_list pl;
    cado_poly pol;
    mpz_poly_ptr F[NB_POLYS_MAX];

    mpz_t ell;
    double t0;

    /* read params */
    param_list_init(pl);
    declare_usage(pl);

    if (argc == 1)
        usage (argv[0], NULL, pl);

    argc--,argv++;
    for ( ; argc ; ) {
        if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
        usage (argv0, NULL, pl);
    }

    /* Read poly filename from command line */
    if ((polyfile = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: parameter -poly is mandatory\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    /* Read ell from command line (assuming radix 10) */
    mpz_init (ell);
    if (!param_list_parse_mpz(pl, "ell", ell)) {
        fprintf(stderr, "Error: parameter -ell is mandatory\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    /* Init polynomial */
    cado_poly_init (pol);
    cado_poly_read(pol, polyfile);
    for(int side = 0; side < pol->nb_polys; side++)
        F[side] = pol->pols[side];

    int nsm_arg[NB_POLYS_MAX];
    for(int side = 0; side < pol->nb_polys; side++)
        nsm_arg[side]=-1;

    param_list_parse_int_list (pl, "nsm", nsm_arg, pol->nb_polys, ",");

    FILE * in = rank ? NULL : stdin;
    FILE * out = rank ? NULL: stdout;
    const char * infilename = param_list_lookup_string(pl, "in");
    const char * outfilename = param_list_lookup_string(pl, "out");

    if (!rank && infilename) {
        in = fopen_maybe_compressed(infilename, "r");
        ASSERT_ALWAYS(in != NULL);
    }
    if (!rank && outfilename) {
        out = fopen_maybe_compressed(outfilename, "w");
        ASSERT_ALWAYS(out != NULL);
    }

    if (param_list_warn_unused(pl))
        usage (argv0, NULL, pl);

    verbose_interpret_parameters(pl);

    if (!rank)
    param_list_print_command_line (stdout, pl);

    sm_side_info sm_info[NB_POLYS_MAX];

    for(int side = 0 ; side < pol->nb_polys; side++) {
        sm_side_info_init(sm_info[side], F[side], ell);
        if (nsm_arg[side] >= 0)
            sm_info[side]->nsm = nsm_arg[side]; /* command line wins */
        if (!rank)
            printf("# Using %d SMs on side %d\n", sm_info[side]->nsm, side);
    }

    /*
       if (!rank) {
           for (int side = 0; side < pol->nb_polys; side++) {
               printf("\n# Polynomial on side %d:\nF[%d] = ", side, side);
               mpz_poly_fprintf(stdout, F[side]);

               printf("# SM info on side %d:\n", side);
               sm_side_info_print(stdout, sm_info[side]);

               fflush(stdout);
           }
       }
       */

    t0 = wct_seconds();

    sm_append(in, out, sm_info, pol->nb_polys);

    if (!rank) {
        printf("\n# sm completed in %2.2lf seconds (WCT)\n", wct_seconds() - t0);
        fflush(stdout);
    }

    if (!rank && infilename) fclose_maybe_compressed(in, infilename);
    if (!rank && out != stdout) fclose_maybe_compressed(out, outfilename);

    for(int side = 0 ; side < pol->nb_polys ; side++) {
        sm_side_info_clear(sm_info[side]);
    }

    mpz_clear(ell);
    cado_poly_clear(pol);
    param_list_clear(pl);

    MPI_Finalize();

    return 0;
}
