#define _POSIX_C_SOURCE 200112L

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "bw-common.h"
#include "params.h"
#include "filenames.h"

struct bw_params bw[1];

const char * dirtext[] = { "left", "right" };

/* Has to be defined by the program */
extern void usage();

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

    bw->checkpoints = 1;

    return 0;
}

int bw_common_init_shared(struct bw_params * bw, param_list pl, int * p_argc, char *** p_argv)
{
    if (bw->can_print) {
        /* print command line */
        fprintf (stderr, "# (%s) %s", REV, (*p_argv)[0]);
        for (int i = 1; i < (*p_argc); i++)
            fprintf (stderr, " %s", (*p_argv)[i]);
        fprintf (stderr, "\n");
    }

    (*p_argv)++, (*p_argc)--;
    param_list_configure_knob(pl, "-v", &bw->verbose);
    for( ; (*p_argc) ; ) {
        if (param_list_update_cmdline(pl, p_argc, p_argv)) { continue; }
        if (strcmp((*p_argv)[0],"--") == 0) {
            (*p_argv)++, (*p_argc)--;
            break;
        }
        fprintf(stderr, "Unhandled parameter %s\n", (*p_argv)[0]);
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
        cfg = BW_CONFIG_FILE;
        if (access(cfg, R_OK) == 0) {
            param_list_read_file(pl, cfg);
        }
    }


    // Here, the matrix filename really ends up in some heap data that will
    // survive after freeing the param_list (if ever we do so early -- as a
    // matter of fact, in prep.c we don't)
    param_list_parse_intxint(pl, "mpi", bw->mpi_split);
    param_list_parse_intxint(pl, "thr", bw->thr_split);
    param_list_parse_int(pl, "seed", &bw->seed);
    param_list_parse_int(pl, "interval", &bw->interval);
    param_list_parse_uint(pl, "nx", &bw->nx);
    param_list_parse_int_and_int(pl, "ys", bw->ys, "..");
    param_list_parse_int(pl, "start", &bw->start);
    param_list_parse_int(pl, "end", &bw->end);
    param_list_parse_int(pl, "checkpoints", &bw->checkpoints);


    mpz_init_set_ui(bw->p, 2);
    param_list_parse_mpz(pl, "p", bw->p);

    if ((tmp = param_list_lookup_string(pl, "nullspace")) != NULL) {
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
        param_list_add_key(pl, "nullspace", dirtext[bw->dir], PARAMETER_FROM_FILE);
    }

    
    int okm=0, okn=0;
    int mn;
    if (param_list_parse_int(pl, "mn", &mn)) {
        bw->m=mn;
        bw->n=mn;
        okm++;
        okn++;
    }
    okm += param_list_parse_int(pl, "m", &bw->m);
    okn += param_list_parse_int(pl, "n", &bw->n);
    if (!okm || !okn)
        usage();

    if (bw->verbose && bw->can_print)
        param_list_display (pl, stderr);

    /* Force lookup of parameters that are used late in the process. This
     * has the effect of eliminating possible warnings. */
    param_list_lookup_string(pl, "bwc_mm_impl");

    param_list_lookup_string(pl, "bwc_mm_l1_cache_size");
    param_list_lookup_string(pl, "bwc_mm_cache_line_size");

    param_list_lookup_string(pl, "bwc_mm_threaded_nthreads");
    param_list_lookup_string(pl, "bwc_mm_threaded_sgroup_size");
    param_list_lookup_string(pl, "bwc_mm_threaded_offset1");
    param_list_lookup_string(pl, "bwc_mm_threaded_offset2");
    param_list_lookup_string(pl, "bwc_mm_threaded_offset3");
    param_list_lookup_string(pl, "bwc_mm_threaded_densify_tolerance");
    param_list_lookup_string(pl, "bwc_mm_store_transposed");

    return 0;
}

int bw_common_init(struct bw_params * bw, param_list pl, int * p_argc, char *** p_argv)
{
    bw_common_init_defaults(bw);
    return bw_common_init_shared(bw, pl, p_argc, p_argv);
}

int bw_common_clear(struct bw_params * bw)
{
    mpz_clear(bw->p);
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
        "\tcheckpoints=<bool>\tsave checkpoints.\n"
        "\tmatrix=<filename>\tset matrix\n"
        "\tinterval=<int>\tset checking bw->interval\n"
        "\tseed=<int>\tseed value for picking random numbers\n"
        "\tys=<int>..<int>\tcoordinate range for krylov/mksol task\n"
        ;
    return t;
}
