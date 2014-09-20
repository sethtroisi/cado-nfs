#include "cado.h"

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include "bwc_config.h"
#include "cado_config.h"
#include "bw-common.h"
#include "params.h"
#include "filenames.h"
#include "portability.h"
#include "utils.h"

struct bw_params bw[1];

const char * dirtext[] = { "left", "right" };

/* Has to be defined by the program */
extern void usage();

int bw_common_init_defaults(struct bw_params * bw)
{
    /*** defaults ***/
    memset(bw, 0, sizeof(*bw));
    bw->interval = 1000;
    bw->can_print = 1;
    bw->ys[0] = bw->ys[1] = -1;
    bw->dir = 1;
    bw->mpi_split[0] = bw->mpi_split[1] = 1;
    bw->thr_split[0] = bw->thr_split[1] = 1;

    bw->checkpoints = 1;

    return 0;
}

typedef int (*sortfunc_t) (const void*, const void*);

static int intcmp(const int * a, const int * b)
{
    return (*a>*b)-(*b>*a);
}


int bw_common_init_shared(struct bw_params * bw, param_list pl, int * p_argc, char *** p_argv)
{
    bw->original_argc = *p_argc;
    bw->original_argv = *p_argv;
    bw->wct_base = wct_seconds();

    (*p_argv)++, (*p_argc)--;
    param_list_configure_switch(pl, "-v", &bw->verbose);
    for( ; (*p_argc) ; ) {
        if (param_list_update_cmdline(pl, p_argc, p_argv)) { continue; }
        if (strcmp((*p_argv)[0],"--") == 0) {
            (*p_argv)++, (*p_argc)--;
            break;
        }
        fprintf(stderr, "Unhandled parameter %s\n", (*p_argv)[0]);
        usage();
    }

    verbose_set_enabled_flags(pl);

    if (bw->can_print) {
        param_list_print_command_line(stderr, pl);
        param_list_print_command_line(stdout, pl);
    }

    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "wdir")) != NULL) {
        /* We now do mkdir -p beforehand on all jobs. Note that at the
         * point where the current function is being called, we're not
         * multithreaded yet */
        mkdir_with_parents(tmp, 1);
        if (chdir(tmp) < 0) {
            fprintf(stderr, "chdir(%s): %s\n", tmp, strerror(errno));
            exit(1);
        }
    }

    const char * cfg;

    if ((cfg = param_list_lookup_string(pl, "cfg"))) {
        param_list_read_file(pl, cfg);
    }


    param_list_parse_intxint(pl, "mpi", bw->mpi_split);
    param_list_parse_intxint(pl, "thr", bw->thr_split);
    param_list_parse_int(pl, "seed", &bw->seed);
    param_list_parse_int(pl, "interval", &bw->interval);
    param_list_parse_int_and_int(pl, "ys", bw->ys, "..");
    param_list_parse_int(pl, "start", &bw->start);
    param_list_parse_int(pl, "end", &bw->end);
    param_list_parse_int(pl, "checkpoints", &bw->checkpoints);
    param_list_parse_int(pl, "skip_online_checks", &bw->skip_online_checks);
    param_list_parse_int(pl, "keep_rolling_checkpoints", &bw->keep_rolling_checkpoints);
    param_list_parse_int(pl, "keep_checkpoints_younger_than", &bw->keep_checkpoints_younger_than);
    param_list_parse_int(pl, "checkpoint_precious", &bw->checkpoint_precious);

    int yes_i_insist = 0;
    param_list_parse_int(pl, "yes_i_insist", &yes_i_insist);

    if (bw->skip_online_checks && bw->keep_rolling_checkpoints) {
        fprintf(stderr, "The combination of skip_online_checks and keep_rolling_checkpointsis a dangerous match.");
        if (!yes_i_insist) {
            printf("\n");
            exit(1);
        }
        fprintf(stderr, " Proceeding anyway\n");
    }

    param_list_lookup_string(pl, "matrix");
    param_list_lookup_string(pl, "balancing");


    mpz_init_set_ui(bw->p, 2);
    param_list_parse_mpz(pl, "prime", bw->p);
    int nullspace_forced = 0;

    if ((tmp = param_list_lookup_string(pl, "nullspace")) != NULL) {
        char * tmp_l = strdup(tmp);
        for(unsigned int i = 0 ; i < strlen(tmp_l) ; i++) {
            char c = tmp[i];
            char cl = tolower(c);
            tmp_l[i] = cl;
            nullspace_forced |= c != cl;
        }
        if (strcmp(tmp_l, dirtext[0]) == 0) {
            bw->dir = 0;
        } else if (strcmp(tmp, dirtext[1]) == 0) {
            bw->dir = 1;
        } else {
            fprintf(stderr, "Parameter nullspace may only be %s|%s\n",
                    dirtext[0], dirtext[1]);
            exit(1);
        }
        free(tmp_l);
    } else {
        /* Default is right nullspace for p>2, and left for p==2 */
        bw->dir = mpz_cmp_ui(bw->p, 2) != 0;
        param_list_add_key(pl, "nullspace", dirtext[bw->dir], PARAMETER_FROM_FILE);
    }

    if ((mpz_cmp_ui(bw->p, 2) == 0) != (bw->dir == 0)) {
        if (mpz_cmp_ui(bw->p, 2) == 0) {
            fprintf(stderr, "p==2 seems appropriate for factoring. Yet, the nullspace parameter has been passed as nullspace=right. This looks odd.\n");
            if (nullspace_forced) {
                fprintf(stderr, "Proceeding anyway (uppercase nullspace argument)\n");
            } else {
                fprintf(stderr, "Aborting. Pass nullspace=LEFT if this is really intended.\n");
                exit(1);
            }
        } else {
            fprintf(stderr, "p>2 seems appropriate for discrete logarithm. Yet, the nullspace parameter has been passed as nullspace=left. This looks odd.\n");
            if (nullspace_forced) {
                fprintf(stderr, "Proceeding anyway (uppercase nullspace argument)\n");
            } else {
                fprintf(stderr, "Aborting. Pass nullspace=RIGHT if this is really intended.\n");
                exit(1);
            }
        }
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
    if (!okm || !okn) {
        fprintf(stderr, "parameter m and/or n is missing\n");
        usage();
    }
    bw->nsolvecs = bw->n;
    param_list_parse_int(pl, "nsolvecs", &bw->nsolvecs);

    param_list_parse_int(pl, "nrhs", &bw->nrhs);

    bw->number_of_check_stops = param_list_parse_int_list(pl, "check_stops", bw->check_stops, MAX_NUMBER_OF_CHECK_STOPS, ",");
    int interval_already_in_check_stops = 0;
    for(int i = 0 ; i < bw->number_of_check_stops ; i++) {
        if (bw->check_stops[i] == bw->interval) {
            interval_already_in_check_stops = 1;
            break;
        }
    }
    if (!interval_already_in_check_stops) {
        ASSERT_ALWAYS(bw->number_of_check_stops < MAX_NUMBER_OF_CHECK_STOPS - 1);
        bw->check_stops[bw->number_of_check_stops++] = bw->interval;
    }
    qsort(bw->check_stops, bw->number_of_check_stops, sizeof(int), (sortfunc_t) &intcmp);



    if (bw->verbose && bw->can_print)
        param_list_display (pl, stderr);

    /* Force lookup of parameters that are used late in the process. This
     * has the effect of eliminating possible warnings. */
    param_list_lookup_string(pl, "mm_impl");

    param_list_lookup_string(pl, "l1_cache_size");
    param_list_lookup_string(pl, "cache_line_size");

    param_list_lookup_string(pl, "mm_threaded_nthreads");
    param_list_lookup_string(pl, "mm_threaded_sgroup_size");
    param_list_lookup_string(pl, "mm_threaded_offset1");
    param_list_lookup_string(pl, "mm_threaded_offset2");
    param_list_lookup_string(pl, "mm_threaded_offset3");
    param_list_lookup_string(pl, "mm_threaded_densify_tolerance");
    param_list_lookup_string(pl, "mm_store_transposed");
    param_list_lookup_string(pl, "interleaving");

    param_list_lookup_string(pl, "rebuild_cache");
    param_list_lookup_string(pl, "sequential_cache_build");
    param_list_lookup_string(pl, "sequential_cache_load");
    param_list_lookup_string(pl, "local_cache_copy_dir");
    param_list_lookup_string(pl, "matmul_bucket_methods");
    param_list_lookup_string(pl, "sequence");   // for lingen
    param_list_lookup_string(pl, "rhs");  // for gather
    param_list_lookup_string(pl, "save_submatrices");
    param_list_lookup_string(pl, "export_cachelist");
    param_list_lookup_string(pl, "sanity_check_vector");
    param_list_lookup_string(pl, "cpubinding");
    param_list_lookup_string(pl, "only_mpi");
    param_list_lookup_string(pl, "balancing_queue_size");
    param_list_lookup_string(pl, "random_matrix");

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
