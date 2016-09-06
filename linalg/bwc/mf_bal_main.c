#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include "portability.h"
#include "utils.h"
#include "mf.h"
#include "balancing.h"
#include "rowset_heap.h"
#include "cheating_vec_init.h"
#include "mf_bal.h"

/* This program computes how a matrix would have to be balanced for
 * fitting a job grid of given size. This does _not_ read the matrix,
 * only the row-weight and col-weight files are read.
 */

/* TODO: Now that all files input and output by bwc are permutation
 * independent, there is room for computing this permutation on the fly
 * if the balancing argument is not provided (and only in this case).
 * Essentially we would hook onto the mmt_fill_fields_from_balancing()
 * function called from mmt_finish_init().
 */
void usage(int rc) {
    fprintf(stderr,
            "Usage: ./mf_bal [options,flags] <nh> <nv> <matrix file>\n");
    fprintf(stderr,
            "Recognized options"
                " (<option_name>=<value>, or --<option_name> <value>:\n"
            " mfile       matrix file (can also be given freeform)\n"
            " rwfile      row weight file (defaults to <mfile>.rw)\n"
            " cwfile      col weight file (defaults to <mfile>.cw)\n"
            " like        balance in a way compatible with this file\n"
            " out         output file name (defaults to stdout)\n"
            "Recognized flags:"
            " --ascii     output in ascii\n"
            " --quiet     be quiet\n"
            " --rowperm   permute rows in priority (defaults to auto)\n"
            " --colperm   permute rows in priority (defaults to auto)\n"
           );
    exit(rc);
}

int main(int argc, char * argv[])
{
    param_list pl;
    unsigned int wild =  0;
    struct mf_bal_args mba[1];
    memset(mba, 0, sizeof(struct mf_bal_args));
    mba->do_perm[0] = MF_BAL_PERM_AUTO;
    mba->do_perm[1] = MF_BAL_PERM_AUTO;
    // int display_correlation = 0;

    param_list_init(pl);
    argv++,argc--;
    param_list_configure_switch(pl, "--quiet", &mba->quiet);
    // param_list_configure_switch(pl, "--display-correlation", &display_correlation);
    param_list_configure_switch(pl, "--rectangular", &mba->rectangular);
    param_list_configure_switch(pl, "--withcoeffs", &mba->withcoeffs);

    for(;argc;) {
        char * q;
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;

        if (argv[0][0] != '-' && wild == 0 && (q = strchr(argv[0],'x')) != NULL) {
            mba->nh = atoi(argv[0]);
            mba->nv = atoi(q+1);
            wild+=2;
            argv++,argc--;
            continue;
        }

        if (argv[0][0] != '-' && wild == 0) { mba->nh = atoi(argv[0]); wild++,argv++,argc--; continue; }
        if (argv[0][0] != '-' && wild == 1) { mba->nv = atoi(argv[0]); wild++,argv++,argc--; continue; }
        if (argv[0][0] != '-' && wild == 2) {
            mba->mfile = argv[0];
            wild++;
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "unknown option %s\n", argv[0]);
        exit(1);
    }

    if (!mba->nh || !mba->nv)
        usage(1);

    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "reorder")) != NULL) {
        if (strcmp(tmp, "auto") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_AUTO;
            mba->do_perm[0] = MF_BAL_PERM_AUTO;
        } else if (strcmp(tmp, "rows") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_NO;
            mba->do_perm[0] = MF_BAL_PERM_YES;
        } else if (strcmp(tmp, "columns") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_YES;
            mba->do_perm[0] = MF_BAL_PERM_NO;
        } else if (strcmp(tmp, "rows,columns") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_YES;
            mba->do_perm[0] = MF_BAL_PERM_YES;
        } else if (strcmp(tmp, "columns,rows") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_YES;
            mba->do_perm[0] = MF_BAL_PERM_YES;
        } else if (strcmp(tmp, "both") == 0) {
            mba->do_perm[1] = MF_BAL_PERM_YES;
            mba->do_perm[0] = MF_BAL_PERM_YES;
        } else {
            fprintf(stderr, "Argument \"%s\" to the \"reorder\" parameter not understood\n"
                    "Supported values are:\n"
                    "\tauto (default)\n"
                    "\trows\n"
                    "\tcolumns\n"
                    "\tboth (equivalent forms: \"rows,columns\" or \"columns,rows\"\n",
                    tmp);
            exit(EXIT_FAILURE);
        }
    }

    if ((tmp = param_list_lookup_string(pl, "mfile")) != NULL) {
        mba->mfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "rwfile")) != NULL) {
        mba->rwfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "cwfile")) != NULL) {
        mba->cwfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "out")) != NULL) {
        mba->bfile = tmp;
    }
    param_list_parse_int(pl, "skip_decorrelating_permutation", &mba->skip_decorrelating_permutation);

    mf_bal(mba);

    /* mba holds pointers to inside the pl struct, so beware ! */
    param_list_clear(pl);

    return 0;
}
