#ifndef BW_COMMON_H_
#define BW_COMMON_H_

/* This file gathers the parameters that are relevant to all BW
 * programs -- most programs parse the corresponding values -- for some,
 * the values are unused, which is not a lot to worry about.
 */

#include <gmp.h>
#include <limits.h>
#include "bwc_config.h"
#include "params.h"

/* relevant to secure.c as well ; really, this is almost a hard-coded
 * constant.  Don't imagine it's easily tunable. If one really insists,
 * perhaps a u64k/secure or u64n/secure or u128/secure would work, but
 * this has not been checked. */
#define NCHECKS_CHECK_VECTOR    64

struct bw_params {

    /* m,n blocking factors as in the textbook BW description */
    int m,n;

    /* modulus -- has to be 2 for now. Depending on abase and some stuff here
     * and there, it could be expanded. */
    mpz_t p;

    /* _at the moment_ this corresponds to the checking & checkpointing
     * interval. Quite clearly, the two could be separated.
     */
    int interval;

    /* defined, but unused */
    int verbose;

    /* Whether the current job/thread may print to stdout */
    int can_print;

    /* This indicates the starting iteration -- only for krylov and mksol */
    int start;
    int end;

    /* for a given krylov/mksol task, indicates the coordinate range in
     * [0..n[ that is relevant for us. ys[1]-ys[0] defines the ``local''
     * blocking factor n'.
     */
    int ys[2];

    /* Number of coordinates in the X vectors */
    unsigned int nx;

    /* dir is a boolean flag equal to 1 if we are looking for the right
     * nullspace of the matrix. In matmul_top speak, it indicates where the
     * source vector is.
     */
    int dir;

    /* How many mpi jobs are defined. First array member indicates the number
     * of horizontal stripes in the mpi job grid, second is vertical
     */
    int mpi_split[2];

    /* Similar, but for threads. Not exclusive of the above.
    */
    int thr_split[2];

    /* Only prep and lingen are not deterministic. */
    int seed;

    /* Save checkpoints or not */
    int checkpoints;
};

extern struct bw_params bw[1];

#ifdef __cplusplus
extern "C" {
#endif

extern int bw_common_init(struct bw_params * bw, param_list pl, int * p_argc, char *** p_argv);
extern int bw_common_clear(struct bw_params * bw);
extern const char * bw_common_usage_string();

/* Some stuff which is shared as well by bw-common-mpi.h, but not within
 * the public interface nevertheless */
extern int bw_common_init_shared(struct bw_params * bw, param_list pl, int * p_argc, char *** p_argv);
extern int bw_common_init_defaults(struct bw_params * bw);

#ifdef __cplusplus
}
#endif

#endif	/* BW_COMMON_H_ */
