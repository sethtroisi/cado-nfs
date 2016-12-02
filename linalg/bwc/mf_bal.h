#ifndef MF_BAL_H_
#define MF_BAL_H_

#ifdef __cplusplus
extern "C" {
#endif

struct mf_bal_args {
    const char * rwfile;
    const char * cwfile;
    /* mf_prepare_matrix_u32 sets mfile with strdup, and there is no
     * associated dtor... So we must allow inspection, otherwise we'll
     * leak that */
    char * mfile;
    const char * bfile;
    int quiet;
    int nh;
    int nv;
    int withcoeffs;
    int rectangular;
    int skip_decorrelating_permutation;
    enum { MF_BAL_PERM_YES, MF_BAL_PERM_NO, MF_BAL_PERM_AUTO } do_perm[2];
};

void mf_bal(struct mf_bal_args * mba);
void mf_bal_adjust_from_option_string(struct mf_bal_args * mba, const char * opts);

#ifdef __cplusplus
}
#endif

#endif	/* MF_BAL_H_ */
