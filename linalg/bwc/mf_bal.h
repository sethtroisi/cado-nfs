#ifndef MF_BAL_H_
#define MF_BAL_H_

#ifdef __cplusplus
extern "C" {
#endif

struct mf_bal_args {
    const char * rwfile;
    const char * cwfile;
    const char * mfile;
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

#ifdef __cplusplus
}
#endif

#endif	/* MF_BAL_H_ */
