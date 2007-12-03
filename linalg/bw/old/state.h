#ifndef STATE_H_
#define STATE_H_

#ifdef	__cplusplus
extern "C" {
#endif

extern void commit_more_values(bw_vector_block);
extern void state_checkpoint(bw_vector_block, int);
extern int recover_state(bw_vector_block, struct mksol_info_block *);

#ifdef	__cplusplus
}
#endif

#endif /* STATE_H_ */
