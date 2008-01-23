#ifndef RIGHT_ACTION_H_
#define RIGHT_ACTION_H_

#ifdef	__cplusplus
extern "C" {
#endif

/* structure containing the relevant function to the matrices we update
 * in the course of a gaussian elimination work */

typedef void (*perm_applier_t)  (void *, unsigned int *);
typedef void (*col_exchanger_t) (void *, int, int);
typedef void (*transvecter_t)   (void *, int, int, int, bw_scalar);
typedef void (*x_multiplier_t)  (void *, int);

struct right_matrix_action {
	perm_applier_t	apply_perm;
	col_exchanger_t	exchange_cols;
	transvecter_t	transvection;
	x_multiplier_t	x_multiply;
	void * data;
};

#ifdef	__cplusplus
}
#endif

#endif	/* RIGHT_ACTION_H_ */
