#ifndef SCALAR_MATRIX_H_
#define SCALAR_MATRIX_H_

#ifdef	__cplusplus
extern "C" {
#endif

struct e_coeff {
	bw_mbmat mat;
	unsigned int * clist;
};

extern void ec_apply_perm(struct e_coeff *, unsigned int *);
extern void ec_exchange_cols(struct e_coeff *, int, int);
extern void ec_transvec(struct e_coeff *, int, int, int, bw_scalar);
extern void ec_x_multiply(struct e_coeff *, int);

#ifdef	__cplusplus
}
#endif

#endif	/* SCALAR_MATRIX_H_ */
