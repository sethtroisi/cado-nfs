#ifndef SCALAR_MATRIX_H_
#define SCALAR_MATRIX_H_

#ifdef	__cplusplus
extern "C" {
#endif

struct e_coeff {
	bw_mbpoly original_p;
	bw_mbpoly p;
	unsigned int * clist;
	int degree;
};

extern struct e_coeff * ec_encapsulate(bw_mbpoly, unsigned int *, int);
extern void ec_unlink(struct e_coeff *);
extern void ec_apply_perm(struct e_coeff *, unsigned int *);
extern void ec_exchange_cols(struct e_coeff *, int, int);
extern void ec_transvec(struct e_coeff *, int, int, int, bw_scalar);
extern void ec_x_multiply(struct e_coeff *, int);
extern void ec_advance(struct e_coeff *, int);
extern void ec_park(struct e_coeff *);
extern void ec_untwist(struct e_coeff *);
extern int ec_is_twisted(struct e_coeff *);
extern void ec_print(struct e_coeff *);

#ifdef	__cplusplus
}
#endif

#endif	/* SCALAR_MATRIX_H_ */
