#ifndef TWISTING_POLYNOMIALS_H_
#define TWISTING_POLYNOMIALS_H_

#ifdef	__cplusplus
extern "C" {
#endif

struct t_poly {
	bw_bbpoly p;
	int degree;		/* this is the maximum degree of the entries
				   (actually the maximum nominal degree) */
	unsigned int alloc;	/* space allocated (must be >= degree+1) */
	unsigned int * clist;
	int * degnom;
	int * degree_table;	/* Degree of each entry (-1 for 0 polynomial) */
};

extern struct t_poly * tp_alloc(unsigned int);
extern struct t_poly * tp_comp_alloc(struct t_poly *,struct t_poly *);
extern void tp_set_ident(struct t_poly *);

extern void tp_apply_perm(struct t_poly *, unsigned int *);
extern void tp_exchange_cols(struct t_poly *, int, int);
extern void tp_transvec(struct t_poly *, int, int, int, bw_scalar);
extern void tp_x_multiply(struct t_poly *, int);
extern void tp_act_on_delta(struct t_poly *, int *);
extern void tp_print(struct t_poly *);
extern void tp_free(struct t_poly *);
extern int tp_write(FILE *, struct t_poly *);
extern struct t_poly * tp_read(FILE *, unsigned int);

#ifdef	__cplusplus
}
#endif

#endif	/* TWISTING_POLYNOMIALS_H_ */
