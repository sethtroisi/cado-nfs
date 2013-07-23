#ifndef ABASE_VBASE_H_
#define ABASE_VBASE_H_

/* MPFQ generated file -- do not edit */

#include <stddef.h>
#include <stdio.h>
#include <gmp.h>
#include "select_mpi.h"
struct abase_vbase_s;
typedef struct abase_vbase_s * abase_vbase_ptr;
typedef struct abase_vbase_s const * abase_vbase_srcptr;

struct abase_vbase_tmpl_s;
typedef struct abase_vbase_tmpl_s * abase_vbase_tmpl_ptr;
typedef struct abase_vbase_tmpl_s const * abase_vbase_tmpl_srcptr;

struct abase_vbase_s {
    void * obj; /* pointer to global implementation private fields */
    void (*field_characteristic)(abase_vbase_ptr, mpz_t);
    int (*field_degree)(abase_vbase_ptr);
    void (*field_init)(abase_vbase_ptr);
    void (*field_clear)(abase_vbase_ptr);
    void (*field_specify)(abase_vbase_ptr, unsigned long, void *);
    void (*field_setopt)(abase_vbase_ptr, unsigned long, void *);
    void (*init)(abase_vbase_ptr, void *);
    void (*clear)(abase_vbase_ptr, void *);
    void (*set)(abase_vbase_ptr, void *, const void *);
    void (*set_ui)(abase_vbase_ptr, void *, unsigned long);
    void (*set_zero)(abase_vbase_ptr, void *);
    unsigned long (*get_ui)(abase_vbase_ptr, const void *);
    void (*set_mpn)(abase_vbase_ptr, void *, mp_limb_t *, size_t);
    void (*set_mpz)(abase_vbase_ptr, void *, mpz_t);
    void (*get_mpn)(abase_vbase_ptr, mp_limb_t *, const void *);
    void (*get_mpz)(abase_vbase_ptr, mpz_t, const void *);
    void (*random)(abase_vbase_ptr, void *, gmp_randstate_t);
    void (*random2)(abase_vbase_ptr, void *, gmp_randstate_t);
    void (*add)(abase_vbase_ptr, void *, const void *, const void *);
    void (*sub)(abase_vbase_ptr, void *, const void *, const void *);
    void (*neg)(abase_vbase_ptr, void *, const void *);
    void (*mul)(abase_vbase_ptr, void *, const void *, const void *);
    void (*sqr)(abase_vbase_ptr, void *, const void *);
    int (*is_sqr)(abase_vbase_ptr, const void *);
    int (*sqrt)(abase_vbase_ptr, void *, const void *);
    void (*pow)(abase_vbase_ptr, void *, const void *, unsigned long *, size_t);
    void (*frobenius)(abase_vbase_ptr, void *, const void *);
    void (*add_ui)(abase_vbase_ptr, void *, const void *, unsigned long);
    void (*sub_ui)(abase_vbase_ptr, void *, const void *, unsigned long);
    void (*mul_ui)(abase_vbase_ptr, void *, const void *, unsigned long);
    int (*inv)(abase_vbase_ptr, void *, const void *);
    void (*hadamard)(abase_vbase_ptr, void *, void *, void *, void *);
    void (*elt_ur_init)(abase_vbase_ptr, void *);
    void (*elt_ur_clear)(abase_vbase_ptr, void *);
    void (*elt_ur_set)(abase_vbase_ptr, void *, const void *);
    void (*elt_ur_set_elt)(abase_vbase_ptr, void *, const void *);
    void (*elt_ur_set_zero)(abase_vbase_ptr, void *);
    void (*elt_ur_set_ui)(abase_vbase_ptr, void *, unsigned long);
    void (*elt_ur_add)(abase_vbase_ptr, void *, const void *, const void *);
    void (*elt_ur_neg)(abase_vbase_ptr, void *, const void *);
    void (*elt_ur_sub)(abase_vbase_ptr, void *, const void *, const void *);
    void (*mul_ur)(abase_vbase_ptr, void *, const void *, const void *);
    void (*sqr_ur)(abase_vbase_ptr, void *, const void *);
    void (*reduce)(abase_vbase_ptr, void *, void *);
    void (*normalize)(abase_vbase_ptr, void *);
    void (*addmul_si_ur)(abase_vbase_ptr, void *, const void *, long);
    int (*cmp)(abase_vbase_ptr, const void *, const void *);
    int (*cmp_ui)(abase_vbase_ptr, const void *, unsigned long);
    int (*is_zero)(abase_vbase_ptr, const void *);
    void (*asprint)(abase_vbase_ptr, char * *, const void *);
    void (*fprint)(abase_vbase_ptr, FILE *, const void *);
    void (*print)(abase_vbase_ptr, const void *);
    int (*sscan)(abase_vbase_ptr, void *, const char *);
    int (*fscan)(abase_vbase_ptr, FILE *, void *);
    int (*scan)(abase_vbase_ptr, void *);
    void (*vec_init)(abase_vbase_ptr, void *, unsigned int);
    void (*vec_reinit)(abase_vbase_ptr, void *, unsigned int, unsigned int);
    void (*vec_clear)(abase_vbase_ptr, void *, unsigned int);
    void (*vec_set)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_set_zero)(abase_vbase_ptr, void *, unsigned int);
    void (*vec_setcoef)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_setcoef_ui)(abase_vbase_ptr, void *, unsigned long, unsigned int);
    void (*vec_getcoef)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_add)(abase_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_neg)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_rev)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_sub)(abase_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_scal_mul)(abase_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_conv)(abase_vbase_ptr, void *, const void *, unsigned int, const void *, unsigned int);
    void (*vec_random)(abase_vbase_ptr, void *, unsigned int, gmp_randstate_t);
    void (*vec_random2)(abase_vbase_ptr, void *, unsigned int, gmp_randstate_t);
    int (*vec_cmp)(abase_vbase_ptr, const void *, const void *, unsigned int);
    int (*vec_is_zero)(abase_vbase_ptr, const void *, unsigned int);
    void (*vec_asprint)(abase_vbase_ptr, char * *, const void *, unsigned int);
    void (*vec_fprint)(abase_vbase_ptr, FILE *, const void *, unsigned int);
    void (*vec_print)(abase_vbase_ptr, const void *, unsigned int);
    int (*vec_sscan)(abase_vbase_ptr, void *, unsigned int *, const char *);
    int (*vec_fscan)(abase_vbase_ptr, FILE *, void *, unsigned int *);
    int (*vec_scan)(abase_vbase_ptr, void *, unsigned int *);
    void (*vec_ur_init)(abase_vbase_ptr, void *, unsigned int);
    void (*vec_ur_set_zero)(abase_vbase_ptr, void *, unsigned int);
    void (*vec_ur_set_vec)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_ur_reinit)(abase_vbase_ptr, void *, unsigned int, unsigned int);
    void (*vec_ur_clear)(abase_vbase_ptr, void *, unsigned int);
    void (*vec_ur_set)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_ur_setcoef)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_ur_getcoef)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_ur_add)(abase_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_ur_sub)(abase_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_ur_neg)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_ur_rev)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_scal_mul_ur)(abase_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_conv_ur)(abase_vbase_ptr, void *, const void *, unsigned int, const void *, unsigned int);
    void (*vec_reduce)(abase_vbase_ptr, void *, void *, unsigned int);
    ptrdiff_t (*vec_elt_stride)(abase_vbase_ptr, int);
    void (*poly_init)(abase_vbase_ptr, void *, unsigned int);
    void (*poly_clear)(abase_vbase_ptr, void *);
    void (*poly_set)(abase_vbase_ptr, void *, const void *);
    void (*poly_setmonic)(abase_vbase_ptr, void *, const void *);
    void (*poly_setcoef)(abase_vbase_ptr, void *, const void *, unsigned int);
    void (*poly_setcoef_ui)(abase_vbase_ptr, void *, unsigned long, unsigned int);
    void (*poly_getcoef)(abase_vbase_ptr, void *, const void *, unsigned int);
    int (*poly_deg)(abase_vbase_ptr, const void *);
    void (*poly_add)(abase_vbase_ptr, void *, const void *, const void *);
    void (*poly_sub)(abase_vbase_ptr, void *, const void *, const void *);
    void (*poly_add_ui)(abase_vbase_ptr, void *, const void *, unsigned long);
    void (*poly_sub_ui)(abase_vbase_ptr, void *, const void *, unsigned long);
    void (*poly_neg)(abase_vbase_ptr, void *, const void *);
    void (*poly_scal_mul)(abase_vbase_ptr, void *, const void *, const void *);
    void (*poly_mul)(abase_vbase_ptr, void *, const void *, const void *);
    void (*poly_divmod)(abase_vbase_ptr, void *, void *, const void *, const void *);
    void (*poly_precomp_mod)(abase_vbase_ptr, void *, const void *);
    void (*poly_mod_pre)(abase_vbase_ptr, void *, const void *, const void *, const void *);
    void (*poly_gcd)(abase_vbase_ptr, void *, const void *, const void *);
    void (*poly_xgcd)(abase_vbase_ptr, void *, void *, void *, const void *, const void *);
    void (*poly_random)(abase_vbase_ptr, void *, unsigned int, gmp_randstate_t);
    void (*poly_random2)(abase_vbase_ptr, void *, unsigned int, gmp_randstate_t);
    int (*poly_cmp)(abase_vbase_ptr, const void *, const void *);
    void (*poly_asprint)(abase_vbase_ptr, char * *, const void *);
    void (*poly_fprint)(abase_vbase_ptr, FILE *, const void *);
    void (*poly_print)(abase_vbase_ptr, const void *);
    int (*poly_sscan)(abase_vbase_ptr, void *, const char *);
    int (*poly_fscan)(abase_vbase_ptr, FILE *, void *);
    int (*poly_scan)(abase_vbase_ptr, void *);
    int (*groupsize)(abase_vbase_ptr);
    int (*offset)(abase_vbase_ptr, int);
    int (*stride)(abase_vbase_ptr);
    void (*set_ui_at)(abase_vbase_ptr, void *, int, unsigned long);
    void (*set_ui_all)(abase_vbase_ptr, void *, unsigned long);
    void (*elt_ur_set_ui_at)(abase_vbase_ptr, void *, int, unsigned long);
    void (*elt_ur_set_ui_all)(abase_vbase_ptr, void *, unsigned long);
    void (*dotprod)(abase_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*mul_constant_ui)(abase_vbase_ptr, void *, const void *, unsigned long);
    void (*mpi_ops_init)(abase_vbase_ptr);
    MPI_Datatype (*mpi_datatype)(abase_vbase_ptr);
    MPI_Datatype (*mpi_datatype_ur)(abase_vbase_ptr);
    MPI_Op (*mpi_addition_op)(abase_vbase_ptr);
    MPI_Op (*mpi_addition_op_ur)(abase_vbase_ptr);
    void (*mpi_ops_clear)(abase_vbase_ptr);
    const char * (*oo_impl_name)(abase_vbase_ptr);
    void (*oo_field_init)(abase_vbase_ptr);
    void (*oo_field_clear)(abase_vbase_ptr);
};

struct abase_vbase_tmpl_s {
    void (*dotprod)(abase_vbase_ptr, abase_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*addmul_tiny)(abase_vbase_ptr, abase_vbase_ptr, void *, const void *, void *, unsigned int);
    void (*transpose)(abase_vbase_ptr, abase_vbase_ptr, void *, const void *);
};
typedef struct abase_vbase_s abase_vbase[1];
typedef struct abase_vbase_tmpl_s abase_vbase_tmpl[1];

void abase_vbase_oo_field_init_byfeatures(abase_vbase_ptr, ...);
void abase_vbase_oo_init_templates(abase_vbase_tmpl_ptr, abase_vbase_ptr, abase_vbase_ptr);

#endif  /* ABASE_VBASE_H_ */

/* vim:set ft=cpp: */
