#ifndef MPFQ_VBASE_H_
#define MPFQ_VBASE_H_

/* MPFQ generated file -- do not edit */

#include <stddef.h>
#include <stdio.h>
#include <gmp.h>
#include "select_mpi.h"
struct mpfq_vbase_s;
typedef struct mpfq_vbase_s * mpfq_vbase_ptr;
typedef struct mpfq_vbase_s const * mpfq_vbase_srcptr;

struct mpfq_vbase_tmpl_s;
typedef struct mpfq_vbase_tmpl_s * mpfq_vbase_tmpl_ptr;
typedef struct mpfq_vbase_tmpl_s const * mpfq_vbase_tmpl_srcptr;

struct mpfq_vbase_s {
    void * obj; /* pointer to global implementation private fields */
    const char * (*impl_name)();
    unsigned long (*impl_max_characteristic_bits)();
    unsigned long (*impl_max_degree)();
    void (*field_characteristic)(mpfq_vbase_ptr, mpz_t);
    unsigned long (*field_characteristic_bits)(mpfq_vbase_ptr);
    int (*field_degree)(mpfq_vbase_ptr);
    void (*field_init)(mpfq_vbase_ptr);
    void (*field_clear)(mpfq_vbase_ptr);
    void (*field_specify)(mpfq_vbase_ptr, unsigned long, void *);
    void (*field_setopt)(mpfq_vbase_ptr, unsigned long, void *);
    void (*init)(mpfq_vbase_ptr, void *);
    void (*clear)(mpfq_vbase_ptr, void *);
    void (*set)(mpfq_vbase_ptr, void *, const void *);
    void (*set_ui)(mpfq_vbase_ptr, void *, unsigned long);
    void (*set_zero)(mpfq_vbase_ptr, void *);
    unsigned long (*get_ui)(mpfq_vbase_ptr, const void *);
    void (*set_mpn)(mpfq_vbase_ptr, void *, mp_limb_t *, size_t);
    void (*set_mpz)(mpfq_vbase_ptr, void *, mpz_t);
    void (*get_mpn)(mpfq_vbase_ptr, mp_limb_t *, const void *);
    void (*get_mpz)(mpfq_vbase_ptr, mpz_t, const void *);
    void (*random)(mpfq_vbase_ptr, void *, gmp_randstate_t);
    void (*random2)(mpfq_vbase_ptr, void *, gmp_randstate_t);
    void (*add)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*sub)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*neg)(mpfq_vbase_ptr, void *, const void *);
    void (*mul)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*sqr)(mpfq_vbase_ptr, void *, const void *);
    int (*is_sqr)(mpfq_vbase_ptr, const void *);
    int (*sqrt)(mpfq_vbase_ptr, void *, const void *);
    void (*pow)(mpfq_vbase_ptr, void *, const void *, unsigned long *, size_t);
    void (*powz)(mpfq_vbase_ptr, void *, const void *, mpz_srcptr);
    void (*frobenius)(mpfq_vbase_ptr, void *, const void *);
    void (*add_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*sub_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*mul_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    int (*inv)(mpfq_vbase_ptr, void *, const void *);
    void (*hadamard)(mpfq_vbase_ptr, void *, void *, void *, void *);
    void (*elt_ur_init)(mpfq_vbase_ptr, void *);
    void (*elt_ur_clear)(mpfq_vbase_ptr, void *);
    void (*elt_ur_set)(mpfq_vbase_ptr, void *, const void *);
    void (*elt_ur_set_elt)(mpfq_vbase_ptr, void *, const void *);
    void (*elt_ur_set_zero)(mpfq_vbase_ptr, void *);
    void (*elt_ur_set_ui)(mpfq_vbase_ptr, void *, unsigned long);
    void (*elt_ur_add)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*elt_ur_neg)(mpfq_vbase_ptr, void *, const void *);
    void (*elt_ur_sub)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*mul_ur)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*sqr_ur)(mpfq_vbase_ptr, void *, const void *);
    void (*reduce)(mpfq_vbase_ptr, void *, void *);
    void (*normalize)(mpfq_vbase_ptr, void *);
    void (*addmul_si_ur)(mpfq_vbase_ptr, void *, const void *, long);
    int (*cmp)(mpfq_vbase_ptr, const void *, const void *);
    int (*cmp_ui)(mpfq_vbase_ptr, const void *, unsigned long);
    int (*is_zero)(mpfq_vbase_ptr, const void *);
    int (*asprint)(mpfq_vbase_ptr, char * *, const void *);
    int (*fprint)(mpfq_vbase_ptr, FILE *, const void *);
    int (*print)(mpfq_vbase_ptr, const void *);
    int (*sscan)(mpfq_vbase_ptr, void *, const char *);
    int (*fscan)(mpfq_vbase_ptr, FILE *, void *);
    int (*scan)(mpfq_vbase_ptr, void *);
    void (*vec_init)(mpfq_vbase_ptr, void *, unsigned int);
    void (*vec_reinit)(mpfq_vbase_ptr, void *, unsigned int, unsigned int);
    void (*vec_clear)(mpfq_vbase_ptr, void *, unsigned int);
    void (*vec_set)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_set_zero)(mpfq_vbase_ptr, void *, unsigned int);
    void (*vec_setcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_setcoeff_ui)(mpfq_vbase_ptr, void *, unsigned long, unsigned int);
    void (*vec_getcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_add)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_neg)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_rev)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_sub)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_scal_mul)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_conv)(mpfq_vbase_ptr, void *, const void *, unsigned int, const void *, unsigned int);
    void (*vec_random)(mpfq_vbase_ptr, void *, unsigned int, gmp_randstate_t);
    void (*vec_random2)(mpfq_vbase_ptr, void *, unsigned int, gmp_randstate_t);
    int (*vec_cmp)(mpfq_vbase_ptr, const void *, const void *, unsigned int);
    int (*vec_is_zero)(mpfq_vbase_ptr, const void *, unsigned int);
    void * (*vec_subvec)(mpfq_vbase_ptr, void *, int);
    const void * (*vec_subvec_const)(mpfq_vbase_ptr, const void *, int);
    void * (*vec_coeff_ptr)(mpfq_vbase_ptr, void *, int);
    const void * (*vec_coeff_ptr_const)(mpfq_vbase_ptr, const void *, int);
    int (*vec_asprint)(mpfq_vbase_ptr, char * *, const void *, unsigned int);
    int (*vec_fprint)(mpfq_vbase_ptr, FILE *, const void *, unsigned int);
    int (*vec_print)(mpfq_vbase_ptr, const void *, unsigned int);
    int (*vec_sscan)(mpfq_vbase_ptr, void *, unsigned int *, const char *);
    int (*vec_fscan)(mpfq_vbase_ptr, FILE *, void *, unsigned int *);
    int (*vec_scan)(mpfq_vbase_ptr, void *, unsigned int *);
    void (*vec_ur_init)(mpfq_vbase_ptr, void *, unsigned int);
    void (*vec_ur_set_zero)(mpfq_vbase_ptr, void *, unsigned int);
    void (*vec_ur_set_vec)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_ur_reinit)(mpfq_vbase_ptr, void *, unsigned int, unsigned int);
    void (*vec_ur_clear)(mpfq_vbase_ptr, void *, unsigned int);
    void (*vec_ur_set)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_ur_setcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_ur_getcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_ur_add)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_ur_sub)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_ur_neg)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_ur_rev)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*vec_scal_mul_ur)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*vec_conv_ur)(mpfq_vbase_ptr, void *, const void *, unsigned int, const void *, unsigned int);
    void (*vec_reduce)(mpfq_vbase_ptr, void *, void *, unsigned int);
    void * (*vec_ur_subvec)(mpfq_vbase_ptr, void *, int);
    const void * (*vec_ur_subvec_const)(mpfq_vbase_ptr, const void *, int);
    void * (*vec_ur_coeff_ptr)(mpfq_vbase_ptr, void *, int);
    const void * (*vec_ur_coeff_ptr_const)(mpfq_vbase_ptr, const void *, int);
    ptrdiff_t (*vec_elt_stride)(mpfq_vbase_ptr, int);
    ptrdiff_t (*vec_ur_elt_stride)(mpfq_vbase_ptr, int);
    void (*poly_init)(mpfq_vbase_ptr, void *, unsigned int);
    void (*poly_clear)(mpfq_vbase_ptr, void *);
    void (*poly_set)(mpfq_vbase_ptr, void *, const void *);
    void (*poly_setmonic)(mpfq_vbase_ptr, void *, const void *);
    void (*poly_setcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*poly_setcoeff_ui)(mpfq_vbase_ptr, void *, unsigned long, unsigned int);
    void (*poly_getcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    int (*poly_deg)(mpfq_vbase_ptr, const void *);
    void (*poly_add)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*poly_sub)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*poly_add_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*poly_sub_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*poly_neg)(mpfq_vbase_ptr, void *, const void *);
    void (*poly_scal_mul)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*poly_mul)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*poly_divmod)(mpfq_vbase_ptr, void *, void *, const void *, const void *);
    void (*poly_precomp_mod)(mpfq_vbase_ptr, void *, const void *);
    void (*poly_mod_pre)(mpfq_vbase_ptr, void *, const void *, const void *, const void *);
    void (*poly_gcd)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*poly_xgcd)(mpfq_vbase_ptr, void *, void *, void *, const void *, const void *);
    void (*poly_random)(mpfq_vbase_ptr, void *, unsigned int, gmp_randstate_t);
    void (*poly_random2)(mpfq_vbase_ptr, void *, unsigned int, gmp_randstate_t);
    int (*poly_cmp)(mpfq_vbase_ptr, const void *, const void *);
    int (*poly_asprint)(mpfq_vbase_ptr, char * *, const void *);
    int (*poly_fprint)(mpfq_vbase_ptr, FILE *, const void *);
    int (*poly_print)(mpfq_vbase_ptr, const void *);
    int (*poly_sscan)(mpfq_vbase_ptr, void *, const char *);
    int (*poly_fscan)(mpfq_vbase_ptr, FILE *, void *);
    int (*poly_scan)(mpfq_vbase_ptr, void *);
    int (*groupsize)(mpfq_vbase_ptr);
    int (*offset)(mpfq_vbase_ptr, int);
    int (*stride)(mpfq_vbase_ptr);
    void (*set_ui_at)(mpfq_vbase_ptr, void *, int, unsigned long);
    void (*set_ui_all)(mpfq_vbase_ptr, void *, unsigned long);
    void (*elt_ur_set_ui_at)(mpfq_vbase_ptr, void *, int, unsigned long);
    void (*elt_ur_set_ui_all)(mpfq_vbase_ptr, void *, unsigned long);
    void (*dotprod)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*mul_constant_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*mpi_ops_init)(mpfq_vbase_ptr);
    MPI_Datatype (*mpi_datatype)(mpfq_vbase_ptr);
    MPI_Datatype (*mpi_datatype_ur)(mpfq_vbase_ptr);
    MPI_Op (*mpi_addition_op)(mpfq_vbase_ptr);
    MPI_Op (*mpi_addition_op_ur)(mpfq_vbase_ptr);
    void (*mpi_ops_clear)(mpfq_vbase_ptr);
    void (*oo_field_init)(mpfq_vbase_ptr);
    void (*oo_field_clear)(mpfq_vbase_ptr);
};

struct mpfq_vbase_tmpl_s {
    void (*dotprod)(mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*addmul_tiny)(mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, void *, unsigned int);
    void (*transpose)(mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *);
};
typedef struct mpfq_vbase_s mpfq_vbase[1];
typedef struct mpfq_vbase_tmpl_s mpfq_vbase_tmpl[1];

void mpfq_vbase_oo_field_init_byfeatures(mpfq_vbase_ptr, ...);
void mpfq_vbase_oo_init_templates(mpfq_vbase_tmpl_ptr, mpfq_vbase_ptr, mpfq_vbase_ptr);

#endif  /* MPFQ_VBASE_H_ */

/* vim:set ft=cpp: */
