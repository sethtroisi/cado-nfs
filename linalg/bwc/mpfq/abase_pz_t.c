/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include "abase_pz_t.h"

/* Active handler: simd_pz */
/* Automatically generated code  */
/* Active handler: pz */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_pz, tag=pz, }, ],
   fieldtype=prime,
   n=mpz_size(k->p),
   nn=(2*mpz_size(k->p) + 1),
   tag=pz,
   type=plain,
   vbase_stuff={
    choose_byfeatures=<code>,
    families=[
     [ u64k1, u64k2, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_1, tag=p_1, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_2, tag=p_2, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_3, tag=p_3, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_4, tag=p_4, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_8, tag=p_8, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_pz, tag=pz, }, ],
     ],
    member_templates_restrict={
     p_1=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_1, tag=p_1, }, ],
     p_2=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_2, tag=p_2, }, ],
     p_3=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_3, tag=p_3, }, ],
     p_4=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_4, tag=p_4, }, ],
     p_8=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_8, tag=p_8, }, ],
     pz=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_pz, tag=pz, }, ],
     u64k1=[ u64k1, u64k2, ],
     u64k2=[ u64k1, u64k2, ],
     },
    vc:includes=[ <stdarg.h>, ],
    },
   virtual_base={
    filebase=abase_vbase,
    global_prefix=abase_,
    name=abase_vbase,
    substitutions=[
     [ (?^:abase_pz_elt \*), void *, ],
     [ (?^:abase_pz_src_elt\b), const void *, ],
     [ (?^:abase_pz_elt\b), void *, ],
     [ (?^:abase_pz_dst_elt\b), void *, ],
     [ (?^:abase_pz_elt_ur \*), void *, ],
     [ (?^:abase_pz_src_elt_ur\b), const void *, ],
     [ (?^:abase_pz_elt_ur\b), void *, ],
     [ (?^:abase_pz_dst_elt_ur\b), void *, ],
     [ (?^:abase_pz_vec \*), void *, ],
     [ (?^:abase_pz_src_vec\b), const void *, ],
     [ (?^:abase_pz_vec\b), void *, ],
     [ (?^:abase_pz_dst_vec\b), void *, ],
     [ (?^:abase_pz_vec_ur \*), void *, ],
     [ (?^:abase_pz_src_vec_ur\b), const void *, ],
     [ (?^:abase_pz_vec_ur\b), void *, ],
     [ (?^:abase_pz_dst_vec_ur\b), void *, ],
     [ (?^:abase_pz_poly \*), void *, ],
     [ (?^:abase_pz_src_poly\b), const void *, ],
     [ (?^:abase_pz_poly\b), void *, ],
     [ (?^:abase_pz_dst_poly\b), void *, ],
     ],
    },
   vtag=pz,
   w=64,
   } */


/* Functions operating on the field structure */

/* Element allocation functions */

/* Elementary assignment functions */

/* Assignment of random values */

/* Arithmetic operations on elements */

/* Operations involving unreduced elements */

/* Comparison functions */

/* Input/output functions */

/* Vector functions */

/* Polynomial functions */

/* Functions related to SIMD operation */

/* Member templates related to SIMD operation */

/* MPI interface */

/* Object-oriented interface */
/* *simd_pz::code_for_member_template_dotprod */
void abase_pz_pz_dotprod(abase_pz_dst_field K0 MAYBE_UNUSED, abase_pz_dst_field K1 MAYBE_UNUSED, abase_pz_dst_vec xw, abase_pz_src_vec xu1, abase_pz_src_vec xu0, unsigned int n)
{
        abase_pz_elt_ur s,t;
        abase_pz_elt_ur_init(K0, &s);
        abase_pz_elt_ur_init(K0, &t);
        abase_pz_elt_ur_set_zero(K0, s);
        for(unsigned int i = 0 ; i < n ; i++) {
            abase_pz_mul_ur(K0, t, abase_pz_vec_coeff_ptr_const(K0, xu0, i), abase_pz_vec_coeff_ptr_const(K1, xu1, i));
            abase_pz_elt_ur_add(K0, s, s, t);
        }
        abase_pz_reduce(K0, abase_pz_vec_coeff_ptr(K0, xw, 0), s);
        abase_pz_elt_ur_clear(K0, &s);
        abase_pz_elt_ur_clear(K0, &t);
}

/* *simd_pz::code_for_member_template_addmul_tiny */
void abase_pz_pz_addmul_tiny(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_field L MAYBE_UNUSED, abase_pz_dst_vec w, abase_pz_src_vec u, abase_pz_dst_vec v, unsigned int n)
{
        abase_pz_elt s;
        abase_pz_init(K, &s);
        for(unsigned int i = 0 ; i < n ; i++) {
            abase_pz_mul(K, s, abase_pz_vec_coeff_ptr_const(K, u, i), abase_pz_vec_coeff_ptr_const(K, v, 0));
            abase_pz_add(K, abase_pz_vec_coeff_ptr(K, w, i), abase_pz_vec_coeff_ptr_const(K, w, i), s);
        }
        abase_pz_clear(K, &s);
}

/* *simd_pz::code_for_member_template_transpose */
void abase_pz_pz_transpose(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_field L MAYBE_UNUSED, abase_pz_dst_vec w, abase_pz_src_vec u)
{
    abase_pz_set(K, abase_pz_vec_coeff_ptr(K, w, 0), abase_pz_vec_coeff_ptr_const(K, u, 0));
}


/* vim:set ft=cpp: */
