/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include "mpfq_p_1_t.h"

/* Active handler: simd_gfp */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::gfp::elt */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_1, tag=p_1, }, ],
   fieldtype=prime,
   n=1,
   nn=3,
   opthw=,
   tag=p_1,
   type=plain,
   vbase_stuff={
    choose_byfeatures=<code>,
    families=[
     [ u64k1, u64k2, u64k4, ],
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
     u64k1=[ u64k1, u64k2, u64k4, ],
     u64k2=[ u64k1, u64k2, u64k4, ],
     u64k4=[ u64k1, u64k2, u64k4, ],
     },
    vc:includes=[ <stdarg.h>, ],
    },
   virtual_base={
    filebase=mpfq_vbase,
    global_prefix=mpfq_,
    name=mpfq_vbase,
    substitutions=[
     [ (?^:mpfq_p_1_elt \*), void *, ],
     [ (?^:mpfq_p_1_src_elt\b), const void *, ],
     [ (?^:mpfq_p_1_elt\b), void *, ],
     [ (?^:mpfq_p_1_dst_elt\b), void *, ],
     [ (?^:mpfq_p_1_elt_ur \*), void *, ],
     [ (?^:mpfq_p_1_src_elt_ur\b), const void *, ],
     [ (?^:mpfq_p_1_elt_ur\b), void *, ],
     [ (?^:mpfq_p_1_dst_elt_ur\b), void *, ],
     [ (?^:mpfq_p_1_vec \*), void *, ],
     [ (?^:mpfq_p_1_src_vec\b), const void *, ],
     [ (?^:mpfq_p_1_vec\b), void *, ],
     [ (?^:mpfq_p_1_dst_vec\b), void *, ],
     [ (?^:mpfq_p_1_vec_ur \*), void *, ],
     [ (?^:mpfq_p_1_src_vec_ur\b), const void *, ],
     [ (?^:mpfq_p_1_vec_ur\b), void *, ],
     [ (?^:mpfq_p_1_dst_vec_ur\b), void *, ],
     [ (?^:mpfq_p_1_poly \*), void *, ],
     [ (?^:mpfq_p_1_src_poly\b), const void *, ],
     [ (?^:mpfq_p_1_poly\b), void *, ],
     [ (?^:mpfq_p_1_dst_poly\b), void *, ],
     ],
    },
   vtag=p_1,
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
/* *simd_gfp::code_for_member_template_dotprod */
void mpfq_p_1_p_1_dotprod(mpfq_p_1_dst_field K0 MAYBE_UNUSED, mpfq_p_1_dst_field K1 MAYBE_UNUSED, mpfq_p_1_dst_vec xw, mpfq_p_1_src_vec xu1, mpfq_p_1_src_vec xu0, unsigned int n)
{
        mpfq_p_1_elt_ur s,t;
        mpfq_p_1_elt_ur_init(K0, &s);
        mpfq_p_1_elt_ur_init(K0, &t);
        mpfq_p_1_elt_ur_set_zero(K0, s);
        for(unsigned int i = 0 ; i < n ; i++) {
            mpfq_p_1_mul_ur(K0, t, xu0[i], xu1[i]);
            mpfq_p_1_elt_ur_add(K0, s, s, t);
        }
        mpfq_p_1_reduce(K0, xw[0], s);
        mpfq_p_1_elt_ur_clear(K0, &s);
        mpfq_p_1_elt_ur_clear(K0, &t);
}

/* *simd_gfp::code_for_member_template_addmul_tiny */
void mpfq_p_1_p_1_addmul_tiny(mpfq_p_1_dst_field K MAYBE_UNUSED, mpfq_p_1_dst_field L MAYBE_UNUSED, mpfq_p_1_dst_vec w, mpfq_p_1_src_vec u, mpfq_p_1_dst_vec v, unsigned int n)
{
        mpfq_p_1_elt s;
        mpfq_p_1_init(K, &s);
        for(unsigned int i = 0 ; i < n ; i++) {
            mpfq_p_1_mul(K, s, u[i], v[0]);
            mpfq_p_1_add(K, w[i], w[i], s);
        }
        mpfq_p_1_clear(K, &s);
}

/* *simd_gfp::code_for_member_template_transpose */
void mpfq_p_1_p_1_transpose(mpfq_p_1_dst_field K MAYBE_UNUSED, mpfq_p_1_dst_field L MAYBE_UNUSED, mpfq_p_1_dst_vec w, mpfq_p_1_src_vec u)
{
    mpfq_p_1_set(K, w[0], u[0]);
}


/* vim:set ft=cpp: */
