/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include "abase_p_8_t.h"

/* Active handler: simd_gfp */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::gfp::elt */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_8, tag=p_8, }, ],
   fieldtype=prime,
   n=8,
   nn=17,
   opthw=,
   tag=p_8,
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
     [ (?^:abase_p_8_elt \*), void *, ],
     [ (?^:abase_p_8_src_elt\b), const void *, ],
     [ (?^:abase_p_8_elt\b), void *, ],
     [ (?^:abase_p_8_dst_elt\b), void *, ],
     [ (?^:abase_p_8_elt_ur \*), void *, ],
     [ (?^:abase_p_8_src_elt_ur\b), const void *, ],
     [ (?^:abase_p_8_elt_ur\b), void *, ],
     [ (?^:abase_p_8_dst_elt_ur\b), void *, ],
     [ (?^:abase_p_8_vec \*), void *, ],
     [ (?^:abase_p_8_src_vec\b), const void *, ],
     [ (?^:abase_p_8_vec\b), void *, ],
     [ (?^:abase_p_8_dst_vec\b), void *, ],
     [ (?^:abase_p_8_vec_ur \*), void *, ],
     [ (?^:abase_p_8_src_vec_ur\b), const void *, ],
     [ (?^:abase_p_8_vec_ur\b), void *, ],
     [ (?^:abase_p_8_dst_vec_ur\b), void *, ],
     [ (?^:abase_p_8_poly \*), void *, ],
     [ (?^:abase_p_8_src_poly\b), const void *, ],
     [ (?^:abase_p_8_poly\b), void *, ],
     [ (?^:abase_p_8_dst_poly\b), void *, ],
     ],
    },
   vtag=p_8,
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
void abase_p_8_p_8_dotprod(abase_p_8_dst_field K0 MAYBE_UNUSED, abase_p_8_dst_field K1 MAYBE_UNUSED, abase_p_8_dst_vec xw, abase_p_8_src_vec xu1, abase_p_8_src_vec xu0, unsigned int n)
{
        abase_p_8_elt_ur s,t;
        abase_p_8_elt_ur_init(K0, &s);
        abase_p_8_elt_ur_init(K0, &t);
        abase_p_8_elt_ur_set_zero(K0, s);
        for(unsigned int i = 0 ; i < n ; i++) {
            abase_p_8_mul_ur(K0, t, xu0[i], xu1[i]);
            abase_p_8_elt_ur_add(K0, s, s, t);
        }
        abase_p_8_reduce(K0, xw[0], s);
        abase_p_8_elt_ur_clear(K0, &s);
        abase_p_8_elt_ur_clear(K0, &t);
}

/* *simd_gfp::code_for_member_template_addmul_tiny */
void abase_p_8_p_8_addmul_tiny(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_field L MAYBE_UNUSED, abase_p_8_dst_vec w, abase_p_8_src_vec u, abase_p_8_dst_vec v, unsigned int n)
{
        abase_p_8_elt s;
        abase_p_8_init(K, &s);
        for(unsigned int i = 0 ; i < n ; i++) {
            abase_p_8_mul(K, s, u[i], v[0]);
            abase_p_8_add(K, w[i], w[i], s);
        }
        abase_p_8_clear(K, &s);
}

/* *simd_gfp::code_for_member_template_transpose */
void abase_p_8_p_8_transpose(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_field L MAYBE_UNUSED, abase_p_8_dst_vec w, abase_p_8_src_vec u)
{
    abase_p_8_set(K, w[0], u[0]);
}


/* vim:set ft=cpp: */
