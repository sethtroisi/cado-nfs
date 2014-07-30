#ifndef MPFQ_P_4_T_H_
#define MPFQ_P_4_T_H_

/* MPFQ generated file -- do not edit */

#include "mpfq_p_4.h"
/* Active handler: simd_gfp */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::gfp::elt */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELD_p_4, tag=p_4, }, ],
   fieldtype=prime,
   n=4,
   nn=9,
   opthw=,
   tag=p_4,
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
    filebase=mpfq_vbase,
    global_prefix=mpfq_,
    name=mpfq_vbase,
    substitutions=[
     [ (?^:mpfq_p_4_elt \*), void *, ],
     [ (?^:mpfq_p_4_src_elt\b), const void *, ],
     [ (?^:mpfq_p_4_elt\b), void *, ],
     [ (?^:mpfq_p_4_dst_elt\b), void *, ],
     [ (?^:mpfq_p_4_elt_ur \*), void *, ],
     [ (?^:mpfq_p_4_src_elt_ur\b), const void *, ],
     [ (?^:mpfq_p_4_elt_ur\b), void *, ],
     [ (?^:mpfq_p_4_dst_elt_ur\b), void *, ],
     [ (?^:mpfq_p_4_vec \*), void *, ],
     [ (?^:mpfq_p_4_src_vec\b), const void *, ],
     [ (?^:mpfq_p_4_vec\b), void *, ],
     [ (?^:mpfq_p_4_dst_vec\b), void *, ],
     [ (?^:mpfq_p_4_vec_ur \*), void *, ],
     [ (?^:mpfq_p_4_src_vec_ur\b), const void *, ],
     [ (?^:mpfq_p_4_vec_ur\b), void *, ],
     [ (?^:mpfq_p_4_dst_vec_ur\b), void *, ],
     [ (?^:mpfq_p_4_poly \*), void *, ],
     [ (?^:mpfq_p_4_src_poly\b), const void *, ],
     [ (?^:mpfq_p_4_poly\b), void *, ],
     [ (?^:mpfq_p_4_dst_poly\b), void *, ],
     ],
    },
   vtag=p_4,
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
void mpfq_p_4_p_4_dotprod(mpfq_p_4_dst_field, mpfq_p_4_dst_field, mpfq_p_4_dst_vec, mpfq_p_4_src_vec, mpfq_p_4_src_vec, unsigned int);
void mpfq_p_4_p_4_addmul_tiny(mpfq_p_4_dst_field, mpfq_p_4_dst_field, mpfq_p_4_dst_vec, mpfq_p_4_src_vec, mpfq_p_4_dst_vec, unsigned int);
void mpfq_p_4_p_4_transpose(mpfq_p_4_dst_field, mpfq_p_4_dst_field, mpfq_p_4_dst_vec, mpfq_p_4_src_vec);

#endif  /* MPFQ_P_4_T_H_ */

/* vim:set ft=cpp: */
