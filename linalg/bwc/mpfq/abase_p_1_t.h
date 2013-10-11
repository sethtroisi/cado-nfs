#ifndef ABASE_P_1_T_H_
#define ABASE_P_1_T_H_

/* MPFQ generated file -- do not edit */

#include "abase_p_1.h"
/* Active handler: simd_gfp */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::gfp::elt */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Options used:{
   type=plain,
   virtual_base={
    filebase=abase_vbase,
    name=abase_vbase,
    global_prefix=abase_,
    substitutions=[
     [ (?^:abase_p_1_elt \*), void *, ],
     [ (?^:abase_p_1_src_elt\b), const void *, ],
     [ (?^:abase_p_1_elt\b), void *, ],
     [ (?^:abase_p_1_dst_elt\b), void *, ],
     [ (?^:abase_p_1_elt_ur \*), void *, ],
     [ (?^:abase_p_1_src_elt_ur\b), const void *, ],
     [ (?^:abase_p_1_elt_ur\b), void *, ],
     [ (?^:abase_p_1_dst_elt_ur\b), void *, ],
     [ (?^:abase_p_1_vec \*), void *, ],
     [ (?^:abase_p_1_src_vec\b), const void *, ],
     [ (?^:abase_p_1_vec\b), void *, ],
     [ (?^:abase_p_1_dst_vec\b), void *, ],
     [ (?^:abase_p_1_vec_ur \*), void *, ],
     [ (?^:abase_p_1_src_vec_ur\b), const void *, ],
     [ (?^:abase_p_1_vec_ur\b), void *, ],
     [ (?^:abase_p_1_dst_vec_ur\b), void *, ],
     [ (?^:abase_p_1_poly \*), void *, ],
     [ (?^:abase_p_1_src_poly\b), const void *, ],
     [ (?^:abase_p_1_poly\b), void *, ],
     [ (?^:abase_p_1_dst_poly\b), void *, ],
     ],
    },
   w=64,
   fieldtype=prime,
   nn=3,
   tag=p_1,
   vbase_stuff={
    choose_byfeatures=<code>,
    vc:includes=[ <stdarg.h>, ],
    families=[
     [ { tag=p_4, cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, }, ],
     [ { tag=p_3, cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_1, }, ],
     [ u64k1, u64k2, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_8, }, ],
     [ { tag=p_2, cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, }, ],
     ],
    member_templates_restrict={
     p_8=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_8, }, ],
     p_2=[ { tag=p_2, cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, }, ],
     u64k2=[ u64k1, u64k2, ],
     p_3=[ { tag=p_3, cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, }, ],
     p_1=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_1, }, ],
     u64k1=[ u64k1, u64k2, ],
     p_4=[ { tag=p_4, cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, }, ],
     },
    },
   opthw=,
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_1, }, ],
   n=1,
   vtag=p_1,
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
void abase_p_1_p_1_dotprod(abase_p_1_dst_field, abase_p_1_dst_field, abase_p_1_dst_vec, abase_p_1_src_vec, abase_p_1_src_vec, unsigned int);
void abase_p_1_p_1_addmul_tiny(abase_p_1_dst_field, abase_p_1_dst_field, abase_p_1_dst_vec, abase_p_1_src_vec, abase_p_1_dst_vec, unsigned int);
void abase_p_1_p_1_transpose(abase_p_1_dst_field, abase_p_1_dst_field, abase_p_1_dst_vec, abase_p_1_src_vec);

#endif  /* ABASE_P_1_T_H_ */

/* vim:set ft=cpp: */
