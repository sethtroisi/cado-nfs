#ifndef ABASE_P_8_T_H_
#define ABASE_P_8_T_H_

/* MPFQ generated file -- do not edit */

#include "abase_p_8.h"
/* Active handler: simd_gfp */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::gfp::elt */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Options used:{
   n=8,
   opthw=,
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_8, }, ],
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
   tag=p_8,
   vtag=p_8,
   fieldtype=prime,
   w=64,
   virtual_base={
    filebase=abase_vbase,
    name=abase_vbase,
    global_prefix=abase_,
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
   type=plain,
   nn=17,
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
void abase_p_8_p_8_dotprod(abase_p_8_dst_field, abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec, abase_p_8_src_vec, unsigned int);
void abase_p_8_p_8_addmul_tiny(abase_p_8_dst_field, abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec, abase_p_8_dst_vec, unsigned int);
void abase_p_8_p_8_transpose(abase_p_8_dst_field, abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec);

#endif  /* ABASE_P_8_T_H_ */

/* vim:set ft=cpp: */
