#ifndef ABASE_U64K1_T_H_
#define ABASE_U64K1_T_H_

/* MPFQ generated file -- do not edit */

#include "abase_u64k1.h"
#include "abase_u64k2.h"
/* Active handler: simd_u64k */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Active handler: simd_dotprod */
/* Active handler: io */
/* Active handler: trivialities */
/* Active handler: simd_char2 */
/* Options used:{
   tag=u64k1,
   virtual_base={
    filebase=abase_vbase,
    name=abase_vbase,
    global_prefix=abase_,
    substitutions=[
     [ (?^:abase_u64k1_elt \*), void *, ],
     [ (?^:abase_u64k1_src_elt\b), const void *, ],
     [ (?^:abase_u64k1_elt\b), void *, ],
     [ (?^:abase_u64k1_dst_elt\b), void *, ],
     [ (?^:abase_u64k1_elt_ur \*), void *, ],
     [ (?^:abase_u64k1_src_elt_ur\b), const void *, ],
     [ (?^:abase_u64k1_elt_ur\b), void *, ],
     [ (?^:abase_u64k1_dst_elt_ur\b), void *, ],
     [ (?^:abase_u64k1_vec \*), void *, ],
     [ (?^:abase_u64k1_src_vec\b), const void *, ],
     [ (?^:abase_u64k1_vec\b), void *, ],
     [ (?^:abase_u64k1_dst_vec\b), void *, ],
     [ (?^:abase_u64k1_vec_ur \*), void *, ],
     [ (?^:abase_u64k1_src_vec_ur\b), const void *, ],
     [ (?^:abase_u64k1_vec_ur\b), void *, ],
     [ (?^:abase_u64k1_dst_vec_ur\b), void *, ],
     [ (?^:abase_u64k1_poly \*), void *, ],
     [ (?^:abase_u64k1_src_poly\b), const void *, ],
     [ (?^:abase_u64k1_poly\b), void *, ],
     [ (?^:abase_u64k1_dst_poly\b), void *, ],
     ],
    },
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
   family=[ u64k1, u64k2, ],
   w=64,
   k=1,
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
void abase_u64k1_u64k1_dotprod(abase_u64k1_dst_field, abase_u64k1_dst_field, abase_u64k1_dst_vec, abase_u64k1_src_vec, abase_u64k1_src_vec, unsigned int);
void abase_u64k1_u64k2_dotprod(abase_u64k1_dst_field, abase_u64k2_dst_field, abase_u64k1_dst_vec, abase_u64k2_src_vec, abase_u64k1_src_vec, unsigned int);
void abase_u64k1_u64k1_addmul_tiny(abase_u64k1_dst_field, abase_u64k1_dst_field, abase_u64k1_dst_vec, abase_u64k1_src_vec, abase_u64k1_dst_vec, unsigned int);
void abase_u64k1_u64k2_addmul_tiny(abase_u64k1_dst_field, abase_u64k2_dst_field, abase_u64k2_dst_vec, abase_u64k1_src_vec, abase_u64k2_dst_vec, unsigned int);
void abase_u64k1_u64k1_transpose(abase_u64k1_dst_field, abase_u64k1_dst_field, abase_u64k1_dst_vec, abase_u64k1_src_vec);
void abase_u64k1_u64k2_transpose(abase_u64k1_dst_field, abase_u64k2_dst_field, abase_u64k1_dst_vec, abase_u64k2_src_vec);

#endif  /* ABASE_U64K1_T_H_ */

/* vim:set ft=cpp: */
