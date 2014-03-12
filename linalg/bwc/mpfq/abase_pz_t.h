#ifndef ABASE_PZ_T_H_
#define ABASE_PZ_T_H_

/* MPFQ generated file -- do not edit */

#include "abase_pz.h"
/* Active handler: simd_pz */
/* Automatically generated code  */
/* Active handler: pz */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=pz, }, ],
   fieldtype=prime,
   tag=pz,
   type=plain,
   vbase_stuff={
    choose_byfeatures=<code>,
    families=[
     [ u64k1, u64k2, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_1, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_2, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_3, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_4, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_8, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=pz, }, ],
     ],
    member_templates_restrict={
     p_1=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_1, }, ],
     p_2=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_2, }, ],
     p_3=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_3, }, ],
     p_4=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_4, }, ],
     p_8=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_8, }, ],
     pz=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=pz, }, ],
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
void abase_pz_pz_dotprod(abase_pz_dst_field, abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec, abase_pz_src_vec, unsigned int);
void abase_pz_pz_addmul_tiny(abase_pz_dst_field, abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec, abase_pz_dst_vec, unsigned int);
void abase_pz_pz_transpose(abase_pz_dst_field, abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec);

#endif  /* ABASE_PZ_T_H_ */

/* vim:set ft=cpp: */
