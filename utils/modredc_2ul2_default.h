#include "modredc_2ul2.h"

/*
   Here are typedef's that rename all functions to mod_* instead of 
   modredc2ul2_*, which one day might become an automatic renaming scheme so 
   that different modulus sizes can be used simply by #including different 
   mod*.h files, but without changing anything else in the source code. 
 */

#ifdef  mod_init
#warning "modredc_2ul2_default.h included after another mod*_default.h ; bindings overwritten"
#endif

/* Define the default names and types for arithmetic with these functions */
#undef residue_t
#undef modulus_t
#undef modint_t
#undef MOD_SIZE
#undef MOD_MINBITS
#undef MOD_MAXBITS
#undef mod_intset
#undef mod_intset_ul
#undef mod_intequal
#undef mod_intequal_ul
#undef mod_intcmp
#undef mod_intcmp_ul
#undef mod_intfits_ul
#undef mod_intadd
#undef mod_intsub
#undef mod_intbits
#undef mod_intdivexact
#undef mod_init
#undef mod_init_noset0
#undef mod_clear
#undef mod_set
#undef mod_set_ul
#undef mod_set_ul_reduced
#undef mod_set_uls
#undef mod_set_uls_reduced
#undef mod_swap
#undef mod_initmod_ul
#undef mod_initmod_uls
#undef mod_getmod_ul
#undef mod_getmod_uls
#undef mod_clearmod
#undef mod_get_ul
#undef mod_get_uls
#undef mod_equal
#undef mod_is0
#undef mod_is1
#undef mod_add
#undef mod_add_ul
#undef mod_sub
#undef mod_sub_ul
#undef mod_neg
#undef mod_mul
#undef mod_div2
#undef mod_div3
#undef mod_div7
#undef mod_div13
#undef mod_pow_ul
#undef mod_2pow_ul
#undef mod_pow_mp
#undef mod_2pow_mp
#undef mod_V_ul
#undef mod_V_mp
#undef mod_sprp
#undef mod_sprp2
#undef mod_gcd
#undef mod_inv
#undef mod_jacobi
#undef mod_set0
#undef mod_set1
#undef mod_next
#undef mod_finished

#define residue_t            residueredc2ul2_t
#define modulus_t            modulusredc2ul2_t
#define modint_t             modintredc2ul2_t
#define MOD_SIZE             MODREDC2UL2_SIZE
#define MOD_MINBITS          MODREDC2UL2_MINBITS
#define MOD_MAXBITS          MODREDC2UL2_MAXBITS
#define mod_intset           modredc2ul2_intset
#define mod_intset_ul        modredc2ul2_intset_ul
#define mod_intequal         modredc2ul2_intequal
#define mod_intequal_ul      modredc2ul2_intequal_ul
#define mod_intcmp           modredc2ul2_intcmp
#define mod_intcmp_ul        modredc2ul2_intcmp_ul
#define mod_intfits_ul       modredc2ul2_intfits_ul
#define mod_intadd           modredc2ul2_intadd
#define mod_intsub           modredc2ul2_intsub
#define mod_intbits          modredc2ul2_intbits
#define mod_intdivexact      modredc2ul2_intdivexact
#define mod_init             modredc2ul2_init
#define mod_init_noset0      modredc2ul2_init_noset0
#define mod_clear            modredc2ul2_clear
#define mod_set              modredc2ul2_set
#define mod_set_ul           modredc2ul2_set_ul
#define mod_set_ul_reduced   modredc2ul2_set_ul_reduced
#define mod_set_uls          modredc2ul2_set_uls
#define mod_set_uls_reduced  modredc2ul2_set_uls_reduced
#define mod_swap             modredc2ul2_swap
#define mod_initmod_uls      modredc2ul2_initmod_uls
#define mod_getmod_ul        modredc2ul2_getmod_ul
#define mod_getmod_uls       modredc2ul2_getmod_uls
#define mod_clearmod         modredc2ul2_clearmod
#define mod_get_ul           modredc2ul2_get_ul
#define mod_get_uls          modredc2ul2_get_uls
#define mod_equal            modredc2ul2_equal
#define mod_is0              modredc2ul2_is0
#define mod_is1              modredc2ul2_is1
#define mod_add              modredc2ul2_add
#define mod_add_ul           modredc2ul2_add_ul
#define mod_sub              modredc2ul2_sub
#define mod_sub_ul           modredc2ul2_sub_ul
#define mod_neg              modredc2ul2_neg
#define mod_mul              modredc2ul2_mul
#define mod_div2             modredc2ul2_div2
#define mod_div3             modredc2ul2_div3
#define mod_div5             modredc2ul2_div5
#define mod_div7             modredc2ul2_div7
#define mod_div13            modredc2ul2_div13
#define mod_pow_ul           modredc2ul2_pow_ul
#define mod_2pow_ul          modredc2ul2_2pow_ul
#define mod_pow_mp           modredc2ul2_pow_mp
#define mod_2pow_mp          modredc2ul2_2pow_mp
#define mod_V_ul             modredc2ul2_V_ul
#define mod_V_mp             modredc2ul2_V_mp
#define mod_sprp             modredc2ul2_sprp
#define mod_sprp2            modredc2ul2_sprp2
#define mod_gcd              modredc2ul2_gcd
#define mod_inv              modredc2ul2_inv
#define mod_jacobi           modredc2ul2_jacobi
#define mod_set0             modredc2ul2_set0
#define mod_set1             modredc2ul2_set1
#define mod_next             modredc2ul2_next
#define mod_finished         modredc2ul2_finished
