#include "modredc_ul.h"

/*
   Here are typedef's that rename all functions to mod_* instead of 
   modredcul_*, which one day might become an automatic renaming scheme so 
   that different modulus sizes can be used simply by #including different 
   mod*.h files, but without changing anything else in the source code. 
 */

#ifdef  mod_init
#warning "modredc_ul_default.h included after another mod*_default.h ; bindings overwritten"
#endif

/* Define the default names and types for arithmetic with these functions */
#undef residue_t
#undef modulus_t
#undef modint_t
#undef MOD_SIZE
#undef MOD_MAXBITS
#undef mod_intset
#undef mod_intset_ul
#undef mod_intequal
#undef mod_intequal_ul
#undef mod_intcmp
#undef mod_intcmp_ul
#undef mod_intfits_ul
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
#undef mod_pow_ul
#undef mod_2pow_ul
#undef mod_pow_mp
#undef mod_2pow_mp
#undef mod_V_ul
#undef mod_V_mp
#undef mod_sprp
#undef mod_sprp2
#undef mod_isprime
#undef mod_gcd
#undef mod_inv
#undef mod_jacobi
#undef mod_set0
#undef mod_set1
#undef mod_next
#undef mod_finished

#define residue_t            residueredcul_t
#define modulus_t            modulusredcul_t
#define modint_t             modintredcul_t
#define MOD_SIZE             MODREDCUL_SIZE
#define MOD_MAXBITS          MODREDCUL_MAXBITS
#define mod_intset           modredcul_intset
#define mod_intset_ul        modredcul_intset_ul
#define mod_intequal         modredcul_intequal
#define mod_intequal_ul      modredcul_intequal_ul
#define mod_intcmp           modredcul_intcmp
#define mod_intcmp_ul        modredcul_intcmp_ul
#define mod_intfits_ul       modredcul_intfits_ul
#define mod_intbits          modredcul_intbits
#define mod_intdivexact      modredcul_intdivexact
#define mod_init             modredcul_init
#define mod_init_noset0      modredcul_init_noset0
#define mod_clear            modredcul_clear
#define mod_set              modredcul_set
#define mod_set_ul           modredcul_set_ul
#define mod_set_ul_reduced   modredcul_set_ul_reduced
#define mod_set_uls          modredcul_set_uls
#define mod_set_uls_reduced  modredcul_set_uls_reduced
#define mod_swap             modredcul_swap
#define mod_initmod_ul       modredcul_initmod_ul
#define mod_initmod_uls      modredcul_initmod_uls
#define mod_getmod_ul        modredcul_getmod_ul
#define mod_getmod_uls       modredcul_getmod_uls
#define mod_clearmod         modredcul_clearmod
#define mod_get_ul           modredcul_get_ul
#define mod_get_uls          modredcul_get_uls
#define mod_equal            modredcul_equal
#define mod_is0              modredcul_is0
#define mod_is1              modredcul_is1
#define mod_add              modredcul_add
#define mod_add_ul           modredcul_add_ul
#define mod_sub              modredcul_sub
#define mod_sub_ul           modredcul_sub_ul
#define mod_neg              modredcul_neg
#define mod_mul              modredcul_mul
#define mod_div2             modredcul_div2
#define mod_div3             modredcul_div3
#define mod_div7             modredcul_div7
#define mod_pow_ul           modredcul_pow_ul
#define mod_2pow_ul          modredcul_2pow_ul
#define mod_pow_mp           modredcul_pow_mp
#define mod_2pow_mp          modredcul_2pow_mp
#define mod_V_ul             modredcul_V_ul
#define mod_V_mp             modredcul_V_mp
#define mod_sprp             modredcul_sprp
#define mod_sprp2            modredcul_sprp2
#define mod_isprime          modredcul_isprime
#define mod_gcd              modredcul_gcd
#define mod_inv              modredcul_inv
#define mod_jacobi           modredcul_jacobi
#define mod_set0             modredcul_set0
#define mod_set1             modredcul_set1
#define mod_next             modredcul_next
#define mod_finished         modredcul_finished

