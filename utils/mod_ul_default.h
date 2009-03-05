#include "mod_ul.h"

/*
   Here are typedef's that rename all functions to mod_* instead of 
   modredcul_*, which one day might become an automatic renaming scheme so 
   that different modulus sizes can be used simply by #including different 
   mod*.h files, but without changing anything else in the source code. 
 */

#ifdef  mod_init
#warning "mod_ul_default.h included after another mod*_default.h ; bindings overwritten"
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
#undef mod_pow_mp
#undef mod_2pow_mp
#undef mod_V_ul
#undef mod_V_mp
#undef mod_sprp
#undef mod_gcd
#undef mod_inv
#undef mod_jacobi
#undef mod_set0
#undef mod_set1
#undef mod_next
#undef mod_finished

#define residue_t            residueul_t
#define modulus_t            modulusul_t
#define modint_t             modintul_t
#define MOD_SIZE             MODUL_SIZE
#define MOD_MAXBITS          MODUL_MAXBITS
#define mod_intset           modul_intset
#define mod_intset_ul        modul_intset_ul
#define mod_intequal         modul_intequal
#define mod_intequal_ul      modul_intequal_ul
#define mod_intcmp           modul_intcmp
#define mod_intcmp_ul        modul_intcmp_ul
#define mod_intfits_ul       modul_intfits_ul
#define mod_intbits          modul_intbits
#define mod_intdivexact      modul_intdivexact
#define mod_init             modul_init
#define mod_init_noset0      modul_init_noset0
#define mod_clear            modul_clear
#define mod_set              modul_set
#define mod_set_ul           modul_set_ul
#define mod_set_ul_reduced   modul_set_ul_reduced
#define mod_set_uls          modul_set_uls
#define mod_set_uls_reduced  modul_set_uls_reduced
#define mod_swap             modul_swap
#define mod_initmod_ul       modul_initmod_ul
#define mod_initmod_uls      modul_initmod_uls
#define mod_getmod_ul        modul_getmod_ul
#define mod_getmod_uls       modul_getmod_uls
#define mod_clearmod         modul_clearmod
#define mod_get_ul           modul_get_ul
#define mod_get_uls          modul_get_uls
#define mod_equal            modul_equal
#define mod_is0              modul_is0
#define mod_is1              modul_is1
#define mod_add              modul_add
#define mod_add_ul           modul_add_ul
#define mod_sub              modul_sub
#define mod_sub_ul           modul_sub_ul
#define mod_neg              modul_neg
#define mod_mul              modul_mul
#define mod_div2             modul_div2
#define mod_div3             modul_div3
#define mod_div7             modul_div7
#define mod_pow_ul           modul_pow_ul
#define mod_pow_mp           modul_pow_mp
#define mod_2pow_mp          modul_2pow_mp
#define mod_V_ul             modul_V_ul
#define mod_V_mp             modul_V_mp
#define mod_sprp             modul_sprp
#define mod_gcd              modul_gcd
#define mod_inv              modul_inv
#define mod_jacobi           modul_jacobi
#define mod_set0             modul_set0
#define mod_set1             modul_set1
#define mod_next             modul_next
#define mod_finished         modul_finished

