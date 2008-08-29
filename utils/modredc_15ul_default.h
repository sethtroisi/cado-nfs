/*
   Here are typedef's that rename all functions to mod_* instead of 
   modredc15ul_*, which one day might become an automatic renaming scheme so 
   that different modulus sizes can be used simply by #including different 
   mod*.h files, but without changing anything else in the source code. 
 */

#ifdef  mod_init
#warning "modredc_15ul_default.h included after another mod*_default.h ; bindings overwritten"
#endif

/* Define the default names and types for arithmetic with these functions */
#undef residue_t
#undef modulus_t
#undef modint_t
#undef MODUL_SIZE
#undef MODUL_MAXBITS
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

#define residue_t            residueredc15ul_t
#define modulus_t            modulusredc15ul_t
#define modint_t             modintredc15ul_t
#define MODUL_SIZE           MODREDC15UL_SIZE
#define MODUL_MAXBITS        MODREDC15UL_MAXBITS
#define mod_intset           modredc15ul_intset
#define mod_intset_ul        modredc15ul_intset_ul
#define mod_intequal         modredc15ul_intequal
#define mod_intequal_ul      modredc15ul_intequal_ul
#define mod_intcmp           modredc15ul_intcmp
#define mod_intcmp_ul        modredc15ul_intcmp_ul
#define mod_intfits_ul       modredc15ul_intfits_ul
#define mod_intbits          modredc15ul_intbits
#define mod_intdivexact      modredc15ul_intdivexact
#define mod_init             modredc15ul_init
#define mod_init_noset0      modredc15ul_init_noset0
#define mod_clear            modredc15ul_clear
#define mod_set              modredc15ul_set
#define mod_set_ul           modredc15ul_set_ul
#define mod_set_ul_reduced   modredc15ul_set_ul_reduced
#define mod_set_uls          modredc15ul_set_uls
#define mod_set_uls_reduced  modredc15ul_set_uls_reduced
#define mod_swap             modredc15ul_swap
#define mod_initmod_ul       modredc15ul_initmod_ul
#define mod_initmod_uls      modredc15ul_initmod_uls
#define mod_getmod_ul        modredc15ul_getmod_ul
#define mod_getmod_uls       modredc15ul_getmod_uls
#define mod_clearmod         modredc15ul_clearmod
#define mod_get_ul           modredc15ul_get_ul
#define mod_get_uls          modredc15ul_get_uls
#define mod_equal            modredc15ul_equal
#define mod_is0              modredc15ul_is0
#define mod_is1              modredc15ul_is1
#define mod_add              modredc15ul_add
#define mod_add_ul           modredc15ul_add_ul
#define mod_sub              modredc15ul_sub
#define mod_sub_ul           modredc15ul_sub_ul
#define mod_neg              modredc15ul_neg
#define mod_mul              modredc15ul_mul
#define mod_div2             modredc15ul_div2
#define mod_div3             modredc15ul_div3
#define mod_div7             modredc15ul_div7
#define mod_pow_ul           modredc15ul_pow_ul
#define mod_pow_mp           modredc15ul_pow_mp
#define mod_2pow_mp          modredc15ul_2pow_mp
#define mod_V_ul             modredc15ul_V_ul
#define mod_V_mp             modredc15ul_V_mp
#define mod_sprp             modredc15ul_sprp
#define mod_gcd              modredc15ul_gcd
#define mod_inv              modredc15ul_inv
#define mod_jacobi           modredc15ul_jacobi
#define mod_set0             modredc15ul_set0
#define mod_set1             modredc15ul_set1
#define mod_next             modredc15ul_next
#define mod_finished         modredc15ul_finished
