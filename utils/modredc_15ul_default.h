#include "modredc_15ul.h"

/*
   Here are typedef's that rename all functions to mod_* instead of 
   modredc15ul_*, which one day might become an automatic renaming scheme so 
   that different modulus sizes can be used simply by #including different 
   mod*.h files, but without changing anything else in the source code. 
 */

#ifdef  mod_init
#warning "modredc_15ul_default.h included after another mod*_default.h ; bindings overwritten"
#endif

#undef MOD_RENAME
#define MOD_RENAME(x) MODREDC15UL_RENAME(x)
#include "mod_rename.h"

#undef residue_t
#undef modulus_t
#undef modint_t
#undef MOD_SIZE
#undef MOD_MINBITS
#undef MOD_MAXBITS
#define residue_t            residueredc15ul_t
#define modulus_t            modulusredc15ul_t
#define modint_t             modintredc15ul_t
#define MOD_SIZE             MODREDC15UL_SIZE
#define MOD_MAXBITS          MODREDC15UL_MAXBITS

/* This function is used in mod_2ul_common.h with different 
   implementations for 15ul and 2ul2. Maybe there's a prettier way */
#undef mod_divn
#define mod_divn modredc15ul_divn
