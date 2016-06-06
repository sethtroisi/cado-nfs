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

#undef MOD_RENAME
#define MOD_RENAME(x) MODREDC2UL2_RENAME(x)
#include "mod_rename.h"

#undef residue_t
#undef modulus_t
#undef modint_t
#undef MOD_SIZE
#undef MOD_MINBITS
#undef MOD_MAXBITS
#define residue_t            residueredc2ul2_t
#define modulus_t            modulusredc2ul2_t
#define modint_t             modintredc2ul2_t
#define MOD_SIZE             MODREDC2UL2_SIZE
#define MOD_MINBITS          MODREDC2UL2_MINBITS
#define MOD_MAXBITS          MODREDC2UL2_MAXBITS

/* This function is used in mod_2ul_common.h with different 
   implementations for 15ul and 2ul2. Maybe there's a prettier way */
#undef mod_divn
#define mod_divn modredc2ul2_divn
