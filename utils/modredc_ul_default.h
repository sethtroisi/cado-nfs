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

#undef MOD_RENAME
#define MOD_RENAME(x) MODREDCUL_RENAME(x)
#include "mod_rename.h"

#undef residue_t
#undef modulus_t
#undef modint_t
#undef MOD_SIZE
#undef MOD_MINBITS
#undef MOD_MAXBITS
#define residue_t            residueredcul_t
#define modulus_t            modulusredcul_t
#define modint_t             modintredcul_t
#define MOD_SIZE             MODREDCUL_SIZE
#define MOD_MAXBITS          MODREDCUL_MAXBITS
