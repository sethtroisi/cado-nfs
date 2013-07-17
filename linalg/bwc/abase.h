#ifndef ABASE_H_
#define ABASE_H_

#if defined(SELECT_ABASE_u64k1) || defined(SELECT_ABASE_u64)
#include "mpfq/abase_u64k1.h"
#elif defined(SELECT_ABASE_u64k2)
#include "mpfq/abase_u64k2.h"
#elif defined(SELECT_ABASE_u64k3)
#include "mpfq/abase_u64k3.h"
#elif defined(SELECT_ABASE_u64k4)
#include "mpfq/abase_u64k4.h"
#elif defined(SELECT_ABASE_u64n)
#error "argh"
#include "mpfq/abase_u64n.h"
#elif defined(SELECT_ABASE_u128)
#error "argh"
#include "mpfq/abase_u128.h"
#elif defined(SELECT_ABASE_p16) /* This is really the first non-gf2 try */
#define NOT_OVER_GF2
#include "mpfq/abase_p16.h"

/* In reality we don't have everything configured for the moment. Update
 * CMakeLists.txt when the need for other moduli arises */

#elif defined(SELECT_ABASE_p_1)
#define NOT_OVER_GF2
#include "mpfq/abase_p_1.h"
#elif defined(SELECT_ABASE_p_2)
#define NOT_OVER_GF2
#include "mpfq/abase_p_2.h"
#elif defined(SELECT_ABASE_p_3)
#define NOT_OVER_GF2
#include "mpfq/abase_p_3.h"
#elif defined(SELECT_ABASE_p_4)
#define NOT_OVER_GF2
#include "mpfq/abase_p_4.h"
#elif defined(SELECT_ABASE_p_5)
#define NOT_OVER_GF2
#include "mpfq/abase_p_5.h"
#elif defined(SELECT_ABASE_p_6)
#define NOT_OVER_GF2
#include "mpfq/abase_p_6.h"
#elif defined(SELECT_ABASE_p_7)
#define NOT_OVER_GF2
#include "mpfq/abase_p_7.h"
#elif defined(SELECT_ABASE_p_8)
#define NOT_OVER_GF2
#include "mpfq/abase_p_8.h"
#elif defined(SELECT_ABASE_p_9)
#define NOT_OVER_GF2
#include "mpfq/abase_p_9.h"
#else
// #warning "Using default selection for abase"
#error "argh. This code must be compiled with some SELECT_ABASE_ macro defined"
// #include "mpfq/abase_u64.h"
#endif

/* This is used as a shorthand throughout in order to ease the access to
 * the _primary_ abase. Other abases have to be accessed via the OO
 * interface.
 */
#include "mpfq/mpfq_name_ab.h"

#endif	/* ABASE_H_ */
