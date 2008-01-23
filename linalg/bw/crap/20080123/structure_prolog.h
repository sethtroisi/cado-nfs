#ifndef STRUCTURE_PROLOG_H_
#define STRUCTURE_PROLOG_H_

#ifndef NO_STRICT_TYPECHECKING
#define STRICT_TYPECHECKING
#endif

#define PREFER_INLINES

/* That's really hackish. <<extern>> is mandatory here for gcc, as well
 * as static is for deccc.
 */

#ifdef __GNUC__
# define FASTFUNC extern __inline__
# define INCLUDING_STRUCTURE_INLINES_H_
#elif defined(__DECC)
# define COMPILER_KNOWS_PRAGMA
# define FASTFUNC static
# define INCLUDING_STRUCTURE_INLINES_H_
#else
# ifdef STRICT_ANSI
#  define FASTFUNC
# else
#  define FASTFUNC extern inline
# endif
#endif

/* There are some tricky declarations to ensure that the compiler will do
 * a good work of typechecking if we instruct him to do so. The code
 * looks a bit dirty. There is no influence on the machine code
 * generated. */

#ifdef STRICT_TYPECHECKING
#define STRICTTYPE_DECL(orig, alias)	\
struct alias ## _s { orig _x; };		\
typedef struct alias ## _s alias
#define STRICTTYPE_VAL(x)	((x)._x)
#ifdef __GNUC__
#define STRICTTYPE_CAST(name,x)	((name) { _x : x})
#else
#define STRICTTYPE_CAST(name,x)	((name) { x})
#endif
#else
#define STRICTTYPE_DECL(orig, alias)	typedef orig alias
#define STRICTTYPE_VAL(x)	x
#define STRICTTYPE_CAST(name,x)	((name) x)
#endif

/* XXX Duh, does this get used at any place ??? */

#if 0
STRICTTYPE_DECL(bw_scalar, bw_lvector);
STRICTTYPE_DECL(bw_scalar, bw_rvector);

#define bw_lvector_alloc(x)	bw_vector_alloc(STRICTTYPE_VAL(x),	\
					m_param,BW_SCALAR_SHORT)
#define bw_rvector_alloc(y)	bw_vector_alloc(STRICTTYPE_VAL(y),	\
					n_param,BW_SCALAR_SHORT)
#define bw_lvector_free(x)	free(STRICTTYPE_VAL(x))
#define bw_rvector_free(y)	free(STRICTTYPE_VAL(y))
#endif

#include "structure_automatic.h"

#ifdef INCLUDING_STRUCTURE_INLINES_H_
# include "structure_inlines.h"
#endif

#endif /* STRUCTURE_PROLOG_H_ */
