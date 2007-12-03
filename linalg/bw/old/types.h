#ifndef TYPES_H_
#define TYPES_H_

#ifdef	__cplusplus
extern "C" {
#endif

/* For C++ code, this is best handled using boost/cstdint.hpp */

/*****************************************************************************/
/* - type32 et al. - */
/* For the time being, we cannot use ISO C types, since they are not yet
 * widely available */
#ifndef ISO_C99_AVAILABLE
#if     defined(i386) || \
        defined(__i386__) || \
        defined(alpha) || \
        defined(__alpha__) || \
        defined(sparc) || \
        defined(__sparc__) || \
	defined(__x86_64__)
typedef unsigned int    type32;
typedef   signed int   stype32;
typedef unsigned short  type16;
typedef   signed short stype16;
typedef unsigned char   type8;
typedef   signed char  stype8;
#else
# error "Try to determine your builtin type sizes, or use ISO_C99_AVAILABLE"
#endif	/* known archs */
/* Find a 64bit integer type */
#if defined(i386) || defined(__i386__)
#ifdef __GNUC__
typedef unsigned long long int    type64;
typedef   signed long long int   stype64;
#else	/* __GNUC__ */
# error "Does your compiler (presumably not gcc) support a 64bit integer type?"
#endif	/* __GNUC__ */
#elif defined(alpha) || defined(__alpha__) || defined(__x86_64__)
typedef unsigned long int    type64;
typedef   signed long int   stype64;
#elif (defined(__sparc) || defined(sparc) ) && \
	(defined(__SUNPRO_C) || defined(__GNUC__))
typedef unsigned long int    type64;
typedef   signed long int   stype64;
#elif defined(__x86_64__) || defined(__x86_64)
typedef unsigned long int    type64;
typedef   signed long int   stype64;
#else	/* i386 / alpha / sparc */
# error "Does your machine (not i386/alpha/sparc) have a 64bit integer type?"
#endif	/* i386 / alpha /sparc */
#else	/* ISO_C99_AVAILABLE */
#include <inttypes.h>
typedef uint64_t  type64;
typedef  int64_t stype64;
typedef uint32_t  type32;
typedef  int32_t stype32;
typedef uint16_t  type16;
typedef  int16_t stype16;
typedef uint8_t   type8;
typedef  int8_t  stype8;
#endif	/* ISO_C99_AVAILABLE */

#define UC64(x)  ((type64)(x))
#define SC64(x) ((stype64)(x))
#define UC32(x)  ((type32)(x))
#define SC32(x) ((stype32)(x))
#define UC16(x)  ((type16)(x))
#define SC16(x) ((stype16)(x))
#define UC8(x)    ((type8)(x))
#define SC8(x)   ((stype8)(x))

#endif	/* TYPES_H_ */

/*****************************************************************************/
/* coord_t */

typedef type32 coord_t;

#ifdef	__cplusplus
}
#endif

/* vim:set sw=8: */
