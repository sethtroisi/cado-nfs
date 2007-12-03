#ifndef ENDIAN_H_
#define ENDIAN_H_

#ifdef	__cplusplus
extern "C" {
#endif

/*****************************************************************************/
/* Endianness concerns */

#if defined(i386) || defined(__i386__) || defined(alpha) || defined(__alpha__) || defined(__alpha) || defined(ia64) || defined(__ia64__)
#define KNOWN_LITTLE_ENDIAN
#endif

extern type32 builtin_endianness;
#define LITTLE_ENDIAN_RESULT    0x03020100U
#define BIG_ENDIAN_RESULT       0x00010203U
/* We don't include PDP-11 definition, as this would probably break at
 * many other places. */

#ifdef KNOWN_LITTLE_ENDIAN
#define M_LITTLE_ENDIAN (1)
#define M_BIG_ENDIAN    (0)
#define DO_LITTLE_ENDIAN(i)	{ i; }
#define DO_BIG_ENDIAN(i)
#else	/* KNOWN_LITTLE_ENDIAN */
#define M_LITTLE_ENDIAN (builtin_endianness==LITTLE_ENDIAN_RESULT)
#define M_BIG_ENDIAN    (builtin_endianness==BIG_ENDIAN_RESULT)
#define DO_LITTLE_ENDIAN(i)	if (M_LITTLE_ENDIAN) { i; }
#define DO_BIG_ENDIAN(i)	if (M_BIG_ENDIAN) { i; }
#endif	/* KNOWN_LITTLE_ENDIAN */

#define mswap32(x)      x=((x)>>24)|((x)<<24)|(((x)>>8)&UC32(0x0000FF00U))|(((x)<<8)&UC32(0x00FF0000))
#define mswap16(x)      x=((x)>>8)|((x)<<8)

void reverse(char *, size_t, size_t);
void check_endianness(FILE *);

#ifdef KNOWN_LITTLE_ENDIAN
#define WRITE32(d,s) *(type32*)(d)=(type32)(s)
#define READ32(d,s)  d=*(type32*)(s)

#else	/* KNOWN_LITTLE_ENDIAN */
void f_WRITE32(type32*,type32);
type32 f_READ32(type32*);
#define WRITE32(d,s) f_WRITE32((type32*)d,(type32)s)
#define READ32(d,s) d=f_READ32((type32*)s)

#endif	/* KNOWN_LITTLE_ENDIAN */

#endif	/* ENDIAN_H_ */

#ifdef	__cplusplus
}
#endif

/* vim:set sw=8: */
