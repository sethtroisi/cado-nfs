#ifndef AUXFUNCS_H_
#define AUXFUNCS_H_

#include <stdlib.h>
#include <stdio.h>

#ifdef	__cplusplus
extern "C" {
#endif

#if defined(__GNUC__)
#define DOES_NOT_RETURN __attribute__ ((noreturn))
#else
#define DOES_NOT_RETURN /**/
#endif

unsigned int ceil_log2(unsigned int d);
void die(const char *, int, ... ) DOES_NOT_RETURN;
void eternal_sleep(void);
void * _my_malloc(size_t,const char *,int);
#define mymalloc(s) _my_malloc(s,__FILE__,__LINE__)
// #define BUG()	_BUG(__FILE__,__LINE__)
// void _BUG(const char *, const int);
int exist(const char *);
void coredump_limit(int, int);
/* That's a bad place to put this... */
void check_limit_requirements(int);
int mkfname(char *, const char *, ...);
long fcopy(FILE *, FILE *, size_t);
/* If the machine lacks snprintf, provide a dummy workaround */
#ifdef HAS_NOT_SNPRINTF
#include <stdarg.h>
#include <stdio.h>
extern int vsnprintf(char *, size_t, const char *, va_list);
extern int snprintf(char *, size_t, const char *, ...);
#endif	/* HAS_NOT_SNPRINTF */

#ifdef	__cplusplus
}
#endif

#endif /* AUXFUNCS_H_ */
