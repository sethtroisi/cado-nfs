#ifndef CADO_POPEN_H_
#define CADO_POPEN_H_

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_GETRUSAGE
#include <sys/resource.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif


FILE * cado_popen(const char * command, const char * mode);
#ifdef HAVE_GETRUSAGE
void cado_pclose2(FILE * stream, struct rusage * r);
#else
void cado_pclose2(FILE * stream, void * r);
#endif

static inline void cado_pclose(FILE * stream) { cado_pclose2(stream, NULL); }

#ifdef __cplusplus
}
#endif

#endif	/* CADO_POPEN_H_ */
