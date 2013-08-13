#ifndef CADO_POPEN_H_
#define CADO_POPEN_H_

#include <stdio.h>
#include <unistd.h>
#include <sys/resource.h>

#ifdef __cplusplus
extern "C" {
#endif


FILE * cado_popen(const char * command, const char * mode);
void cado_pclose2(FILE * stream, struct rusage * r);
static inline void cado_pclose(FILE * stream) { cado_pclose2(stream, NULL); }

#ifdef __cplusplus
}
#endif

#endif	/* CADO_POPEN_H_ */
