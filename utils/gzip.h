#ifndef GZIP_H_
#define GZIP_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

extern const char * path_basename(const char * path);
extern int is_supported_compression_format(const char * s);
extern FILE * fopen_compressed_r(const char * name, int* p_pipeflag, char const ** suf);
extern FILE * fopen_compressed_w(const char * name, int* p_pipeflag, char const ** suf);

/* for compatibility only */
FILE * gzip_open(const char * name, const char * mode);
void gzip_close(FILE * f, const char * name);

#ifdef __cplusplus
}
#endif

#endif	/* GZIP_H_ */
