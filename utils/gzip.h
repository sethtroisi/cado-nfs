#ifndef GZIP_H_
#define GZIP_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

extern const char * path_basename(const char * path);
extern int is_supported_compression_format(const char * s);

/* Takes a filename, possibly ending with any recognized compression
 * extension, and returns the associated file stream. The stream may have
 * been opened by either fopen() of popen(), depending on whether an
 * external decompression program has to be called. If the caller is
 * intersted in knowing, the integer *p_pipeflag is filled with 1 (for
 * popen) or 0 (for popen). In fact, the caller should care, because this
 * can be used to decide whether to close the stream with pclose or
 * fclose (even if fclose works, it's almost guaranteed to create
 * zombies).
 * If non-NULL, suf is the location of a pointer which is position to the
 * recognized suffix, which has been used to decide on which compression
 * method.
 */
extern FILE * fopen_compressed_r(const char * name, int* p_pipeflag, char const ** suf);

/* analogous */
extern FILE * fopen_compressed_w(const char * name, int* p_pipeflag, char const ** suf);

/* for compatibility only */
FILE * gzip_open(const char * name, const char * mode);
void gzip_close(FILE * f, const char * name);

#ifdef __cplusplus
}
#endif

#endif	/* GZIP_H_ */
