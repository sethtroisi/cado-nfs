#ifndef GZIP_H_
#define GZIP_H_

#include <stdio.h>
#include <libgen.h>
#include "prempt.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Put the directory of cado in rep_cado */
void set_rep_cado(const char *argv0, char *rep_cado);

/* Search the executable in PATH and, if found, return in real_path the
   complete path WITH the executable in the end */
char * search_real_exec_in_path(const char *executable, char *real_path);

/* Search the path of antebuffer and put the complete path + name of
   the executable ("antebuffer" or, if not found, cat) in antebuffer.
   Needs in rep_cado the cado directory, and, if it's given, in path_antebuffer
   the argument -path_antebuffer */
void search_antebuffer (const char *rep_cado, const char *path_antebuffer, char *antebuffer);

/* There are of course scores of existing basename() codes accessible,
 * starting with POSIX basename. However we fear posible inconsistencies
 * here, so we stick to a simple-and-stupid version, whose specification
 * meets our needs. Here, the returned string is always a substring of
 * the input string, and the latter never undergoes any modification.
 */
extern const char * path_basename(const char * path);

extern int is_supported_compression_format(const char * s);

/* Takes a filename, possibly ending with any recognized compression
 * extension, and returns the associated file stream. The stream may have
 * been opened by either fopen() of popen(), depending on whether an
 * external decompression program has to be called. If the caller is
 * intersted in knowing, the integer *p_pipeflag is filled with 1 (for
 * popen) or 0 (for fopen). In fact, the caller should care, because this
 * can be used to decide whether to close the stream with pclose or
 * fclose (even if fclose works, it's almost guaranteed to create
 * zombies).
 * If non-NULL, suf is the location of a pointer which is position to the
 * recognized suffix, which has been used to decide on which compression
 * method.
 */
extern FILE * fopen_maybe_compressed2(const char * name, const char * mode, int* p_pipeflag, char const ** suf);
extern FILE * fopen_maybe_compressed(const char * name, const char * mode);

/* This one just looks at the file name, and guesses again whether popen() or
 * fopen() was used. The file stream is then closed with pclose() or
 * fclose() accordingly.  */
extern void fclose_maybe_compressed(FILE *, const char * name);

#ifdef __cplusplus
}
#endif

#endif	/* GZIP_H_ */
