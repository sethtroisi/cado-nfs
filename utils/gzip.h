#ifndef GZIP_H_
#define GZIP_H_

#ifdef __cplusplus
extern "C" {
#endif

extern int is_gzip(const char * s);
extern FILE *gzip_open(const char * s, const char * t);
extern void gzip_close(FILE *ofile, const char * s);

#ifdef __cplusplus
}
#endif

#endif	/* GZIP_H_ */
