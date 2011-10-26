#ifndef LAS_COORDINATES_H_
#define LAS_COORDINATES_H_

#include <stdint.h>
#include "las-types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  Forward declarations of conversion functions */
void xToIJ(int *i, unsigned int *j, const unsigned int X, sieve_info_srcptr si);
void NxToIJ(int *i, unsigned int *j, const unsigned int N, const unsigned int x, sieve_info_srcptr si);

void IJTox(unsigned int * x, int i, unsigned int j, sieve_info_srcptr si);
void IJToNx(unsigned int *N, unsigned int * x, int i, unsigned int j, sieve_info_srcptr si);
void IJToAB(int64_t *a, uint64_t *b, const int i, const unsigned int j, sieve_info_srcptr si);
void xToAB(int64_t *a, uint64_t *b, const unsigned int x, sieve_info_srcptr si);
void NxToAB(int64_t *a, uint64_t *b, const unsigned int N, const unsigned int x, sieve_info_srcptr si);
int ABToIJ(int *i, unsigned int *j, const int64_t a, const uint64_t b, sieve_info_srcptr si);
int ABTox(unsigned int *x, const int64_t a, const uint64_t b, sieve_info_srcptr si);
int ABToNx(unsigned int * N, unsigned int *x, const int64_t a, const uint64_t b, sieve_info_srcptr si);
/*  */


#ifdef __cplusplus
}
#endif

#endif	/* LAS_COORDINATES_H_ */
